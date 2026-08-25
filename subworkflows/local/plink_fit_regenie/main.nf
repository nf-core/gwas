// Fit REGENIE Step 1 prediction models in standard mode or through the split L0/L1 route.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Installed directly from nf-core/modules
include { REGENIE_STEP1   } from '../../../modules/nf-core/regenie/step1/main'
include { REGENIE_SPLITL0 } from '../../../modules/nf-core/regenie/splitl0/main'
include { REGENIE_RUNL0   } from '../../../modules/nf-core/regenie/runl0/main'
include { REGENIE_RUNL1   } from '../../../modules/nf-core/regenie/runl1/main'

workflow PLINK_FIT_REGENIE {
    take:
    ch_genotypes // channel: [ val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file) ]
    ch_pheno // channel: [ val(meta2), path(pheno) ]
    ch_covar // channel: [ val(meta3), path(covar) ], use [] when absent
    ch_bsize // channel: [ val(meta4), val(bsize) ]
    ch_step1_mode // channel: [ val(meta5), val(step1_mode) ], 'standard' or 'chunked'
    ch_n_l0_jobs // channel: [ val(meta6), val(n_l0_jobs) ], use [] for standard mode

    main:
    ch_modes = ch_step1_mode.map { meta, step1_mode ->
        if (!['standard', 'chunked'].contains(step1_mode)) {
            error("[nf-core/gwas] ERROR: PLINK_FIT_REGENIE step1_mode must be 'standard' or 'chunked', got '${step1_mode}'")
        }
        tuple(meta.id, step1_mode)
    }

    ch_step1_inputs = ch_genotypes
        .map { meta, plink_genotype_file, plink_variant_file, plink_sample_file -> tuple(meta.id, tuple(meta, plink_genotype_file, plink_variant_file, plink_sample_file)) }
        .join(ch_pheno.map { meta, pheno -> tuple(meta.id, tuple(meta, pheno)) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_covar.map { meta, covar -> tuple(meta.id, tuple(meta, covar)) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_bsize.map { meta, bsize -> tuple(meta.id, bsize) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_modes, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_n_l0_jobs.map { meta, n_l0_jobs -> tuple(meta.id, n_l0_jobs) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .branch { analysis_id, genotypes, pheno, covar, bsize, step1_mode, n_l0_jobs ->
            standard: step1_mode == 'standard'
            return tuple(analysis_id, genotypes, pheno, covar, bsize)
            chunked: step1_mode == 'chunked'
            return tuple(analysis_id, genotypes, pheno, covar, bsize, n_l0_jobs)
        }

    ch_standard = ch_step1_inputs.standard.multiMap { _analysis_id, genotypes, pheno, covar, bsize ->
        genotypes: genotypes
        pheno: pheno
        covar: covar
        bsize: bsize
    }
    REGENIE_STEP1(ch_standard.genotypes, ch_standard.pheno, ch_standard.covar, ch_standard.bsize)

    ch_chunked_inputs = ch_step1_inputs.chunked
    ch_splitl0 = ch_chunked_inputs.multiMap { _analysis_id, genotypes, pheno, covar, bsize, n_l0_jobs ->
        genotypes: genotypes
        pheno: pheno
        covar: covar
        bsize: bsize
        n_l0_jobs: n_l0_jobs
    }
    REGENIE_SPLITL0(ch_splitl0.genotypes, ch_splitl0.pheno, ch_splitl0.covar, ch_splitl0.bsize, ch_splitl0.n_l0_jobs)

    ch_split_jobs = REGENIE_SPLITL0.out.master
        .map { meta, master -> tuple(meta.id, master) }
        .join(REGENIE_SPLITL0.out.snplists.map { meta, snplists -> tuple(meta.id, snplists) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { analysis_id, master, snplists ->
            (snplists instanceof List ? snplists : [snplists]).collect { snplist ->
                def job_match = snplist.name =~ /_job(\d+)\.snplist$/
                tuple(analysis_id, master, snplist, job_match[0][1] as Integer)
            }
        }

    ch_runl0_inputs = ch_chunked_inputs
        .combine(ch_split_jobs, by: 0)
        .map { analysis_id, genotypes, pheno, covar, bsize, _n_l0_jobs, master, snplist, job_number ->
            tuple(analysis_id, job_number, genotypes, tuple(genotypes[0], master, snplist, job_number), pheno, covar, bsize)
        }
    ch_runl0 = ch_runl0_inputs.multiMap { _analysis_id, _job_number, genotypes, l0_inputs, pheno, covar, bsize ->
        genotypes: genotypes
        l0: l0_inputs
        pheno: pheno
        covar: covar
        bsize: bsize
    }
    REGENIE_RUNL0(ch_runl0.genotypes, ch_runl0.l0, ch_runl0.pheno, ch_runl0.covar, ch_runl0.bsize)

    ch_l0_bundle = ch_runl0_inputs
        .map { analysis_id, job_number, _genotypes, l0_inputs, _pheno, _covar, _bsize -> tuple(analysis_id, job_number, l0_inputs[1], l0_inputs[2]) }
        .join(
            REGENIE_RUNL0.out.l0_predictions.map { meta, l0_predictions ->
                def prediction_files = l0_predictions instanceof List ? l0_predictions : [l0_predictions]
                def job_match = prediction_files[0].name =~ /_job(\d+)_l0_Y/
                tuple(meta.id, job_match[0][1] as Integer, prediction_files)
            },
            by: [0, 1],
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .groupTuple(by: 0)
        .map { analysis_id, job_numbers, masters, snplists, l0_predictions ->
            def indices = (0..<job_numbers.size()).toList().sort { index -> job_numbers[index] }
            tuple(analysis_id, masters[indices[0]], indices.collect { index -> snplists[index] }, indices.collectMany { index -> l0_predictions[index] })
        }

    ch_runl1 = ch_chunked_inputs
        .join(ch_l0_bundle, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .multiMap { _analysis_id, genotypes, pheno, covar, bsize, _n_l0_jobs, master, snplists, l0_predictions ->
            genotypes: genotypes
            l1: tuple(genotypes[0], master, snplists, l0_predictions)
            pheno: pheno
            covar: covar
            bsize: bsize
        }
    REGENIE_RUNL1(ch_runl1.genotypes, ch_runl1.l1, ch_runl1.pheno, ch_runl1.covar, ch_runl1.bsize)

    def ch_logs = REGENIE_STEP1.out.log
    ch_logs = ch_logs.mix(REGENIE_SPLITL0.out.log)
    ch_logs = ch_logs.mix(REGENIE_RUNL0.out.log)
    ch_logs = ch_logs.mix(REGENIE_RUNL1.out.log)

    emit:
    predictions = REGENIE_STEP1.out.predictions.mix(REGENIE_RUNL1.out.predictions) // channel: [ val(meta), path(predictions) ]
    loco        = REGENIE_STEP1.out.loco.mix(REGENIE_RUNL1.out.loco) // channel: [ val(meta), path(loco) ]
    logs        = ch_logs // channel: [ val(meta), path(log) ]
}
