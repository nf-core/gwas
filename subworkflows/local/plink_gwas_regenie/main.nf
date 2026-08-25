// Run REGENIE Step 1 fitting and Step 2 association across PLINK-format genotype shards.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.
//
// The pipeline uses a separate route that reuses one fitted prediction bundle across multiple
// analysis requests. This reusable composition instead fits and associates per request.

// SUBWORKFLOW: Consisting entirely of nf-core/modules
include { PLINK_FIT_REGENIE } from '../plink_fit_regenie/main'

// MODULE: Installed directly from nf-core/modules
include { REGENIE_STEP2     } from '../../../modules/nf-core/regenie/step2/main'

workflow PLINK_GWAS_REGENIE {
    take:
    ch_step1_genotypes // channel: [ val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file) ], Step 1 focal identity
    ch_step2_genotypes // channel: [ val(meta2), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file) ], Step 2 shards; public outputs retain the Step 1 focal identity
    ch_pheno // channel: [ val(meta3), path(pheno) ]
    ch_covar // channel: [ val(meta4), path(covar) ], use [] when absent
    ch_step1_bsize // channel: [ val(meta5), val(step1_bsize) ]
    ch_step2_bsize // channel: [ val(meta6), val(step2_bsize) ]
    ch_step1_mode // channel: [ val(meta7), val(step1_mode) ], 'standard' or 'chunked'
    ch_n_l0_jobs // channel: [ val(meta8), val(n_l0_jobs) ], use [] for standard mode

    main:
    PLINK_FIT_REGENIE(ch_step1_genotypes, ch_pheno, ch_covar, ch_step1_bsize, ch_step1_mode, ch_n_l0_jobs)

    def ch_step2 = ch_step2_genotypes
        .map { meta, plink_genotype_file, plink_variant_file, plink_sample_file -> [meta.id, [meta, plink_genotype_file, plink_variant_file, plink_sample_file]] }
        .groupTuple(by: 0)
        .join(PLINK_FIT_REGENIE.out.predictions.map { meta, predictions -> [meta.id, [meta, predictions]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(PLINK_FIT_REGENIE.out.loco.map { meta, loco -> [meta.id, [meta, loco]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_pheno.map { meta, pheno -> [meta.id, [meta, pheno]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_covar.map { meta, covar -> [meta.id, [meta, covar]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_step2_bsize.map { meta, step2_bsize -> [meta.id, step2_bsize] }, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { _analysis_id, genotype_shards, predictions, loco, pheno, covar, step2_bsize ->
            genotype_shards.collect { genotype_shard ->
                [[predictions[0], genotype_shard[1], genotype_shard[2], genotype_shard[3]], [predictions[0], predictions[1], loco[1]], pheno, covar, step2_bsize]
            }
        }
        .multiMap { genotypes, predictions, pheno, covar, step2_bsize ->
            genotypes: genotypes
            predictions: predictions
            pheno: pheno
            covar: covar
            bsize: step2_bsize
        }

    REGENIE_STEP2(ch_step2.genotypes, ch_step2.predictions, ch_step2.pheno, ch_step2.covar, ch_step2.bsize)

    emit:
    results     = REGENIE_STEP2.out.results // channel: [ val(meta), path(regenie_results) ]
    logs        = PLINK_FIT_REGENIE.out.logs.mix(REGENIE_STEP2.out.log) // channel: [ val(meta), path(log) ]
    predictions = PLINK_FIT_REGENIE.out.predictions // channel: [ val(meta), path(predictions) ]
    loco        = PLINK_FIT_REGENIE.out.loco // channel: [ val(meta), path(loco) ]
}
