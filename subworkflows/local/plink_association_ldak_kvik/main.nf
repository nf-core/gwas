// Run LDAK-KVIK Step 1 once and Step 2 across PLINK 1 genotype shards for one focal analysis.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.
//
// LOCAL SUBMISSION NOTE:
// nf-core/gwas does not call this generic fit-to-association composition because its
// pipeline-local LDAK-KVIK route reuses one fitted prediction bundle across multiple
// analysis requests. This subworkflow is retained as an independent nf-core/modules
// submission candidate for callers that want Step 1 and Step 2 composed per request.
// Remove this pipeline-specific note from the upstream submission.

// SUBWORKFLOW: Consisting entirely of upstream-ready nf-core/modules
include { LDAK_THINCOMMON } from '../../../modules/local/ldak/thincommon/main'
include { LDAK_KVIKSTEP1  } from '../../../modules/local/ldak/kvikstep1/main'
include { LDAK_KVIKSTEP2  } from '../../../modules/local/ldak/kvikstep2/main'

workflow PLINK_ASSOCIATION_LDAK_KVIK {
    take:
    ch_step1_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], once per analysis
    ch_genotype_shards // channel: [ val(meta2), path(bed), path(bim), path(fam) ], one or more per analysis
    ch_pheno // channel: [ val(meta3), path(phenotype_file), val(is_binary) ], once per analysis
    ch_qcovar // channel: [ val(meta4), path(quant_covariates_file) ], use [] for the optional file
    ch_covar // channel: [ val(meta5), path(cat_covariates_file) ], use [] for the optional file
    ch_step1_extract_policy // channel: [ val(meta6), path(extract_file), val(subset_policy) ], extract is [] for all/thin_common
    ch_keep // channel: [ val(meta7), path(keep_file) ], use [] for the optional file

    main:
    def ch_step1_inputs = ch_step1_genotypes
        .map { meta, bed, bim, fam -> [meta.id, [meta, bed, bim, fam]] }
        .join(ch_pheno.map { meta, phenotype_file, is_binary -> [meta.id, [meta, phenotype_file, is_binary]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_qcovar.map { meta, quant_covariates_file -> [meta.id, [meta, quant_covariates_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_covar.map { meta, cat_covariates_file -> [meta.id, [meta, cat_covariates_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_step1_extract_policy.map { meta, extract_file, subset_policy -> [meta.id, [meta, extract_file, subset_policy]] }, failOnDuplicate: true, failOnMismatch: true)

    def ch_thin_common_inputs = ch_step1_inputs.filter { _analysis_id, _genotypes, _pheno, _qcovar, _covar, extract_policy -> extract_policy[2] == 'thin_common' }
    LDAK_THINCOMMON(ch_thin_common_inputs.map { _analysis_id, genotypes, _pheno, _qcovar, _covar, _extract_policy -> genotypes })

    def ch_direct_step1_inputs = ch_step1_inputs
        .filter { _analysis_id, _genotypes, _pheno, _qcovar, _covar, extract_policy -> extract_policy[2] in ['all', 'provided'] }
        .map { _analysis_id, genotypes, pheno, qcovar, covar, extract_policy -> [genotypes, pheno, qcovar, covar, [extract_policy[0], extract_policy[1]]] }

    def ch_thin_step1_inputs = ch_thin_common_inputs
        .map { analysis_id, genotypes, pheno, qcovar, covar, _extract_policy -> [analysis_id, genotypes, pheno, qcovar, covar] }
        .join(LDAK_THINCOMMON.out.predictors.map { meta, extract_file -> [meta.id, [meta, extract_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .map { _analysis_id, genotypes, pheno, qcovar, covar, extract -> [genotypes, pheno, qcovar, covar, extract] }

    def ch_step1 = ch_direct_step1_inputs
        .mix(ch_thin_step1_inputs)
        .multiMap { genotypes, pheno, qcovar, covar, extract ->
            genotypes: genotypes
            pheno: pheno
            qcovar: qcovar
            covar: covar
            extract: extract
        }

    LDAK_KVIKSTEP1(ch_step1.genotypes, ch_step1.pheno, ch_step1.qcovar, ch_step1.covar, ch_step1.extract)

    def ch_step2 = ch_genotype_shards
        .map { meta, bed, bim, fam -> [meta.id, [meta, bed, bim, fam]] }
        .groupTuple(by: 0)
        .join(LDAK_KVIKSTEP1.out.predictions.map { meta, root, loco_details, loco_prs -> [meta.id, [meta, root, loco_details, loco_prs]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_pheno.map { meta, phenotype_file, _is_binary -> [meta.id, [meta, phenotype_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_qcovar.map { meta, quant_covariates_file -> [meta.id, [meta, quant_covariates_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_covar.map { meta, cat_covariates_file -> [meta.id, [meta, cat_covariates_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_keep.map { meta, keep_file -> [meta.id, [meta, keep_file]] }, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { _analysis_id, genotype_shards, predictions, pheno, qcovar, covar, keep ->
            genotype_shards.collect { genotype_shard ->
                [[predictions[0], genotype_shard[1], genotype_shard[2], genotype_shard[3]], pheno, predictions, qcovar, covar, keep]
            }
        }
        .multiMap { genotypes, pheno, predictions, qcovar, covar, keep ->
            genotypes: genotypes
            pheno: pheno
            predictions: predictions
            qcovar: qcovar
            covar: covar
            keep: keep
        }

    LDAK_KVIKSTEP2(ch_step2.genotypes, ch_step2.pheno, ch_step2.predictions, ch_step2.qcovar, ch_step2.covar, ch_step2.keep)

    emit:
    results             = LDAK_KVIKSTEP2.out.results // channel: [ val(meta), path(assoc) ], once per Step 2 shard
    harmonisation_input = LDAK_KVIKSTEP2.out.harmonisation_input // channel: [ val(meta), path(tsv) ], once per Step 2 shard
    summaries           = LDAK_KVIKSTEP2.out.summaries // channel: [ val(meta), path(summaries) ], once per Step 2 shard
    pvalues             = LDAK_KVIKSTEP2.out.pvalues // channel: [ val(meta), path(pvalues) ], optional per shard
    predictions         = LDAK_KVIKSTEP1.out.predictions // channel: [ val(meta), path(root), path(loco_details), path(loco_prs) ]
    effects             = LDAK_KVIKSTEP1.out.effects // channel: [ val(meta), path(effects) ], optional per analysis
    progress            = LDAK_THINCOMMON.out.progress // channel: [ val(meta), path(progress) ], only for thin_common analyses
    logs                = LDAK_KVIKSTEP1.out.log.mix(LDAK_KVIKSTEP2.out.log) // channel: [ val(meta), path(log) ]
}
