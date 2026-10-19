// Prepare IID-keyed trait/design CSVs, then fit one or more ordered native MPH GRMs.
// Callers with native-ready CSVs invoke MPH_REML directly.

include { PREPARE_MPH_INPUTS } from '../../../modules/local/prepare_mph_inputs/main'
include { MPH_REML           } from '../../../modules/local/mph/reml/main'

workflow FIT_MPH_REML {
    take:
    ch_inputs // channel: [ val(meta), path(grm_files), path(phenotypes), path(quant_covariates), path(cat_covariates), path(fam), val(trait_names), val(missing_tokens) ]

    main:
    def ch_preparation = ch_inputs.multiMap { meta, grm_files, phenotypes, quant_covariates, cat_covariates, fam, trait_names, missing_tokens ->
        tables: [meta, phenotypes, quant_covariates, cat_covariates, trait_names, missing_tokens]
        samples: [meta, grm_files.find { grm_file -> grm_file.name.endsWith('.grm.iid') }, fam]
    }
    PREPARE_MPH_INPUTS(ch_preparation.tables, ch_preparation.samples)

    // MPH fits exactly the covariate columns it is named, and ignores a covariate file given without them,
    // so the fit needs the encoded column names as a value. Treatment-coded names exist only once the
    // serializer has observed the factor levels, and the header of the file it wrote is where it states them.
    // The design and its names travel as one element, because the optional join below represents an absent
    // right-hand side with a single `null`.
    def ch_covariate_design = PREPARE_MPH_INPUTS.out.covariates
        .map { meta, covariate_csv -> [meta, covariate_csv, covariate_csv] }
        .splitCsv(elem: 2, limit: 1)
        .map { meta, covariate_csv, header -> [meta, [covariate_csv, header.tail()]] }

    def ch_fit = ch_inputs
        .map { meta, grm_files, _phenotypes, _quant_covariates, _cat_covariates, _fam, trait_names, _missing_tokens -> [meta, grm_files, trait_names] }
        .join(PREPARE_MPH_INPUTS.out.phenotype, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_covariate_design, failOnDuplicate: true, remainder: true)
        .multiMap { meta, grm_files, trait_names, phenotype_csv, covariate_design ->
            grm: [meta, grm_files]
            phenotype: [meta, phenotype_csv, trait_names]
            covariates: [meta, covariate_design ? covariate_design.first() : [], covariate_design ? covariate_design.last() : []]
        }
    MPH_REML(ch_fit.grm, ch_fit.phenotype, ch_fit.covariates)

    emit:
    phenotype            = PREPARE_MPH_INPUTS.out.phenotype // channel: [ val(meta), path(phenotype_csv) ]
    covariates           = PREPARE_MPH_INPUTS.out.covariates // channel: [ val(meta), path(covariate_csv) ], absent when none supplied
    variance_components  = MPH_REML.out.variance_components // channel: [ val(meta), path(variance_components) ]
    fixed_effects        = MPH_REML.out.fixed_effects // channel: [ val(meta), path(fixed_effects) ]
    iterations           = MPH_REML.out.iterations // channel: [ val(meta), path(iterations) ]
    projected_phenotypes = MPH_REML.out.projected_phenotypes // channel: [ val(meta), path(projected_phenotypes) ]
    correlations         = MPH_REML.out.correlations // channel: [ val(meta), path(correlations) ], multi-trait fits only
    native_log           = MPH_REML.out.log // channel: [ val(meta), path(log) ]
}
