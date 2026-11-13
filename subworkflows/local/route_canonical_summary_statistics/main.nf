// Converge every summary-statistics origin at the GWASLab boundary. Internal association results carry producer
// metadata and a producer-specific GWASLab column mapping; every external row supplies its explicit GWASLab
// format, including a pre-harmonised table declared with the `gwaslab` format. The GWASLab table is the pipeline
// summary-statistics representation and is published once under its stable scientific identity.
//
// This is pipeline routing, producer identity and publication policy, not an nf-core/modules submission
// candidate. This controller holds no meta-analysis eligibility rule and no harmonisation policy beyond the
// mappings the validated relational contract already registers. Every
// constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions,
// and it reads no params, no workflow and no projectDir — the build-keyed GWASLab resources arrive resolved.

// MODULE: Local to the pipeline
include { GWASLAB_HARMONIZE               } from '../../../modules/local/gwaslab/harmonize/main'

// FUNCTION: Local to the pipeline
include { getAssociationColumnMappingJson } from '../validate_gwas_input/method_registry'
include { getInternalSummaryMetadata      } from '../validate_gwas_input/identity_helpers'

workflow ROUTE_CANONICAL_SUMMARY_STATISTICS {
    take:
    ch_association_results // channel: [ val(meta), path(sumstats) ], one per analysis per association method actually exercised, `meta.method` named by the producing route
    ch_external_summary_statistics // channel: [ val(meta), path(source) ], every declared external summary row with an explicit GWASLab source format
    gwaslab_references // value: map of genome build to the resolved GWASLab reference resources [fasta, fasta_index, rsid_vcf, rsid_vcf_index, strand_vcf, strand_vcf_index], each `[]` when unconfigured

    main:

    // Every internal and external source converges before GWASLab. The producer-specific internal mappings
    // remain explicit and unchanged; an external row supplies a named GWASLab format.
    def ch_harmonise_records = ch_association_results
        .map { meta, source ->
            def summary_meta = getInternalSummaryMetadata(meta, meta.method) + [
                method: meta.method,
                source_name: source.name,
                gwaslab_input_format: getAssociationColumnMappingJson(meta.method, meta.is_binary),
            ]
            [summary_meta, source]
        }
        .mix(
            ch_external_summary_statistics.map { meta, source -> [meta + [method: 'external', gwaslab_input_format: meta.source_format], source] }
        )

    // The optional GWASLab resources remain build-keyed pipeline parameters, resolved once by the spine and
    // passed in whole. This is independent from the request-owned LDSC/LDAK reference catalog and does not
    // infer a scientific analysis reference.
    def ch_harmonise_input = ch_harmonise_records.multiMap { meta, source ->
        def references = gwaslab_references[meta.build]
        sumstats: [meta, source, meta.gwaslab_input_format, meta.build]
        reference_fasta: [[id: meta.build], references.fasta, references.fasta_index]
        rsid_reference: [[id: meta.build], references.rsid_vcf, references.rsid_vcf_index]
        strand_reference: [[id: meta.build], references.strand_vcf, references.strand_vcf_index]
    }

    GWASLAB_HARMONIZE(
        ch_harmonise_input.sumstats,
        ch_harmonise_input.reference_fasta,
        ch_harmonise_input.rsid_reference,
        ch_harmonise_input.strand_reference,
    )

    emit:
    summary_statistics = GWASLAB_HARMONIZE.out.sumstats // channel: [ val(meta), path(gwaslab_summary_statistics) ], exactly one per summary_statistics_id
    harmonization_log  = GWASLAB_HARMONIZE.out.log // channel: [ val(meta), path(gwaslab_log) ], one-to-one with summary_statistics on the same task
}
