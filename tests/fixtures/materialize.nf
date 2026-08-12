nextflow.enable.dsl = 2

// Test-only prelaunch workflow. The production preparation workflow owns the same conversions;
// this small composition exists because route-profile manifests are resolved before a pipeline
// workflow can run and therefore need their literal PLINK paths materialized in advance.
include { PLINK2_VCF as PLINK2_VCF_ALL  } from '../../modules/local/plink2/vcf/main'
include { PLINK2_VCF as PLINK2_VCF_CHR1 } from '../../modules/local/plink2/vcf/main'
include { PLINK2_MAKEBED                  } from '../../modules/local/plink2/makebed/main'

workflow MATERIALIZE_GWAS_TEST_FIXTURES {
    def source_vcf = file(params.fixture_source_vcf, checkIfExists: true)

    PLINK2_VCF_ALL(Channel.value([[id: 'example_all'], source_vcf]))
    PLINK2_VCF_CHR1(Channel.value([[id: 'example_chr1'], source_vcf]))
    PLINK2_MAKEBED(PLINK2_VCF_ALL.out.pgen)
}

workflow {
    MATERIALIZE_GWAS_TEST_FIXTURES()
}
