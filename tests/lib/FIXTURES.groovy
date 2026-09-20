// Canonical nf-core/test-datasets fixture bundle used by executable pipeline tests.
//
// Every fixture-backed test runs against a materialized, checksum-verified bundle rooted at
// GWAS_TEST_FIXTURES, rather than reading the published branch over the network per case.
// tests/fixtures/materialize.sh builds that bundle from the published branch and prints its root;
// tests/fixtures/nf-test.sh does it for a focused run. GWAS_FIXTURE_SOURCE selects the canonical
// source materialize.sh reads from, without committing a machine-specific path.
class FIXTURES {

    // The canonical bundle's published address, and the same prefix the published manifests embed.
    // Inside the tests it is a substitution token rather than a fetch target: RELATIONAL rewrites it
    // to the materialized root, so a declared fixture path stays portable in the test sources while
    // every read stays local and checksum-verified.
    static final String UPSTREAM = 'https://raw.githubusercontent.com/nf-core/test-datasets/gwas/'
    // The 23 files published on the nf-core/test-datasets `gwas` branch, plus the bivariate options
    // document materialize.sh copies in from assets/examples/relational. Mirrors the canonical list in
    // tests/fixtures/materialize.sh; the nine PLINK derivatives resolve over the network only once the
    // follow-up nf-core/test-datasets PR adding them has merged.
    private static final List<String> REQUIRED = [
        'results/fixtures/genotypes/example_all.vcf.gz',
        'results/fixtures/genotypes/example_all.pgen',
        'results/fixtures/genotypes/example_all.psam',
        'results/fixtures/genotypes/example_all.pvar',
        'results/fixtures/genotypes/example_all.bed',
        'results/fixtures/genotypes/example_all.bim',
        'results/fixtures/genotypes/example_all.fam',
        'results/fixtures/genotypes/example_chr1.pgen',
        'results/fixtures/genotypes/example_chr1.psam',
        'results/fixtures/genotypes/example_chr1.pvar',
        'results/fixtures/pheno_cov/example.pheno',
        'results/fixtures/pheno_cov/example.qcovar',
        'results/fixtures/pheno_cov/example.catcovar',
        'results/fixtures/relational/cohort_manifest.csv',
        'results/fixtures/relational/analysis_manifest_quantitative.csv',
        'results/fixtures/relational/analysis_manifest_binary.csv',
        'results/fixtures/relational/analysis_manifest_association_only.csv',
        'results/fixtures/relational/analysis_manifest_heritability_only.csv',
        'results/fixtures/relational/analysis_manifest_heterogeneous.csv',
        'results/fixtures/relational/method_options_heterogeneous.json',
        'results/fixtures/relational/resources/gcta_grm_extract.txt',
        'results/fixtures/relational/resources/ldak_predictor_extract.txt',
        'results/fixtures/relational/resources/ldak_weights.txt',
        'results/fixtures/relational/method_options_heterogeneous_bivariate.json',
    ].asImmutable()

    static String base(Object projectDir = null) {
        def declared = System.getenv('GWAS_TEST_FIXTURES')
        if (!declared) {
            // The suite deliberately does not fall back to UPSTREAM per case: hundreds of cases each
            // resolving the published branch over the network is slow and non-hermetic, and when it
            // fails it fails minutes in on schema validation with nothing naming the cause. Fail here
            // instead, in seconds, and let materialize.sh fetch the branch once.
            throw new IllegalStateException(
                'GWAS_TEST_FIXTURES is not set, and the fixture-backed tests do not read the published ' +
                'branch per case. Materialize the bundle and export the root it prints:\n' +
                '  export GWAS_TEST_FIXTURES=$(tests/fixtures/materialize.sh --profile docker)\n' +
                'or launch a focused run through tests/fixtures/nf-test.sh, which does that for you. ' +
                'Set GWAS_FIXTURE_SOURCE to the canonical fixture source before materializing the bundle.')
        }

        def directory = new File(declared).absoluteFile
        validateMaterializedRoot(directory, declared)
        def resolved = withTrailingSlash(directory.absolutePath)
        return resolved
    }

    private static void validateMaterializedRoot(File root, String declared) {
        def manifest = new File(root, '.complete.sha256')
        if (!manifest.isFile()) {
            throw new IllegalStateException("GWAS_TEST_FIXTURES is set to '${declared}' but ${manifest} is not a file")
        }
        def manifestLines = manifest.readLines()
        def entries = manifestLines.collectEntries { line ->
            def fields = line.trim().split(/\s+/, 2)
            if (fields.size() != 2 || !(fields[0] ==~ /[0-9a-fA-F]{64}/)) {
                throw new IllegalStateException("Invalid fixture checksum entry in ${manifest}: '${line}'")
            }
            [(fields[1]): fields[0].toLowerCase()]
        }
        if (manifestLines.size() != REQUIRED.size() || entries.keySet() != REQUIRED.toSet()) {
            throw new IllegalStateException("Fixture checksum manifest ${manifest} does not describe the complete runtime bundle")
        }
        REQUIRED.each { relativePath ->
            def fixture = new File(root, relativePath)
            if (!fixture.isFile()) {
                throw new IllegalStateException("GWAS_TEST_FIXTURES is set to '${declared}' but ${fixture} is not a file")
            }
            if (sha256(fixture) != entries[relativePath]) {
                throw new IllegalStateException("Fixture checksum mismatch for ${fixture}")
            }
        }
    }

    private static String sha256(File file) {
        def digest = java.security.MessageDigest.getInstance('SHA-256')
        file.withInputStream { input ->
            byte[] buffer = new byte[8192]
            int count
            while ((count = input.read(buffer)) != -1) {
                digest.update(buffer, 0, count)
            }
        }
        return digest.digest().collect { value -> String.format('%02x', value & 0xff) }.join()
    }

    private static String withTrailingSlash(String path) {
        return path.endsWith('/') ? path : path + '/'
    }
}
