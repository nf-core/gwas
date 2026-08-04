// Canonical nf-core/test-datasets fixture bundle used by executable pipeline tests.
class FIXTURES {

    static final String UPSTREAM = 'https://raw.githubusercontent.com/nf-core/test-datasets/gwas/'

    static String base(Object projectDir = null) {
        return UPSTREAM
    }
}
