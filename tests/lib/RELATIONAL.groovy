// Linked-manifest builders for pipeline-level tests.
//
// Fixture paths name the portable nf-core/test-datasets GWAS bundle and are resolved through FIXTURES.base,
// which reads the materialized, checksum-verified bundle rooted at GWAS_TEST_FIXTURES.
class RELATIONAL {

    static String cohorts(Object projectDir, Object outputDir, String name, Closure mutate = null) {
        def header = ['cohort_id', 'genome_build', 'ancestry', 'pgen', 'psam', 'pvar', 'bed', 'bim', 'fam', 'genotype_view_id']
        def rows = [cohort('example_pgen')]
        if (mutate) {
            mutate(rows)
        }
        return materialise(projectDir, outputDir, name, 'cohorts', header, rows, ['pgen', 'psam', 'pvar', 'bed', 'bim', 'fam'])
    }

    static Map cohort(String cohortId) {
        def fixture = { path -> "${FIXTURES.UPSTREAM}results/fixtures/${path}" }
        def cohorts = [
            example_pgen: [
                cohort_id: 'example_pgen',
                genome_build: 'GRCh37',
                ancestry: 'EUR',
                pgen: fixture('genotypes/example_all.pgen'),
                psam: fixture('genotypes/example_all.psam'),
                pvar: fixture('genotypes/example_all.pvar'),
                bed: '',
                bim: '',
                fam: '',
                genotype_view_id: '',
            ],
            example_bfile: [
                cohort_id: 'example_bfile',
                genome_build: 'GRCh37',
                ancestry: 'EUR',
                pgen: '',
                psam: '',
                pvar: '',
                bed: fixture('genotypes/example_all.bed'),
                bim: fixture('genotypes/example_all.bim'),
                fam: fixture('genotypes/example_all.fam'),
                genotype_view_id: '',
            ],
        ]
        if (!cohorts.containsKey(cohortId)) {
            throw new IllegalArgumentException("No shipped cohort fixture '${cohortId}'")
        }
        return new LinkedHashMap(cohorts[cohortId])
    }

    static String appendColumn(String manifestPath, String column, Object value) {
        def manifest = new File(manifestPath)
        def lines = manifest.readLines()
        manifest.text = (["${lines.first()},${column}"] + lines.tail().collect { line -> "${line},${quote(value)}" }).join('\n') + '\n'
        return manifestPath
    }

    static String retainColumns(String manifestPath, List<String> columns) {
        def manifest = new File(manifestPath)
        def rows = manifest.readLines().collect { line -> line.split(',', -1).toList() }
        def indices = columns.collect { column -> rows.first().indexOf(column) }
        manifest.text = rows.collect { row -> indices.collect { index -> row[index] }.join(',') }.join('\n') + '\n'
        return manifestPath
    }

    static String quoteHeaders(String manifestPath) {
        def manifest = new File(manifestPath)
        def lines = manifest.readLines()
        manifest.text = ([lines.first().split(',', -1).collect { column -> quote(column) }.join(',')] + lines.tail()).join('\n') + '\n'
        return manifestPath
    }

    static String analyses(Object projectDir, Object outputDir, String name, Closure mutate = null) {
        def header = [
            'analysis_id',
            'cohort_id',
            'trait_id',
            'trait_type',
            'phenotype',
            'phenotype_column',
            'control_value',
            'case_value',
            'quant_covariates',
            'cat_covariates',
            'association_methods',
            'heritability_methods',
            'population_prevalence',
            'sample_prevalence',
        ]
        def rows = [analysis(projectDir, 'example_pgen_qt')]
        if (mutate) {
            mutate(rows)
        }
        return materialise(
            projectDir,
            outputDir,
            name,
            'analyses',
            header,
            rows,
            ['phenotype', 'quant_covariates', 'cat_covariates'],
        )
    }

    static Map analysis(Object projectDir, String analysisId) {
        def fixture = { path -> "${FIXTURES.UPSTREAM}results/fixtures/${path}" }
        def common = [
            cohort_id: 'example_pgen',
            phenotype: fixture('pheno_cov/example.pheno'),
            quant_covariates: fixture('pheno_cov/example.qcovar'),
            cat_covariates: fixture('pheno_cov/example.catcovar'),
            association_methods: '',
            heritability_methods: '',
            population_prevalence: '',
            sample_prevalence: '',
        ]
        def analyses = [
            example_pgen_qt: common + [
                analysis_id: 'example_pgen_qt',
                trait_id: 'QT',
                trait_type: 'quantitative',
                phenotype_column: 'QT',
                control_value: '',
                case_value: '',
            ],
            example_pgen_bt: common + [
                analysis_id: 'example_pgen_bt',
                trait_id: 'BT',
                trait_type: 'binary',
                phenotype_column: 'BT',
                control_value: '1',
                case_value: '2',
            ],
            example_bfile_qt: common + [
                analysis_id: 'example_bfile_qt',
                cohort_id: 'example_bfile',
                trait_id: 'QT',
                trait_type: 'quantitative',
                phenotype_column: 'QT',
                control_value: '',
                case_value: '',
            ],
        ]
        if (!analyses.containsKey(analysisId)) {
            throw new IllegalArgumentException("No shipped analysis fixture '${analysisId}'")
        }
        return new LinkedHashMap(analyses[analysisId])
    }

    // A byte-for-byte copy of a fixture file under a new name. `resource` writes text, which corrupts a
    // binary member such as a `.pgen`, so a test that needs a renamed-but-identical bundle uses this.
    static String copyResource(Object outputDir, String name, Object source) {
        def directory = new File(new File(outputDir.toString()).parentFile, 'resources')
        directory.mkdirs()
        def target = new File(directory, name)
        target.parentFile.mkdirs()
        target.bytes = new File(source.toString()).bytes
        return target.absolutePath
    }

    static String resource(Object outputDir, String name, String content) {
        def directory = new File(new File(outputDir.toString()).parentFile, 'resources')
        directory.mkdirs()
        def resource = new File(directory, name)
        resource.parentFile.mkdirs()
        resource.text = content
        return resource.absolutePath
    }

    static String methodOptions(Object outputDir, String name, Object options) {
        def directory = new File(new File(outputDir.toString()).parentFile, "manifests/${name}")
        directory.mkdirs()
        def document = new File(directory, 'method_options.json')
        document.text = options instanceof String
            ? options
            : groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(options)) + '\n'
        return document.absolutePath
    }

    private static String materialise(
        Object projectDir,
        Object outputDir,
        String name,
        String role,
        List<String> header,
        List<Map> sourceRows,
        List<String> fileColumns
    ) {
        def rows = sourceRows.collect { row -> new LinkedHashMap(row) }
        def fixtureBase = FIXTURES.base(projectDir)
        rows.each { row ->
            fileColumns.each { column ->
                if (row[column]) {
                    row[column] = row[column].toString().replace(FIXTURES.UPSTREAM, fixtureBase)
                }
            }
        }

        def directory = new File(new File(outputDir.toString()).parentFile, "manifests/${name}")
        directory.mkdirs()
        def manifest = new File(directory, "${role}.csv")
        def body = rows.collect { row -> header.collect { column -> quote(row[column]) }.join(',') }
        manifest.text = ([header.join(',')] + body).join('\n') + '\n'
        return manifest.absolutePath
    }

    private static String quote(Object value) {
        def text = value == null ? '' : value.toString()
        return text.contains(',') ? "\"${text}\"" : text
    }
}
