// Linked-manifest builders for pipeline-level tests.
//
// Fixture paths name the portable nf-core/test-datasets GWAS bundle and are resolved through FIXTURES.base,
// so executable tests consistently use the public nf-core/test-datasets fixture bundle.
class RELATIONAL {

    static String cohorts(Object projectDir, Object outputDir, String name, Closure mutate = null) {
        def header = ['cohort_id', 'genome_build', 'ancestry', 'pgen', 'psam', 'pvar', 'bed', 'bim', 'fam', 'vcf']
        def rows = [cohort('example_pgen')]
        if (mutate) {
            mutate(rows)
        }
        return materialise(projectDir, outputDir, name, 'cohorts', header, rows, ['pgen', 'psam', 'pvar', 'bed', 'bim', 'fam', 'vcf'])
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
                vcf: '',
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
                vcf: '',
            ],
            example_vcf: [
                cohort_id: 'example_vcf',
                genome_build: 'GRCh37',
                ancestry: 'EUR',
                pgen: '',
                psam: '',
                pvar: '',
                bed: '',
                bim: '',
                fam: '',
                vcf: fixture('genotypes/example_all.vcf.gz'),
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
            association_methods: 'plink2',
            heritability_methods: '',
            population_prevalence: '',
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
            example_vcf_qt: common + [
                analysis_id: 'example_vcf_qt',
                cohort_id: 'example_vcf',
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

    static Map twoChromosomePlink1(Object projectDir, Object outputDir, String name) {
        def source = cohort('example_bfile')
        def fixtureBase = FIXTURES.base(projectDir)
        def directory = new File(new File(outputDir.toString()).parentFile, "resources/${name}")
        directory.mkdirs()
        def sourcePath = { value -> value.toString().replace(FIXTURES.UPSTREAM, fixtureBase) }
        def sourceBytes = { value ->
            def path = sourcePath(value)
            return path.startsWith('http') ? new URL(path).bytes : new File(path).bytes
        }
        def sourceLines = { value ->
            def path = sourcePath(value)
            return path.startsWith('http') ? new URL(path).readLines() : new File(path).readLines()
        }
        def bed = new File(directory, 'example_two_chromosomes.bed')
        def bim = new File(directory, 'example_two_chromosomes.bim')
        def fam = new File(directory, 'example_two_chromosomes.fam')
        bed.bytes = sourceBytes(source.bed)
        fam.bytes = sourceBytes(source.fam)
        def bimLines = sourceLines(source.bim)
        def split = Math.max(1, (int) (bimLines.size() / 2))
        bim.text = bimLines.withIndex().collect { line, index ->
            def fields = line.split('\t', -1)
            fields[0] = index < split ? '1' : '2'
            fields[3] = ((index + 1) * 100).toString()
            return fields.join('\t')
        }.join('\n') + '\n'
        return [bed: bed.absolutePath, bim: bim.absolutePath, fam: fam.absolutePath]
    }

    static String predictors(Object projectDir, Object outputDir, String name, int count) {
        def source = cohort('example_pgen').pvar
        def fixture = source.toString().replace(FIXTURES.UPSTREAM, FIXTURES.base(projectDir))
        def lines = fixture.startsWith('http') ? new URL(fixture).readLines() : new File(fixture).readLines()
        def records = lines.findAll { line -> line && !line.startsWith('#') }
        def firstCount = Math.max(1, (int) (count / 2))
        def predictors = (records.take(firstCount) + records.takeRight(count - firstCount))
            .collect { line -> line.tokenize()[2] }
        return resource(outputDir, name, predictors.join('\n') + '\n')
    }

    static String weights(Object projectDir, Object outputDir, String name, int value) {
        def source = cohort('example_pgen').pvar
        def fixture = source.toString().replace(FIXTURES.UPSTREAM, FIXTURES.base(projectDir))
        def lines = fixture.startsWith('http') ? new URL(fixture).readLines() : new File(fixture).readLines()
        def content = lines
            .findAll { line -> line && !line.startsWith('#') }
            .collect { line -> "${line.tokenize()[2]} ${value}" }
            .join('\n') + '\n'
        return resource(outputDir, name, content)
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

    static String legacy35Header(Object outputDir, String name) {
        def header = [
            'analysis_id', 'cohort_id', 'trait_id', 'trait_type', 'genome_build', 'ancestry',
            'pgen', 'psam', 'pvar', 'bed', 'bim', 'fam', 'vcf', 'vcf_index',
            'phenotype', 'phenotype_column', 'control_value', 'case_value', 'quant_covariates',
            'cat_covariates', 'association_methods', 'heritability_methods', 'population_prevalence',
            'sample_prevalence', 'ldak_model', 'ldak_power', 'ldak_weights',
            'ldak_relatedness_filter', 'ldak_kvik_step1_subset', 'ldak_kvik_step1_extract',
            'gcta_grm_parts', 'gcta_sparse_cutoff', 'gcta_ld_score_region_kb', 'gcta_ld_bins',
            'gcta_ldms_maf_edges',
        ]
        def directory = new File(new File(outputDir.toString()).parentFile, "manifests/${name}")
        directory.mkdirs()
        def manifest = new File(directory, 'removed_monolithic_input.csv')
        manifest.text = header.join(',') + '\n'
        return manifest.absolutePath
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
