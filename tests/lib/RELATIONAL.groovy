// Linked-manifest builders for pipeline-level tests.
//
// Fixture paths name the portable nf-core/test-datasets GWAS bundle and are resolved through FIXTURES.base,
// so executable tests consistently use the public nf-core/test-datasets fixture bundle.
class RELATIONAL {

    // Gzip member metadata is not scientific content. Snapshot decompressed text with
    // normalised line endings while retaining the relative filename for attribution.
    static List<String> gzipTextHashes(Object directory) {
        def root = new File(directory.toString())
        def hashes = []
        root.eachFileRecurse { entry ->
            if (entry.isFile() && entry.name.endsWith('.gz')) {
                def relative = root.toPath().relativize(entry.toPath()).toString()
                def text = new java.util.zip.GZIPInputStream(new FileInputStream(entry))
                    .withReader('UTF-8') { reader -> reader.readLines().join('\n') }
                def digest = java.security.MessageDigest.getInstance('MD5').digest(text.getBytes('UTF-8'))
                def md5 = digest.collect { value -> String.format('%02x', value & 0xff) }.join()
                hashes << "${relative}:md5,${md5}"
            }
        }
        return hashes.sort()
    }

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

    static String relationships(Object projectDir, Object outputDir, String name, Closure mutate = null) {
        def header = [
            'relationship_id',
            'left_analysis_id',
            'right_analysis_id',
            'left_summary_statistics_id',
            'right_summary_statistics_id',
            'relationship_methods',
            'pair_quant_covariates',
            'pair_cat_covariates',
        ]
        def rows = [relationship('qt_bt')]
        if (mutate) {
            mutate(rows)
        }
        return materialise(
            projectDir,
            outputDir,
            name,
            'relationships',
            header,
            rows,
            ['pair_quant_covariates', 'pair_cat_covariates'],
        )
    }

    static String summaryStatistics(Object projectDir, Object outputDir, String name, Closure mutate = null) {
        def header = [
            'summary_statistics_id',
            'trait_id',
            'trait_type',
            'source',
            'source_mode',
            'source_format',
            'producer_analysis_id',
            'producer_association_method',
            'genome_build',
            'ancestry',
            'source_method',
            'source_release',
            'heritability_methods',
            'population_prevalence',
            'sample_prevalence',
            'access_constraints',
        ]
        def rows = [summary(outputDir, 'external_qt')]
        if (mutate) {
            mutate(rows)
        }
        return materialise(projectDir, outputDir, name, 'summary_statistics', header, rows, ['source'])
    }

    static Map summary(Object outputDir, String summaryStatisticsId) {
        def summaries = [
            external_qt: [
                summary_statistics_id: 'external_qt',
                trait_id: 'QT_external',
                trait_type: 'quantitative',
                source: canonicalSummary(outputDir, 'external_qt'),
                source_mode: 'canonical',
                source_format: 'nfcore_gwas_canonical_v1',
                producer_analysis_id: '',
                producer_association_method: '',
                genome_build: 'GRCh37',
                ancestry: 'EUR',
                source_method: 'published_gwas',
                source_release: 'v1',
                heritability_methods: 'ldsc_h2',
                population_prevalence: '',
                sample_prevalence: '',
                access_constraints: 'controlled_access',
            ],
            external_bt: [
                summary_statistics_id: 'external_bt',
                trait_id: 'BT_external',
                trait_type: 'binary',
                source: canonicalSummary(outputDir, 'external_bt'),
                source_mode: 'canonical',
                source_format: 'nfcore_gwas_canonical_v1',
                producer_analysis_id: '',
                producer_association_method: '',
                genome_build: 'GRCh37',
                ancestry: 'EUR',
                source_method: 'published_gwas',
                source_release: 'v2',
                heritability_methods: '',
                population_prevalence: 0.1,
                sample_prevalence: 0.2,
                access_constraints: '',
            ],
        ]
        if (!summaries.containsKey(summaryStatisticsId)) {
            throw new IllegalArgumentException("No shipped summary fixture '${summaryStatisticsId}'")
        }
        return new LinkedHashMap(summaries[summaryStatisticsId])
    }

    static Map internalSummary(String analysisId, String associationMethod, String heritabilityMethods = '') {
        return [
            summary_statistics_id: "${analysisId}--${associationMethod}",
            trait_id: '',
            trait_type: '',
            source: '',
            source_mode: '',
            source_format: '',
            producer_analysis_id: analysisId,
            producer_association_method: associationMethod,
            genome_build: '',
            ancestry: '',
            source_method: '',
            source_release: '',
            heritability_methods: heritabilityMethods,
            population_prevalence: '',
            sample_prevalence: '',
            access_constraints: '',
        ]
    }

    static String summaryRelationships(Object projectDir, Object outputDir, String name, Closure mutate = null) {
        def header = [
            'relationship_id',
            'left_analysis_id',
            'right_analysis_id',
            'left_summary_statistics_id',
            'right_summary_statistics_id',
            'relationship_methods',
            'pair_quant_covariates',
            'pair_cat_covariates',
        ]
        def rows = [[
            relationship_id: 'external_pair',
            left_analysis_id: '',
            right_analysis_id: '',
            left_summary_statistics_id: 'external_qt',
            right_summary_statistics_id: 'external_bt',
            relationship_methods: 'ldsc_rg,ldak_sumcors',
            pair_quant_covariates: '',
            pair_cat_covariates: '',
        ]]
        if (mutate) {
            mutate(rows)
        }
        return materialise(
            projectDir,
            outputDir,
            name,
            'summary_relationships',
            header,
            rows,
            ['pair_quant_covariates', 'pair_cat_covariates'],
        )
    }

    static String referenceCatalog(Object outputDir, String name) {
        def hapmap3 = resource(outputDir, "${name}/hapmap3.snp", "rs1\n")
        def resourceRoot = new File(new File(outputDir.toString()).parentFile, "resources/${name}")
        def referenceLd = new File(resourceRoot, 'reference')
        def regressionWeights = new File(resourceRoot, 'weights')
        referenceLd.mkdirs()
        regressionWeights.mkdirs()
        new File(referenceLd, '1.l2.ldscore').text = "CHR SNP BP L2\n1 rs1 1 1\n"
        new File(regressionWeights, '1.l2.ldscore').text = "CHR SNP BP L2\n1 rs1 1 1\n"
        def tagging = resource(outputDir, "${name}/reference.tagging", "Predictor Tagging\nrs1 1\n")
        def document = [
            ldsc: [
                ldsc_eur: [
                    genome_build: 'GRCh37',
                    ancestry: 'EUR',
                    variant_id_system: 'rsid',
                    hapmap3_snplist: hapmap3,
                    reference_ld_scores: referenceLd.absolutePath,
                    regression_weights: regressionWeights.absolutePath,
                ],
            ],
            ldak: [
                ldak_thin_eur: [
                    genome_build: 'GRCh37',
                    ancestry: 'EUR',
                    variant_id_system: 'rsid',
                    model: 'LDAK-Thin',
                    tagging_file: tagging,
                ],
            ],
        ]
        def directory = new File(new File(outputDir.toString()).parentFile, "manifests/${name}")
        directory.mkdirs()
        def catalog = new File(directory, 'reference_catalog.json')
        catalog.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(document)) + '\n'
        return catalog.absolutePath
    }

    // Build one deterministic, self-contained LDSC fixture family large enough for the native regression
    // smoke tests. These files are generated inside nf-test output state and are not pipeline dependencies.
    static Map ldscResources(Object outputDir, String name) {
        def root = new File(new File(outputDir.toString()).parentFile, "resources/${name}")
        def reference = new File(root, 'reference')
        def weights = new File(root, 'weights')
        reference.mkdirs()
        weights.mkdirs()
        def gzipText = { File target, String body ->
            target.withOutputStream { output ->
                def gzip = new java.util.zip.GZIPOutputStream(output)
                gzip.write(body.getBytes('UTF-8'))
                gzip.close()
            }
        }
        def hapmap3 = new File(root, 'w_hm3.snplist')
        def alleleRows = []
        (1..2).each { chromosome ->
            def referenceRows = []
            def weightRows = []
            (1..500).each { index ->
                def globalIndex = (chromosome - 1) * 500 + index
                def ldScore = 1.0 + (globalIndex % 47) / 10.0
                def weight = 1.0 + (globalIndex % 13) / 20.0
                referenceRows << "${chromosome}\trs${globalIndex}\t${globalIndex}\t${ldScore}"
                weightRows << "${chromosome}\trs${globalIndex}\t${globalIndex}\t${weight}"
                alleleRows << "rs${globalIndex}\tA\tG"
            }
            gzipText(new File(reference, "${chromosome}.l2.ldscore.gz"), "CHR\tSNP\tBP\tL2\n${referenceRows.join('\n')}\n")
            new File(reference, "${chromosome}.l2.M_5_50").text = '500\n'
            gzipText(new File(weights, "${chromosome}.l2.ldscore.gz"), "CHR\tSNP\tBP\tL2\n${weightRows.join('\n')}\n")
        }
        hapmap3.text = "SNP\tA1\tA2\n${alleleRows.join('\n')}\n"

        def summary = { String id, double coefficient, int sampleSize ->
            def source = new File(root, "${id}.canonical.tsv")
            def rows = (1..1000).collect { index ->
                def chromosome = index <= 500 ? 1 : 2
                def localIndex = index <= 500 ? index : index - 500
                def ldScore = 1.0 + (index % 47) / 10.0
                def sign = index % 2 == 0 ? 1.0 : -1.0
                def z = sign * Math.sqrt(1.0 + coefficient * ldScore) + ((index % 7) - 3) / 20.0
                def se = 0.01
                def beta = z * se
                "rs${index}\t${chromosome}\t${localIndex}\tA\tG\t1900000\t0.25\t${String.format(java.util.Locale.ROOT, '%.8f', beta)}\t${se}\t0.001\t${sampleSize}"
            }
            source.text = "SNPID\tCHR\tPOS\tEA\tNEA\tSTATUS\tEAF\tBETA\tSE\tP\tN\n${rows.join('\n')}\n"
            source.absolutePath
        }

        def document = [
            ldsc: [
                ldsc_eur: [
                    genome_build: 'GRCh37',
                    ancestry: 'EUR',
                    variant_id_system: 'rsid',
                    hapmap3_snplist: hapmap3.absolutePath,
                    reference_ld_scores: reference.absolutePath,
                    regression_weights: weights.absolutePath,
                ],
            ],
        ]
        def catalog = new File(root, 'reference_catalog.json')
        catalog.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(document)) + '\n'
        return [
            catalog: catalog.absolutePath,
            quantitative: summary('external_qt', 0.5, 10000),
            binary: summary('external_bt', 0.7, 12000),
            hapmap3: hapmap3.absolutePath,
            reference: reference.absolutePath,
            weights: weights.absolutePath,
        ]
    }

    static String summaryMethodOptions(Object outputDir, String name) {
        return methodOptions(outputDir, name, [
            unary_requests: [
                'ldsc_h2--external_qt': [reference_bundle_id: 'ldsc_eur'],
            ],
            pair_requests: [
                'ldsc_rg--external_pair': [reference_bundle_id: 'ldsc_eur'],
                'ldak_sumcors--external_pair': [reference_bundle_id: 'ldak_thin_eur'],
            ],
        ])
    }

    static String canonicalSummary(Object outputDir, String name) {
        return resource(
            outputDir,
            "${name}.canonical.tsv",
            'SNPID\tCHR\tPOS\tEA\tNEA\tSTATUS\tEAF\tBETA\tSE\tP\tN\n' +
                'rs1\t1\t1\tA\tG\t1900000\t0.25\t0.1\t0.01\t1e-4\t1000\n',
        )
    }

    static Map relationship(String relationshipId) {
        def fixture = { path -> "${FIXTURES.UPSTREAM}results/fixtures/${path}" }
        def relationships = [
            qt_bt: [
                relationship_id: 'qt_bt',
                left_analysis_id: 'example_pgen_qt',
                right_analysis_id: 'example_pgen_bt',
                left_summary_statistics_id: '',
                right_summary_statistics_id: '',
                relationship_methods: 'gcta_bivariate_reml',
                pair_quant_covariates: fixture('pheno_cov/example.qcovar'),
                pair_cat_covariates: fixture('pheno_cov/example.catcovar'),
            ],
        ]
        if (!relationships.containsKey(relationshipId)) {
            throw new IllegalArgumentException("No shipped relationship fixture '${relationshipId}'")
        }
        return new LinkedHashMap(relationships[relationshipId])
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
