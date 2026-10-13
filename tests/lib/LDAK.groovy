// LDAK HE and PCGC tables, snapshotted without the two columns that differ between CPU vendors.
//
// LDAK 6.1 links Intel MKL, which takes a different code path on Intel and non-Intel CPUs, so the same run
// does not give byte-identical tables on both (lyh970817/gwas#131). In the HE and PCGC component tables only
// `Mega_Intensity` and its `SE` differ. They restate `Heritability` and its `SE` per megabase of `Size`
// (heritability / Size * 1e6) with six decimals, which exposes digits that the six-decimal `Heritability`
// column rounds away. Dropping them loses nothing the table does not already state, and every other line
// is snapshotted exactly.
//
// A test passes HERITABILITY_TABLES as `ignore:` to the `getAllFilesFromDir` call that hashes file content,
// so those files keep only their names there, and snapshots `heritabilityTables` beside it.
class LDAK {

    static final List<String> HERITABILITY_TABLES = [
        '**/*.ldak_he.he',
        '**/*.ldak_he.he.liab',
        '**/*.ldak_pcgc.pcgc',
        '**/*.ldak_pcgc.pcgc.marginal',
    ]

    // Every HE or PCGC table under `directory`, keyed by its path relative to `directory`, as its lines
    // without the `Mega_Intensity` column and the `SE` column that follows it.
    static Map<String, List<String>> heritabilityTables(Object directory) {
        def root = new File(directory.toString()).toPath()
        def matchers = HERITABILITY_TABLES.collect { glob -> root.fileSystem.getPathMatcher("glob:${glob}") }
        def tables = [:]
        root.toFile().eachFileRecurse(groovy.io.FileType.FILES) { table ->
            def relative = root.relativize(table.toPath())
            if (matchers.any { matcher -> matcher.matches(relative) }) {
                tables[relative.toString()] = withoutMegaIntensity(table)
            }
        }
        return tables.sort()
    }

    private static List<String> withoutMegaIntensity(File table) {
        def dropped = -1
        return table.readLines().collect { line ->
            def fields = line.tokenize()
            if (fields && fields.first() == 'Component') {
                dropped = fields.indexOf('Mega_Intensity')
            }
            return dropped < 0 ? line : (fields.take(dropped) + fields.drop(dropped + 2)).join(' ')
        }
    }
}
