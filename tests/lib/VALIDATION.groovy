// Collects everything a failed run said, so an input-validation test can assert on the message a
// researcher actually sees.
//
// The two validation passes surface their errors in different places. The Groovy preflight raises
// from the workflow body and its message reaches the console, while nf-schema rejects the samplesheet
// during parameter validation and prints only a summary line to the console, leaving the per-row
// detail in the run log. Both are joined here so a test does not have to know which pass caught it.
class VALIDATION {

    // Nextflow 26.10 prints errors and log messages on stderr rather than stdout, so read both streams.
    static String console(Object run) {
        return ((run.stdout ?: []) + (run.stderr ?: [])).join('\n')
    }

    static String report(Object outputDir, Object run) {
        def log = new File(new File(outputDir.toString()).parentFile, 'meta/nextflow.log')
        return [console(run), log.exists() ? log.text : ''].join('\n')
    }
}
