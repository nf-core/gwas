// Collects everything a failed run said, so an input-validation test can assert on the message a
// researcher actually sees.
//
// The two validation passes surface their errors in different places. The Groovy preflight raises
// from the workflow body and its message reaches the console, while nf-schema rejects the samplesheet
// during parameter validation and prints only a summary line to the console, leaving the per-row
// detail in the run log. Both are joined here so a test does not have to know which pass caught it.
class VALIDATION {

    static String report(Object outputDir, List<String> consoleOutput) {
        def log = new File(new File(outputDir.toString()).parentFile, 'meta/nextflow.log')
        return [consoleOutput ? consoleOutput.join('\n') : '', log.exists() ? log.text : ''].join('\n')
    }
}
