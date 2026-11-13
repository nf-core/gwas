// Re-runs a finished pipeline test with `-resume`, so a test can assert what a resumed run reused.
//
// "This matrix was not rebuilt when a trait was added to the samplesheet" is a claim about the task
// cache, and the task cache is only observable across two runs that share a work directory. nf-test
// drives exactly one `nextflow run` per test and passes its own `-w`, which Nextflow refuses to see
// twice (`Can only specify option -w once`), so a second test can never be pointed at the first
// one's work directory and there is no two-stage resume mode to reach for.
//
// What nf-test does leave behind is everything a resume needs. It launches Nextflow from the test
// directory, one level above `params.outdir`, so that run's `.nextflow/history` and `.nextflow/cache`
// sit there, and its work directory does too. Re-launching from there resumes that same run.
//
// The command is reconstructed from the launch line Nextflow wrote into the run log rather than
// rebuilt by hand, so the profile, the config files, the params file and the work directory are
// whatever nf-test actually used and cannot drift from how the suite was invoked. Two consequences
// worth knowing: the line is tokenised on spaces, so a repository path containing a space would break
// it, and the resumed run is a second full pipeline execution, so a test using this costs roughly
// double.
class RESUME {

    // Re-run the finished test with `-resume` and the given parameter overrides, and return the rows
    // of the resumed run's execution trace as maps keyed by the trace header. Command-line `--input`
    // and `--outdir` override the `-params-file` nf-test passes, which is what lets the second run
    // use a samplesheet the first one did not have.
    static List<Map> rerun(Object outputDir, Map overrides) {
        def test_dir = new File(outputDir.toString()).parentFile
        def log = new File(test_dir, 'meta/nextflow.log')
        def launch = log.readLines().find { line -> line.contains('$> nextflow ') }
        if (!launch) {
            throw new IllegalStateException("no nextflow launch command found in ${log}")
        }
        def argv = launch.substring(launch.indexOf('$> ') + 3).tokenize(' ')

        // The first run's log and trace destinations are the only two options that must not be
        // reused; everything else — the config files, the params file, the profile, the work
        // directory — is exactly what makes this a resume of that run rather than a fresh one.
        def trace = new File(test_dir, "resume_trace_${java.util.UUID.randomUUID()}.txt")
        def command = []
        for (int index = 0; index < argv.size(); index++) {
            if (argv[index] in ['-log', '-with-trace']) {
                index++
                continue
            }
            command << argv[index]
        }
        command += ['-resume', '-with-trace', trace.absolutePath]
        overrides.each { name, value -> command += ["--${name}", value.toString()] }

        def process = new ProcessBuilder(command.collect { argument -> argument.toString() })
            .directory(test_dir)
            .redirectErrorStream(true)
            .start()
        def output = process.inputStream.text
        if (process.waitFor() != 0) {
            throw new IllegalStateException("resumed run failed:\n${output}")
        }

        def rows = trace.readLines().findAll { line -> line.trim() }
        def header = rows.first().split('\t').toList()
        return rows.tail().collect { row -> [header, row.split('\t').toList()].transpose().collectEntries() }
    }

    // Run the same cache probe while also recording the canonical identity order accepted by a terminal
    // one-task-per-output consumer. The sort makes the assembled-result contract independent of task
    // scheduling while retaining duplicates, so both cardinality and identity order remain observable.
    static Map rerunWithOutputs(Object outputDir, Map overrides, String outputProcess) {
        def trace = rerun(outputDir, overrides)
        return [trace: trace, outputIdentities: taskTags(trace, outputProcess)]
    }

    // The trace rows of one process, by simple process name, matching how TRACE reads nf-test's own.
    static List<Map> tasks(List<Map> trace, String process) {
        return trace.findAll { row -> row.name.split(/ \(/).first().tokenize(':').last().trim() == process }
    }

    // Canonically ordered task tags are output identities when the observed process consumes exactly one task
    // per assembled result, as GWASLAB_HARMONIZE does for association outputs.
    static List<String> taskTags(List<Map> trace, String process) {
        return tasks(trace, process)
            .collect { row -> (row.name =~ /\(([^()]*)\)$/)[0][1] }
            .sort()
    }

    // The task statuses of one process, sorted, so a test can say "every one of these was reused" and
    // "exactly one of those was not" without caring which task is which.
    static List<String> statuses(List<Map> trace, String process) {
        return tasks(trace, process).collect { row -> row.status }.sort()
    }

    // Exact task hashes let cache-boundary tests prove identity directly rather than treating a
    // CACHED/COMPLETED status alone as a proxy for which task Nextflow matched.
    static List<String> hashes(List<Map> trace, String process) {
        return tasks(trace, process).collect { row -> row.hash }.sort()
    }
}
