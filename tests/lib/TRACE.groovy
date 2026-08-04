// Reads the execution trace of a finished run.
//
// "This expensive step ran exactly once, however many analysis units asked for it" is a claim about
// what was executed, not about what was published, so it is answered from the trace rather than by
// counting output files — a cached or reused task publishes the same files as a re-executed one.
//
// nf-test exposes each trace row as a task whose `name` is the fully qualified process name followed
// by the process tag in parentheses, for example
// `NFCORE_GWAS:GWAS:PREPARE_COHORT_GENOTYPES:PLINK2_VCF (example_vcf)`.
class TRACE {

    // How many tasks of `process` the run executed.
    static int count(Object workflow, String process) {
        return tasks(workflow, process).size()
    }

    // The process tags of every task of `process` the run executed, sorted. For a per-cohort step the
    // tag is the cohort identifier, so this says which cohorts were prepared as well as how many.
    static List<String> tags(Object workflow, String process) {
        return tasks(workflow, process).collect { task -> tag(task.name) }.sort()
    }

    private static List tasks(Object workflow, String process) {
        return workflow.trace.succeeded().findAll { task -> name(task.name) == process }
    }

    private static String name(String task_name) {
        return task_name.split(/ \(/).first().tokenize(':').last().trim()
    }

    private static String tag(String task_name) {
        def opened = task_name.indexOf('(')
        return opened < 0 ? '' : task_name.substring(opened + 1, task_name.lastIndexOf(')'))
    }
}
