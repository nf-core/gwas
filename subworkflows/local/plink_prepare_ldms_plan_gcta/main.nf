// Build one LD-by-MAF component plan: the LD scores of a focal PLINK 1 bundle and the ordered, disjoint SNP
// groups derived from them. The plan is the shared input of every LD- and MAF-stratified matrix family, so it
// is built here once rather than inside any one tool's matrix builder.
//
// GCTA's LD-score pass runs once per chromosome, each a task of its own, and the per-chromosome tables are
// joined in the BIM's chromosome order into the one table the stratifier reads. The scatter does not change
// what the stratifier reads: `ldscore_SNP` and `freq` are the same per chromosome as genome-wide. GCTA's
// segment-based `ldscore_region` becomes a within-chromosome segment mean, where a genome-wide pass lets the
// segments at a chromosome end draw on the next chromosome.
// Every constituent process reports directly to the run-wide versions topic, so this subworkflow emits no versions.

// MODULE: Local to the pipeline
include { CUSTOM_PLINKBIMCHROMOSOMES  } from '../../../modules/local/custom/plinkbimchromosomes/main'
include { GCTA_CALCULATELDSCORES      } from '../../../modules/local/gcta/calculateldscores/main'
include { CUSTOM_GCTACONCATLDSCORES   } from '../../../modules/local/custom/gctaconcatldscores/main'
include { CUSTOM_GCTASTRATIFYLDSCORES } from '../../../modules/local/custom/gctastratifyldscores/main'

workflow PLINK_PREPARE_LDMS_PLAN_GCTA {
    take:
    ch_genotypes // channel: [ val(meta), path(bed), path(bim), path(fam) ], focal genotypes
    ch_ld_score_region_kb // channel: [ val(meta2), val(ld_score_region_kb) ], LD-score region size
    ch_ld_bins // channel: [ val(meta3), val(ld_bins) ], LD-score bin count
    ch_maf_edges // channel: [ val(meta4), val(maf_edges) ], MAF edges

    main:
    ch_plan_inputs = ch_genotypes
        .map { meta, bed, bim, fam -> tuple(meta.id, meta, bed, bim, fam) }
        .join(ch_ld_score_region_kb.map { meta2, ld_score_region_kb -> tuple(meta2.id, ld_score_region_kb) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_ld_bins.map { meta3, ld_bins -> tuple(meta3.id, ld_bins) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .join(ch_maf_edges.map { meta4, maf_edges -> tuple(meta4.id, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)

    // The chromosome axis is the BIM's own first column, listed in a task because a BIM can hold millions of
    // rows. The listing it returns has one line per chromosome, so reading it here is cheap.
    CUSTOM_PLINKBIMCHROMOSOMES(
        ch_plan_inputs.map { _plan_id, meta, _bed, bim, _fam, _ld_score_region_kb, _ld_bins, _maf_edges -> tuple(meta, bim) }
    )
    ch_chromosomes = CUSTOM_PLINKBIMCHROMOSOMES.out.chromosomes
        .map { meta, chromosome_list -> tuple(meta.id, chromosome_list.readLines().findAll { line -> line.trim() }) }

    // Every chromosome task stages the same bundle and differs only in the chromosome it passes to `--chr`.
    // The chromosome also rides in a private metadata key so the gather can restore BIM order; the key is
    // dropped again when the plan's own metadata is restored below.
    ch_ld_score_inputs = ch_plan_inputs
        .map { plan_id, meta, bed, bim, fam, ld_score_region_kb, _ld_bins, _maf_edges -> tuple(plan_id, meta, bed, bim, fam, ld_score_region_kb) }
        .join(ch_chromosomes, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .flatMap { _plan_id, meta, bed, bim, fam, ld_score_region_kb, chromosomes ->
            chromosomes.collect { chromosome -> tuple(meta + [ld_score_chromosome: chromosome], bed, bim, fam, ld_score_region_kb, chromosome) }
        }
    GCTA_CALCULATELDSCORES(ch_ld_score_inputs)

    // The gather releases a plan as soon as all of its chromosomes have arrived, and orders the tables by
    // the chromosome's position in the BIM rather than by completion order or file name.
    ch_concat_inputs = GCTA_CALCULATELDSCORES.out.ld_scores
        .map { chromosome_meta, ld_scores -> tuple(chromosome_meta.id, chromosome_meta.ld_score_chromosome, ld_scores) }
        .combine(ch_chromosomes, by: 0)
        .map { plan_id, chromosome, ld_scores, chromosomes -> tuple(groupKey(plan_id, chromosomes.size()), chromosomes.indexOf(chromosome), ld_scores) }
        .groupTuple()
        .map { plan_key, positions, ld_scores ->
            tuple(plan_key.getGroupTarget(), [positions, ld_scores].transpose().sort(false) { entry -> entry[0] }.collect { entry -> entry[1] })
        }
        .join(ch_plan_inputs.map { plan_id, meta, _bed, _bim, _fam, _ld_score_region_kb, _ld_bins, _maf_edges -> tuple(plan_id, meta) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _plan_id, ld_scores, meta -> tuple(meta, ld_scores) }
    CUSTOM_GCTACONCATLDSCORES(ch_concat_inputs)

    ch_stratify_inputs = CUSTOM_GCTACONCATLDSCORES.out.ld_scores
        .map { meta, ld_scores -> tuple(meta.id, meta, ld_scores) }
        .join(ch_plan_inputs.map { plan_id, _meta, _bed, _bim, _fam, _ld_score_region_kb, ld_bins, maf_edges -> tuple(plan_id, ld_bins, maf_edges) }, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _plan_id, meta, ld_scores, ld_bins, maf_edges -> tuple(meta, ld_scores, ld_bins, maf_edges) }
    CUSTOM_GCTASTRATIFYLDSCORES(ch_stratify_inputs)

    emit:
    plan      = CUSTOM_GCTASTRATIFYLDSCORES.out.strata_bundle // channel: [ val(meta), path(strata_manifest), [ path(snp_group_file), ... ] ], LD-major MAF-minor row order
    groups    = CUSTOM_GCTASTRATIFYLDSCORES.out.groups // channel: [ val(meta), path(group_manifest), path(snp_group_files) ], same row order as plan
    ld_scores = CUSTOM_GCTACONCATLDSCORES.out.ld_scores // channel: [ val(meta), path(ld_scores) ], per-chromosome tables joined in BIM chromosome order
}
