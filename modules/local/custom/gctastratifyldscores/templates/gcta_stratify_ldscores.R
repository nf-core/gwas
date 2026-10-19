#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT

arguments <- c('${ld_scores}', '${prefix}', '${ld_bins}', '${maf_edges.join(",")}')
process_name <- ${task_process_literal}

input_file <- arguments[[1]]
output_prefix <- arguments[[2]]
ld_bin_value <- suppressWarnings(as.numeric(arguments[[3]]))
maf_boundaries <- suppressWarnings(as.numeric(strsplit(arguments[[4]], ",", fixed = TRUE)[[1]]))

fail <- function(message) stop(message, call. = FALSE)
format_bound <- function(value) format(signif(value, 8), scientific = FALSE, trim = TRUE)

if (length(ld_bin_value) != 1L || !is.finite(ld_bin_value) || ld_bin_value < 1L || ld_bin_value != floor(ld_bin_value)) {
  fail("ld_bins must be one positive integer")
}
ld_bin_count <- as.integer(ld_bin_value)

if (length(maf_boundaries) < 2L || any(!is.finite(maf_boundaries))) {
  fail("maf_edges must contain at least two finite numeric boundaries")
}
if (any(diff(maf_boundaries) <= 0)) fail("maf_edges must be strictly increasing")
if (any(maf_boundaries < 0 | maf_boundaries > 0.5)) fail("maf_edges must be between 0 and 0.5")

scores <- read.table(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

ld_scores <- scores[["ldscore_SNP"]]
maf <- pmin(scores[["freq"]], 1 - scores[["freq"]])

ld_boundaries <- as.numeric(quantile(
  ld_scores,
  probs = seq(0, 1, length.out = ld_bin_count + 1L),
  names = FALSE,
  type = 7
))

manifest_rows <- list()
row_index <- 1L
assignment_count <- integer(nrow(scores))

for (ld_index in seq_len(ld_bin_count)) {
  ld_lower <- ld_boundaries[ld_index]
  ld_upper <- ld_boundaries[ld_index + 1L]
  ld_member <- ld_scores >= ld_lower & if (ld_index == ld_bin_count) ld_scores <= ld_upper else ld_scores < ld_upper

  for (maf_index in seq_len(length(maf_boundaries) - 1L)) {
    maf_lower <- maf_boundaries[maf_index]
    maf_upper <- maf_boundaries[maf_index + 1L]
    maf_member <- maf >= maf_lower & if (maf_index == length(maf_boundaries) - 1L) maf <= maf_upper else maf < maf_upper
    member <- ld_member & maf_member
    assignment_count <- assignment_count + as.integer(member)
    predictor_count <- sum(member)
    if (predictor_count == 0L) next

    stratum_key <- sprintf("ld%02d_maf%02d", ld_index, maf_index)
    group_filename <- sprintf("%s_snp_group_%s.txt", output_prefix, stratum_key)
    writeLines(as.character(scores[["SNP"]][member]), group_filename)

    manifest_rows[[row_index]] <- data.frame(
      model_key = "ldms",
      stratum_key = stratum_key,
      ld_lower = format_bound(ld_lower),
      ld_upper = format_bound(ld_upper),
      maf_lower = format_bound(maf_lower),
      maf_upper = format_bound(maf_upper),
      predictor_count = predictor_count,
      group_filename = group_filename,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}

if (any(assignment_count != 1L)) fail("Each predictor must be assigned to exactly one stratum")
if (length(manifest_rows) == 0L) fail("No non-empty strata produced; check ld_bins and maf_edges")

manifest <- do.call(rbind, manifest_rows)
write.table(manifest, paste0(output_prefix, ".strata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
groups <- data.frame(group_id = manifest[["stratum_key"]], group_filename = manifest[["group_filename"]])
write.table(groups, paste0(output_prefix, ".groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Written here rather than captured by an `eval` output, which Nextflow allows only on a Bash script.
writeLines(
  c(
    sprintf('"%s":', process_name),
    "    gcta_stratify_ldscores: 1.0.0",
    sprintf("    r-base: %s", as.character(getRversion()))
  ),
  "versions.yml"
)
