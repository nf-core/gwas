#!/usr/bin/env Rscript

ldak_reml_file <- "$ldak_reml_file"
quarter_files <- c($quarter_reml_files_r)

parse_reml_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Required REML file does not exist: ", file_path, call. = FALSE)
  }

  lines <- readLines(file_path)
  her_all_line <- grep("^Her_All[[:space:]]", lines, value = TRUE)

  if (length(her_all_line) != 1L) {
    stop(
      "Expected exactly one Her_All row in required REML file ",
      file_path,
      "; found ",
      length(her_all_line),
      call. = FALSE
    )
  }

  parts <- strsplit(her_all_line, "\\\\s+")[[1]]
  if (length(parts) < 3L) {
    stop("Malformed Her_All row in required REML file: ", file_path, call. = FALSE)
  }

  estimates <- type.convert(parts[2:3], as.is = TRUE)
  if (!is.numeric(estimates) || anyNA(estimates) || any(!is.finite(estimates))) {
    stop(
      "Her_All heritability and standard error must be finite numbers in required REML file: ",
      file_path,
      call. = FALSE
    )
  }

  data.frame(
    file = basename(file_path),
    heritability = estimates[1],
    se = estimates[2],
    stringsAsFactors = FALSE
  )
}

if (length(quarter_files) < 2L) {
  stop("At least two quarter REML files are required for the n - 1 inflation formula", call. = FALSE)
}

quarter_results <- lapply(quarter_files, parse_reml_file)
quarter_data <- do.call(rbind, quarter_results)
quarter_data\$type <- "quarter"

ldak_data <- parse_reml_file(ldak_reml_file)
ldak_data\$type <- "ldak"

all_results <- rbind(quarter_data, ldak_data)

n_quarters <- nrow(quarter_data)
quarter_mean_h2 <- mean(quarter_data\$heritability)
quarter_mean_se <- mean(quarter_data\$se)
ldak_h2 <- ldak_data\$heritability[1]
ldak_se <- ldak_data\$se[1]

quarter_sum_h2 <- sum(quarter_data\$heritability)
inflation_T1 <- (quarter_sum_h2 - ldak_h2) / (n_quarters - 1)
inflation_factor <- quarter_mean_h2 / ldak_h2

mean_t1 <- inflation_T1
sd_t1 <- sqrt(sum(quarter_data\$se^2) + ldak_se^2) / (n_quarters - 1)

statistical_test_results <- list(
  n_quarters = n_quarters,
  pvalue = pnorm(0, mean = mean_t1, sd = sd_t1),
  mean_T1samp = mean_t1,
  sd_T1samp = sd_t1
)

output_lines <- c(
  "LDAK Inflation Analysis Results",
  "================================",
  "",
  paste("Number of quarter files processed:", length(quarter_files)),
  paste("LDAK file processed:", basename(ldak_reml_file)),
  "",
  "Quarter Results:",
  paste("  Mean Heritability:", round(quarter_mean_h2, 6)),
  paste("  Mean SE:", round(quarter_mean_se, 6)),
  "",
  "LDAK Results:",
  paste("  Heritability:", round(ldak_h2, 6)),
  paste("  SE:", round(ldak_se, 6)),
  "",
  "Inflation Analysis:",
  paste("  T1 Statistic (Documentation Formula):", round(inflation_T1, 6)),
  paste("  Inflation Ratio (Legacy, Quarter/LDAK):", round(inflation_factor, 6)),
  "  Interpretation:",
  "    - T1 ~= 0: No inflation (good quality control)",
  "    - T1 > 0.05: Possible population structure or relatedness not fully captured",
  "    - T1 < -0.05: Possible over-correction",
  "",
  "Statistical Test Results:",
  paste("  Number of quarters used:", statistical_test_results\$n_quarters),
  paste("  P-value:", round(statistical_test_results\$pvalue, 6)),
  paste("  Mean T1 (analytical):", round(statistical_test_results\$mean_T1samp, 6)),
  paste("  SD T1 (analytical):", round(statistical_test_results\$sd_T1samp, 6)),
  "",
  "Individual Results:"
)

for (i in seq_len(nrow(all_results))) {
  row <- all_results[i, ]
  output_lines <- c(
    output_lines,
    paste(
      " ", row\$file, "(", row\$type, "):",
      "H2 =", round(row\$heritability, 6),
      "SE =", round(row\$se, 6)
    )
  )
}

writeLines(output_lines, "${prefix}.txt")

r_version <- system("R --version | sed -n '1s/.*version //;s/ .*//p'", intern = TRUE)
writeLines(c(
  "\"${task.process}\":",
  paste0("    r-base: ", r_version[1])
), "versions.yml")
