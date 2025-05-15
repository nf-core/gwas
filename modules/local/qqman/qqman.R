
# Load required packages
suppressPackageStartupMessages({
  library(argparse)
  library(qqman)
})

# ---- ARGUMENT PARSING ----
parser <- ArgumentParser(description = "Run qqman on a PLINK .assoc file")

parser$add_argument("-i", "--input", required = TRUE,
                    help = "Input PLINK .assoc file")
parser$add_argument("-m", "--manhattan", default = "manhattan_plot.png",
                    help = "Output filename for Manhattan plot [default: manhattan_plot.png]")
parser$add_argument("-q", "--qq", default = "qq_plot.png",
                    help = "Output filename for Q-Q plot [default: qq_plot.png]")

args <- parser$parse_args()

# ---- LOAD DATA ----
assoc_data <- read.table(args$input, header = TRUE)

# Validate required columns
required_cols <- c("CHR", "BP", "SNP", "P")
missing_cols <- setdiff(required_cols, colnames(assoc_data))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in .assoc file:", paste(missing_cols, collapse = ", ")))
}

# ---- PLOTTING ----

# Manhattan plot
png(args$manhattan, width = 1000, height = 600)
manhattan(assoc_data, main = "Manhattan Plot")
dev.off()

# Q-Q plot
png(args$qq, width = 600, height = 600)
qq(assoc_data$P, main = "Q-Q Plot")
dev.off()

cat("Plots saved:\n")
cat("  Manhattan plot:", args$manhattan, "\n")
cat("  Q-Q plot:", args$qq, "\n")
