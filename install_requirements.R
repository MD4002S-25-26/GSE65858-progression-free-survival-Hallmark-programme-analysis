# Installs any missing packages needed by the pipeline. Run once, explicitly:
#   Rscript install_requirements.R
# Versions used for the manuscript are listed in README.md.
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("dplyr", "tidyr", "tibble", "ggplot2", "ggrepel", "survival", "msigdbr",
          "gtsummary", "flextable", "officer")
bioc <- c("limma", "fgsea")
miss <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss(cran))) install.packages(miss(cran))
if (length(miss(bioc))) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss(bioc), ask = FALSE, update = FALSE)
}
cat("Missing packages:", if (length(miss(c(cran, bioc)))) paste(miss(c(cran, bioc)), collapse = ", ") else "none", "\n")
