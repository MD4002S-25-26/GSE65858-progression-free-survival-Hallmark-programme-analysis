# ============================================================
# 00a_prepare_inputs_from_GEO.R
#
# Purpose:
#   Build the two input objects used by the analysis scripts
#   (pheno3_matched.rds and expr2_probe_matrix.rds), and fetch the
#   GPL10558 platform annotation, directly from the public GEO files.
#
# Run from the repository root:  source("scripts/00a_prepare_inputs_from_GEO.R")
#
# Inputs (downloaded if absent, checksum-verified):
#   - GSE65858_series_matrix.txt.gz   (NCBI GEO, series GSE65858)
#   - GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz (NCBI GEO, platform GPL10558)
# Outputs (repository root):
#   - GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt
#   - pheno3_matched.rds
#   - expr2_probe_matrix.rds
# ============================================================

suppressPackageStartupMessages(library(dplyr))
options(timeout = max(1200, getOption("timeout")))

SERIES_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65858/matrix/GSE65858_series_matrix.txt.gz"
GPL_URL    <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/suppl/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"

# MD5 of the UNCOMPRESSED files used for the manuscript (downloaded 17 June 2026,
# re-verified identical on 9 November 2026).
SERIES_MD5 <- "8cd09a480c1b4e904d250461f613ed20"
GPL_MD5    <- "3ce00b33f57be39902e5387116b52da5"

dir.create("data_raw", showWarnings = FALSE)

gunzip_file <- function(gz, out) {
  con_in  <- gzfile(gz, "rb")
  con_out <- file(out, "wb")
  on.exit({ close(con_in); close(con_out) })
  repeat {
    chunk <- readBin(con_in, "raw", n = 1e7)
    if (!length(chunk)) break
    writeBin(chunk, con_out)
  }
  invisible(out)
}

fetch_and_check <- function(url, gz, txt, expected_md5) {
  if (!file.exists(txt)) {
    if (!file.exists(gz)) {
      message("Downloading ", basename(gz))
      download.file(url, gz, mode = "wb", quiet = TRUE)
    }
    gunzip_file(gz, txt)
  }
  observed <- unname(tools::md5sum(txt))
  if (!identical(observed, expected_md5)) {
    stop("Checksum mismatch for ", txt, "\n  expected ", expected_md5,
         "\n  observed ", observed,
         "\nGEO may have updated this file; the published results were produced with the expected version.")
  }
  message("Checksum OK: ", basename(txt))
}

series_txt <- "data_raw/GSE65858_series_matrix.txt"
fetch_and_check(SERIES_URL, "data_raw/GSE65858_series_matrix.txt.gz", series_txt, SERIES_MD5)
gpl_txt <- "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt"
fetch_and_check(GPL_URL, "data_raw/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz", gpl_txt, GPL_MD5)

# ------------------------------------------------------------
# Parse the series matrix
# ------------------------------------------------------------
gse_txt <- readLines(series_txt)
strip_q <- function(x) gsub('^"|"$', "", x)

sample_line <- function(prefix) {
  ln <- grep(paste0("^", prefix, "\t"), gse_txt, value = TRUE)
  strip_q(strsplit(ln, "\t", fixed = TRUE)[[1]][-1])
}
geo_accession <- sample_line("!Sample_geo_accession")

# "key: value" characteristics rows -> named columns
char_lines <- grep("^!Sample_characteristics_ch1\t", gse_txt, value = TRUE)
char_vals <- lapply(char_lines, function(ln) strip_q(strsplit(ln, "\t", fixed = TRUE)[[1]][-1]))
keys <- vapply(char_vals, function(v) sub(":.*$", "", v[1]), character(1))
chr  <- setNames(lapply(char_vals, function(v) sub("^[^:]*: ?", "", v)), keys)
stopifnot(length(geo_accession) == 270, all(lengths(chr) == 270))

# GEO stores missing values in text fields as the literal string "NA". These are
# kept as-is for character variables (downstream scripts convert them to missing);
# numeric and logical variables are converted to numeric/logical with NA.
to_na   <- function(x) { x[x %in% c("NA", "")] <- NA; x }
to_lgl  <- function(x) as.logical(to_na(x))
to_num  <- function(x) suppressWarnings(as.numeric(to_na(x)))
prefix  <- function(p, x) paste0(p, x)

pheno <- tibble::tibble(
  sample_geo_accession      = geo_accession,
  sample_id                 = chr$ID,
  sample_age                = to_num(chr$age),
  sample_gender             = chr$gender,
  sample_smoking            = chr$smoking,
  sample_packyears          = chr$packyears,
  sample_alcohol            = chr$alcohol,
  sample_tumour_type        = chr$tumor_type,
  sample_tumour_site        = chr$tumor_site,
  sample_uicc_stage         = chr$uicc_stage,
  sample_t_category         = prefix("T", chr$t_category),
  sample_n_category         = prefix("N", chr$n_category),
  sample_distant_metastasis = to_num(chr$distant_metastasis),
  sample_treatment          = chr$treatment,
  sample_hpv_dna            = chr$hpv_dna,
  sample_hpv16_dna_rna      = chr$hpv16_dna_rna,
  os                        = to_num(chr$os),
  os_event                  = to_lgl(chr$os_event),
  sample_pfs                = to_num(chr$pfs),
  sample_pfs_event          = to_lgl(chr$pfs_event),
  sample_consensus_cluster  = chr$consensus_cluster
)

# ------------------------------------------------------------
# Parse the probe-level expression table and match samples
# ------------------------------------------------------------
start <- grep("^!series_matrix_table_begin", gse_txt)
end   <- grep("^!series_matrix_table_end", gse_txt)
expr <- read.delim(
  text = gse_txt[(start + 1):(end - 1)],
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

common_ids <- intersect(colnames(expr)[-1], pheno$sample_geo_accession)
expr2 <- expr[, c(colnames(expr)[1], common_ids)]
pheno3 <- pheno |>
  filter(sample_geo_accession %in% common_ids) |>
  arrange(match(sample_geo_accession, common_ids))

stopifnot(
  nrow(pheno3) == 270,
  all(pheno3$sample_geo_accession == colnames(expr2)[-1]),
  !anyDuplicated(expr2$ID_REF)
)

saveRDS(pheno3, "pheno3_matched.rds")
saveRDS(expr2, "expr2_probe_matrix.rds")
message("Wrote pheno3_matched.rds (", nrow(pheno3), " samples) and expr2_probe_matrix.rds (",
        nrow(expr2), " probes)")
