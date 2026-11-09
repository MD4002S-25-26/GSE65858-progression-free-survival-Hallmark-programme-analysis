file.exists("GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt")
# ============================================================
# 02_probe_to_gene_mapping.R
#
# Purpose:
# Map the GSE65858 Illumina probe-level expression matrix to
# human gene symbols using GPL10558 platform annotation.
#
# Multiple probes mapping to the same gene symbol are collapsed
# by taking the mean expression value across probes.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Define input files
# ------------------------------------------------------------
expression_file <- "expr2_probe_matrix.rds"

annotation_file <-
  "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt"

# Confirm that both required files exist
stopifnot(file.exists(expression_file))
stopifnot(file.exists(annotation_file))

# ------------------------------------------------------------
# 2. Load probe-level expression data
# ------------------------------------------------------------
expr2 <- readRDS(expression_file)

cat(
  "Original probe-level matrix dimensions:",
  nrow(expr2), "rows x", ncol(expr2), "columns\n"
)

stopifnot("ID_REF" %in% colnames(expr2))

# Expression sample columns are all columns other than ID_REF
sample_columns <- setdiff(
  colnames(expr2),
  "ID_REF"
)

cat(
  "Number of expression samples:",
  length(sample_columns),
  "\n"
)

# ------------------------------------------------------------
# 3. Read GPL10558 annotation file
# ------------------------------------------------------------
annot_txt <- readLines(
  annotation_file,
  warn = FALSE
)

# Locate the start of the probe annotation table
probe_start <- grep(
  "^\\[Probes\\]$",
  annot_txt
)

stopifnot(length(probe_start) == 1)

annot_lines <- annot_txt[
  (probe_start + 1):length(annot_txt)
]

annot <- read.delim(
  text = annot_lines,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

cat(
  "Full annotation dimensions:",
  nrow(annot), "rows x", ncol(annot), "columns\n"
)

# Confirm required annotation columns exist
required_annotation_columns <- c(
  "Species",
  "Probe_Id",
  "Symbol",
  "Entrez_Gene_ID",
  "ILMN_Gene"
)

stopifnot(
  all(
    required_annotation_columns %in%
      colnames(annot)
  )
)

# ------------------------------------------------------------
# 4. Retain human probes
# ------------------------------------------------------------
annot_human <- annot |>
  filter(Species == "Homo sapiens")

cat(
  "Human annotation probes:",
  nrow(annot_human),
  "\n"
)

# Check how many expression probes match the annotation
matched_probe_count <- sum(
  expr2$ID_REF %in% annot_human$Probe_Id
)

cat(
  "Expression probes matched to GPL10558:",
  matched_probe_count,
  "of",
  nrow(expr2),
  "\n"
)

stopifnot(matched_probe_count > 0)

# ------------------------------------------------------------
# 5. Join probe expression to gene annotation
# ------------------------------------------------------------
expr_annot <- expr2 |>
  rename(Probe_Id = ID_REF) |>
  left_join(
    annot_human |>
      select(
        Probe_Id,
        Symbol,
        Entrez_Gene_ID,
        ILMN_Gene
      ),
    by = "Probe_Id"
  )

cat(
  "Annotated probe matrix dimensions:",
  nrow(expr_annot), "rows x",
  ncol(expr_annot), "columns\n"
)

# ------------------------------------------------------------
# 6. Remove probes without a valid gene symbol
# ------------------------------------------------------------
expr_annot_valid <- expr_annot |>
  filter(
    !is.na(Symbol),
    trimws(Symbol) != ""
  )

cat(
  "Probes with valid gene symbols:",
  nrow(expr_annot_valid),
  "\n"
)

cat(
  "Unique valid gene symbols:",
  n_distinct(expr_annot_valid$Symbol),
  "\n"
)

# ------------------------------------------------------------
# 7. Collapse multiple probes per gene
#
# For genes represented by more than one probe, mean expression
# is calculated separately for every sample.
# ------------------------------------------------------------
expr_gene <- expr_annot_valid |>
  select(
    Symbol,
    all_of(sample_columns)
  ) |>
  group_by(Symbol) |>
  summarise(
    across(
      all_of(sample_columns),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat(
  "Final gene-level matrix dimensions:",
  nrow(expr_gene), "genes x",
  length(sample_columns), "samples\n"
)

# ------------------------------------------------------------
# 8. Quality-control checks
# ------------------------------------------------------------
stopifnot(!anyDuplicated(expr_gene$Symbol))

stopifnot(
  identical(
    colnames(expr_gene)[-1],
    sample_columns
  )
)

stopifnot(
  nrow(expr_gene) > 0
)

print(head(expr_gene[, 1:5]))

# ------------------------------------------------------------
# 9. Save final gene-level expression matrix
# ------------------------------------------------------------
saveRDS(
  expr_gene,
  "expr_gene_level.rds"
)

write.csv(
  expr_gene,
  "results_v2/GSE65858_gene_level_expression.csv",
  row.names = FALSE
)

cat(
  "\nProbe-to-gene mapping completed successfully.\n"
)

cat(
  "Saved: expr_gene_level.rds\n"
)

cat(
  "Saved: results_v2/GSE65858_gene_level_expression.csv\n"
)
expr_gene_check <- readRDS(
  "expr_gene_level.rds"
)

dim(expr_gene_check)
