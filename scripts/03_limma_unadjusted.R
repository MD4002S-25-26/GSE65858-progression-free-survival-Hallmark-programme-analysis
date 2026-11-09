# ============================================================
# 03_limma_unadjusted.R
#
# Purpose:
# Compare gene expression between:
#   1. PFS event within 2 years
#   2. PFS-event-free at 2 years
#
# Patients censored before 2 years are excluded.
# ============================================================

library(dplyr)
library(limma)

# ------------------------------------------------------------
# 1. Load corrected phenotype and verified gene matrix
# ------------------------------------------------------------
pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

expr_gene <- readRDS(
  "expr_gene_level.rds"
)

# ------------------------------------------------------------
# 2. Create the corrected binary comparison cohort
# ------------------------------------------------------------
pheno_binary <- pheno_analysis |>
  filter(
    pfs_2yr_group %in% c(
      "PFS event within 2 years",
      "PFS-event-free at 2 years"
    )
  ) |>
  mutate(
    group = case_when(
      pfs_2yr_group == "PFS event within 2 years" ~
        "Early_PFS_event",
      
      pfs_2yr_group == "PFS-event-free at 2 years" ~
        "Event_free_2y",
      
      TRUE ~ NA_character_
    ),
    
    group = factor(
      group,
      levels = c(
        "Event_free_2y",
        "Early_PFS_event"
      )
    )
  )

# Expected:
# Event_free_2y = 135
# Early_PFS_event = 106
print(table(pheno_binary$group, useNA = "ifany"))

# ------------------------------------------------------------
# 3. Build gene-expression matrix
# ------------------------------------------------------------
expr_mat <- as.matrix(expr_gene[, -1])
rownames(expr_mat) <- expr_gene$Symbol

# Select and order samples to exactly match phenotype
expr_binary <- expr_mat[
  ,
  pheno_binary$sample_geo_accession,
  drop = FALSE
]

# ------------------------------------------------------------
# 4. Essential alignment checks
# ------------------------------------------------------------
print(dim(expr_binary))

alignment_correct <- all(
  colnames(expr_binary) ==
    pheno_binary$sample_geo_accession
)

cat(
  "Expression and phenotype alignment:",
  alignment_correct,
  "\n"
)

stopifnot(alignment_correct)
stopifnot(ncol(expr_binary) == nrow(pheno_binary))

# ------------------------------------------------------------
# 5. Construct unadjusted design matrix
#
# Event_free_2y is the reference group.
# Therefore:
# positive logFC = higher in Early_PFS_event
# negative logFC = higher in Event_free_2y
# ------------------------------------------------------------
design_unadj <- model.matrix(
  ~ group,
  data = pheno_binary
)

print(colnames(design_unadj))

# ------------------------------------------------------------
# 6. Run limma
# ------------------------------------------------------------
fit_unadj <- lmFit(
  expr_binary,
  design_unadj
)

fit_unadj <- eBayes(fit_unadj)

de_unadj <- topTable(
  fit_unadj,
  coef = "groupEarly_PFS_event",
  number = Inf,
  sort.by = "P"
)

# Add gene symbols as a normal column
de_unadj$Symbol <- rownames(de_unadj)

de_unadj <- de_unadj |>
  select(
    Symbol,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B
  )

# ------------------------------------------------------------
# 7. Inspect results
# ------------------------------------------------------------
print(head(de_unadj, 20))

cat(
  "Genes with BH-adjusted p < 0.05:",
  sum(de_unadj$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Genes with absolute log2 fold change > 1:",
  sum(abs(de_unadj$logFC) > 1, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------
# 8. Save corrected outputs
# ------------------------------------------------------------
write.csv(
  de_unadj,
  "results_v2/GSE65858_v2_limma_unadjusted.csv",
  row.names = FALSE
)

saveRDS(
  de_unadj,
  "results_v2/GSE65858_v2_limma_unadjusted.rds"
)

saveRDS(
  pheno_binary,
  "results_v2/pheno_binary_corrected_2yr.rds"
)

saveRDS(
  expr_binary,
  "results_v2/expr_binary_corrected_2yr.rds"
)