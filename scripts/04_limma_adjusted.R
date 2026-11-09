# ============================================================
# 04_limma_adjusted.R
#
# Purpose:
# Compare gene expression between:
#   1. PFS event within 2 years
#   2. PFS-event-free at 2 years
#
# Adjusted for:
#   - UICC stage
#   - HPV16 DNA/RNA status
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
# 2. Create corrected binary cohort
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
    ),
    
    sample_uicc_stage = factor(
      sample_uicc_stage,
      levels = c("I", "II", "III", "IVA", "IVB", "IVC")
    ),
    
    sample_hpv16_dna_rna = factor(
      sample_hpv16_dna_rna
    )
  )

# ------------------------------------------------------------
# 3. Remove samples missing adjustment variables
# ------------------------------------------------------------
pheno_adj <- pheno_binary |>
  filter(
    !is.na(sample_uicc_stage),
    !is.na(sample_hpv16_dna_rna),
    as.character(sample_hpv16_dna_rna) != "NA"
  ) |>
  droplevels() |>
  mutate(
    group = relevel(
      group,
      ref = "Event_free_2y"
    ),
    
    sample_hpv16_dna_rna = relevel(
      sample_hpv16_dna_rna,
      ref = "DNA-"
    )
  )

# Check sample counts retained
print(table(pheno_adj$group, useNA = "ifany"))
print(table(pheno_adj$sample_uicc_stage, useNA = "ifany"))
print(table(pheno_adj$sample_hpv16_dna_rna, useNA = "ifany"))

cat(
  "Adjusted model sample size:",
  nrow(pheno_adj),
  "\n"
)

cat(
  "Samples removed for missing covariates:",
  nrow(pheno_binary) - nrow(pheno_adj),
  "\n"
)

# ------------------------------------------------------------
# 4. Build and align expression matrix
# ------------------------------------------------------------
expr_mat <- as.matrix(expr_gene[, -1])
rownames(expr_mat) <- expr_gene$Symbol

expr_adj <- expr_mat[
  ,
  pheno_adj$sample_geo_accession,
  drop = FALSE
]

alignment_correct <- all(
  colnames(expr_adj) ==
    pheno_adj$sample_geo_accession
)

cat(
  "Expression and phenotype alignment:",
  alignment_correct,
  "\n"
)

print(dim(expr_adj))

stopifnot(alignment_correct)
stopifnot(ncol(expr_adj) == nrow(pheno_adj))

# ------------------------------------------------------------
# 5. Construct adjusted design matrix
#
# Event_free_2y is the progression-group reference.
# DNA- is the HPV reference.
#
# Positive group logFC:
#   higher expression in Early_PFS_event
#
# Negative group logFC:
#   higher expression in Event_free_2y
# ------------------------------------------------------------
design_adj <- model.matrix(
  ~ sample_uicc_stage +
    sample_hpv16_dna_rna +
    group,
  data = pheno_adj
)

print(colnames(design_adj))
print(dim(design_adj))

# ------------------------------------------------------------
# 6. Run adjusted limma
# ------------------------------------------------------------
fit_adj <- lmFit(
  expr_adj,
  design_adj
)

fit_adj <- eBayes(fit_adj)

de_adj <- topTable(
  fit_adj,
  coef = "groupEarly_PFS_event",
  number = Inf,
  sort.by = "P"
)

de_adj$Symbol <- rownames(de_adj)

de_adj <- de_adj |>
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
# 7. Inspect corrected adjusted results
# ------------------------------------------------------------
print(head(de_adj, 20))

cat(
  "Genes with BH-adjusted p < 0.05:",
  sum(de_adj$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Genes with absolute log2 fold change > 1:",
  sum(abs(de_adj$logFC) > 1, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------
# 8. Save corrected adjusted outputs
# ------------------------------------------------------------
write.csv(
  de_adj,
  "results_v2/GSE65858_v2_limma_adjusted.csv",
  row.names = FALSE
)

saveRDS(
  de_adj,
  "results_v2/GSE65858_v2_limma_adjusted.rds"
)

saveRDS(
  pheno_adj,
  "results_v2/pheno_adjusted_corrected_2yr.rds"
)

saveRDS(
  expr_adj,
  "results_v2/expr_adjusted_corrected_2yr.rds"
)

cat(
  "\nCorrected adjusted limma analysis completed successfully.\n"
)