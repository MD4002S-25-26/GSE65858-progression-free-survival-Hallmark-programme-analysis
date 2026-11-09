# ============================================================
# 08_pipeline_audit.R
#
# Purpose:
# Summarise the corrected GSE65858 analysis and identify
# proportional-hazards diagnostic flags before final reporting.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Load final corrected result objects
# ------------------------------------------------------------
pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

de_unadj <- readRDS(
  "results_v2/GSE65858_v2_limma_unadjusted.rds"
)

de_adj <- readRDS(
  "results_v2/GSE65858_v2_limma_adjusted.rds"
)

fgsea_results <- readRDS(
  "results_v2/GSE65858_v2_fgsea_hallmark.rds"
)

adjusted_score_results <- readRDS(
  "results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.rds"
)

cox_results <- readRDS(
  "results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.rds"
)

# ------------------------------------------------------------
# 2. Endpoint summary
# ------------------------------------------------------------
cat("\nCORRECTED TWO-YEAR PFS GROUPS\n")

print(
  table(
    pheno_analysis$pfs_2yr_group,
    useNA = "ifany"
  )
)

cat(
  "Total cohort:",
  nrow(pheno_analysis),
  "\n"
)

cat(
  "Total PFS events:",
  sum(
    pheno_analysis$sample_pfs_event == 1,
    na.rm = TRUE
  ),
  "\n"
)

# ------------------------------------------------------------
# 3. Gene-level result summary
# ------------------------------------------------------------
cat("\nGENE-LEVEL ANALYSIS\n")

cat(
  "Genes tested:",
  nrow(de_adj),
  "\n"
)

cat(
  "Unadjusted genes at BH FDR < 0.05:",
  sum(de_unadj$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Adjusted genes at BH FDR < 0.05:",
  sum(de_adj$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------
# 4. Ranked enrichment summary
# ------------------------------------------------------------
cat("\nHALLMARK FGSEA\n")

cat(
  "Pathways tested:",
  nrow(fgsea_results),
  "\n"
)

cat(
  "Pathways at BH FDR < 0.05:",
  sum(fgsea_results$padj < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Pathways at BH FDR < 0.25:",
  sum(fgsea_results$padj < 0.25, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------
# 5. Sample-level score summary
# ------------------------------------------------------------
cat("\nSAMPLE-LEVEL PATHWAY SCORES\n")

cat(
  "Adjusted pathways at raw p < 0.05:",
  sum(
    adjusted_score_results$p_value < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Adjusted pathways at BH FDR < 0.05:",
  sum(
    adjusted_score_results$p_adj_BH < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

# ------------------------------------------------------------
# 6. Cox summary
# ------------------------------------------------------------
cat("\nADJUSTED COX MODELS\n")

cat(
  "Models fitted:",
  nrow(cox_results),
  "\n"
)

cat(
  "Complete-case sample size:",
  unique(cox_results$n),
  "\n"
)

cat(
  "Complete-case PFS events:",
  unique(cox_results$events),
  "\n"
)

cat(
  "Pathways at raw p < 0.05:",
  sum(cox_results$p_value < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Pathways at BH FDR < 0.05:",
  sum(cox_results$p_adj_BH < 0.05, na.rm = TRUE),
  "\n"
)

# ------------------------------------------------------------
# 7. Proportional-hazards diagnostic flags
# ------------------------------------------------------------
ph_flags <- cox_results |>
  mutate(
    pathway_PH_q_BH = p.adjust(
      pathway_PH_p,
      method = "BH"
    )
  ) |>
  filter(
    pathway_PH_p < 0.05
  ) |>
  arrange(
    pathway_PH_p
  ) |>
  select(
    pathway,
    hazard_ratio,
    p_value,
    p_adj_BH,
    pathway_PH_p,
    pathway_PH_q_BH,
    global_PH_p
  )

cat("\nPATHWAY-SPECIFIC PH FLAGS\n")

print(
  ph_flags,
  n = Inf,
  width = Inf
)

write.csv(
  ph_flags,
  "results_v2/GSE65858_v2_cox_PH_flags.csv",
  row.names = FALSE
)

cat(
  "\nPipeline audit completed successfully.\n"
)