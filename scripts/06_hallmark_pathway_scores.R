# ============================================================
# 06_hallmark_pathway_scores.R
#
# Purpose:
# Calculate sample-level scores for all 50 MSigDB Hallmark
# pathways and test their association with corrected two-year
# PFS landmark status.
#
# Scores are calculated by:
#   1. z-standardising each gene across all 270 samples
#   2. averaging z-scores for genes belonging to each pathway
#
# All 50 pathways are tested to avoid selecting pathways only
# because they appeared prominent in fgsea.
# ============================================================

library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Create output folder if required
# ------------------------------------------------------------
dir.create(
  "results_v2",
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------
expr_gene <- readRDS(
  "expr_gene_level.rds"
)

pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

hallmark_list <- readRDS(
  "results_v2/GSE65858_v2_hallmark_gene_sets.rds"
)

stopifnot(length(hallmark_list) == 50)

# ------------------------------------------------------------
# 3. Construct gene-expression matrix
# ------------------------------------------------------------
expr_mat <- as.matrix(
  expr_gene[, -1]
)

rownames(expr_mat) <- expr_gene$Symbol

storage.mode(expr_mat) <- "numeric"

cat(
  "Expression matrix dimensions:",
  nrow(expr_mat), "genes x",
  ncol(expr_mat), "samples\n"
)

# Confirm sample IDs match phenotype IDs
stopifnot(
  all(
    colnames(expr_mat) %in%
      pheno_analysis$sample_geo_accession
  )
)

# ------------------------------------------------------------
# 4. Z-standardise every gene across all 270 samples
#
# Each gene therefore has:
#   mean approximately 0
#   standard deviation approximately 1
# ------------------------------------------------------------
gene_means <- rowMeans(
  expr_mat,
  na.rm = TRUE
)

gene_sds <- apply(
  expr_mat,
  1,
  sd,
  na.rm = TRUE
)

valid_genes <- (
  is.finite(gene_means) &
    is.finite(gene_sds) &
    gene_sds > 0
)

cat(
  "Genes retained for z-standardisation:",
  sum(valid_genes),
  "\n"
)

cat(
  "Genes excluded because of zero or invalid variance:",
  sum(!valid_genes),
  "\n"
)

expr_valid <- expr_mat[
  valid_genes,
  ,
  drop = FALSE
]

z_expr <- sweep(
  expr_valid,
  1,
  gene_means[valid_genes],
  FUN = "-"
)

z_expr <- sweep(
  z_expr,
  1,
  gene_sds[valid_genes],
  FUN = "/"
)

stopifnot(
  all(is.finite(z_expr))
)

# ------------------------------------------------------------
# 5. Record pathway gene coverage
# ------------------------------------------------------------
pathway_coverage <- data.frame(
  pathway = names(hallmark_list),
  genes_in_msigdb = lengths(hallmark_list),
  genes_present = vapply(
    hallmark_list,
    function(genes) {
      length(
        intersect(
          genes,
          rownames(z_expr)
        )
      )
    },
    FUN.VALUE = integer(1)
  ),
  stringsAsFactors = FALSE
) |>
  mutate(
    coverage_percent =
      100 * genes_present / genes_in_msigdb
  ) |>
  arrange(genes_present)

print(pathway_coverage)

cat(
  "Minimum pathway gene coverage:",
  min(pathway_coverage$genes_present),
  "\n"
)

stopifnot(
  all(pathway_coverage$genes_present >= 15)
)

write.csv(
  pathway_coverage,
  "results_v2/GSE65858_v2_hallmark_pathway_coverage.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Calculate all 50 sample-level pathway scores
#
# Pathway score =
# mean z-standardised expression across available pathway genes
# ------------------------------------------------------------
score_list <- lapply(
  names(hallmark_list),
  function(pathway_name) {
    
    pathway_genes <- intersect(
      hallmark_list[[pathway_name]],
      rownames(z_expr)
    )
    
    colMeans(
      z_expr[
        pathway_genes,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
)

score_matrix <- do.call(
  cbind,
  score_list
)

colnames(score_matrix) <- names(hallmark_list)

rownames(score_matrix) <- colnames(z_expr)

stopifnot(
  nrow(score_matrix) == 270,
  ncol(score_matrix) == 50
)

cat(
  "Pathway score matrix dimensions:",
  nrow(score_matrix), "samples x",
  ncol(score_matrix), "pathways\n"
)

# ------------------------------------------------------------
# 7. Save sample-level scores
# ------------------------------------------------------------
score_wide <- data.frame(
  sample_geo_accession = rownames(score_matrix),
  score_matrix,
  check.names = FALSE,
  row.names = NULL
)

score_long <- score_wide |>
  pivot_longer(
    cols = -sample_geo_accession,
    names_to = "pathway",
    values_to = "pathway_score"
  )

saveRDS(
  score_wide,
  "results_v2/GSE65858_v2_hallmark_pathway_scores.rds"
)

write.csv(
  score_wide,
  "results_v2/GSE65858_v2_hallmark_pathway_scores_wide.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Add corrected phenotype information
# ------------------------------------------------------------
score_analysis <- score_long |>
  left_join(
    pheno_analysis |>
      select(
        sample_geo_accession,
        pfs_2yr_group,
        sample_pfs,
        sample_pfs_event,
        sample_uicc_stage,
        sample_hpv16_dna_rna
      ),
    by = "sample_geo_accession"
  ) |>
  mutate(
    group = case_when(
      pfs_2yr_group ==
        "PFS event within 2 years" ~
        "Early_PFS_event",
      
      pfs_2yr_group ==
        "PFS-event-free at 2 years" ~
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

stopifnot(
  !any(
    is.na(score_analysis$pfs_2yr_group)
  )
)

# ------------------------------------------------------------
# 9. Unadjusted landmark-group comparisons
#
# Wilcoxon rank-sum tests are used because pathway-score
# distributions are not assumed to be normally distributed.
#
# BH correction is applied across all 50 pathways.
# ------------------------------------------------------------
score_binary <- score_analysis |>
  filter(
    !is.na(group)
  )

wilcoxon_results <- score_binary |>
  group_by(pathway) |>
  group_modify(
    ~ {
      
      early_scores <- .x$pathway_score[
        .x$group == "Early_PFS_event"
      ]
      
      event_free_scores <- .x$pathway_score[
        .x$group == "Event_free_2y"
      ]
      
      test_result <- wilcox.test(
        pathway_score ~ group,
        data = .x,
        exact = FALSE
      )
      
      data.frame(
        n_early = length(early_scores),
        n_event_free = length(event_free_scores),
        
        median_early =
          median(early_scores, na.rm = TRUE),
        
        median_event_free =
          median(event_free_scores, na.rm = TRUE),
        
        median_difference =
          median(early_scores, na.rm = TRUE) -
          median(event_free_scores, na.rm = TRUE),
        
        p_value = test_result$p.value
      )
    }
  ) |>
  ungroup() |>
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    score_direction = case_when(
      median_difference > 0 ~
        "Higher in early PFS-event group",
      
      median_difference < 0 ~
        "Higher in PFS-event-free group",
      
      TRUE ~
        "No difference"
    )
  ) |>
  arrange(
    p_adj_BH,
    p_value
  )

write.csv(
  wilcoxon_results,
  "results_v2/GSE65858_v2_hallmark_scores_wilcoxon_all50.csv",
  row.names = FALSE
)

cat(
  "\nPathways significant in unadjusted Wilcoxon analysis:",
  sum(
    wilcoxon_results$p_adj_BH < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

print(
  wilcoxon_results |>
    select(
      pathway,
      median_early,
      median_event_free,
      median_difference,
      p_value,
      p_adj_BH,
      score_direction
    ) |>
    head(15)
)

# ------------------------------------------------------------
# 10. Prepare adjusted analysis cohort
#
# Matches the covariate handling used in adjusted limma:
#   - UICC stage
#   - HPV16 DNA/RNA status
# ------------------------------------------------------------
score_adjusted <- score_binary |>
  filter(
    !is.na(sample_uicc_stage),
    !is.na(sample_hpv16_dna_rna),
    as.character(sample_hpv16_dna_rna) != "NA"
  ) |>
  mutate(
    sample_uicc_stage = factor(
      sample_uicc_stage,
      levels = c(
        "I",
        "II",
        "III",
        "IVA",
        "IVB",
        "IVC"
      )
    ),
    
    sample_hpv16_dna_rna = factor(
      sample_hpv16_dna_rna
    ),
    
    sample_hpv16_dna_rna = relevel(
      sample_hpv16_dna_rna,
      ref = "DNA-"
    ),
    
    group = factor(
      group,
      levels = c(
        "Event_free_2y",
        "Early_PFS_event"
      )
    )
  )

cat(
  "\nAdjusted pathway-score model sample size per pathway:",
  n_distinct(score_adjusted$sample_geo_accession),
  "\n"
)

# ------------------------------------------------------------
# 11. Adjusted linear models for all 50 pathways
#
# Positive beta:
#   higher score in early PFS-event group
#
# Negative beta:
#   higher score in PFS-event-free group
# ------------------------------------------------------------
adjusted_score_results <- score_adjusted |>
  group_by(pathway) |>
  group_modify(
    ~ {
      
      model_data <- droplevels(.x)
      
      fit <- lm(
        pathway_score ~
          sample_uicc_stage +
          sample_hpv16_dna_rna +
          group,
        data = model_data
      )
      
      coefficient_table <-
        summary(fit)$coefficients
      
      term_name <- "groupEarly_PFS_event"
      
      estimate <-
        coefficient_table[
          term_name,
          "Estimate"
        ]
      
      standard_error <-
        coefficient_table[
          term_name,
          "Std. Error"
        ]
      
      degrees_freedom <-
        df.residual(fit)
      
      critical_value <- qt(
        0.975,
        df = degrees_freedom
      )
      
      data.frame(
        n = nobs(fit),
        beta_early_vs_event_free = estimate,
        standard_error = standard_error,
        confidence_interval_low =
          estimate -
          critical_value * standard_error,
        confidence_interval_high =
          estimate +
          critical_value * standard_error,
        t_statistic =
          coefficient_table[
            term_name,
            "t value"
          ],
        p_value =
          coefficient_table[
            term_name,
            "Pr(>|t|)"
          ]
      )
    }
  ) |>
  ungroup() |>
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    score_direction = case_when(
      beta_early_vs_event_free > 0 ~
        "Higher in early PFS-event group",
      
      beta_early_vs_event_free < 0 ~
        "Higher in PFS-event-free group",
      
      TRUE ~
        "No difference"
    )
  ) |>
  arrange(
    p_adj_BH,
    p_value
  )

write.csv(
  adjusted_score_results,
  "results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.csv",
  row.names = FALSE
)

saveRDS(
  adjusted_score_results,
  "results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.rds"
)

cat(
  "\nPathways significant after adjusted-model BH correction:",
  sum(
    adjusted_score_results$p_adj_BH < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nTop adjusted sample-level pathway results:\n"
)

print(
  adjusted_score_results |>
    select(
      pathway,
      n,
      beta_early_vs_event_free,
      confidence_interval_low,
      confidence_interval_high,
      p_value,
      p_adj_BH,
      score_direction
    ) |>
    head(15)
)

cat(
  "\nAll-50 Hallmark pathway scoring completed successfully.\n"
)
pathways_of_initial_interest <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

adjusted_score_results |>
  filter(
    pathway %in%
      pathways_of_initial_interest
  ) |>
  select(
    pathway,
    beta_early_vs_event_free,
    confidence_interval_low,
    confidence_interval_high,
    p_value,
    p_adj_BH,
    score_direction
  )
options(tibble.width = Inf)

adjusted_score_results |>
  filter(
    pathway %in% pathways_of_initial_interest
  ) |>
  arrange(p_value) |>
  select(
    pathway,
    beta_early_vs_event_free,
    confidence_interval_low,
    confidence_interval_high,
    p_value,
    p_adj_BH,
    score_direction
  ) |>
  print(
    n = Inf,
    width = Inf
  )