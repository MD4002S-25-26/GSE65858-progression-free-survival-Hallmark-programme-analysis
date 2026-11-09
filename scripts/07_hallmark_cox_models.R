# ============================================================
# 07_hallmark_cox_models.R
#
# Purpose:
# Test associations between all 50 sample-level Hallmark scores
# and progression-free survival using adjusted Cox models.
#
# Each pathway is tested in a separate model adjusted for:
#   - UICC stage
#   - HPV16 DNA/RNA status
#
# Hazard ratios are reported per 1 SD increase in pathway score.
# ============================================================

library(dplyr)
library(tidyr)
library(survival)

# ------------------------------------------------------------
# 1. Load phenotype and pathway scores
# ------------------------------------------------------------
pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

score_wide <- readRDS(
  "results_v2/GSE65858_v2_hallmark_pathway_scores.rds"
)

# ------------------------------------------------------------
# 2. Convert pathway-score table to long format
# ------------------------------------------------------------
score_long <- score_wide |>
  pivot_longer(
    cols = -sample_geo_accession,
    names_to = "pathway",
    values_to = "pathway_score"
  )

# ------------------------------------------------------------
# 3. Add survival and clinical variables
# ------------------------------------------------------------
cox_data <- score_long |>
  left_join(
    pheno_analysis |>
      select(
        sample_geo_accession,
        sample_pfs,
        sample_pfs_event,
        sample_uicc_stage,
        sample_hpv16_dna_rna
      ),
    by = "sample_geo_accession"
  ) |>
  mutate(
    sample_pfs = as.numeric(sample_pfs),
    sample_pfs_event = as.numeric(sample_pfs_event),
    
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
    )
  )

# Basic outcome checks
cat(
  "Total cohort samples:",
  n_distinct(cox_data$sample_geo_accession),
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

stopifnot(
  all(
    cox_data$sample_pfs > 0,
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 4. Fit one adjusted Cox model per Hallmark pathway
# ------------------------------------------------------------
cox_results <- cox_data |>
  group_by(pathway) |>
  group_modify(
    ~ {
      
      model_data <- .x |>
        filter(
          !is.na(sample_pfs),
          !is.na(sample_pfs_event),
          !is.na(sample_uicc_stage),
          !is.na(sample_hpv16_dna_rna),
          as.character(sample_hpv16_dna_rna) != "NA",
          !is.na(pathway_score)
        ) |>
        droplevels()
      
      model_data$sample_hpv16_dna_rna <-
        relevel(
          model_data$sample_hpv16_dna_rna,
          ref = "DNA-"
        )
      
      # Standardise the pathway score within the model cohort
      model_data$pathway_score_z <-
        as.numeric(
          scale(model_data$pathway_score)
        )
      
      stopifnot(
        all(is.finite(model_data$pathway_score_z))
      )
      
      fit <- coxph(
        Surv(
          sample_pfs,
          sample_pfs_event
        ) ~
          pathway_score_z +
          sample_uicc_stage +
          sample_hpv16_dna_rna,
        data = model_data,
        ties = "efron",
        x = TRUE,
        model = TRUE
      )
      
      fit_summary <- summary(fit)
      
      coefficient_row <-
        fit_summary$coefficients[
          "pathway_score_z",
          ,
          drop = FALSE
        ]
      
      confidence_row <-
        fit_summary$conf.int[
          "pathway_score_z",
          ,
          drop = FALSE
        ]
      
      # Test proportional-hazards assumptions
      ph_test <- cox.zph(
        fit,
        transform = "km"
      )
      
      ph_table <- ph_test$table
      
      pathway_ph_p <- if (
        "pathway_score_z" %in% rownames(ph_table)
      ) {
        ph_table[
          "pathway_score_z",
          "p"
        ]
      } else {
        NA_real_
      }
      
      global_ph_p <- if (
        "GLOBAL" %in% rownames(ph_table)
      ) {
        ph_table[
          "GLOBAL",
          "p"
        ]
      } else {
        NA_real_
      }
      
      data.frame(
        n = nrow(model_data),
        
        events = sum(
          model_data$sample_pfs_event == 1
        ),
        
        beta = coefficient_row[
          1,
          "coef"
        ],
        
        standard_error = coefficient_row[
          1,
          "se(coef)"
        ],
        
        hazard_ratio = confidence_row[
          1,
          "exp(coef)"
        ],
        
        confidence_interval_low = confidence_row[
          1,
          "lower .95"
        ],
        
        confidence_interval_high = confidence_row[
          1,
          "upper .95"
        ],
        
        z_statistic = coefficient_row[
          1,
          "z"
        ],
        
        p_value = coefficient_row[
          1,
          "Pr(>|z|)"
        ],
        
        pathway_PH_p = pathway_ph_p,
        
        global_PH_p = global_ph_p
      )
    }
  ) |>
  ungroup() |>
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    
    hazard_direction = case_when(
      hazard_ratio > 1 ~
        "Higher pathway score associated with higher PFS-event hazard",
      
      hazard_ratio < 1 ~
        "Higher pathway score associated with lower PFS-event hazard",
      
      TRUE ~
        "No directional association"
    )
  ) |>
  arrange(
    p_adj_BH,
    p_value
  )

# ------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------
write.csv(
  cox_results,
  "results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.csv",
  row.names = FALSE
)

saveRDS(
  cox_results,
  "results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.rds"
)

# ------------------------------------------------------------
# 6. Print audit summary
# ------------------------------------------------------------
cat(
  "\nCox models fitted:",
  nrow(cox_results),
  "\n"
)

cat(
  "Pathways with raw Cox p < 0.05:",
  sum(
    cox_results$p_value < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Pathways with BH-adjusted Cox p < 0.05:",
  sum(
    cox_results$p_adj_BH < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Pathway-score PH tests with p < 0.05:",
  sum(
    cox_results$pathway_PH_p < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Models with global PH-test p < 0.05:",
  sum(
    cox_results$global_PH_p < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nTop adjusted Cox results:\n"
)

cox_results |>
  select(
    pathway,
    n,
    events,
    hazard_ratio,
    confidence_interval_low,
    confidence_interval_high,
    p_value,
    p_adj_BH,
    pathway_PH_p,
    global_PH_p
  ) |>
  head(15) |>
  print(
    n = 15,
    width = Inf
  )

cat(
  "\nAll-50 adjusted Cox analysis completed successfully.\n"
)