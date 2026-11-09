# ============================================================
# 10_table2_integrated_pathway_summary.R
#
# Purpose:
# Create a concise main-text table integrating selected
# Hallmark results across fgsea, sample-level score models,
# adjusted Cox models and PH diagnostics.
# ============================================================

library(dplyr)
library(flextable)
library(officer)

dir.create(
  "tables_v2",
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load result objects
# ------------------------------------------------------------
fgsea_results <- readRDS(
  "results_v2/GSE65858_v2_fgsea_hallmark.rds"
)

score_results <- readRDS(
  "results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.rds"
)

cox_results <- readRDS(
  "results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.rds"
)

# ------------------------------------------------------------
# 2. Select representative pathways
# ------------------------------------------------------------
selected_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

pathway_labels <- c(
  "HALLMARK_E2F_TARGETS" =
    "E2F targets",
  
  "HALLMARK_MYC_TARGETS_V1" =
    "MYC targets V1",
  
  "HALLMARK_G2M_CHECKPOINT" =
    "G2M checkpoint",
  
  "HALLMARK_ALLOGRAFT_REJECTION" =
    "Allograft rejection",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" =
    "TNFα signalling via NF-κB",
  
  "HALLMARK_INFLAMMATORY_RESPONSE" =
    "Inflammatory response",
  
  "HALLMARK_KRAS_SIGNALING_UP" =
    "KRAS signalling up",
  
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" =
    "IL6–JAK–STAT3 signalling"
)

# ------------------------------------------------------------
# 3. Formatting functions
# ------------------------------------------------------------
format_probability <- function(x) {
  vapply(
    x,
    function(value) {
      if (is.na(value)) {
        "—"
      } else if (value < 0.001) {
        "<0.001"
      } else {
        sprintf("%.3f", value)
      }
    },
    FUN.VALUE = character(1)
  )
}

# ------------------------------------------------------------
# 4. Join results
# ------------------------------------------------------------
table2_data <- fgsea_results |>
  filter(
    pathway %in% selected_pathways
  ) |>
  select(
    pathway,
    NES,
    fgsea_q = padj
  ) |>
  left_join(
    score_results |>
      select(
        pathway,
        score_beta =
          beta_early_vs_event_free,
        score_ci_low =
          confidence_interval_low,
        score_ci_high =
          confidence_interval_high,
        score_q =
          p_adj_BH
      ),
    by = "pathway"
  ) |>
  left_join(
    cox_results |>
      select(
        pathway,
        cox_hr =
          hazard_ratio,
        cox_ci_low =
          confidence_interval_low,
        cox_ci_high =
          confidence_interval_high,
        cox_q =
          p_adj_BH,
        ph_p =
          pathway_PH_p
      ),
    by = "pathway"
  ) |>
  mutate(
    display_order =
      match(
        pathway,
        selected_pathways
      ),
    
    Pathway =
      unname(
        pathway_labels[pathway]
      ),
    
    `Ranked enrichment\nNES (BH q)` =
      paste0(
        sprintf("%.2f", NES),
        " (",
        format_probability(fgsea_q),
        ")"
      ),
    
    `Adjusted score estimate\n(95% CI; BH q)` =
      paste0(
        sprintf("%.3f", score_beta),
        " (",
        sprintf("%.3f", score_ci_low),
        " to ",
        sprintf("%.3f", score_ci_high),
        "); ",
        format_probability(score_q)
      ),
    
    `Adjusted Cox HR\n(95% CI; BH q)` =
      paste0(
        sprintf("%.2f", cox_hr),
        " (",
        sprintf("%.2f", cox_ci_low),
        " to ",
        sprintf("%.2f", cox_ci_high),
        "); ",
        format_probability(cox_q)
      ),
    
    `PH p` =
      format_probability(ph_p)
  ) |>
  arrange(
    display_order
  ) |>
  select(
    Pathway,
    `Ranked enrichment\nNES (BH q)`,
    `Adjusted score estimate\n(95% CI; BH q)`,
    `Adjusted Cox HR\n(95% CI; BH q)`,
    `PH p`
  )

print(table2_data)

write.csv(
  table2_data,
  "tables_v2/Table2_integrated_pathway_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Create formatted table
# ------------------------------------------------------------
table2_flex <- flextable(
  table2_data
) |>
  theme_booktabs() |>
  font(
    fontname = "Arial",
    part = "all"
  ) |>
  fontsize(
    size = 8.5,
    part = "all"
  ) |>
  bold(
    part = "header"
  ) |>
  align(
    j = 1,
    align = "left",
    part = "all"
  ) |>
  align(
    j = 2:5,
    align = "center",
    part = "all"
  ) |>
  valign(
    valign = "center",
    part = "all"
  ) |>
  width(
    j = 1,
    width = 1.45
  ) |>
  width(
    j = 2,
    width = 0.95
  ) |>
  width(
    j = 3,
    width = 1.70
  ) |>
  width(
    j = 4,
    width = 1.65
  ) |>
  width(
    j = 5,
    width = 0.55
  ) |>
  add_footer_lines(
    values = paste0(
      "Negative NES or score estimates indicate higher values ",
      "in the PFS-event-free group; positive values indicate ",
      "higher values in the early PFS-event group. Cox hazard ",
      "ratios are reported per one-standard-deviation increase ",
      "in pathway score. BH q-values were calculated across all ",
      "50 Hallmark pathways within each analysis. PH p-values ",
      "are from pathway-specific proportional-hazards diagnostics."
    )
  ) |>
  fontsize(
    i = 1,
    size = 7.5,
    part = "footer"
  )

# ------------------------------------------------------------
# 6. Save Word table
# ------------------------------------------------------------
save_as_docx(
  "Table 2" = table2_flex,
  path =
    "tables_v2/Table2_integrated_pathway_summary.docx"
)

cat(
  "\nTable 2 completed successfully.\n"
)

cat(
  "Saved:",
  "tables_v2/Table2_integrated_pathway_summary.docx\n"
)