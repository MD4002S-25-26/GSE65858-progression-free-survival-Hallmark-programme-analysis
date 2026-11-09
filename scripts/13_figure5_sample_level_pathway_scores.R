# ============================================================
# 13_figure5_sample_level_pathway_scores.R
#
# Purpose:
# Display sample-level pathway-score distributions for six
# representative Hallmark programmes.
#
# The plotted cohort matches the adjusted pathway-score models:
# 222 patients with complete UICC stage and HPV16 DNA/RNA data.
#
# Adjusted beta:
#   Positive = higher score in the early PFS-event group.
#   Negative = higher score in the PFS-event-free group.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

dir.create(
  "figures_v2",
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Helper: replace hidden Unicode dash characters with "-"
# ------------------------------------------------------------
ascii_hyphens <- function(x) {
  gsub(
    "[\u2010\u2011\u2012\u2013\u2014\u2212]",
    "-",
    x,
    perl = TRUE
  )
}

event_free_label <- "PFS-event-free\nat 2 years"
early_event_label <- "PFS event\nwithin 2 years"

# ------------------------------------------------------------
# 1. Load pathway scores, phenotype and adjusted results
# ------------------------------------------------------------
score_wide <- readRDS(
  "results_v2/GSE65858_v2_hallmark_pathway_scores.rds"
)

pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

adjusted_score_results <- readRDS(
  "results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.rds"
)

# ------------------------------------------------------------
# 2. Select representative pathways
# ------------------------------------------------------------
selected_pathways <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

pathway_labels <- c(
  "HALLMARK_MYC_TARGETS_V1" = "MYC targets V1",
  "HALLMARK_E2F_TARGETS" = "E2F targets",
  "HALLMARK_G2M_CHECKPOINT" = "G2M checkpoint",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammatory response",
  "HALLMARK_KRAS_SIGNALING_UP" = "KRAS signalling up",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6-JAK-STAT3 signalling"
)

pathway_labels <- ascii_hyphens(pathway_labels)

# ------------------------------------------------------------
# 3. Convert pathway scores to long format and add phenotype
# ------------------------------------------------------------
plot_data <- score_wide |>
  pivot_longer(
    cols = -sample_geo_accession,
    names_to = "pathway",
    values_to = "pathway_score"
  ) |>
  filter(
    pathway %in% selected_pathways
  ) |>
  left_join(
    pheno_analysis |>
      select(
        sample_geo_accession,
        pfs_2yr_group,
        sample_uicc_stage,
        sample_hpv16_dna_rna
      ),
    by = "sample_geo_accession"
  ) |>
  filter(
    pfs_2yr_group %in% c(
      "PFS-event-free at 2 years",
      "PFS event within 2 years"
    ),
    !is.na(sample_uicc_stage),
    !is.na(sample_hpv16_dna_rna),
    as.character(sample_hpv16_dna_rna) != "NA"
  ) |>
  mutate(
    group = case_when(
      pfs_2yr_group == "PFS-event-free at 2 years" ~
        event_free_label,
      pfs_2yr_group == "PFS event within 2 years" ~
        early_event_label,
      TRUE ~ NA_character_
    ),
    group = ascii_hyphens(group),
    group = factor(
      group,
      levels = c(
        event_free_label,
        early_event_label
      )
    ),
    pathway_label = unname(pathway_labels[pathway]),
    pathway_label = ascii_hyphens(pathway_label),
    pathway_label = factor(
      pathway_label,
      levels = unname(pathway_labels[selected_pathways])
    )
  )

cat(
  "Unique patients plotted:",
  n_distinct(plot_data$sample_geo_accession),
  "\n"
)

cat(
  "Pathways plotted:",
  n_distinct(plot_data$pathway),
  "\n"
)

print(
  with(
    plot_data |>
      distinct(sample_geo_accession, group),
    table(group)
  )
)

stopifnot(
  n_distinct(plot_data$sample_geo_accession) == 222,
  n_distinct(plot_data$pathway) == 6
)

# ------------------------------------------------------------
# 4. Prepare aligned statistical annotations
# ------------------------------------------------------------
annotation_data <- adjusted_score_results |>
  filter(
    pathway %in% selected_pathways
  ) |>
  mutate(
    pathway_label = unname(pathway_labels[pathway]),
    pathway_label = ascii_hyphens(pathway_label),
    pathway_label = factor(
      pathway_label,
      levels = unname(pathway_labels[selected_pathways])
    ),
    annotation_label = paste0(
      "Adjusted beta = ",
      sprintf("%.3f", beta_early_vs_event_free),
      "\n95% CI ",
      sprintf("%.3f", confidence_interval_low),
      " to ",
      sprintf("%.3f", confidence_interval_high),
      "\nBH q = ",
      sprintf("%.3f", p_adj_BH)
    ),
    annotation_label = ascii_hyphens(annotation_label),
    annotation_x = 1.5,
    annotation_y = 1.55
  ) |>
  select(
    pathway,
    pathway_label,
    annotation_label,
    annotation_x,
    annotation_y
  )

# ------------------------------------------------------------
# 5. Create the six-panel distribution figure
# ------------------------------------------------------------
pathway_plot <- ggplot(
  plot_data,
  aes(
    x = group,
    y = pathway_score,
    fill = group
  )
) +
  geom_violin(
    alpha = 0.25,
    width = 0.85,
    trim = FALSE,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.20,
    alpha = 0.75,
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  geom_jitter(
    aes(colour = group),
    width = 0.12,
    height = 0,
    alpha = 0.30,
    size = 0.75,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_text(
    data = annotation_data,
    aes(
      x = annotation_x,
      y = annotation_y,
      label = annotation_label
    ),
    inherit.aes = FALSE,
    size = 2.65,
    lineheight = 0.95,
    hjust = 0.5,
    vjust = 0.5
  ) +
  facet_wrap(
    ~ pathway_label,
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "PFS-event-free\nat 2 years" = "#1769AA",
      "PFS event\nwithin 2 years" = "#E69F00"
    ),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(
      "PFS-event-free\nat 2 years" = "#1769AA",
      "PFS event\nwithin 2 years" = "#E69F00"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.25)
    )
  ) +
  labs(
    title = "Sample-level Hallmark pathway scores by two-year PFS status",
    subtitle = paste0(
      "Adjusted for UICC stage and HPV16 DNA/RNA status; ",
      "no pathway survived BH correction\n",
      "across all 50 Hallmark tests"
    ),
    x = NULL,
    y = "Mean pathway gene z-score"
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 9.5,
      margin = margin(b = 10)
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    strip.background = element_blank(),
    axis.title.y = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 8.2
    ),
    panel.spacing = grid::unit(
      1.1,
      "lines"
    ),
    plot.margin = margin(
      10,
      12,
      10,
      10
    )
  )

print(pathway_plot)

# ------------------------------------------------------------
# 6. Save high-resolution outputs
# ------------------------------------------------------------
ggsave(
  filename = "figures_v2/Figure5_sample_level_pathway_scores.png",
  plot = pathway_plot,
  width = 8.2,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "figures_v2/Figure5_sample_level_pathway_scores.pdf",
  plot = pathway_plot,
  width = 8.2,
  height = 7.2,
  units = "in"
)

cat("\nSample-level pathway-score figure completed successfully.\n")