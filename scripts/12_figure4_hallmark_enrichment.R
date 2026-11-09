# ============================================================
# 12_figure4_hallmark_enrichment.R
#
# Purpose:
# Create a diverging horizontal bar plot of Hallmark pathways
# reaching BH-adjusted q < 0.05 in the corrected fgsea analysis.
#
# Negative NES:
#   Enriched toward PFS-event-free status at two years.
#
# Positive NES:
#   Enriched toward a PFS event within two years.
# ============================================================

library(dplyr)
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

event_free_label <- "PFS-event-free at two years"
early_event_label <- "PFS event within two years"

# ------------------------------------------------------------
# 1. Load corrected fgsea results
# ------------------------------------------------------------
fgsea_results <- readRDS(
  "results_v2/GSE65858_v2_fgsea_hallmark.rds"
)

required_columns <- c(
  "pathway",
  "NES",
  "pval",
  "padj"
)

stopifnot(
  all(required_columns %in% colnames(fgsea_results))
)

# ------------------------------------------------------------
# 2. Function to produce readable pathway names
# ------------------------------------------------------------
clean_pathway_name <- function(x) {
  cleaned <- sub("^HALLMARK_", "", x)
  cleaned <- gsub("_", " ", cleaned)
  cleaned <- tools::toTitleCase(tolower(cleaned))

  replacements <- c(
    "E2f" = "E2F",
    "Myc" = "MYC",
    "G2m" = "G2M",
    "Kras" = "KRAS",
    "Dna" = "DNA",
    "Mtorc1" = "mTORC1",
    "Il2" = "IL2",
    "Il6" = "IL6",
    "Jak" = "JAK",
    "Stat3" = "STAT3",
    "Stat5" = "STAT5"
  )

  for (old_text in names(replacements)) {
    cleaned <- gsub(
      old_text,
      replacements[[old_text]],
      cleaned,
      fixed = TRUE
    )
  }

  ascii_hyphens(cleaned)
}

# ------------------------------------------------------------
# 3. Retain significant Hallmark pathways
# ------------------------------------------------------------
plot_data <- fgsea_results |>
  filter(
    !is.na(padj),
    padj < 0.05,
    is.finite(NES)
  ) |>
  arrange(NES) |>
  mutate(
    pathway_label = clean_pathway_name(pathway),
    direction = case_when(
      NES < 0 ~ event_free_label,
      NES > 0 ~ early_event_label,
      TRUE ~ "No direction"
    )
  ) |>
  mutate(
    pathway_label = case_when(
      pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~
        "TNF-alpha signalling via NF-kappaB",
      pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING" ~
        "IL6-JAK-STAT3 signalling",
      pathway == "HALLMARK_IL2_STAT5_SIGNALING" ~
        "IL2-STAT5 signalling",
      pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE" ~
        "Interferon gamma response",
      pathway == "HALLMARK_MYC_TARGETS_V1" ~
        "MYC targets V1",
      pathway == "HALLMARK_MYC_TARGETS_V2" ~
        "MYC targets V2",
      pathway == "HALLMARK_E2F_TARGETS" ~
        "E2F targets",
      pathway == "HALLMARK_G2M_CHECKPOINT" ~
        "G2M checkpoint",
      pathway == "HALLMARK_KRAS_SIGNALING_UP" ~
        "KRAS signalling up",
      TRUE ~ pathway_label
    ),
    pathway_label = ascii_hyphens(pathway_label),
    direction = ascii_hyphens(direction),
    NES_label = sprintf("%.2f", NES),
    label_x = case_when(
      NES < 0 ~ NES - 0.07,
      NES > 0 ~ NES + 0.07,
      TRUE ~ NES
    ),
    label_hjust = case_when(
      NES < 0 ~ 1,
      NES > 0 ~ 0,
      TRUE ~ 0.5
    )
  )

# ------------------------------------------------------------
# 4. Audit results
# ------------------------------------------------------------
cat("Hallmark pathways displayed:", nrow(plot_data), "\n")
cat("Negative NES pathways:", sum(plot_data$NES < 0), "\n")
cat("Positive NES pathways:", sum(plot_data$NES > 0), "\n")

stopifnot(nrow(plot_data) == 23)

plot_data |>
  select(pathway_label, NES, padj, direction) |>
  tibble::as_tibble() |>
  print(n = Inf, width = Inf)

# ------------------------------------------------------------
# 5. Set pathway ordering
# ------------------------------------------------------------
plot_data$pathway_label <- factor(
  plot_data$pathway_label,
  levels = rev(plot_data$pathway_label)
)

# ------------------------------------------------------------
# 6. Create diverging Hallmark enrichment plot
# ------------------------------------------------------------
hallmark_plot <- ggplot(
  plot_data,
  aes(
    x = NES,
    y = pathway_label,
    fill = direction
  )
) +
  geom_col(
    width = 0.74
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      x = label_x,
      label = NES_label,
      hjust = label_hjust
    ),
    size = 3.1,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "PFS-event-free at two years" = "#2C73B3",
      "PFS event within two years" = "#E69F00"
    ),
    breaks = c(
      event_free_label,
      early_event_label
    )
  ) +
  scale_x_continuous(
    limits = c(-3.6, 3.1),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    title = "Hallmark enrichment by two-year PFS status",
    subtitle = "Twenty-three of 50 Hallmark gene sets reached BH q < 0.05",
    x = paste0(
      "Normalised enrichment score\n",
      "Negative = PFS-event-free; ",
      "positive = PFS event within two years"
    ),
    y = NULL,
    fill = "Enrichment direction"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0
    ),
    plot.title.position = "plot",
    plot.subtitle = element_text(
      size = 10.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 9)
    ),
    axis.text.y = element_text(
      size = 9
    ),
    axis.text.x = element_text(
      size = 9
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    ),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(
      t = 12,
      r = 35,
      b = 12,
      l = 12
    )
  )

print(hallmark_plot)

# ------------------------------------------------------------
# 7. Save high-resolution outputs
# ------------------------------------------------------------
ggsave(
  filename = "figures_v2/Figure4_hallmark_fgsea_significant.png",
  plot = hallmark_plot,
  width = 7.6,
  height = 8.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "figures_v2/Figure4_hallmark_fgsea_significant.pdf",
  plot = hallmark_plot,
  width = 7.6,
  height = 8.8,
  units = "in"
)

cat("\nHallmark enrichment figure completed successfully.\n")