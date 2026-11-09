# ============================================================
# 11_figure3_adjusted_volcano.R
#
# Purpose:
# Create a volcano plot for the corrected, stage- and
# HPV-adjusted limma comparison.
#
# Positive log2FC:
#   Higher expression among patients with a PFS event
#   within two years.
#
# Negative log2FC:
#   Higher expression among patients who remained
#   PFS-event-free at two years.
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

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

# ------------------------------------------------------------
# 1. Load adjusted limma results
# ------------------------------------------------------------
de_adj <- readRDS(
  "results_v2/GSE65858_v2_limma_adjusted.rds"
)

required_columns <- c(
  "Symbol",
  "logFC",
  "P.Value",
  "adj.P.Val"
)

stopifnot(
  all(required_columns %in% colnames(de_adj))
)

# ------------------------------------------------------------
# 2. Prepare plotting data
# ------------------------------------------------------------
volcano_data <- de_adj |>
  filter(
    !is.na(Symbol),
    Symbol != "",
    is.finite(logFC),
    is.finite(P.Value),
    P.Value > 0,
    is.finite(adj.P.Val)
  ) |>
  mutate(
    minus_log10_p = -log10(P.Value)
  )

cat("Genes plotted:", nrow(volcano_data), "\n")
cat(
  "Genes at BH q < 0.05:",
  sum(volcano_data$adj.P.Val < 0.05),
  "\n"
)
cat(
  "Genes with |log2FC| > 1:",
  sum(abs(volcano_data$logFC) > 1),
  "\n"
)

# ------------------------------------------------------------
# 3. Label the six genes with the smallest raw p-values
# ------------------------------------------------------------
label_data <- volcano_data |>
  arrange(P.Value) |>
  slice_head(n = 6)

print(
  label_data |>
    select(Symbol, logFC, P.Value, adj.P.Val)
)

# ------------------------------------------------------------
# 4. Calculate plotting limits
# ------------------------------------------------------------
x_limit <- max(
  abs(volcano_data$logFC),
  na.rm = TRUE
)

x_limit <- ceiling(x_limit * 12) / 10

y_limit <- max(
  volcano_data$minus_log10_p,
  na.rm = TRUE
)

figure_subtitle <- ascii_hyphens(
  "PFS event within two years versus PFS-event-free at two years"
)

summary_label <- ascii_hyphens(
  paste0(
    "0 genes at BH q < 0.05\n",
    "0 genes with |log2FC| > 1"
  )
)

# ------------------------------------------------------------
# 5. Create plot
# ------------------------------------------------------------
volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = logFC,
    y = minus_log10_p
  )
) +
  geom_point(
    size = 0.75,
    alpha = 0.22
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = Symbol),
    size = 3.2,
    seed = 123,
    max.overlaps = Inf,
    box.padding = 0.55,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.35
  ) +
  annotate(
    geom = "label",
    x = -0.96 * x_limit,
    y = 0.97 * y_limit,
    label = summary_label,
    hjust = 0,
    vjust = 1,
    size = 3.35,
    linewidth = 0.25
  ) +
  coord_cartesian(
    xlim = c(-x_limit, x_limit),
    clip = "off"
  ) +
  labs(
    title = "Adjusted gene-level differential expression",
    subtitle = figure_subtitle,
    x = expression(
      atop(
        log[2]~"fold change",
        "Positive values indicate higher expression with a PFS event within two years"
      )
    ),
    y = expression(
      -log[10]~"raw p-value"
    )
  ) +
  theme_classic(
    base_size = 11.5
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10.5,
      margin = margin(b = 8)
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )

print(volcano_plot)

# ------------------------------------------------------------
# 6. Save high-resolution outputs
# ------------------------------------------------------------
ggsave(
  filename = "figures_v2/Figure3_adjusted_limma_volcano.png",
  plot = volcano_plot,
  width = 7.2,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "figures_v2/Figure3_adjusted_limma_volcano.pdf",
  plot = volcano_plot,
  width = 7.2,
  height = 5.2,
  units = "in"
)

cat("\nAdjusted volcano plot completed successfully.\n")