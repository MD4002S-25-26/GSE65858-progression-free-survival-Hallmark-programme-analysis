source("scripts/01_prepare_phenotype_and_endpoint.R")
source("scripts/02_probe_to_gene_mapping.R")
source("scripts/03_limma_unadjusted.R")
source("scripts/04_limma_adjusted.R")
parse(
  file = "scripts/04_limma_adjusted.R"
)
file.exists(
  "results_v2/GSE65858_v2_limma_adjusted.csv"
)

file.exists(
  "results_v2/GSE65858_v2_limma_adjusted.rds"
)

file.exists(
  "results_v2/pheno_adjusted_corrected_2yr.rds"
)

file.exists(
  "results_v2/expr_adjusted_corrected_2yr.rds"
)
# ============================================================
# 05_fgsea_hallmark.R
#
# Purpose:
# Perform Hallmark ranked gene-set enrichment analysis using
# moderated t-statistics from the corrected adjusted limma model.
#
# Ranking direction:
#   Positive t / NES = enriched toward PFS event within 2 years
#   Negative t / NES = enriched toward PFS-event-free at 2 years
# ============================================================

library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)

# ------------------------------------------------------------
# 1. Create output folders if required
# ------------------------------------------------------------
dir.create(
  "results_v2",
  showWarnings = FALSE
)

dir.create(
  "figures_v2",
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load corrected adjusted limma results
# ------------------------------------------------------------
de_adj <- readRDS(
  "results_v2/GSE65858_v2_limma_adjusted.rds"
)

required_columns <- c(
  "Symbol",
  "t",
  "P.Value",
  "adj.P.Val"
)

stopifnot(
  all(required_columns %in% colnames(de_adj))
)

# ------------------------------------------------------------
# 3. Construct ranked gene statistic
#
# The moderated t-statistic is used rather than p-value because
# it retains both effect direction and relative evidence.
# ------------------------------------------------------------
rank_data <- de_adj |>
  filter(
    !is.na(Symbol),
    Symbol != "",
    is.finite(t)
  ) |>
  distinct(
    Symbol,
    .keep_all = TRUE
  )

ranked_stats <- rank_data$t
names(ranked_stats) <- rank_data$Symbol

ranked_stats <- sort(
  ranked_stats,
  decreasing = TRUE
)

cat(
  "Number of ranked genes:",
  length(ranked_stats),
  "\n"
)

stopifnot(
  length(ranked_stats) > 0,
  !anyDuplicated(names(ranked_stats)),
  all(is.finite(ranked_stats))
)

# Save ranked statistic for reproducibility
saveRDS(
  ranked_stats,
  "results_v2/GSE65858_v2_adjusted_ranked_t_statistics.rds"
)

# ------------------------------------------------------------
# 4. Obtain the MSigDB Hallmark gene sets
# ------------------------------------------------------------
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
) |>
  select(
    gs_name,
    gene_symbol
  ) |>
  distinct() |>
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  )

hallmark_list <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

cat(
  "Number of Hallmark pathways:",
  length(hallmark_list),
  "\n"
)

stopifnot(length(hallmark_list) == 50)

saveRDS(
  hallmark_list,
  "results_v2/GSE65858_v2_hallmark_gene_sets.rds"
)

# ------------------------------------------------------------
# 5. Run ranked enrichment analysis
# ------------------------------------------------------------
set.seed(123)

fgsea_results <- fgseaMultilevel(
  pathways = hallmark_list,
  stats = ranked_stats,
  minSize = 15,
  maxSize = 500,
  eps = 0
) |>
  as.data.frame() |>
  arrange(
    padj,
    desc(abs(NES))
  )

# Add an interpretable direction label
fgsea_results <- fgsea_results |>
  mutate(
    enrichment_direction = case_when(
      NES > 0 ~ "PFS event within 2 years",
      NES < 0 ~ "PFS-event-free at 2 years",
      TRUE ~ "No direction"
    )
  )

# ------------------------------------------------------------
# 6. Inspect results
# ------------------------------------------------------------
cat(
  "Hallmark pathways tested:",
  nrow(fgsea_results),
  "\n"
)

cat(
  "Pathways with BH-adjusted p < 0.05:",
  sum(fgsea_results$padj < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Pathways with BH-adjusted p < 0.25:",
  sum(fgsea_results$padj < 0.25, na.rm = TRUE),
  "\n"
)

cat("\nTop corrected Hallmark results:\n")

print(
  fgsea_results |>
    select(
      pathway,
      size,
      NES,
      pval,
      padj,
      enrichment_direction
    ) |>
    head(20)
)

cat("\nTop pathways enriched toward early PFS events:\n")

print(
  fgsea_results |>
    filter(NES > 0) |>
    arrange(padj, desc(NES)) |>
    select(pathway, size, NES, pval, padj) |>
    head(10)
)

cat("\nTop pathways enriched toward PFS-event-free status:\n")

print(
  fgsea_results |>
    filter(NES < 0) |>
    arrange(padj, NES) |>
    select(pathway, size, NES, pval, padj) |>
    head(10)
)

# ------------------------------------------------------------
# 7. Save complete results
#
# Convert the leading-edge list column to text for the CSV.
# Retain the original list column in the RDS file.
# ------------------------------------------------------------
saveRDS(
  fgsea_results,
  "results_v2/GSE65858_v2_fgsea_hallmark.rds"
)

fgsea_csv <- fgsea_results |>
  mutate(
    leadingEdge = vapply(
      leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  )

write.csv(
  fgsea_csv,
  "results_v2/GSE65858_v2_fgsea_hallmark.csv",
  row.names = FALSE
)

write.csv(
  fgsea_csv |>
    filter(padj < 0.05),
  "results_v2/GSE65858_v2_fgsea_hallmark_FDR05.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Create a provisional top-pathway plot
#
# This is an audit plot, not necessarily the final dissertation
# figure. Final presentation will be refined later.
# ------------------------------------------------------------
top_positive <- fgsea_results |>
  filter(NES > 0) |>
  arrange(padj, desc(NES)) |>
  slice_head(n = 5)

top_negative <- fgsea_results |>
  filter(NES < 0) |>
  arrange(padj, NES) |>
  slice_head(n = 5)

plot_data <- bind_rows(
  top_positive,
  top_negative
) |>
  mutate(
    pathway_label = gsub(
      "^HALLMARK_",
      "",
      pathway
    ),
    pathway_label = gsub(
      "_",
      " ",
      pathway_label
    ),
    pathway_label = factor(
      pathway_label,
      levels = rev(pathway_label)
    )
  )

fgsea_plot <- ggplot(
  plot_data,
  aes(
    x = NES,
    y = pathway_label,
    fill = enrichment_direction
  )
) +
  geom_col() +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  labs(
    title = "Hallmark programmes associated with two-year PFS status",
    x = "Normalised enrichment score",
    y = NULL,
    fill = "Enrichment direction"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(fgsea_plot)

ggsave(
  filename = "figures_v2/GSE65858_v2_fgsea_hallmark_top10.png",
  plot = fgsea_plot,
  width = 9,
  height = 6,
  units = "in",
  dpi = 300
)

cat(
  "\nCorrected Hallmark enrichment completed successfully.\n"
)
expr_mat <- as.matrix(
  readRDS("expr_gene_level.rds")[, -1]
)

gene_sd <- apply(
  expr_mat,
  1,
  sd,
  na.rm = TRUE
)

cat(
  "Genes with valid non-zero SD:",
  sum(is.finite(gene_sd) & gene_sd > 0),
  "\n"
)

cat(
  "Genes excluded for zero/invalid SD:",
  sum(!is.finite(gene_sd) | gene_sd <= 0),
  "\n"
)