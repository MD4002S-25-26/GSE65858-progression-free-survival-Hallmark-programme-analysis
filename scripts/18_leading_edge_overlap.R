# ============================================================
# 18_leading_edge_overlap.R
#
# Post-dissertation analysis for the manuscript.
# Question: are the 23 significant Hallmark sets (BH q < 0.05 in the
# adjusted ranked fgsea analysis) independent, or do they share
# leading-edge genes?
#
# Pre-specified rules (fixed before inspecting any overlap results):
#   - significant sets: padj < 0.05 in the saved dissertation fgsea output
#   - overlap metric: Jaccard index between leading-edge gene sets
#   - clustering: average-linkage hierarchical clustering on 1 - Jaccard,
#     tree cut at height 0.8 (i.e. mean Jaccard >= 0.2 between merged groups)
#   - recurrent genes: leading-edge genes present in >= 50% of a cluster's sets
#
# Run from the repository root, after scripts 01-05 have created
# results_v2/GSE65858_v2_fgsea_hallmark.rds and ..._hallmark_gene_sets.rds.
# Outputs are written to results_v2/leading_edge/.
# No statistical re-fitting is performed; the dissertation fgsea
# object is used unchanged.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

IN_DIR  <- "results_v2"
OUT_DIR <- file.path("results_v2", "leading_edge")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

H_CUT <- 0.8

fg <- readRDS(file.path(IN_DIR, "GSE65858_v2_fgsea_hallmark.rds"))
hallmark <- readRDS(file.path(IN_DIR, "GSE65858_v2_hallmark_gene_sets.rds"))

sig <- fg |> filter(padj < 0.05) |> arrange(NES)
stopifnot(nrow(sig) == 23, sum(sig$NES < 0) == 13, sum(sig$NES > 0) == 10)

le <- setNames(sig$leadingEdge, sig$pathway)
short <- function(x) tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", x))))
labs <- setNames(short(names(le)), names(le))

# ---- 1. leading-edge summary ----------------------------------
summ <- sig |>
  transmute(
    pathway = short(pathway),
    direction = ifelse(NES < 0, "PFS-event-free at 2 years", "PFS event within 2 years"),
    NES = round(NES, 2), padj = signif(padj, 3),
    set_size_in_ranking = size,
    leading_edge_n = lengths(leadingEdge),
    leading_edge_pct_of_set = round(100 * lengths(leadingEdge) / size, 1)
  )
write.csv(summ, file.path(OUT_DIR, "leading_edge_summary.csv"), row.names = FALSE)

# ---- 2. pairwise Jaccard (leading edge and full set) ----------
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
pairwise <- function(lst) {
  n <- length(lst)
  m <- matrix(NA_real_, n, n, dimnames = list(names(lst), names(lst)))
  for (i in seq_len(n)) for (j in seq_len(n)) m[i, j] <- jac(lst[[i]], lst[[j]])
  m
}
J_le   <- pairwise(le)
J_full <- pairwise(hallmark[names(le)])
write.csv(round(J_le, 3), file.path(OUT_DIR, "pairwise_jaccard_leading_edge.csv"))
write.csv(round(J_full, 3), file.path(OUT_DIR, "pairwise_jaccard_full_sets.csv"))

# ---- 3. clustering --------------------------------------------
hc <- hclust(as.dist(1 - J_le), method = "average")
cl <- cutree(hc, h = H_CUT)
# renumber clusters by size (largest first), singletons last
tab <- sort(table(cl), decreasing = TRUE)
cl <- setNames(match(cl, as.integer(names(tab))), names(cl))

cluster_tbl <- data.frame(
  pathway = short(names(cl)),
  cluster = unname(cl),
  direction = ifelse(sig$NES[match(names(cl), sig$pathway)] < 0,
                     "PFS-event-free at 2 years", "PFS event within 2 years"),
  NES = round(sig$NES[match(names(cl), sig$pathway)], 2)
) |> arrange(cluster, NES)
write.csv(cluster_tbl, file.path(OUT_DIR, "clusters.csv"), row.names = FALSE)

# ---- 4. per-cluster redundancy and recurrent genes -------------
cluster_stats <- lapply(sort(unique(cl)), function(k) {
  members <- names(cl)[cl == k]
  sets <- le[members]
  allg <- unlist(sets)
  freq <- table(allg)
  within <- if (length(members) > 1) mean(J_le[members, members][upper.tri(diag(length(members)))]) else NA
  data.frame(
    cluster = k, n_sets = length(members),
    members = paste(short(members), collapse = "; "),
    directions = paste(unique(ifelse(sig$NES[match(members, sig$pathway)] < 0, "event-free", "early-event")), collapse = "/"),
    union_leading_edge_genes = length(unique(allg)),
    genes_in_ge_half_of_sets = sum(freq >= max(2, ceiling(length(members) / 2))),
    mean_within_cluster_jaccard = round(within, 3)
  )
}) |> bind_rows()
write.csv(cluster_stats, file.path(OUT_DIR, "cluster_summary.csv"), row.names = FALSE)

rec <- lapply(sort(unique(cl)), function(k) {
  members <- names(cl)[cl == k]
  if (length(members) < 2) return(NULL)
  freq <- sort(table(unlist(le[members])), decreasing = TRUE)
  keep <- freq[freq >= max(2, ceiling(length(members) / 2))]
  if (!length(keep)) return(NULL)
  data.frame(cluster = k, gene = names(keep), n_sets_with_gene = as.integer(keep),
             n_sets_in_cluster = length(members))
}) |> bind_rows()
write.csv(rec, file.path(OUT_DIR, "recurrent_core_genes.csv"), row.names = FALSE)

# ---- 5. overall redundancy ------------------------------------
all_inc <- unlist(le)
overall <- data.frame(
  metric = c("significant sets", "leading-edge gene incidences (sum over sets)",
             "unique leading-edge genes", "mean incidences per unique gene",
             "genes in >=3 sets", "genes in >=5 sets",
             "median pairwise Jaccard, leading edge (all pairs)",
             "median pairwise Jaccard, leading edge (same direction)",
             "median pairwise Jaccard, leading edge (opposite directions)",
             "number of clusters at cut h = 0.8"),
  value = c(length(le), length(all_inc), length(unique(all_inc)),
            round(length(all_inc) / length(unique(all_inc)), 2),
            sum(table(all_inc) >= 3), sum(table(all_inc) >= 5),
            round(median(J_le[upper.tri(J_le)]), 3),
            round(median(J_le[upper.tri(J_le) & outer(sig$NES < 0, sig$NES < 0, "==")]), 3),
            round(median(J_le[outer(sig$NES < 0, sig$NES > 0, "&")]), 3),
            length(unique(cl)))
)
write.csv(overall, file.path(OUT_DIR, "overall_redundancy.csv"), row.names = FALSE)
print(overall); print(cluster_stats[, c("cluster","n_sets","directions","union_leading_edge_genes","mean_within_cluster_jaccard")])

# ---- 6. figures -----------------------------------------------
ord <- hc$order
lv <- names(le)[ord]
long <- expand.grid(a = lv, b = lv, stringsAsFactors = FALSE)
long$J <- mapply(function(x, y) J_le[x, y], long$a, long$b)
long$a <- factor(long$a, levels = lv); long$b <- factor(long$b, levels = rev(lv))
dir_col <- ifelse(sig$NES[match(lv, sig$pathway)] < 0, "#1b6ca8", "#b5472f")

p <- ggplot(long, aes(a, b, fill = J)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#222222", limits = c(0, 1), name = "Jaccard\n(leading edge)") +
  scale_x_discrete(labels = labs[lv]) + scale_y_discrete(labels = labs[rev(lv)]) +
  labs(x = NULL, y = NULL,
       caption = "Axis label colour: blue = enriched toward event-free at 2 years; red = enriched toward PFS event within 2 years") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = dir_col),
        axis.text.y = element_text(colour = rev(dir_col)),
        panel.grid = element_blank(), plot.caption = element_text(size = 8, hjust = 0))
ggsave(file.path(OUT_DIR, "leading_edge_jaccard_heatmap.png"), p, width = 9, height = 8, dpi = 300)

png(file.path(OUT_DIR, "leading_edge_dendrogram.png"), width = 2400, height = 1600, res = 300)
par(mar = c(1, 4, 2, 1))
plot(hc, labels = labs[hc$labels], hang = -1, cex = 0.7, main = "Average-linkage clustering of leading-edge sets (1 - Jaccard)", xlab = "", sub = "")
abline(h = H_CUT, lty = 2, col = "grey40")
dev.off()

# ---- 7. descriptive extras (no new tests) ----------------------
ix <- which(upper.tri(J_le), arr.ind = TRUE)
pairs <- data.frame(
  set_a = short(rownames(J_le)[ix[, 1]]), set_b = short(colnames(J_le)[ix[, 2]]),
  jaccard_leading_edge = round(J_le[ix], 3),
  shared_genes = mapply(function(i, j) length(intersect(le[[i]], le[[j]])), ix[, 1], ix[, 2]),
  jaccard_full_set = round(J_full[ix], 3),
  same_direction = (sig$NES[ix[, 1]] < 0) == (sig$NES[ix[, 2]] < 0)
) |> arrange(desc(jaccard_leading_edge))
write.csv(pairs, file.path(OUT_DIR, "all_pairs_ranked.csv"), row.names = FALSE)
print(head(pairs, 10))

gene_freq <- sort(table(unlist(le)), decreasing = TRUE)
g5 <- data.frame(gene = names(gene_freq[gene_freq >= 5]), n_sets = as.integer(gene_freq[gene_freq >= 5]))
g5$sets <- vapply(g5$gene, function(g) paste(short(names(le)[vapply(le, function(x) g %in% x, logical(1))]), collapse = "; "), character(1))
write.csv(g5, file.path(OUT_DIR, "genes_in_5plus_sets.csv"), row.names = FALSE)
print(g5)
cat("median Jaccard, full Hallmark sets (same 23):", round(median(J_full[upper.tri(J_full)]), 3), "\n")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nDone. Outputs in", normalizePath(OUT_DIR), "\n")
