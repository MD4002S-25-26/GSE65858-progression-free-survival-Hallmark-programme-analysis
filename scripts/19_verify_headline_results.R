# ============================================================
# 19_verify_headline_results.R
#
# Checks that a fresh run reproduces the headline numbers reported in the
# manuscript. Run from the repository root after the pipeline has finished.
# Numeric comparisons use a small tolerance (floating-point differences
# between platforms are expected in the last digits).
# ============================================================

fail <- character()
check <- function(label, ok) { cat(sprintf("%-70s %s\n", label, if (isTRUE(ok)) "OK" else "FAILED")); if (!isTRUE(ok)) fail <<- c(fail, label) }
near <- function(a, b, tol = 0.01) isTRUE(abs(a - b) <= tol)

de_adj <- read.csv("results_v2/GSE65858_v2_limma_adjusted.csv")
check("22,001 genes ranked", nrow(de_adj) == 22001)
check("no gene with BH q < 0.05 (adjusted)", sum(de_adj$adj.P.Val < 0.05) == 0)
check("largest |log2FC| below 1", max(abs(de_adj$logFC)) < 1)
check("smallest q = 0.628 (AMIGO2)", near(min(de_adj$adj.P.Val), 0.628, 0.001) && de_adj$Symbol[which.min(de_adj$P.Value)] == "AMIGO2")

fg <- readRDS("results_v2/GSE65858_v2_fgsea_hallmark.rds")
check("23 Hallmark sets with q < 0.05", sum(fg$padj < 0.05) == 23)
check("13 negative and 10 positive NES among significant", sum(fg$padj < 0.05 & fg$NES < 0) == 13 && sum(fg$padj < 0.05 & fg$NES > 0) == 10)
check("E2F Targets NES = -2.89", near(fg$NES[fg$pathway == "HALLMARK_E2F_TARGETS"], -2.89))
check("MYC Targets V1 NES = -2.73", near(fg$NES[fg$pathway == "HALLMARK_MYC_TARGETS_V1"], -2.73))
check("Allograft Rejection NES = 2.40", near(fg$NES[fg$pathway == "HALLMARK_ALLOGRAFT_REJECTION"], 2.40))

sc <- read.csv("results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.csv")
check("adjusted score models: n = 222", all(sc$n == 222))
check("adjusted score models: 9 nominal, none with q < 0.05", sum(sc$p_value < 0.05) == 9 && sum(sc$p_adj_BH < 0.05) == 0)

cx <- read.csv("results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.csv")
check("Cox: 250 patients, 120 events", all(cx$n == 250) && all(cx$events == 120))
check("Cox: 4 nominal, none with q < 0.05", sum(cx$p_value < 0.05) == 4 && sum(cx$p_adj_BH < 0.05) == 0)
check("Cox: pathway-specific PH p < 0.05 in 11 models, no global test < 0.05", sum(cx$pathway_PH_p < 0.05) == 11 && sum(cx$global_PH_p < 0.05) == 0)

ov <- read.csv("results_v2/leading_edge/overall_redundancy.csv")
val <- function(m) as.numeric(ov$value[ov$metric == m])
check("leading edge: 23 sets, 1,533 entries, 1,128 unique genes", val("significant sets") == 23 && val("leading-edge gene incidences (sum over sets)") == 1533 && val("unique leading-edge genes") == 1128)
check("leading edge: 21 clusters at cut height 0.8", val("number of clusters at cut h = 0.8") == 21)

if (length(fail)) stop(length(fail), " check(s) failed: ", paste(fail, collapse = "; ")) else cat("\nAll headline checks passed.\n")
