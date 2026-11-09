# ============================================================
# 21_treatment_sensitivity.R
#
# Post hoc sensitivity analysis: do the main results change when treatment is
# added to the adjustment set (UICC stage and HPV16 DNA/RNA status)?
#
# The analysis was specified before it was run; see
# verification/treatment_sensitivity_plan_2026-11-09.md. Settings are copied from
# scripts 04 (limma), 05 (fgsea), 06 (adjusted score models) and 07 (Cox models),
# with one added treatment term. Those scripts and their outputs are not changed.
#
# Treatment coding (fixed in advance):
#   primary   : multi-modality (reference) vs less than multi-modality
#               (single-modality and palliative combined)
#   secondary : the same models with palliative patients excluded
#
# Run from the repository root after scripts 01 to 07 have finished.
# Outputs: results_v2/treatment_sensitivity/
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(limma); library(fgsea); library(survival)
})

OUT <- "results_v2/treatment_sensitivity"
dir.create(OUT, showWarnings = FALSE)

ph        <- readRDS("pheno_analysis_corrected_landmark.rds")
pheno_adj <- readRDS("results_v2/pheno_adjusted_corrected_2yr.rds")   # the 222-patient cohort
expr_adj  <- readRDS("results_v2/expr_adjusted_corrected_2yr.rds")
hallmark  <- readRDS("results_v2/GSE65858_v2_hallmark_gene_sets.rds")
scores    <- readRDS("results_v2/GSE65858_v2_hallmark_pathway_scores.rds")
main_de   <- readRDS("results_v2/GSE65858_v2_limma_adjusted.rds")
main_fg   <- readRDS("results_v2/GSE65858_v2_fgsea_hallmark.rds")
main_sc   <- read.csv("results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.csv")
main_cx   <- read.csv("results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.csv")

# ---- treatment variable ---------------------------------------------------------------
stopifnot(!anyNA(ph$sample_treatment), all(ph$sample_treatment %in% c("mono", "multi", "palliative")))
trt <- ph |> select(sample_geo_accession, sample_treatment)
add_trt <- function(d, exclude_palliative) {
  d <- d |> select(-any_of("sample_treatment")) |> left_join(trt, by = "sample_geo_accession")
  stopifnot(!anyNA(d$sample_treatment))
  if (exclude_palliative) d <- d |> filter(sample_treatment != "palliative")
  d$trt <- factor(ifelse(d$sample_treatment == "multi", "multi", "less_than_multi"),
                  levels = c("multi", "less_than_multi"))
  d
}

# ---- (a) limma and ranked enrichment ---------------------------------------------------
run_limma_fgsea <- function(exclude_palliative, tag, with_trt = TRUE) {
  d <- add_trt(pheno_adj, exclude_palliative)
  d <- d |> mutate(group = relevel(factor(group), ref = "Event_free_2y"),
                   sample_hpv16_dna_rna = relevel(factor(sample_hpv16_dna_rna), ref = "DNA-")) |>
    droplevels()
  x <- expr_adj[, d$sample_geo_accession, drop = FALSE]
  stopifnot(all(colnames(x) == d$sample_geo_accession))
  design <- if (with_trt) model.matrix(~ sample_uicc_stage + sample_hpv16_dna_rna + trt + group, data = d) else
    model.matrix(~ sample_uicc_stage + sample_hpv16_dna_rna + group, data = d)
  fit <- eBayes(lmFit(x, design))
  de <- topTable(fit, coef = "groupEarly_PFS_event", number = Inf, sort.by = "P")
  de$Symbol <- rownames(de)
  de <- de |> select(Symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
  if (with_trt) write.csv(de, file.path(OUT, paste0("limma_adjusted_plus_treatment_", tag, ".csv")), row.names = FALSE)
  rk <- de |> filter(!is.na(Symbol), Symbol != "", is.finite(t)) |> distinct(Symbol, .keep_all = TRUE)
  ranked <- sort(setNames(rk$t, rk$Symbol), decreasing = TRUE)
  set.seed(123)
  fg <- fgseaMultilevel(pathways = hallmark, stats = ranked, minSize = 15, maxSize = 500, eps = 0) |>
    as.data.frame() |> arrange(padj, desc(abs(NES)))
  fgc <- fg |> mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))
  if (with_trt) write.csv(fgc, file.path(OUT, paste0("fgsea_plus_treatment_", tag, ".csv")), row.names = FALSE)
  list(n = nrow(d), de = de, fg = fg)
}

# ---- (b) sample-level score models ---------------------------------------------------
run_scores <- function(exclude_palliative, tag, with_trt = TRUE) {
  d <- add_trt(pheno_adj, exclude_palliative) |>
    mutate(group = relevel(factor(group), ref = "Event_free_2y"),
           sample_hpv16_dna_rna = relevel(factor(sample_hpv16_dna_rna), ref = "DNA-")) |>
    droplevels()
  paths <- setdiff(colnames(scores), "sample_geo_accession")
  res <- lapply(paths, function(p) {
    md <- d |> select(sample_geo_accession, sample_uicc_stage, sample_hpv16_dna_rna, trt, group) |>
      left_join(scores[, c("sample_geo_accession", p)] |> setNames(c("sample_geo_accession", "pathway_score")),
                by = "sample_geo_accession")
    fit <- if (with_trt) lm(pathway_score ~ sample_uicc_stage + sample_hpv16_dna_rna + trt + group, data = md) else
      lm(pathway_score ~ sample_uicc_stage + sample_hpv16_dna_rna + group, data = md)
    cf <- summary(fit)$coefficients["groupEarly_PFS_event", ]
    crit <- qt(0.975, df.residual(fit))
    data.frame(pathway = p, n = nobs(fit), beta_early_vs_event_free = cf[["Estimate"]],
               standard_error = cf[["Std. Error"]],
               confidence_interval_low = cf[["Estimate"]] - crit * cf[["Std. Error"]],
               confidence_interval_high = cf[["Estimate"]] + crit * cf[["Std. Error"]],
               t_statistic = cf[["t value"]], p_value = cf[["Pr(>|t|)"]])
  }) |> bind_rows() |> mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |> arrange(p_adj_BH, p_value)
  if (with_trt) write.csv(res, file.path(OUT, paste0("scores_adjusted_plus_treatment_", tag, ".csv")), row.names = FALSE)
  res
}

# ---- (c) Cox models --------------------------------------------------------------------
run_cox <- function(exclude_palliative, tag, with_trt = TRUE) {
  base <- ph |> mutate(sample_pfs = as.numeric(sample_pfs), sample_pfs_event = as.numeric(sample_pfs_event),
                       sample_uicc_stage = factor(sample_uicc_stage, levels = c("I", "II", "III", "IVA", "IVB", "IVC")),
                       sample_hpv16_dna_rna = factor(sample_hpv16_dna_rna)) |>
    filter(!is.na(sample_pfs), !is.na(sample_pfs_event), !is.na(sample_uicc_stage),
           !is.na(sample_hpv16_dna_rna), as.character(sample_hpv16_dna_rna) != "NA",
           sample_geo_accession %in% scores$sample_geo_accession) |>
    droplevels()
  base$sample_hpv16_dna_rna <- relevel(base$sample_hpv16_dna_rna, ref = "DNA-")
  base <- add_trt(base, exclude_palliative)
  paths <- setdiff(colnames(scores), "sample_geo_accession")
  res <- lapply(paths, function(p) {
    md <- base |> left_join(scores[, c("sample_geo_accession", p)] |> setNames(c("sample_geo_accession", "pathway_score")),
                            by = "sample_geo_accession") |> filter(!is.na(pathway_score))
    md$z <- as.numeric(scale(md$pathway_score))
    f <- if (with_trt) Surv(sample_pfs, sample_pfs_event) ~ z + sample_uicc_stage + sample_hpv16_dna_rna + trt else
      Surv(sample_pfs, sample_pfs_event) ~ z + sample_uicc_stage + sample_hpv16_dna_rna
    fit <- coxph(f, data = md, ties = "efron", x = TRUE, model = TRUE)
    s <- summary(fit); zph <- cox.zph(fit, transform = "km")$table
    data.frame(pathway = p, n = nrow(md), events = sum(md$sample_pfs_event == 1),
               beta = s$coefficients["z", "coef"], standard_error = s$coefficients["z", "se(coef)"],
               hazard_ratio = s$conf.int["z", "exp(coef)"], confidence_interval_low = s$conf.int["z", "lower .95"],
               confidence_interval_high = s$conf.int["z", "upper .95"], p_value = s$coefficients["z", "Pr(>|z|)"],
               pathway_PH_p = zph["z", "p"], global_PH_p = zph["GLOBAL", "p"])
  }) |> bind_rows() |> mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |> arrange(p_adj_BH, p_value)
  if (with_trt) write.csv(res, file.path(OUT, paste0("cox_plus_treatment_", tag, ".csv")), row.names = FALSE)
  res
}


# ---- internal control: without the treatment term these functions must reproduce the main results ----
ctl_lf <- run_limma_fgsea(FALSE, "control", with_trt = FALSE)
ctl_sc <- run_scores(FALSE, "control", with_trt = FALSE)
ctl_cx <- run_cox(FALSE, "control", with_trt = FALSE)
j <- merge(main_de, ctl_lf$de, by = "Symbol", suffixes = c("_m", "_c"))
stopifnot(nrow(j) == nrow(main_de), max(abs(j$t_m - j$t_c)) < 1e-8, max(abs(j$P.Value_m - j$P.Value_c)) < 1e-8)
jf <- merge(main_fg[, c("pathway", "NES", "padj")], ctl_lf$fg[, c("pathway", "NES", "padj")], by = "pathway", suffixes = c("_m", "_c"))
stopifnot(nrow(jf) == 50, max(abs(jf$NES_m - jf$NES_c)) < 1e-6, max(abs(jf$padj_m - jf$padj_c)) < 1e-6)
js <- merge(main_sc, ctl_sc, by = "pathway", suffixes = c("_m", "_c"))
stopifnot(nrow(js) == 50, max(abs(js$beta_early_vs_event_free_m - js$beta_early_vs_event_free_c)) < 1e-8, max(abs(js$p_value_m - js$p_value_c)) < 1e-8)
jc <- merge(main_cx, ctl_cx, by = "pathway", suffixes = c("_m", "_c"))
stopifnot(nrow(jc) == 50, max(abs(jc$beta_m - jc$beta_c)) < 1e-8, max(abs(jc$p_value_m - jc$p_value_c)) < 1e-8)
cat("Control passed: models without the treatment term reproduce the main results.\n")

# ---- run the primary coding and the palliative-excluded check -----------------------------
runs <- list(primary = FALSE, exclude_palliative = TRUE)
res <- lapply(names(runs), function(tag) {
  ex <- runs[[tag]]
  lf <- run_limma_fgsea(ex, tag); sc <- run_scores(ex, tag); cx <- run_cox(ex, tag)
  list(lf = lf, sc = sc, cx = cx)
})
names(res) <- names(runs)

# ---- comparison with the main analysis: outcomes P1 to P5 -----------------------------------
key4 <- c("E2F_TARGETS", "MYC_TARGETS_V1", "G2M_CHECKPOINT", "ESTROGEN_RESPONSE_LATE")
norm <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", sub("^HALLMARK_", "", x)))
summarise_run <- function(r, tag) {
  fg <- r$lf$fg; sig <- main_fg |> filter(padj < 0.05)
  m <- merge(sig[, c("pathway", "NES", "padj")], fg[, c("pathway", "NES", "padj")], by = "pathway",
             suffixes = c("_main", "_sens"))
  stopifnot(nrow(m) == 23)
  m$ratio <- m$NES_sens / m$NES_main
  newly <- fg |> filter(padj < 0.05, !pathway %in% sig$pathway)
  de <- merge(main_de[, c("Symbol", "t")], r$lf$de[, c("Symbol", "t", "adj.P.Val")], by = "Symbol", suffixes = c("_main", "_sens"))
  # score and Cox models
  cmp <- function(main, sens, est, lab) {
    j <- merge(main[, c("pathway", est, "p_value", "p_adj_BH")], sens[, c("pathway", est, "p_value", "p_adj_BH")],
               by = "pathway", suffixes = c("_main", "_sens"))
    nom <- j[j$p_value_main < 0.05, ]
    list(n_nominal_main = nrow(nom),
         direction_kept = sum(sign(nom[[paste0(est, "_main")]]) == sign(nom[[paste0(est, "_sens")]])),
         still_p05 = sum(nom$p_value_sens < 0.05),
         any_q05_sens = sum(j$p_adj_BH_sens < 0.05),
         key = j[norm(j$pathway) %in% key4, ] |>
           transmute(model = lab, pathway = pathway, estimate_main = .data[[paste0(est, "_main")]],
                     estimate_sens = .data[[paste0(est, "_sens")]],
                     pct_change = 100 * (estimate_sens - estimate_main) / abs(estimate_main),
                     p_main = p_value_main, p_sens = p_value_sens))
  }
  s <- cmp(main_sc, r$sc, "beta_early_vs_event_free", "sample-level score (beta)")
  cx <- cmp(main_cx, r$cx, "beta", "Cox (log hazard ratio)")
  n_sign <- sum(sign(m$NES_main) == sign(m$NES_sens)); n_q <- sum(m$padj_sens < 0.05, na.rm = TRUE)
  med <- median(m$ratio)
  verdict <- if ((23 - n_sign) > 2 || n_q < 12) "not robust"
             else if (n_sign == 23 && n_q >= 19 && med >= 0.8 && med <= 1.2) "robust" else "partly attenuated"
  list(summary = data.frame(
    coding = tag, n_limma_cohort = r$lf$n,
    P1_sets_keeping_NES_sign_of_23 = n_sign,
    P2_sets_still_q_below_0.05_of_23 = n_q, P2_other_sets_newly_q_below_0.05 = nrow(newly),
    P3_median_NES_ratio = round(med, 3), P3_min_NES_ratio = round(min(m$ratio), 3), P3_max_NES_ratio = round(max(m$ratio), 3),
    P4_spearman_t_main_vs_sens = round(cor(de$t_main, de$t_sens, method = "spearman"), 4),
    P4_genes_q_below_0.05_sens = sum(de$adj.P.Val < 0.05),
    P5_score_nominal_main = s$n_nominal_main, P5_score_direction_kept = s$direction_kept,
    P5_score_still_p_below_0.05 = s$still_p05, P5_score_any_q_below_0.05 = s$any_q05_sens,
    P5_cox_nominal_main = cx$n_nominal_main, P5_cox_direction_kept = cx$direction_kept,
    P5_cox_still_p_below_0.05 = cx$still_p05, P5_cox_any_q_below_0.05 = cx$any_q05_sens,
    prespecified_description = verdict),
    sets = m |> arrange(NES_main) |> mutate(coding = tag), key = bind_rows(s$key, cx$key) |> mutate(coding = tag))
}
sm <- lapply(names(res), function(tag) summarise_run(res[[tag]], tag))
write.csv(bind_rows(lapply(sm, `[[`, "summary")), file.path(OUT, "summary_P1_to_P5.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(sm, `[[`, "sets")), file.path(OUT, "comparison_23_sets.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(sm, `[[`, "key")), file.path(OUT, "comparison_key_programmes.csv"), row.names = FALSE)
print(t(bind_rows(lapply(sm, `[[`, "summary"))))
