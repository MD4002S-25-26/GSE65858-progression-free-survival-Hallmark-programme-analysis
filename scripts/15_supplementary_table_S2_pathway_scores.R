# ============================================================
# 15_supplementary_table_S2_pathway_scores.R
#
# Creates:
#   Supplementary Table S2: all 50 sample-level pathway-score analyses
#   - unadjusted Wilcoxon comparison
#   - adjusted linear-regression coefficient and 95% CI
#   - editable Word document
#   - clean CSV
# ============================================================

source("scripts/supplementary_table_helpers.R")

wilcoxon_path <-
  "results_v2/GSE65858_v2_hallmark_scores_wilcoxon_all50.csv"

adjusted_path <-
  "results_v2/GSE65858_v2_hallmark_scores_adjusted_all50.csv"

wilcoxon_results <- read_result_file(wilcoxon_path)
adjusted_results <- read_result_file(adjusted_path)

assert_columns(
  wilcoxon_results,
  c(
    "pathway",
    "n_early",
    "n_event_free",
    "median_early",
    "median_event_free",
    "median_difference",
    "p_value",
    "p_adj_BH"
  ),
  "Wilcoxon pathway-score output"
)

assert_columns(
  adjusted_results,
  c(
    "pathway",
    "n",
    "beta_early_vs_event_free",
    "confidence_interval_low",
    "confidence_interval_high",
    "p_value",
    "p_adj_BH"
  ),
  "Adjusted pathway-score output"
)

assert_50_pathways(wilcoxon_results)
assert_50_pathways(adjusted_results)

wilcoxon_clean <- wilcoxon_results |>
  transmute(
    pathway,
    n_event_free,
    n_early,
    median_event_free,
    median_early,
    median_difference_early_minus_event_free = median_difference,
    wilcoxon_raw_p = p_value,
    wilcoxon_BH_q = p_adj_BH
  )

adjusted_clean <- adjusted_results |>
  transmute(
    pathway,
    adjusted_n = n,
    adjusted_beta_early_vs_event_free = beta_early_vs_event_free,
    adjusted_ci_low = confidence_interval_low,
    adjusted_ci_high = confidence_interval_high,
    adjusted_raw_p = p_value,
    adjusted_BH_q = p_adj_BH
  )

table_s2_numeric <- full_join(
  wilcoxon_clean,
  adjusted_clean,
  by = "pathway"
) |>
  mutate(pathway_label = pathway_label(pathway)) |>
  select(pathway_label, everything(), -pathway) |>
  arrange(adjusted_BH_q, adjusted_raw_p)

stopifnot(nrow(table_s2_numeric) == 50L)
stopifnot(!anyDuplicated(table_s2_numeric$pathway_label))

write.csv(
  table_s2_numeric,
  "supplementary_tables/Supplementary_Table_S2_all50_pathway_scores.csv",
  row.names = FALSE
)

table_s2a_display <- table_s2_numeric |>
  transmute(
    Pathway = pathway_label,
    `Event-free median` = format_num(median_event_free, 3),
    `Early-event median` = format_num(median_early, 3),
    `Median difference` =
      format_num(median_difference_early_minus_event_free, 3),
    `Raw p` = format_p(wilcoxon_raw_p),
    `BH q` = format_p(wilcoxon_BH_q)
  )

table_s2b_display <- table_s2_numeric |>
  transmute(
    Pathway = pathway_label,
    N = as.character(as.integer(adjusted_n)),
    `Adjusted beta (95% CI)` = format_ci(
      adjusted_beta_early_vs_event_free,
      adjusted_ci_low,
      adjusted_ci_high,
      digits = 3
    ),
    `Raw p` = format_p(adjusted_raw_p),
    `BH q` = format_p(adjusted_BH_q)
  )

ft_s2a <- make_flextable(table_s2a_display, font_size = 11)
ft_s2a <- align(
  ft_s2a,
  j = setdiff(names(table_s2a_display), "Pathway"),
  align = "center",
  part = "all"
)

ft_s2b <- make_flextable(table_s2b_display, font_size = 11)
ft_s2b <- align(
  ft_s2b,
  j = setdiff(names(table_s2b_display), "Pathway"),
  align = "center",
  part = "all"
)

footnote_s2a <- paste(
  "All 50 MSigDB Hallmark gene sets are shown.",
  "Median difference is calculated as the PFS-event-within-two-years group minus the PFS-event-free group;",
  "negative values therefore indicate higher pathway scores in the PFS-event-free group.",
  "Unadjusted comparisons used Wilcoxon rank-sum tests in 241 patients.",
  "BH correction was applied across all 50 pathways.",
  "BH, Benjamini-Hochberg; PFS, progression-free survival."
)

footnote_s2b <- paste(
  "Adjusted beta is calculated as the PFS-event-within-two-years group minus the PFS-event-free group;",
  "negative values therefore indicate higher pathway scores in the PFS-event-free group.",
  "Linear-regression models included UICC stage and HPV16 DNA/RNA status and used 222 complete cases.",
  "Ninety-five per cent confidence intervals are model-based and were not adjusted for multiple testing.",
  "BH correction was applied across all 50 adjusted pathway-score models.",
  "BH, Benjamini-Hochberg; CI, confidence interval; HPV, human papillomavirus;",
  "PFS, progression-free survival; UICC, Union for International Cancer Control."
)

doc <- read_docx()
doc <- body_add_par(
  doc,
  "Supplementary Table S2A. Complete unadjusted sample-level Hallmark pathway-score comparisons.",
  style = "heading 1"
)
doc <- body_add_flextable(doc, ft_s2a)
doc <- body_add_par(doc, footnote_s2a, style = "Normal")
doc <- body_add_break(doc)
doc <- body_add_par(
  doc,
  "Supplementary Table S2B. Complete adjusted sample-level Hallmark pathway-score models.",
  style = "heading 1"
)
doc <- body_add_flextable(doc, ft_s2b)
doc <- body_add_par(doc, footnote_s2b, style = "Normal")
doc <- body_set_default_section(
  doc,
  prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(
      top = 0.6,
      bottom = 0.6,
      left = 0.6,
      right = 0.6
    ),
    type = "continuous"
  )
)

print(
  doc,
  target =
    "supplementary_tables/Supplementary_Table_S2_all50_pathway_scores.docx"
)
message(
  "Created: supplementary_tables/Supplementary_Table_S2_all50_pathway_scores.docx"
)
