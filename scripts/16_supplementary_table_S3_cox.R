# ============================================================
# 16_supplementary_table_S3_cox.R
#
# Creates:
#   Supplementary Table S3: all 50 adjusted Cox models and
#   proportional-hazards diagnostics.
#
# The script first tries likely output filenames. If your Cox file has
# a different name, set cox_input_path manually below.
# ============================================================

source("scripts/supplementary_table_helpers.R")

# Set this manually only when automatic discovery does not find the file.
cox_input_path <- NA_character_

candidate_paths <- c(
  "results_v2/GSE65858_v2_hallmark_cox_all50.csv",
  "results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.csv",
  "results_v2/GSE65858_v2_cox_hallmark_all50.csv",
  "results_v2/GSE65858_v2_hallmark_scores_cox_all50.csv",
  "results_v2/GSE65858_v2_hallmark_cox_all50.rds",
  "results_v2/GSE65858_v2_hallmark_cox_adjusted_all50.rds",
  "results_v2/GSE65858_v2_cox_hallmark_all50.rds"
)

pick_column <- function(data, exact = character(), pattern = NULL) {
  exact_match <- exact[exact %in% names(data)]
  if (length(exact_match) > 0) {
    return(exact_match[[1]])
  }

  if (!is.null(pattern)) {
    regex_match <- grep(
      pattern,
      names(data),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(regex_match) > 0) {
      return(regex_match[[1]])
    }
  }

  NA_character_
}

looks_like_cox_output <- function(path) {
  object <- tryCatch(read_result_file(path), error = function(e) NULL)
  if (is.null(object) || !is.data.frame(object)) return(FALSE)

  pathway_col <- pick_column(
    object,
    exact = c("pathway", "Pathway", "gs_name"),
    pattern = "^pathway$|gs_name"
  )

  hr_col <- pick_column(
    object,
    exact = c("hazard_ratio", "HR", "hr"),
    pattern = "hazard.*ratio|^hr$|exp.*coef"
  )

  !is.na(pathway_col) && !is.na(hr_col) && nrow(object) >= 50
}

if (is.na(cox_input_path)) {
  existing_candidates <- candidate_paths[file.exists(candidate_paths)]

  if (length(existing_candidates) > 0) {
    valid_candidates <- existing_candidates[
      vapply(existing_candidates, looks_like_cox_output, logical(1))
    ]
  } else {
    valid_candidates <- character()
  }

  if (length(valid_candidates) == 0) {
    discovered <- c(
      list.files(
        "results_v2",
        pattern = "\\.(csv|rds)$",
        full.names = TRUE,
        ignore.case = TRUE
      )
    )

    valid_candidates <- discovered[
      vapply(discovered, looks_like_cox_output, logical(1))
    ]
  }

  if (length(valid_candidates) == 0) {
    stop(
      "No all-50 Cox output could be identified automatically.\n",
      "Set cox_input_path near the top of this script to the saved ",
      "Cox CSV or RDS file, then rerun it."
    )
  }

  cox_input_path <- valid_candidates[[1]]
}

message("Using Cox input: ", cox_input_path)
cox_results <- read_result_file(cox_input_path)

pathway_col <- pick_column(
  cox_results,
  exact = c("pathway", "Pathway", "gs_name"),
  pattern = "^pathway$|gs_name"
)
n_col <- pick_column(
  cox_results,
  exact = c("n", "N", "sample_size"),
  pattern = "^n$|sample.*size"
)
events_col <- pick_column(
  cox_results,
  exact = c("events", "n_events", "event_count"),
  pattern = "event.*count|n.*event|^events$"
)
hr_col <- pick_column(
  cox_results,
  exact = c("hazard_ratio", "HR", "hr"),
  pattern = "hazard.*ratio|^hr$|exp.*coef"
)
low_col <- pick_column(
  cox_results,
  exact = c(
    "confidence_interval_low",
    "ci_low",
    "conf_low",
    "lower_95",
    "hr_ci_low"
  ),
  pattern = "ci.*low|conf.*low|lower.*95|hr.*low"
)
high_col <- pick_column(
  cox_results,
  exact = c(
    "confidence_interval_high",
    "ci_high",
    "conf_high",
    "upper_95",
    "hr_ci_high"
  ),
  pattern = "ci.*high|conf.*high|upper.*95|hr.*high"
)
p_col <- pick_column(
  cox_results,
  exact = c(
    "p_value",
    "cox_p_value",
    "association_p_value",
    "p"
  ),
  pattern = "^(p|p_value)$|cox.*p|association.*p"
)
q_col <- pick_column(
  cox_results,
  exact = c(
    "p_adj_BH",
    "q_value",
    "BH_q",
    "padj",
    "cox_q_value"
  ),
  pattern = "adj.*bh|q.*value|bh.*q|^padj$"
)
pathway_ph_col <- pick_column(
  cox_results,
  exact = c(
    "pathway_ph_p",
    "ph_pathway_p",
    "ph_p_value",
    "cox_zph_pathway_p",
    "ph_diagnostic_p"
  ),
  pattern = "pathway.*ph.*p|ph.*pathway.*p|zph.*pathway|ph.*diagnostic"
)
global_ph_col <- pick_column(
  cox_results,
  exact = c(
    "global_ph_p",
    "ph_global_p",
    "cox_zph_global_p",
    "global_p"
  ),
  pattern = "global.*ph.*p|ph.*global.*p|zph.*global"
)

required_mapped <- c(
  pathway = pathway_col,
  HR = hr_col,
  CI_low = low_col,
  CI_high = high_col,
  p_value = p_col,
  BH_q = q_col,
  pathway_PH_p = pathway_ph_col
)

if (anyNA(required_mapped)) {
  stop(
    "Could not map all required Cox columns.\n",
    paste(
      names(required_mapped),
      required_mapped,
      sep = " = ",
      collapse = "\n"
    ),
    "\n\nAvailable columns:\n",
    paste(names(cox_results), collapse = ", "),
    "\n\nEdit the candidate names in pick_column() if your output uses ",
    "different column names."
  )
}

value_or_na <- function(data, column_name) {
  if (is.na(column_name)) {
    rep(NA_real_, nrow(data))
  } else {
    data[[column_name]]
  }
}

table_s3_numeric <- data.frame(
  pathway = pathway_label(cox_results[[pathway_col]]),
  n = value_or_na(cox_results, n_col),
  events = value_or_na(cox_results, events_col),
  HR = cox_results[[hr_col]],
  ci_low = cox_results[[low_col]],
  ci_high = cox_results[[high_col]],
  raw_p = cox_results[[p_col]],
  BH_q = cox_results[[q_col]],
  pathway_PH_p = cox_results[[pathway_ph_col]],
  global_PH_p = value_or_na(cox_results, global_ph_col),
  stringsAsFactors = FALSE
) |>
  arrange(BH_q, raw_p)

assert_50_pathways(table_s3_numeric, "pathway")
stopifnot(nrow(table_s3_numeric) == 50L)

write.csv(
  table_s3_numeric,
  "supplementary_tables/Supplementary_Table_S3_all50_cox.csv",
  row.names = FALSE
)

table_s3_display <- table_s3_numeric |>
  transmute(
    Pathway = pathway,
    N = ifelse(is.na(n), "", as.character(as.integer(n))),
    Events = ifelse(is.na(events), "", as.character(as.integer(events))),
    `Adjusted HR (95% CI)` = format_ci(HR, ci_low, ci_high, digits = 2),
    `Raw p` = format_p(raw_p),
    `BH q` = format_p(BH_q),
    `Pathway PH p` = format_p(pathway_PH_p),
    `Global PH p` = format_p(global_PH_p)
  )

# Omit entirely unavailable optional columns.
empty_columns <- vapply(
  table_s3_display,
  function(x) all(is.na(x) | x == ""),
  logical(1)
)
table_s3_display <- table_s3_display[, !empty_columns, drop = FALSE]

ft_s3 <- make_flextable(table_s3_display, font_size = 11)
ft_s3 <- align(
  ft_s3,
  j = setdiff(names(table_s3_display), "Pathway"),
  align = "center",
  part = "all"
)

footnote_s3 <- paste(
  "All 50 MSigDB Hallmark pathway scores were examined in separate Cox proportional-hazards models adjusted for UICC stage and HPV16 DNA/RNA status.",
  "Scores were standardised within the complete-case survival cohort; HRs are therefore reported per one-standard-deviation increase in pathway score.",
  "Ninety-five per cent confidence intervals are model-based and were not adjusted for multiple testing.",
  "BH q-values were calculated across all 50 pathway-association tests.",
  "Pathway and global PH p-values are proportional-hazards diagnostics and were not treated as biological association tests.",
  "BH, Benjamini-Hochberg; CI, confidence interval; HPV, human papillomavirus;",
  "HR, hazard ratio; PFS, progression-free survival; PH, proportional hazards;",
  "UICC, Union for International Cancer Control."
)

write_table_docx(
  ft = ft_s3,
  title = paste(
    "Supplementary Table S3.",
    "Complete adjusted Cox models and proportional-hazards diagnostics."
  ),
  footnote = footnote_s3,
  output_path =
    "supplementary_tables/Supplementary_Table_S3_all50_cox.docx"
)
