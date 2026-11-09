# ============================================================
# 20_median_followup.R
#
# Median follow-up for STROBE item 14(c), by the reverse Kaplan-Meier method
# (censoring and event indicators swapped; Schemper and Smith 1996).
# Run from the repository root after scripts 01 and 06 have finished.
#
# Assumption: the deposited PFS times are in days. GEO does not label the unit; days is
# consistent with the 730-day landmark used in this analysis. Months = days / 30.4375.
# The 241-patient fixed-time cohort is deliberately not summarised here: it excludes
# patients censored before 730 days, so reverse Kaplan-Meier would be biased.
# ============================================================
suppressPackageStartupMessages({ library(survival); library(dplyr) })

ph <- readRDS("pheno_analysis_corrected_landmark.rds")
sc <- readRDS("results_v2/GSE65858_v2_hallmark_pathway_scores.rds")
ph$t <- as.numeric(ph$sample_pfs)
ph$e <- as.numeric(ph$sample_pfs_event)

# Cox cohort: the same inclusion rules as script 07
cox <- ph |>
  filter(sample_geo_accession %in% sc$sample_geo_accession,
         !is.na(t), !is.na(e), !is.na(sample_uicc_stage),
         !is.na(sample_hpv16_dna_rna), as.character(sample_hpv16_dna_rna) != "NA")
all_pts <- ph |> filter(!is.na(t), !is.na(e))

summarise_followup <- function(d, label) {
  q <- quantile(survfit(Surv(t, 1 - e) ~ 1, data = d), probs = 0.5, conf.int = TRUE)
  data.frame(
    cohort = label, n = nrow(d), pfs_events = sum(d$e),
    median_followup_days = unname(q$quantile),
    ci_lower_days = unname(q$lower), ci_upper_days = unname(q$upper),
    median_followup_months = round(unname(q$quantile) / 30.4375, 1),
    ci_lower_months = round(unname(q$lower) / 30.4375, 1),
    ci_upper_months = round(unname(q$upper) / 30.4375, 1),
    person_years = round(sum(d$t) / 365.25, 0), max_followup_days = max(d$t)
  )
}
out <- bind_rows(summarise_followup(all_pts, "All patients with PFS data"),
                 summarise_followup(cox, "Cox analysis cohort"))
write.csv(out, "results_v2/median_followup.csv", row.names = FALSE)
print(out)
