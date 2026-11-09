# ============================================================
# 00_run_complete_pipeline.R
#
# Runs the whole analysis from the public GEO files to the manuscript outputs.
# Run from the repository root:
#   Rscript -e 'source("scripts/00_run_complete_pipeline.R")'
# See README.md for requirements and expected outputs.
# ============================================================

# Output folders expected by the analysis scripts
for (d in c("results_v2", "figures_v2", "tables_v2", "supplementary_tables")) {
  dir.create(d, showWarnings = FALSE)
}

# Step 0: build pheno3_matched.rds and expr2_probe_matrix.rds from GEO,
# and fetch the GPL10558 annotation (checksum-verified)
source("scripts/00a_prepare_inputs_from_GEO.R")

# Steps 1-17: the dissertation analysis (unchanged scripts)
source("scripts/01_prepare_phenotype_and_endpoint.R")
source("scripts/02_probe_to_gene_mapping.R")
source("scripts/03_limma_unadjusted.R")
source("scripts/04_limma_adjusted.R")
source("scripts/05_fgsea_hallmark.R")
source("scripts/06_hallmark_pathway_scores.R")
source("scripts/07_hallmark_cox_models.R")
source("scripts/08_pipeline_audit.R")
source("scripts/09_table1_cohort_characteristics.R")
source("scripts/10_table2_integrated_pathway_summary.R")
source("scripts/11_figure3_adjusted_volcano.R")
source("scripts/12_figure4_hallmark_enrichment.R")
source("scripts/13_figure5_sample_level_pathway_scores.R")
source("scripts/14_supplementary_table_S1_fgsea.R")
source("scripts/15_supplementary_table_S2_pathway_scores.R")
source("scripts/16_supplementary_table_S3_cox.R")
source("scripts/17_run_all_supplementary_tables.R")

# Step 18: leading-edge overlap analysis added for the manuscript
source("scripts/18_leading_edge_overlap.R")

# Step 20: median follow-up (STROBE 14c); run before the final verification step
source("scripts/20_median_followup.R")

# Step 21: post hoc treatment-adjusted sensitivity analysis (plan: verification/treatment_sensitivity_plan_2026-11-09.md)
source("scripts/21_treatment_sensitivity.R")

# Step 19: confirm the headline numbers reported in the manuscript
source("scripts/19_verify_headline_results.R")

# Record the software environment of this run
writeLines(capture.output(sessionInfo()), "results_v2/sessionInfo_run.txt")
