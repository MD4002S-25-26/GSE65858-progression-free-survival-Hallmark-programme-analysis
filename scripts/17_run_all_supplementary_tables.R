# ============================================================
# 17_run_all_supplementary_tables.R
#
# Run from the project root after the corrected analytical pipeline.
# Expected project structure:
#   scripts/
#   results_v2/
#   supplementary_tables/
# ============================================================

scripts_to_run <- c(
  "scripts/14_supplementary_table_S1_fgsea.R",
  "scripts/15_supplementary_table_S2_pathway_scores.R",
  "scripts/16_supplementary_table_S3_cox.R"
)

for (script in scripts_to_run) {
  message("\n--- Running ", script, " ---")
  source(script, echo = FALSE)
}

message(
  "\nAll supplementary-table scripts completed. ",
  "Check the supplementary_tables folder."
)
