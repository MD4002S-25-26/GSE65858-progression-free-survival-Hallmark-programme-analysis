# Reproduction test (9 November 2026)

A clean folder containing only `scripts/` was used. The pipeline was run with
`Rscript -e 'source("scripts/00_run_complete_pipeline.R")'` on macOS, R 4.6.0
(see `sessionInfo_run_2026-11-09.txt`). Total run time was about 65 seconds, plus the
download of the GEO files.

| Check | Result |
|---|---|
| Fresh GEO series matrix and GPL10558 annotation vs the files used for the dissertation (MD5) | identical |
| `pheno3_matched.rds` rebuilt from GEO vs the dissertation object | `identical()` = TRUE |
| `expr2_probe_matrix.rds` rebuilt from GEO vs the dissertation object | `identical()` = TRUE |
| 16 result and table CSV files vs the dissertation outputs (tolerance 1e-8) | all identical |
| 9 leading-edge CSV files vs the earlier manuscript analysis | byte-identical |
| `19_verify_headline_results.R` (headline numbers in the manuscript) | all checks passed |

Defects found and fixed while testing (scripts 01-17 otherwise unchanged):
1. Nothing created the `results_v2/` folder before script 02 wrote to it. The master script now creates all output folders.
2. Script 09 called `install.packages()` unconditionally. Replaced with an availability check.
3. The phenotype object was originally built by hand in Excel from the GEO metadata. Script `00a_prepare_inputs_from_GEO.R` now builds it directly from the GEO file and reproduces the original object exactly.
