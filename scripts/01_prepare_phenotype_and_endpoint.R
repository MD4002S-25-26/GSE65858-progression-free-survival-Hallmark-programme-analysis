# ============================================================
# 01_prepare_phenotype_and_endpoint.R
# Purpose:
#   Clean the matched GSE65858 phenotype data and define the
#   corrected 2-year PFS landmark groups.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Load matched phenotype
# ------------------------------------------------------------
pheno3 <- readRDS("pheno3_matched.rds")

# ------------------------------------------------------------
# 2. Helper function for binary event variables
# ------------------------------------------------------------
to_binary <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  
  x <- as.character(x)
  
  dplyr::case_when(
    x %in% c("TRUE", "True", "1") ~ 1,
    x %in% c("FALSE", "False", "0") ~ 0,
    TRUE ~ suppressWarnings(as.numeric(x))
  )
}

# ------------------------------------------------------------
# 3. Clean outcome and covariate variables
# ------------------------------------------------------------
pheno_analysis <- pheno3 |>
  mutate(
    sample_pfs = as.numeric(sample_pfs),
    sample_pfs_event = to_binary(sample_pfs_event),
    os = as.numeric(os),
    os_event = to_binary(os_event),
    sample_age = as.numeric(sample_age),
    
    sample_uicc_stage = factor(sample_uicc_stage),
    sample_hpv16_dna_rna = factor(sample_hpv16_dna_rna)
  )

# ------------------------------------------------------------
# 4. Define corrected 2-year PFS landmark grouping
#
# PFS event within 2 years:
#   Event occurred on or before day 730.
#
# PFS-event-free at 2 years:
#   Patient was still under observation beyond day 730,
#   whether subsequently experiencing an event or being censored.
#
# Censored before 2 years:
#   No event occurred, but observation ended on or before day 730.
# ------------------------------------------------------------
landmark_days <- 730

pheno_analysis <- pheno_analysis |>
  mutate(
    pfs_2yr_group = case_when(
      sample_pfs_event == 1 &
        sample_pfs <= landmark_days ~
        "PFS event within 2 years",
      
      sample_pfs > landmark_days ~
        "PFS-event-free at 2 years",
      
      sample_pfs_event == 0 &
        sample_pfs <= landmark_days ~
        "Censored before 2 years",
      
      TRUE ~ NA_character_
    ),
    
    pfs_2yr_group = factor(
      pfs_2yr_group,
      levels = c(
        "PFS-event-free at 2 years",
        "PFS event within 2 years",
        "Censored before 2 years"
      )
    )
  )

# ------------------------------------------------------------
# 5. Audit checks
# ------------------------------------------------------------

# Overall group counts
print(table(pheno_analysis$pfs_2yr_group, useNA = "ifany"))

# Group by event status
print(
  with(
    pheno_analysis,
    table(pfs_2yr_group, sample_pfs_event, useNA = "ifany")
  )
)

# Number of events occurring after the landmark
late_events <- sum(
  pheno_analysis$sample_pfs_event == 1 &
    pheno_analysis$sample_pfs > landmark_days,
  na.rm = TRUE
)

cat("Events occurring after 2 years:", late_events, "\n")

# Number genuinely censored before or at the landmark
early_censored <- sum(
  pheno_analysis$sample_pfs_event == 0 &
    pheno_analysis$sample_pfs <= landmark_days,
  na.rm = TRUE
)

cat("Censored before or at 2 years:", early_censored, "\n")

# Logical checks
stopifnot(
  all(
    pheno_analysis$sample_pfs[
      pheno_analysis$pfs_2yr_group == "PFS event within 2 years"
    ] <= landmark_days,
    na.rm = TRUE
  )
)

stopifnot(
  all(
    pheno_analysis$sample_pfs[
      pheno_analysis$pfs_2yr_group == "PFS-event-free at 2 years"
    ] > landmark_days,
    na.rm = TRUE
  )
)

stopifnot(
  all(
    pheno_analysis$sample_pfs_event[
      pheno_analysis$pfs_2yr_group == "Censored before 2 years"
    ] == 0,
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 6. Save definitive corrected phenotype
# ------------------------------------------------------------
saveRDS(
  pheno_analysis,
  "pheno_analysis_corrected_landmark.rds"
)

write.csv(
  pheno_analysis,
  "pheno_analysis_corrected_landmark.csv",
  row.names = FALSE
)
table(pheno_analysis$pfs_2yr_group, useNA = "ifany")
with(
  pheno_analysis,
  table(pfs_2yr_group, sample_pfs_event, useNA = "ifany")
)
late_events
early_censored
file.exists("pheno_analysis_corrected_landmark.rds")
file.exists("pheno_analysis_corrected_landmark.csv")