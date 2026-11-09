pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

names(pheno_analysis)
data.frame(
  variable = names(pheno_analysis),
  missing_n = vapply(
    pheno_analysis,
    function(x) sum(is.na(x)),
    integer(1)
  ),
  unique_values = vapply(
    pheno_analysis,
    function(x) length(unique(x)),
    integer(1)
  )
)
# (The original script called install.packages() here. Replaced by an availability
# check so that the pipeline never modifies the user's R library.)
stopifnot(
  requireNamespace("gtsummary", quietly = TRUE),
  requireNamespace("flextable", quietly = TRUE),
  requireNamespace("officer", quietly = TRUE)
)
dir.create(
  "tables_v2",
  showWarnings = FALSE
)
# ============================================================
# 09_table1_cohort_characteristics.R
#
# Purpose:
# Create Table 1 describing the clinical characteristics of
# the corrected two-year PFS landmark analysis cohort.
#
# The landmark cohort includes:
#   - PFS-event-free at 2 years
#   - PFS event within 2 years
#
# Patients censored before 2 years are not included in this
# table because they were excluded from landmark analyses.
# ============================================================

library(dplyr)
library(gtsummary)
library(flextable)
library(officer)

# ------------------------------------------------------------
# 1. Create output folder
# ------------------------------------------------------------
dir.create(
  "tables_v2",
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load corrected phenotype data
# ------------------------------------------------------------
pheno_analysis <- readRDS(
  "pheno_analysis_corrected_landmark.rds"
)

# ------------------------------------------------------------
# 3. Helper function for coded missing values
#
# The phenotype data contain some missing information stored as
# text such as "NA" rather than as a genuine R NA value.
# ------------------------------------------------------------
clean_character_missing <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "na",
      "Na",
      "Unknown",
      "UNKNOWN",
      "unknown"
    )
  ] <- NA_character_
  
  x
}

# ------------------------------------------------------------
# 4. Create the two-year landmark Table 1 dataset
# ------------------------------------------------------------
table1_data <- pheno_analysis |>
  filter(
    pfs_2yr_group %in% c(
      "PFS-event-free at 2 years",
      "PFS event within 2 years"
    )
  ) |>
  mutate(
    landmark_group = case_when(
      pfs_2yr_group ==
        "PFS-event-free at 2 years" ~
        "PFS-event-free at 2 years",
      
      pfs_2yr_group ==
        "PFS event within 2 years" ~
        "PFS event within 2 years",
      
      TRUE ~ NA_character_
    ),
    
    landmark_group = factor(
      landmark_group,
      levels = c(
        "PFS-event-free at 2 years",
        "PFS event within 2 years"
      )
    ),
    
    # Continuous variables
    sample_age = as.numeric(
      sample_age
    ),
    
    sample_packyears = suppressWarnings(
      as.numeric(
        clean_character_missing(
          sample_packyears
        )
      )
    ),
    
    # Categorical variables
    sample_gender = factor(
      clean_character_missing(
        sample_gender
      )
    ),
    
    sample_smoking = factor(
      clean_character_missing(
        sample_smoking
      )
    ),
    
    sample_alcohol = factor(
      clean_character_missing(
        sample_alcohol
      )
    ),
    
    sample_tumour_type = factor(
      clean_character_missing(
        sample_tumour_type
      )
    ),
    
    sample_tumour_site = factor(
      clean_character_missing(
        sample_tumour_site
      )
    ),
    
    sample_uicc_stage = factor(
      clean_character_missing(
        sample_uicc_stage
      ),
      levels = c(
        "I",
        "II",
        "III",
        "IVA",
        "IVB",
        "IVC"
      )
    ),
    
    sample_t_category = factor(
      clean_character_missing(
        sample_t_category
      )
    ),
    
    sample_n_category = factor(
      clean_character_missing(
        sample_n_category
      )
    ),
    
    sample_distant_metastasis = factor(
      clean_character_missing(
        sample_distant_metastasis
      )
    ),
    
    sample_treatment = factor(
      clean_character_missing(
        sample_treatment
      )
    ),
    
    sample_hpv16_dna_rna = factor(
      clean_character_missing(
        sample_hpv16_dna_rna
      ),
      levels = c(
        "DNA-",
        "DNA+RNA-",
        "DNA+RNA+"
      )
    )
  ) |>
  select(
    landmark_group,
    sample_age,
    sample_gender,
    sample_smoking,
    sample_packyears,
    sample_alcohol,
    sample_tumour_type,
    sample_tumour_site,
    sample_uicc_stage,
    sample_t_category,
    sample_n_category,
    sample_distant_metastasis,
    sample_treatment,
    sample_hpv16_dna_rna
  )

# ------------------------------------------------------------
# 5. Audit the Table 1 dataset
# ------------------------------------------------------------
cat(
  "\nTABLE 1 GROUP COUNTS\n"
)

print(
  table(
    table1_data$landmark_group,
    useNA = "ifany"
  )
)

cat(
  "\nTotal Table 1 sample size:",
  nrow(table1_data),
  "\n"
)

stopifnot(
  nrow(table1_data) == 241
)

stopifnot(
  sum(
    table1_data$landmark_group ==
      "PFS-event-free at 2 years"
  ) == 135
)

stopifnot(
  sum(
    table1_data$landmark_group ==
      "PFS event within 2 years"
  ) == 106
)

# Display levels so they can be checked before final reporting
categorical_variables <- c(
  "sample_gender",
  "sample_smoking",
  "sample_alcohol",
  "sample_tumour_type",
  "sample_tumour_site",
  "sample_uicc_stage",
  "sample_t_category",
  "sample_n_category",
  "sample_distant_metastasis",
  "sample_treatment",
  "sample_hpv16_dna_rna"
)

cat(
  "\nCATEGORICAL VARIABLE LEVELS\n"
)

for (variable_name in categorical_variables) {
  
  cat(
    "\n",
    variable_name,
    "\n",
    sep = ""
  )
  
  print(
    table(
      table1_data[[variable_name]],
      useNA = "ifany"
    )
  )
}

# ------------------------------------------------------------
# 6. Save the underlying Table 1 dataset
# ------------------------------------------------------------
write.csv(
  table1_data,
  "tables_v2/Table1_landmark_analysis_dataset.csv",
  row.names = FALSE
)

saveRDS(
  table1_data,
  "tables_v2/Table1_landmark_analysis_dataset.rds"
)

# ------------------------------------------------------------
# 7. Generate Table 1
#
# Continuous variables:
#   median (Q1, Q3)
#
# Categorical variables:
#   n (%)
#
# No baseline-comparison p-values are included.
# ------------------------------------------------------------
table1 <- table1_data |>
  tbl_summary(
    by = landmark_group,
    
    include = c(
      sample_age,
      sample_gender,
      sample_smoking,
      sample_packyears,
      sample_alcohol,
      sample_tumour_type,
      sample_tumour_site,
      sample_uicc_stage,
      sample_t_category,
      sample_n_category,
      sample_distant_metastasis,
      sample_treatment,
      sample_hpv16_dna_rna
    ),
    
    label = list(
      sample_age ~
        "Age, years",
      
      sample_gender ~
        "Sex",
      
      sample_smoking ~
        "Smoking status",
      
      sample_packyears ~
        "Smoking exposure, pack-years",
      
      sample_alcohol ~
        "Alcohol consumption",
      
      sample_tumour_type ~
        "Tumour type",
      
      sample_tumour_site ~
        "Primary tumour site",
      
      sample_uicc_stage ~
        "UICC stage",
      
      sample_t_category ~
        "T category",
      
      sample_n_category ~
        "N category",
      
      sample_distant_metastasis ~
        "Distant metastasis",
      
      sample_treatment ~
        "Treatment",
      
      sample_hpv16_dna_rna ~
        "HPV16 DNA/RNA group"
    ),
    
    type = list(
      sample_age ~ "continuous",
      sample_packyears ~ "continuous",
      all_categorical() ~ "categorical"
    ),
    
    statistic = list(
      all_continuous() ~
        "{median} ({p25}, {p75})",
      
      all_categorical() ~
        "{n} ({p}%)"
    ),
    
    digits = list(
      sample_age ~ 1,
      sample_packyears ~ 1
    ),
    
    missing = "ifany",
    missing_text = "Missing/unknown",
    
    percent = "column"
  ) |>
  add_overall(
    last = FALSE
  ) |>
  modify_header(
    label ~ "**Characteristic**"
  ) |>
  modify_caption(
    paste0(
      "**Table 1. Clinical characteristics of the ",
      "two-year PFS landmark analysis cohort.**  \n",
      "Values are median (Q1, Q3) or n (%). ",
      "Twenty-nine patients censored before two years ",
      "were excluded from the landmark comparison."
    )
  ) |>
  bold_labels()

# Display the table in RStudio
table1

# ------------------------------------------------------------
# 8. Save the table object
# ------------------------------------------------------------
saveRDS(
  table1,
  "tables_v2/Table1_cohort_characteristics_gtsummary.rds"
)

# ------------------------------------------------------------
# 9. Convert to a Word-compatible flextable
# ------------------------------------------------------------
table1_flex <- as_flex_table(
  table1
)

table1_flex <- table1_flex |>
  autofit() |>
  fontsize(
    size = 9,
    part = "all"
  ) |>
  font(
    fontname = "Arial",
    part = "all"
  ) |>
  valign(
    valign = "top",
    part = "all"
  )

# ------------------------------------------------------------
# 10. Save as a Word document
# ------------------------------------------------------------
save_as_docx(
  "Table 1" = table1_flex,
  path =
    "tables_v2/Table1_cohort_characteristics.docx"
)

cat(
  "\nTable 1 completed successfully.\n"
)

cat(
  "Saved Word table:",
  "tables_v2/Table1_cohort_characteristics.docx\n"
)

cat(
  "Saved analysis dataset:",
  "tables_v2/Table1_landmark_analysis_dataset.csv\n"
)
file.exists(
  "tables_v2/Table1_cohort_characteristics.docx"
)

file.exists(
  "tables_v2/Table1_landmark_analysis_dataset.csv"
)
cat(
  "\nTable 1 completed successfully.\n"
)

cat(
  "Saved Word table:",
  "tables_v2/Table1_cohort_characteristics.docx\n"
)

cat(
  "Saved analysis dataset:",
  "tables_v2/Table1_landmark_analysis_dataset.csv\n"
)

stopifnot(
  file.exists("tables_v2/Table1_cohort_characteristics.docx"),
  file.exists("tables_v2/Table1_landmark_analysis_dataset.csv")
)