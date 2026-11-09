# ============================================================
# supplementary_table_helpers.R
# Shared helpers for MD4002 supplementary Hallmark tables
# ============================================================

required_packages <- c("dplyr", "flextable", "officer")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(flextable)
library(officer)

dir.create("supplementary_tables", showWarnings = FALSE, recursive = TRUE)

ascii_text <- function(x) {
  x <- as.character(x)
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x, perl = TRUE)
  x <- gsub("\u03B1", "alpha", x, fixed = TRUE)
  x <- gsub("\u03BA", "kappa", x, fixed = TRUE)
  x
}

pathway_label <- function(x) {
  key <- gsub("^HALLMARK_", "", as.character(x))

  special <- c(
    "E2F_TARGETS" = "E2F Targets",
    "MYC_TARGETS_V1" = "MYC Targets V1",
    "MYC_TARGETS_V2" = "MYC Targets V2",
    "G2M_CHECKPOINT" = "G2M Checkpoint",
    "TNFA_SIGNALING_VIA_NFKB" = "TNF-alpha Signaling via NF-kappaB",
    "IL6_JAK_STAT3_SIGNALING" = "IL6-JAK-STAT3 Signaling",
    "IL2_STAT5_SIGNALING" = "IL2-STAT5 Signaling",
    "KRAS_SIGNALING_UP" = "KRAS Signaling Up",
    "KRAS_SIGNALING_DN" = "KRAS Signaling Down",
    "MTORC1_SIGNALING" = "mTORC1 Signaling",
    "TGF_BETA_SIGNALING" = "TGF-beta Signaling",
    "PI3K_AKT_MTOR_SIGNALING" = "PI3K-AKT-mTOR Signaling",
    "P53_PATHWAY" = "p53 Pathway",
    "UV_RESPONSE_UP" = "UV Response Up",
    "UV_RESPONSE_DN" = "UV Response Down",
    "INTERFERON_ALPHA_RESPONSE" = "Interferon Alpha Response",
    "INTERFERON_GAMMA_RESPONSE" = "Interferon Gamma Response",
    "OXIDATIVE_PHOSPHORYLATION" = "Oxidative Phosphorylation",
    "UNFOLDED_PROTEIN_RESPONSE" = "Unfolded Protein Response",
    "ESTROGEN_RESPONSE_EARLY" = "Estrogen Response Early",
    "ESTROGEN_RESPONSE_LATE" = "Estrogen Response Late",
    "ALLOGRAFT_REJECTION" = "Allograft Rejection",
    "INFLAMMATORY_RESPONSE" = "Inflammatory Response",
    "HEME_METABOLISM" = "Heme Metabolism",
    "MITOTIC_SPINDLE" = "Mitotic Spindle",
    "DNA_REPAIR" = "DNA Repair",
    "PROTEIN_SECRETION" = "Protein Secretion",
    "APICAL_JUNCTION" = "Apical Junction",
    "APICAL_SURFACE" = "Apical Surface",
    "EPITHELIAL_MESENCHYMAL_TRANSITION" =
      "Epithelial-Mesenchymal Transition"
  )

  out <- unname(special[key])
  missing <- is.na(out)

  fallback <- gsub("_", " ", key[missing])
  fallback <- tools::toTitleCase(tolower(fallback))
  out[missing] <- fallback

  ascii_text(out)
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

format_ci <- function(estimate, low, high, digits = 3) {
  ifelse(
    is.na(estimate) | is.na(low) | is.na(high),
    "",
    paste0(
      sprintf(paste0("%.", digits, "f"), estimate),
      " (",
      sprintf(paste0("%.", digits, "f"), low),
      " to ",
      sprintf(paste0("%.", digits, "f"), high),
      ")"
    )
  )
}

read_result_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    object <- readRDS(path)
    if (!is.data.frame(object)) {
      object <- as.data.frame(object)
    }
    return(object)
  }

  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

assert_columns <- function(data, columns, object_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      "\nAvailable columns: ",
      paste(names(data), collapse = ", ")
    )
  }
}

assert_50_pathways <- function(data, pathway_column = "pathway") {
  n_pathways <- dplyr::n_distinct(data[[pathway_column]])
  if (n_pathways != 50L) {
    stop(
      "Expected exactly 50 Hallmark pathways, but found ",
      n_pathways,
      ". Check that the corrected all-50 output was loaded."
    )
  }
}

make_flextable <- function(data, font_size = 11, landscape = TRUE) {
  ft <- flextable(data)
  ft <- theme_booktabs(ft)
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "left", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")
  ft <- font(ft, fontname = "Arial", part = "all")
  ft <- padding(ft, padding = 2, part = "all")
  ft <- set_table_properties(
    ft,
    layout = "autofit",
    opts_word = list(split = TRUE, repeat_headers = TRUE)
  )
  ft <- autofit(ft)
  ft
}

write_table_docx <- function(
    ft,
    title,
    footnote,
    output_path,
    landscape = TRUE
) {
  doc <- read_docx()
  doc <- body_add_par(
    doc,
    value = title,
    style = "heading 1"
  )
  doc <- body_add_flextable(doc, value = ft)
  doc <- body_add_par(
    doc,
    value = footnote,
    style = "Normal"
  )

  if (landscape) {
    section <- prop_section(
      page_size = page_size(orient = "landscape"),
      page_margins = page_mar(
        top = 0.6,
        bottom = 0.6,
        left = 0.6,
        right = 0.6
      ),
      type = "continuous"
    )
    doc <- body_set_default_section(doc, section)
  }

  print(doc, target = output_path)
  message("Created: ", output_path)
}
