# ============================================================
# 14_supplementary_table_S1_fgsea.R
#
# Creates:
#   Supplementary Table S1: all 50 ranked Hallmark enrichment results
#   - editable Word document
#   - clean CSV
# ============================================================

source("scripts/supplementary_table_helpers.R")

input_path <- "results_v2/GSE65858_v2_fgsea_hallmark.csv"

fgsea_results <- read_result_file(input_path)

assert_columns(
  fgsea_results,
  c("pathway", "size", "NES", "pval", "padj"),
  "Ranked Hallmark enrichment output"
)
assert_50_pathways(fgsea_results)

if (!"enrichment_direction" %in% names(fgsea_results)) {
  fgsea_results <- fgsea_results |>
    mutate(
      enrichment_direction = case_when(
        NES > 0 ~ "PFS event within two years",
        NES < 0 ~ "PFS-event-free at two years",
        TRUE ~ "No direction"
      )
    )
}

table_s1_numeric <- fgsea_results |>
  transmute(
    pathway = pathway_label(pathway),
    genes_in_analysis = as.integer(size),
    NES = NES,
    raw_p = pval,
    BH_q = padj,
    direction = ascii_text(enrichment_direction)
  ) |>
  arrange(BH_q, desc(abs(NES)))

stopifnot(nrow(table_s1_numeric) == 50L)

write.csv(
  table_s1_numeric,
  "supplementary_tables/Supplementary_Table_S1_all50_fgsea.csv",
  row.names = FALSE
)

table_s1_display <- table_s1_numeric |>
  transmute(
    Pathway = pathway,
    `Genes in analysis` = genes_in_analysis,
    NES = format_num(NES, 2),
    `Raw p` = format_p(raw_p),
    `BH q` = format_p(BH_q),
    Direction = direction
  )

ft_s1 <- make_flextable(table_s1_display, font_size = 11)
ft_s1 <- align(
  ft_s1,
  j = c("Genes in analysis", "NES", "Raw p", "BH q"),
  align = "center",
  part = "all"
)

footnote_s1 <- paste(
  "All 50 MSigDB Hallmark gene sets are shown.",
  "Genes in analysis denotes the number of pathway genes represented in the ranked expression dataset.",
  "Positive NES values indicate enrichment toward patients experiencing a PFS event within two years;",
  "negative values indicate enrichment toward the PFS-event-free group.",
  "BH q-values were calculated across all 50 Hallmark tests.",
  "BH, Benjamini-Hochberg; NES, normalised enrichment score; PFS, progression-free survival."
)

write_table_docx(
  ft = ft_s1,
  title = paste(
    "Supplementary Table S1.",
    "Complete ranked Hallmark gene-set enrichment results."
  ),
  footnote = footnote_s1,
  output_path =
    "supplementary_tables/Supplementary_Table_S1_all50_fgsea.docx"
)
