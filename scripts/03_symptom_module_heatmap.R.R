############################
# 0) INSTALL PACKAGES (RUN ONCE IF NEEDED)
############################
# install.packages(c("dplyr", "ggplot2", "tibble", "tidyr", "purrr",
#                    "stringr", "readxl", "openxlsx", "BiocManager"))
# BiocManager::install(c("DESeq2", "fgsea", "AnnotationDbi", "org.Hs.eg.db"))

############################
# 1) LOAD PACKAGES
############################
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)
library(openxlsx)
library(DESeq2)
library(fgsea)
library(AnnotationDbi)
library(org.Hs.eg.db)

set.seed(123)

############################
# 2) FILE PATHS
############################
metadata_file     <- "BRAVE_RNASeq_metadata.xlsx"
counts_file       <- "GSE231409_BRAVE_RNASeq_counts.csv"
gmt_file          <- "modules_61.gmt"
limma_filter_file <- "annotated_limma_results_under21_oldLogic.xlsx"

############################
# 3) READ BRAVE METADATA ONLY
############################
meta_all <- readxl::read_excel(metadata_file, sheet = "brave") %>%
  dplyr::mutate(
    sample_id      = alias_sequencing_id,
    age            = as.numeric(age),
    corona         = as.character(corona),
    SampleTiming   = as.character(SampleTiming),
    tissue_from_id = toupper(sub(".*\\.([A-Za-z]+)$", "\\1", sample_id)),
    timepoint      = sub(".*\\.(M[0-9]+)\\..*$", "\\1", sample_id)
  )

cat("BRAVE metadata rows:", nrow(meta_all), "\n")

############################
# 4) READ COUNTS
############################
counts_df <- read.csv(counts_file, check.names = FALSE)

gene_ids <- counts_df[[1]]
count_mat <- as.matrix(counts_df[, -1])
rownames(count_mat) <- gene_ids
storage.mode(count_mat) <- "integer"

############################
# 5) MATCH COUNTS TO METADATA
############################
common_samples <- intersect(colnames(count_mat), meta_all$sample_id)

count_mat <- count_mat[, common_samples, drop = FALSE]

meta_all <- meta_all %>%
  dplyr::filter(sample_id %in% common_samples) %>%
  dplyr::arrange(match(sample_id, colnames(count_mat)))

stopifnot(all(meta_all$sample_id == colnames(count_mat)))

cat("Matched BRAVE samples:", ncol(count_mat), "\n")

############################
# 6) DEFINE FIGURE 4 COHORT
############################
# Final Figure 4 cohort:
# - BRAVE only
# - age < 21
# - acute samples only
# - NSB and PAX only
# - SARS-CoV-2 positive only

meta_fig4 <- meta_all %>%
  dplyr::filter(
    age < 21,
    !is.na(SampleTiming),
    tolower(SampleTiming) == "acute",
    tissue_from_id %in% c("NSB", "PAX"),
    corona == "Positive"
  ) %>%
  dplyr::mutate(
    tissue = tissue_from_id
  )

count_mat_fig4 <- count_mat[, meta_fig4$sample_id, drop = FALSE]

cat("\nFigure 4 cohort sizes:\n")
print(table(meta_fig4$tissue))

cat("\nAge summary:\n")
print(summary(meta_fig4$age))

############################
# 7) DEFINE SYMPTOMS
############################
symptoms_to_use <- c(
  "fever",
  "cough",
  "headache",
  "congestion",
  "rhinorrhea",
  "anosmia",
  "dysgeusia",
  "sorethroat"
)

cat("\nSymptom counts by tissue:\n")
for (tt in c("NSB", "PAX")) {
  cat("\n", tt, "\n", sep = "")
  tmp <- meta_fig4 %>% dplyr::filter(tissue == tt)
  for (sx in symptoms_to_use) {
    cat(sx, ":\n")
    print(table(tmp[[sx]], useNA = "ifany"))
  }
}

############################
# 8) MAP ENSEMBL IDS TO GENE SYMBOLS
############################
ensembl_ids <- sub("\\..*$", "", rownames(count_mat_fig4))

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = unique(ensembl_ids),
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

gene_info <- data.frame(
  ensembl_full = rownames(count_mat_fig4),
  ENSEMBL      = ensembl_ids,
  SYMBOL       = unname(symbol_map[ensembl_ids]),
  stringsAsFactors = FALSE
)

count_df_annot <- cbind(gene_info, as.data.frame(count_mat_fig4, check.names = FALSE))

count_df_annot <- count_df_annot %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "")

sample_cols <- meta_fig4$sample_id

count_df_annot <- count_df_annot %>%
  dplyr::mutate(mean_count = rowMeans(dplyr::across(dplyr::all_of(sample_cols)))) %>%
  dplyr::arrange(dplyr::desc(mean_count)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

count_mat_symbol <- as.matrix(count_df_annot[, sample_cols, drop = FALSE])
rownames(count_mat_symbol) <- count_df_annot$SYMBOL
storage.mode(count_mat_symbol) <- "integer"

cat("\nMapped unique gene symbols:", nrow(count_mat_symbol), "\n")

############################
# 9) READ MODULE 61 GMT + USE FIXED FILTERED GENE UNIVERSE
############################
read_gmt_simple <- function(gmt_path) {
  x <- readLines(gmt_path)
  pathways <- lapply(x, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    genes <- unique(parts[-c(1, 2)])
    genes[genes != ""]
  })
  names(pathways) <- sapply(strsplit(x, "\t"), `[`, 1)
  pathways
}

pathways_raw <- read_gmt_simple(gmt_file)
module61_genes <- sort(unique(unlist(pathways_raw)))

pb_limma <- readxl::read_excel(limma_filter_file, sheet = "Peripheral_Blood_annotated")
ur_limma <- readxl::read_excel(limma_filter_file, sheet = "Upper_Respiratory_annotated")

pb_universe <- pb_limma %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::pull(SYMBOL) %>%
  unique()

ur_universe <- ur_limma %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::pull(SYMBOL) %>%
  unique()

cat("Peripheral blood fixed filtered universe:", length(pb_universe), "\n")
cat("Upper respiratory fixed filtered universe:", length(ur_universe), "\n")

pb_module61_genes <- intersect(pb_universe, module61_genes)
ur_module61_genes <- intersect(ur_universe, module61_genes)

cat("Peripheral blood Module 61 genes in filtered universe:", length(pb_module61_genes), "\n")
cat("Upper respiratory Module 61 genes in filtered universe:", length(ur_module61_genes), "\n")

count_mat_module61_pax <- count_mat_symbol[
  rownames(count_mat_symbol) %in% pb_module61_genes,
  ,
  drop = FALSE
]

count_mat_module61_nsb <- count_mat_symbol[
  rownames(count_mat_symbol) %in% ur_module61_genes,
  ,
  drop = FALSE
]

pathways_module61_pax <- lapply(pathways_raw, function(gs) {
  intersect(gs, rownames(count_mat_module61_pax))
})
pathways_module61_pax <- pathways_module61_pax[lengths(pathways_module61_pax) >= 5]

pathways_module61_nsb <- lapply(pathways_raw, function(gs) {
  intersect(gs, rownames(count_mat_module61_nsb))
})
pathways_module61_nsb <- pathways_module61_nsb[lengths(pathways_module61_nsb) >= 5]

cat("Peripheral blood pathways retained (>=5 genes):", length(pathways_module61_pax), "\n")
cat("Upper respiratory pathways retained (>=5 genes):", length(pathways_module61_nsb), "\n")

############################
# 10) DEFINE MODULE PANELS
############################
innate_modules <- c(
  "Innate Immune Cell Activation",
  "Interferon Response",
  "Type I Interferon Signaling",
  "Type II Interferon Signaling",
  "TNF Signaling",
  "TLR Signaling",
  "NLR Signaling",
  "Chemokine Signaling",
  "RNA Sensing",
  "Phagocytosis",
  "Myeloid Activation",
  "NK Activity"
)

adaptive_modules <- c(
  "Adaptive Immune Response",
  "MHC Class I Antigen Presentation",
  "TCR Signaling",
  "BCR Signaling",
  "NF-kappaB Signaling",
  "Mononuclear Cell Migration",
  "Lymphocyte Trafficking",
  "Immune Memory"
)

selected_modules <- c(innate_modules, adaptive_modules)

############################
# 11) RUN DESEQ2 + FGSEA USING FIXED TISSUE-SPECIFIC GENE UNIVERSE
############################
run_symptom_fgsea <- function(tissue_code, symptom_var, meta_df) {
  
  if (tissue_code == "NSB") {
    count_matrix  <- count_mat_module61_nsb
    pathways_list <- pathways_module61_nsb
  } else if (tissue_code == "PAX") {
    count_matrix  <- count_mat_module61_pax
    pathways_list <- pathways_module61_pax
  } else {
    stop("Unknown tissue_code")
  }
  
  sub_meta <- meta_df %>%
    dplyr::filter(
      tissue == tissue_code,
      .data[[symptom_var]] %in% c("Y", "N")
    ) %>%
    dplyr::mutate(
      symptom_status = factor(.data[[symptom_var]], levels = c("N", "Y"))
    )
  
  group_counts <- table(sub_meta$symptom_status)
  cat("\n", tissue_code, " | ", symptom_var,
      " | N=", group_counts["N"], " Y=", group_counts["Y"], "\n", sep = "")
  
  if (length(group_counts) < 2 || min(group_counts) < 3) {
    cat("Skipped: too few samples in one group\n")
    return(NULL)
  }
  
  sub_counts <- count_matrix[, sub_meta$sample_id, drop = FALSE]
  
  keep_genes <- rowSums(sub_counts) > 0
  sub_counts <- sub_counts[keep_genes, , drop = FALSE]
  
  pathways_sub <- lapply(pathways_list, function(gs) intersect(gs, rownames(sub_counts)))
  pathways_sub <- pathways_sub[lengths(pathways_sub) >= 5]
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData   = as.data.frame(sub_meta),
    design    = ~ symptom_status
  )
  
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  res <- DESeq2::results(
    dds,
    contrast = c("symptom_status", "Y", "N"),
    independentFiltering = FALSE
  )
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("SYMBOL")
  
  ranks <- res_df$stat
  names(ranks) <- res_df$SYMBOL
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  
  fg <- fgsea::fgseaMultilevel(
    pathways = pathways_sub,
    stats    = ranks,
    eps      = 0
  )
  
  fg <- tibble::as_tibble(fg) %>%
    dplyr::mutate(
      effective_size = purrr::map_int(pathway, ~ length(pathways_sub[[.x]])),
      tissue         = tissue_code,
      symptom        = symptom_var,
      n_absent       = sum(sub_meta$symptom_status == "N"),
      n_present      = sum(sub_meta$symptom_status == "Y")
    ) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES)))
  
  list(
    deseq = res_df,
    fgsea = fg
  )
}

all_fgsea <- list()
all_deseq <- list()

for (tt in c("NSB", "PAX")) {
  for (sx in symptoms_to_use) {
    out <- run_symptom_fgsea(
      tissue_code = tt,
      symptom_var = sx,
      meta_df     = meta_fig4
    )
    
    if (!is.null(out)) {
      all_fgsea[[paste(tt, sx, sep = "_")]] <- out$fgsea
      all_deseq[[paste(tt, sx, sep = "_")]] <- out$deseq
    }
  }
}

fgsea_all_df <- dplyr::bind_rows(all_fgsea)

cat("\nTotal FGSEA result rows:", nrow(fgsea_all_df), "\n")

############################
# 12) BUILD PLOT TABLE
############################
symptom_labels <- c(
  fever      = "Fever",
  cough      = "Cough",
  headache   = "Headache",
  congestion = "Congestion",
  rhinorrhea = "Rhinorrhea",
  anosmia    = "Loss of smell",
  dysgeusia  = "Loss of taste",
  sorethroat = "Sore throat"
)

symptom_order <- c(
  "Fever",
  "Cough",
  "Headache",
  "Congestion",
  "Rhinorrhea",
  "Loss of smell",
  "Loss of taste",
  "Sore throat"
)

plot_df <- tidyr::expand_grid(
  tissue  = c("NSB", "PAX"),
  symptom = symptoms_to_use,
  pathway = selected_modules
) %>%
  dplyr::left_join(
    fgsea_all_df %>%
      dplyr::select(tissue, symptom, pathway, NES, padj, effective_size),
    by = c("tissue", "symptom", "pathway")
  ) %>%
  dplyr::mutate(
    module_group = dplyr::case_when(
      pathway %in% innate_modules   ~ "Innate immunity",
      pathway %in% adaptive_modules ~ "Adaptive immunity",
      TRUE ~ NA_character_
    ),
    tissue_label = dplyr::recode(
      tissue,
      "NSB" = "Upper respiratory",
      "PAX" = "Peripheral blood"
    ),
    symptom_label = dplyr::recode(symptom, !!!symptom_labels),
    label = ifelse(!is.na(padj) & padj < 0.05, sprintf("%.2f", NES), "")
  ) %>%
  dplyr::mutate(
    tissue_label = factor(
      tissue_label,
      levels = c("Upper respiratory", "Peripheral blood")
    ),
    module_group = factor(
      module_group,
      levels = c("Innate immunity", "Adaptive immunity")
    ),
    symptom_label = factor(
      symptom_label,
      levels = rev(symptom_order)
    ),
    pathway = factor(
      pathway,
      levels = c(innate_modules, adaptive_modules)
    )
  )

cat("\nPlot table dimensions:\n")
print(dim(plot_df))

############################
# 13) PLOT HEATMAP
############################
p_fig4 <- ggplot(plot_df, aes(x = pathway, y = symptom_label, fill = NES)) +
  geom_tile(width = 0.98, height = 0.98) +
  geom_text(
    aes(label = label),
    size = 3.8,
    color = "black",
    na.rm = TRUE
  ) +
  facet_grid(
    rows = vars(tissue_label),
    cols = vars(module_group),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low      = "#0C6291",
    mid      = "#FBFEF9",
    high     = "#A63446",
    midpoint = 0,
    limits   = c(-4, 4),
    oob      = scales::squish,
    na.value = "grey85",
    name     = "NES"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "grey92", color = NA),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.x   = unit(1.2, "lines"),
    panel.spacing.y   = unit(1.2, "lines"),
    strip.background  = element_blank(),
    strip.placement   = "outside",
    strip.text.x      = element_text(face = "bold", size = 13),
    strip.text.y.left = element_text(face = "bold", size = 13, angle = 90),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    axis.text.y       = element_text(size = 12),
    axis.title        = element_blank(),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 11)
  )

print(p_fig4)

############################
# 14) SAVE FIGURE
############################
ggsave(
  filename = "Figure3_under21_BRAVEonly_acute_symptom_heatmap_Module61_final.png",
  plot     = p_fig4,
  width    = 18,
  height   = 9,
  dpi      = 300,
  bg       = "grey92"
)

cat("Saved figure: Figure3_under21_BRAVEonly_acute_symptom_heatmap_Module61_final.png\n")

############################
# 15) SAVE EXCEL OUTPUTS
############################
cohort_summary <- meta_fig4 %>%
  dplyr::count(tissue, name = "n_samples")

symptom_summary <- purrr::map_dfr(c("NSB", "PAX"), function(tt) {
  tmp <- meta_fig4 %>% dplyr::filter(tissue == tt)
  
  purrr::map_dfr(symptoms_to_use, function(sx) {
    tibble::tibble(
      tissue    = tt,
      symptom   = sx,
      absent_n  = sum(tmp[[sx]] == "N", na.rm = TRUE),
      present_n = sum(tmp[[sx]] == "Y", na.rm = TRUE)
    )
  })
})

fgsea_export <- fgsea_all_df %>%
  dplyr::mutate(
    leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
  )

plot_export <- plot_df %>%
  dplyr::arrange(tissue_label, module_group, symptom_label, pathway)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Cohort_Summary")
openxlsx::writeData(wb, "Cohort_Summary", cohort_summary)

openxlsx::addWorksheet(wb, "Symptom_Summary")
openxlsx::writeData(wb, "Symptom_Summary", symptom_summary)

openxlsx::addWorksheet(wb, "FGSEA_All")
openxlsx::writeData(wb, "FGSEA_All", fgsea_export)

openxlsx::addWorksheet(wb, "Plot_Table")
openxlsx::writeData(wb, "Plot_Table", plot_export)

openxlsx::saveWorkbook(
  wb,
  file = "Figure3_under21_BRAVEonly_acute_symptom_Module61_results_final.xlsx",
  overwrite = TRUE
)

cat("Saved workbook: Figure3_under21_BRAVEonly_acute_symptom_Module61_results_final.xlsx\n")

############################
# 16) OPTIONAL QUICK CHECKS
############################
cat("\nSignificant tiles (padj < 0.05):\n")
print(sum(!is.na(plot_df$padj) & plot_df$padj < 0.05))

cat("\nPer tissue significant tiles:\n")
print(
  plot_df %>%
    dplyr::group_by(tissue_label) %>%
    dplyr::summarise(
      sig_tiles = sum(!is.na(padj) & padj < 0.05),
      .groups = "drop"
    )
)

cat("\nDone.\n")
