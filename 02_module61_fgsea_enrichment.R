# =========================
# 0. Load packages
# =========================
cran_pkgs <- c(
  "readxl",
  "writexl",
  "dplyr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "stringr"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE)
}

library(readxl)
library(writexl)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(fgsea)

# =========================
# 1. File paths
# =========================
limma_file <- "annotated_limma_results_under21_oldLogic.xlsx"
gmt_file   <- "modules_61.gmt"

outdir <- "GSEA_under21_dissertation_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(limma_file)) stop("Missing file: ", limma_file)
if (!file.exists(gmt_file)) stop("Missing file: ", gmt_file)

# =========================
# 2. Read under-21 annotated limma sheets
# =========================
pb_df <- read_excel(limma_file, sheet = "Peripheral_Blood_annotated")
urt_df <- read_excel(limma_file, sheet = "Upper_Respiratory_annotated")

cat("Peripheral blood rows:", nrow(pb_df), "\n")
cat("Upper respiratory rows:", nrow(urt_df), "\n")

# =========================
# 3. Basic checks
# =========================
required_cols <- c("SYMBOL", "t", "logFC", "adj.P.Val")

for (cc in required_cols) {
  if (!cc %in% colnames(pb_df)) {
    stop("Peripheral_Blood_annotated is missing column: ", cc)
  }
  if (!cc %in% colnames(urt_df)) {
    stop("Upper_Respiratory_annotated is missing column: ", cc)
  }
}

# =========================
# 4. Read Module 61 pathways
# =========================
gene_sets_raw <- fgsea::gmtPathways(gmt_file)

cat("Number of Module 61 pathways:", length(gene_sets_raw), "\n")

module61_genes <- unique(unlist(gene_sets_raw))
module61_genes <- module61_genes[!is.na(module61_genes) & module61_genes != ""]
cat("Unique Module 61 genes:", length(module61_genes), "\n")

# =========================
# 5. Helper: make ranked vector
#    Rank genes by limma moderated t-statistic
# =========================
make_ranks <- function(df, symbol_col = "SYMBOL", stat_col = "t") {
  
  tmp <- df %>%
    mutate(
      SYMBOL_CLEAN = toupper(trimws(as.character(.data[[symbol_col]]))),
      STAT_CLEAN   = suppressWarnings(as.numeric(.data[[stat_col]]))
    ) %>%
    filter(
      !is.na(SYMBOL_CLEAN),
      SYMBOL_CLEAN != "",
      is.finite(STAT_CLEAN)
    )
  
  # If a gene symbol appears more than once, keep the row
  # with the strongest absolute t-statistic
  rank_df <- tmp %>%
    group_by(SYMBOL_CLEAN) %>%
    slice_max(order_by = abs(STAT_CLEAN), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(SYMBOL_CLEAN, STAT_CLEAN)
  
  ranks <- rank_df$STAT_CLEAN
  names(ranks) <- rank_df$SYMBOL_CLEAN
  ranks <- sort(ranks, decreasing = TRUE)
  
  return(ranks)
}

ranks_pb <- make_ranks(pb_df, symbol_col = "SYMBOL", stat_col = "t")
ranks_urt <- make_ranks(urt_df, symbol_col = "SYMBOL", stat_col = "t")

cat("Peripheral blood ranked genes:", length(ranks_pb), "\n")
cat("Upper respiratory ranked genes:", length(ranks_urt), "\n")

# =========================
# 6. Restrict Module 61 pathways
#    to genes present in each ranked list
# =========================
prepare_pathways_for_ranks <- function(gene_sets_raw, ranks, minSize = 5, maxSize = 5000) {
  
  gs <- lapply(gene_sets_raw, function(x) {
    unique(intersect(toupper(x), names(ranks)))
  })
  
  gs <- gs[lengths(gs) >= minSize]
  gs <- gs[lengths(gs) <= maxSize]
  
  return(gs)
}

gene_sets_pb <- prepare_pathways_for_ranks(gene_sets_raw, ranks_pb, minSize = 5, maxSize = 5000)
gene_sets_urt <- prepare_pathways_for_ranks(gene_sets_raw, ranks_urt, minSize = 5, maxSize = 5000)

cat("Peripheral blood pathways retained:", length(gene_sets_pb), "\n")
cat("Upper respiratory pathways retained:", length(gene_sets_urt), "\n")

# =========================
# 7. Run FGSEA
# =========================
run_fgsea_module61 <- function(ranks, gene_sets, tissue_name, minSize = 5, maxSize = 5000) {
  
  fg <- fgseaMultilevel(
    pathways = gene_sets,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize,
    eps      = 0
  ) %>%
    as_tibble() %>%
    mutate(
      tissue = tissue_name,
      FDR = padj,
      neglogFDR = -log10(pmax(FDR, 1e-300)),
      effective_size = lengths(gene_sets[pathway]),
      sig_group = case_when(
        !is.na(FDR) & FDR < 0.05 & NES > 0 ~ "Up_sig",
        !is.na(FDR) & FDR < 0.05 & NES < 0 ~ "Down_sig",
        TRUE                               ~ "Not_sig"
      )
    ) %>%
    arrange(desc(NES), FDR)
  
  return(fg)
}

gsea_pb_u21 <- run_fgsea_module61(
  ranks = ranks_pb,
  gene_sets = gene_sets_pb,
  tissue_name = "Peripheral_Blood"
)

gsea_urt_u21 <- run_fgsea_module61(
  ranks = ranks_urt,
  gene_sets = gene_sets_urt,
  tissue_name = "Upper_Respiratory"
)

cat("Peripheral blood pathways tested:", nrow(gsea_pb_u21), "\n")
cat("Upper respiratory pathways tested:", nrow(gsea_urt_u21), "\n")

# =========================
# 8. Make label tables
#    top 10 positive + top 10 negative significant pathways
# =========================
make_label_df <- function(fg_df, label_n_up = 10, label_n_down = 10) {
  
  lab_up <- fg_df %>%
    filter(sig_group == "Up_sig") %>%
    arrange(desc(NES), FDR) %>%
    slice_head(n = label_n_up)
  
  lab_down <- fg_df %>%
    filter(sig_group == "Down_sig") %>%
    arrange(NES, FDR) %>%
    slice_head(n = label_n_down)
  
  bind_rows(lab_up, lab_down) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    mutate(pathway_label_wrap = str_wrap(pathway, width = 34))
}

labels_pb_u21 <- make_label_df(gsea_pb_u21, label_n_up = 10, label_n_down = 10)
labels_urt_u21 <- make_label_df(gsea_urt_u21, label_n_up = 10, label_n_down = 10)

# =========================
# 9. Plot function
# =========================
plot_gsea_all_pathways <- function(plot_df, lab_df, title_text, fdr_cut = 0.05) {
  
  sig_fill <- grDevices::adjustcolor("#D7AAA5", alpha.f = 0.55)
  sig_line <- "#D26C6C"
  ns_fill  <- grDevices::adjustcolor("grey70", alpha.f = 0.35)
  ns_line  <- grDevices::adjustcolor("grey60", alpha.f = 0.60)
  
  p <- ggplot(plot_df, aes(x = NES, y = neglogFDR)) +
    
    # Non-significant pathways first
    geom_point(
      data = filter(plot_df, sig_group == "Not_sig"),
      aes(size = effective_size),
      shape = 21,
      fill = ns_fill,
      color = ns_line,
      stroke = 0.6
    ) +
    
    # Significant pathways on top
    geom_point(
      data = filter(plot_df, sig_group != "Not_sig"),
      aes(size = effective_size),
      shape = 21,
      fill = sig_fill,
      color = sig_line,
      stroke = 0.9
    ) +
    
    # Reference lines
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dotdash",
      linewidth = 0.5
    ) +
    
    geom_hline(
      yintercept = -log10(fdr_cut),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    # Labels
    geom_label_repel(
      data = lab_df,
      aes(label = pathway_label_wrap),
      size = 3,
      label.size = 0.25,
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = Inf,
      min.segment.length = 0,
      seed = 1
    ) +
    
    scale_size_continuous(range = c(2, 12)) +
    
    labs(
      title = title_text,
      x = "Normalized enrichment score (NES)",
      y = expression(-log[10]~"FDR"),
      size = "Effective size"
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

# =========================
# 10. Make plots
# =========================
p_pb_u21 <- plot_gsea_all_pathways(
  plot_df = gsea_pb_u21,
  lab_df = labels_pb_u21,
  title_text = "Peripheral Blood"
)

p_urt_u21 <- plot_gsea_all_pathways(
  plot_df = gsea_urt_u21,
  lab_df = labels_urt_u21,
  title_text = "Upper Respiratory"
)

print(p_pb_u21)
print(p_urt_u21)

# =========================
# 11. Save plots
# =========================
ggsave(
  filename = file.path(outdir, "GSEA_Peripheral_Blood_under21_dissertation.png"),
  plot = p_pb_u21,
  width = 11,
  height = 6.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(outdir, "GSEA_Upper_Respiratory_under21_dissertation.png"),
  plot = p_urt_u21,
  width = 11,
  height = 6.8,
  dpi = 300,
  bg = "white"
)

# =========================
# 12. Helper to export FGSEA tables cleanly
# =========================
fgsea_export <- function(df) {
  df %>%
    mutate(
      leadingEdge = vapply(
        leadingEdge,
        paste,
        collapse = ", ",
        FUN.VALUE = character(1)
      )
    )
}

# =========================
# 13. Summary table
# =========================
summary_table <- bind_rows(
  gsea_pb_u21 %>%
    summarise(
      tissue = "Peripheral_Blood",
      pathways_tested = n(),
      up_sig = sum(sig_group == "Up_sig", na.rm = TRUE),
      down_sig = sum(sig_group == "Down_sig", na.rm = TRUE),
      not_sig = sum(sig_group == "Not_sig", na.rm = TRUE)
    ),
  gsea_urt_u21 %>%
    summarise(
      tissue = "Upper_Respiratory",
      pathways_tested = n(),
      up_sig = sum(sig_group == "Up_sig", na.rm = TRUE),
      down_sig = sum(sig_group == "Down_sig", na.rm = TRUE),
      not_sig = sum(sig_group == "Not_sig", na.rm = TRUE)
    )
)

# =========================
# 14. Save dissertation Excel workbook
# =========================
write_xlsx(
  list(
    Summary = summary_table,
    Peripheral_Blood_All = fgsea_export(gsea_pb_u21),
    Upper_Respiratory_All = fgsea_export(gsea_urt_u21),
    Peripheral_Blood_Labelled = labels_pb_u21,
    Upper_Respiratory_Labelled = labels_urt_u21
  ),
  path = file.path(outdir, "GSEA_under21_Module61_dissertation.xlsx")
)

# =========================
# 15. Optional: save csv files too
# =========================
write.csv(
  fgsea_export(gsea_pb_u21),
  file = file.path(outdir, "GSEA_Peripheral_Blood_under21_all.csv"),
  row.names = FALSE
)

write.csv(
  fgsea_export(gsea_urt_u21),
  file = file.path(outdir, "GSEA_Upper_Respiratory_under21_all.csv"),
  row.names = FALSE
)

write.csv(
  labels_pb_u21,
  file = file.path(outdir, "GSEA_Peripheral_Blood_under21_labelled.csv"),
  row.names = FALSE
)

write.csv(
  labels_urt_u21,
  file = file.path(outdir, "GSEA_Upper_Respiratory_under21_labelled.csv"),
  row.names = FALSE
)

# =========================
# 16. Done
# =========================
cat("\nDONE\n")
print(list.files(outdir))
