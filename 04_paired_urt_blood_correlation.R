############################
# STEP 0: INSTALL + LOAD PACKAGES
############################

cran_packages <- c(
  "readxl",
  "dplyr",
  "tibble",
  "tidyr",
  "purrr",
  "stringr",
  "ggplot2",
  "writexl",
  "scales"
)

bioc_packages <- c(
  "edgeR",
  "GSVA"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(readxl)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(writexl)
library(scales)
library(edgeR)
library(GSVA)


############################
# STEP 1: INPUT FILE NAMES
############################

metadata_file <- "BRAVE_RNASeq_metadata.xlsx"
counts_file   <- "GSE231409_BRAVE_RNASeq_counts.csv"
gmt_file      <- "modules_61.gmt"
limma_file    <- "annotated_limma_results_under21_oldLogic.xlsx"


############################
# STEP 2: SELECT MAIN FIGURE PATHWAYS
############################

selected_modules <- c(
  "Innate Immune Cell Activation",
  "Interferon Response",
  "Type I Interferon Signaling",
  "Type II Interferon Signaling",
  "RNA Sensing",
  "Phagocytosis",
  "Myeloid Activation",
  "Adaptive Immune Response",
  "TCR Signaling",
  "BCR Signaling",
  "NK Activity",
  "Mononuclear Cell Migration",
  "Lymphocyte Trafficking",
  "Immune Memory"
)

module_label_map <- c(
  "Innate Immune Cell Activation" = "Innate immune cell activation",
  "Interferon Response" = "Interferon response",
  "Type I Interferon Signaling" = "Type I interferon signaling",
  "Type II Interferon Signaling" = "Type II interferon signaling",
  "RNA Sensing" = "RNA sensing",
  "Phagocytosis" = "Phagocytosis",
  "Myeloid Activation" = "Myeloid activation",
  "Adaptive Immune Response" = "Adaptive immune response",
  "TCR Signaling" = "T cell receptor signaling",
  "BCR Signaling" = "B cell receptor signaling",
  "NK Activity" = "NK cell activity",
  "Mononuclear Cell Migration" = "Mononuclear cell migration",
  "Lymphocyte Trafficking" = "Lymphocyte trafficking",
  "Immune Memory" = "Immune memory"
)

plot_fdr_cutoff <- 0.05


############################
# STEP 3: HELPER FUNCTIONS
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

prepare_tissue_expression <- function(count_mat,
                                      sample_metadata,
                                      annot_df,
                                      tissue_label,
                                      module_genes) {
  
  tissue_meta <- sample_metadata %>%
    filter(Compartment == tissue_label) %>%
    arrange(match(alias_study_id, unique(sample_metadata$alias_study_id)))
  
  sample_ids <- tissue_meta$alias_sequencing_id
  
  if (!all(sample_ids %in% colnames(count_mat))) {
    stop(paste("Some", tissue_label, "sample IDs are missing from count matrix columns."))
  }
  
  annot_map <- annot_df %>%
    transmute(
      target_id = as.character(target_id),
      SYMBOL = as.character(SYMBOL)
    ) %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(target_id, .keep_all = TRUE)
  
  keep_ids <- intersect(annot_map$target_id, rownames(count_mat))
  
  expr_counts <- count_mat[keep_ids, sample_ids, drop = FALSE]
  annot_map2  <- annot_map[match(rownames(expr_counts), annot_map$target_id), , drop = FALSE]
  
  collapse_tbl <- data.frame(
    target_id = rownames(expr_counts),
    SYMBOL = annot_map2$SYMBOL,
    mean_count = rowMeans(expr_counts, na.rm = TRUE),
    stringsAsFactors = FALSE
  ) %>%
    arrange(SYMBOL, desc(mean_count)) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  
  expr_counts <- expr_counts[collapse_tbl$target_id, , drop = FALSE]
  rownames(expr_counts) <- collapse_tbl$SYMBOL
  
  expr_counts <- expr_counts[intersect(rownames(expr_counts), module_genes), , drop = FALSE]
  
  dge <- DGEList(counts = expr_counts)
  dge <- calcNormFactors(dge)
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  list(
    expr = logcpm,
    meta = tissue_meta
  )
}

run_ssgsea_robust <- function(expr_mat, pathways_list, min_size = 5) {
  pathways_use <- lapply(pathways_list, function(gs) intersect(gs, rownames(expr_mat)))
  pathways_use <- pathways_use[sapply(pathways_use, length) >= min_size]
  
  if (length(pathways_use) == 0) {
    stop("No pathways passed the minimum overlap threshold.")
  }
  
  out <- tryCatch(
    {
      param_obj <- GSVA::ssgseaParam(
        as.matrix(expr_mat),
        pathways_use,
        minSize = min_size,
        maxSize = Inf,
        alpha = 0.25,
        normalize = TRUE
      )
      GSVA::gsva(param_obj, verbose = FALSE)
    },
    error = function(e) {
      message("Modern GSVA API failed. Trying legacy gsva() syntax...")
      GSVA::gsva(
        as.matrix(expr_mat),
        pathways_use,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = FALSE,
        min.sz = min_size,
        max.sz = Inf,
        verbose = FALSE
      )
    }
  )
  
  as.matrix(out)
}

cor_test_safe <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(c(r = NA, p_value = NA, n = length(x)))
  }
  
  if (sd(x) == 0 || sd(y) == 0) {
    return(c(r = NA, p_value = NA, n = length(x)))
  }
  
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  
  c(
    r = unname(ct$estimate),
    p_value = ct$p.value,
    n = length(x)
  )
}


############################
# STEP 4: READ INPUT DATA
############################

metadata_brave <- read_excel(metadata_file, sheet = "brave")
metadata_messi <- read_excel(metadata_file, sheet = "messi")
metadata <- bind_rows(metadata_brave, metadata_messi)

counts_df <- read.csv(counts_file, check.names = FALSE)
count_mat <- as.matrix(counts_df[, -1])
rownames(count_mat) <- counts_df$target_id
storage.mode(count_mat) <- "numeric"

module61 <- read_gmt_simple(gmt_file)

pb_annot  <- read_excel(limma_file, sheet = "Peripheral_Blood_annotated")
urt_annot <- read_excel(limma_file, sheet = "Upper_Respiratory_annotated")


############################
# STEP 5: CHECK METADATA AND COUNTS MATCH
############################

cat("Metadata rows:", nrow(metadata), "\n")
cat("Count matrix samples:", ncol(count_mat), "\n")

missing_in_counts <- setdiff(metadata$alias_sequencing_id, colnames(count_mat))
missing_in_meta   <- setdiff(colnames(count_mat), metadata$alias_sequencing_id)

cat("Metadata sample IDs missing from count matrix:", length(missing_in_counts), "\n")
cat("Count matrix sample IDs missing from metadata:", length(missing_in_meta), "\n")

if (length(missing_in_counts) > 0 || length(missing_in_meta) > 0) {
  stop("Metadata and counts do not match cleanly. Check sample IDs.")
}


############################
# STEP 6: FILTER TO FIGURE 5 COHORT
############################

metadata_f5 <- metadata %>%
  mutate(
    SampleType   = tolower(SampleType),
    SampleTiming = tolower(SampleTiming),
    corona       = tolower(corona)
  ) %>%
  filter(
    age < 21,
    SampleTiming == "acute",
    corona == "positive",
    SampleType %in% c("pax", "np", "nasal")
  ) %>%
  mutate(
    Compartment = if_else(SampleType == "pax", "PAX", "NSB")
  )

cat("\nFiltered Figure 5 sample counts:\n")
print(table(metadata_f5$SampleType))
print(table(metadata_f5$Compartment))


############################
# STEP 7: KEEP ONLY PAIRED NSB–PAX PARTICIPANTS
############################

pair_summary <- metadata_f5 %>%
  group_by(alias_study_id) %>%
  summarise(
    has_PAX = any(Compartment == "PAX"),
    has_NSB = any(Compartment == "NSB"),
    n_PAX = sum(Compartment == "PAX"),
    n_NSB = sum(Compartment == "NSB"),
    URT_type = paste(sort(unique(SampleType[Compartment == "NSB"])), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(has_PAX, has_NSB)

cat("\nNumber of paired participants:", nrow(pair_summary), "\n")

if (any(pair_summary$n_PAX != 1 | pair_summary$n_NSB != 1)) {
  print(pair_summary %>% filter(n_PAX != 1 | n_NSB != 1))
  stop("Some participants do not have exactly one PAX and one NSB sample.")
}

participant_order <- pair_summary$alias_study_id

metadata_paired <- metadata_f5 %>%
  filter(alias_study_id %in% participant_order) %>%
  arrange(match(alias_study_id, participant_order), Compartment)

cat("\nPaired sample counts by raw upper respiratory type:\n")
print(table(metadata_paired$SampleType))


############################
# STEP 8: READ MODULE 61 AND KEEP SELECTED PANEL
############################

module61_all_genes <- unique(unlist(module61))

missing_modules <- setdiff(selected_modules, names(module61))
if (length(missing_modules) > 0) {
  stop(
    paste(
      "These selected modules were not found in modules_61.gmt:",
      paste(missing_modules, collapse = ", ")
    )
  )
}

module61_plot <- module61[selected_modules]

cat("\nSelected plotting modules:\n")
print(names(module61_plot))


############################
# STEP 9: BUILD TISSUE-SPECIFIC EXPRESSION MATRICES
# Uses the fixed tissue-specific under-21 limma gene universe
############################

urt_obj <- prepare_tissue_expression(
  count_mat = count_mat,
  sample_metadata = metadata_paired,
  annot_df = urt_annot,
  tissue_label = "NSB",
  module_genes = module61_all_genes
)

pb_obj <- prepare_tissue_expression(
  count_mat = count_mat,
  sample_metadata = metadata_paired,
  annot_df = pb_annot,
  tissue_label = "PAX",
  module_genes = module61_all_genes
)

cat("\nExpression matrix dimensions after tissue-specific filtering:\n")
cat("NSB:", dim(urt_obj$expr)[1], "genes x", dim(urt_obj$expr)[2], "samples\n")
cat("PAX:", dim(pb_obj$expr)[1], "genes x", dim(pb_obj$expr)[2], "samples\n")


############################
# STEP 10: RUN ssGSEA IN NSB AND PAX
############################

ssgsea_nsb <- run_ssgsea_robust(
  expr_mat = urt_obj$expr,
  pathways_list = module61_plot,
  min_size = 5
)

ssgsea_pax <- run_ssgsea_robust(
  expr_mat = pb_obj$expr,
  pathways_list = module61_plot,
  min_size = 5
)

common_modules <- selected_modules[
  selected_modules %in% rownames(ssgsea_nsb) &
    selected_modules %in% rownames(ssgsea_pax)
]

ssgsea_nsb <- ssgsea_nsb[common_modules, , drop = FALSE]
ssgsea_pax <- ssgsea_pax[common_modules, , drop = FALSE]

colnames(ssgsea_nsb) <- urt_obj$meta$alias_study_id
colnames(ssgsea_pax) <- pb_obj$meta$alias_study_id

common_pairs <- participant_order[
  participant_order %in% colnames(ssgsea_nsb) &
    participant_order %in% colnames(ssgsea_pax)
]

ssgsea_nsb <- ssgsea_nsb[, common_pairs, drop = FALSE]
ssgsea_pax <- ssgsea_pax[, common_pairs, drop = FALSE]

cat("\nFinal paired participants used in correlations:", length(common_pairs), "\n")


############################
# STEP 11: CORRELATE NSB MODULES AGAINST PAX MODULES
############################

cor_list <- list()
k <- 1

for (urt_mod in common_modules) {
  for (pax_mod in common_modules) {
    
    tmp <- cor_test_safe(
      x = as.numeric(ssgsea_nsb[urt_mod, common_pairs]),
      y = as.numeric(ssgsea_pax[pax_mod, common_pairs])
    )
    
    cor_list[[k]] <- data.frame(
      URT_Module = urt_mod,
      PAX_Module = pax_mod,
      r = as.numeric(tmp["r"]),
      p_value = as.numeric(tmp["p_value"]),
      n = as.integer(tmp["n"]),
      stringsAsFactors = FALSE
    )
    
    k <- k + 1
  }
}

cor_results <- bind_rows(cor_list) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    abs_r = abs(r),
    significant = if_else(!is.na(padj) & padj < plot_fdr_cutoff, TRUE, FALSE)
  )

cat("\nNumber of BH-significant correlations:", sum(cor_results$significant, na.rm = TRUE), "\n")


############################
# STEP 12: BUILD PLOT TABLE
############################

x_labels <- module_label_map[common_modules]
y_labels <- rev(module_label_map[common_modules])

plot_df <- cor_results %>%
  mutate(
    PAX_Label = factor(module_label_map[PAX_Module], levels = x_labels),
    URT_Label = factor(module_label_map[URT_Module], levels = y_labels)
  )

sig_points <- plot_df %>%
  filter(significant, !is.na(r))


############################
# STEP 13: PLOT MAIN FIGURE
############################

figure5_plot <- ggplot(plot_df, aes(x = PAX_Label, y = URT_Label)) +
  geom_tile(
    fill = NA,
    color = "#CDCDCD",
    linewidth = 0.4
  ) +
  geom_point(
    data = sig_points,
    aes(size = abs_r, fill = r),
    shape = 21,
    color = "#4A4A4A",
    stroke = 0.25
  ) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(
    range = c(2.6, 8.4),
    limits = c(0, 1),
    guide = "none"
  ) +
  scale_fill_stepsn(
    colours = c("#B5641D", "#E4AE57", "#FAFAFA", "#B8B0D9", "#6D57AF"),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.2),
    labels = c("-1", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1"),
    name = NULL
  ) +
  coord_fixed() +
  labs(
    title = "Peripheral Blood",
    x = NULL,
    y = "Upper respiratory"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.background  = element_rect(fill = "#FAFAFA", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 18,
      margin = margin(r = 24),
      color = "black"
    ),
    
    axis.text.x.top = element_text(
      angle = 42,
      hjust = 0,
      vjust = 0.15,
      size = 10.5,
      color = "#3F3F3F",
      margin = margin(b = 2)
    ),
    axis.text.y = element_text(
      size = 10.5,
      color = "#3F3F3F"
    ),
    axis.ticks = element_blank(),
    
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "plain",
      margin = margin(b = 20)
    ),
    
    legend.position = "right",
    legend.text = element_text(size = 13),
    
    plot.margin = margin(t = 18, r = 18, b = 18, l = 18)
  ) +
  guides(
    fill = guide_coloursteps(
      frame.colour = "black",
      ticks.colour = "black",
      show.limits = TRUE,
      even.steps = TRUE,
      barwidth = grid::unit(10, "mm"),
      barheight = grid::unit(118, "mm")
    )
  )

print(figure5_plot)


############################
# STEP 14: SAVE FIGURE
############################

ggsave(
  filename = "Figure5_under21_NSB_vs_PAX_module_correlation_reducedPanel_main.png",
  plot = figure5_plot,
  width = 12,
  height = 8.6,
  dpi = 300,
  bg = "#FAFAFA"
)


############################
# STEP 15: SAVE OUTPUT TABLES
############################

cohort_summary <- tibble(
  Metric = c(
    "Age filter",
    "Timing filter",
    "COVID filter",
    "Upper respiratory definition",
    "Final upper respiratory label",
    "Final peripheral blood label",
    "Paired participants",
    "NSB samples",
    "PAX samples",
    "NSB raw sample type: np",
    "NSB raw sample type: nasal",
    "Reduced module panel size",
    "Correlation method",
    "Multiple testing correction",
    "Plotting cutoff"
  ),
  Value = c(
    "< 21",
    "acute only",
    "Positive only",
    "np + nasal",
    "NSB",
    "PAX",
    length(common_pairs),
    sum(metadata_paired$Compartment == "NSB"),
    sum(metadata_paired$Compartment == "PAX"),
    sum(metadata_paired$SampleType == "np"),
    sum(metadata_paired$SampleType == "nasal"),
    length(common_modules),
    "Pearson correlation",
    "BH / FDR",
    paste0("padj < ", plot_fdr_cutoff)
  )
)

pair_metadata_sheet <- metadata_paired %>%
  select(
    alias_study_id,
    alias_sequencing_id,
    SampleType,
    Compartment,
    age,
    corona,
    SampleTiming
  )

ssgsea_nsb_df <- as.data.frame(ssgsea_nsb) %>%
  rownames_to_column("Module")

ssgsea_pax_df <- as.data.frame(ssgsea_pax) %>%
  rownames_to_column("Module")

plot_table <- plot_df %>%
  select(
    URT_Module,
    PAX_Module,
    URT_Label,
    PAX_Label,
    r,
    abs_r,
    p_value,
    padj,
    n,
    significant
  )

write_xlsx(
  list(
    Cohort_Summary = cohort_summary,
    Pair_Metadata = pair_metadata_sheet,
    ssGSEA_NSB = ssgsea_nsb_df,
    ssGSEA_PAX = ssgsea_pax_df,
    Correlation_Results = cor_results,
    Plot_Table = plot_table
  ),
  path = "Figure5_under21_NSB_vs_PAX_module_correlation_reducedPanel_main.xlsx"
)


############################
# STEP 16: OPTIONAL CONSOLE OUTPUT
############################

cat("\nTop positive correlations:\n")
print(
  cor_results %>%
    arrange(desc(r)) %>%
    select(URT_Module, PAX_Module, r, padj, n) %>%
    head(15)
)

cat("\nTop negative correlations:\n")
print(
  cor_results %>%
    arrange(r) %>%
    select(URT_Module, PAX_Module, r, padj, n) %>%
    head(15)
)

cat("\nDone.\n")
cat("Saved figure: Figure5_under21_NSB_vs_PAX_module_correlation_reducedPanel_main.png\n")
cat("Saved workbook: Figure5_under21_NSB_vs_PAX_module_correlation_reducedPanel_main.xlsx\n")
