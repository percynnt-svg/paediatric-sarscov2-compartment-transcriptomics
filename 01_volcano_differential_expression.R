############################
# 0) INSTALL PACKAGES
############################
# Run once if needed

# install.packages(c(
#   "dplyr", "ggplot2", "ggrepel", "msigdbr",
#   "conflicted", "readxl", "writexl"
# ))
#
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
#
# BiocManager::install(c(
#   "GEOquery", "edgeR", "limma",
#   "AnnotationDbi", "org.Hs.eg.db"
# ))


############################
# 1) LOAD PACKAGES
############################
library(GEOquery)
library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(conflicted)
library(readxl)
library(writexl)
library(AnnotationDbi)
library(org.Hs.eg.db)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("distinct", "dplyr")
conflict_prefer("left_join", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("slice_head", "dplyr")
conflict_prefer("transmute", "dplyr")


############################
# 2) SETTINGS
############################
lfc_cut <- 1
fdr_cut <- 0.05

N_LABEL_UP <- 10
N_LABEL_DOWN <- 5

counts_file <- "GSE231409_BRAVE_RNASeq_counts.csv"
meta_file   <- "BRAVE_RNASeq_metadata.xlsx"

outdir <- file.path(getwd(), "volcano_under21_oldlogic_outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


############################
# 3) READ COUNTS
############################
counts_df <- read.csv(counts_file, check.names = FALSE)
stopifnot("target_id" %in% colnames(counts_df))

counts_mat <- as.matrix(counts_df[, -1])
mode(counts_mat) <- "numeric"
rownames(counts_mat) <- counts_df$target_id

sample_names <- colnames(counts_mat)

cat("Counts data frame dimensions:\n")
print(dim(counts_df))

cat("Counts matrix dimensions:\n")
print(dim(counts_mat))

cat("First 10 sample names:\n")
print(head(sample_names, 10))


############################
# 4) PARSE SAMPLE NAMES
#    Old logic: SUBJECT.TIMEPOINT.TISSUE
############################
sp <- strsplit(sample_names, "\\.")

sample_info <- data.frame(
  sample = sample_names,
  subject = vapply(sp, `[`, "", 1),
  timepoint = vapply(sp, `[`, "", 2),
  tissue_code = vapply(sp, `[`, "", 3),
  stringsAsFactors = FALSE
)

sample_info$tissue <- ifelse(
  sample_info$tissue_code == "NSB", "NSB",
  ifelse(sample_info$tissue_code == "PAX", "PAX", NA)
)

cat("\nTissue assignment from sample names:\n")
print(table(sample_info$tissue, useNA = "ifany"))


############################
# 5) PULL COVID STATUS FROM GEO
#    Old logic
############################
gse <- getGEO("GSE231409", GSEMatrix = TRUE)
pd <- pData(gse[[1]])

char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)

parse_characteristics <- function(chars) {
  chars <- chars[!is.na(chars)]
  out <- list()
  for (x in chars) {
    if (grepl(":", x, fixed = TRUE)) {
      key <- tolower(trimws(sub(":.*$", "", x)))
      val <- trimws(sub("^.*?:", "", x))
      out[[key]] <- val
    }
  }
  out
}

idx <- match(sample_info$sample, pd$title)

if (any(is.na(idx))) {
  cat("\nUnmatched samples against GEO pd$title:\n")
  print(head(sample_info$sample[is.na(idx)], 20))
  stop("Some samples did not match GEO pd$title.")
}

covid <- character(nrow(sample_info))

for (i in seq_len(nrow(sample_info))) {
  kv <- parse_characteristics(unlist(pd[idx[i], char_cols], use.names = FALSE))
  
  covid_raw <- kv[["covid"]]
  if (is.null(covid_raw)) covid_raw <- kv[["covid cat"]]
  if (is.null(covid_raw)) covid_raw <- kv[["covid_cat"]]
  
  if (is.null(covid_raw) || is.na(covid_raw)) {
    covid[i] <- NA
  } else if (grepl("pos|positive|yes|true", covid_raw, ignore.case = TRUE)) {
    covid[i] <- "pos"
  } else if (grepl("neg|negative|no|false", covid_raw, ignore.case = TRUE)) {
    covid[i] <- "neg"
  } else {
    covid[i] <- NA
  }
}

sample_info$covid <- covid

cat("\nCOVID labels from GEO metadata:\n")
print(table(sample_info$covid, useNA = "ifany"))

cat("\nTissue x COVID table:\n")
print(table(sample_info$tissue, sample_info$covid, useNA = "ifany"))


############################
# 6) READ LOCAL METADATA
#    For age filter
############################
meta <- read_excel(meta_file)

cat("\nMetadata columns:\n")
print(colnames(meta))

meta_age <- meta %>%
  transmute(
    sample = alias_sequencing_id,
    age_num = suppressWarnings(as.numeric(age)),
    age_cat = age_cat
  )

sample_info <- sample_info %>%
  left_join(meta_age, by = "sample")

cat("\nMissing age counts:\n")
print(table(is.na(sample_info$age_num)))

cat("\nAge summary:\n")
print(summary(sample_info$age_num))


############################
# 7) FINAL UNDER-21 COHORT
############################
sample_info_u21 <- sample_info %>%
  filter(
    !is.na(tissue),
    tissue %in% c("NSB", "PAX"),
    covid %in% c("pos", "neg"),
    !is.na(age_num),
    age_num < 21
  )

cat("\nFinal under-21 cohort used:\n")
print(with(sample_info_u21, table(tissue, covid)))

cat("\nTimepoint distribution:\n")
print(with(sample_info_u21, table(tissue, covid, timepoint)))


############################
# 8) SPLIT INTO PAX AND NSB
############################
meta_pax_u21 <- sample_info_u21 %>%
  filter(tissue == "PAX") %>%
  mutate(covid = factor(covid, levels = c("neg", "pos")))

meta_nsb_u21 <- sample_info_u21 %>%
  filter(tissue == "NSB") %>%
  mutate(covid = factor(covid, levels = c("neg", "pos")))

cat("\nPAX dimensions:\n")
print(dim(meta_pax_u21))

cat("\nNSB dimensions:\n")
print(dim(meta_nsb_u21))

cat("\nPAX covid counts:\n")
print(table(meta_pax_u21$covid))

cat("\nNSB covid counts:\n")
print(table(meta_nsb_u21$covid))

cat("\nDo PAX and NSB sample IDs match exactly?\n")
print(identical(meta_pax_u21$sample, meta_nsb_u21$sample))


############################
# 9) RUN LIMMA-VOOM FOR PAX
############################
y_pax <- edgeR::DGEList(
  counts = counts_mat[, meta_pax_u21$sample, drop = FALSE]
)

keep_pax <- edgeR::filterByExpr(y_pax, group = meta_pax_u21$covid)
y_pax <- y_pax[keep_pax, , keep.lib.sizes = FALSE]
y_pax <- edgeR::calcNormFactors(y_pax)

design_pax <- model.matrix(~0 + meta_pax_u21$covid)
colnames(design_pax) <- levels(meta_pax_u21$covid)

v_pax <- limma::voom(y_pax, design_pax, plot = FALSE)
fit_pax <- limma::lmFit(v_pax, design_pax)
cont_pax <- limma::makeContrasts(pos_vs_neg = pos - neg, levels = design_pax)
fit2_pax <- limma::eBayes(limma::contrasts.fit(fit_pax, cont_pax))

pax_all <- limma::topTable(
  fit2_pax,
  coef = "pos_vs_neg",
  number = Inf,
  sort.by = "P"
)

pax_all$target_id <- rownames(pax_all)
rownames(pax_all) <- NULL

cat("\nPAX counts after filterByExpr:\n")
print(dim(y_pax$counts))

cat("\nPAX limma result dimensions:\n")
print(dim(pax_all))


############################
# 10) RUN LIMMA-VOOM FOR NSB
############################
y_nsb <- edgeR::DGEList(
  counts = counts_mat[, meta_nsb_u21$sample, drop = FALSE]
)

keep_nsb <- edgeR::filterByExpr(y_nsb, group = meta_nsb_u21$covid)
y_nsb <- y_nsb[keep_nsb, , keep.lib.sizes = FALSE]
y_nsb <- edgeR::calcNormFactors(y_nsb)

design_nsb <- model.matrix(~0 + meta_nsb_u21$covid)
colnames(design_nsb) <- levels(meta_nsb_u21$covid)

v_nsb <- limma::voom(y_nsb, design_nsb, plot = FALSE)
fit_nsb <- limma::lmFit(v_nsb, design_nsb)
cont_nsb <- limma::makeContrasts(pos_vs_neg = pos - neg, levels = design_nsb)
fit2_nsb <- limma::eBayes(limma::contrasts.fit(fit_nsb, cont_nsb))

nsb_all <- limma::topTable(
  fit2_nsb,
  coef = "pos_vs_neg",
  number = Inf,
  sort.by = "P"
)

nsb_all$target_id <- rownames(nsb_all)
rownames(nsb_all) <- NULL

cat("\nNSB counts after filterByExpr:\n")
print(dim(y_nsb$counts))

cat("\nNSB limma result dimensions:\n")
print(dim(nsb_all))


############################
# 11) CHECK PAX AND NSB ARE DIFFERENT
############################
cat("\nAre PAX and NSB logFC identical?\n")
print(identical(pax_all$logFC, nsb_all$logFC))

cat("\nAre PAX and NSB adj.P.Val identical?\n")
print(identical(pax_all$adj.P.Val, nsb_all$adj.P.Val))


############################
# 12) BUILD IMMUNE GENE LIST
#     Old logic
############################
ms_h <- msigdbr(species = "Homo sapiens", category = "H")

immune_hallmark <- ms_h %>%
  filter(gs_name %in% c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  )) %>%
  pull(gene_symbol) %>%
  unique()

ms_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

immune_bp <- ms_bp %>%
  filter(grepl(
    "IMMUNE|IMMUNITY|INNATE|ADAPTIVE|LEUKOCYTE|LYMPHOCYTE|T_CELL|B_CELL|INTERFERON|CYTOKINE|INFLAM",
    gs_name,
    ignore.case = TRUE
  )) %>%
  pull(gene_symbol) %>%
  unique()

immune_genes <- unique(c(immune_hallmark, immune_bp))

cat("\nNumber of immune genes:\n")
print(length(immune_genes))

cat("\nIFI44L in immune gene list?\n")
print("IFI44L" %in% immune_genes)


############################
# 13) PREP VOLCANO DATA
############################
prep_volcano <- function(res_all, immune_genes_vec, lfc_cut = 1, fdr_cut = 0.05) {
  
  res_all <- res_all %>%
    mutate(
      ensembl_nover = sub("\\..*$", "", target_id),
      adjP_safe = pmax(adj.P.Val, 1e-300),
      neglogFDR = -log10(adjP_safe)
    )
  
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(res_all$ensembl_nover),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  ) %>%
    distinct(ENSEMBL, .keep_all = TRUE)
  
  res_plot <- res_all %>%
    left_join(map, by = c("ensembl_nover" = "ENSEMBL")) %>%
    mutate(
      label_name = ifelse(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, ensembl_nover),
      is_immune = !is.na(SYMBOL) & (SYMBOL %in% immune_genes_vec),
      status = case_when(
        is_immune & adj.P.Val < fdr_cut & logFC >=  lfc_cut ~ "Up",
        is_immune & adj.P.Val < fdr_cut & logFC <= -lfc_cut ~ "Down",
        TRUE ~ "Other"
      )
    )
  
  return(res_plot)
}

pax_plotdf <- prep_volcano(pax_all, immune_genes, lfc_cut = lfc_cut, fdr_cut = fdr_cut)
nsb_plotdf <- prep_volcano(nsb_all, immune_genes, lfc_cut = lfc_cut, fdr_cut = fdr_cut)

cat("\nPAX immune-pass status counts:\n")
print(table(pax_plotdf$status))

cat("\nNSB immune-pass status counts:\n")
print(table(nsb_plotdf$status))


############################
# 14) FINAL VOLCANO PLOT FUNCTION
############################
plot_volcano_final <- function(df, title_text, out_png,
                               outdir = "volcano_under21_oldlogic_outputs",
                               lfc_cut = 1, fdr_cut = 0.05,
                               n_up = 10, n_down = 5,
                               x_limits = NULL,
                               y_limits = NULL) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  top_up <- df %>%
    filter(status == "Up") %>%
    arrange(desc(logFC), adj.P.Val) %>%
    slice_head(n = n_up)
  
  top_down <- df %>%
    filter(status == "Down") %>%
    arrange(logFC, adj.P.Val) %>%
    slice_head(n = n_down)
  
  lab_df <- bind_rows(top_up, top_down) %>%
    distinct(label_name, .keep_all = TRUE)
  
  p <- ggplot(df, aes(x = logFC, y = neglogFDR, color = status)) +
    geom_point(
      data = df %>% filter(status == "Other"),
      alpha = 0.35,
      size = 1.5
    ) +
    geom_point(
      data = df %>% filter(status == "Up"),
      alpha = 0.95,
      size = 2.0
    ) +
    geom_point(
      data = df %>% filter(status == "Down"),
      alpha = 0.95,
      size = 2.0
    ) +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut),
      linetype = "dashed",
      linewidth = 0.7,
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(fdr_cut),
      linetype = "dashed",
      linewidth = 0.7,
      color = "black"
    ) +
    ggrepel::geom_label_repel(
      data = lab_df,
      aes(label = label_name),
      size = 4.0,
      fontface = "plain",
      label.size = 0.3,
      fill = "white",
      box.padding = 0.35,
      point.padding = 0.25,
      segment.size = 0.5,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c(
        "Down" = "blue",
        "Other" = "grey75",
        "Up" = "red"
      )
    ) +
    labs(
      title = title_text,
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ "FDR"),
      color = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 13, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85")
    )
  
  if (!is.null(x_limits) & !is.null(y_limits)) {
    p <- p + coord_cartesian(xlim = x_limits, ylim = y_limits)
  } else if (!is.null(x_limits)) {
    p <- p + coord_cartesian(xlim = x_limits)
  } else if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  ggsave(
    filename = file.path(outdir, out_png),
    plot = p,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  return(list(
    plot = p,
    top_up = top_up,
    top_down = top_down
  ))
}


############################
# 15) MAKE FINAL VOLCANO PLOTS
############################
blood_volcano_final <- plot_volcano_final(
  df = pax_plotdf,
  title_text = "Peripheral Blood",
  out_png = "volcano_peripheral_blood_under21_oldLogic_final.png",
  outdir = outdir,
  lfc_cut = lfc_cut,
  fdr_cut = fdr_cut,
  n_up = 10,
  n_down = 5,
  x_limits = c(-2.2, 4.5),
  y_limits = c(0, 9.5)
)

urt_volcano_final <- plot_volcano_final(
  df = nsb_plotdf,
  title_text = "Upper Respiratory",
  out_png = "volcano_upper_respiratory_under21_oldLogic_final.png",
  outdir = outdir,
  lfc_cut = lfc_cut,
  fdr_cut = fdr_cut,
  n_up = 10,
  n_down = 5,
  x_limits = c(-2.7, 5.0),
  y_limits = c(0, 10.2)
)

blood_volcano_final$plot
urt_volcano_final$plot


############################
# 16) SAVE FULL LIMMA RESULTS
############################
write_xlsx(
  list(
    Peripheral_Blood_limma = pax_all,
    Upper_Respiratory_limma = nsb_all
  ),
  path = file.path(outdir, "full_limma_results_under21_oldLogic.xlsx")
)


############################
# 17) SAVE ANNOTATED LIMMA RESULTS
############################
write_xlsx(
  list(
    Peripheral_Blood_annotated = pax_plotdf,
    Upper_Respiratory_annotated = nsb_plotdf
  ),
  path = file.path(outdir, "annotated_limma_results_under21_oldLogic.xlsx")
)


############################
# 18) SAVE ALL PASSING IMMUNE GENES
############################
passing_immune_pax <- pax_plotdf %>%
  filter(
    is_immune,
    adj.P.Val < fdr_cut,
    abs(logFC) >= lfc_cut
  ) %>%
  arrange(desc(logFC))

passing_immune_nsb <- nsb_plotdf %>%
  filter(
    is_immune,
    adj.P.Val < fdr_cut,
    abs(logFC) >= lfc_cut
  ) %>%
  arrange(desc(logFC))

write_xlsx(
  list(
    Peripheral_Blood = passing_immune_pax,
    Upper_Respiratory = passing_immune_nsb
  ),
  path = file.path(outdir, "all_passing_immune_genes_under21_oldLogic.xlsx")
)


############################
# 19) SAVE TOP LABEL TABLES
############################
write.csv(
  blood_volcano_final$top_up,
  file.path(outdir, "top10_up_peripheral_blood_under21_oldLogic_final.csv"),
  row.names = FALSE
)

write.csv(
  blood_volcano_final$top_down,
  file.path(outdir, "top5_down_peripheral_blood_under21_oldLogic_final.csv"),
  row.names = FALSE
)

write.csv(
  urt_volcano_final$top_up,
  file.path(outdir, "top10_up_upper_respiratory_under21_oldLogic_final.csv"),
  row.names = FALSE
)

write.csv(
  urt_volcano_final$top_down,
  file.path(outdir, "top5_down_upper_respiratory_under21_oldLogic_final.csv"),
  row.names = FALSE
)


############################
# 20) FINAL FILE CHECK
############################
cat("\nFiles saved in output folder:\n")
print(list.files(outdir))
