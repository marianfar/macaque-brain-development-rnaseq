#!/usr/bin/env Rscript
# ==============================================================================
# 01_preprocessing_and_qc.R
# Steps: TPM import, biomaRt annotation, protein-coding filtering, PCA & clustering
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(biomaRt)
  library(ggplot2)
  library(pheatmap)
  library(viridis)
})

cat("=== 1. Loading TPM Expression & Metadata ===\n")

# Load metadata
meta_file <- "data/metadata/SraRunTable_new.csv"
if (!file.exists(meta_file)) {
  stop("Metadata file not found at: ", meta_file)
}
meta_raw <- read.csv(meta_file, header = TRUE)
meta <- dplyr::select(meta_raw, any_of(c("Run", "AGE", "Age", "sex", "Sample_Name", "tissue", "Instrument")))
if ("AGE" %in% colnames(meta) && !"Age" %in% colnames(meta)) {
  meta <- dplyr::rename(meta, Age = AGE)
}

# Order age factor
age_levels <- c("E60", "E81", "E82", "E110", "P0", "P2", "7M", "1Y", "2Y", "3.5Y", "4Y", "5Y", "7Y", "11Y")
meta$Age <- factor(meta$Age, levels = age_levels)
meta$AgeGroup <- factor(ifelse(grepl("^E", meta$Age), "Prenatal", "Postnatal"), levels = c("Prenatal", "Postnatal"))

# Load RSEM TPM files if directory exists
rsem_dir <- "results/rsem_genes_results"
if (dir.exists(rsem_dir) && length(list.files(rsem_dir)) > 0) {
  samples <- list.files(rsem_dir)
  expr <- sapply(samples, function(s) {
    quant <- read.csv(file.path(rsem_dir, s), sep = "\t", header = TRUE)
    setNames(quant$TPM, quant$gene_id)
  })
  colnames(expr) <- sub("\\.genes\\.results$", "", colnames(expr))
  expr <- expr[, meta$Run]
} else {
  cat("RSEM directory not found or empty. Using pre-saved expression matrix if available.\n")
}

cat("=== 2. Gene Annotation with biomaRt (Macaque Ensembl) ===\n")
# Query Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "mmulatta_gene_ensembl")
meta_genes_raw <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"),
  filters    = "external_gene_name",
  values     = rownames(expr),
  mart       = ensembl
)

meta_genes <- data.frame(external_gene_name = rownames(expr), stringsAsFactors = FALSE) %>%
  left_join(meta_genes_raw, by = "external_gene_name") %>%
  distinct(external_gene_name, .keep_all = TRUE)

cat("=== 3. Biotype Filtering & Expression Thresholding ===\n")
# Filter to annotated protein-coding genes
keep_pc <- !is.na(meta_genes$gene_biotype) & meta_genes$gene_biotype == "protein_coding"
expr_pc <- expr[keep_pc, ]

# Filter lowly expressed genes (mean >= 1 TPM or detected in >= 50% samples)
expressed_pc <- rowMeans(expr_pc > 0) >= 0.5 | rowMeans(expr_pc) >= 1
expr_pc <- expr_pc[which(expressed_pc), ]

cat(sprintf("Retained %d expressed protein-coding genes across %d samples.\n", nrow(expr_pc), ncol(expr_pc)))

cat("=== 4. Sample Correlations & PCA ===\n")
corr_pearson  <- cor(log1p(expr_pc))
corr_spearman <- cor(expr_pc, method = "spearman")

pca <- prcomp(log1p(t(expr_pc)), center = TRUE, scale. = TRUE)
prop_var <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 1)

# Export cleaned data for downstream scripts
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(
  expr_pc = expr_pc,
  meta = meta,
  corr_pearson = corr_pearson,
  corr_spearman = corr_spearman,
  pca = pca,
  prop_var = prop_var
), file = "results/tables/01_preprocessed_data.rds")

cat("Preprocessed data successfully saved to results/tables/01_preprocessed_data.rds\n")
