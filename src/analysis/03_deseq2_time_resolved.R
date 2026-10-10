#!/usr/bin/env Rscript
# ==============================================================================
# 03_deseq2_time_resolved.R
# Dynamic Developmental Modeling: DESeq2 LRT & K-means Clustering (k=5)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(viridis)
})

cat("=== Loading Preprocessed Data ===\n")
prep <- readRDS("results/tables/01_preprocessed_data.rds")
expr_pc <- prep$expr_pc
meta <- prep$meta

cat("=== Running DESeq2 (LRT across developmental ages) ===\n")
files <- file.path("results/rsem_genes_results", paste0(meta$Run, ".genes.results"))
names(files) <- meta$Run

txi <- tximport(files, type = "rsem", geneIdCol = "gene_id")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Age)
dds_filtered <- dds[intersect(rownames(expr_pc), rownames(dds)), ]

dds_lrt <- DESeq(dds_filtered, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
res_df <- as.data.frame(res_lrt) %>% filter(!is.na(padj))

write.csv(res_df, "results/tables/deg_time_resolved_lrt.csv", row.names = TRUE)

cat("=== K-means Clustering on Top 1000 LRT DEGs ===\n")
sig_genes <- res_df %>% filter(padj < 0.01) %>% arrange(padj)
top1000_genes <- rownames(sig_genes)[1:min(1000, nrow(sig_genes))]

expr_sub <- expr_pc[top1000_genes, ]
expr_scaled <- t(scale(t(expr_sub)))

set.seed(123)
k <- 5
km <- kmeans(expr_scaled, centers = k, nstart = 50, iter.max = 20)

cat("K-means cluster sizes:\n")
print(table(km$cluster))

# Save clustering object and genes
saveRDS(list(res_df = res_df, km = km, top1000_genes = top1000_genes, expr_scaled = expr_scaled),
        file = "results/tables/03_kmeans_lrt_results.rds")

cat("Time-resolved analysis complete. Results saved in results/tables/\n")
