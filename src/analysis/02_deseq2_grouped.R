#!/usr/bin/env Rscript
# ==============================================================================
# 02_deseq2_grouped.R
# Differential Expression: Binary Stage Grouping (Wald Test ~ AgeGroup)
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

cat("=== Running DESeq2 (Wald Test: Prenatal vs Postnatal) ===\n")
files <- file.path("results/rsem_genes_results", paste0(meta$Run, ".genes.results"))
names(files) <- meta$Run

txi <- tximport(files, type = "rsem", geneIdCol = "gene_id")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ AgeGroup)
dds_filtered <- dds[intersect(rownames(expr_pc), rownames(dds)), ]

dds_filtered <- DESeq(dds_filtered)
res_DESeq2 <- results(dds_filtered, contrast = c("AgeGroup", "Postnatal", "Prenatal"))
res_df <- as.data.frame(res_DESeq2) %>% filter(!is.na(padj))

# Classify significance (|log2FC| > 2, FDR < 0.05)
res_df$group <- "NS"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange >  2] <- "Up_Postnatal"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < -2] <- "Up_Prenatal"
res_df$group <- factor(res_df$group, levels = c("NS", "Up_Postnatal", "Up_Prenatal"))

cat("Summary of DEGs (FDR < 0.05, |log2FC| > 2):\n")
print(table(res_df$group))

# Save DEG results table
write.csv(res_df, "results/tables/deg_grouped_wald_test.csv", row.names = TRUE)

# Volcano Plot
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red3", "dodgerblue3")) +
  labs(title = "Volcano Plot: Postnatal vs Prenatal (Wald Test)",
       x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_minimal()

ggsave("results/figures/volcano_prenatal_vs_postnatal.png", p_volcano, width = 7, height = 5)

# Save results for GO analysis
saveRDS(res_df, "results/tables/02_res_grouped_df.rds")
cat("Grouped analysis complete. Figures saved to results/figures/\n")
