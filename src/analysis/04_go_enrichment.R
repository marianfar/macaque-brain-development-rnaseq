#!/usr/bin/env Rscript
# ==============================================================================
# 04_go_enrichment.R
# Pathway Enrichment & Functional Networks (clusterProfiler & enrichplot)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

cat("=== Loading Data for GO Enrichment ===\n")
prep <- readRDS("results/tables/01_preprocessed_data.rds")
km_res <- readRDS("results/tables/03_kmeans_lrt_results.rds")

universe <- rownames(prep$expr_pc)
clusters <- km_res$km$cluster
k <- 5

cat("=== Running GO Biological Process for K-means Clusters ===\n")
ego_list <- vector("list", k)

for (cl in 1:k) {
  cluster_genes <- names(clusters[clusters == cl])
  cat(sprintf("Enriching Cluster %d (%d genes)...\n", cl, length(cluster_genes)))
  
  ego <- enrichGO(
    gene          = cluster_genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  ego_list[[cl]] <- ego
  
  if (!is.null(ego) && nrow(ego) > 0) {
    p_dot <- dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP: Cluster", cl))
    ggsave(sprintf("results/figures/go_cluster_%d_dotplot.png", cl), p_dot, width = 8, height = 6)
  }
}

saveRDS(ego_list, "results/tables/04_go_cluster_results.rds")
cat("GO enrichment complete. All figures and tables exported.\n")
