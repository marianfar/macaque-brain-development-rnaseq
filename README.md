# Bulk RNA-seq Analysis of Developing vs. Adult Macaque Frontal Cortex

[![Language: R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)
[![Tool: DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2-brightgreen.svg)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
[![Course: ETH Systems Genomics](https://img.shields.io/badge/ETH%20Zürich-Systems%20Genomics%20Fall%202025-blue.svg)](https://ethz.ch/)

End-to-end bulk RNA-seq analysis of *Macaca mulatta* dorsolateral prefrontal cortex (DFC) across lifespan (E60 to 11 years). This project compares two analytical frameworks: **Binary Stage Grouping** (Prenatal vs. Postnatal) and **Time-Resolved Dynamic Modeling** (DESeq2 LRT across ages and *k*-means trajectory clustering).

Completed as part of the ETH Zürich Systems Genomics course (Fall 2025). The complete presentation slide deck is available in [`docs/SystemsGenomics_Presentation.pdf`](docs/SystemsGenomics_Presentation.pdf).


---

## Biological Context & Objectives

Primate cortical development is orchestrated by coordinated transcriptional waves regulating neurogenesis, synaptogenesis, and myelination. This project addresses three primary questions:
* What key genes distinguish prenatal from postnatal cortical development?
* How does gene expression evolve dynamically across continuous developmental ages?
* Does standard binary stage grouping mask stage-specific transition programs that dynamic modeling recovers?

---

## Data Availability & Reference

The raw RNA-seq data are not included in this repository due to file size and licensing restrictions.

Sample metadata and accession numbers are available in [`data/metadata/SraRunTable_new.csv`](data/metadata/SraRunTable_new.csv).

The dataset originates from:
> Zhu Y, Sousa AMM, Gao T, Skarica M, Li M, Santpere G, Esteller-Cucala P, Juan D, Ferrández-Peral L, Gulden FO, Yang M, Miller DJ, Marques-Bonet T, Imamura Kawasawa Y, Zhao H, Sestan N. **Spatiotemporal transcriptomic divergence across human and macaque brain development.** *Science*. 2018 Dec 14;362(6420):eaat8077. doi: [10.1126/science.aat8077](https://doi.org/10.1126/science.aat8077). Epub 2018 Dec 13. PMID: 30545855; PMCID: PMC6900982.

---

## Pipeline & Analysis Strategies

### 1. Upstream Processing
* Quality Control & Standardization: Assessed with FastQC/MultiQC ([`docs/multiqc_report.html`](docs/multiqc_report.html)). Older Illumina 1.5 runs (Phred+64) converted to Phred+33 using `seqtk`.
* Trimming: Adapter trimming and quality filtering via `fastp`.
* Alignment & Quantification: Reads mapped to `rheMac10` with `STAR` and quantified with `RSEM`.
* Annotation: Ensembl gene annotations retrieved with `biomaRt`, filtering uncharacterized `LOC*` genes and non-protein-coding features.

### 2. Grouped Analysis (Prenatal vs. Postnatal)
* Design: `DESeq2` Wald test (`~ AgeGroup`).
* Outcome: Identifies 2,206 DEGs (|log2FC| > 2, FDR < 0.05).
* Biological Recovery: Captures broad transitions (prenatal chromatin/RNA processing vs. postnatal circulation/immune signaling).
* Limitation: Collapses 14 distinct age classes, obscuring intermediate perinatal transitions.

### 3. Time-Resolved Analysis (Dynamic LRT & Clustering)
* Design: `DESeq2` Likelihood Ratio Test (`~ Age` vs. `~ 1`).
* Trajectory Modeling: Unsupervised *k*-means clustering (*k* = 5) on top LRT genes.
* Biological Recovery: Recovers temporal developmental waves: early neurogenesis (Clusters 2 & 3), perinatal synaptogenesis (Cluster 1), synaptic maturation (Cluster 4), and postnatal myelination (Cluster 5).

---

## Repository Structure
```macaque-brain-development-rnaseq/
├── README.md                          # Project documentation
├── .gitignore                         # Excludes raw data and local caches
├── data/
│   └── metadata/
│       └── SraRunTable_new.csv        # Macaque DFC SRA sample metadata
├── docs/
│   ├── SystemsGenomics_Presentation.pdf # Project slides and methodology
│   └── multiqc_report.html            # Post-trim MultiQC report
└── src/
    ├── upstream/
    │   ├── 01_qc_and_trim.sh          # fastp, seqtk, and FastQC pipeline
    │   └── 02_alignment_and_quant.sh  # STAR alignment and RSEM quantification
    └── analysis/
        ├── 01_preprocessing_and_qc.R  # biomaRt annotation, filtering, PCA
        ├── 02_deseq2_grouped.R        # Wald test (Prenatal vs Postnatal)
        ├── 03_deseq2_time_resolved.R  # DESeq2 LRT across ages & k-means clustering
        └── 04_go_enrichment.R         # clusterProfiler GO over-representation
```
---

## Key Tools & Dependencies

- Read Processing: `FastQC`, `MultiQC`, `seqtk`, `fastp`
- Alignment & Quantification: `STAR`, `RSEM`
- Differential Expression & Trajectories: `DESeq2`
- Annotation: `biomaRt`, `org.Mmu.eg.db`
- Functional Enrichment: `clusterProfiler`
- Visualization: `ggplot2`, `pheatmap`

---

## Contributors

- Marian Farrera
- Dea Karam
- Remi Tregan

Systems Genomics Course, Fall 2025 -- ETH Zurich / D-BSSE Basel.
'''

with open('README.md', 'w') as f:
    f.write(content.strip() + '\n')
print('README.md updated successfully.')
"
