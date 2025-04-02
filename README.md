# Aneuploid Embryo scRNA-seq Analysis Pipeline

This repository contains the complete single-cell RNA sequencing (scRNA-seq) analysis workflow for **in vitro time-lapse cultured aneuploid embryos**. The pipeline integrates multiple analytical modules to comprehensively characterize transcriptional and regulatory changes across different embryonic states.

---

## 🧬 Overview of the Analysis Pipeline

1. **Seurat Pipeline**
   - Quality control
   - Normalization and scaling
   - Dimensionality reduction (PCA, UMAP)
   - Clustering and cell type annotation

2. **Differential Expression Analysis**
   - Identification of differential expression genes
   - Comparison between euploid and aneuploid conditions

3. **Repeat Element Analysis**
   - Quantification and visualization of transposable element expression

4. **SCENIC Analysis** 
   - Gene regulatory network inference
   - Transcription factor activity scoring

5. **Developmental Trajectory Inference**  
   - Pseudotime ordering
   - Branch-specific gene expression patterns

6. **Alternative Splicing Analysis**  
   - Identification of differential splicing events
   - Integration with gene expression data

7. **LinGe Aneuploid Blastocyst Dataset**  
   - Identification of differential expression genes
   - Comparative analysis between our data and LinGe's data

---
