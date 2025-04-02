# ========================
# - Title: XuYanWen aneuploid project
# - Dataset: GSA-Human: HRA005378
# - Article: Abnormal lineage differentiation of peri-implantation aneuploid embryos 
# -          revealed by single-cell RNA sequencing
# - Content:
# -         1. Global setting;
# -         2. Library;
# -         3. Seurat pipeline;
# -         4. Differential expression analysis;
# -         5. Repeat analysis;
# -         6. SCENIC analysiswo wowo
# -         7. Infer developmental trajectory;
# -         8. Alternative splicing analysis;
# -         9. LinGe aneuploid blastocysts;
# ========================


# ========================
# 1st part: Global setting ----
# ========================

### >>> 1. Setting workding directory
setwd("/home/yhw/bioinfo/project-xuyanwen/aneuploid")
dir.create("R/CodeData", recursive = T)
dir.create("R/Graphs", recursive = T)
dir.create("R/Table", recursive = T)


### >>> 2. Setting library
pkg.lib <- "/home/laborer/software/anaconda3/envs/rs-4.2.3/lib/R/library"
.libPaths(pkg.lib)



# =================
# 2nd part: Library ----
# =================

### >>> 1. CRAN packages
cran.pks <- c("SCP", "forcats", "dplyr", "tidyr", "tidyverse", "stringr", "circlize", 
              "ggplot2", "ggalluvial", "ggrepel", "RColorBrewer", "cowplot", "Seurat")
for (pks in cran.pks) {
  library(pks, character.only = T)
}


### >>> 2. Bioconductor packages
bioc.pks <- c()
for (pks in bioc.pks) {
  library(pks, character.only = T)
}


### >>> 3. GitHub packages
git.pks <- c("scCustomize", "SCP", "harmony")
for (pks in git.pks) {
  library(pks, character.only = T)
}


### >>> 4. Functions
theme_dp <- theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 13, angle = 0),
                  axis.title.y = element_text(face = "plain", colour = "#000000", size = 13, angle = 90),
                  axis.text.x  = element_text(face = "plain", colour = "#000000", size = 11, angle = 45, hjust = 1),
                  axis.text.y  = element_text(face = "plain", colour = "#000000", size = 11, angle = 0))
WriteLoom.new <- function (sc.obj, loom.file) 
{
  if (!(class(sc.obj)[1] %in% c("Seurat", "SingleCellExperiment"))) {
    stop("[Error] The data type is not supported by this function !!!")
  }
  if (class(sc.obj)[1] == "Seurat") {
    exprMat <- GetAssayData(sc.obj, slot = "count")
    cellInfo <- sc.obj@meta.data
  }
  else if (class(sc.obj)[1] == "SingleCellExperiment") {
    exprMat <- assay(sc.obj, "counts")
    cellInfo <- colData(sc.obj)
  }
  loom <- build_loom(loom.file, dgem = exprMat)
  #loom <- add_cell_annotation(loom, cellInfo)
  close_loom(loom)
}



# =========================
# 3rd part: Seurat pipeline ----
# =========================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Seurat_output")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Load count data
ge.count <- read.table("analysis/scrna/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt", 
                       header = T, sep = "\t", row.names = 1)
colnames(ge.count)[-1:-5] <- gsub("\\.", "_", 
                                  gsub("d", "D", 
                                       gsub("_R1.*", "", 
                                            gsub("X.*gene.", "", colnames(ge.count)[-1:-5]))))
cell.meta.cxy <- read.csv("R/CodeData/meta.csv")
cell.meta <- data.frame(CellID = colnames(ge.count)[-1:-5],
                        Stage = str_split_fixed(colnames(ge.count)[-1:-5], "_", 3)[, 1],
                        Embryo = str_split_fixed(colnames(ge.count)[-1:-5], "_", 3)[, 2],
                        row.names = colnames(ge.count)[-1:-5]) %>% 
  mutate(Karyotype = case_when(Embryo %in% c("A038", "A039", "A040", "A043", "A044", "A045") ~ "T16",
                               Embryo %in% c("A060", "A061", "A062") ~ "M22",
                               Embryo %in% c("A097", "A101") ~ "M16",
                               Embryo %in% c("A095") ~ "M16?",
                               Embryo %in% c("E16", "E18", "E19", "E21", "E22", "E23") ~ "Euploidy"))
table(cell.meta.cxy$embryo, cell.meta.cxy$karyotype)
table(cell.meta$Embryo, cell.meta$Karyotype, cell.meta$Stage)
# 
ge.tpm <- CountToTpm(count = ge.count[, -1:-5], length = ge.count$Length)
if (all(colnames(ge.tpm) == colnames(sr.yhw.no22))) {
  write.csv(ge.tpm, file.path(res.out, "All_sample_gene_expression_level_TPM.csv"), 
            row.names = T, col.names = T, quote = F)
  write.csv(ge.count, file.path(res.out, "All_sample_gene_expression_level_count.csv"), 
            row.names = T, col.names = T, quote = F)
  saveRDS(ge.tpm, file.path(res.out, "All_sample_gene_expression_level_TPM.rds"))
  saveRDS(ge.count, file.path(res.out, "All_sample_gene_expression_level_count.rds"))
}
ge.cpm <- GetAssayData(sr.yhw.no22, slot = "data", assay = "RNA")
write.csv(ge.cpm, file.path(res.out, "All_sample_gene_expression_level_log-normalized_count.csv"), 
          row.names = T, col.names = T, quote = F)
saveRDS(ge.cpm, file.path(res.out, "All_sample_gene_expression_level_log-normalized_count.rds"))


### >>> 3. Seurat pipeline
pd.col <- c("#607d8b","#795548","#ff5722","#ffc107","#cddc39","#4caf50","#009688",
            "#00bcd4","#2196f3","#3f51b5","#673ab7","#9c27b0","#e91e63","#f44336",
            "#b0bec5","#bcaaa4","#ffab91","#ffe082","#e6ee9c","#a5d6a7","#80cbc4",
            "#80deea","#90caf9","#9fa8da","#b39ddb","#ce93d8","#f48fb1","#ef9a9a",
            "#37474f","#4e342e","#d84315","#ff8f00","#9e9d24","#2e7d32","#00695c",
            "#00838f","#1565c0","#283593","#4527a0","#6a1b9a","#ad1457","#c62828")
sr.yhw <- CreateSeuratObject(counts = ge.count[, -1:-5], meta.data = cell.meta, project = "Aneuploid")
sr.yhw <- NormalizeData(sr.yhw) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.yhw))
sr.yhw <- RunPCA(sr.yhw, features = VariableFeatures(object = sr.yhw))
ElbowPlot(sr.yhw, ndims = 30)
dev.off()
sr.yhw <- RunUMAP(sr.yhw, dims = 1:8)
pdf(file.path(res.out, "XuYanWen_Aneuploid_clustering_umap_plot_before_annotation.pdf"), height = 4.5, width = 16)
DimPlot(sr.yhw, reduction = "umap", group.by = c("Stage", "Embryo", "Karyotype"), cols = pd.col, pt.size = 1)
dev.off()
markers <- c("SOX2", "POU5F1", "NANOG", "KLF4", # EPI
             "GATA4", "GATA6", "PDGFRA", "SOX17", # PE
             "CDX2", "GATA3", "KRT8", "KRT18", # TE
             "CTNNB1", "ITGA6", "LRP5", "TP63", "TEAD4", "ELF5", "FGFR2", "FZD5", # CT
             "FN1", "HLA-G", "MMP2", "ITGA5", "CD9", "ITGA1", # EVT
             "PSG1", "HSD3B1", "CYP19A1", "SDC1", "INHA", "ERVW-1", "ERVV-1", "CGA", "CGB5" # ST
             )
pdf(file.path(res.out, "XuYanWen_Aneuploid_clustering_umap_plot_with_markers_expr.pdf"), height = 10, width = 10)
FeatureDimPlot(srt = sr.yhw, features = c(markers[1:12]), reduction = "UMAP", theme_use = "theme_blank")
dev.off()


### >>> 4. Remove chr16 genes
chr16.gene <- read.table("analysis/scrna/metadata/chr16.gene.list")
sr.yhw.no16 <- sr.yhw
sr.yhw.no16 <- NormalizeData(sr.yhw.no16) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000 + length(intersect(chr16.gene$V1, VariableFeatures(object = sr.yhw)))) %>%
  ScaleData(features = rownames(sr.yhw.no16))
sr.yhw.no16 <- RunPCA(sr.yhw.no16, features = setdiff(VariableFeatures(object = sr.yhw.no16), chr16.gene$V1))
ElbowPlot(sr.yhw.no16, ndims = 30)
dev.off()
sr.yhw.no16 <- RunUMAP(sr.yhw.no16, dims = 1:8)
sr.yhw.no16 <- FindNeighbors(sr.yhw.no16, dims = 1:8) %>%
  FindClusters(resolution = seq(0, 3, 0.5))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_clustering_umap_plot_before_annotation.pdf"), height = 9, width = 11)
DimPlot(sr.yhw.no16, reduction = "umap", group.by = c("Stage", "Embryo", "Karyotype", "RNA_snn_res.2"), 
        cols = pd.col, pt.size = 1)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_clustering_umap_plot_with_markers_expr.pdf"), height = 10, width = 10)
FeatureDimPlot(srt = sr.yhw.no16, features = c(markers), reduction = "UMAP", theme_use = "theme_blank")
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_marker_gene_Dot_plot_before_annotation.pdf"), height = 6, width = 15)
DotPlot(sr.yhw.no16, features = unique(markers), cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2", assay = "RNA",
        scale = T, scale.min = 0, scale.max = 100) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(markers)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_clustering_umap_plot_with_markers_expr.pdf"), height = 20, width = 20)
FeatureDimPlot(srt = sr.yhw.no16, features = c(markers), reduction = "UMAP", theme_use = "theme_blank")
dev.off()


### >>> 5. Remove chr22 genes
chr22.gene <- read.table("analysis/scrna/metadata/chr22.gene.list")
sr.yhw.no22 <- sr.yhw
sr.yhw.no22 <- NormalizeData(sr.yhw.no22) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000 + length(intersect(c(chr16.gene$V1, chr22.gene$V1), 
                                                                                     VariableFeatures(object = sr.yhw)))) %>%
  ScaleData(features = rownames(sr.yhw.no22))
sr.yhw.no22 <- RunPCA(sr.yhw.no22, features = setdiff(VariableFeatures(object = sr.yhw.no22), c(chr16.gene$V1, chr22.gene$V1)))
ElbowPlot(sr.yhw.no22, ndims = 30)
dev.off()
sr.yhw.no22 <- RunUMAP(sr.yhw.no22, dims = 1:7)
sr.yhw.no22 <- FindNeighbors(sr.yhw.no22, dims = 1:7) %>%
  FindClusters(resolution = seq(0, 3, 0.5))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_clustering_umap_plot_before_annotation.pdf"), height = 9, width = 11)
DimPlot(sr.yhw.no22, reduction = "umap", group.by = c("Stage", "Embryo", "Karyotype", "RNA_snn_res.2"), 
        cols = pd.col, pt.size = 1, label = T)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_clustering_umap_plot_with_markers_expr.pdf"), height = 20, width = 20)
FeatureDimPlot(srt = sr.yhw.no22, features = c(markers), reduction = "UMAP", theme_use = "theme_blank")
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_gene_Dot_plot_before_annotation.pdf"), height = 6, width = 15)
DotPlot(sr.yhw.no22, features = unique(markers), cols = c("#eaeaea", "#fc0330"),
        col.min = 0, col.max = 2, group.by = "RNA_snn_res.2", assay = "RNA",
        scale = T, scale.min = 0, scale.max = 100) +
  RotatedAxis() +
  scale_x_discrete(labels = unique(markers)) +
  labs(x = "Marker Gene", y = "Cluster") +
  theme_dp
dev.off()


### >>> 6. Rename cell types
sr.yhw.no22@meta.data <- sr.yhw.no22@meta.data %>% 
  mutate(CellType = case_when(RNA_snn_res.2 %in% c(1,4,13) ~ "EPI",
                              RNA_snn_res.2 %in% c(12) ~ "Hypoblast",
                              RNA_snn_res.2 %in% c(0,3,9) ~ "CTB",
                              RNA_snn_res.2 %in% c(5,8,10,11) ~ "pre-STB",
                              RNA_snn_res.2 %in% c(2,6,7) ~ "STB"))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_clustering_umap_plot_after_annotation.pdf"), height = 9, width = 11)
DimPlot(sr.yhw.no22, reduction = "umap", group.by = c("Stage", "Embryo", "Karyotype", "CellType"), cols = pd.col, pt.size = 1, label = T)
dev.off()
sr.yhw.no22$Stage <- factor(sr.yhw.no22$Stage, levels = c("D8", "D10"))
sr.yhw.no22$Karyotype <- factor(sr.yhw.no22$Karyotype, levels = c("Euploidy", "M16", "M16?", "M22", "T16"))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_clustering_umap_plot_before_annotation_modified.pdf"), height = 12, width = 12)
CellDimPlot(
  srt = sr.yhw.no22, group.by = c("Stage", "Embryo", "Karyotype", "CellType"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 1.5
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_clustering_umap_plot_before_annotation_modified2.pdf"), height = 12, width = 12)
CellDimPlot(
  srt = sr.yhw.no22, group.by = c("Stage", "Embryo", "Karyotype", "CellType"),
  reduction = "UMAP", theme_use = "theme_blank", pt.size = 1.5, palcolor = pd.col
)
dev.off()
# Highlight cells
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_clustering_umap_plot_before_annotation_modified_highlight.pdf"), height = 5, width = 5)
CellDimPlot(sr.yhw.no22, pt.size = 1, sizes.highlight = 3, cols.highlight = "black", pt.alpha = 0.5,
            group.by = "CellType", reduction = "UMAP", 
            cells.highlight = colnames(sr.yhw.no22)[sr.yhw.no22$CellID == "D8_A095_53"]
)
dev.off()


### >>> 7. Cell cycle analysis
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sr.yhw.no22 <- CellCycleScoring(sr.yhw.no22, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# plot marker gene expression
ht <- GroupHeatmap(
  srt = sr.yhw.no22,
  features = c(
    "SOX2", "POU5F1", "NANOG", # EPI
    "GATA4", "GATA6", "PDGFRA", "SOX17", # PE
    "CDX2", "GATA3", "KRT8", "KRT18", # TE
    "CTNNB1", "ITGA6", "TEAD4", "FGFR2", "FZD5", # CT
    "HSD3B1", "CYP19A1", "SDC1", "ERVW-1", "ERVV-1", "CGA", "CGB5" # ST
  ),
  group.by = c("Stage", "Karyotype", "CellType"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_gene_Dot_plot_after_annotation.pdf"), 
    height = 8, width = 10)
print(ht$plot)
dev.off()
# plot cell state (all cells)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_count_stacked_barplot_to_show_cell_cycle_after_annotation.pdf"), 
    height = 6, width = 7)
CellStatPlot(sr.yhw.no22, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "trend")
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_ratio_stacked_barplot_to_show_cell_cycle_after_annotation.pdf"), 
    height = 6, width = 7)
CellStatPlot(sr.yhw.no22, stat.by = "Phase", group.by = "CellType", stat_type = "percent", plot_type = "trend")
dev.off()
# plot cell state (split by Karyotype)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_count_stacked_barplot_split_by_karyotype_to_show_cell_cycle_after_annotation.pdf"), 
    height = 8, width = 10)
CellStatPlot(sr.yhw.no22, stat.by = "Phase", group.by = "CellType", stat_type = "count", plot_type = "trend", split.by = "Karyotype")
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_ratio_stacked_barplot_split_by_karyotype_to_show_cell_cycle_after_annotation.pdf"), 
    height = 8, width = 10)
CellStatPlot(sr.yhw.no22, stat.by = "Phase", group.by = "CellType", stat_type = "percent", plot_type = "trend", split.by = "Karyotype")
dev.off()
meta <- data.frame(g1 = c("CTB_M16", "pre-STB_M16", "pre-STB_T16", "CTB_T16", "EPI_T16", "EPI_M22", "EPI_M16",
                          "STB_T16", "STB_M22", "pre-STB_M22", "CTB_M22", "Hypoblast_M22", "Hypoblast_T16"),
                   g2 = c("CTB_Euploidy", "pre-STB_Euploidy", "pre-STB_Euploidy", "CTB_Euploidy", "EPI_Euploidy", "EPI_Euploidy", "EPI_Euploidy",
                          "STB_Euploidy", "STB_Euploidy", "pre-STB_Euploidy", "CTB_Euploidy", "Hypoblast_Euploidy", "Hypoblast_Euploidy"),
                   chi_square_test_statistic = NA,
                   chi_square_test_parameter = NA,
                   chi_square_test_p.value = NA,
                   fisher_test_estimate = NA,
                   fisher_test_p.value = NA)
cc.sig.g2m <- CellCycle.sig(sr.obj = sr.yhw.no22, group.by = "MergedCellType", state = "G2M", meta = meta, outdir = res.out)
cc.sig.g1 <- CellCycle.sig(sr.obj = sr.yhw.no22, group.by = "MergedCellType", state = "G1", meta = meta, outdir = res.out)
cc.sig.s1 <- CellCycle.sig(sr.obj = sr.yhw.no22, group.by = "MergedCellType", state = "S", meta = meta, outdir = res.out)
# plot cell component
sr.yhw.no22$MergedGroup <- paste(sr.yhw.no22$Karyotype, sr.yhw.no22$Stage, sep = "_")
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_count_stacked_barplot_to_show_cell_component_after_annotation.pdf"), 
    height = 8, width = 10)
CellStatPlot(sr.yhw.no22, stat.by = "CellType", group.by = "MergedGroup", stat_type = "count", plot_type = "trend")
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_ratio_stacked_barplot_to_show_cell_component_after_annotation.pdf"), 
    height = 8, width = 10)
CellStatPlot(sr.yhw.no22, stat.by = "CellType", group.by = "MergedGroup", stat_type = "percent", plot_type = "trend")
dev.off()
# plot cell component (filterd version)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_count_stacked_barplot_to_show_cell_component_after_annotation_filtered.pdf"), 
    height = 8, width = 10)
CellStatPlot(subset(sr.yhw.no22, MergedGroup != "M16?_D8"), stat.by = "CellType", 
             group.by = "MergedGroup", stat_type = "count", plot_type = "trend")
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_ratio_stacked_barplot_to_show_cell_component_after_annotation_filtered.pdf"), 
    height = 8, width = 10)
CellStatPlot(subset(sr.yhw.no22, MergedGroup != "M16?_D8"), stat.by = "CellType", 
             group.by = "MergedGroup", stat_type = "percent", plot_type = "trend")
dev.off()


### >>> 8. Plot cell type specific genes
spe.gene <- read.table("R/Table/TrBgenes_xiang.txt", sep = "\t", header = F) %>% 
  filter(V1 %in% rownames(sr.yhw.no22))
table(spe.gene$V2)
library(SCP)
sr.yhw.no22 <- AnnotateFeatures(sr.yhw.no22, species = "Homo_sapiens", db = c("TF", "CSPA"))
sr.yhw.no22$MergedCellType <- paste(sr.yhw.no22$CellType, sr.yhw.no22$Karyotype, sep = "_")
#sr.yhw.no22$MergedCellType <- factor(sr.yhw.no22$MergedCellType, levels = sort(unique(sr.yhw.no22$MergedCellType)))
ht <- FeatureHeatmap(
  srt = sr.yhw.no22, group.by = "CellType", features = spe.gene$V1, feature_split = spe.gene$V2, split.by = "Karyotype",
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 8
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_cell_type_specific_gene_expression_D10_violin_after_annotation.pdf"), 
    height = 10, width = 10)
print(ht$plot)
dev.off()
# calculate gene score
for (i in unique(spe.gene$V2)) {
  sr.yhw.no22@meta.data[, i] <- colMeans(GetAssayData(sr.yhw.no22, slot = "data")[subset(spe.gene, V2 == i)$V1, ])
}
tmp <- subset(sr.yhw.no22, Stage == "D8")
tmp$MergedCellType <- paste(tmp$CellType, tmp$Karyotype, sep = "_")
tmp$MergedCellType <- factor(tmp$MergedCellType, levels = sort(unique(tmp$MergedCellType)))
tmp$Karyotype <- factor(tmp$Karyotype, levels = c("Euploidy", "M16", "M16?", "M22", "T16"))
pd.col.fix <- pd.col[1:5]
names(pd.col.fix) <- c("T16", "Euploidy", "M22", "M16?", "M16")
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_cell_type_specific_gene_expression_D8_violin_after_annotation.pdf"), 
    height = 10, width = 15)
FeatureStatPlot(tmp, stat.by = unique(spe.gene$V2), split.by = "Karyotype", palcolor = pd.col.fix, bg_palcolor = pd.col.fix,
                group.by = "MergedCellType", bg.by = "Karyotype", add_box = TRUE, stack = T, 
                comparisons = list(c("CTB_M16", "CTB_Euploidy"), 
                                   c("STB_T16", "STB_Euploidy"), 
                                   c("pre-STB_Euploidy", "pre-STB_M16"),
                                   c("CTB_Euploidy", "CTB_M22"),
                                   c("pre-STB_Euploidy", "pre-STB_M22"),
                                   c("STB_Euploidy", "STB_M22")))
dev.off()
tmp <- subset(sr.yhw.no22, Stage == "D10")
tmp$MergedCellType <- paste(tmp$CellType, tmp$Karyotype, sep = "_")
tmp$MergedCellType <- factor(tmp$MergedCellType, levels = sort(unique(tmp$MergedCellType)))
tmp$Karyotype <- factor(tmp$Karyotype, levels = c("Euploidy", "M16", "M16?", "M22", "T16"))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_cell_type_specific_gene_expression_D10_violin_after_annotation.pdf"), 
    height = 10, width = 10)
FeatureStatPlot(tmp, stat.by = unique(spe.gene$V2), split.by = "Karyotype", palcolor = pd.col.fix, bg_palcolor = pd.col.fix,
                group.by = "MergedCellType", bg.by = "Karyotype", add_box = TRUE, stack = T, 
                comparisons = list(c("STB_T16", "STB_Euploidy")))
dev.off()


### >>> 9. Velocyto
DefaultAssay(sr.yhw.no22) <- "RNA"
SeuratToVelocyto(sr.ob = sr.yhw.no22, reduction = "umap", out.dir = file.path(res.out, "velocyto"), 
                 sample.name = "sr_NO22", str.in.barcode = NULL)
tmp@meta.data <- tmp@meta.data %>% 
  mutate(Highlight = case_when(MergedGroup == "Euploidy_D10" ~ "Yes",
                               MergedGroup == "T16_D10" ~ "Yes",
                               MergedGroup != "Euploidy_D10" & MergedGroup != "T16_D10" ~ "Others"))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_cell_type_highlighted_for_velocyto_T16.pdf"), 
    height = 5, width = 5)
CellDimPlot(tmp, group.by = "MergedGroup", 
            cells.highlight = colnames(tmp)[tmp$Highlight == "Yes"], 
            sizes.highlight = 1, reduction = "UMAP",
            theme_use = "theme_blank")
dev.off()
tmp@meta.data <- tmp@meta.data %>% 
  mutate(Highlight = case_when(MergedGroup == "Euploidy_D8" ~ "Yes",
                               MergedGroup == "M16_D8" ~ "Yes",
                               MergedGroup != "Euploidy_D8" & MergedGroup != "M16_D8" ~ "Others"))
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_cell_type_highlighted_for_velocyto_M16.pdf"), 
    height = 5, width = 5)
CellDimPlot(tmp, group.by = "MergedGroup", 
            cells.highlight = colnames(tmp)[tmp$Highlight == "Yes"], 
            sizes.highlight = 1, reduction = "UMAP",
            theme_use = "theme_blank")
dev.off()



# ==========================================
# 4th part: Differential expression analysis ----
# ==========================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/DEGs")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Run edgeR
hs.tfs <- read.table("/home/yhw/bioinfo/project-xuyanwen/aneuploid/R/Table/Homo_sapiens_TF.txt", header = T, sep = "\t")
# T16 vs Euploidy
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType)
deg.t16 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "T16"), sample.n = NULL)
  deg.t16[[paste0("T16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_T16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("T16/", i, "_T16_vs_Euploidy")))
}
for (i in names(deg.t16)) {
  pdf(file.path(res.out, paste0("T16/", i, "_T16_vs_Euploidy_volcano.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.t16[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.01, lfc = 2, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
  pdf(file.path(res.out, paste0("T16/", i, "_T16_vs_Euploidy_volcano_TFS.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.t16[[i]]$all, geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 3, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
for (i in names(deg.t16)) {
  p1 <- VisDEG.volcano.v2(pd = deg.t16[[i]]$all, geneset =c(hs.tfs$Symbol, "CDH1"), 
                          p.col = "PValue", q.col = "FDR", lfc.col = "logFC", sig = 0.05, lfc = log2(1.25),
                          curve.plot = F, nosig.remove = F,
                          pt.size.by = NULL, pt.size = 2, pt.shape = 19, pt.alpha = 1,
                          l.width = 0.8, top.n = NULL, label.gene = c("CDH1"), gene.size = 3,
                          high.col = "#D9212A", low.col = "#045EC3", nosig.col = "#B2B2B2") + ggtitle(paste0(i, ": CDH1"))
  p2 <- VisDEG.volcano.v2(pd = deg.t16[[i]]$all, geneset =c(hs.tfs$Symbol, "CDH1"), 
                          p.col = "PValue", q.col = "FDR", lfc.col = "logFC", sig = 0.05, lfc = log2(1.25),
                          curve.plot = F, nosig.remove = F, 
                          pt.size.by = NULL, pt.size = 2, pt.shape = 19, pt.alpha = 1,
                          l.width = 0.8, top.n = 10, label.gene = NULL, gene.size = 3,
                          high.col = "#D9212A", low.col = "#045EC3", nosig.col = "#B2B2B2") + ggtitle(paste0(i, ": Top10 TFs"))
  pdf(file.path(res.out, paste0("T16/Volcano_plot_to_show_", i, "_with_CDH1.pdf")), height = 5, width = 13)
  print(p1 + p2)
  dev.off()
}
# M22 vs Euploidy
deg.m22 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M22"), sample.n = NULL)
  deg.m22[[paste0("M22_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_M22", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("M22/", i, "_M22_vs_Euploidy")))
}
for (i in names(deg.m22)) {
  pdf(file.path(res.out, paste0("M22/", i, "_M22_vs_Euploidy_volcano.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.m22[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.01, lfc = 2, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
  pdf(file.path(res.out, paste0("M22/", i, "_M22_vs_Euploidy_volcano_TFS.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.m22[[i]]$all, geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 3, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M16? vs Euploidy
deg.m16x <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp.sr$Karyotype <- gsub("M16.", "M16x", tmp.sr$Karyotype)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16x"), sample.n = NULL)
  deg.m16x[[paste0("M16x_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                sample.n = NULL, group.by = "CellType",
                                                g1 = "1_Euploidy", g2 = "2_M16x", lfc = 1, sig = 0.05,
                                                res.out = file.path(res.out, paste0("M16?/", i, "_M16x_vs_Euploidy")))
}
for (i in names(deg.m16x)) {
  pdf(file.path(res.out, paste0("M16?/", i, "_M16x_vs_Euploidy_volcano.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.m16x[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.01, lfc = 2, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
  pdf(file.path(res.out, paste0("M16?/", i, "_M16x_vs_Euploidy_volcano_TFS.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.m16x[[i]]$all, geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 3, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M16 vs Euploidy
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType)
deg.m16 <- list()
for (i in c("CTB", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16"), sample.n = NULL)
  deg.m16[[paste0("M16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_M16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("M16/", i, "_M16_vs_Euploidy")))
}
for (i in names(deg.m16)) {
  pdf(file.path(res.out, paste0("M16/", i, "_M16_vs_Euploidy_volcano.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.m16[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.01, lfc = 2, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
  pdf(file.path(res.out, paste0("M16/", i, "_M16_vs_Euploidy_volcano_TFS.pdf")), height = 6, width = 8)
  print(VisDEG.volcano(deg.data = deg.m16[[i]]$all, geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 3, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# multiple volcano (all, 0.001, 0.005, 0.01, 0.05, 0.1)
pd <- c(deg.m16, deg.m22, deg.t16.d8, deg.t16.d10)
names(pd)[7:14] <- c("T16_D8_CTB", "T16_D8_EPI", "T16_D8_Hypoblast", "T16_D8_pre-STB", "T16_D8_STB",
                     "T16_D10_EPI", "T16_D10_pre-STB", "T16_D10_STB")
for (i in names(pd)) {
  pd[[i]] <- subset(pd[[i]]$all, abs(logFC) >= log2(1.5) & PValue <= 0.05)
  pd[[i]] <- pd[[i]][pd[[i]][, 7] >= 0.005, ]
}
pdf(file.path(res.out, "All_DEG_vs_Ctrl_volcano_all_TFs_0.1.pdf"), height = 6, width = 6)
p1 <- VisDEG.volcano.multi(deg.list = pd, geneset = hs.tfs$Symbol, lfc = log2(1.5), sig = 0.05,
                           top.n = 1, pt.size = 0.25, jitter.width = 0.35, gene.size = 4, 
                           label.gene = NULL,
                           high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9")
table(p1$data$GroupToVolcano, p1$data$GroupToChange)
write.csv(p1$data, file.path(res.out, "All_DEG_vs_Ctrl_volcano_all_TFs_suppl.csv"), row.names = F, quote = F)
p2 <- table(p1$data$GroupToVolcano, p1$data$GroupToChange) %>% 
  as.data.frame() %>% 
  spread(key = "Var1", value = "Freq")
tbl <- tableGrob(p2, rows = NULL)
grid.arrange(p1$plot,
             tbl,
             nrow = 2,
             as.table = TRUE,
             heights = c(4, 1))
dev.off()
pdf(file.path(res.out, "All_DEG_vs_Ctrl_volcano_all_DEGs_0.1.pdf"), height = 6, width = 6)
p1 <- VisDEG.volcano.multi(deg.list = pd, geneset = NULL, lfc = log2(1.5), sig = 0.05,
                           top.n = 1, pt.size = 0.5, jitter.width = 0.35, gene.size = 4, 
                           label.gene = NULL,
                           high.col = "#33A02C", low.col = "#6A3D9A", nosig.col = "#D9D9D9")
table(p1$data$GroupToVolcano, p1$data$GroupToChange)
write.csv(p1$data, file.path(res.out, "All_DEG_vs_Ctrl_volcano_all_DEGs_suppl.csv"), row.names = F, quote = F)
p2 <- table(p1$data$GroupToVolcano, p1$data$GroupToChange) %>% 
  as.data.frame() %>% 
  spread(key = "Var1", value = "Freq")
library(grid)
library(gridExtra)
tbl <- tableGrob(p2, rows = NULL)
grid.arrange(p1$plot,
             tbl,
             nrow = 2,
             as.table = TRUE,
             heights = c(4, 1))
dev.off()
# overlapped DEGs
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  j <- grep(paste0("_", i, "$"), names(pd))
  pd.list <- pd[j]
  for (k in names(pd.list)) {
    pd.list[[k]] <- pd.list[[k]]$SYMBOL
  }
  PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = i, 
                 file.name = paste0("Overlapped_DEGs_of_", i), 
                 res.out = file.path(res.out, "Overlapped_0.1"),
                 pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
}

  
### >>> 3. GO analysis
# T16 vs Euploidy
go.t16 <- list()
for (i in names(deg.t16)) {
  go.t16[[paste0(i, ".up")]] <- Pipe.GO(species = "human", 
                                        genelist = subset(deg.t16[[i]]$sig, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                                        basename = paste0(i, "_T16_vs_Euploidy_up"), genetype = "SYMBOL",
                                        res.out = file.path(res.out, paste0("T16/GO/", i)))
  go.t16[[paste0(i, ".down")]] <- Pipe.GO(species = "human", 
                                          genelist = subset(deg.t16[[i]]$sig, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                                          basename = paste0(i, "_T16_vs_Euploidy_down"), genetype = "SYMBOL",
                                          res.out = file.path(res.out, paste0("T16/GO/", i)))
}
for (i in names(deg.t16)) {
  go.t16[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                              genelist = intersect(chr16.gene$V1,
                                                                   subset(deg.t16[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                              basename = paste0(i, "_T16_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                              res.out = file.path(res.out, paste0("T16/GO_chr16/", i)))
  go.t16[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                genelist = intersect(chr16.gene$V1,
                                                                     subset(deg.t16[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                basename = paste0(i, "_T16_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                res.out = file.path(res.out, paste0("T16/GO_chr16/", i)))
}
gsea.t16 <- list()
for (i in names(deg.t16)) {
  gsea.t16[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.t16[[i]]$all, deg.type = "edger",
                                                         lfc = 1, sig = 0.05,
                                                         reversed = FALSE, species = "human",
                                                         basename = paste0(i, "_lfc1_pvalue0.05"),
                                                         genetype = "SYMBOL", gene.col = "SYMBOL",
                                                         outdir = file.path(res.out, paste0("GSEA/T16_lfc1_pvalue0.05/", i)))
  gsea.t16[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.t16[[i]]$all, deg.type = "edger",
                                                          lfc = 2, sig = 0.05,
                                                          reversed = FALSE, species = "human",
                                                          basename = paste0(i, "_lfc2_pvalue0.005"),
                                                          genetype = "SYMBOL", gene.col = "SYMBOL",
                                                          outdir = file.path(res.out, paste0("GSEA/T16_lfc2_pvalue0.005/", i)))
  gsea.t16[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.t16[[i]]$all, deg.type = "edger", 
                                             lfc = 0, sig = 1, 
                                             reversed = FALSE, species = "human", 
                                             basename = paste0(i, "_all"), 
                                             genetype = "SYMBOL", gene.col = "SYMBOL", 
                                             outdir = file.path(res.out, paste0("GSEA/T16_all/", i)))
}
# M22 vs Euploidy
go.m22 <- list()
for (i in names(deg.m22)) {
  go.m22[[paste0(i, ".up")]] <- Pipe.GO(species = "human", 
                                        genelist = subset(deg.m22[[i]]$sig, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                                        basename = paste0(i, "_M22_vs_Euploidy_up"), genetype = "SYMBOL",
                                        res.out = file.path(res.out, paste0("M22/GO/", i)))
  go.m22[[paste0(i, ".down")]] <- Pipe.GO(species = "human", 
                                          genelist = subset(deg.m22[[i]]$sig, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                                          basename = paste0(i, "_M22_vs_Euploidy_down"), genetype = "SYMBOL",
                                          res.out = file.path(res.out, paste0("M22/GO/", i)))
}
for (i in names(deg.m22)) {
  go.m22[[paste0(i, ".chr22.up")]] <- Pipe.GO(species = "human", 
                                              genelist = intersect(chr22.gene$V1,
                                                                   subset(deg.m22[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                              basename = paste0(i, "_M22_vs_Euploidy_chr22_up"), genetype = "SYMBOL",
                                              res.out = file.path(res.out, paste0("M22/GO_chr22/", i)))
  go.m22[[paste0(i, ".chr22.down")]] <- Pipe.GO(species = "human", 
                                                genelist = intersect(chr22.gene$V1,
                                                                     subset(deg.m22[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                basename = paste0(i, "_M22_vs_Euploidy_chr22_down"), genetype = "SYMBOL",
                                                res.out = file.path(res.out, paste0("M22/GO_chr22/", i)))
}
gsea.m22 <- list()
for (i in names(deg.m22)) {
  gsea.m22[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.m22[[i]]$all, deg.type = "edger",
                                                         lfc = 1, sig = 0.05,
                                                         reversed = FALSE, species = "human",
                                                         basename = paste0(i, "_lfc1_pvalue0.05"),
                                                         genetype = "SYMBOL", gene.col = "SYMBOL",
                                                         outdir = file.path(res.out, paste0("GSEA/M22_lfc1_pvalue0.05/", i)))
  gsea.m22[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.m22[[i]]$all, deg.type = "edger",
                                                          lfc = 2, sig = 0.05,
                                                          reversed = FALSE, species = "human",
                                                          basename = paste0(i, "_lfc2_pvalue0.005"),
                                                          genetype = "SYMBOL", gene.col = "SYMBOL",
                                                          outdir = file.path(res.out, paste0("GSEA/M22_lfc2_pvalue0.005/", i)))
  gsea.m22[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.m22[[i]]$all, deg.type = "edger", 
                                             lfc = 0, sig = 1, 
                                             reversed = FALSE, species = "human", 
                                             basename = paste0(i, "_all"), 
                                             genetype = "SYMBOL", gene.col = "SYMBOL", 
                                             outdir = file.path(res.out, paste0("GSEA/M22_all/", i)))
}
# M16? vs Euploidy
go.m16x <- list()
for (i in names(deg.m16x)) {
  go.m16x[[paste0(i, ".up")]] <- Pipe.GO(species = "human", 
                                         genelist = subset(deg.m16x[[i]]$sig, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                                         basename = paste0(i, "_M16x_vs_Euploidy_up"), genetype = "SYMBOL",
                                         res.out = file.path(res.out, paste0("M16?/GO/", i)))
  go.m16x[[paste0(i, ".down")]] <- Pipe.GO(species = "human", 
                                           genelist = subset(deg.m16x[[i]]$sig, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                                           basename = paste0(i, "_M16x_vs_Euploidy_down"), genetype = "SYMBOL",
                                           res.out = file.path(res.out, paste0("M16?/GO/", i)))
}
for (i in names(deg.m16x)) {
  go.m16x[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                               genelist = intersect(chr16.gene$V1,
                                                                    subset(deg.m16x[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                               basename = paste0(i, "_M16x_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                               res.out = file.path(res.out, paste0("M16?/GO_chr16/", i)))
  go.m16x[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                 genelist = intersect(chr16.gene$V1,
                                                                      subset(deg.m16x[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                 basename = paste0(i, "_M16x_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("M16?/GO_chr16/", i)))
}
gsea.m16x <- list()
for (i in names(deg.m16x)) {
  gsea.m16x[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.m16x[[i]]$all, deg.type = "edger",
                                                          lfc = 1, sig = 0.05,
                                                          reversed = FALSE, species = "human",
                                                          basename = paste0(i, "_lfc1_pvalue0.05"),
                                                          genetype = "SYMBOL", gene.col = "SYMBOL",
                                                          outdir = file.path(res.out, paste0("GSEA/M16x_lfc1_pvalue0.05/", i)))
  gsea.m16x[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.m16x[[i]]$all, deg.type = "edger",
                                                           lfc = 2, sig = 0.05,
                                                           reversed = FALSE, species = "human",
                                                           basename = paste0(i, "_lfc2_pvalue0.005"),
                                                           genetype = "SYMBOL", gene.col = "SYMBOL",
                                                           outdir = file.path(res.out, paste0("GSEA/M16x_lfc2_pvalue0.005/", i)))
  gsea.m16x[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.m16x[[i]]$all, deg.type = "edger", 
                                              lfc = 0, sig = 1, 
                                              reversed = FALSE, species = "human", 
                                              basename = paste0(i, "_all"), 
                                              genetype = "SYMBOL", gene.col = "SYMBOL", 
                                              outdir = file.path(res.out, paste0("GSEA/M16x_all/", i)))
}
# M16 vs Euploidy
go.m16 <- list()
for (i in names(deg.m16)) {
  go.m16[[paste0(i, ".up")]] <- Pipe.GO(species = "human", 
                                        genelist = subset(deg.m16[[i]]$sig, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                                        basename = paste0(i, "_M16_vs_Euploidy_up"), genetype = "SYMBOL",
                                        res.out = file.path(res.out, paste0("M16/GO/", i)))
  go.m16[[paste0(i, ".down")]] <- Pipe.GO(species = "human", 
                                          genelist = subset(deg.m16[[i]]$sig, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                                          basename = paste0(i, "_M16_vs_Euploidy_down"), genetype = "SYMBOL",
                                          res.out = file.path(res.out, paste0("M16/GO/", i)))
}
for (i in names(deg.m16)) {
  go.m16[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                              genelist = intersect(chr16.gene$V1,
                                                                   subset(deg.m16[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                              basename = paste0(i, "_M16_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                              res.out = file.path(res.out, paste0("M16/GO_chr16/", i)))
  go.m16[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                genelist = intersect(chr16.gene$V1,
                                                                     subset(deg.m16[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                basename = paste0(i, "_M16_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                res.out = file.path(res.out, paste0("M16/GO_chr16/", i)))
}
gsea.m16 <- list()
for (i in names(deg.m16)) {
  gsea.m16[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.m16[[i]]$all, deg.type = "edger",
                                                         lfc = 1, sig = 0.05,
                                                         reversed = FALSE, species = "human",
                                                         basename = paste0(i, "_lfc1_pvalue0.05"),
                                                         genetype = "SYMBOL", gene.col = "SYMBOL",
                                                         outdir = file.path(res.out, paste0("GSEA/M16_lfc1_pvalue0.05/", i)))
  gsea.m16[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.m16[[i]]$all, deg.type = "edger",
                                                          lfc = 2, sig = 0.05,
                                                          reversed = FALSE, species = "human",
                                                          basename = paste0(i, "_lfc2_pvalue0.005"),
                                                          genetype = "SYMBOL", gene.col = "SYMBOL",
                                                          outdir = file.path(res.out, paste0("GSEA/M16_lfc2_pvalue0.005/", i)))
  gsea.m16[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.m16[[i]]$all, deg.type = "edger", 
                                             lfc = 0, sig = 1, 
                                             reversed = FALSE, species = "human", 
                                             basename = paste0(i, "_all"), 
                                             genetype = "SYMBOL", gene.col = "SYMBOL", 
                                             outdir = file.path(res.out, paste0("GSEA/M16_all/", i)))
}
# plot GO
for (j in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  go.data <- list(CTB.M16.down = go.m16.d8[[paste0("M16_CTB.down.", j)]],
                  CTB.M22.down = go.m22.d8[[paste0("M22_CTB.down.", j)]],
                  CTB.T16.D8.down = go.t16.d8[[paste0("T16_CTB.down.", j)]],
                  EPI.M22.down = go.m22.d8[[paste0("M22_EPI.down.", j)]],
                  EPI.T16.D10.down = go.t16.d10[[paste0("T16_EPI.down.", j)]],
                  EPI.T16.D8.down = go.t16.d8[[paste0("T16_EPI.down.", j)]],
                  Hypo.M22.down = go.m22.d8[[paste0("M22_Hypoblast.down.", j)]],
                  Hypo.T16.D8.down = go.t16.d8[[paste0("T16_Hypoblast.down.", j)]],
                  preSTB.M16.down = go.m16.d8[[paste0("M16_pre-STB.down.", j)]],
                  preSTB.M22.down = go.m22.d8[[paste0("M22_pre-STB.down.", j)]],
                  preSTB.T16.D10.down = go.t16.d10[[paste0("T16_pre-STB.down.", j)]],
                  preSTB.T16.D8.down = go.t16.d8[[paste0("T16_pre-STB.down.", j)]],
                  STB.T16.D10.down = go.t16.d10[[paste0("T16_STB.down.", j)]],
                  STB.T16.D8.down = go.t16.d8[[paste0("T16_STB.down.", j)]],
                  CTB.M16.up = go.m16.d8[[paste0("M16_CTB.up.", j)]],
                  CTB.M22.up = go.m22.d8[[paste0("M22_CTB.up.", j)]],
                  CTB.T16.D8.up = go.t16.d8[[paste0("T16_CTB.up.", j)]],
                  EPI.M22.up = go.m22.d8[[paste0("M22_EPI.up.", j)]],
                  EPI.T16.D10.up = go.t16.d10[[paste0("T16_EPI.up.", j)]],
                  EPI.T16.D8.up = go.t16.d8[[paste0("T16_EPI.up.", j)]],
                  Hypo.M22.up = go.m22.d8[[paste0("M22_Hypoblast.up.", j)]],
                  Hypo.T16.D8.up = go.t16.d8[[paste0("T16_Hypoblast.up.", j)]],
                  preSTB.M16.up = go.m16.d8[[paste0("M16_pre-STB.up.", j)]],
                  preSTB.M22.up = go.m22.d8[[paste0("M22_pre-STB.up.", j)]],
                  preSTB.T16.D10.up = go.t16.d10[[paste0("T16_pre-STB.up.", j)]],
                  preSTB.T16.D8.up = go.t16.d8[[paste0("T16_pre-STB.up.", j)]],
                  STB.T16.D10.up = go.t16.d10[[paste0("T16_STB.up.", j)]],
                  STB.T16.D8.up = go.t16.d8[[paste0("T16_STB.up.", j)]])
  go.terms <- c("hsa04020", "hsa04148", "hsa04310", "hsa04390", "hsa04510", "hsa04512",
                "hsa04514", "hsa04550", "hsa04916", "hsa04934", "hsa04974", "hsa05032",
                "hsa05034", "hsa05165", "hsa05166", "hsa05202", "hsa05217", "hsa05222",
                "hsa05224", "hsa05225", "hsa05226", "hsa05410", "hsa05412", "hsa05418")
  pd.seq <- names(go.data)
  pd <- GO.replot(go.data = go.data, go.terms = go.terms, 
                  color.by = "GeneRatio", size.by = "pvalue", cols.pal = viridis(12),
                  multi.go = T, multi.cluster = T, seq.x = pd.seq, seq.y = NULL,
                  col.min = 0, col.max = 0.05, pd.title = paste0("GO analysis from DEGs ", j))
  dev.off()
  pdf(file.path(res.out, paste0("Replot_of_GO_analysis_DEG_All_", j, ".pdf")), 
      height = 6, width = 14)
  print(pd$bubble)
  pd$data
  dev.off()
}


### >>> 4. D8 and D10
# T16 vs Euploidy: D8
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType, sr.yhw.no22$Stage)
deg.t16.d8 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "T16"), sample.n = NULL)
  deg.t16.d8[[paste0("T16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_T16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("T16_D8/", i, "_T16_vs_Euploidy")))
}
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.t16.d8)) {
    pdf(file.path(res.out, paste0("T16_D8/", i, "_T16_vs_Euploidy_volcano_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.t16.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.t16.d8[[i]]$all), value = T), low.expr = j,
                         geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.01, lfc = 2, title = i,
                         pt.size = 2, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
    pdf(file.path(res.out, paste0("T16_D8/", i, "_T16_vs_Euploidy_volcano_TFs_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.t16.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.t16.d8[[i]]$all), value = T), low.expr = j,
                         geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.05, lfc = 1, title = i,
                         pt.size = 3, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
  }
}
# T16 vs Euploidy: D10
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType, sr.yhw.no22$Stage)
deg.t16.d10 <- list()
for (i in c("EPI", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i & Stage == "D10")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "T16"), sample.n = NULL)
  deg.t16.d10[[paste0("T16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                  sample.n = NULL, group.by = "CellType",
                                                  g1 = "1_Euploidy", g2 = "2_T16", lfc = 1, sig = 0.05,
                                                  res.out = file.path(res.out, paste0("T16_D10/", i, "_T16_vs_Euploidy")))
}
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.t16.d10)) {
    pdf(file.path(res.out, paste0("T16_D10/", i, "_T16_vs_Euploidy_volcano_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.t16.d10[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.t16.d10[[i]]$all), value = T), low.expr = j,
                         geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.01, lfc = 2, title = i,
                         pt.size = 2, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
    pdf(file.path(res.out, paste0("T16_D10/", i, "_T16_vs_Euploidy_volcano_TFs_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.t16.d10[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.t16.d10[[i]]$all), value = T), low.expr = j,
                         geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.05, lfc = 1, title = i,
                         pt.size = 3, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
  }
}
# M16 vs Euploidy: D8
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType, sr.yhw.no22$Stage)
deg.m16.d8 <- list()
for (i in c("CTB", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16"), sample.n = NULL)
  deg.m16.d8[[paste0("M16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                 sample.n = NULL, group.by = "CellType",
                                                 g1 = "1_Euploidy", g2 = "2_M16", lfc = 1, sig = 0.05,
                                                 res.out = file.path(res.out, paste0("M16_D8/", i, "_M16_vs_Euploidy")))
}
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.m16.d8)) {
    pdf(file.path(res.out, paste0("M16_D8/", i, "_M16_vs_Euploidy_volcano_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.m16.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.m16.d8[[i]]$all), value = T), low.expr = j,
                         geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.01, lfc = 2, title = i,
                         pt.size = 2, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
    pdf(file.path(res.out, paste0("M16_D8/", i, "_M16_vs_Euploidy_volcano_TFs_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.m16.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.m16.d8[[i]]$all), value = T), low.expr = j,
                         geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.05, lfc = 1, title = i,
                         pt.size = 3, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
  }
}
# M16? vs Euploidy: D8
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType, sr.yhw.no22$Stage)
deg.m16x.d8 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i & Stage == "D8")
  tmp.sr$Karyotype <- gsub("M16.", "M16x", tmp.sr$Karyotype)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16x"), sample.n = NULL)
  deg.m16x.d8[[paste0("M16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                 sample.n = NULL, group.by = "CellType",
                                                 g1 = "1_Euploidy", g2 = "2_M16x", lfc = 1, sig = 0.05,
                                                 res.out = file.path(res.out, paste0("M16?_D8/", i, "_M16_vs_Euploidy")))
}
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.m16x.d8)) {
    pdf(file.path(res.out, paste0("M16?_D8/", i, "_M16_vs_Euploidy_volcano_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.m16x.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.m16x.d8[[i]]$all), value = T), low.expr = j,
                         geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.01, lfc = 2, title = i,
                         pt.size = 2, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
    pdf(file.path(res.out, paste0("M16?_D8/", i, "_M16_vs_Euploidy_volcano_TFs_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.m16x.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.m16x.d8[[i]]$all), value = T), low.expr = j,
                         geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.05, lfc = 1, title = i,
                         pt.size = 3, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
  }
}
# M22 vs Euploidy: D8
table(sr.yhw.no22$Karyotype, sr.yhw.no22$CellType, sr.yhw.no22$Stage)
deg.m22.d8 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M22"), sample.n = NULL)
  deg.m22.d8[[paste0("M22_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                 sample.n = NULL, group.by = "CellType",
                                                 g1 = "1_Euploidy", g2 = "2_M22", lfc = 1, sig = 0.05,
                                                 res.out = file.path(res.out, paste0("M22_D8/", i, "_M22_vs_Euploidy")))
}
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.m22.d8)) {
    pdf(file.path(res.out, paste0("M22_D8/", i, "_M22_vs_Euploidy_volcano_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.m22.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.m22.d8[[i]]$all), value = T), low.expr = j,
                         geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.01, lfc = 2, title = i,
                         pt.size = 2, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
    pdf(file.path(res.out, paste0("M22_D8/", i, "_M22_vs_Euploidy_volcano_TFs_", j, ".pdf")), height = 6, width = 9)
    print(VisDEG.volcano(deg.data = deg.m22.d8[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.m22.d8[[i]]$all), value = T), low.expr = j,
                         geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.05, lfc = 1, title = i,
                         pt.size = 3, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
  }
}


### >>> 5. GO analysis
# T16: D8
go.t16.d8 <- list()
for (j in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.t16.d8)) {
    pd <- subset(deg.t16.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.t16.d8[[paste0(i, ".up.", j)]] <- Pipe.GO(species = "human", 
                                                 genelist = pd$SYMBOL,
                                                 basename = paste0(i, "_T16_vs_Euploidy_up_", j), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("T16_D8/GO/", i, "_", j)))
    pd <- subset(deg.t16.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.t16.d8[[paste0(i, ".down.", j)]] <- Pipe.GO(species = "human", 
                                                   genelist = pd$SYMBOL,
                                                   basename = paste0(i, "_T16_vs_Euploidy_down_", j), genetype = "SYMBOL",
                                                   res.out = file.path(res.out, paste0("T16_D8/GO/", i, "_", j)))
  }
}
for (i in names(deg.t16.d8)) {
  go.t16.d8[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                                 genelist = intersect(chr16.gene$V1,
                                                                      subset(deg.t16.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                 basename = paste0(i, "_T16_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("T16_D8/GO_chr16/", i)))
  go.t16.d8[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                   genelist = intersect(chr16.gene$V1,
                                                                        subset(deg.t16.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                   basename = paste0(i, "_T16_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                   res.out = file.path(res.out, paste0("T16_D8/GO_chr16/", i)))
}
gsea.t16.d8 <- list()
for (i in names(deg.t16.d8)) {
  gsea.t16.d8[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.t16.d8[[i]]$all, deg.type = "edger",
                                                            lfc = 1, sig = 0.05,
                                                            reversed = FALSE, species = "human",
                                                            basename = paste0(i, "_lfc1_pvalue0.05"),
                                                            genetype = "SYMBOL", gene.col = "SYMBOL",
                                                            outdir = file.path(res.out, paste0("GSEA/T16_D8_lfc1_pvalue0.05/", i)))
  gsea.t16.d8[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.t16.d8[[i]]$all, deg.type = "edger",
                                                             lfc = 2, sig = 0.05,
                                                             reversed = FALSE, species = "human",
                                                             basename = paste0(i, "_lfc2_pvalue0.005"),
                                                             genetype = "SYMBOL", gene.col = "SYMBOL",
                                                             outdir = file.path(res.out, paste0("GSEA/T16_D8_lfc2_pvalue0.005/", i)))
  gsea.t16.d8[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.t16.d8[[i]]$all, deg.type = "edger", 
                                                lfc = 0, sig = 1, 
                                                reversed = FALSE, species = "human", 
                                                basename = paste0(i, "_all"), 
                                                genetype = "SYMBOL", gene.col = "SYMBOL", 
                                                outdir = file.path(res.out, paste0("GSEA/T16_D8_all/", i)))
}
# M16: D8
go.m16.d8 <- list()
for (j in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.m16.d8)) {
    pd <- subset(deg.m16.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.m16.d8[[paste0(i, ".up.", j)]] <- Pipe.GO(species = "human", 
                                                 genelist = pd$SYMBOL,
                                                 basename = paste0(i, "_M16_vs_Euploidy_up_", j), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("M16_D8/GO/", i, "_", j)))
    pd <- subset(deg.m16.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.m16.d8[[paste0(i, ".down.", j)]] <- Pipe.GO(species = "human", 
                                                   genelist = pd$SYMBOL,
                                                   basename = paste0(i, "_M16_vs_Euploidy_down_", j), genetype = "SYMBOL",
                                                   res.out = file.path(res.out, paste0("M16_D8/GO/", i, "_", j)))
  }
}
for (i in names(deg.m16.d8)) {
  go.m16.d8[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                                 genelist = intersect(chr16.gene$V1,
                                                                      subset(deg.m16.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                 basename = paste0(i, "_M16_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("M16_D8/GO_chr16/", i)))
  go.m16.d8[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                   genelist = intersect(chr16.gene$V1,
                                                                        subset(deg.m16.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                   basename = paste0(i, "_M16_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                   res.out = file.path(res.out, paste0("M16_D8/GO_chr16/", i)))
}
gsea.m16.d8 <- list()
for (i in names(deg.m16.d8)) {
  gsea.m16.d8[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.m16.d8[[i]]$all, deg.type = "edger",
                                                            lfc = 1, sig = 0.05,
                                                            reversed = FALSE, species = "human",
                                                            basename = paste0(i, "_lfc1_pvalue0.05"),
                                                            genetype = "SYMBOL", gene.col = "SYMBOL",
                                                            outdir = file.path(res.out, paste0("GSEA/M16_D8_lfc1_pvalue0.05/", i)))
  gsea.m16.d8[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.m16.d8[[i]]$all, deg.type = "edger",
                                                             lfc = 2, sig = 0.05,
                                                             reversed = FALSE, species = "human",
                                                             basename = paste0(i, "_lfc2_pvalue0.005"),
                                                             genetype = "SYMBOL", gene.col = "SYMBOL",
                                                             outdir = file.path(res.out, paste0("GSEA/M16_D8_lfc2_pvalue0.005/", i)))
  gsea.m16.d8[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.m16.d8[[i]]$all, deg.type = "edger", 
                                                lfc = 0, sig = 1, 
                                                reversed = FALSE, species = "human", 
                                                basename = paste0(i, "_all"), 
                                                genetype = "SYMBOL", gene.col = "SYMBOL", 
                                                outdir = file.path(res.out, paste0("GSEA/M16_D8_all/", i)))
}
# M16?: D8
go.m16x.d8 <- list()
for (j in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.m16x.d8)) {
    pd <- subset(deg.m16x.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.m16x.d8[[paste0(i, ".up.", j)]] <- Pipe.GO(species = "human", 
                                                  genelist = pd$SYMBOL,
                                                  basename = paste0(i, "_M16x_vs_Euploidy_up_", j), genetype = "SYMBOL",
                                                  res.out = file.path(res.out, paste0("M16?_D8/GO/", i, "_", j)))
    pd <- subset(deg.m16x.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.m16x.d8[[paste0(i, ".down.", j)]] <- Pipe.GO(species = "human", 
                                                    genelist = pd$SYMBOL,
                                                    basename = paste0(i, "_M16x_vs_Euploidy_down_", j), genetype = "SYMBOL",
                                                    res.out = file.path(res.out, paste0("M16?_D8/GO/", i, "_", j)))
  }
}
for (i in names(deg.m16x.d8)) {
  go.m16x.d8[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                                  genelist = intersect(chr16.gene$V1,
                                                                       subset(deg.m16x.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                  basename = paste0(i, "_M16x_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                                  res.out = file.path(res.out, paste0("M16?_D8/GO_chr16/", i)))
  go.m16x.d8[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                    genelist = intersect(chr16.gene$V1,
                                                                         subset(deg.m16x.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                    basename = paste0(i, "_M16x_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                    res.out = file.path(res.out, paste0("M16?_D8/GO_chr16/", i)))
}
gsea.m16x.d8 <- list()
for (i in names(deg.m16x.d8)) {
  gsea.m16x.d8[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.m16x.d8[[i]]$all, deg.type = "edger",
                                                             lfc = 1, sig = 0.05,
                                                             reversed = FALSE, species = "human",
                                                             basename = paste0(i, "_lfc1_pvalue0.05"),
                                                             genetype = "SYMBOL", gene.col = "SYMBOL",
                                                             outdir = file.path(res.out, paste0("GSEA/M16x_D8_lfc1_pvalue0.05/", i)))
  gsea.m16x.d8[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.m16x.d8[[i]]$all, deg.type = "edger",
                                                              lfc = 2, sig = 0.05,
                                                              reversed = FALSE, species = "human",
                                                              basename = paste0(i, "_lfc2_pvalue0.005"),
                                                              genetype = "SYMBOL", gene.col = "SYMBOL",
                                                              outdir = file.path(res.out, paste0("GSEA/M16x_D8_lfc2_pvalue0.005/", i)))
  gsea.m16x.d8[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.m16x.d8[[i]]$all, deg.type = "edger", 
                                                 lfc = 0, sig = 1, 
                                                 reversed = FALSE, species = "human", 
                                                 basename = paste0(i, "_all"), 
                                                 genetype = "SYMBOL", gene.col = "SYMBOL", 
                                                 outdir = file.path(res.out, paste0("GSEA/M16x_D8_all/", i)))
}
# M22: D8
go.m22.d8 <- list()
for (j in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.m22.d8)) {
    pd <- subset(deg.m22.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.m22.d8[[paste0(i, ".up.", j)]] <- Pipe.GO(species = "human", 
                                                 genelist = pd$SYMBOL,
                                                 basename = paste0(i, "_M22_vs_Euploidy_up_", j), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("M22_D8/GO/", i, "_", j)))
    pd <- subset(deg.m22.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.m22.d8[[paste0(i, ".down.", j)]] <- Pipe.GO(species = "human", 
                                                   genelist = pd$SYMBOL,
                                                   basename = paste0(i, "_M22_vs_Euploidy_down_", j), genetype = "SYMBOL",
                                                   res.out = file.path(res.out, paste0("M22_D8/GO/", i, "_", j)))
  }
}
for (i in names(deg.m22.d8)[4]) {
  print(length(intersect(chr22.gene$V1,
                         subset(deg.m22.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL)))
  go.m22.d8[[paste0(i, ".chr22.up")]] <- Pipe.GO(species = "human", 
                                                 genelist = intersect(chr22.gene$V1,
                                                                      subset(deg.m22.d8[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                 basename = paste0(i, "_M22_vs_Euploidy_chr22_up"), genetype = "SYMBOL",
                                                 res.out = file.path(res.out, paste0("M22_D8/GO_chr22/", i)))
  go.m22.d8[[paste0(i, ".chr22.down")]] <- Pipe.GO(species = "human", 
                                                   genelist = intersect(chr22.gene$V1,
                                                                        subset(deg.m22.d8[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                   basename = paste0(i, "_M22_vs_Euploidy_chr22_down"), genetype = "SYMBOL",
                                                   res.out = file.path(res.out, paste0("M22_D8/GO_chr22/", i)))
}
gsea.m22.d8 <- list()
for (i in names(deg.m22.d8)) {
  gsea.m22.d8[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.m22.d8[[i]]$all, deg.type = "edger",
                                                            lfc = 1, sig = 0.05,
                                                            reversed = FALSE, species = "human",
                                                            basename = paste0(i, "_lfc1_pvalue0.05"),
                                                            genetype = "SYMBOL", gene.col = "SYMBOL",
                                                            outdir = file.path(res.out, paste0("GSEA/M22_D8_lfc1_pvalue0.05/", i)))
  gsea.m22.d8[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.m22.d8[[i]]$all, deg.type = "edger",
                                                             lfc = 2, sig = 0.05,
                                                             reversed = FALSE, species = "human",
                                                             basename = paste0(i, "_lfc2_pvalue0.005"),
                                                             genetype = "SYMBOL", gene.col = "SYMBOL",
                                                             outdir = file.path(res.out, paste0("GSEA/M22_D8_lfc2_pvalue0.005/", i)))
  gsea.m22.d8[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.m22.d8[[i]]$all, deg.type = "edger", 
                                                lfc = 0, sig = 1, 
                                                reversed = FALSE, species = "human", 
                                                basename = paste0(i, "_all"), 
                                                genetype = "SYMBOL", gene.col = "SYMBOL", 
                                                outdir = file.path(res.out, paste0("GSEA/M22_D8_all/", i)))
}
# D10
go.t16.d10 <- list()
for (j in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.t16.d10)) {
    pd <- subset(deg.t16.d10[[i]]$sig, PValue <= 0.05 & logFC >= log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.t16.d10[[paste0(i, ".up.", j)]] <- Pipe.GO(species = "human", 
                                                  genelist = pd$SYMBOL,
                                                  basename = paste0(i, "_T16_vs_Euploidy_up_", j), genetype = "SYMBOL",
                                                  res.out = file.path(res.out, paste0("T16_D10/GO/", i, "_", j)))
    pd <- subset(deg.t16.d10[[i]]$sig, PValue <= 0.05 & logFC <= -log2(1.5))
    pd <- pd[rowSums(pd[, grep("^[1|2]_", colnames(pd), value = T)] >= j) >= 1, ]
    go.t16.d10[[paste0(i, ".down.", j)]] <- Pipe.GO(species = "human", 
                                                    genelist = pd$SYMBOL,
                                                    basename = paste0(i, "_T16_vs_Euploidy_down_", j), genetype = "SYMBOL",
                                                    res.out = file.path(res.out, paste0("T16_D10/GO/", i, "_", j)))
  }
}
for (i in names(deg.t16.d10)) {
  go.t16.d10[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                                  genelist = intersect(chr16.gene$V1,
                                                                       subset(deg.t16.d10[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                  basename = paste0(i, "_T16_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                                  res.out = file.path(res.out, paste0("T16_D10/GO_chr16/", i)))
  go.t16.d10[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                    genelist = intersect(chr16.gene$V1,
                                                                         subset(deg.t16.d10[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                    basename = paste0(i, "_T16_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                    res.out = file.path(res.out, paste0("T16_D10/GO_chr16/", i)))
}
gsea.t16.d10 <- list()
for (i in names(deg.t16.d10)) {
  gsea.t16.d10[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.t16.d10[[i]]$all, deg.type = "edger",
                                                             lfc = 1, sig = 0.05,
                                                             reversed = FALSE, species = "human",
                                                             basename = paste0(i, "_lfc1_pvalue0.05"),
                                                             genetype = "SYMBOL", gene.col = "SYMBOL",
                                                             outdir = file.path(res.out, paste0("GSEA/T16_D10_lfc1_pvalue0.05/", i)))
  gsea.t16.d10[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.t16.d10[[i]]$all, deg.type = "edger",
                                                              lfc = 2, sig = 0.05,
                                                              reversed = FALSE, species = "human",
                                                              basename = paste0(i, "_lfc2_pvalue0.005"),
                                                              genetype = "SYMBOL", gene.col = "SYMBOL",
                                                              outdir = file.path(res.out, paste0("GSEA/T16_D10_lfc2_pvalue0.005/", i)))
  gsea.t16.d10[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.t16.d10[[i]]$all, deg.type = "edger", 
                                                 lfc = 0, sig = 1, 
                                                 reversed = FALSE, species = "human", 
                                                 basename = paste0(i, "_all"), 
                                                 genetype = "SYMBOL", gene.col = "SYMBOL", 
                                                 outdir = file.path(res.out, paste0("GSEA/T16_D10_all/", i)))
}


### >>> 6. Overlapped DEGs
go.overlap <- list()
# T16 vs M16: D8 CTB
pd.list <- list(M16.up = subset(deg.m16.d8$M16_CTB$up, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                T16.up = subset(deg.t16.d8$T16_CTB$up, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                M16.down = subset(deg.m16.d8$M16_CTB$down, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                T16.down = subset(deg.t16.d8$T16_CTB$down, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                Chr16 = chr16.gene$V1)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), 
               pd.title = "T16 vs M16 (CTB in D8)", file.name = "D8_T16_vs_M16_CTB", 
               res.out = file.path(res.out, "Overlapped_DEGs/D8_T16_vs_M16_CTB"))
go.overlap[["CTB.up"]] <- Pipe.GO(species = "human", 
                                  genelist = intersect(pd.list$T16.up, pd.list$M16.up),
                                  basename = "CTB_M16_and_T16_up", genetype = "SYMBOL",
                                  res.out = file.path(res.out, "Overlapped_DEGs/CTB_GO"))
go.overlap[["CTB.down"]] <- Pipe.GO(species = "human", 
                                    genelist = intersect(pd.list$T16.down, pd.list$M16.down),
                                    basename = "CTB_M16_and_T16_down", genetype = "SYMBOL",
                                    res.out = file.path(res.out, "Overlapped_DEGs/CTB_GO"))
go.overlap[["CTB.up.chr16"]] <- Pipe.GO(species = "human", 
                                        genelist = intersect(chr16.gene$V1,
                                                             intersect(pd.list$T16.up, pd.list$M16.up)),
                                        basename = "CTB_M16_and_T16_chr16_up", genetype = "SYMBOL",
                                        res.out = file.path(res.out, "Overlapped_DEGs/CTB_GO_chr16"))
go.overlap[["CTB.down.chr16"]] <- Pipe.GO(species = "human", 
                                          genelist = intersect(chr16.gene$V1,
                                                               intersect(pd.list$T16.down, pd.list$M16.down)),
                                          basename = "CTB_M16_and_T16_chr16_down", genetype = "SYMBOL",
                                          res.out = file.path(res.out, "Overlapped_DEGs/CTB_GO_chr16"))
pd <- merge(deg.m16.d8$M16_CTB$all %>% mutate(cell = "M16"), 
            deg.t16.d8$T16_CTB$all %>% mutate(cell = "T16"), by = "SYMBOL")
pd %>% 
  ggplot(aes(x = logFC.x, y = logFC.y)) + 
  geom_point() + 
  geom_hline(yintercept = c(-2, 2), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  geom_vline(xintercept = c(-2, 2), colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5)
# T16 vs M16: D8 pre-STB
pd.list <- list(M16.up = subset(deg.m16.d8$`M16_pre-STB`$up, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                T16.up = subset(deg.t16.d8$`T16_pre-STB`$up, PValue <= 0.01 & logFC >= log2(4))$SYMBOL,
                M16.down = subset(deg.m16.d8$`M16_pre-STB`$down, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                T16.down = subset(deg.t16.d8$`T16_pre-STB`$down, PValue <= 0.01 & logFC <= -log2(4))$SYMBOL,
                Chr16 = chr16.gene$V1)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), 
               pd.title = "T16 vs M16 (pre-STB in D8)", file.name = "D8_T16_vs_M16_pre-STB", 
               res.out = file.path(res.out, "Overlapped_DEGs/D8_T16_vs_M16_pre-STB"))
go.overlap[["pre-STB.up"]] <- Pipe.GO(species = "human", 
                                      genelist = intersect(pd.list$T16.up, pd.list$M16.up),
                                      basename = "pre-STB_M16_and_T16_up", genetype = "SYMBOL",
                                      res.out = file.path(res.out, "Overlapped_DEGs/pre-STB_GO"))
go.overlap[["pre-STB.down"]] <- Pipe.GO(species = "human", 
                                        genelist = intersect(pd.list$T16.down, pd.list$M16.down),
                                        basename = "pre-STB_M16_and_T16_down", genetype = "SYMBOL",
                                        res.out = file.path(res.out, "Overlapped_DEGs/pre-STB_GO"))
go.overlap[["pre-STB.up.chr16"]] <- Pipe.GO(species = "human", 
                                            genelist = intersect(chr16.gene$V1,
                                                                 intersect(pd.list$T16.up, pd.list$M16.up)),
                                            basename = "pre-STB_M16_and_T16_chr16_up", genetype = "SYMBOL",
                                            res.out = file.path(res.out, "Overlapped_DEGs/pre-STB_GO_chr16"))
go.overlap[["pre-STB.down.chr16"]] <- Pipe.GO(species = "human", 
                                              genelist = intersect(chr16.gene$V1,
                                                                   intersect(pd.list$T16.down, pd.list$M16.down)),
                                              basename = "pre-STB_M16_and_T16_chr16_down", genetype = "SYMBOL",
                                              res.out = file.path(res.out, "Overlapped_DEGs/pre-STB_GO_chr16"))


### >>> 7. DEGs based on Karyotype
# DEGs
table(sr.yhw.no22$Stage, sr.yhw.no22$Karyotype)
sr.yhw.no22$MergedKaryotype <- gsub("\\?", "x", paste0(sr.yhw.no22$Stage, "-", sr.yhw.no22$Karyotype))
table(sr.yhw.no22$MergedKaryotype)
deg.kary <- list()
deg.meta <- data.frame(g1 = c("D10-Euploidy", "D8-Euploidy", "D8-Euploidy", "D8-Euploidy", "D8-Euploidy"), 
                       g2 = c("D10-T16", "D8-M16", "D8-M16x", "D8-M22", "D8-T16"))
for (i in 1:nrow(deg.meta)) {
  tmp.sr <- subset(sr.yhw.no22, MergedKaryotype %in% deg.meta[i,])
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "MergedKaryotype", 
                       comparison = c(deg.meta[i, "g1"], deg.meta[i, "g2"]), 
                       sample.n = NULL)
  tmp.name <- paste0(deg.meta[i, "g2"], "_vs_", deg.meta[i, "g1"])
  deg.kary[[tmp.name]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                      sample.n = NULL, group.by = "CellType", 
                                      g1 = paste0("1_", deg.meta[i, "g1"]), 
                                      g2 = paste0("2_", deg.meta[i, "g2"]), 
                                      lfc = 1, sig = 0.05,
                                      res.out = file.path(res.out, paste0("DEG/", tmp.name)))
}
# plot Volcano
for (i in names(deg.kary)) {
  pdf(file.path(res.out, paste0("DEG/", i, "_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.kary[[i]]$all, geneset = NULL, p.col = "FDR", lfc.col = "logFC",
                       sig = 0.005, lfc = 2, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
  pdf(file.path(res.out, paste0("DEG/", i, "_volcano_TFs.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.kary[[i]]$all, geneset = hs.tfs$Symbol, p.col = "FDR", lfc.col = "logFC",
                       sig = 0.005, lfc = 2, title = i,
                       pt.size = 3, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
go.kary <- list()
for (i in names(deg.kary)) {
  go.kary[[paste0(i, ".up")]] <- Pipe.GO(species = "human", 
                                         genelist = subset(deg.kary[[i]]$sig, FDR <= 0.005 & logFC >= log2(4))$SYMBOL,
                                         basename = paste0(i, "_up"), 
                                         genetype = "SYMBOL",
                                         res.out = file.path(res.out, paste0("GO/", i)))
  go.kary[[paste0(i, ".down")]] <- Pipe.GO(species = "human", 
                                           genelist = subset(deg.kary[[i]]$sig, FDR <= 0.005 & logFC <= -log2(4))$SYMBOL,
                                           basename = paste0(i, "_down"), 
                                           genetype = "SYMBOL",
                                           res.out = file.path(res.out, paste0("GO/", i)))
}
gsea.kary <- list()
for (i in names(deg.kary)) {
  gsea.kary[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.kary[[i]]$all, deg.type = "edger",
                                                          lfc = 0.05, sig = 1,
                                                          reversed = FALSE, species = "human",
                                                          basename = paste0(i, "_lfc1_pvalue0.05"),
                                                          genetype = "SYMBOL",
                                                          gene.col = "SYMBOL",
                                                          outdir = file.path(res.out, paste0("GSEA_karyotype/", i)))
  gsea.kary[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.kary[[i]]$all, deg.type = "edger",
                                                           lfc = 0.005, sig = 2,
                                                           reversed = FALSE, species = "human",
                                                           basename = paste0(i, "_lfc2_pvalue0.005"),
                                                           genetype = "SYMBOL",
                                                           gene.col = "SYMBOL",
                                                           outdir = file.path(res.out, paste0("GSEA_karyotype/", i)))
  gsea.kary[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.kary[[i]]$all, deg.type = "edger", 
                                              lfc = 0, sig = 1, 
                                              reversed = FALSE, species = "human", 
                                              basename = paste0(i, "_all"), 
                                              genetype = "SYMBOL", 
                                              gene.col = "SYMBOL", 
                                              outdir = file.path(res.out, paste0("GSEA_karyotype_all/", i)))
}
# DEGs stacked plot
gene.anno <- read.table("/home/laborer/refgenome/embl/rs110/homo_sapiens/anno/Homo_sapiens.GRCh38.110.chr.no.contig_genebody.bed",
                        sep = "\t") %>% 
  separate(col = "V4", into = c("ID", "SYMBOL", "Type"), sep = ":") %>% 
  unique()
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.kary)) {
    pdf(file.path(res.out, paste0(i, "_LFC_Stacked_Plot_", j, ".pdf")), height = 6, width = 12)
    print(VisDEG.StackBarByLFC(deg.data = deg.kary[[i]]$all, 
                               expr.by = grep("^[1|2]_", colnames(deg.kary[[i]]$all), value = T), low.expr = j,
                               gene.anno = gene.anno, chr.seq = c(paste0("chr", 1:22), "chrX", "chrY"), 
                               remove.nonsig = T, sig.col = "PValue", sig = 0.05, 
                               lfc.col = "logFC", split.by = NULL,
                               remove.chr = TRUE, chr.list = c("chrM")))
    dev.off()
  }
}
VisDEG.StackBarByLFC <- function(deg.data, gene.anno, expr.by = NULL, low.expr = 0.01,
                                 chr.seq = c(paste0("chr", 1:22), "chrX", "chrY"),
                                 remove.nonsig = TRUE, sig.col = "PValue", sig = 0.05,
                                 lfc.col = "logFC", split.by = "GroupLFC",
                                 remove.chr = TRUE, chr.list = c("chrM")) {
  # process data
  deg.data$log2FC <- deg.data[, lfc.col]
  # filter gene by expression level
  if (!is.null(expr.by)) {
    if (length(expr.by) >= 2) {
      index <- rowMeans(deg.data[, expr.by]) >= low.expr
    } else {
      index <- deg.data[, expr.by] >= low.expr
    }
    deg.data <- deg.data[index, ]
  }
  # filter gene by significance
  if (isTRUE(remove.nonsig)) {
    pd <- deg.data[deg.data[, sig.col] <= 0.05, ]
  }
  # merge data
  pd <- merge(pd, gene.anno, by = "SYMBOL") %>% 
    mutate(Grade = case_when(log2FC > log2(1) & log2FC <= log2(1.5) ~ "Up.FC:0-1.5",
                             log2FC > log2(1.5) & log2FC <= log2(3) ~ "Up.FC:1.5-3",
                             log2FC > log2(3) & log2FC <= log2(4.5) ~ "Up.FC:3-4.5",
                             log2FC > log2(4.5) & log2FC <= log2(6) ~ "Up.FC:4.5-6",
                             log2FC > log2(6) ~ "Up.FC:6-inf",
                             log2FC < -log2(1) & log2FC >= -log2(1.5) ~ "Down.FC:0-1.5",
                             log2FC < -log2(1.5) & log2FC >= -log2(3) ~ "Down.FC:1.5-3",
                             log2FC < -log2(3) & log2FC >= -log2(4.5) ~ "Down.FC:3-4.5",
                             log2FC < -log2(4.5) & log2FC >= -log2(6) ~ "Down.FC:4.5-6",
                             log2FC < -log2(6) ~ "Down.FC:6-inf"),
           GroupLFC = case_when(log2FC == 0 ~ "No.changed",
                                log2FC > 0 ~ "Up.regulated",
                                log2FC < 0 ~ "Down.regulated"),
           Value = 1) %>%
    dplyr::group_by(V1, GroupLFC, Grade) %>% 
    dplyr::summarise(Count = length(unique(SYMBOL)))
  # filter chromosome
  if (isTRUE(remove.chr)) {
    pd <- pd[!(pd$V1 %in% chr.list), ]
  }
  pd.gene <- gene.anno %>% 
    dplyr::group_by(V1) %>% 
    dplyr::summarise(TotalCount = length(unique(SYMBOL)))
  # make a stat
  pd <- merge(pd, pd.gene, by = "V1") %>% 
    dplyr::mutate(Ratio = Count/TotalCount) %>% 
    dplyr::mutate(Label = round(Ratio, 2))
  # plot
  if (!is.null(split.by)) {
    pd$split.col <- pd[, split.by]
    p <- pd %>% 
      ggplot(aes(x = V1, y = Ratio, fill = Grade)) +
      geom_bar(position = "stack", stat = "identity") +
      geom_text(aes(label = Label), position = position_stack(vjust= 0.5), colour = "white", size = 3) +
      scale_fill_brewer(palette = "Paired") + 
      scale_x_discrete(limits = chr.seq) +
      facet_wrap(. ~ split.col) + 
      theme_bw() +
      theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
            axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
            axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
            axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
            panel.grid = element_blank())
  } else {
    p <- pd %>% 
      ggplot(aes(x = V1, y = Ratio, fill = Grade)) +
      geom_bar(position = "stack", stat = "identity") +
      geom_text(aes(label = Label), position = position_stack(vjust= 0.5), colour = "white", size = 3) +
      scale_fill_brewer(palette = "Paired") + 
      scale_x_discrete(limits = chr.seq) +
      theme_bw() +
      theme(axis.title.x = element_text(face = "plain", colour = "#000000", size = 14, angle = 0),
            axis.title.y = element_text(face = "plain", colour = "#000000", size = 14, angle = 90),
            axis.text.x  = element_text(face = "plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
            axis.text.y  = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
            panel.grid = element_blank())
  }
  return(p)
}



# =========================
# 5th part: Repeat analysis ----
# =========================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Repeat")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Load count data
# repeat locus
repeat.count <- read.table("analysis/scrna/results/featurecounts/repeat/all_samples_repeat_count_locus_unique_matrix.txt", 
                           header = T, sep = "\t", row.names = 1)
colnames(repeat.count)[-1:-5] <- gsub("\\.", "_", 
                                      gsub("d", "D", 
                                           gsub("_R1.*", "", 
                                                gsub("X.*gene.", "", colnames(repeat.count)[-1:-5]))))
repeat.count <- repeat.count[rowSums(repeat.count[, -1:-5]) >= 50, ]
if (all(colnames(repeat.count)[-1:-5] == colnames(sr.yhw.no22))) {
  sr.yhw.no22[["Locus"]] <- CreateAssayObject(count = repeat.count[, -1:-5])
}
DefaultAssay(sr.yhw.no22) <- "Locus"
sr.yhw.no22 <- NormalizeData(sr.yhw.no22) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.yhw.no22))
rm(repeat.count)
# repeat family
repeat.count <- read.table("analysis/scrna/results/featurecounts/repeat/all_samples_repeat_count_family_multi_matrix.txt", 
                           header = T, sep = "\t", row.names = 1)
dim(repeat.count)
colnames(repeat.count)[-1:-5] <- gsub("\\.", "_", 
                                      gsub("d", "D", 
                                           gsub("_R1.*", "", 
                                                gsub("X.*gene.", "", colnames(repeat.count)[-1:-5]))))
repeat.tpm <- CountToTpm(repeat.count[, -1:-5], repeat.count$Length)
dim(repeat.tpm)
if (all(colnames(repeat.tpm) == colnames(sr.yhw.no22))) {
  write.csv(repeat.tpm, file.path(res.out, "All_sample_repeat_expression_level_TPM.csv"), 
            row.names = T, col.names = T, quote = F)
  saveRDS(repeat.tpm, file.path(res.out, "All_sample_repeat_expression_level_TPM.rds"))
}
repeat.cpm <- GetAssayData(sr.yhw.no22, slot = "data", assay = "Family")
dim(repeat.cpm)
write.csv(repeat.cpm, file.path(res.out, "All_sample_repeat_expression_level_log-normalized_count.csv"), 
          row.names = T, col.names = T, quote = F)
saveRDS(repeat.cpm, file.path(res.out, "All_sample_repeat_expression_level_log-normalized_count.rds"))
if (all(colnames(repeat.count)[-1:-5] == colnames(sr.yhw.no22))) {
  sr.yhw.no22[["Family"]] <- CreateAssayObject(count = repeat.count[, -1:-5])
}
DefaultAssay(sr.yhw.no22) <- "Family"
sr.yhw.no22 <- NormalizeData(sr.yhw.no22) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 100) %>%
  ScaleData(features = rownames(sr.yhw.no22))
rm(repeat.count)


### >>> 3. Differential expression analysis (family level)
# find differential families (cell type)
DefaultAssay(sr.yhw.no22) <- "Family"
Idents(sr.yhw.no22) <- sr.yhw.no22$CellType
sr.yhw.no22.family.markers <- FindAllMarkers(sr.yhw.no22, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
pd.gene <- sr.yhw.no22.family.markers %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "Family",
  group.by = c("CellType"), exp_cutoff = "15",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_Family_expression_Dot_plot_after_annotation.pdf"), 
    height = 14, width = 6)
print(ht$plot)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_Family_expression_Scatter_plot_after_annotation.pdf"), 
    height = 30, width = 30)
FeatureDimPlot(srt = sr.yhw.no22.scenic, features = pd.gene$gene, assay = "Family", 
               reduction = "UMAP", theme_use = "theme_blank")
dev.off()
# find differential families (Karyotype)
DefaultAssay(sr.yhw.no22) <- "Family"
Idents(sr.yhw.no22) <- sr.yhw.no22$Karyotype
sr.yhw.no22.family.markers2 <- FindAllMarkers(sr.yhw.no22, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
pd.gene <- sr.yhw.no22.family.markers2 %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "Family",
  group.by = c("Karyotype"), exp_cutoff = "15",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_Family_expression_Dot_plot_after_annotation_Karyotype.pdf"), 
    height = 14, width = 6)
print(ht$plot)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_Family_expression_Scatter_plot_after_annotation_Karyotype.pdf"), 
    height = 30, width = 30)
FeatureDimPlot(srt = sr.yhw.no22.scenic, features = pd.gene$gene, assay = "Family", 
               reduction = "UMAP", theme_use = "theme_blank")
dev.off()
# T16 vs Euploidy
def.t16 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "T16"), sample.n = NULL, assay = "Family")
  def.t16[[paste0("T16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_T16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("DEFs/T16/", i, "_T16_vs_Euploidy")))
}
for (i in names(def.t16)) {
  pdf(file.path(res.out, paste0("DEFs/T16/", i, "_T16_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = def.t16[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M22 vs Euploidy
def.m22 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M22"), sample.n = NULL, assay = "Family")
  def.m22[[paste0("M22_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_M22", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("DEFs/M22/", i, "_M22_vs_Euploidy")))
}
for (i in names(def.m22)) {
  pdf(file.path(res.out, paste0("DEFs/M22/", i, "_M22_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = def.m22[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M16? vs Euploidy
def.m16x <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp.sr$Karyotype <- gsub("M16.", "M16x", tmp.sr$Karyotype)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16x"), sample.n = NULL, assay = "Family")
  def.m16x[[paste0("M16x_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                sample.n = NULL, group.by = "CellType",
                                                g1 = "1_Euploidy", g2 = "2_M16x", lfc = 1, sig = 0.05,
                                                res.out = file.path(res.out, paste0("DEFs/M16?/", i, "_M16x_vs_Euploidy")))
}
for (i in names(def.m16x)) {
  pdf(file.path(res.out, paste0("DEFs/M16?/", i, "_M16x_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = def.m16x[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M16 vs Euploidy
def.m16 <- list()
for (i in c("CTB", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16"), sample.n = NULL, assay = "Family")
  def.m16[[paste0("M16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_M16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("DEFs/M16/", i, "_M16_vs_Euploidy")))
}
for (i in names(def.m16)) {
  pdf(file.path(res.out, paste0("DEFs/M16/", i, "_M16_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = def.m16[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}


### >>> 4. Differential expression analysis (locus level)
# T16 vs Euploidy
del.t16 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "T16"), sample.n = NULL, assay = "Locus")
  del.t16[[paste0("T16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_T16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("DELs/T16/", i, "_T16_vs_Euploidy")))
}
for (i in names(del.t16)) {
  pdf(file.path(res.out, paste0("DELs/T16/", i, "_T16_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = del.t16[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 1, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M22 vs Euploidy
del.m22 <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M22"), sample.n = NULL, assay = "Locus")
  del.m22[[paste0("M22_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_M22", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("DELs/M22/", i, "_M22_vs_Euploidy")))
}
for (i in names(del.m22)) {
  pdf(file.path(res.out, paste0("DELs/M22/", i, "_M22_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = del.m22[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 1, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M16? vs Euploidy
del.m16x <- list()
for (i in c("CTB", "EPI", "Hypoblast", "pre-STB", "STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp.sr$Karyotype <- gsub("M16.", "M16x", tmp.sr$Karyotype)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16x"), sample.n = NULL, assay = "Locus")
  del.m16x[[paste0("M16x_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                sample.n = NULL, group.by = "CellType",
                                                g1 = "1_Euploidy", g2 = "2_M16x", lfc = 1, sig = 0.05,
                                                res.out = file.path(res.out, paste0("DELs/M16?/", i, "_M16x_vs_Euploidy")))
}
for (i in names(del.m16x)) {
  pdf(file.path(res.out, paste0("DELs/M16?/", i, "_M16x_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = del.m16x[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 1, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# M16 vs Euploidy
del.m16 <- list()
for (i in c("CTB", "pre-STB")) {
  tmp.sr <- subset(sr.yhw.no22, CellType == i)
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("Euploidy", "M16"), sample.n = NULL, assay = "Locus")
  del.m16[[paste0("M16_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                              sample.n = NULL, group.by = "CellType",
                                              g1 = "1_Euploidy", g2 = "2_M16", lfc = 1, sig = 0.05,
                                              res.out = file.path(res.out, paste0("DELs/M16/", i, "_M16_vs_Euploidy")))
}
for (i in names(del.m16)) {
  pdf(file.path(res.out, paste0("DELs/M16/", i, "_M16_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = del.m16[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 1, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}



# =========================
# 6th part: SCENIC analysis ----
# =========================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/SCENIC")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Output loom file
library("AfterChat")
library("SCopeLoomR")
DefaultAssay(sr.yhw.no22) <- "RNA"
WriteLoom.new(sc.obj = sr.yhw.no22, loom.file = file.path(res.out, "sr.yhw.no22.loom"))


### >>> 3. Load SCENIC results
sr.yhw.no22@meta.data
sr.yhw.no22.scenic <- ScenicAUC(sr.obj = sr.yhw.no22, workdir = file.path(getwd(), "R/Graphs/SCENIC"), 
                                gp.name = "CellType", sp.name = "XuYanWen_Aneuploid_no_chr22", 
                                pd.gene = markers, pt.size = 2, pd.order = T)
tmp <- ScenicAUC(sr.obj = sr.yhw.no22, workdir = file.path(getwd(), "R/Graphs/SCENIC_new"), 
                 gp.name = "MergedCellType", sp.name = "XuYanWen_Aneuploid_no_chr22", 
                 pd.gene = markers, pt.size = 2, pd.order = T)
rm(tmp)
DefaultAssay(sr.yhw.no22.scenic) <- "Regulon_AUC"
Idents(sr.yhw.no22.scenic) <- sr.yhw.no22.scenic$CellType
# find differential TFs
sr.yhw.no22.scenic.markers <- FindAllMarkers(sr.yhw.no22.scenic, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
# plot marker TFs activity
pd.gene <- sr.yhw.no22.scenic.markers %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "Regulon_AUC",
  group.by = c("CellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_activity_Dot_plot_after_annotation.pdf"), 
    height = 12, width = 6)
print(ht$plot)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_activity_Scatter_plot_after_annotation.pdf"), 
    height = 20, width = 20)
FeatureDimPlot(srt = sr.yhw.no22.scenic, features = pd.gene$gene, assay = "Regulon_AUC", 
               reduction = "UMAP", theme_use = "theme_blank")
dev.off()
# plot marker TFs expression
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "RNA",
  group.by = c("CellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_expression_Dot_plot_after_annotation.pdf"), 
    height = 12, width = 6)
print(ht$plot)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_expression_Scatter_plot_after_annotation.pdf"), 
    height = 20, width = 20)
FeatureDimPlot(srt = sr.yhw.no22.scenic, features = pd.gene$gene, assay = "RNA", 
               reduction = "UMAP", theme_use = "theme_blank")
dev.off()


### >>> 4. Show marker TFs
sr.yhw.no22.scenic$MergedCellType <- paste(sr.yhw.no22.scenic$CellType,
                                           sr.yhw.no22.scenic$Stage,
                                           sr.yhw.no22.scenic$Karyotype, sep = "_")
Idents(sr.yhw.no22.scenic) <- sr.yhw.no22.scenic$MergedCellType
sr.yhw.no22.scenic$MergedCellType <- factor(sr.yhw.no22.scenic$MergedCellType)
# find differential TFs
sr.yhw.no22.scenic.markers2 <- FindAllMarkers(sr.yhw.no22.scenic, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
# plot marker TFs activity (top15)
pd.gene <- sr.yhw.no22.scenic.markers %>% 
  group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "Regulon_AUC",
  group.by = c("MergedCellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm"), cluster_rows = T, cluster_columns = T
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_activity_Dot_plot_more_groups_Top15_after_annotation.pdf"), 
    height = 20, width = 20)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "RNA",
  group.by = c("MergedCellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_expression_Dot_plot_more_groups_Top15_after_annotation.pdf"), 
    height = 20, width = 20)
print(ht$plot)
dev.off()
# plot marker TFs activity (top30)
pd.gene <- sr.yhw.no22.scenic.markers %>% 
  group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "Regulon_AUC",
  group.by = c("MergedCellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm"), cluster_rows = T, cluster_columns = T
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_activity_Dot_plot_more_groups_Top30_after_annotation.pdf"), 
    height = 35, width = 20)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "RNA",
  group.by = c("MergedCellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_expression_Dot_plot_more_groups_Top30_after_annotation.pdf"), 
    height = 20, width = 20)
print(ht$plot)
dev.off()
# plot marker TFs activity (top50)
pd.gene <- sr.yhw.no22.scenic.markers %>% 
  group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "Regulon_AUC",
  group.by = c("MergedCellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, show_column_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm"), cluster_rows = T, cluster_columns = T
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_activity_Dot_plot_more_groups_Top50_after_annotation.pdf"), 
    height = 42, width = 20)
print(ht$plot)
dev.off()
ht <- GroupHeatmap(
  srt = sr.yhw.no22.scenic, 
  features = pd.gene$gene, assay = "RNA",
  group.by = c("MergedCellType"), exp_cutoff = "0.1",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "S.Score"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE, dot_size = unit(8, "mm")
)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr22_marker_TFs_expression_Dot_plot_more_groups_Top50_after_annotation.pdf"), 
    height = 20, width = 20)
print(ht$plot)
dev.off()


### >>> 5. Differential activity analysis
dir.create(file.path(res.out, "DE_analysis"), recursive = T)
sr.yhw.no22.scenic$MergedCellType <- gsub("\\?", "x", sr.yhw.no22.scenic$MergedCellType)
table(gsub("\\?", "x", sr.yhw.no22.scenic$MergedCellType)) %>% as.data.frame()
for (i in grep("Euploidy", unique(sr.yhw.no22.scenic$MergedCellType), invert = T, value = T)) {
  group1.cell <- length(grep(paste0("^", i, "$"), sr.yhw.no22.scenic$MergedCellType))
  group2.cell <- length(grep(paste0("^", gsub("_[T|M].*", "_Euploidy", i), "$"), sr.yhw.no22.scenic$MergedCellType))
  if (group1.cell >= 3 & group2.cell >= 3) {
    sr.yhw.no22.scenic <- RunDEtest(srt = sr.yhw.no22.scenic, group_by = "MergedCellType", 
                                    fc.threshold = 1, only.pos = FALSE, assay = "Regulon_AUC",
                                    group1 = i, group2 = gsub("_[T|M].*", "_Euploidy", i))
    pdf(file.path(res.out, paste0("DE_analysis/TFs_activity_Volcano_of_", i, ".pdf")), 
        height = 7, width = 8)
    print(VolcanoPlot(srt = sr.yhw.no22.scenic, nlabel = 15, pt.size = 1))
    dev.off()
    write.csv(sr.yhw.no22.scenic@tools$DEtest_custom$AllMarkers_wilcox, 
              file.path(res.out, paste0("DE_analysis/TFs_activity_Volcano_of_", i, ".csv")),
              row.names = F, col.names = T)
  }
}


### >>> 6. Regulatory network
dir.create(file.path(res.out, "Network"), recursive = T)
# define regulators
pd.gene <- sr.yhw.no22.scenic.markers %>% 
  group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
#
library("GENIE3")
set.seed(123)
exprMatr <- GetAssayData(sr.yhw.no22, assay = "RNA", slot = "data") %>% as.matrix()
exprMatr[1:5,1:5]
tf.exprMatr <- exprMatr[intersect(pd.gene$gene, rownames(exprMatr)),]
weightMat <- GENIE3(exprMatr, regulators = rownames(tf.exprMatr)[rowSums(tf.exprMatr) > 0],
                    nCores = 16, verbose = TRUE)
# Get all the regulatory links
linkList <- getLinkList(weightMat)
dim(linkList)
write.csv(linkList, file.path(res.out, "Network/Regulon_by_TFs_all_links.csv"), row.names = F)
# Get only the top-ranked links
linkList <- getLinkList(weightMat, reportMax = 1000)
dim(linkList)
write.csv(linkList, file.path(res.out, "Network/Regulon_by_TFs_top1000_links.csv"), row.names = F)
# Get only the links with a weight higher than some threshold
linkList <- getLinkList(weightMat, threshold = 0.1)
dim(linkList)
write.csv(linkList, file.path(res.out, "Network/Regulon_by_Weight_0.1_links.csv"), row.names = F)
# load TFs
net.tfs <- read.table("R/Table/DEG_TFsActivity_regulators.txt") %>% unique()
net.final <- subset(linkList, regulatoryGene %in% net.tfs$V1) %>% 
  mutate(regulatoryGene = as.character(regulatoryGene),
         targetGene = as.character(targetGene))
# add targets LFC
net.final <- merge(net.final, deg.t16.d8$T16_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_T16_CTB", "Targets_Expr.D8_T16_CTB")
net.final <- merge(net.final, deg.t16.d8$`T16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_T16_pre-STB", "Targets_Expr.D8_T16_pre-STB")
net.final <- merge(net.final, deg.t16.d8$T16_STB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_T16_STB", "Targets_Expr.D8_T16_STB")
net.final <- merge(net.final, deg.t16.d10$`T16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D10_T16_pre-STB", "Targets_Expr.D10_T16_pre-STB")
net.final <- merge(net.final, deg.t16.d10$T16_STB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D10_T16_STB", "Targets_Expr.D10_T16_STB")
net.final <- merge(net.final, deg.m16.d8$M16_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M16_CTB", "Targets_Expr.D8_M16_CTB")
net.final <- merge(net.final, deg.m16.d8$`M16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M16_pre-STB", "Targets_Expr.D8_M16_pre-STB")
net.final <- merge(net.final, deg.m22.d8$M22_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M22_CTB", "Targets_Expr.D8_M22_CTB")
net.final <- merge(net.final, deg.m22.d8$`M22_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M22_pre-STB", "Targets_Expr.D8_M22_pre-STB")
net.final <- merge(net.final, deg.m16x.d8$M16_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M16x_CTB", "Targets_Expr.D8_M16x_CTB")
net.final <- merge(net.final, deg.m16x.d8$`M16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M16x_pre-STB", "Targets_Expr.D8_M16x_pre-STB")
net.final <- merge(net.final, deg.m16x.d8$M16_STB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(net.final)[(ncol(net.final)-1):ncol(net.final)] <- c("Targets_LFC.D8_M16x_STB", "Targets_Expr.D8_M16x_STB")
# add TF activity
tfs.act <- AverageExpression(sr.yhw.no22.scenic, slot = "data", group.by = "MergedCellType")
tfs.act$RNA <- as.data.frame(tfs.act$RNA)
tfs.act$RNA$SYMBOL <- rownames(tfs.act$RNA)
tfs.act$Regulon_AUC <- as.data.frame(tfs.act$Regulon_AUC)
tfs.act$Regulon_AUC$SYMBOL <- rownames(tfs.act$Regulon_AUC)
net.final <- merge(net.final, tfs.act$RNA, by.x = "regulatoryGene", by.y = "SYMBOL")
col.index <- (ncol(net.final)-31):ncol(net.final)
colnames(net.final)[col.index] <- paste0("TFs_Expr_", colnames(net.final)[col.index])
net.final <- merge(net.final, tfs.act$Regulon_AUC, by.x = "regulatoryGene", by.y = "SYMBOL")
col.index <- (ncol(net.final)-31):ncol(net.final)
colnames(net.final)[col.index] <- paste0("TFs_Act_", colnames(net.final)[col.index])
# output network
write.csv(net.final, file.path(res.out, "Network/Regulon_by_Weight_0.1_network.csv"), row.names = F, quote = F)


### >>> 7. SCENIC regulatory network
# - "FOXO1", "ELF3", "CLOCK", "KMT2A", "KMT2B", "CEBPA", "CEBPB"
# FOXO1 in M16 pre-STB和T16_D8 pre-STB
# ELF3 in M16 CTB
# CLOCK in T16 D10 STB
# KMT2A in T16_D8 pre-STB
# KMT2B in T16_D8 CTB
# load data
grep("CEBP", scenic.net$regulatoryGene, value = T)
scenic.net <- read.table("/home/yhw/bioinfo/project-xuyanwen/aneuploid/R/Graphs/SCENIC/pySCENIC_expr_mat.adjacencies.cor.tfs.tsv", sep = "\t")
colnames(scenic.net) <- c("regulatoryGene", "targetGene", "importance", "regulation", "rho")
scenic.net <- scenic.net %>% 
  filter(regulatoryGene %in% c("FOXO1", "ELF3", "CLOCK", "KMT2A", "KMT2B") & regulation == 1) %>% 
  filter(rho >= mean(rho) & importance >= mean(importance))
# GO analysis
scenic.go <- list()
for (i in unique(scenic.net$regulatoryGene)) {
  scenic.go[[i]] <- Pipe.GO(species = "human", 
                            genelist = subset(scenic.net, regulatoryGene == i)$targetGene,
                            basename = paste0(i, "_regulated_genes"), genetype = "SYMBOL",
                            res.out = file.path(res.out, paste0("SCENIC_GO/", i)))
}
go.data <- list(FOXO1 = scenic.go$FOXO1,
                ELF3 = scenic.go$ELF3,
                KMT2A = scenic.go$KMT2A,
                KMT2B = scenic.go$KMT2B,
                CLOCK = scenic.go$CLOCK)
go.terms <- c("GO:0001890", "GO:0001892", "GO:0001893", "GO:0060669", "GO:0060674", "GO:0060706",
              "GO:0060707", "GO:0061450", "GO:1901163",
              "GO:0001829")
pd <- GO.replot(go.data = go.data, go.terms = go.terms, 
                multi.go = T, col.max = 3, pd.title = "TFs-regulated genes")
pdf(file.path(res.out, "SCENIC_GO/Replot_of_GO_analysis_TFs-regulated_genes.pdf"), 
    height = 6, width = 8)
pd$plot
dev.off()
# build network
scenic.net <- read.table("/home/yhw/bioinfo/project-xuyanwen/aneuploid/R/Graphs/SCENIC/pySCENIC_expr_mat.adjacencies.cor.tfs.tsv", sep = "\t")
colnames(scenic.net) <- c("regulatoryGene", "targetGene", "importance", "regulation", "rho")
scenic.net <- rbind(scenic.net %>% 
                      filter(regulatoryGene %in% c("FOXO1", "ELF3", "CLOCK", "KMT2A", "KMT2B") & regulation == 1) %>% 
                      filter(rho >= mean(rho) & importance >= mean(importance)),
                    scenic.net %>% 
                      filter(regulatoryGene %in% c("CEBPA", "CEBPB") & regulation == 1))
genie3 <- linkList %>% 
  mutate(regulation = 1, rho = weight) %>% 
  filter(regulatoryGene %in% c("CEBPA", "CEBPB"))
colnames(genie3) <- colnames(scenic.net)
scenic.net <- rbind(scenic.net, genie3)
scenic.net <- merge(scenic.net, deg.t16.d8$T16_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_T16_CTB", "Targets_Expr.D8_T16_CTB")
scenic.net <- merge(scenic.net, deg.t16.d8$`T16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_T16_pre-STB", "Targets_Expr.D8_T16_pre-STB")
scenic.net <- merge(scenic.net, deg.t16.d8$T16_STB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_T16_STB", "Targets_Expr.D8_T16_STB")
scenic.net <- merge(scenic.net, deg.t16.d10$`T16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D10_T16_pre-STB", "Targets_Expr.D10_T16_pre-STB")
scenic.net <- merge(scenic.net, deg.t16.d10$T16_STB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D10_T16_STB", "Targets_Expr.D10_T16_STB")
scenic.net <- merge(scenic.net, deg.m16.d8$M16_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M16_CTB", "Targets_Expr.D8_M16_CTB")
scenic.net <- merge(scenic.net, deg.m16.d8$`M16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M16_pre-STB", "Targets_Expr.D8_M16_pre-STB")
scenic.net <- merge(scenic.net, deg.m22.d8$M22_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M22_CTB", "Targets_Expr.D8_M22_CTB")
scenic.net <- merge(scenic.net, deg.m22.d8$`M22_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M22_pre-STB", "Targets_Expr.D8_M22_pre-STB")
scenic.net <- merge(scenic.net, deg.m16x.d8$M16_CTB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M16x_CTB", "Targets_Expr.D8_M16x_CTB")
scenic.net <- merge(scenic.net, deg.m16x.d8$`M16_pre-STB`$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M16x_pre-STB", "Targets_Expr.D8_M16x_pre-STB")
scenic.net <- merge(scenic.net, deg.m16x.d8$M16_STB$all[, c(1, 2, 8)], by.x = "targetGene", by.y = "SYMBOL")
colnames(scenic.net)[(ncol(scenic.net)-1):ncol(scenic.net)] <- c("Targets_LFC.D8_M16x_STB", "Targets_Expr.D8_M16x_STB")
tfs.act <- AverageExpression(sr.yhw.no22.scenic, slot = "data", group.by = "MergedCellType")
tfs.act$RNA <- as.data.frame(tfs.act$RNA)
tfs.act$RNA$SYMBOL <- rownames(tfs.act$RNA)
tfs.act$Regulon_AUC <- as.data.frame(tfs.act$Regulon_AUC)
tfs.act$Regulon_AUC$SYMBOL <- rownames(tfs.act$Regulon_AUC)
scenic.net <- merge(scenic.net, tfs.act$RNA, by.x = "regulatoryGene", by.y = "SYMBOL")
col.index <- (ncol(scenic.net)-31):ncol(scenic.net)
colnames(scenic.net)[col.index] <- paste0("TFs_Expr_", colnames(scenic.net)[col.index])
scenic.net <- merge(scenic.net, tfs.act$Regulon_AUC, by.x = "regulatoryGene", by.y = "SYMBOL")
col.index <- (ncol(scenic.net)-31):ncol(scenic.net)
colnames(scenic.net)[col.index] <- paste0("TFs_Act_", colnames(scenic.net)[col.index])
write.csv(subset(scenic.net, regulatoryGene == "ELF3"), 
          file.path(res.out, "SCENIC_GO/Regulon_by_ELF3_network.csv"), row.names = F, quote = F)
write.csv(subset(scenic.net, regulatoryGene == "FOXO1"), 
          file.path(res.out, "SCENIC_GO/Regulon_by_FOXO1_network.csv"), row.names = F, quote = F)
write.csv(subset(scenic.net, regulatoryGene %in% c("ELF3", "FOXO1")), 
          file.path(res.out, "SCENIC_GO/Regulon_by_ELF3_and_FOXO1_network.csv"), row.names = F, quote = F)
write.csv(subset(scenic.net, regulatoryGene %in% c("CEBPA", "CEBPB", "ELF3", "FOXO1")), 
          file.path(res.out, "SCENIC_GO/Regulon_by_CEBPA_CEBPB_ELF3_and_FOXO1_network.csv"), row.names = F, quote = F)



# ========================================
# 7th part: Infer developmental trajectory ----
# ========================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Dev_Trajectory")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Load data
hs.embryo <- readRDS("/home/yhw/document/public_data/E-MTAB-3929_human_pre-implantation/R/CodeData/CMQ_GSE136447_seurat_object.rds")
DimPlot(hs.embryo, reduction = "umap", pt.size = 2, cols = pd.col, group.by = "CellType", label = T)
hs.embryo@meta.data


### >>> 3. Integration by Harmony
# integration
DefaultAssay(sr.yhw.no22) <- "RNA"
sr.merge <- merge(sr.yhw.no22, y = hs.embryo, add.cell.ids = c("ThisStudy", "GSE136447"), project = "Merged")
sr.merge@meta.data %>% 
  mutate(DataSet = case_when(orig.ident %in% c("D8", "D10") ~ "ThisStudy",
                             orig.ident %in% c("GSE136447") ~ "GSE136447"),
         CellType = gsub("hs", "", CellType),
         Karyotype = case_when(orig.ident %in% c("D8", "D10") ~ Karyotype,
                               orig.ident %in% c("GSE136447") ~ "GSE136447"),
         Stage = case_when(orig.ident %in% c("D8", "D10") ~ Stage,
                           orig.ident %in% c("GSE136447") ~ Day)) -> sr.merge@meta.data
sr.merge@meta.data <- sr.merge@meta.data[, c("DataSet", "CellType", "Karyotype", "Stage")]
sr.merge <- NormalizeData(sr.merge) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sr.merge))
pd.gene <- setdiff(VariableFeatures(object = sr.merge), c(chr16.gene$V1, chr22.gene$V1))
sr.merge <- RunPCA(sr.merge, features = pd.gene)
ElbowPlot(sr.merge, ndims = 30)
dev.off()
sr.merge <- RunUMAP(sr.merge, dims = 1:20)
DimPlot(sr.merge, reduction = "umap", group.by = c("DataSet", "CellType"), cols = pd.col, label = T)
sr.merge <- RunHarmony(sr.merge, group.by.vars = "DataSet")
sr.merge <- RunUMAP(sr.merge, reduction = "harmony", dims = 1:20)
# correlation
table(sr.merge$Stage, sr.merge$CellType)
sr.merge$MergedCellType <- paste(sr.merge$Stage, sr.merge$Karyotype, sr.merge$CellType, sep = "_")
sr.merge.list <- SplitObject(sr.merge, split.by = "DataSet")
ht <- CellCorHeatmap(
  srt_query = sr.merge.list$ThisStudy, srt_ref = sr.merge.list$GSE136447,
  query_group = "MergedCellType", ref_group = "MergedCellType",
  nlabel = 3, label_by = "row", cluster_rows = T, cluster_columns = T,
  show_row_names = TRUE, show_column_names = TRUE
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_correlation_plot.pdf"), height = 10, width = 25)
print(ht$plot)
dev.off()
# plot cells
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot_highlighted_by_T16.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T,
            cells.highlight = colnames(sr.merge)[sr.merge$Karyotype == "T16"], sizes.highlight = 0.75)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot_highlighted_by_M16.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T,
            cells.highlight = colnames(sr.merge)[sr.merge$Karyotype == "M16"], sizes.highlight = 0.75)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot_highlighted_by_M16?.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T,
            cells.highlight = colnames(sr.merge)[sr.merge$Karyotype == "M16?"], sizes.highlight = 0.75)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot_highlighted_by_M22.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T,
            cells.highlight = colnames(sr.merge)[sr.merge$Karyotype == "M22"], sizes.highlight = 0.75)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot_highlighted_by_Euploidy.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T,
            cells.highlight = colnames(sr.merge)[sr.merge$Karyotype == "Euploidy"], sizes.highlight = 0.75)
dev.off()
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_clustering_umap_plot_highlighted_by_GSE136447.pdf"), height = 11, width = 11)
CellDimPlot(srt = sr.merge, group.by = c("DataSet", "CellType", "Karyotype", "Stage"),
            reduction = "umap", theme_use = "theme_blank", label = T, pt.size = 0.5, label_repel = T,
            cells.highlight = colnames(sr.merge)[sr.merge$Karyotype == "GSE136447"], sizes.highlight = 0.75)
dev.off()


### >>> 4. Monocle3
sr.merge$nCount_RNA = colSums(x = sr.merge, slot = "counts")
sr.merge$nFeature_RNA = colSums(x = GetAssayData(object =sr.merge, slot = "counts") > 0)
sr.merge[["percent.mt"]] <- PercentageFeatureSet(sr.merge, pattern = "^MT-")
DimPlot(sr.merge)
tmp.cor <- data.frame(xmin = c(0, -20), xmax = c(15, -10), ymin = c(-10, -10), ymax = c(10, 0),
                      group = c("TE-lineage", "ICM-lineage"))
test <- ExtractCellByPos(object = sr.merge, object.type = "seurat", dim.name = "umap",
                         group.coor = tmp.cor, group.name = "test", pt.size = 1)
# TE-lineage
tmp <- sr.merge[,colnames(sr.merge) %in% test$id$`TE-lineage`]
mn3 <- list()
mn3$TE <- TI.Mncl3(sr.ob = tmp, btc.column = "DataSet", gp.column = "CellType",
                   mt.column = "percent.mt", use.sr.umap = T,
                   cluster.res = NULL, root.gp = "TE", sp.name = "TE_lineage",
                   pt.size = 1, out.dir = file.path(res.out, "Monocle3"))
mn3$TE.gene <- GeneModule.Mncl3(mncl3.cds = mn3$TE, sp.name = "TE_lineage_gene_module", 
                                gp.column = "CellType", ncores = 16, qvalue = 0.05, 
                                out.dir = file.path(res.out, "Monocle3"))
if (all(names(pseudotime(mn3$TE)) == colnames(tmp))) {
  tmp$Pseudotime <- log2(pseudotime(mn3$TE)+1)
}
tmp$MergedCellType <- paste(tmp$Stage, tmp$Karyotype, tmp$CellType, sep = "_")
tmp <- subset(tmp, DataSet == "ThisStudy")
p1 <- FeatureStatPlot(
  srt = tmp[,grep("^STB$", tmp$CellType)], group.by = "MergedCellType",
  stat.by = c("Pseudotime"), add_box = TRUE,
  comparisons = list(
    c("D10_T16_STB", "D10_Euploidy_STB"),
    c("D8_T16_STB", "D8_Euploidy_STB"),
    c("D8_M22_STB", "D8_Euploidy_STB"),
    c("D8_M16?_STB", "D8_Euploidy_STB")
  )
)
p2 <- FeatureStatPlot(
  srt = tmp[,grep("^pre-STB$", tmp$CellType)], group.by = "MergedCellType",
  stat.by = c("Pseudotime"), add_box = TRUE,
  comparisons = list(
    c("D10_T16_pre-STB", "D10_Euploidy_pre-STB"),
    c("D8_T16_pre-STB", "D8_Euploidy_pre-STB"),
    c("D8_M22_pre-STB", "D8_Euploidy_pre-STB"),
    c("D8_M16_pre-STB", "D8_Euploidy_pre-STB"),
    c("D8_M16?_pre-STB", "D8_Euploidy_pre-STB")
  )
)
p3 <- FeatureStatPlot(
  srt = tmp[,grep("^CTB$", tmp$CellType)], group.by = "MergedCellType",
  stat.by = c("Pseudotime"), add_box = TRUE,
  comparisons = list(
    c("D10_T16_CTB", "D10_Euploidy_CTB"),
    c("D8_T16_CTB", "D8_Euploidy_CTB"),
    c("D8_M22_CTB", "D8_Euploidy_CTB"),
    c("D8_M16_CTB", "D8_Euploidy_CTB"),
    c("D8_M16?_CTB", "D8_Euploidy_CTB")
  )
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_pseudotime_plot_TE_lineage.pdf"), height = 6, width = 12)
p1 + p2 + p3
dev.off()
# ICM-lineage
tmp <- sr.merge[,colnames(sr.merge) %in% test$id$`ICM-lineage`]
mn3$ICM <- TI.Mncl3(sr.ob = tmp, btc.column = "DataSet", gp.column = "CellType",
                    mt.column = "percent.mt", use.sr.umap = T,
                    cluster.res = NULL, root.gp = "PreEPI", sp.name = "ICM_lineage",
                    pt.size = 1, out.dir = file.path(res.out, "Monocle3"))
mn3$ICM.gene <- GeneModule.Mncl3(mncl3.cds = mn3$ICM, sp.name = "ICM_lineage_gene_module", 
                                 gp.column = "CellType", ncores = 16, qvalue = 0.05, 
                                 out.dir = file.path(res.out, "Monocle3"))
if (all(names(pseudotime(mn3$ICM)) == colnames(tmp))) {
  tmp$Pseudotime <- log2(pseudotime(mn3$ICM)+1)
}
tmp$MergedCellType <- paste(tmp$Stage, tmp$Karyotype, tmp$CellType, sep = "_")
tmp <- subset(tmp, DataSet == "ThisStudy")
p1 <- FeatureStatPlot(
  srt = tmp[,grep("EPI", tmp$CellType)], group.by = "MergedCellType",
  stat.by = c("Pseudotime"), add_box = TRUE,
  comparisons = list(
    c("D10_T16_EPI", "D10_Euploidy_EPI"),
    c("D8_T16_EPI", "D8_Euploidy_EPI"),
    c("D8_M22_EPI", "D8_Euploidy_EPI"),
    c("D8_M16_EPI", "D8_Euploidy_EPI"),
    c("D8_M16?_EPI", "D8_Euploidy_EPI")
  )
)
pdf(file.path(res.out, "XuYanWen_Aneuploid_no_chr16_chr22_integrated_pseudotime_plot_ICM_lineage.pdf"), height = 6, width = 8)
p1
dev.off()


### >>> 5. Differentiation state
cyto <- list()
# TE-lineage
tmp <- sr.merge[,colnames(sr.merge) %in% test$id$`TE-lineage`]
tmp$Karyotype <- factor(tmp$Karyotype)
cyto$TE <- Pipe.CytoTrace(sr.obj = tmp, sample.n = NULL, 
                          group.by = "Karyotype", 
                          gene.num = 10, 
                          geneset = markers, 
                          prefix = "TE-lineage", 
                          outdir = file.path(res.out, "CytoTrace/TE-lineage"))
plotCytoGenes(cyto$TE, numOfGenes = 10, outputDir = file.path(res.out, "CytoTrace/TE-lineage/"))
# ICM-lineage
tmp <- sr.merge[,colnames(sr.merge) %in% test$id$`ICM-lineage`]
tmp$Karyotype <- factor(tmp$Karyotype)
cyto$TE <- Pipe.CytoTrace(sr.obj = tmp, sample.n = NULL, 
                          group.by = "Karyotype", 
                          gene.num = 10, 
                          geneset = markers, 
                          prefix = "ICM-lineage", 
                          outdir = file.path(res.out, "CytoTrace/ICM-lineage"))
plotCytoGenes(cyto$TE, numOfGenes = 10, outputDir = file.path(res.out, "CytoTrace/ICM-lineage/"))


### >>> 6. Differential expression analysis
table(sr.yhw.no22$CellType, sr.yhw.no22$Karyotype, sr.yhw.no22$Stage)
table(sr.merge$Stage, sr.merge$DataSet, sr.merge$CellType)
table(sr.merge$CellType, sr.merge$DataSet)
sr.merge.deg <- sr.merge
sr.merge.deg@meta.data <- sr.merge.deg@meta.data %>% 
  mutate(CellType = gsub("PostEPI.E1", "EPI", gsub("PrE", "Hypoblast", gsub("Early.STB", "pre-STB", CellType))),
         Karyotype = gsub("\\?", "x", Karyotype))
sr.merge.deg <- sr.merge.deg[, sr.merge.deg$CellType %in% c("EPI", "Hypoblast", "CTB", "pre-STB", "STB")]
table(sr.merge.deg$CellType, sr.merge.deg$DataSet)
sr.merge.deg <- sr.merge.deg[, !(sr.merge.deg$Karyotype %in% c("Euploidy"))]
sr.merge.deg <- sr.merge.deg[, sr.merge.deg$Stage %in% c("D8", "D10")]
table(sr.merge.deg$Karyotype, sr.merge.deg$CellType, sr.merge.deg$Stage)
table(sr.merge.deg$Karyotype)
deg.vs.li <- list()
go.vs.li <- list()
# EPI D8: M16? M22 T16 vs GSE136447-PostEPI.E1
for (i in c("M16x", "M22", "T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "EPI" & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("EPI_D8_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                   sample.n = NULL, group.by = "CellType",
                                                   g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                   res.out = file.path(res.out, paste0("DEG/EPI_D8/", i, "_vs_Euploidy")))
}
# EPI D10: T16 vs GSE136447-PostEPI.E1
for (i in c("T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "EPI" & Stage == "D10")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("EPI_D10_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                    sample.n = NULL, group.by = "CellType",
                                                    g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                    res.out = file.path(res.out, paste0("DEG/EPI_D10/", i, "_vs_Euploidy")))
}
# Hypoblast D8: M16? M22 T16 vs GSE136447-PrE
for (i in c("M16x", "M22", "T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "Hypoblast" & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("Hypoblast_D8_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                         sample.n = NULL, group.by = "CellType",
                                                         g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                         res.out = file.path(res.out, paste0("DEG/Hypoblast_D8/", i, "_vs_Euploidy")))
}

# CTB D8: M16 M16? M22 T16 vs GSE136447-CTB
for (i in c("M16", "M16x", "M22", "T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "CTB" & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("CTB_D8_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                   sample.n = NULL, group.by = "CellType",
                                                   g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                   res.out = file.path(res.out, paste0("DEG/CTB_D8/", i, "_vs_Euploidy")))
}
# CTB D10: T16 vs GSE136447-CTB
for (i in c("T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "CTB" & Stage == "D10")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("CTB_D10_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                    sample.n = NULL, group.by = "CellType",
                                                    g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                    res.out = file.path(res.out, paste0("DEG/CTB_D10/", i, "_vs_Euploidy")))
}

# pre-STB D8: T16 vs GSE136447-STB
for (i in c("M16", "M16x", "M22", "T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "pre-STB" & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("pre-STB_D10_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                        sample.n = NULL, group.by = "CellType",
                                                        g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                        res.out = file.path(res.out, paste0("DEG/pre-STB_D8/", i, "_vs_Euploidy")))
}
# pre-STB D10: T16 vs GSE136447-STB
for (i in c("T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "pre-STB" & Stage == "D10")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("pre-STB_D10_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                        sample.n = NULL, group.by = "CellType",
                                                        g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                        res.out = file.path(res.out, paste0("DEG/pre-STB_D10/", i, "_vs_Euploidy")))
}

# STB D8: T16 vs GSE136447-STB
for (i in c("M16x", "T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "STB" & Stage == "D8")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("STB_D10_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                    sample.n = NULL, group.by = "CellType",
                                                    g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                    res.out = file.path(res.out, paste0("DEG/STB_D8/", i, "_vs_Euploidy")))
}
# STB D10: T16 vs GSE136447-STB
for (i in c("T16")) {
  tmp.sr <- subset(sr.merge.deg, CellType == "STB" & Stage == "D10")
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype", comparison = c("GSE136447", i), sample.n = NULL)
  deg.vs.li[[paste0("STB_D10_", i)]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                                    sample.n = NULL, group.by = "CellType",
                                                    g1 = "1_GSE136447", g2 = paste0("2_", i), lfc = 1, sig = 0.05,
                                                    res.out = file.path(res.out, paste0("DEG/STB_D10/", i, "_vs_Euploidy")))
}
# Volcano
for (i in names(deg.vs.li)) {
  pdf(file.path(res.out, paste0("DEG/", i, "_vs_Euploidy_volcano.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.vs.li[[i]]$all, geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.01, lfc = 2, title = i,
                       pt.size = 2, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
  pdf(file.path(res.out, paste0("DEG/", i, "_vs_Euploidy_volcano_TFS.pdf")), height = 6, width = 9)
  print(VisDEG.volcano(deg.data = deg.vs.li[[i]]$all, geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                       sig = 0.05, lfc = 1, title = i,
                       pt.size = 3, pt.shape = 16,
                       label.gene = NULL, gene.size = 3, num.size = 5,
                       up.col = c("#D9212A"), down.col = c("#045EC3"),
                       nosig.col = "#B2B2B2"))
  dev.off()
}
# GO
for (i in names(deg.vs.li)) {
  go.vs.li[[paste0(i, ".up")]] <- Pipe.GO(species = "human", 
                                          genelist = subset(deg.vs.li[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL,
                                          basename = paste0(i, "_vs_Euploidy_up"), genetype = "SYMBOL",
                                          res.out = file.path(res.out, paste0("DEG/GO/", i)))
  go.vs.li[[paste0(i, ".down")]] <- Pipe.GO(species = "human", 
                                            genelist = subset(deg.vs.li[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL,
                                            basename = paste0(i, "_vs_Euploidy_down"), genetype = "SYMBOL",
                                            res.out = file.path(res.out, paste0("DEG/GO/", i)))
}
# GO chr16
for (i in names(deg.vs.li)) {
  go.vs.li[[paste0(i, ".chr16.up")]] <- Pipe.GO(species = "human", 
                                                genelist = intersect(chr16.gene$V1,
                                                                     subset(deg.vs.li[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                basename = paste0(i, "_vs_Euploidy_chr16_up"), genetype = "SYMBOL",
                                                res.out = file.path(res.out, paste0("DEG/GO_chr16/", i)))
  go.vs.li[[paste0(i, ".chr16.down")]] <- Pipe.GO(species = "human", 
                                                  genelist = intersect(chr16.gene$V1,
                                                                       subset(deg.vs.li[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                  basename = paste0(i, "_vs_Euploidy_chr16_down"), genetype = "SYMBOL",
                                                  res.out = file.path(res.out, paste0("DEG/GO_chr16/", i)))
}
# GO chr22
for (i in names(deg.vs.li)) {
  go.vs.li[[paste0(i, ".chr22.up")]] <- Pipe.GO(species = "human", 
                                                genelist = intersect(chr22.gene$V1,
                                                                     subset(deg.vs.li[[i]]$sig, PValue <= 0.05 & logFC >= log2(2))$SYMBOL),
                                                basename = paste0(i, "_vs_Euploidy_chr22_up"), genetype = "SYMBOL",
                                                res.out = file.path(res.out, paste0("DEG/GO_chr22/", i)))
  go.vs.li[[paste0(i, ".chr22.down")]] <- Pipe.GO(species = "human", 
                                                  genelist = intersect(chr22.gene$V1,
                                                                       subset(deg.vs.li[[i]]$sig, PValue <= 0.05 & logFC <= -log2(2))$SYMBOL),
                                                  basename = paste0(i, "_vs_Euploidy_chr22_down"), genetype = "SYMBOL",
                                                  res.out = file.path(res.out, paste0("DEG/GO_chr22/", i)))
}



# =======================================
# 8th part: Alternative splicing analysis ----
# =======================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/Alternative_Splicing")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Load data
# splice junction counts
sj <- readRDS("analysis/scrna/results/marvel/star/Splice_Junction_Counts.rds")
sj[!is.na(sj[,2]), ][1:5, 1:5]
colnames(sj)[-1] <- gsub("-", "_", gsub("\\.", "_", gsub("d", "D", gsub("_R1.*", "", colnames(sj)[-1]))))
if (all(colnames(sj)[-1] %in% colnames(sr.yhw.no22))) {
  sj <- sj[, c("coord.intron", colnames(sr.yhw.no22))]
}
all(colnames(sj)[-1] == colnames(sr.yhw.no22))
# cell meta
df.pheno <- sr.yhw.no22@meta.data
all(colnames(sj)[-1] == rownames(df.pheno))
. <- data.frame("sample.id" = rownames(df.pheno), stringsAsFactors = FALSE)
df.pheno <- cbind.data.frame(., df.pheno)
row.names(df.pheno) <- NULL
all(colnames(sj)[-1] == df.pheno$sample.id)
# feature.list
df.feature.list <- readRDS("analysis/scrna/results/marvel/rMATs/df.feature.list.rds")
df.feature.list
# intron counts
df.intron.counts <- as.data.frame(fread("analysis/scrna/results/marvel/intron/Counts_by_Region.txt", 
                                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA"))


colnames(df.intron.counts)[-1] <- gsub("-", "_", gsub("\\.", "_", gsub("d", "D", gsub("_R1.*", "", colnames(df.intron.counts)[-1]))))
if (all(colnames(df.intron.counts)[-1] %in% colnames(sr.yhw.no22))) {
  df.intron.counts <- df.intron.counts[, c("coord.intron", colnames(sr.yhw.no22))]
}
all(colnames(df.intron.counts)[-1] == colnames(sr.yhw.no22))
# gene feature
df.tpm.feature <- read.table("/home/laborer/refgenome/embl/homo_sapiens/anno/Homo_sapiens.GRCh38.110.chr.no.contig_gene_anno.txt",
                             sep = "\t")
df.tpm.feature <- df.tpm.feature[, -1:-4]
colnames(df.tpm.feature) <- c("gene_id", "gene_short_name", "gene_type")
df.tpm.feature <- df.tpm.feature[!duplicated(df.tpm.feature$gene_short_name), ]
# gene tpm
df.tpm <- CountToTpm(ge.count[, -1:-5], ge.count$Length)
all(colnames(df.tpm) == rownames(df.pheno))
keep.gene <- unique(intersect(df.tpm.feature$gene_short_name, rownames(df.tpm)))
df.tpm.feature <- df.tpm.feature[df.tpm.feature$gene_short_name %in% keep.gene, ]
if (all(df.tpm.feature$gene_short_name %in% rownames(df.tpm))) {
  df.tpm <- df.tpm[df.tpm.feature$gene_short_name, ]
}
if (all(rownames(df.tpm) == df.tpm.feature$gene_short_name)) {
  . <- data.frame("gene_id" = df.tpm.feature$gene_id, stringsAsFactors = FALSE)
  df.tpm <- cbind.data.frame(., df.tpm)
  row.names(df.tpm) <- NULL
}
# gtf
gtf <- as.data.frame(fread("analysis/scrna/metadata/gene_chr.gtf", 
                           sep = "\t", header = FALSE, stringsAsFactors = FALSE, na.strings = "NA", quote = "\""))


### >>> 3. Create MARVEL object
library("MARVEL")
marvel <- CreateMarvelObject(SpliceJunction = sj,
                             SplicePheno = df.pheno,
                             SpliceFeature = df.feature.list,
                             IntronCounts = df.intron.counts,
                             GeneFeature = df.tpm.feature,
                             Exp = df.tpm,
                             GTF = gtf)


### >>> 4. Detect additional events
# AFE
marvel <- DetectEvents(MarvelObject = marvel,
                       min.cells = 15,
                       min.expr = 1,
                       track.progress = FALSE,
                       EventType = "AFE"
)
# ALE
marvel <- DetectEvents(MarvelObject = marvel,
                       min.cells = 15,
                       min.expr = 1,
                       track.progress = FALSE,
                       EventType = "ALE"
)


### >>> 5. Check splicing junction data
marvel <- CheckAlignment(MarvelObject = marvel, level = "SJ")


### >>> 6. Validate, filter, compute splicing events
# SE
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="SE"
)
# MXE
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     UnevenCoverageMultiplier=10,
                     EventType="MXE"
)
# RI
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="RI",
                     thread=4
)
# A5SS
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A5SS"
)
# A3SS
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="A3SS"
)
# AFE
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="AFE"
)
# ALE
marvel <- ComputePSI(MarvelObject=marvel,
                     CoverageThreshold=10,
                     EventType="ALE"
)


### >>> 7. Transform expression values
marvel <- TransformExpValues(MarvelObject=marvel,
                             offset=1,
                             transformation="log2",
                             threshold.lower=1
)


### >>> 8. Check matrices and metadata
# Check splicing data
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing")
# Check gene data
marvel <- CheckAlignment(MarvelObject=marvel, level="gene")
# Cross-check splicing and gene data
marvel <- CheckAlignment(MarvelObject=marvel, level="splicing and gene")


### >>> 9. Overview of splicing events
if (! dir.exists(file.path(res.out, "Splicing_Event"))) { 
  dir.create(file.path(res.out, "Splicing_Event"), recursive = T)
}
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
table(df.pheno$Karyotype, df.pheno$CellType, df.pheno$Stage)
for (i in unique(df.pheno$Stage)) {
  for (j in unique(df.pheno$CellType)) {
    for (k in unique(df.pheno$Karyotype)) {
      # Define sample ids
      sample.ids <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == k), 
                             "sample.id"]
      if (length(sample.ids) > 0) {
        # Tabulate expressed events
        tmp <- CountEvents(MarvelObject=marvel,
                           sample.ids=sample.ids,
                           min.cells=3
        )
        # Output Plot
        pdf(file.path(res.out, paste0("Splicing_Event/", i, "_", j, "_", k, ".pdf")), height = 6, width = 7)
        print(tmp$N.Events$Plot)
        dev.off()
        write.table(tmp$N.Events$Table, file.path(res.out, paste0("Splicing_Event/", i, "_", j, "_", k, ".txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      }
    }
  }
}


### >>> 10. Modality analysis
if (! dir.exists(file.path(res.out, "Modality_Analysis"))) { 
  dir.create(file.path(res.out, "Modality_Analysis"), recursive = T)
}
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
table(df.pheno$Karyotype, df.pheno$CellType, df.pheno$Stage)
for (i in unique(df.pheno$Stage)) {
  for (j in unique(df.pheno$CellType)) {
    for (k in unique(df.pheno$Karyotype)) {
      # Define sample ids
      sample.ids <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == k), 
                             "sample.id"]
      if (length(sample.ids) >= 3) {
        # Assign modality
        tmp <- AssignModality(MarvelObject=marvel,
                              sample.ids=sample.ids,
                              min.cells=3,
                              seed=1
        )
        tmp$Modality$Results[1:5, c("tran_id", "event_type", "gene_id", "gene_short_name", "modality.bimodal.adj")]
        # Tabulate modality proportion (overall)
        tmp <- PropModality(MarvelObject=tmp,
                            modality.column="modality.bimodal.adj",
                            modality.type="extended",
                            event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                            across.event.type=FALSE
        )
        # Output Plot
        pdf(file.path(res.out, paste0("Modality_Analysis/", i, "_", j, "_", k, "_donut.pdf")), height = 6, width = 7)
        print(tmp$Modality$Prop$DoughnutChart$Plot)
        dev.off()
        write.table(tmp$Modality$Prop$DoughnutChart$Table, 
                    file.path(res.out, paste0("Modality_Analysis/", i, "_", j, "_", k, "_donut.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        # Tabulate modality proportion (by event type)
        tmp <- PropModality(MarvelObject=tmp,
                            modality.column="modality.bimodal.adj",
                            modality.type="extended",
                            event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                            across.event.type=TRUE,
                            prop.test="chisq",
                            prop.adj="fdr",
                            xlabels.size=8
        )
        # Output Plot
        pdf(file.path(res.out, paste0("Modality_Analysis/", i, "_", j, "_", k, "_barplot.pdf")), height = 6, width = 7)
        print(tmp$Modality$Prop$BarChart$Plot)
        dev.off()
        write.table(tmp$Modality$Prop$BarChart$Table, 
                    file.path(res.out, paste0("Modality_Analysis/", i, "_", j, "_", k, "_barplot.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      }
    }
  }
}


### >>> 11. Differential splicing analysis
if (! dir.exists(file.path(res.out, "Differential_Splicing_Gene"))) { 
  dir.create(file.path(res.out, "Differential_Splicing_Gene"), recursive = T)
}
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
table(df.pheno$Karyotype, df.pheno$CellType, df.pheno$Stage)
for (g in c("T16", "M16", "M22", "M16?")) {
  for (i in unique(df.pheno$Stage)) {
    for (j in unique(df.pheno$CellType)) {
      # Cell group 1 (reference)
      cell.group.g1 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == "Euploidy"), "sample.id"]
      # Cell group 2
      cell.group.g2 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == g), "sample.id"]
      if (length(cell.group.g1) >= 2 & length(cell.group.g2) >= 2) {
        # Differential gene expression analysis
        tmp <- CompareValues(MarvelObject=marvel,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             min.cells=3,
                             method="wilcox",
                             method.adjust="fdr",
                             level="gene",
                             show.progress=FALSE
        )
        results <- tmp$DE$Exp$Table
        write.table(results, file.path(res.out, paste0("Differential_Splicing_Gene/", i, "_", j, "_", g, "_vs_Euploidy_DEGs.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        pd.label <- rbind(
          results %>% filter(p.val <= 0.05) %>% arrange(desc(log2fc)) %>% head(n = 5),
          results %>% filter(p.val <= 0.05) %>% arrange(log2fc) %>% head(n = 5),
          results %>% filter(p.val <= 0.05) %>% arrange(p.val) %>% head(n = 10)
        ) %>% unique()
        gene_short_names <- pd.label[, "gene_short_name"]
        tmp <- PlotDEValues(MarvelObject=tmp,
                            pval=0.10,
                            log2fc=0.5,
                            point.size=0.1,
                            xlabel.size=10,
                            level="gene.global",
                            anno=TRUE,
                            anno.gene_short_name=gene_short_names
        )
        pdf(file.path(res.out, paste0("Differential_Splicing_Gene/", i, "_", j, "_", g, "_vs_Euploidy_DEGs.pdf")), height = 6, width = 7)
        print(tmp$DE$Exp.Global$Plot)
        dev.off()
        # Differential splicing analysis
        tmp <- CompareValues(MarvelObject=tmp,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             min.cells=2,
                             method=c("ad", "dts"),
                             method.adjust="fdr",
                             level="splicing",
                             event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE"),
                             show.progress=FALSE
        )
        write.table(tmp$DE$PSI$Table[["ad"]], 
                    file.path(res.out, paste0("Differential_Splicing_Gene/", i, "_", j, "_", g, "_vs_Euploidy_DS_ad.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        write.table(tmp$DE$PSI$Table[["dts"]], 
                    file.path(res.out, paste0("Differential_Splicing_Gene/", i, "_", j, "_", g, "_vs_Euploidy_DS_dts.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        # Differential (spliced) gene analysis
        tmp <- CompareValues(MarvelObject=tmp,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             psi.method=c("ad", "dts"),
                             psi.pval=c(0.10, 0.10),
                             psi.delta=0,
                             method.de.gene="wilcox",
                             method.adjust.de.gene="fdr",
                             downsample=FALSE,
                             show.progress=FALSE,
                             level="gene.spliced"
        )
        results <- tmp$DE$Exp.Spliced$Table
        write.table(results, file.path(res.out, paste0("Differential_Splicing_Gene/", i, "_", j, "_", g, "_vs_Euploidy_DSGs.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        pd.label <- rbind(
          results %>% filter(p.val <= 0.05) %>% arrange(desc(log2fc)) %>% head(n = 5),
          results %>% filter(p.val <= 0.05) %>% arrange(log2fc) %>% head(n = 5),
          results %>% filter(p.val <= 0.05) %>% arrange(p.val) %>% head(n = 10)
        ) %>% unique()
        gene_short_names <- pd.label[, "gene_short_name"]
        tmp <- PlotDEValues(MarvelObject=tmp,
                            method=c("ad", "dts"),
                            psi.pval=c(0.10, 0.10),
                            psi.delta=0,
                            gene.pval=0.10,
                            gene.log2fc=0.5,
                            point.size=0.1,
                            xlabel.size=8,
                            level="gene.spliced",
                            anno=TRUE,
                            anno.gene_short_name=gene_short_names
        )
        pdf(file.path(res.out, paste0("Differential_Splicing_Gene/", i, "_", j, "_", g, "_vs_Euploidy_DSGs.pdf")), height = 6, width = 7)
        print(tmp$DE$Exp.Spliced$Plot)
        dev.off()
      }
    }
  }
}


### >>> 12. Modality dynamics
if (! dir.exists(file.path(res.out, "Modality_Dynamics"))) { 
  dir.create(file.path(res.out, "Modality_Dynamics"), recursive = T)
}
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
table(df.pheno$Karyotype, df.pheno$CellType, df.pheno$Stage)
for (g in c("T16", "M16", "M22", "M16?")) {
  for (i in unique(df.pheno$Stage)) {
    for (j in unique(df.pheno$CellType)) {
      # Cell group 1 (reference)
      cell.group.g1 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == "Euploidy"), "sample.id"]
      # Cell group 2
      cell.group.g2 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == g), "sample.id"]
      # Merge
      cell.group.list <- list("Euploidy" = cell.group.g1,
                              g = cell.group.g2)
      names(cell.group.list)[2] <- g
      if (length(cell.group.g1) >= 2 & length(cell.group.g2) >= 2) {
        # Assign modality dynamics
        tmp <- ModalityChange(MarvelObject=tmp,
                              method=c("ad", "dts"),
                              psi.pval=c(0.10, 0.10)
        )
        pdf(file.path(res.out, paste0("Modality_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Modality_plot.pdf")), height = 6, width = 7)
        print(tmp$DE$Modality$Plot)
        dev.off()
        write.table(tmp$DE$Modality$Table, 
                    file.path(res.out, paste0("Modality_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Modality_table.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        write.table(tmp$DE$Modality$Plot.Stats, 
                    file.path(res.out, paste0("Modality_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Modality_stats.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        # Plot
        for (l in 1:nrow(tmp$DE$Modality$Table)) {
          # gene
          tmp <- PlotValues(MarvelObject=tmp,
                            cell.group.list=cell.group.list,
                            feature=tmp$DE$Modality$Table$gene_id[l],
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
          )
          plot.1_gene <- tmp$adhocPlot$Exp
          # splicing
          tmp <- PlotValues(MarvelObject=tmp,
                            cell.group.list=cell.group.list,
                            feature=tmp$DE$Modality$Table$tran_id[l],
                            xlabels.size=5,
                            level="splicing",
                            min.cells=2
          )
          plot.1_splicing <- tmp$adhocPlot$PSI
          file.name <- paste(tmp$DE$Modality$Table$gene_short_name[l],
                             tmp$DE$Modality$Table$event_type[l],
                             tmp$DE$Modality$Table$modality.change[l],
                             tmp$DE$Modality$Table$tran_id[l], sep = "_")
          pdf(file.path(res.out, paste0("Modality_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Modality_PSI_plot_", file.name, ".pdf")), 
              height = 6, width = 7)
          print(plot.1_gene + plot.1_splicing)
          dev.off()
        }
      }
    }
  }
}


### >>> 13. Gene-splicing dynamics
# - Coordinated gene-splicing relationship refers to the change in mean gene expression 
#   is in the same direction with the corresponding splicing event(s).

# - Opposing gene-splicing relationship refers to the change in mean gene expression 
#   is in the opposite direction to the corresponding splicing event(s).

# - Isoform-switching refers to genes that are differentially spliced 
#   without being differentially expressed.

# - Complex gene-splicing relationship refers to genes with both coordinated and 
#   opposing relationships with the corresponding splicing events.
if (! dir.exists(file.path(res.out, "Gene-splicing_Dynamics"))) { 
  dir.create(file.path(res.out, "Gene-splicing_Dynamics"), recursive = T)
}
# Retrieve sample metadata
df.pheno <- marvel$SplicePheno
table(df.pheno$Karyotype, df.pheno$CellType, df.pheno$Stage)
for (g in c("T16", "M16", "M22", "M16?")) {
  for (i in unique(df.pheno$Stage)) {
    for (j in unique(df.pheno$CellType)) {
      # Cell group 1 (reference)
      cell.group.g1 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == "Euploidy"), "sample.id"]
      # Cell group 2
      cell.group.g2 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == g), "sample.id"]
      # Merge
      cell.group.list <- list("Euploidy" = cell.group.g1,
                              g = cell.group.g2)
      names(cell.group.list)[2] <- g
      if (length(cell.group.g1) >= 2 & length(cell.group.g2) >= 2) {
        # Assign modality dynamics
        tmp <- ModalityChange(MarvelObject=tmp,
                              method=c("ad", "dts"),
                              psi.pval=c(0.10, 0.10)
        )
        tmp <- IsoSwitch(MarvelObject=tmp,
                         method=c("ad", "dts"),
                         psi.pval=c(0.10, 0.10),
                         psi.delta=0,
                         gene.pval=0.10,
                         gene.log2fc=0.5
        )
        pdf(file.path(res.out, paste0("Gene-splicing_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Splicing_plot.pdf")), height = 6, width = 7)
        print(tmp$DE$Cor$Plot)
        dev.off()
        write.table(tmp$DE$Cor$Table, 
                    file.path(res.out, paste0("Gene-splicing_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Splicing_table.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        write.table(tmp$DE$Cor$Plot.Stats, 
                    file.path(res.out, paste0("Gene-splicing_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Splicing_stats.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        # Plot
        for (type in unique(tmp$DE$Cor$Table$cor)) {
          results <- subset(tmp$DE$Cor$Table, cor == type)
          for (l in 1:nrow(results)) {
            # gene
            tmp <- PlotValues(MarvelObject=tmp,
                              cell.group.list=cell.group.list,
                              feature=results$gene_id[l],
                              maintitle="gene_short_name",
                              xlabels.size=7,
                              level="gene"
            )
            plot.1_gene <- tmp$adhocPlot$Exp
            # splicing
            tmp <- PlotValues(MarvelObject=tmp,
                              cell.group.list=cell.group.list,
                              feature=subset(tmp$DE$Modality$Table, gene_id == results$gene_id[l])$tran_id,
                              xlabels.size=5,
                              level="splicing",
                              min.cells=2
            )
            plot.1_splicing <- tmp$adhocPlot$PSI
            file.name <- paste(type,
                               tmp$DE$Modality$Table$gene_short_name[l],
                               tmp$DE$Modality$Table$event_type[l],
                               tmp$DE$Modality$Table$modality.change[l],
                               tmp$DE$Modality$Table$tran_id[l], sep = "_")
            pdf(file.path(res.out, paste0("Gene-splicing_Dynamics/", i, "_", j, "_", g, "_vs_Euploidy_Modality_PSI_plot_", file.name, ".pdf")), 
                height = 6, width = 7)
            print(plot.1_gene + plot.1_splicing)
            dev.off()
          }
        }
      }
    }
  }
}


### >>> 14. NMD analysis
if (! dir.exists(file.path(res.out, "NMD_Analysis"))) { 
  dir.create(file.path(res.out, "NMD_Analysis"), recursive = T)
}
# Parse GTF
tmp <- ParseGTF(MarvelObject=tmp)
df.pheno <- marvel$SplicePheno
table(df.pheno$Karyotype, df.pheno$CellType, df.pheno$Stage)
for (g in c("T16", "M16", "M22", "M16?")) {
  for (i in unique(df.pheno$Stage)) {
    for (j in unique(df.pheno$CellType)) {
      # Cell group 1 (reference)
      cell.group.g1 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == "Euploidy"), "sample.id"]
      # Cell group 2
      cell.group.g2 <- df.pheno[which(df.pheno$Stage == i & df.pheno$CellType == j & df.pheno$Karyotype == g), "sample.id"]
      if (length(cell.group.g1) >= 2 & length(cell.group.g2) >= 2) {
        # Differential splicing analysis
        tmp <- CompareValues(MarvelObject=tmp,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             min.cells=2,
                             method=c("ad", "dts"),
                             method.adjust="fdr",
                             level="splicing",
                             event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE"),
                             show.progress=FALSE)
        # Find PTCs
        tmp <- FindPTC(MarvelObject=tmp,
                       method=c("ad", "dts"),
                       pval=c(0.10, 0.10),
                       delta=5)
        # Show novel SJ and non-CDS transcripts
        tmp <- PropPTC(MarvelObject=tmp,
                       xlabels.size=8,
                       show.NovelSJ.NoCDS=TRUE,
                       prop.test="chisq")
        pdf(file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_PTC_Prop.pdf")), 
            height = 6, width = 7)
        print(tmp$NMD$PTC.Prop$Plot)
        dev.off()
        write.table(tmp$NMD$PTC.Prop$Table, 
                    file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_PTC_Prop_table.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        write.table(tmp$NMD$PTC.Prop$Plot.Stats, 
                    file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_PTC_Prop_stats.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        # NMD prediction
        tmp <- CompareExpr(MarvelObject=tmp,
                           xlabels.size=8)
        pdf(file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_NMD_Expr.pdf")), 
            height = 6, width = 7)
        print(tmp$NMD$NMD.Expr$Plot)
        dev.off()
        write.table(tmp$NMD$NMD.Expr$Table, 
                    file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_NMD_Expr_table.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        write.table(tmp$NMD$NMD.Expr$Plot.Stats, 
                    file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_NMD_Expr_stats.txt")),
                    col.names = T, row.names = F, sep = "\t", quote = F)
        if (sum(tmp$NMD$PTC.Prop$Table$NMD == "PTC") > 0) {
          # NMD gene candidates
          tmp <- AnnoVolcanoPlot(MarvelObject=tmp,
                                 anno=FALSE,
                                 xlabel.size=7)
          # Volcano plot: Annotate candidate NMD genes
          results <- tmp$NMD$AnnoVolcanoPlot$Table
          write.table(results, file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_NMD_Candidates.txt")),
                      col.names = T, row.names = F, sep = "\t", quote = F)
          pd.label <- rbind(
            results %>% filter(p.val <= 0.05) %>% arrange(desc(log2fc)) %>% head(n = 5),
            results %>% filter(p.val <= 0.05) %>% arrange(log2fc) %>% head(n = 5),
            results %>% filter(p.val <= 0.05) %>% arrange(p.val) %>% head(n = 10)
          ) %>% unique()
          gene_short_names <- pd.label[, "gene_short_name"]
          tmp <- AnnoVolcanoPlot(MarvelObject=tmp,
                                 anno=TRUE,
                                 anno.gene_short_name=gene_short_names,
                                 point.size=1.0,
                                 label.size=2.5,
                                 xlabel.size=7
          )
          pdf(file.path(res.out, paste0("NMD_Analysis/", i, "_", j, "_", g, "_vs_Euploidy_NMD_Candidates_Volcano.pdf")), 
              height = 6, width = 7)
          print(tmp$NMD$AnnoVolcanoPlot$Plot)
          dev.off()
        }
      }
    }
  }
}



# =====================================
# 9th part: LinGe aneuploid blastocysts ----
# =====================================

### >>> 1. Setting output directory
res.out <- file.path(getwd(), "R/Graphs/LinGe_Aneuploid")
if (! dir.exists(res.out)) { dir.create(res.out, recursive = T) }


### >>> 2. Load data
sr.linge <- readRDS("/home/yhw/bioinfo/project-xuyanwen/aneuploid/R/CodeData/embryo.14908cells.rawcounts.rds")
sr.linge <- sr.linge[!duplicated(gsub("-ENSG.*", "", rownames(sr.linge))), ]
sr.linge[1:10, 1:10]
rownames(sr.linge) <- gsub("-ENSG.*", "", rownames(sr.linge))
meta <- read.csv("/home/yhw/bioinfo/project-xuyanwen/aneuploid/R/CodeData/Metadata.csv", header = T, row.names = 1)
if (all(rownames(meta) %in% colnames(sr.linge)) & all(colnames(sr.linge) %in% rownames(meta))) {
  meta <- meta[colnames(sr.linge), ]
}
if (all(rownames(meta) == colnames(sr.linge))) {
  sr.linge <- CreateSeuratObject(counts = sr.linge, meta.data = meta)
}
rm(meta)


### >>> 3. Differential expression analysis
table(sr.linge$Karyotype.By.PGT.A)
# DEGs between karyotypes (without cell types)
deg.linge <- list()
deg.meta <- data.frame(g1 = c("E", "E", "E"), 
                       g2 = c("16T", "16M", "22M"))
for (i in 1:nrow(deg.meta)) {
  tmp.sr <- subset(sr.linge, Karyotype.By.PGT.A %in% deg.meta[i,])
  tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype.By.PGT.A", 
                       comparison = c(deg.meta[i, "g1"], deg.meta[i, "g2"]), 
                       sample.n = 200)
  tmp.name <- paste0(deg.meta[i, "g2"], "_vs_", deg.meta[i, "g1"])
  deg.linge[[tmp.name]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                       sample.n = NULL, group.by = "CellType", 
                                       g1 = paste0("1_", deg.meta[i, "g1"]), 
                                       g2 = paste0("2_", deg.meta[i, "g2"]), 
                                       lfc = 1, sig = 0.05,
                                       res.out = file.path(res.out, paste0("DEG/", tmp.name)))
}
# DEGs between karyotypes (with cell types)
for (j in unique(sr.linge$Lineage)) {
  for (i in 1:nrow(deg.meta)) {
    tmp.sr <- subset(sr.linge, Karyotype.By.PGT.A %in% deg.meta[i,] & Lineage == j)
    tmp <- SeuratToEdger(sr.obj = tmp.sr, group.by = "Karyotype.By.PGT.A", 
                         comparison = c(deg.meta[i, "g1"], deg.meta[i, "g2"]), 
                         sample.n = 300)
    tmp.name <- paste0(deg.meta[i, "g2"], "_vs_", deg.meta[i, "g1"], "_", j)
    deg.linge[[tmp.name]] <- edgeR.scRNA(count = tmp$count, meta = tmp$meta, expr = tmp$expr,
                                         sample.n = NULL, group.by = "CellType", 
                                         g1 = paste0("1_", deg.meta[i, "g1"]), 
                                         g2 = paste0("2_", deg.meta[i, "g2"]), 
                                         lfc = 1, sig = 0.05,
                                         res.out = file.path(res.out, paste0("DEG/", tmp.name)))
  }
}
# draw volcano plots
for (j in c(0, 0.001, 0.005, 0.01, 0.05, 0.1)) {
  for (i in names(deg.linge)) {
    pdf(file.path(res.out, paste0(i, "_volcano_", j, ".pdf")), height = 6, width = 8)
    print(VisDEG.volcano(deg.data = deg.linge[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.linge[[i]]$all), value = T), low.expr = j,
                         geneset = NULL, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.01, lfc = 2, title = i,
                         pt.size = 2, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
    pdf(file.path(res.out, paste0(i, "_volcano_TFs_", j, ".pdf")), height = 6, width = 8)
    print(VisDEG.volcano(deg.data = deg.linge[[i]]$all, 
                         expr.by = grep("^[1|2]_", colnames(deg.linge[[i]]$all), value = T), low.expr = j,
                         geneset = hs.tfs$Symbol, p.col = "PValue", lfc.col = "logFC",
                         sig = 0.05, lfc = 1, title = i,
                         pt.size = 3, pt.shape = 16,
                         label.gene = NULL, gene.size = 3, num.size = 5,
                         up.col = c("#D9212A"), down.col = c("#045EC3"),
                         nosig.col = "#B2B2B2"))
    dev.off()
  }
}
# GO analysis
go.linge <- list()
for (i in names(deg.linge)) {
  deg.data <- deg.linge[[i]]$sig
  expr.by <- grep("^[1|2]_", colnames(deg.data), value = T)
  index <- rowSums(deg.data[, expr.by] >= 0.005) >= 1
  deg.data <- deg.data[index, ]
  print(nrow(subset(deg.data, FDR <= 0.005 & logFC >= log2(3))))
  print(nrow(subset(deg.data, FDR <= 0.005 & logFC <= -log2(3))))
  go.linge[[paste0(i, ".up")]] <- Pipe.GO(species = "human",
                                         genelist = subset(deg.data, FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                                         basename = paste0(i, "_up"), genetype = "SYMBOL",
                                         res.out = file.path(res.out, paste0("GO/", i)))
  go.linge[[paste0(i, ".down")]] <- Pipe.GO(species = "human",
                                           genelist = subset(deg.data, FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                                           basename = paste0(i, "_down"), genetype = "SYMBOL",
                                           res.out = file.path(res.out, paste0("GO/", i)))
}
# GSEA analysis
gsea.linge <- list()
for (i in names(deg.linge)) {
  gsea.linge[[paste0(i, "_lfc1_pvalue0.05")]] <- Pipe.GSEA(deg.obj = deg.linge[[i]]$all, deg.type = "edger",
                                                           lfc = 1, sig = 0.05,
                                                           reversed = FALSE, species = "human",
                                                           basename = paste0(i, "_lfc1_pvalue0.05"),
                                                           genetype = "SYMBOL", gene.col = "SYMBOL",
                                                           outdir = file.path(res.out, paste0("GSEA_lfc1_pvalue0.05/", i)))
  gsea.linge[[paste0(i, "_lfc2_pvalue0.005")]] <- Pipe.GSEA(deg.obj = deg.linge[[i]]$all, deg.type = "edger",
                                                            lfc = 2, sig = 0.05,
                                                            reversed = FALSE, species = "human",
                                                            basename = paste0(i, "_lfc2_pvalue0.005"),
                                                            genetype = "SYMBOL", gene.col = "SYMBOL",
                                                            outdir = file.path(res.out, paste0("GSEA_lfc2_pvalue0.005/", i)))
  gsea.linge[[paste0(i, "_all")]] <- Pipe.GSEA(deg.obj = deg.linge[[i]]$all, deg.type = "edger", 
                                               lfc = 0, sig = 1, 
                                               reversed = FALSE, species = "human", 
                                               basename = paste0(i, "_all"), 
                                               genetype = "SYMBOL", gene.col = "SYMBOL", 
                                               outdir = file.path(res.out, paste0("GSEA_all/", i)))
}
# Overlapped genes (T16)
pd.list <- list(Lin.up = subset(deg.linge$`16T_vs_E`$up, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`16T_vs_E`$down, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.kary$`D8-T16_vs_D8-Euploidy`$up, (`1_D8-Euploidy` >= 0.005 | `2_D8-T16` >= 0.005) & 
                                   FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.kary$`D8-T16_vs_D8-Euploidy`$down, (`1_D8-Euploidy` >= 0.005 | `2_D8-T16` >= 0.005) & 
                                     FDR <= 0.005 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_T16_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_of_T16_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_T16_D8_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
pd.list <- list(Lin.up = subset(deg.linge$`16T_vs_E_EPI`$up, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`16T_vs_E_EPI`$down, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.t16.d8$T16_EPI$up, (`1_Euploidy` >= 0.005 | `2_T16` >= 0.005) & 
                                   PValue <= 0.005 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.t16.d8$T16_EPI$down, (`1_Euploidy` >= 0.005 | `2_T16` >= 0.005) & 
                                     PValue <= 0.005 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_T16_EPI_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_of_T16_EPI_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_T16_EPI_D8_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
pd.list <- list(Lin.up = subset(deg.linge$`16T_vs_E_HYP`$up, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`16T_vs_E_HYP`$down, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.t16.d8$T16_Hypoblast$up, (`1_Euploidy` >= 0.005 | `2_T16` >= 0.005) & 
                                   PValue <= 0.05 & logFC >= log2(2))$SYMBOL,
                Ours.down = subset(deg.t16.d8$T16_Hypoblast$down, (`1_Euploidy` >= 0.005 | `2_T16` >= 0.005) & 
                                     PValue <= 0.05 & logFC <= -log2(2))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_T16_HYP_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_of_T16_HYP_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_T16_HYP_D8_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
# Overlapped genes (T16)
pd.list <- list(Lin.up = subset(deg.linge$`16T_vs_E`$up, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`16T_vs_E`$down, (`1_E` >= 0.005 | `2_16T` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.kary$`D10-T16_vs_D10-Euploidy`$up, (`1_D10-Euploidy` >= 0.005 | `2_D10-T16` >= 0.005) & 
                                   FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.kary$`D10-T16_vs_D10-Euploidy`$down, (`1_D10-Euploidy` >= 0.005 | `2_D10-T16` >= 0.005) & 
                                     FDR <= 0.005 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_T16_between_LinGe_and_Ours_D10", 
               file.name = paste0("Overlapped_DEGs_of_T16_between_LinGe_and_Ours_D10"), 
               res.out = file.path(res.out, "Overlapped_T16_D10_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
# Overlapped genes (M16)
pd.list <- list(Lin.up = subset(deg.linge$`16M_vs_E`$up, (`1_E` >= 0.005 | `2_16M` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`16M_vs_E`$down, (`1_E` >= 0.005 | `2_16M` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.kary$`D8-M16_vs_D8-Euploidy`$up, (`1_D8-Euploidy` >= 0.005 | `2_D8-M16` >= 0.005) & 
                                   FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.kary$`D8-M16_vs_D8-Euploidy`$down, (`1_D8-Euploidy` >= 0.005 | `2_D8-M16` >= 0.005) & 
                                     FDR <= 0.005 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_M16_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_of_M16_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_M16_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
# Overlapped genes (M22)
pd.list <- list(Lin.up = subset(deg.linge$`22M_vs_E`$up, (`1_E` >= 0.005 | `2_22M` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`22M_vs_E`$down, (`1_E` >= 0.005 | `2_22M` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.kary$`D8-M16_vs_D8-Euploidy`$up, (`1_D8-Euploidy` >= 0.005 | `2_D8-M16` >= 0.005) & 
                                   FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.kary$`D8-M16_vs_D8-Euploidy`$down, (`1_D8-Euploidy` >= 0.005 | `2_D8-M16` >= 0.005) & 
                                     FDR <= 0.005 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_M22_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_of_M22_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_M22_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
pd.list <- list(Lin.up = subset(deg.linge$`22M_vs_E_EPI`$up, (`1_E` >= 0.005 | `2_22M` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`22M_vs_E_EPI`$down, (`1_E` >= 0.005 | `2_22M` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.m22.d8$`M22_EPI`$up, (`1_Euploidy` >= 0.005 | `2_M22` >= 0.005) & 
                                   PValue <= 0.05 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.m22.d8$`M22_EPI`$down, (`1_Euploidy` >= 0.005 | `2_M22` >= 0.005) & 
                                     PValue <= 0.05 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_M22_EPI_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_EPI_of_M22_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_M22_EPI_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))
pd.list <- list(Lin.up = subset(deg.linge$`22M_vs_E_HYP`$up, (`1_E` >= 0.005 | `2_22M` >= 0.005) & 
                                  FDR <= 0.005 & logFC >= log2(3))$SYMBOL,
                Lin.down = subset(deg.linge$`22M_vs_E_HYP`$down, (`1_E` >= 0.005 | `2_22M` >= 0.005) & 
                                    FDR <= 0.005 & logFC <= -log2(3))$SYMBOL,
                Ours.up = subset(deg.m22.d8$`M22_Hypoblast`$up, (`1_Euploidy` >= 0.005 | `2_M22` >= 0.005) & 
                                   PValue <= 0.05 & logFC >= log2(3))$SYMBOL,
                Ours.down = subset(deg.m22.d8$`M22_Hypoblast`$down, (`1_Euploidy` >= 0.005 | `2_M22` >= 0.005) & 
                                     PValue <= 0.05 & logFC <= -log2(3))$SYMBOL)
PlotOverlapped(pd.list = pd.list, pd.label = names(pd.list), pd.title = "DEGs_of_M22_Hypoblast_between_LinGe_and_Ours_D8", 
               file.name = paste0("Overlapped_DEGs_Hypoblast_of_M22_between_LinGe_and_Ours_D8"), 
               res.out = file.path(res.out, "Overlapped_M22_Hypoblast_0.005"),
               pd.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))



# ====================
# Last part: Save Data ----
# ====================
saveRDS(sr.yhw.no22, "R/CodeData/project-xuyanwen_aneuploid.rds")
save.image("R/CodeData/project-xuyanwen_aneuploid.RData")
# load("R/CodeData/project-xuyanwen_aneuploid.RData")
# load("/home/yhw/bioinfo/project-xuyanwen/aneuploid/R/CodeData/project-xuyanwen_aneuploid.RData")
# rm(list = ls())
