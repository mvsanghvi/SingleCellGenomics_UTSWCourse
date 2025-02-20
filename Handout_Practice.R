.libPaths("~/.conda/envs/sc_nanocourse_feb20_morning/lib/R/library")
dir.create("~/feb20_morning", showWarnings = FALSE, recursive = TRUE)
setwd("~/feb20_morning")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  library(celldex)
  library(SingleR)
  library(scDblFinder)
  library(SingleCellExperiment)
})
if(!file.exists("~/feb20_morning/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5")) {
  system("wget 'https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5'")
}
h5 <- Read10X_h5("~/feb20_morning/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
obj <- CreateSeuratObject(counts = h5, project = "hands-on")

#Introduction to the Seurat object
#1. Get the number of genes/features and number of cells/barcodes for this dataset.
dim(obj)
#2. Display list of genes (first 20)
rownames(obj) %>% head(20)
#3. Get the cell metadata table from the object, save it to a variable called meta, and inspect its content.
meta <- obj@meta.data
head(meta)
#4. Quantile statistics for number of UMIs
summary(meta$nCount_RNA)

#Basic quality control
#1. Run normalization
obj <- NormalizeData(obj)
#2. Feature selection, keep 2000 most variable features.
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
#3. Scale the data
obj <- ScaleData(obj)
#4. Run PCA. Create Elbow plot to determine the number of PCs to keep for next steps.
obj <- RunPCA(obj)
ElbowPlot(obj, ndims = 50)
#5. Calculate clustering with Louvain algorithm.
obj <- FindNeighbors(obj)
obj <- FindClusters(obj)
#6. Compute UMAP
obj <- RunUMAP(obj, dims=1:10)
DimPlot(obj)
#7. We can also compute tSNE projection with RunTSNE
obj <- RunTSNE(obj, dims = 1:10)
DimPlot(obj, reduction = "tsne")

#Cell type annotation
#1. Use SingleR to perform automated cell type annotation
humanprim.ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(obj)
humanprim.main <- SingleR(test = sce, assay.type.test = 1, ref = humanprim.ref, labels = humanprim.ref$label.main)
obj@meta.data$humanprim.main <- humanprim.main$pruned.labels
DimPlot(obj, group.by = "humanprim.main", label = T , repel = T)
#2. Manual cell type annotation with marker genes

gene_list <- list(
  `CD4 T` = c('IL7R', 'MAL', 'LTB', 'CD4', 'LDHB', 'TPT1', 'TRAC', 'TMSB10', 'CD3D', 'CD3G'),
  `CD8 T` = c('CD8B', 'CD8A', 'CD3D', 'TMSB10', 'HCST', 'CD3G', 'LINC02446', 'CTSW', 'CD3E', 'TRAC'),
  `NK` = c('NKG7', 'KLRD1', 'TYROBP', 'GNLY', 'FCER1G', 'PRF1', 'CD247', 'KLRF1', 'CST7', 'GZMB'),
  `other T` = c('CD3D', 'TRDC', 'GZMK', 'KLRB1', 'NKG7', 'TRGC2', 'CST7', 'LYAR', 'KLRG1', 'GZMA'),
  `B` = c('CD79A', 'RALGPS2', 'CD79B', 'MS4A1', 'BANK1', 'CD74', 'TNFRSF13C', 'HLA-DQA1', 'IGHM', 'MEF2C'),
  `DC` = c('CD74', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'CCDC88A', 'HLA-DRA', 'HLA-DMA', 'CST3', 'HLA-DQB1', 'HLA-DRB1'),
  `Monocyte` = c('CTSS', 'FCN1', 'NEAT1', 'LYZ', 'PSAP', 'S100A9', 'AIF1', 'MNDA', 'SERPINA1', 'TYROBP')
)
obj <- AddModuleScore(
  object = obj,
  features = gene_list,
  name = names(gene_list)
)
ft_ids <- tail(colnames(obj@meta.data),length(gene_list))
FeaturePlot(obj, features = ft_ids, ncol = 3)
DotPlot(obj, features = ft_ids, group.by = "seurat_clusters")
VlnPlot(obj, features = ft_ids, group.by = "seurat_clusters",
        stack = T, flip = T)
##Manually add label to each cluster based on expression level. 
obj$`manual_labels` <- case_when(
  obj$seurat_clusters %in% c(2,3,4,5,7,10) ~ "T",
  obj$seurat_clusters %in% c(1,12) ~ "NK",
  obj$seurat_clusters %in% c(0,8) ~ "Mono",
  obj$seurat_clusters %in% c(11) ~ "DC",
  obj$seurat_clusters %in% c(6,9) ~ "B"
)
DimPlot(obj, group.by = "manual_labels", label = TRUE, repel = TRUE)
#3. Identify marker genes for each cluster
all.markers <- FindAllMarkers(obj, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
DoHeatmap(obj, features = top5_markers$gene, slot = "scale.data", size = 3) +
  theme(axis.text = element_text(size = 8))