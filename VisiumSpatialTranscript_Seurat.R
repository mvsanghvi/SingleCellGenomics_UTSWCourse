# Seurat workflow for Visium spatial transcriptomics
# https://satijalab.org/seurat/articles/spatial_vignette.html#x-visium

# Load libraries ----------------------------------------------------------

# if needed, install this recently updated github version of the SeuratData package
#remotes::install_github("satijalab/seurat-data")

# install this newer version of Seurat for Visium analysis
#install.packages("Seurat", version="5.1.0")

library(Matrix)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(shiny)
library(future)

options(future.globals.maxSize = 8000 * 1024^2)

# Load spatial data -------------------------------------------------------

# Using Seurat shortcut
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

# Load data without shortcut if you have your own data or if you download from 10X website
# datadir <- '/project/InternalMedicine/Chan_lab/shared/for_Shawn/SCG_Nanocourse/Spatial_Practice_Mouse_Brain'
# brain <- Load10X_Spatial(datadir, filename = "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5")

# Quality Control and Pre-Processing --------------------------------------

# Visualize molecular count distribution
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Calculate QC metrics
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-") # For mitochondrial genes

# Visualize basic QC metrics
VlnPlot(brain, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# If desired, filter out low gene expression and transcript count and high mitochondrial gene content.
# For this example, we will keep all spots except those with 0 features (genes) and counts (transcripts).
brain <- subset(brain, subset = nFeature_Spatial > 0 & nCount_Spatial > 0)

# Normalization with SCTransform
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# Unsupervised Analysis ---------------------------------------------------

# Visualize expression of genes of interest (Hpca, Ttr) on spatial map
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

# Adjust transparency and size of spots
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

# Dimensional Reduction ------------------------------------------------
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# Cluster visualization
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)

# Interactive Plots -
# SpatialDimPlot(brain, interactive = TRUE)
# SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
# LinkedDimPlot(brain)

# Identify spatially variable features

# First method: FindMarkers()
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

# Second method: FindSpatiallyVariableFeatures()
SpatiallyVariableFeatures_workaround <- function(object, assay="SCT", selection.method = "moransi") {
  #' This is a workaround function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[
    data[[paste0(selection.method, ".spatially.variable")]] & (!is.na(data[[paste0(selection.method, ".spatially.variable")]])), 
    moransi_cols
  ]
  
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  
  # Return row names of the sorted data frame
  rownames(sorted_data)
}

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures_workaround(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

# Subset the Data -----------------------------------------------------

# subset cortex based on desired clusters
cortex <- subset(brain, idents = c(2, 3, 4, 6, 7))

# pull spatial coordinates of subset and store in meta.data
coords <- GetTissueCoordinates(cortex)
rownames(coords) <- colnames(cortex)
cortex$imagerow <- coords$x
cortex$imagecol <- coords$y

# check to make sure that coordinates have been added
head(cortex@meta.data) 

# subset based on specific spatial condition representing anatomical ROI
# there may be some trial-and-error here with the cropping window
condition <- cortex$imagerow < 8000 & cortex$imagecol < 9750
cortex <- subset(cortex, cells = WhichCells(cortex, cells=which(condition)))

# use these plots to check your subsetting
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

# Cell Type Annotation ----------------------------------------------------
# load scRNAseq reference dataset
allen_reference <- readRDS("/project/nanocourse/SCGenomics/shared/Feb21_Spatial/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells. This speeds up SCTransform dramatically with no loss in performance.
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# after subsetting, re-normalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

# perform cell type annotation
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

# visualize prediction scores on spatial map
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

# add cell type labels to the meta.data of the original Seurat object
labels = GetAssayData(cortex,assay = "predictions")
labels = t(labels)
labels = cbind(labels, apply(labels, 1, function(x) colnames(labels)[which.max(x)]))
colnames(labels)[25] = "FTA_Spot_Annotation"

# visualize cell type labels on the spatial map
cortex <- AddMetaData(cortex, metadata = labels)
cortex@meta.data <- cortex@meta.data[, -c(13:36)]
Idents(cortex) <- cortex$FTA_Spot_Annotation
SpatialDimPlot(cortex, label = TRUE, label.size = 3, stroke = NA)

# Supervised Analysis -----------------------------------------------------
# use FindSpatiallyVariableFeatures() to predict cell types on spatial map
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), slot = 'data')
top.clusters <- head(SpatiallyVariableFeatures_workaround(cortex, assay = 'predictions', selection.method = "moransi"), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
SpatialFeaturePlot(cortex, features = c("Astro", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

# Integrate Multiple Datasets ----------------------------------------------
# Load second dataset
brain2 <- LoadData("stxBrain", type = "posterior1")

# Calculate QC metrics
brain2[["percent.mt"]] <- PercentageFeatureSet(brain2, pattern = "^mt-") # For mitochondrial genes

# Visualize basic QC metrics
VlnPlot(brain2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)

# If desired, filter out low gene expression and transcript count and high mitochondrial gene content.
# For this example, we will keep all spots except those with 0 features (genes) and counts (transcripts).
brain2 <- subset(brain2, subset = nFeature_Spatial > 0 & nCount_Spatial > 0)

brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

# merge the two datasets
brain.merge <- merge(brain, brain2)

# set default assay to the SCTransform normalized one
DefaultAssay(brain.merge) <- "SCT"

# select variable features to focus on features with highest variation for more informative downstream analysis
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))

# Unintegrated dimensional reduction and clustering
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

# Visualize unintegrated dataset
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))

# perform integration for batch correction. We use Harmony Integration, but feel free to experiment with other methods.
library(harmony)

brain.merge <- IntegrateLayers(
  object = brain.merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# Integrated clustering
brain.merge <- FindNeighbors(brain.merge, reduction = "harmony", dims = 1:30)
brain.merge <- FindClusters(brain.merge, cluster.name = "harmony_clusters")
brain.merge <- RunUMAP(brain.merge, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# Visualize integrated dataset
DimPlot(brain.merge, reduction = "umap.harmony", group.by = c("ident", "orig.ident"))
SpatialDimPlot(brain.merge)
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))