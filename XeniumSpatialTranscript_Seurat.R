# Seurat workflow for Xenium spatial transcriptomics
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ

# Load libraries ----------------------------------------------------------
library(Seurat)
library(future)
library(ggplot2)
#library(sf)
library(spacexr)

plan("multisession", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)

# Use the "Packages" interface (right) to search "Seurat" and check which
# version got loaded. Switch to Seurat 5.0.3 and unselect SeuratObject
# if previously running a newer version of Seurat (uncheck other version
# check this version)

# Load data ---------------------------------------------------------------
# Set path to Xenium data
path <- "/project/nanocourse/SCGenomics/shared/Feb21_Spatial/Xenium_Practice_Mouse_Brain"
# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov")

# Quality Control and Pre-Processing --------------------------------------

# Visualize basic QC metrics
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 3, pt.size = 0)

# Remove cells with 0 counts and/or 0 features
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0 & nFeature_Xenium > 0)

# Normalize data
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")

# Unsupervised Analysis ---------------------------------------------------

# Visualize gene expression on spatial map
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)
ImageFeaturePlot(xenium.obj, features = c("Cux2", "Rorb", "Bcl11b", "Foxp2"), max.cutoff = c(25,
                                                                                             35, 12, 10), size = 0.75, cols = c("white", "red"))

# Crop Xenium data
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords

# Visualize cropped area with cell segmentations & selected genes
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Gad1", "Sst", "Npy2r", "Pvalb", "Nrn1"), nmols = 10000)

# Dimensional Reduction and Clustering ------------------------------------
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)

# Visualize clusters on UMAP
DimPlot(xenium.obj)

# Visualize gene expression on UMAP
FeaturePlot(xenium.obj, features = c("Cux2", "Bcl11b", "Foxp2", "Gad1", "Sst", "Gfap"))

# Visualize clusters on spatial map
ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)

# Cropping to Region of Interest ------------------------------------------

# Visualize expression of region-specific gene
p1 <- ImageFeaturePlot(xenium.obj, features = "Slc17a7", axes = TRUE, max.cutoff = "q90")
p1

# Crop Xenium sample based on anatomical region of interest
crop <- Crop(xenium.obj[["fov"]], x = c(600, 2100), y = c(900, 4700))
xenium.obj[["crop"]] <- crop
p2 <- ImageFeaturePlot(xenium.obj, fov = "crop", features = "Slc17a7", size = 1, axes = TRUE, max.cutoff = "q90")
p2

# Cell Type Annotation ----------------------------------------------------

# Create query object based on Xenium Seurat object
query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")[, Cells(xenium.obj[["crop"]])]
coords <- GetTissueCoordinates(xenium.obj[["crop"]], which = "centroids")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# allen.cortex.ref can be downloaded here:
# https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1

# Load scRNAseq reference dataset
allen.cortex.ref <- readRDS("/project/nanocourse/SCGenomics/shared/Feb21_Spatial/allen_cortex.rds")
allen.cortex.ref <- UpdateSeuratObject(allen.cortex.ref)

# Create reference object based on scRNAseq reference dataset
Idents(allen.cortex.ref) <- "subclass"
# remove CR cells because there aren't enough of them for annotation
allen.cortex.ref <- subset(allen.cortex.ref, subset = subclass != "CR")
counts <- GetAssayData(allen.cortex.ref, assay = "RNA", slot = "counts")
cluster <- as.factor(allen.cortex.ref$subclass)
names(cluster) <- colnames(allen.cortex.ref)
nUMI <- allen.cortex.ref$nCount_RNA
names(nUMI) <- colnames(allen.cortex.ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)

# run RCTD with doublet mode for Xenium
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# store RCTD results in Xenium Seurat object
annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

# Supervised Analysis - Niche/Neighborhood Analysis -----------------------

# Run niche analysis
xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "crop", group.by = "predicted.celltype",
                              niches.k = 5, neighbors.k = 30)

# Plot cell types and niches on spatial map
celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
                              dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot

# View composition of each niche based on cell type
table(xenium.obj$predicted.celltype, xenium.obj$niches)