library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scDblFinder)
library(harmony)
library(monocle3)

loc <- "/project/nanocourse/SCGenomics/shared/Feb20_afternoon/5-timepoints/"

folders <- c("GSM4994960_E18-ME","GSM4994962_Pre-D5-BL6","GSM4994963_Pre-BL6",
             "GSM2759554_5wk-1","GSM4994967_Adult-BL6")
samples <- sapply(strsplit(folders,"_"),'[',1)
tps <- sapply(strsplit(folders,"_"),'[',2)

commonGenes <- read.table(paste0(loc,"commonGenes.tsv"), sep = "\t", header = F)$V1

data.list <- list()
obj.list <- list()
epi.list <- list()
epi.cells <- list(c(0,1,2,6,7),c(2,5,10),c(1,2,5),c(0,1,4,7),c(0,2,7))
for (i in 1:length(folders)){
  data.list[[i]] <- Read10X(data.dir = paste0(loc,folders[i]))
  obj.list[[i]] <- CreateSeuratObject(counts = data.list[[i]], project = tps[i], min.cells = 3, min.features = 200)
  obj.list[[i]]$sample <- samples[i]
  obj.list[[i]]$tps <- tps[i]
  obj.list[[i]][["percent.mt"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^mt-")
  obj.list[[i]][["percent.ribo"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^Rp[sl]")
  
  obj.list[[i]] <- subset(obj.list[[i]], features = commonGenes) #common genes between 2 different GEO accession
  obj.list[[i]] <- subset(obj.list[[i]], subset = percent.mt < 10 & nFeature_RNA > 500)
  obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE) #LogNormalize
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj.list[[i]] <- ScaleData(obj.list[[i]], features = VariableFeatures(object = obj.list[[i]]))
  obj.list[[i]] <- RunPCA(obj.list[[i]], features = VariableFeatures(object = obj.list[[i]]))
  num_PC <- 30
  obj.list[[i]] <- obj.list[[i]] %>%
    RunUMAP(reduction = "pca", dims = 1:num_PC) %>%
    FindNeighbors(reduction = "pca", dims = 1:num_PC) %>%
    FindClusters(resolution = 0.2) %>%
    identity()
  
  #for doublet detection
  sce <- as.SingleCellExperiment(DietSeurat(obj.list[[i]], graphs = c("pca","umap")))
  set.seed(42)
  sce <- scDblFinder(sce)
  obj.list[[i]]$db_score <- sce$scDblFinder.score
  obj.list[[i]]$db_type <- factor(sce$scDblFinder.class, levels = c("singlet","doublet"))
  
  epi.list[[i]] <- subset(obj.list[[i]], subset = seurat_clusters %in% epi.cells[[i]] )
}

# max_nFeatures <- c(5000,6000,6000,3000,4000)
# i <- 1
# FeaturePlot(obj.list[[i]], features = c("nFeature_RNA")) +
#   scale_color_gradient2(low = "blue", mid = "whitesmoke", high = "red", midpoint = max_nFeatures[i]) +
#   ggtitle(tps[i])
# 
# plot.list <- list()
# for (i in 1:length(folders)){
#   # plot.list[[i]] <- DimPlot(obj.list[[i]], reduction = "umap",label = TRUE, pt.size = .1) + 
#   #   NoLegend() + ggtitle(tps[i])
#   # plot.list[[i]] <- FeaturePlot(obj.list[[i]], features = c("Epcam")) +
#   #    ggtitle("") + theme(legend.position = "top", legend.text = element_text(angle = 90, vjust = 0.5, hjust=1))
#   # plot.list[[i]] <- DimPlot(obj.list[[i]], reduction = "umap", group.by = "db_type") + 
#   #   NoLegend() + ggtitle(tps[i])
#   # plot.list[[i]] <- DotPlot(obj.list[[i]], features = c("Epcam")) +
#   #   scale_color_gradient2(low = "blue", mid = "whitesmoke", high = "red", midpoint = 0)
# }
# plot_grid(plotlist = plot.list, ncol = 5)


epi.all <- merge(x = epi.list[[1]], y = epi.list[2:length(folders)], add.cell.ids = tps)
epi.all <- NormalizeData(epi.all, verbose = FALSE) #LogNormalize
epi.all <- FindVariableFeatures(epi.all, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
epi.all <- ScaleData(epi.all, features = VariableFeatures(object = epi.all))
epi.all <- RunPCA(epi.all, features = VariableFeatures(object = epi.all))
epi.all$sorted_tps <- factor(epi.all$tps, levels = tps)

ElbowPlot(epi.all, ndims=50)

#before harmony
p1 <- DimPlot(object = epi.all, reduction = "pca", pt.size = .1, group.by = "sorted_tps")
p2 <- VlnPlot(object = epi.all, features = "PC_1", group.by = "sorted_tps",  pt.size = .1)
plot_grid(p1,p2)

epi.all <- RunHarmony(epi.all, "sorted_tps", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
#after harmony
p1 <- DimPlot(object = epi.all, reduction = "harmony", pt.size = .1, group.by = "sorted_tps")
p2 <- VlnPlot(object = epi.all, features = "harmony_1", group.by = "sorted_tps",  pt.size = .1)
plot_grid(p1,p2)

# t-SNE and Clustering
num_PC <- 30
epi.all <- epi.all %>%
  RunUMAP(reduction = "harmony", dims = 1:num_PC) %>%
  FindNeighbors(reduction = "harmony", dims = 1:num_PC) %>%
  FindClusters(resolution = 0.2) %>%
  identity()

DimPlot(epi.all, reduction = "umap",label = TRUE, pt.size = .1)
markers <- c("Krt14","Acta2","Csn3","Elf5","Prlr","Areg","Hmgb2","Mki67","Igfbp7","Fabp4")
FeaturePlot(epi.all, features = markers)
DotPlot(epi.all, features = markers) +
  scale_color_gradient2(low = "blue", mid = "whitesmoke", high = "red", midpoint = 0)

# epi.all$celltype <- epi.all$tps
# epi.all$celltype[which(epi.all$py_louvainID %in% c(6,9))] <- "basal"
# epi.all$celltype[which(epi.all$py_louvainID %in% c(0,2,11,13))] <- "LP"
# epi.all$celltype[which(epi.all$py_louvainID %in% c(1,4,7,8,10))] <- "ML"
# epi.all$celltype[which(epi.all$py_louvainID %in% c(3))] <- "cycling"
# epi.all$celltype[which(epi.all$py_louvainID %in% c(12))] <- "stromal"


seu <- readRDS(paste0(loc,"all_epithelial.rds"))
seu <- JoinLayers(seu)
##set up for FDL
library(reticulate)
use_condaenv("nanocourse-trajectory")
sc <- import("scanpy", convert = FALSE)
sce <- import("scanpy.external", convert = FALSE)
ad <- import("anndata", convert = FALSE)
adata <- sc$AnnData(
  X = np_array(t(seu@assays[["RNA"]]@layers[["scale.data"]]), dtype="float32"),
  var = data.frame(geneName = VariableFeatures(seu)))
adata$obs['tps'] = seu@meta.data[, c("tps")]
adata$var_names <- VariableFeatures(seu)
adata$obs_names <- colnames(seu)
adata$obsm$update(X_pca = Embeddings(seu, "pca"))
sce$pp$harmony_integrate(adata, key = "tps", adjusted_basis = "X_pca_harmony")
sc$pp$neighbors(adata, use_rep = "X_pca_harmony")
sc$tl$umap(adata)
sc$tl$louvain(adata)

# FDL
sc$tl$draw_graph(adata, layout = "fa", init_pos = "X_umap")
oupDR <- py_to_r(adata$obsm['X_draw_graph_fa'])
rownames(oupDR) <- colnames(seu)
colnames(oupDR) <- c("FDL_1","FDL_2")
oupDR = oupDR / 10^(floor(log10(diff(range(oupDR))))-1)
seu[["fdl"]] <- CreateDimReducObject(embeddings = oupDR, key = "FDL_",
                                     assay = DefaultAssay(seu))
oupObs <- py_to_r(adata$obs)
seu$py_louvainID <- oupObs$louvain

# exprs <- data.frame(FetchData(epi.all, vars = c("py_louvainID","sorted_tps")))
# df <- exprs %>%
#   dplyr::group_by(sorted_tps, py_louvainID) %>%
#   dplyr::count(name = "Frequency")
# total_df <- df %>%
#   group_by(sorted_tps) %>%
#   summarise(Total = sum(Frequency))
# ggplot() +
#   geom_bar(data = df, aes(x = sorted_tps, y = Frequency, fill =  py_louvainID), width = 0.5, stat = "identity", position = "fill") +
#   geom_text(data = total_df, aes(y = 100, x = sorted_tps, label = Total), size = 4, position = position_fill(vjust = 1.02)) 

# Monocle on UMAP
cds <- SeuratWrappers::as.cell_data_set(epi.all)
set.seed(42)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE) # for single trajectory
plot_cells(cds, color_cells_by="seurat_clusters", group_label_size=4, graph_label_size=3,
           label_cell_groups=FALSE, label_principal_points=TRUE, label_groups_by_cluster=FALSE) #+
#  scale_color_manual(values = my_colors_15)
cds <- order_cells(cds, root_pr_nodes="Y_34")
plot_cells(cds, color_cells_by="pseudotime", label_groups_by_cluster=FALSE, 
           label_leaves=FALSE, label_branch_points=FALSE)

# Monocle on FDL
cds2 <- SeuratWrappers::as.cell_data_set(epi.all)
cds2@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- 
  cds2@int_colData@listData[["reducedDims"]]@listData[["FDL"]]
set.seed(42)
cds2 <- cluster_cells(cds2, reduction_method = "UMAP")
cds2 <- learn_graph(cds2, use_partition=FALSE, close_loop=FALSE)
plot_cells(cds2, color_cells_by="py_louvainID", group_label_size=4, graph_label_size=3,
           label_cell_groups=FALSE, label_principal_points=TRUE, label_groups_by_cluster=FALSE)
cds2 <- order_cells(cds2, root_pr_nodes="Y_1")
plot_cells(cds2, color_cells_by="pseudotime", label_groups_by_cluster=FALSE, 
           label_leaves=FALSE, label_branch_points=FALSE)


#saveRDS(epi.all, file = paste0(loc,"all_epithelial.rds"))
#saveRDS(obj.list, file = paste0(loc,"ind_samples.rds"))

###############################################################################
#from Cheng et al.
epi.all <- readRDS(paste0(loc,"all_epithelial.rds"))

seurat_epi <- SplitObject(epi.all, split.by = "sorted_tps")
seurat_epi <- lapply(seurat_epi, NormalizeData)
seurat_epi <- lapply(seurat_epi, FindVariableFeatures, nfeatures = 2000)
options(future.globals.maxSize = 8000 * 1024^2)
anchor_features <- SelectIntegrationFeatures(seurat_epi, nfeatures = 2000, verbose = FALSE)
anchors <- FindIntegrationAnchors(seurat_epi, verbose = FALSE,
                                  anchor.features = anchor_features)
seurat_int <- IntegrateData(anchors, verbose = FALSE)
DefaultAssay(seurat_int) <- "integrated"
seurat_int <- ScaleData(seurat_int, verbose = FALSE)
seurat_int <- RunPCA(seurat_int, npcs = 30, verbose = FALSE)
seurat_int <- RunUMAP(seurat_int, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat_int <- FindNeighbors(seurat_int, dims = 1:30, verbose = FALSE)
seurat_int <- FindClusters(seurat_int, resolution = 0.2, verbose = FALSE)

DimPlot(seurat_int, reduction = "umap",label = TRUE, pt.size = .1)
##############################################################################

#using IntegrateLayers()
library(SeuratWrappers)

epi.all <- readRDS(paste0(loc,"all_epithelial.rds"))
epi.all[["RNA"]] <- split(epi.all[["RNA"]], f = epi.all$sorted_tps)
# run standard analysis workflow
epi.all <- NormalizeData(epi.all)
epi.all <- FindVariableFeatures(epi.all)
epi.all <- ScaleData(epi.all)
epi.all <- RunPCA(epi.all)
#integration
epi.all <- IntegrateLayers(object = epi.all, method = HarmonyIntegration, orig.reduction = "pca", 
                           new.reduction = "harmony", verbose = FALSE)
# re-join layers after integration
epi.all[["RNA"]] <- JoinLayers(epi.all[["RNA"]])
epi.all <- FindNeighbors(epi.all, reduction = "harmony", dims = 1:30)
epi.all <- FindClusters(epi.all, resolution = 0.2, cluster.name = "harmony_clusters")
epi.all <- RunUMAP(epi.all, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")