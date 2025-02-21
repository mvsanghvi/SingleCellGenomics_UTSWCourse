setwd("~/BICF Nanocourse/Single Cell Genomics")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
set.seed(1234)

#Load Signac input files

# Unnormalized data (raw counts)
counts <- Read10X_h5(filename = "Signac Data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")

# Additional cell-level metadata
metadata <- read.csv(
  file = "Signac Data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
  header = TRUE,
  row.names = 1
)

#Create ChromatinAssay and Seurat ojbects

# Create ChromatinAssay object
## The expected format of the input matrix is features x cells. 
## A set of genomic ranges must be supplied along with the matrix
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "Signac Data/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

# Creat a Seurat object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc

## An object of class Seurat 
## 165434 features across 10246 samples within 1 assay 
## Active assay: peaks (165434 features, 0 variable features)
##  2 layers present: counts, data

pbmc[['peaks']]

## ChromatinAssay data with 165434 features for 10246 cells
## Variable features: 0 
## Genome: 
## Annotation present: FALSE 
## Motifs present: FALSE 
## Fragment files: 1

granges(pbmc)

## GRanges object with 165434 ranges and 0 metadata columns:
##              seqnames        ranges strand
##                 <Rle>     <IRanges>  <Rle>
##        [1]       chr1    9772-10660      *
##        [2]       chr1 180712-181178      *
##        [3]       chr1 181200-181607      *
##        [4]       chr1 191183-192084      *
##        [5]       chr1 267576-268461      *
##        ...        ...           ...    ...
##   [165430] KI270713.1   13054-13909      *
##   [165431] KI270713.1   15212-15933      *
##   [165432] KI270713.1   21459-22358      *
##   [165433] KI270713.1   29676-30535      *
##   [165434] KI270713.1   36913-37813      *
##   -------
##   seqinfo: 35 sequences from an unspecified genome; no seqlengths



#Remove features

# Remove features that correspond to chromosome scaffolds or other non-standard chromosomes
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]

pbmc

## An object of class Seurat 
## 165376 features across 10246 samples within 1 assay 
## Active assay: peaks (165376 features, 0 variable features)
##  2 layers present: counts, data



#Add gene annotations to the pbmc objec

library(AnnotationHub)
ah <- AnnotationHub()

# The reference package 10x Genomics used to perform the mapping was “GRCh38-2020-A”,
# which corresponds to the Ensembl v98 patch release.
# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")

## AnnotationHub with 1 record
## # snapshotDate(): 2024-04-30
## # names(): AH75011
## # $dataprovider: Ensembl
## # $species: Homo sapiens
## # $rdataclass: EnsDb
## # $rdatadateadded: 2019-05-02
## # $title: Ensembl 98 EnsDb for Homo sapiens
## # $description: Gene and protein annotations for Homo sapiens based on Ensembl version 98.
## # $taxonomyid: 9606
## # $genome: GRCh38
## # $sourcetype: ensembl
## # $sourceurl: http://www.ensembl.org
## # $sourcesize: NA
## # $tags: c("98", "AHEnsDbs", "Annotation", "EnsDb", "Ensembl", "Gene", "Protein", "Transcript") 
## # retrieve record with 'object[["AH75011"]]'

ensdb_v98 <- ah[["AH75011"]] # retrieve record'

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
head(annotations)

## GRanges object with 6 ranges and 5 metadata columns:
##                   seqnames        ranges strand |           tx_id   gene_name         gene_id   gene_biotype
##                      <Rle>     <IRanges>  <Rle> |     <character> <character>     <character>    <character>
##   ENSE00001489430        X 276322-276394      + | ENST00000399012      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001536003        X 276324-276394      + | ENST00000484611      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00002160563        X 276353-276394      + | ENST00000430923      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001750899        X 281055-281121      + | ENST00000445062      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001719251        X 281194-281256      + | ENST00000429181      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001489388        X 281381-281684      + | ENST00000381657      PLCXD1 ENSG00000182378 protein_coding
##                       type
##                   <factor>
##   ENSE00001489430     exon
##   ENSE00001536003     exon
##   ENSE00002160563     exon
##   ENSE00001750899     exon
##   ENSE00001719251     exon
##   ENSE00001489388     exon
##   -------
##   seqinfo: 25 sequences (1 circular) from GRCh38 genome

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

head(annotations)

## GRanges object with 6 ranges and 5 metadata columns:
##                   seqnames        ranges strand |           tx_id   gene_name         gene_id   gene_biotype
##                      <Rle>     <IRanges>  <Rle> |     <character> <character>     <character>    <character>
##   ENSE00001489430     chrX 276322-276394      + | ENST00000399012      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001536003     chrX 276324-276394      + | ENST00000484611      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00002160563     chrX 276353-276394      + | ENST00000430923      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001750899     chrX 281055-281121      + | ENST00000445062      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001719251     chrX 281194-281256      + | ENST00000429181      PLCXD1 ENSG00000182378 protein_coding
##   ENSE00001489388     chrX 281381-281684      + | ENST00000381657      PLCXD1 ENSG00000182378 protein_coding
##                       type
##                   <factor>
##   ENSE00001489430     exon
##   ENSE00001536003     exon
##   ENSE00002160563     exon
##   ENSE00001750899     exon
##   ENSE00001719251     exon
##   ENSE00001489388     exon
##   -------
##   seqinfo: 25 sequences (1 circular) from hg38 genome

# add the gene information to the object
Annotation(pbmc) <- annotations



#Compute QC Metrics

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)
head( pbmc$nucleosome_signal)

## AAACGAAAGAGAGGTA-1 AAACGAAAGCAGGAGG-1 AAACGAAAGGAAGAAC-1 AAACGAAAGTCGACCC-1 AAACGAACAAGCACTT-1 AAACGAACAAGCGGTA-1 
##          0.5054022          0.4906667          0.6323430          0.6112150          0.4175417          0.8787062

## compute TSS enrichment score per cell (it may take ~20minutes)
#pbmc <- TSSEnrichment(object = pbmc)
pbmc$TSS.enrichment <- readRDS("TSS.enrichment.rds")

head(pbmc$TSS.enrichment)

## AAACGAAAGAGAGGTA-1 AAACGAAAGCAGGAGG-1 AAACGAAAGGAAGAAC-1 AAACGAAAGTCGACCC-1 AAACGAACAAGCACTT-1 AAACGAACAAGCGGTA-1 
##           6.403815           8.516617           4.847291           5.718537           6.887816           5.187235

# Add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

head(pbmc$pct_reads_in_peaks)

## AAACGAAAGAGAGGTA-1 AAACGAAAGCAGGAGG-1 AAACGAAAGGAAGAAC-1 AAACGAAAGTCGACCC-1 AAACGAACAAGCACTT-1 AAACGAACAAGCGGTA-1 
##           75.77924           76.34671           67.93552           72.71074           74.96191           61.05524

# add blacklist ratio
pbmc$blacklist_ratio <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

head(pbmc$blacklist_ratio)

## AAACGAAAGAGAGGTA-1 AAACGAAAGCAGGAGG-1 AAACGAAAGGAAGAAC-1 AAACGAAAGTCGACCC-1 AAACGAACAAGCACTT-1 AAACGAACAAGCGGTA-1 
##       0.0001839701       0.0001505231       0.0006087606       0.0002740340       0.0002106150       0.0002040816



# Density scatter plot visualize the relationship between variables stored in the object metadata
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Outlier cells: nucleosome signal score (mononucleosomal / nucleosome-free ratio, NS) > 4
# Typical cells: NS < 4
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# Outlier cells have different nucleosomal banding patterns
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 3
)



# Remove cells that are outliers for these QC metrics
pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

pbmc

## An object of class Seurat 
## 165376 features across 9651 samples within 1 assay 
## Active assay: peaks (165376 features, 0 variable features)
##  2 layers present: counts, data



#Normalization and linear dimensional reduction

# TF-IDF normalization
pbmc <- RunTFIDF(pbmc)

# Select n% top features only because of the low dynamic range of scATAC-seq data
# q0, all features (peaks); q75, top 25% all peaks
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')

# Run singular value decomposition (SVD) on the TD-IDF matrix, using the features selected
# This returns a reduced dimension representation of the object (similar to PCA of scRNA-seq data)
# This combined steps is known as latent semantic indexing (LSI) 
pbmc <- RunSVD(pbmc)



#Non-linear dimension reduction and clustering

# The first LSI component often captures sequencing depth rather than biological variation. 
# If this is the case, the component should be removed from downstream analysis.
# DepthCor assess the correlation between each LSI component and sequencing depth
DepthCor(pbmc)

#  A very strong correlation between the first LSI component and the sequencing depth
#  We need to perform downstream steps without this component

# Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)

# Computes the k nearest neighbors
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30) # reduction='lsi'

# Identify clusters of cells by a shared NN modularity optimization
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3) # Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm)

DimPlot(object = pbmc, label = TRUE) + NoLegend()



#Gene activity matrix

# Quantify the gene activity by assessing the chromatin accessibility associated with the gene
# Here we will use a simple approach of summing the fragments intersecting the gene body and promoter region
#gene.activities <- GeneActivity(pbmc) #it may take 20-30 minutes
gene.activities <- readRDS("gene.activities.rds")

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)



# Visualize the activities of canonical marker genes 
DefaultAssay(pbmc) <- 'RNA'
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# We can discern populations of monocytes, B, T, and NK cells based on these gene activity profiles. 
# However, further subdivision of these cell types is challenging





#Integrating with scRNA-seq data

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("Signac Data/pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

# We aim to identify shared correlation patterns in the gene activity matrix and scRNA-seq dataset 
# to identify matched biological states across the two modalities

## Find a set of anchors between a reference and query object
#transfer.anchors <- FindTransferAnchors(
#  reference = pbmc_rna,
#  query = pbmc,
#  reduction = 'cca'
#)
transfer.anchors <- readRDS("transfer.anchors.rds")

# Transfer categorical or continuous data across single-cell datasets
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)
head(predicted.labels)

##                         predicted.id prediction.score.CD4.Memory prediction.score.CD14..Monocytes
## AAACGAAAGAGAGGTA-1   CD14+ Monocytes                  0.00000000                       1.00000000
## AAACGAAAGCAGGAGG-1         CD4 Naive                  0.06129442                       0.03582392
## AAACGAAAGGAAGAAC-1 B cell progenitor                  0.00000000                       0.00000000
## AAACGAAAGTCGACCC-1        CD4 Memory                  0.95809377                       0.00000000
## AAACGAACAAGCACTT-1        CD4 Memory                  0.74456137                       0.00000000
## AAACGAACAAGCGGTA-1            NK dim                  0.00000000                       0.00000000
##                    prediction.score.NK.dim prediction.score.pre.B.cell prediction.score.NK.bright
## AAACGAAAGAGAGGTA-1              0.00000000                 0.000000000                 0.00000000
## AAACGAAAGCAGGAGG-1              0.00000000                 0.012526765                 0.00000000
## AAACGAAAGGAAGAAC-1              0.00000000                 0.412073293                 0.00000000
## AAACGAAAGTCGACCC-1              0.00000000                 0.000000000                 0.00000000
## AAACGAACAAGCACTT-1              0.04237848                 0.006261266                 0.00000000
## AAACGAACAAGCGGTA-1              0.94911555                 0.000000000                 0.05088445
##                    prediction.score.CD4.Naive prediction.score.CD8.Naive prediction.score.pDC
## AAACGAAAGAGAGGTA-1                 0.00000000                          0                    0
## AAACGAAAGCAGGAGG-1                 0.86950867                          0                    0
## AAACGAAAGGAAGAAC-1                 0.00000000                          0                    0
## AAACGAAAGTCGACCC-1                 0.02746352                          0                    0
## AAACGAACAAGCACTT-1                 0.03790109                          0                    0
## AAACGAACAAGCGGTA-1                 0.00000000                          0                    0
##                    prediction.score.Double.negative.T.cell prediction.score.CD16..Monocytes
## AAACGAAAGAGAGGTA-1                               0.0000000                                0
## AAACGAAAGCAGGAGG-1                               0.0000000                                0
## AAACGAAAGGAAGAAC-1                               0.0000000                                0
## AAACGAAAGTCGACCC-1                               0.0144427                                0
## AAACGAACAAGCACTT-1                               0.1333011                                0
## AAACGAACAAGCGGTA-1                               0.0000000                                0
##                    prediction.score.Platelet prediction.score.CD8.effector prediction.score.B.cell.progenitor
## AAACGAAAGAGAGGTA-1                0.00000000                    0.00000000                          0.0000000
## AAACGAAAGCAGGAGG-1                0.02084622                    0.00000000                          0.0000000
## AAACGAAAGGAAGAAC-1                0.00000000                    0.00000000                          0.5879267
## AAACGAAAGTCGACCC-1                0.00000000                    0.00000000                          0.0000000
## AAACGAACAAGCACTT-1                0.00000000                    0.03559667                          0.0000000
## AAACGAACAAGCGGTA-1                0.00000000                    0.00000000                          0.0000000
##                    prediction.score.Dendritic.cell prediction.score.max
## AAACGAAAGAGAGGTA-1                               0            1.0000000
## AAACGAAAGCAGGAGG-1                               0            0.8695087
## AAACGAAAGGAAGAAC-1                               0            0.5879267
## AAACGAAAGTCGACCC-1                               0            0.9580938
## AAACGAACAAGCACTT-1                               0            0.7445614
## AAACGAACAAGCGGTA-1                               0            0.9491156

# Add the transferred cell types to the pbmc object
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

# for cells with a small number, it is difficult to figure out their precise cellular identity
predicted_id_counts <- table(pbmc$predicted.id)
predicted_id_counts

## 
##      B cell progenitor        CD14+ Monocytes        CD16+ Monocytes             CD4 Memory              CD4 Naive 
##                    294                   2541                    241                   2596                   1119 
##           CD8 effector              CD8 Naive         Dendritic cell Double negative T cell              NK bright 
##                    972                    626                     67                    172                     44 
##                 NK dim                    pDC               Platelet             pre-B cell 
##                    693                     31                     10                    245

# Identify the predicted.id values that have more than 20 cells
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
pbmc <- pbmc[, pbmc$predicted.id %in% major_predicted_ids]

# Change cell identities to the per-cell predicted labels for downstream analyses 
Idents(pbmc) <- pbmc$predicted.id



#Find differentially accessible peaks between cell types

# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

### the presto package has implemented an extremely fast Wilcoxon test
#if (!requireNamespace("remotes", quietly = TRUE))
#  install.packages('remotes')
#remotes::install_github('immunogenomics/presto')
library(presto) 

# comparison between Naive CD4 cells and CD14 monocytes
# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)

##                           p_val avg_log2FC pct.1 pct.2 p_val_adj
## chr6-13302533-13303459        0  -5.247815 0.026 0.772         0
## chr19-54207815-54208728       0  -4.384378 0.049 0.793         0
## chr17-78198651-78199583       0  -5.505552 0.024 0.760         0
## chr12-119988511-119989430     0   4.153581 0.782 0.089         0
## chr7-142808530-142809435      0   3.663154 0.758 0.087         0
## chr17-82126458-82127377       0   4.997442 0.704 0.043         0

# We can visualize these marker peaks on a violin plot, feature plot, dot plot, heat map, etc.
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

#Find the closest gene to each of these peaks

open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

# Because peak coordinates are difficult to interpret alone, 
# we need to find the closest gene to each of these peaks
closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
head(closest_genes_cd4naive)

##                           tx_id gene_name         gene_id   gene_biotype type            closest_region
## ENSE00002206071 ENST00000397558    BICDL1 ENSG00000135127 protein_coding exon chr12-119989869-119990297
## ENST00000632998 ENST00000632998     PRSS2 ENSG00000275896 protein_coding  utr  chr7-142774509-142774564
## ENST00000665763 ENST00000665763    CCDC57 ENSG00000176155 protein_coding  gap   chr17-82101867-82127691
## ENST00000645515 ENST00000645515  ATP6V0A4 ENSG00000105929 protein_coding  cds  chr7-138752625-138752837
## ENST00000545320 ENST00000545320       CD6 ENSG00000013725 protein_coding  gap   chr11-60971915-60987907
## ENST00000603548 ENST00000603548     FBXW7 ENSG00000109670 protein_coding  utr  chr4-152320544-152322880
##                              query_region distance
## ENSE00002206071 chr12-119988511-119989430      438
## ENST00000632998  chr7-142808530-142809435    33965
## ENST00000665763   chr17-82126458-82127377        0
## ENST00000645515  chr7-138752283-138753197        0
## ENST00000545320   chr11-60985909-60986801        0
## ENST00000603548  chr4-152100248-152101142   219401

closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)
head(closest_genes_cd14mono)

##                           tx_id gene_name         gene_id   gene_biotype type           closest_region
## ENST00000606214 ENST00000606214    TBC1D7 ENSG00000145979 protein_coding  gap   chr6-13267836-13305061
## ENST00000448962 ENST00000448962      RPS9 ENSG00000170889 protein_coding  gap  chr19-54201610-54231740
## ENST00000592988 ENST00000592988     AFMID ENSG00000183077 protein_coding  gap  chr17-78191061-78202498
## ENST00000336600 ENST00000336600  C6orf223 ENSG00000181577 protein_coding  utr   chr6-44003127-44007612
## ENSE00002618192 ENST00000569518    VCPKMT ENSG00000100483 protein_coding exon  chr14-50108632-50109713
## ENSE00001389095 ENST00000340607     PTGES ENSG00000148344 protein_coding exon chr9-129752887-129753042
##                             query_region distance
## ENST00000606214   chr6-13302533-13303459        0
## ENST00000448962  chr19-54207815-54208728        0
## ENST00000592988  chr17-78198651-78199583        0
## ENST00000336600   chr6-44058439-44059230    50826
## ENSE00002618192  chr14-50038381-50039286    69345
## ENSE00001389095 chr9-129776928-129777838    23885





#Plotting genomic regions

pbmc <- SortIdents(pbmc) # sort the plotting order according to similarities across the annotated cell types

# find DA peaks overlapping gene of interest
# gene of interest: CD4
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(pbmc, "CD4"))
CoveragePlot(
  object = pbmc,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)

# gene of interest: LYZ
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd14mono), LookupGeneCoords(pbmc, "LYZ"))
CoveragePlot(
  object = pbmc,
  region = "LYZ",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 5000
)