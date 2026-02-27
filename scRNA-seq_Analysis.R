library(tidyverse)
library(Seurat)
library(SeuratData)



############################################
# Data Prep
############################################


# reading scRNA-seq data
# read the counts matrix from 10X HDF5 files.
Control_mat <- Read10X_h5("~/CSTAT/C2000_Sjogren/Data/GSE175649/GSM5343233_CTRL_filtered_feature_bc_matrix.h5")
Treatment_mat <- Read10X_h5("~/CSTAT/C2000_Sjogren/Data/GSE175649/GSM5343234_NOD_filtered_feature_bc_matrix.h5")
# dgCMatrix is a compact matrix format in R.
# using float32, a dense matrix of the control data would need 27998*21386*4/1e9 = 2.4GB space on your RAM.
# Instead we are takig 250MB using dgCMatrix format.
# @i: stores the row number of the non-zero elements.
# @p: stores the number of non-zero elements in each column. p[k+1]-p[k] is the number of non-zero elements in column k.
# @x: stores the non-zero elements.
# How many genes and how many cells in each dataset measured?




# Create Seurat objects from the dgCMatrix data
Control_obj <- CreateSeuratObject(Control_mat)
Treatment_obj <- CreateSeuratObject(Treatment_mat)
# take a look at the Seurat objects
# Where are the counts stored?
# Try rownames and colnames on the objects. Are the genes measured the same? length(intersect(rownames(Control_obj), rownames(Treatment_obj)))
# Are there duplicated cell names? intersect(colnames(Treatment_obj), colnames(Control_obj))
# Duplicated cell names will be fixed by Seurat when merging objects. We will see this later.




# Later when we merge these two objects into ONE Seurat object, we would like to have a label to distinguish cells from each sample.
# Let's add labels to the Seurat object metadata to distinguish Control and Treatment
Control_obj <- AddMetaData(object = Control_obj, metadata = rep("Control", ncol(Control_obj)),
                                 col.name = "Status")
Treatment_obj <- AddMetaData(object = Treatment_obj, metadata = rep("Treatment", ncol(Treatment_obj)),
                             col.name = "Status")
# check the Seurat objects and confirm that the MetaData has been added.
# What are the other info in the MetaData? nCount_RNA and nFeature_RNA?


# quality control of the data
FeatureScatter(Control_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(Control_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)


summary(Control_obj@meta.data$nCount_RNA)
summary(Control_obj@meta.data$nFeature_RNA)
# There are cells with more than 2500 genes expressed. Likely more than one cell trapped in a droplet.
# There are cells with total RNA expression of more than 10000 and up to 200k! Certainly more than one cell in a droplet.
# There are cells with less than 100 genes expressed. Likely empty droplets.
# We would like to remove the low quality observations.

# similar plots for the Treatment sample.
FeatureScatter(Treatment_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(Treatment_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)





# Before removing low quality observations, let us add one more MetaData, percentage of Mitochondrial RNA in each observation.
# The mitochondrial RNA are named starting with "mt-" for mice. For homo-sapiens, they start with "MT-"
# For example run: rownames(Control_obj)[Control_obj@assays[["RNA"]]@layers[["counts"]]@i]
# Seurat can calculate the percentage of mitochondrial RNA in each observation and add it to the object as MetaData.
Control_obj[["percent.mt"]] <- PercentageFeatureSet(Control_obj, pattern = "^mt-")
Treatment_obj[["percent.mt"]] <- PercentageFeatureSet(Treatment_obj, pattern = "^mt-")


# Violin plots for all three important measures of observation quality:
VlnPlot(Control_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Treatment_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




# Remove low quality observations: empty droplets, high mitochondrial contamination, multiple cells per droplet.
# check how many cells are there before filtering
Control_obj <- subset(Control_obj, 
                            subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                              nCount_RNA < 10000 &
                              percent.mt < 35)
# mitochondrial content in salivary glands is estimated to be 30-35%
# How many cells are remaining after filtering?

FeatureScatter(Control_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(Control_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Repeat for the Treatment sample
# Remove low quality observations: empty droplets, high mitochondrial contamination, multiple cells per droplet.
Treatment_obj <- subset(Treatment_obj, 
                      subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                        nCount_RNA < 10000 &
                        percent.mt < 35)

FeatureScatter(Treatment_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(Treatment_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)





############################################
# Data Integration
############################################

# merge the two scRNA-seq objects
scRNA_seq_merged <- merge(Control_obj, y = Treatment_obj, project = "Sjogren")

# log normalize the counts
# similar to Bulk RNA-seq data, the SD increases with the Mean.
# try running: plot(rowMeans(Control_mat), rowSds(Control_mat), xlim = c(0,50), ylim = c(0,50))
scRNA_seq_merged <- NormalizeData(scRNA_seq_merged, normalization.method = "LogNormalize")
# two new layers added: data.1 and data.2


# How does the Mean-SD plot look like now?
plot(rowMeans(scRNA_seq_merged@assays$RNA@layers$data.1), rowSds(scRNA_seq_merged@assays$RNA@layers$data.1))
plot(rowMeans(scRNA_seq_merged@assays$RNA@layers$data.2), rowSds(scRNA_seq_merged@assays$RNA@layers$data.2))
# Note the curve is caused by the log transformation of the counts.
# Before logtransformation there was a linear relationship between the mean and sd of the genes.


scRNA_seq_merged <- FindVariableFeatures(scRNA_seq_merged, selection.method = "vst", nfeatures = 3000)
# Look at the documentation for FindVariableFeatures and read the "vst" selection method details.
# The variable feature info is stored in the Assays -> RNA -> MetaData.

# Store the Highly Variable Features from the control sample
HVG_Sample1 <- scRNA_seq_merged@assays[["RNA"]]@meta.data[["vf_vst_counts.1_variable"]]
plot(rowMeans(scRNA_seq_merged@assays$RNA@layers$data.1[HVG_Sample1,]),
     rowSds(scRNA_seq_merged@assays$RNA@layers$data.1[HVG_Sample1,]))
# It has selected a subset of the genes (3000 features) that had a large residual compared to the fittet polynomial on the mean-sd plot.
# In this way we are guaranteed to keep highly variable genes from low-expressed, medium-expressed, and highly-expressed genes
VariableFeaturePlot(scRNA_seq_merged)
# Standardized Variance is the variance remaining after fitting a polynomial to the LogNormalized mean vs sd plot.



# center and scale the data before running PCA
scRNA_seq_merged <- ScaleData(scRNA_seq_merged)

# run PCA (Principal Component Analysis, a dimensional reduction method)
scRNA_seq_merged <- RunPCA(scRNA_seq_merged)

# How many PCs should we use in downstream analysis?
ElbowPlot(scRNA_seq_merged, ndims = 50)
# 17 PCs are enough


# Integrate the data
# we could save the integrated object back into the scRNA_seq_merged or make a new object.
scRNA_seq_merged <- IntegrateLayers(object = scRNA_seq_merged, method = CCAIntegration,
                                      orig.reduction = "pca", new.reduction = "integrated.cca", verbose = TRUE)
# what can happen if we don't integrate the data? 
# Look up batch effect in scRNA-seq clustering results.

# Finding nearest neighbors: which cells are similar to each other
scRNA_seq_merged <- FindNeighbors(scRNA_seq_merged, reduction = "integrated.cca", dims = 1:20)

# Clustering cells together, based on how similar they are.
scRNA_seq_merged <- FindClusters(scRNA_seq_merged, resolution = 0.5, cluster.name = "cca_clusters")
# What will happen if we use different resolutions?
# How to decide a good resolution?

# Now we can visualize the clusters.
scRNA_seq_merged <- RunUMAP(scRNA_seq_merged, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
# Uniform Manifold Approximation and Projection (UMAP) is a fast, non-linear dimension reduction technique
# used for visualizing high-dimensional data in 2D or 3D. Older methods, tSNE


DimPlot(scRNA_seq_merged, reduction = "umap.cca",
        group.by = c("cca_clusters"),
        combine = FALSE, label.size = 2)

# try the UMAP plot with splitting by control vs treatment (split.by = "Status")
DimPlot(scRNA_seq_merged, reduction = "umap.cca",
        group.by = c("cca_clusters"),
        combine = FALSE, label.size = 2,
        split.by = "Status")
# Has clustering being successful? Are the clusters overlapping correctly?





#################################################
# Annotations
################################################

# There are many annotated data sets from SeuratData package to use.
AvailableData()
# No reference dataset for salivary glands
# For more common data, e.g. PBMC, you can use a reference annotated data set to help you annotate.

# annotation using markers
# joinig the raw count layers
scRNA_seq_merged <- JoinLayers(scRNA_seq_merged)

# What does the word "Marker" here mean? It is not exactly the same as what pathologists call marker.
Markers <- FindAllMarkers(scRNA_seq_merged, logfc.threshold = 1, min.pct = 0.25, min.diff.pct = 0.10, only.pos = TRUE)
# what are the different options?
# why "only.pos" option? When should we put it to FALSE?

Markers <- Markers %>%
  filter(p_val_adj < 0.01)
writexl::write_xlsx(Markers, path = "~/CSTAT/C2000_Sjogren/Output/scRNA_seq_markers.xlsx", col_names = TRUE)

# annotation using marker genes
# acinar cells
VlnPlot(scRNA_seq_merged, features = c("Bhlha15", "Aqp5", "Slc12a2", "Smgc", "Pip", "Amy1"), group.by = "seurat_clusters")
# Mist1 (Bhlha15) is a crucial transcription factor for differentiation and function.
# Aqp5 (water channel) and Nkcc1 (or Slc12a2, an ion transporter) are commonly used to identify all major acinar types.
FeaturePlot(scRNA_seq_merged, features = c("Bhlha15", "Aqp5", "Slc12a2", "Smgc", "Pip"))

# immune cells
VlnPlot(scRNA_seq_merged, features = c("Cd79a", "Cd3d", "Cd3e", "Cd3g", "Adgre1", "S100a8", "S100a9"), group.by = "seurat_clusters")
# Cluster 9 is B cells and cluster 7 is T cells
# Cluster 6 is Macrophages

# ductal cells
VlnPlot(scRNA_seq_merged, features = c("Dcpp1", "Dcpp2", "Dcpp3", "Fxyd2"), group.by = "seurat_clusters")

# myoepithelial cells
VlnPlot(scRNA_seq_merged, features = c("Acta2", "Krt14", "Krt5"), group.by = "seurat_clusters")
# clusters 5, 12, and 16 are potentially myoepithelial cells

# stromal and fibroblasts
VlnPlot(scRNA_seq_merged, features = c("Col1a1", "Vim", "Pdgfra"), group.by = "seurat_clusters")
# Cluster 13 is most probably stromal, other clusters need pathologist feedback

# endothelial cells
VlnPlot(scRNA_seq_merged, features = c("Pecam1", "Cd31", "Ackr1"), group.by = "seurat_clusters")
# clusters 2 and 17 are potentially endothelial cells

# We can summarize the above results in a dotplot for all the biological markers
DotPlot(scRNA_seq_merged, features = c("Epcam", "Smgc", "Aqp5", "Bhlha15", "Prol1", "Pip", "Lpo", "Amy1", "Alcam", "Krt19",
                                         "Klk1", "Cftr", "Ascl3", "Krt14", "Krt5", "Acta2", "Cnn1", "Cd79a", "Icos", "Cd68",
                                         "Pecam1", "Tubb3", "Ncam1", "Vim", "Col1a1", "Twist1", "Adgre1")) + RotatedAxis()




###########################################
# Pathway Analysis
###########################################

NOTCH_Pathway_Mouse <- c(
  # Receptors
  "Notch1", "Notch2", "Notch3", "Notch4",
  # Ligands
  "Dll1", "Dll3", "Dll4", "Jag1", "Jag2",
  # Proteolytic Processing Components (γ-secretase complex & ADAMs)
  "Adam10", "Adam17", "Psen1", "Psen2", 
  "Ncstn", "Aph1a", "Aph1b", "Psenen",
  # Transcription Complex
  "Rbpj", "Maml1", "Maml2", "Maml3",
  # Canonical Target Genes
  "Hes1", "Hes5", "Hes6", "Hes7",
  "Hey1", "Hey2", "Heyl",
  # Additional Targets
  "Myc", "Ccnd1", "Nrarp",
  "Dtx1", "Dtx3l",
  "Il7r", "Bcl2",
  # Regulators / Modulators
  "Lfng", "Mfng", "Rfng",
  "Numb", "Numbl",
  "Fbxw7", "Itch")

# which NOTCH pathway genes are markers?
Markers[rownames(Markers) %in% NOTCH_Pathway_Mouse,]
# plot the NOTCH pathway genes and split by status
DotPlot(object = scRNA_seq_merged, features = BMP_Pathway_Mouse, split.by = "Status")


IFN_Pathway_Mouse <- c(
  # Interferon ligands
  "Ifna1", "Ifna2", "Ifna4", "Ifna5", "Ifna6",
  "Ifna9", "Ifna11", "Ifna12", "Ifna13",
  "Ifnb1",
  "Ifng",
  # Receptors
  "Ifnar1", "Ifnar2",
  "Ifngr1", "Ifngr2",
  # JAK-STAT core signaling
  "Jak1", "Jak2", "Tyk2",
  "Stat1", "Stat2", "Stat3",
  "Irf9",
  # Interferon Regulatory Factors
  "Irf1", "Irf2", "Irf3", "Irf5", "Irf7", "Irf8",
  # Canonical ISGs (Interferon-Stimulated Genes)
  "Isg15", "Isg20",
  "Mx1", "Mx2",
  "Oas1a", "Oas1b", "Oas2", "Oas3",
  "Ifit1", "Ifit2", "Ifit3",
  "Rsad2",        # Viperin
  "Ddx58",        # RIG-I
  "Ifih1",        # MDA5
  "Gbp2", "Gbp5",
  "Cxcl10", "Cxcl9",
  "Bst2",
  "Tap1", "Tap2",
  "Psmb8", "Psmb9",
  "Socs1", "Socs3")

BMP_Pathway_Mouse <- c(
  # BMP Ligands
  "Bmp2", "Bmp3", "Bmp4", "Bmp5", "Bmp6", "Bmp7", "Bmp8a", "Bmp8b",
  "Gdf2",   # BMP9
  "Gdf5", "Gdf6", "Gdf7",
  # Type I Receptors
  "Bmpr1a",   # Alk3
  "Bmpr1b",   # Alk6
  "Acvr1",    # Alk2
  # Type II Receptors
  "Bmpr2",
  "Acvr2a",
  "Acvr2b",
  # SMAD Signaling Core
  "Smad1",
  "Smad5",
  "Smad9",    # Smad8
  "Smad4",
  # Inhibitory SMADs
  "Smad6",
  "Smad7",
  # BMP Antagonists (Extracellular Regulators)
  "Nog",
  "Chrd",
  "Chrdl1",
  "Grem1",
  "Grem2",
  "Fst",
  "Fstl1",
  "Sost",
  # Canonical BMP Target Genes
  "Id1", "Id2", "Id3",
  "Runx2",
  "Msx1", "Msx2",
  "Dlx5",
  "Hey1", "Hey2",
  "Junb")


# Can we find the markers (DEGs) between Control and Treatment for each cluster?
table(scRNA_seq_merged@active.ident, scRNA_seq_merged$Status)
# cluster 17 has less than 3 cells for the treatment status. Remove it from the below analysis to avoid errors.

Status_Markers <- lapply(levels(scRNA_seq_merged$seurat_clusters)[1:17], function(clust){
  cluster_obj <- subset(scRNA_seq_merged, idents = clust)
  Idents(cluster_obj) <- "Status"
  markers <- FindMarkers(cluster_obj, ident.1 = "Treatment", ident.2 = "Control",
                         logfc.threshold = 0.25, min.pct = 0.05)
  markers <- markers %>%
    filter(p_val_adj < 0.05)
  notch_markers <- markers[rownames(markers) %in% NOTCH_Pathway_Mouse,]
  bmp_markers <- markers[rownames(markers) %in% BMP_Pathway_Mouse,]
  return(list(markers=markers, NOTCH = notch_markers, BMP=bmp_markers))
})





