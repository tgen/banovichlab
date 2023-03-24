# ==============================================================================
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: Dec 2021
# Description: Cell type annotation for the ILD data
# ==============================================================================

# ======================================
# Load libraries
# ======================================

library(data.table)
library(googlesheets4)
library(RCurl)
library(qdap)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ade4)
library(Matrix)
library(DoubletFinder)
library(tidyverse)
library(harmony)
library(findPC)

# ======================================
# Environment parameters
# ======================================

getwd()
#setwd("/scratch/hnatri/ILD/Seurat_objects")

options(future.globals.maxSize = 4096*1024^2 )
set.seed(2811)

# ======================================
# Helper functions
# ======================================

source("utilities.R")
source("lung_celltype_markers.R")
source("lung_celltype_colors.R")

# ======================================
# Import data
# ======================================

ild_all <- readRDS("/scratch/hnatri/ILD/Seurat_objects/ild_all.6_12_18_24_integrated.clusters.4pop.210704.rds")

# =====================================
# Epithelial cell population
# =====================================

# Split out the epithelial pop
epi <- subset(ild_all,
              cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Epithelial",]))

DefaultAssay(epi) <- "integrated"
# Finding variable features, rescaling and running PCA
# nfeatures = 3000
epi <- FindVariableFeatures(epi, assay = "integrated")
epi.features <- VariableFeatures(epi, assay = "integrated")
epi <- ScaleData(epi, assay = "integrated", verbose = T)
epi <- RunPCA(epi,
              reduction.name = "integrated_pca",
              assay = "integrated",
              features = epi.integrated.features,
              verbose = T)
best_pcs_epi <- get_pcs(epi.features, reduction_name="integrated_pca")
ElbowPlot(epi, reduction="integrated_pca", ndims=50)

# n.neighbors = 50, min.dist = 0.5, spread = 1, 
# Constructing the UMAP
epi <- RunUMAP(epi,
               dims = 1:18,
               reduction = "integrated_pca",
               reduction.name = "integrated_umap",
               assay = "integrated",
               verbose = T)
epi <- FindNeighbors(epi,
                     dims = 1:18,
                     reduction = "integrated_pca",
                     assay = "integrated")

epi <- FindClusters(
    object = epi,
    resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
    graph.name = "integrated_snn"
)

epi@meta.data$seurat_clusters <- epi@meta.data$integrated_snn_res.1.5
Idents(epi) <- "seurat_clusters"

# Add in the label for doublets
# Looking for clusters that express more than one cell population marker.
# PTPRC+ (immune cells), EPCAM+ (epithelial cells), PECAM1+/PTPRC− 
# (endothelial cells), and PTPRC−/EPCAM−/PECAM1− (mesenchymal cells)
DimPlot(epi, reduction = "integrated_umap", group.by = "seurat_clusters")
DotPlot(epi, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN"), assay="RNA")
FeaturePlot(epi, reduction = "integrated_umap", features=c("EPCAM","PECAM1","PTPRC","LUM","DCN"), assay="RNA")

doublets <- ifelse(epi@meta.data$seurat_clusters %in% c(21,24,25,26,29,30,32,35,41), "Doublets", "Singlets")
epi@meta.data$doublet <- doublets
DimPlot(epi, reduction="integrated_umap", group.by = "doublet")

# Remove doublets, rescaling and reconstructing UMAP
epi2 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$doublet != "Doublets", ]))
DefaultAssay(epi2) <- "integrated"
epi2 <- FindVariableFeatures(epi2, assay = "integrated")
epi2.features <- VariableFeatures(epi2, assay = "integrated")
epi2 <- ScaleData(epi2, assay = "integrated", verbose = T)
epi2 <- RunPCA(epi2,
               reduction.name = "integrated_pca",
               assay = "integrated",
               verbose = T,
               features = epi2.features)
best_pcs_epi2 <- get_pcs(epi2, reduction_name="integrated_pca")

epi2 <- RunUMAP(epi2,
                dims = 1:14,
                reduction = "integrated_pca",
                reduction.name = "integrated_umap",
                assay = "integrated",
                verbose = T)
epi2 <- FindNeighbors(epi2,
                      dims = 1:14,
                      reduction = "integrated_pca")
epi2 <- FindClusters(epi2,
                     resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                     graph.name = "integrated_snn")
epi2@meta.data$seurat_clusters <- epi2@meta.data$integrated_snn_res.1.5
Idents(epi2) <- "seurat_clusters"

epi <- epi2

# Run FindtransferAnchors
epi.ref <- readRDS("/scratch/hnatri/ILD/Seurat_objects/Sci_Adv/sci_adv_epi_reannotate.rds")

best_pcs_epi.ref <- get_pcs(epi.ref)

# Using the dimensionality of the reference object
epi.anchors <- FindTransferAnchors(reference = epi.ref,
                                   query = epi,
                                   dims = 1:16,
                                   normalization.method = "SCT",
                                   reference.assay = "SCT",
                                   query.assay = "RNA",
                                   reduction = "pcaproject")

epi.predictions <- TransferData(anchorset = epi.anchors,
                                refdata = epi.ref@meta.data$celltype,
                                dims = 1:16,
                                verbose = T)

epi <- AddMetaData(epi, metadata = epi.predictions)

# NoLegend() + 
p1 <- DimPlot(epi.ref, group.by = "celltype", label = T, repel = T, cols = epi_col) + 
    ggtitle("Reference") +
    theme(legend.position = "none")
p2 <- DimPlot(epi,  reduction="integrated_umap", group.by = "predicted.id", label = T, repel = T, cols = epi_col) + 
    ggtitle("Query") +
    theme(legend.position = "none")

p1 + p2

# Manual annotation of cell types based on marker gene expression
DefaultAssay(epi) <- "SCT"

# AT1 markers
DotPlot(epi, features = c("EPCAM", at1_markers))
FeaturePlot(epi, features = c("EPCAM", at1_markers))

# AT2 markers
DotPlot(epi, features = c("EPCAM", at2_markers))
FeaturePlot(epi, features = c("EPCAM", at2_markers))

# Secretory markers
DotPlot(epi, features = c("EPCAM", secretory_markers))
FeaturePlot(epi, features = c("EPCAM", secretory_markers))

# Basal markers
DotPlot(epi, features = c("EPCAM", basal_markers))
FeaturePlot(epi, features = c("EPCAM", basal_markers))

# Ciliated markers
DotPlot(epi, features = c("EPCAM", ciliated_markers))
FeaturePlot(epi, features = c("EPCAM", ciliated_markers))

# Proliferation markers
DotPlot(epi, features = c("EPCAM", proliferating_markers))
FeaturePlot(epi, features = c("EPCAM", proliferating_markers))
# Proliferating: 34, 23. 40: proliferating or transitional AT2?

# PNEC: SYP, CHGA
FeaturePlot(epi, features = c("EPCAM", "SYP", "CHGA"))

# Annotations
manual_annotation_1 <- as.numeric(as.character(epi@meta.data$seurat_clusters))
manual_annotation_1[manual_annotation_1 %in% c(16)] <- "AT1"
manual_annotation_1[manual_annotation_1 %in% c(23, 29, 21, 1, 7, 3, 30, 32)] <- "AT2" # Remove 32, 30. 29? 23? 21?
manual_annotation_1[manual_annotation_1 %in% c(17)] <- "Transitional AT2" # 0 is secretory or transitional AT2
manual_annotation_1[manual_annotation_1 %in% c(0, 27)] <- "Secretory - SCGB3A2+"
manual_annotation_1[manual_annotation_1 %in% c(10)] <- "Secretory - SCGB1A1+/SCGB3A2+"
manual_annotation_1[manual_annotation_1 %in% c(11, 8, 20)] <- "Secretory - SCGB1A1+/MUC5B+"
manual_annotation_1[manual_annotation_1 %in% c(14, 22)] <- "Basal"
manual_annotation_1[manual_annotation_1 %in% c(2, 4, 25, 13, 6, 19, 12, 5, 15, 18, 33)] <- "Ciliated"
manual_annotation_1[manual_annotation_1 %in% c(9, 31)] <- "Differentiating Ciliated"
manual_annotation_1[manual_annotation_1 %in% c(26)] <- "KRT5-/KRT17+"
manual_annotation_1[manual_annotation_1 %in% c(24, 28)] <- "Proliferating"
manual_annotation_1[manual_annotation_1 %in% c(34)] <- "PNEC"

epi@meta.data$manual_annotation_1 <- manual_annotation_1

DimPlot(epi, group.by = "predicted.id", label = T, repel = T, cols = epi_col) +
    theme(legend.position = "none")
DimPlot(epi, group.by = "manual_annotation_1", label = T, repel = T, cols = epi_col) +
    theme(legend.position = "none")

saveRDS(epi, file = "/scratch/hnatri/ILD/Seurat_objects/epithelial.rds")

# =====================================
# Immune cell population
# ====================================

# Split out the immune cell population
immune <- subset(ild_all, 
                 cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Immune",]))

DefaultAssay(immune) <- "integrated"

# Finding variable features, rescaling and running PCA
immune <- FindVariableFeatures(immune, assay = "integrated")
immune.features <- VariableFeatures(immune,  assay = "integrated")
immune <- ScaleData(immune, assay = "integrated", verbose = T)
immune <- RunPCA(immune,
                 reduction.name = "integrated_pca",
                 assay = "integrated",
                 features = immune.features,
                 verbose = T)
best_pcs_immune.integrated <- get_pcs(immune, reduction_name="integrated_pca")

# How many PCs to include?
ElbowPlot(immune.integrated, ndims=50)

# Constructing the UMAP
immune <- RunUMAP(immune, dims = 1:15,
                  reduction = "integrated_pca",
                  reduction.name = "integrated_umap",
                  assay = "integrated",
                  verbose = T)
immune <- FindNeighbors(immune,
                        dims = 1:15,
                        reduction = "integrated_pca",
                        assay = "integrated")
immune <- FindClusters(object = immune,
                       resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                       graph.name = "integrated_snn")

immune@meta.data$seurat_clusters <- immune@meta.data$integrated_snn_res.2
Idents(immune) <- "integrated_snn_res.2"

# Look for doublets
DimPlot(immune, reduction="integrated_umap", group.by="integrated_snn_res.2")
DotPlot(immune, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN", "PCNA", "TOP2A", "MCM6", "MKI67"), assay="RNA")
FeaturePlot(immune, c("EPCAM","PECAM1","PTPRC","LUM","DCN", "PCNA", "TOP2A", "MCM6", "MKI67"))

doublets <- ifelse(immune@meta.data$seurat_clusters %in% c(27,36,44,46), "Doublets", "Singlets")
immune@meta.data$doublet <- doublets

# Remove doublets 
immune2 <- subset(immune, cells = rownames(immune@meta.data[immune@meta.data$doublet != "Doublets", ]))

DefaultAssay(immune2) <- "integrated"
immune2 <- FindVariableFeatures(immune2, assay = "integrated")
immune2.features <- VariableFeatures(immune2, assay = "integrated")
immune2 <- ScaleData(immune2, assay = "integrated", verbose = T)
immune2 <- RunPCA(immune2,
                  reduction.name = "integrated_pca",
                  assay = "integrated",
                  verbose = T,
                  features = immune2.features)

best_pcs_immune2 <- get_pcs(immune2, reduction_name="integrated_pca")

immune2 <- RunUMAP(immune2,
                   dims = 1:14,
                   reduction = "integrated_pca",
                   reduction.name = "integrated_umap",
                   assay = "integrated",
                   verbose = T)
immune2 <- FindNeighbors(immune2,
                         dims = 1:14,
                         reduction = "integrated_pca")
immune2 <- FindClusters(immune2,
                        resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                        graph.name = "integrated_snn")

immune2@meta.data$seurat_clusters <- immune2@meta.data$integrated_snn_res.2
Idents(immune2) <- "seurat_clusters"

DotPlot(immune2, features=c("EPCAM","PECAM1","PTPRC","LUM","DCN", "PCNA", "TOP2A", "MCM6", "MKI67"), assay="RNA")
FeaturePlot(immune2, reduction="integrated_umap", c("EPCAM","PECAM1","PTPRC","LUM","DCN"))

immune <- immune2

immune.ref <- readRDS("/scratch/hnatri/ILD/Seurat_objects/Sci_Adv/sci_adv_immune_reannotate.rds")
immune.ref <- UpdateSeuratObject(immune.ref)
immune.ref <- UpdateSCTAssays(immune.ref)

# How many PCs to include?
ElbowPlot(immune.ref, ndims=50)
best_pcs_immune.ref <- get_pcs(immune.ref)

# Need the SCT assay for FindTransferAnchors
DefaultAssay(immune) <- "RNA"
immune <- SCTransform(immune)

DefaultAssay(immune.ref) <- "RNA"
immune.ref[["SCT"]] <- NULL
immune.ref.list <- SplitObject(object = immune.ref, split.by = "orig.ident")
for (i in 1:length(immune.ref.list)){
    immune.ref.list[[i]] <- SCTransform(immune.ref.list[[i]], method = "glmGamPoi" )
}

immune.ref.features <- SelectIntegrationFeatures(object.list = immune.ref.list, nfeatures = 3000)
immune.ref.merge <- merge(x = immune.ref.list[[1]], y = immune.ref.list[2:length(immune.ref.list)])
VariableFeatures(immune.ref.merge[["SCT"]]) <- immune.ref.features
immune.ref.merge <- RunPCA(immune.ref.merge) %>% RunUMAP(dims = 1:15)

# dims should be the dimensionality of the reference dataset.
# Note: unable to use integrated or RNA assays as query.assay.
immune.anchors <- FindTransferAnchors(reference = immune.ref.merge,
                                      query = immune,
                                      dims = 1:15,
                                      normalization.method = "SCT",
                                      reference.assay = "SCT",
                                      query.assay = "SCT",
                                      reduction = "pcaproject",
                                      #features = common_features,
                                      recompute.residuals = FALSE)

immune.predictions <- TransferData(anchorset = immune.anchors,
                                   refdata = immune.ref.merge@meta.data$celltype,
                                   dims = 1:15,
                                   verbose = T)

immune <- AddMetaData(immune, metadata = immune.predictions)

p1 <- DimPlot(immune.ref, group.by = "celltype", label = T, repel = T, cols = immune_col) + 
    NoLegend() + ggtitle ("Reference")
p2 <- DimPlot(immune, reduction="integrated_umap", group.by = "predicted.id", label = T, repel = T, cols = immune_col) + 
    NoLegend() + ggtitle ("Query")

p1 + p2

# Manual annotations based on marker gene expression

# B cell markers
DotPlot(immune, features=c("PTPRC", bcell_markers))
FeaturePlot(immune, features=bcell_markers)

# Plasma cell markers
DotPlot(immune, features=c("PTPRC", plasma_markers))
FeaturePlot(immune, features=plasma_markers)

# Mast cell markers
DotPlot(immune, features=c("PTPRC", mast_markers))
FeaturePlot(immune, features=mast_markers)

# Monocyte markers
# Inflammatory monocytes: IL-6, IL-8, CCL2, CCL3, and CCL5. CCR2, GR1
DotPlot(immune, features=c("PTPRC", monocyte_markers))
FeaturePlot(immune, features=monocyte_markers)

# cDC markers
# cDC1: CD11C (ITGAX)
# cDC2: IRF4
DotPlot(immune, features=c("PTPRC", cdc_markers))
FeaturePlot(immune, features=cdc_markers)

# pDC markers
DotPlot(immune, features=c("PTPRC", pdc_markers))
FeaturePlot(immune, features=pdc_markers)

# Macrophage markers
DotPlot(immune, features=c("PTPRC", macrophage_markers))
FeaturePlot(immune, features=macrophage_markers)

# Proliferating Macrophages
DotPlot(immune, features=c("PTPRC", prolif_macrophage_markers))
FeaturePlot(immune, features=prolif_macrophage_markers)

# NK cell markers. NKG7 (high), CD8A-
DotPlot(immune, features=c("PTPRC", nk_markers))
FeaturePlot(immune, features=nk_markers)

# T cell markers
DotPlot(immune, features=c("PTPRC", tcell_markers))
FeaturePlot(immune, features=tcell_markers)

# Proliferating T Cells
DotPlot(immune, features=c("PTPRC", prolif_tcell_markers))
FeaturePlot(immune, features=prolif_tcell_markers)

# TReg markers
DotPlot(immune, features=c("PTPRC", treg_markers))
FeaturePlot(immune, features=treg_markers)

# CD8 T cells
DotPlot(immune, features=c("PTPRC", cd8_tcell_markers))
FeaturePlot(immune, features=cd8_tcell_markers)

# CD4 T cells, CD8A-
DotPlot(immune, features=c("PTPRC", cd4_tcell_markers))
FeaturePlot(immune, features=cd4_tcell_markers)

# Annotations
manual_annotation_1 <- as.numeric(as.character(immune@meta.data$seurat_clusters))
manual_annotation_1[manual_annotation_1 %in% c(41, 35)] <- "B cells"
manual_annotation_1[manual_annotation_1 %in% c(40, 46)] <- "Plasma" # 46?
manual_annotation_1[manual_annotation_1 %in% c(32, 45)] <- "Mast"
manual_annotation_1[manual_annotation_1 %in% c(6)] <- "Monocyte"
manual_annotation_1[manual_annotation_1 %in% c(23, 22, 10, 20, 29, 34, 39, 19)] <- "Inflammatory monocyte"
manual_annotation_1[manual_annotation_1 %in% c(11)] <- "cDC1" # cDC markers
manual_annotation_1[manual_annotation_1 %in% c(43)] <- "cDC2" #?
manual_annotation_1[manual_annotation_1 %in% c(42)] <- "pDC"
manual_annotation_1[manual_annotation_1 %in% c(8, 31)] <- "moDC"
manual_annotation_1[manual_annotation_1 %in% c(2, 26, 17, 14, 1, 4)] <- "Alveolar macrophage"
manual_annotation_1[manual_annotation_1 %in% c(28, 3, 18, 25, 36, 5, 0, 16, 21)] <- "Monocyte-derived macrophage"
manual_annotation_1[manual_annotation_1 %in% c(9, 44)] <- "Macrophage - SPP1+"
manual_annotation_1[manual_annotation_1 %in% c(27, 38, 37)] <- "Proliferating" # 39?
manual_annotation_1[manual_annotation_1 %in% c(7)] <- "NK" # !
manual_annotation_1[manual_annotation_1 %in% c(13, 24)] <- "CD8/NKT"
manual_annotation_1[manual_annotation_1 %in% c(12, 15, 33, 30)] <- "CD4"
#manual_annotation_1[manual_annotation_1 %in% c(29)] <- "unclear"

immune@meta.data$manual_annotation_1 <- manual_annotation_1

# Saving the object
saveRDS(immune, file = "/scratch/hnatri/ILD/Seurat_objects/immune.rds")

# =====================================
# Mesenchymal cell population
# ====================================

# Split out the mesenchymal population
mesen <- subset(ild_all,
                cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Mesenchymal",]))

DefaultAssay(mesen) <- "integrated"

# Finding variable features, rescaling and running PCA
mesen <- FindVariableFeatures(mesen, assay = "integrated")
mesen.features <- VariableFeatures(mesen.integrated, assay = "integrated")
mesen <- ScaleData(mesen, assay = "integrated", verbose = T)
mesen <- RunPCA(mesen,
                reduction.name = "integrated_pca",
                assay = "integrated",
                features = mesen.integrated.features,
                verbose = T)

# How many PCs to include?
ElbowPlot(mesen.integrated, ndims=50)
get_pcs(mesen, reduction_name="integrated_pca")

# Constructing the UMAP
mesen <- RunUMAP(mesen, dims = 1:25,
                 reduction = "integrated_pca",
                 reduction.name = "integrated_umap",
                 assay = "integrated",
                 verbose = T)
mesen <- FindNeighbors(mesen,
                       dims = 1:25,
                       reduction = "integrated_pca",
                       assay = "integrated")
mesen <- FindClusters(object = mesen,
                      resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                      graph.name = "integrated_snn")

mesen@meta.data$seurat_clusters <- mesen@meta.data$integrated_snn_res.2
Idents(mesen) <- "integrated_snn_res.2"

# Look for doublets
DefaultAssay(mesen) <- "SCT"
DotPlot(mesen, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN","COL1A1","ACTA2","COL8A1","WNT2","MYLK"))
FeaturePlot(mesen, c("EPCAM","PECAM1","PTPRC","LUM","DCN","COL1A1","ACTA2","COL8A1","WNT2","MYLK"))

DimPlot(mesen, group.by = "seurat_clusters", label=T, repel=T) + NoLegend()
doublets <- ifelse(mesen@meta.data$seurat_clusters %in% c(11, 13, 15, 16, 17, 18, 21, 22, 25, 26, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44), "Doublets", "Singlets")
mesen@meta.data$doublet <- doublets

# Remove doublets
DefaultAssay(mesen) <- "integrated"
mesen2 <- subset(mesen, cells = rownames(mesen@meta.data[mesen@meta.data$doublet != "Doublets", ]))
DimPlot(mesen2, group.by = "seurat_clusters", label=T, repel=T) + NoLegend()
mesen2 <- FindVariableFeatures(mesen2, assay = "integrated")
mesen2.features <- VariableFeatures(mesen2)
mesen2 <- ScaleData(mesen2, assay = "integrated")
mesen2 <- RunPCA(mesen2,
                 reduction.name = "integrated_pca",
                 assay = "integrated",
                 features = mesen2.features,
                 verbose = T)

# How many PCs to include?
ElbowPlot(mesen2, ndims=50)
get_pcs(mesen2, reduction_name="integrated_pca")

# Constructing the UMAP
mesen2 <- RunUMAP(mesen2,
                  dims = 1:16,
                  reduction = "integrated_pca",
                  reduction.name = "integrated_umap",
                  assay = "integrated",
                  verbose = T)
mesen2 <- FindNeighbors(mesen2,
                        dims = 1:16,
                        reduction = "integrated_pca",
                        assay = "integrated")
mesen2 <- FindClusters(object = mesen2,
                       resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                       graph.name = "integrated_snn")

mesen2@meta.data$seurat_clusters <- mesen2@meta.data$integrated_snn_res.1
Idents(mesen2) <- "integrated_snn_res.1"

mesen <- mesen2

# Run FindtransferAnchors
mesen.ref <- readRDS("/scratch/hnatri/ILD/Seurat_objects/Sci_Adv/mesenchymal_cleaned_annotated_reclustered_20200220.rds")
mesen.ref <- UpdateSeuratObject(mesen.ref)
get_pcs(mesen.ref)

mesen.anchors <- FindTransferAnchors(reference = mesen.ref,
                                     query = mesen,
                                     dims = 1:19,
                                     normalization.method = "SCT",
                                     reference.assay = "SCT",
                                     query.assay = "RNA",
                                     reduction = "pcaproject")

mesen.predictions <- TransferData(anchorset = mesen.anchors,
                                  refdata = mesen.ref@meta.data$celltypes,
                                  dims = 1:19,
                                  verbose = T)

mesen <- AddMetaData(mesen, metadata = mesen.predictions)

p1 <- DimPlot(mesen.ref, group.by = "celltypes", label = T, repel = T, cols = mesen_col) + 
    NoLegend() + ggtitle ("Reference")
p2 <- DimPlot(mesen, group.by = "predicted.id", label = T, repel = T, cols = mesen_col) + 
    NoLegend() + ggtitle ("Query")

p1 + p2

# Manual annotations
DefaultAssay(mesen) <- "SCT"
FeaturePlot(mesen, features = mesothelial_markers)
FeaturePlot(mesen, features = smc_markers)
FeaturePlot(mesen, features = pericyte_markers)
FeaturePlot(mesen, features = fibroblast_markers)
FeaturePlot(mesen, features = myofb_markers)
FeaturePlot(mesen, features = activated_myofb_markers)
FeaturePlot(mesen, features = WNT2_fibro_markers)
FeaturePlot(mesen, features = matrixfb_markers)
FeaturePlot(mesen, features = PLIN2_fibro_markers)
FeaturePlot(mesen, features = HAS1_fibro_markers)
FeaturePlot(mesen, features = act_HAS1_fibro_markers)

DimPlot(mesen, group.by = "predicted.id", repel = T, label = T) + NoLegend()

manual_annotation_1 <- as.numeric(as.character(mesen@meta.data$seurat_clusters))

manual_annotation_1[manual_annotation_1 %in% c(3, 7, 16)] <- "Mesothelial"
manual_annotation_1[manual_annotation_1 %in% c(13)] <- "MyoFB"
manual_annotation_1[manual_annotation_1 %in% c()] <- "HAS1 High Activated FB"
manual_annotation_1[manual_annotation_1 %in% c(6, 15)] <- "PLIN2+ FB"
manual_annotation_1[manual_annotation_1 %in% c(2, 14, 17, 11)] <- "WNT2+ FB"
manual_annotation_1[manual_annotation_1 %in% c(5, 9)] <- "Pericyte"
manual_annotation_1[manual_annotation_1 %in% c(1)] <- "SMC"
manual_annotation_1[manual_annotation_1 %in% c(0, 8, 4)] <- "Matrix FB"
manual_annotation_1[manual_annotation_1 %in% c(12, 10, 18)] <- "MyoFB - Activated"

mesen@meta.data$manual_annotation_1 <- manual_annotation_1

saveRDS(mesen, "/scratch/hnatri/ILD/Seurat_objects/mesenchymal.rds")

# =====================================
# Endothelial cell population
# =====================================

# Split out the endothelial population
endo <- subset(ild_all,
               cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Endothelial",]))

DefaultAssay(endo) <- "integrated"

# Finding variable features, rescaling and running PCA
endo <- FindVariableFeatures(endo, assay = "integrated")
endo.features <- VariableFeatures(endo, assay = "integrated")
endo <- ScaleData(endo, assay = "integrated", verbose = T)
endo <- RunPCA(endo,
               reduction.name = "integrated_pca",
               assay = "integrated",
               features = endo.integrated.features,
               verbose = T)

# How many PCs to include?
ElbowPlot(endo.integrated, ndims=50)
get_pcs(endo.integrated, reduction_name="integrated_pca")

# Constructing the UMAP
endo <- RunUMAP(endo, dims = 1:17,
                reduction = "integrated_pca",
                reduction.name = "integrated_umap",
                assay = "integrated",
                verbose = T)
endo <- FindNeighbors(endo,
                      dims = 1:17,
                      reduction = "integrated_pca",
                      assay = "integrated")
endo.clusters <- FindClusters(object = endo,
                              resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                              graph.name = "integrated_snn")

endo@meta.data$seurat_clusters <- endo@meta.data$integrated_snn_res.1
Idents(endo) <- "integrated_snn_res.1"

# Mark doublets
DefaultAssay(endo) <- "SCT"
DotPlot(endo, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN"))
FeaturePlot(endo, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN"))
doublets <- as.character(endo@meta.data$seurat_clusters)
doublets <- ifelse(endo@meta.data$seurat_clusters %in% c(8,10,11,14,15,16,17,20,21,23,24), "Doublets", "Singlets")
endo@meta.data$doublet <- doublets

# Remove doublets
DefaultAssay(endo) <- "integrated"
endo2 <- subset(endo, cells = rownames(endo@meta.data[endo@meta.data$doublet != "Doublets", ]))
DimPlot(endo2, group.by = "seurat_clusters", label=T, repel=T) + NoLegend()
endo2 <- FindVariableFeatures(endo2, assay = "integrated")
endo2.features <- VariableFeatures(endo2)
endo2 <- ScaleData(endo2, assay = "integrated")
endo2 <- RunPCA(endo2,
                reduction.name = "integrated_pca",
                assay = "integrated",
                features = mesen2.features,
                verbose = T)

# How many PCs to include?
ElbowPlot(endo2, ndims=50)
get_pcs(endo2, reduction_name="integrated_pca")

# Constructing the UMAP
endo2 <- RunUMAP(endo2,
                 dims = 1:13,
                 reduction = "integrated_pca",
                 reduction.name = "integrated_umap",
                 assay = "integrated",
                 verbose = T)
endo2 <- FindNeighbors(endo2,
                       dims = 1:13,
                       reduction = "integrated_pca",
                       assay = "integrated")
endo2 <- FindClusters(object = endo2,
                      resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                      graph.name = "integrated_snn")

endo2@meta.data$seurat_clusters <- endo2@meta.data$integrated_snn_res.1
Idents(endo2) <- "integrated_snn_res.1"

# Prepare the reference file
endo.ref <- readRDS("/scratch/hnatri/ILD/Seurat_objects/Sci_Adv/sci_adv_stromal_reannotate.rds")
endo.ref <- UpdateSeuratObject(endo.ref)
endo.ref <- UpdateSCTAssays(endo.ref)

# Removing mesenchymal cell types
endo.cells <- c("Endothelial - capillary","Endothelial - venule","Endothelial - arteriole", 
                "Lymphatic", "Endothelial - inflamed","Endothelial - CA4+ capillary", 
                "Endothelial - peribronchiolar")
endo.ref <- subset(endo.ref, cells=rownames(endo.ref@meta.data[endo.ref@meta.data$celltype %in% endo.cells,]))
endo.ref <- SCTransform(endo.ref)
endo.ref <- RunPCA(endo.ref)

get_pcs(endo.ref)
endo.ref <- RunUMAP(endo.ref, dims = 1:17)
endo <- endo2

# Run FindtransferAnchors
endo.anchors <- FindTransferAnchors(reference = endo.ref,
                                    query = endo,
                                    dims = 1:17,
                                    normalization.method = "SCT",
                                    reference.assay = "SCT",
                                    query.assay = "RNA",
                                    reduction = "pcaproject")

endo.predictions <- TransferData(anchorset = endo.anchors,
                                 refdata = endo.ref@meta.data$celltype,
                                 dims = 1:17,
                                 verbose = T)

endo <- AddMetaData(endo, metadata = endo.predictions)

p1 <- DimPlot(endo.ref, group.by = "celltype", label = T, repel = T, cols = endo_col) + 
    NoLegend() + ggtitle("Reference")
p2 <- DimPlot(endo, group.by = "predicted.id", label = T, repel = T, cols = endo_col) + 
    NoLegend() + ggtitle("Query")

p1 + p2

# Manual annotations
endo_markers <- c("PECAM1", "CLDN5", "PTPRC")
FeaturePlot(endo, features = endo_markers)

vascular_markers <- c("VWF")
FeaturePlot(endo, features = vascular_markers)

capillary_markers <- c("CA4", "VIPR1", "RGCC", "CYB5A", "ADGRL2")
FeaturePlot(endo, features = capillary_markers)

g1_markers <- c("VWF", "HPGD", "EDNRB", "EMCN")
FeaturePlot(endo, features = g1_markers)

g2_markers <- c("VWF", "FCN3", "IL7R", "SLC6A4")
FeaturePlot(endo, features = g2_markers)

arterial_markers <- c("DKK2", "GJA5", "BMX", "HEY1")
FeaturePlot(endo, features = arterial_markers)

venous_markers <- c("PLA1A", "CPE", "PTGDS", "ACKR1", "SPRY1")
FeaturePlot(endo, features = venous_markers)

bronchial_markers <- c("SPRY1", "PLVAP", "COL15A1", "VWA1", "MYC")
FeaturePlot(endo, features = bronchial_markers)

bronchial_g1_markers <- c("POSTN", "ACKR1", "SPRY1")
FeaturePlot(endo, features = bronchial_g1_markers)

bronchial_g2_markers <- c("HBEGF", "ACKR1", "SPRY1")
FeaturePlot(endo, features = bronchial_g2_markers)

lymphatic_markers <- c("CCL21", "PROX1", "PDPN")
FeaturePlot(endo, features = lymphatic_markers)

manual_annotation_1 <- as.numeric(as.character(endo@meta.data$seurat_clusters))

manual_annotation_1[manual_annotation_1 %in% c(4)] <- "Endothelial - capillary"
manual_annotation_1[manual_annotation_1 %in% c(1, 10, 11, 14)] <- "Endothelial - venule"
manual_annotation_1[manual_annotation_1 %in% c(6, 7, 17, 15, 16, 26)] <- "Endothelial - arteriole"
manual_annotation_1[manual_annotation_1 %in% c(12, 2)] <- "Lymphatic"
manual_annotation_1[manual_annotation_1 %in% c(23)] <- "Endothelial - inflamed" # 23?
manual_annotation_1[manual_annotation_1 %in% c(5, 8, 20)] <- "Endothelial - CA4+ capillary"
manual_annotation_1[manual_annotation_1 %in% c(0, 3, 9, 10, 25, 22)] <- "Endothelial - peribronchiolar"

endo@meta.data$manual_annotation_1 <- manual_annotation_1

saveRDS(endo, file = "/scratch/hnatri/ILD/Seurat_objects/endothelial.rds")

# =====================================================
# Merge objects
# ====================================================

epi <- readRDS("/scratch/hnatri/ILD/Seurat_objects/epithelial.rds")
immune <- readRDS("/scratch/hnatri/ILD/Seurat_objects/immune.rds")
mesen <- readRDS("/scratch/hnatri/ILD/Seurat_objects/mesenchymal.rds")
endo <- readRDS("/scratch/hnatri/ILD/Seurat_objects/endothelial.rds")

# Merge all cell population objects
ild <- merge(x=epi, y=c(mesen, endo, immune))

# Some gene_names have extra .1 at the end, dropping those
counts <- GetAssayData(epi, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% setdiff(rownames(epi), gtf_df_genes$gene_name))),]

DefaultAssay(ild) <- "RNA"
ild[["SCT"]] <- NULL
ild <- FindVariableFeatures(ild)
ild.features <- VariableFeatures(ild)

# Rerun SCTransform for ILD object
# SCTransform took about 3h
ild <- SCTransform(ild)
ild <- RunPCA(ild)
get_pcs(ild)
ild <- RunUMAP(ild, dims=1:15)

# Saving the object
saveRDS(ild, file = "/scratch/hnatri/ILD/Seurat_objects/ILD_annotated_210819.rds")

