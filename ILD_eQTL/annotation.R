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
epi.integrated <- subset(ild_all, 
                         cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Epithelial",]))

DefaultAssay(epi.integrated) <- "integrated"
# Finding variable features, rescaling and running PCA
# nfeatures = 3000
epi.integrated <- FindVariableFeatures(epi.integrated,
                                       assay = "integrated")
epi.integrated.features <- VariableFeatures(epi.integrated,
                                            assay = "integrated")
epi.integrated <- ScaleData(epi.integrated,
                            assay = "integrated",
                            verbose = T)
epi.integrated <- RunPCA(epi.integrated,
                         reduction.name = "integrated_pca",
                         assay = "integrated",
                         features = epi.integrated.features,
                         verbose = T)
best_pcs_epi.integrated <- get_pcs(epi.integrated, reduction_name="integrated_pca")

# Saving the ElbowPlot
filename <- "/home/hnatri/ILD_processing_annotation/epithelial_elbowplot.pdf"
pdf(file = filename,
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
ElbowPlot(epi.integrated, reduction="integrated_pca", ndims=50)
dev.off()

# n.neighbors = 50, min.dist = 0.5, spread = 1, 
# Constructing the UMAP
epi.integrated <- RunUMAP(epi.integrated,
                          dims = 1:18,
                          reduction = "integrated_pca",
                          reduction.name = "integrated_umap",
                          assay = "integrated",
                          verbose = T)
epi.integrated <- FindNeighbors(epi.integrated,
                                dims = 1:18,
                                reduction = "integrated_pca",
                                assay = "integrated")

epi.integrated.clusters <- FindClusters(
    object = epi.integrated,
    resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
    graph.name = "integrated_snn"
)

saveRDS(epi.integrated.clusters, file = "/scratch/hnatri/ILD/Seurat_objects/epithelial_integrated_dims18_clusters_210803.rds")
epi <- epi.integrated.clusters

epi@meta.data$seurat_clusters <- epi@meta.data$integrated_snn_res.1.5
Idents(epi) <- "seurat_clusters"

# Proportions of doublets in each cluster:
doublet_summary <- epi@meta.data %>% group_by(seurat_clusters, doublet_finder) %>%
    summarize(n = n())

doublet_summary <- pivot_wider(doublet_summary, names_from = doublet_finder, values_from = n) 
doublet_summary$doublet_prop <- doublet_summary$Doublet/(doublet_summary$Doublet+doublet_summary$Singlet)
doublet_summary <- arrange(doublet_summary, desc(doublet_prop))

# Add in the label for doublets
# Looking for clusters that express more than one cell population marker.
# PTPRC+ (immune cells), EPCAM+ (epithelial cells), PECAM1+/PTPRC− 
# (endothelial cells), and PTPRC−/EPCAM−/PECAM1− (mesenchymal cells)
DimPlot(epi, reduction = "integrated_umap", group.by = "seurat_clusters")
DotPlot(epi, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN"), assay="RNA")
FeaturePlot(epi, reduction = "integrated_umap", features=c("EPCAM","PECAM1","PTPRC","LUM","DCN"), assay="RNA")

doublets <- ifelse(epi@meta.data$seurat_clusters %in% c(21,24,25,26,29,30,32,35,41), "Doublets", "Singlets")
epi@meta.data$doublet <- doublets

filename <- "/home/hnatri/ILD_processing_annotation/epithelial_doublet_dimplot.pdf"
pdf(file = filename,
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches
DimPlot(epi, reduction="integrated_umap", group.by = "doublet")
dev.off()

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

saveRDS(epi2, file = "/scratch/hnatri/ILD/Seurat_objects/epithelial_integrated_dims18_clusters_res1.5_woutdoublets_dims14_nneighbors30_mindist0.3_210804.rds")

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

# Make some plots
epi_col <- c("AT1" = "#548BC5",
             "AT2" = "#EE7342",
             "Proliferating" = "#B874FF",
             "Transitional AT2" = "#F1A5C4",
             "AT2 - low quality" = "#2FC895",
             "Basal" = "#D38402",
             "Ciliated" = "#B19302",
             "PNEC" = "#FFA19D",
             "Differentiating Ciliated" = "#9C9966",
             "KRT5-/KRT17+" = "#03AF21",
             "Secretory - MUC5AC+" = "#003366",
             "Secretory - SCGB1A1+/MUC5B+" = "#00B2DB",
             "Secretory - SCGB3A1+/MUC5B+" = "#F659DD",
             "Secretory - SCGB1A1+/SCGB3A2+" = "#9900CC",
             "Secretory - SCGB3A2+" = "#FF5BC8")

# NoLegend() + 
p1 <- DimPlot(epi.ref, group.by = "celltype", label = T, repel = T, cols = epi_col) + 
    ggtitle("Reference") +
    theme(legend.position = "none")
p2 <- DimPlot(epi,  reduction="integrated_umap", group.by = "predicted.id", label = T, repel = T, cols = epi_col) + 
    ggtitle("Query") +
    theme(legend.position = "none")

p1 + p2

saveRDS(epi, file = "/scratch/hnatri/ILD/Seurat_objects/epithelial_integrated_dims18_clusters_res1.5_woutdoublets_dims14_nneighbors30_mindist0.3_annotated_210804.rds")

# Manual annotation of cell types based on marker gene expression
DefaultAssay(epi) <- "SCT"

# AT1 markers
at1_markers <- c("AGER", "PDPN", "CAV1", "EMP2")
DotPlot(epi, features = c("EPCAM", at1_markers))
FeaturePlot(epi, features = c("EPCAM", at1_markers))

# AT2 markers
at2_markers <- c("SFTPC", "ABCA3", "LAMP3", "AGER")
DotPlot(epi, features = c("EPCAM", at2_markers))
FeaturePlot(epi, features = c("EPCAM", at2_markers))

# Secretory markers
secretory_markers <- c("SCGB1A1", "SCGB3A2", "MGP", "MUC5B", "MUC5AC")
DotPlot(epi, features = c("EPCAM", secretory_markers))
FeaturePlot(epi, features = c("EPCAM", secretory_markers))

# Basal markers
basal_markers <- c("KRT5", "KRT17", "CHGA", "CALCA", "FOXI1", "SERPINB3")
DotPlot(epi, features = c("EPCAM", basal_markers))
FeaturePlot(epi, features = c("EPCAM", basal_markers))

# Ciliated markers
ciliated_markers <- c("FOXJ1", "SFTPB", "OMG", "PPIL6", "NME5", "NWD1", "SLAIN2", "SPAG16")
DotPlot(epi, features = c("EPCAM", ciliated_markers))
FeaturePlot(epi, features = c("EPCAM", ciliated_markers))

# Proliferation markers
proliferating_markers <- c("MKI67", "CDK1")
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

saveRDS(epi, file = "/scratch/hnatri/ILD/Seurat_objects/epithelial_integrated_dims18_clusters_res1.5_woutdoublets_dims14_nneighbors30_mindist0.3_annotated_manualannot_210805.rds")

# =====================================
# Immune cell population
# ====================================

# Split out the immune cell population
immune.integrated <- subset(ild_all, 
                            cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Immune",]))

DefaultAssay(immune.integrated) <- "integrated"

# Finding variable features, rescaling and running PCA
immune.integrated <- FindVariableFeatures(immune.integrated,
                                          assay = "integrated")
immune.integrated.features <- VariableFeatures(immune.integrated,
                                               assay = "integrated")
immune.integrated <- ScaleData(immune.integrated,
                               assay = "integrated",
                               verbose = T)
immune.integrated <- RunPCA(immune.integrated,
                            reduction.name = "integrated_pca",
                            assay = "integrated",
                            features = immune.integrated.features,
                            verbose = T)
best_pcs_immune.integrated <- get_pcs(immune.integrated, reduction_name="integrated_pca")

# How many PCs to include?
ElbowPlot(immune.integrated, ndims=50)

# Constructing the UMAP
immune.integrated <- RunUMAP(immune.integrated, dims = 1:15,
                             reduction = "integrated_pca",
                             reduction.name = "integrated_umap",
                             assay = "integrated",
                             verbose = T)
immune.integrated <- FindNeighbors(immune.integrated,
                                   dims = 1:15,
                                   reduction = "integrated_pca",
                                   assay = "integrated")
immune.integrated.clusters <- FindClusters(object = immune.integrated,
                                           resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                                           graph.name = "integrated_snn")

saveRDS(immune.integrated.clusters, file = "/scratch/hnatri/ILD/Seurat_objects/immune_integrated_dims15_clusters_210804.rds")

immune <- immune.integrated.clusters

immune@meta.data$seurat_clusters <- immune@meta.data$integrated_snn_res.2
Idents(immune) <- "integrated_snn_res.2"

# Proportions of doublets in each cluster:
doublet_summary <- immune@meta.data %>% group_by(seurat_clusters, doublet_finder) %>%
    summarize(n = n())

doublet_summary <- pivot_wider(doublet_summary, names_from = doublet_finder, values_from = n) 
doublet_summary$doublet_prop <- doublet_summary$Doublet/(doublet_summary$Doublet+doublet_summary$Singlet)
doublet_summary <- arrange(doublet_summary, desc(doublet_prop))

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

DefaultAssay(immune.ref)
DefaultAssay(immune)

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

immune_col <- c("Proliferating" = "#009FFA",
                "Monocyte-derived macrophage" = "#FF6B65",
                "moDC" = "#A39922",
                "Inflammatory monocyte" = "#CC8B22",
                "Alveolar macrophage" = "#EB7B27",
                "cDC2" = "#FF57C3",
                "cDC1" = "#00B7BB",
                "Macrophage - SPP1+" = "#00AE21",
                "Monocyte" = "#00AFE0",
                "CD4" = "#6CA421",
                "NK" = "#00B990",
                "Mast" = "#C471FA",
                "pDC" = "#FF5D97",
                "CD8/NKT" = "#00B560",
                "Plasma" = "#F35EE6",
                "B cells" = "#748AFA")

p1 <- DimPlot(immune.ref, group.by = "celltype", label = T, repel = T, cols = immune_col) + 
    NoLegend() + ggtitle ("Reference")
p2 <- DimPlot(immune, reduction="integrated_umap", group.by = "predicted.id", label = T, repel = T, cols = immune_col) + 
    NoLegend() + ggtitle ("Query")

p1 + p2

saveRDS(immune, file = "/scratch/hnatri/ILD/Seurat_objects/immune_sct_ncells5k_nfeat3k_dims15_clusters_res2_woutDoublets_dims14_annotated_210804.rds")

# Manual annotations based on marker gene expression

# B cell markers
bcell_markers <- c("MS4A1", "CD19", "CD79A")
DotPlot(immune, features=c("PTPRC", bcell_markers))
FeaturePlot(immune, features=bcell_markers)

# Plasma cell markers
plasma_markers <- c("JCHAIN", "IGHG1", "IGLL5")
DotPlot(immune, features=c("PTPRC", plasma_markers))
FeaturePlot(immune, features=plasma_markers)

# Mast cell markers
mast_markers <- c("CPA3", "KIT")
DotPlot(immune, features=c("PTPRC", mast_markers))
FeaturePlot(immune, features=mast_markers)

# Monocyte markers
# Inflammatory monocytes: IL-6, IL-8, CCL2, CCL3, and CCL5. CCR2, GR1
monocyte_markers <- c("CD14", "CD16", "S100A12", "FCN1", "S100A9", "LYZ", "CD14", "CCL2", "CCL3", "CCL5", "IL6", "IL8", "CCR2", "GR1")
DotPlot(immune, features=c("PTPRC", monocyte_markers))
FeaturePlot(immune, features=monocyte_markers)

# cDC markers
# cDC1: CD11C (ITGAX)
# cDC2: IRF4
cdc_markers <- c("FCER1A", "CD1C", "CLEC9A", "IRF4", "CD11C")
DotPlot(immune, features=c("PTPRC", cdc_markers))
FeaturePlot(immune, features=cdc_markers)

DimPlot(immune, split.by = "predicted.id", ncol = 4) + NoLegend()

# pDC markers
pdc_markers <- c("LILRA4", "CLEC4C", "JCHAIN")
DotPlot(immune, features=c("PTPRC", pdc_markers))
FeaturePlot(immune, features=pdc_markers)

# Macrophage markers
macrophage_markers <- c("LYZ", "MARCO", "FCGR1A", "C1QA", "APOC1", "SPP1")
DotPlot(immune, features=c("PTPRC", macrophage_markers))
FeaturePlot(immune, features=macrophage_markers)

# Proliferating Macrophages
prolif_macrophage_markers <- c("MKI67", "CD1", "LYZ")
DotPlot(immune, features=c("PTPRC", prolif_macrophage_markers))
FeaturePlot(immune, features=prolif_macrophage_markers)

# NK cell markers. NKG7 (high), CD8A-
nk_markers <- c("NCR1", "KLRB1", "NKG7", "CD8A", "GNLY")
DotPlot(immune, features=c("PTPRC", nk_markers))
FeaturePlot(immune, features=nk_markers)

# T cell markers
tcell_markers <- c("CD3E")
DotPlot(immune, features=c("PTPRC", tcell_markers))
FeaturePlot(immune, features=tcell_markers)

# Proliferating T Cells
prolif_tcell_markers <- c("MKI67", "CD1", "CD3E")
DotPlot(immune, features=c("PTPRC", prolif_tcell_markers))
FeaturePlot(immune, features=prolif_tcell_markers)

# TReg markers
treg_markers <- c("FOXP3")
DotPlot(immune, features=c("PTPRC", treg_markers))
FeaturePlot(immune, features=treg_markers)

# CD8 T cells
cd8_tcell_markers <- c("CD8A")
DotPlot(immune, features=c("PTPRC", cd8_tcell_markers))
FeaturePlot(immune, features=cd8_tcell_markers)

# CD4 T cells, CD8A-
cd4_tcell_markers <- c("CD4", "IL7R", "CD8A")
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
saveRDS(immune, file = "/scratch/hnatri/ILD/Seurat_objects/immune_sct_ncells5k_nfeat3k_dims15_clusters_res2_woutDoublets_dims14_annotated_manualannot_210805.rds")

# =====================================
# Mesenchymal cell population
# ====================================

# Split out the mesenchymal population
mesen.integrated <- subset(ild_all,
                           cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Mesenchymal",]))

DefaultAssay(mesen.integrated) <- "integrated"

# Finding variable features, rescaling and running PCA
mesen.integrated <- FindVariableFeatures(mesen.integrated,
                                         assay = "integrated")
mesen.integrated.features <- VariableFeatures(mesen.integrated,
                                              assay = "integrated")
mesen.integrated <- ScaleData(mesen.integrated,
                              assay = "integrated",
                              verbose = T)
# TODO: redo
mesen.integrated <- RunPCA(mesen.integrated,
                           reduction.name = "integrated_pca",
                           assay = "integrated",
                           features = mesen.integrated.features,
                           verbose = T)

# How many PCs to include?
ElbowPlot(mesen.integrated, ndims=50)
best_pcs_mesen.integrated <- get_pcs(mesen.integrated, reduction_name="integrated_pca")

# Constructing the UMAP
mesen.integrated <- RunUMAP(mesen.integrated, dims = 1:25,
                            reduction = "integrated_pca",
                            reduction.name = "integrated_umap",
                            assay = "integrated",
                            verbose = T)
mesen.integrated <- FindNeighbors(mesen.integrated,
                                  dims = 1:25,
                                  reduction = "integrated_pca",
                                  assay = "integrated")
mesen.integrated.clusters <- FindClusters(object = mesen.integrated,
                                          resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                                          graph.name = "integrated_snn")

mesen <- mesen.integrated.clusters

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
mesen2 <- FindVariableFeatures(mesen2,
                               assay = "integrated")
mesen2.features <- VariableFeatures(mesen2)
mesen2 <- ScaleData(mesen2,
                    assay = "integrated")
mesen2 <- RunPCA(mesen2,
                 reduction.name = "integrated_pca",
                 assay = "integrated",
                 features = mesen2.features,
                 verbose = T)

# How many PCs to include?
ElbowPlot(mesen2, ndims=50)
best_pcs_mesen.integrated <- get_pcs(mesen2, reduction_name="integrated_pca")

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

best_pcs_mesen.ref <- get_pcs(mesen.ref)

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

mesen_col <- c("MyoFB" = "#00B13B",
               "HAS1 High Activated FB" = "#FF6B65",
               "SMC" = "#DC68F6",
               "MyoFB - Activated" = "#00B995",
               "Matrix FB" ="#D38721",
               "PLIN2+ FB" = "#4A92FA",
               "Pericyte" = "#00B0DC",
               "Mesothelial" = "#87A022",
               "WNT2+ FB" = "#FF57B9")

# cols = mesen_col
p1 <- DimPlot(mesen.ref, group.by = "celltypes", label = T, repel = T, cols = mesen_col) + 
    NoLegend() + ggtitle ("Reference")
p2 <- DimPlot(mesen, group.by = "predicted.id", label = T, repel = T, cols = mesen_col) + 
    NoLegend() + ggtitle ("Query")

p1 + p2

saveRDS(mesen, file = "/scratch/hnatri/ILD/Seurat_objects/mesenchymal_sct_ndims25_res2_woutDoublets_ndims16_annotated_210809.rds")

# Manual annotations
DefaultAssay(mesen) <- "SCT"

mesothelial_markers <- c("MSLN", "UPK3B", "HP", "WT1")
FeaturePlot(mesen, features = mesothelial_markers)

smc_markers <- c("ACTA2", "PDGFRB", "MYH11", "TAGLN", "DES", "ACTG2")
FeaturePlot(mesen, features = smc_markers)

pericyte_markers <- c("ACTA2", "PDGFRB", "RGS5", "HIGD1B", "GJA4")
FeaturePlot(mesen, features = pericyte_markers)

fibroblast_markers <- c("LUM", "DCN", "PDGFRA")
FeaturePlot(mesen, features = fibroblast_markers)

myofb_markers <- c("MYLK", "ACTA2", "COL8A1", "COL1A1", "WNT2")
FeaturePlot(mesen, features = myofb_markers)

activated_myofb_markers <- c("COL1A1", "POSTN", "CTHRC1")
FeaturePlot(mesen, features = activated_myofb_markers)

WNT2_fibro_markers <- c("MYLK", "WNT2", "A2M", "GPC3", "MACF1", "CES1", "LIMCH1")
FeaturePlot(mesen, features = WNT2_fibro_markers)

matrixfb_markers <- c("SFRP2", "CLU", "APOD", "FBLN1", "CST3", "IGFBP6", "SCARA5", "CD34")
FeaturePlot(mesen, features = matrixfb_markers)

PLIN2_fibro_markers <- c("PLIN2", "HAS1")
FeaturePlot(mesen, features = PLIN2_fibro_markers)

HAS1_fibro_markers <- c("HAS1", "TWIST1", "PLIN2")
FeaturePlot(mesen, features = HAS1_fibro_markers)

act_HAS1_fibro_markers <- c("LIF", "ICAM1", "CXCL2", "HAS1", "PLIN2")
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

saveRDS(mesen, "/scratch/hnatri/ILD/Seurat_objects/mesenchymal_integrated_ndims17_res2_woutDoublets_ndims16_annotated_210809_manualannot.rds")

# =====================================
# Endothelial cell population
# =====================================

# Split out the endothelial population
endo.integrated <- subset(ild_all,
                          cells = rownames(ild_all@meta.data[ild_all@meta.data$population == "Endothelial",]))

DefaultAssay(endo.integrated) <- "integrated"

# Finding variable features, rescaling and running PCA
endo.integrated <- FindVariableFeatures(endo.integrated,
                                        assay = "integrated")
endo.integrated.features <- VariableFeatures(endo.integrated,
                                             assay = "integrated")
endo.integrated <- ScaleData(endo.integrated,
                             assay = "integrated",
                             verbose = T)
endo.integrated <- RunPCA(endo.integrated,
                          reduction.name = "integrated_pca",
                          assay = "integrated",
                          features = endo.integrated.features,
                          verbose = T)

# How many PCs to include?
ElbowPlot(endo.integrated, ndims=50)
best_pcs_endo.integrated <- get_pcs(endo.integrated, reduction_name="integrated_pca")

# Constructing the UMAP
endo.integrated <- RunUMAP(endo.integrated, dims = 1:17,
                           reduction = "integrated_pca",
                           reduction.name = "integrated_umap",
                           assay = "integrated",
                           verbose = T)
endo.integrated <- FindNeighbors(endo.integrated,
                                 dims = 1:17,
                                 reduction = "integrated_pca",
                                 assay = "integrated")
endo.integrated.clusters <- FindClusters(object = endo.integrated,
                                         resolution = c(0.1, 0.5, 1, 1.5, 2, 2.5),
                                         graph.name = "integrated_snn")

endo <- endo.integrated.clusters

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
endo2 <- FindVariableFeatures(endo2,
                              assay = "integrated")
endo2.features <- VariableFeatures(endo2)
endo2 <- ScaleData(endo2,
                   assay = "integrated")
endo2 <- RunPCA(endo2,
                reduction.name = "integrated_pca",
                assay = "integrated",
                features = mesen2.features,
                verbose = T)

# How many PCs to include?
ElbowPlot(endo2, ndims=50)
best_pcs_endo2 <- get_pcs(endo2, reduction_name="integrated_pca")

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

# Removing mesenchymal populations
endo.cells <- c("Endothelial - capillary","Endothelial - venule","Endothelial - arteriole", 
                "Lymphatic", "Endothelial - inflamed","Endothelial - CA4+ capillary", 
                "Endothelial - peribronchiolar")
endo.ref <- subset(endo.ref, cells=rownames(endo.ref@meta.data[endo.ref@meta.data$celltype %in% endo.cells,]))
endo.ref <- SCTransform(endo.ref)
endo.ref <- RunPCA(endo.ref)

best_pcs_endo.ref <- get_pcs(endo.ref)

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

endo_col <- c("Endothelial - capillary" = "#FF6B65",
              "Endothelial - venule" = "#00B7BB",
              "Endothelial - arteriole" = "#00ADE5",
              "Lymphatic" = "#009AFA", 
              "Endothelial - inflamed" = "#9B7FFA",
              "Endothelial - CA4+ capillary" = "#FF57CF", 
              "Endothelial - peribronchiolar" = "#FF5C9E")

p1 <- DimPlot(endo.ref, group.by = "celltype", label = T, repel = T, cols = endo_col) + 
    NoLegend() + ggtitle("Reference")
p2 <- DimPlot(endo, group.by = "predicted.id", label = T, repel = T, cols = endo_col) + 
    NoLegend() + ggtitle("Query")

p1 + p2

saveRDS(endo, file = "/scratch/hnatri/ILD/Seurat_objects/endothelial_sct_ndims19_res1_woutDoublets_ndims17_annotated_210808.rds")

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

saveRDS(endo, file = "/scratch/hnatri/ILD/Seurat_objects/endothelial_sct_ndims19_res1_woutDoublets_ndims17_annotated_210808_manualannot.rds")

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
epi <- subset(epi, features = rownames(counts))

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

DimPlot(ild, group.by = "manual_annotation_1", label = T, repel = T) + NoLegend()

# Saving the object
saveRDS(ild, file = "/scratch/hnatri/ILD/Seurat_objects/ILD_annotated_210819.rds")

# Less granular annotations
simple_celltypes <- ild@meta.data$predicted.id
simple_celltypes <- multigsub(c("Secretory - SCGB1A1+/MUC5B+",
                                "Secretory - MUC5AC+",
                                "Secretory - SCGB1A1+/SCGB3A2+",
                                "Secretory - SCGB3A1+/MUC5B+",
                                "Secretory - SCGB3A2+"),
                              replicate(5, "Secretory"), simple_celltypes)

simple_celltypes <- multigsub(c("Endothelial - capillary",
                                "Endothelial - CA4+ capillary",
                                "Endothelial - inflamed",
                                "Endothelial - peribronchiolar",
                                "Endothelial - venule",
                                "Endothelial - arteriole"),
                              replicate(6, "Vascular endothelial"), simple_celltypes)

simple_celltypes <- gsub("AT2 - low quality", "AT2", simple_celltypes)
simple_celltypes <- gsub("Lymphatic", "Lymphatic endothelial", simple_celltypes)

simple_celltypes <- gsub("Macrophage - SPP1+", "Macrophage", simple_celltypes)

simple_celltypes <- multigsub(unique(mesen@meta.data$predicted.id),
                              replicate(length(unique(mesen@meta.data$predicted.id)), "Fibroblasts"), simple_celltypes)

ild@meta.data$simple_celltypes <- simple_celltypes


filename <- "/scratch/hnatri/ILD/Processing_annotation/all_pops_umap_simple_celltypes.pdf"
pdf(file = filename,
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
DimPlot(ild, group.by = "simple_celltypes", label=T) +
    theme(legend.position = "none")
dev.off()

ild_lowdims <- RunUMAP(ild, dims=1:5)
ild_highdims <- RunUMAP(ild, dims=1:50)
ild_highdims2 <- ild_highdims

ild_highdims2@meta.data$simple_celltypes <- simple_celltypes

filename <- "/scratch/hnatri/ILD/Processing_annotation/all_pops_umap_simple_celltypes_dims50.pdf"
pdf(file = filename,
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches
DimPlot(ild_highdims2, group.by = "simple_celltypes", label=T, repel=T) +
    theme(legend.position = "none") + ggtitle("")
dev.off()

saveRDS(ild, file = "/scratch/hnatri/ILD/Seurat_objects/ILD_annotated_210819.rds")
