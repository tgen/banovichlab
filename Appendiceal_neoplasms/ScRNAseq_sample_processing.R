# ==============================================================================
# Date: 2023/02/16
# Author: Linh T. Bui (lbui@tgen.org)
# Sample processing with SoupX and DoubletFinder on Appendix cancer scRNAseq
# Cell type annotation all levels
# This file also contains codes used to generate the following figures:
# Supplementary Figure 1
# Supplementary Figure 2
# Supplementary Figure 4
# Supplementary Figure 6
# Supplementary Figure 9
# Supplementary Figure 11
# ==============================================================================
# ==============================================================================
# SET UP THE ENVIRONMENT VARIABLES 
# ==============================================================================
# Set up working directory
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

options(future.globals.maxSize = 4096*1024^2 )
set.seed(12345)

# Load required libraries
library(Matrix)
library(SoupX)
library(Seurat)
library(DoubletFinder)
library(findPC)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(nord)
library(circlize)
library(ComplexHeatmap)

# ==============================================================================
# Sample processing with SoupX and DoubletFinder for all samples 
# These steps were documented here for your preference
# ==============================================================================
#tmpDir = setwd("/scratch/lbui/Appendiceal_data")
#dataset_loc <- "/labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/COHP/CellRangerOuts/GeneExpression/"
#inputfiles_toc <- "outs/filtered_feature_bc_matrix"
#inputfiles_tod <- "outs/raw_feature_bc_matrix"

#batchids <-  c("COHP_41962_1_UN_Whole_C1_X3SC3_L35149_HLWVKDSXY",
#               "COHP_42049_1_UN_Whole_C1_X3SC3_L35152_HLWVKDSXY",
#               "COHP_47035_2_UN_Whole_C1_X3SC4_K28667_HG33KDSX3",
#               "COHP_38494_1_UN_Whole_C1_X3SC3_L35144_HLWVKDSXY",
#               "COHP_39915_1_UN_Whole_C1_X3SC3_L35145_HLWVKDSXY",
#               "COHP_40832_1_UN_Whole_C1_X3SC3_L35146_HLWVKDSXY",
#               "COHP_41423_1_UN_Whole_C1_X3SC3_L35147_HLWVKDSXY",
#               "COHP_41622_1_UN_Whole_C1_X3SC3_L35148_HLWVKDSXY",
#               "COHP_41963_1_UN_Whole_C1_X3SC3_L35150_HLWVKDSXY",
#               "COHP_41993_1_UN_Whole_C1_X3SC3_L35151_HLWVKDSXY",
#               "COHP_45655_2_UN_Whole_C1_X3SC4_K28668_HG33KDSX3",
#               "COHP_46684_1_UN_Whole_C1_X3SC4_L42572_HCY7WDSX3",
#               "COHP_48936_2_UN_Whole_C1_X3SC4_L49293_HCLW7DSX5",
#               "COHP_49354_1_UN_Whole_C1_X3SC4_L49294_HCLW7DSX5",
#               "COHP_49378_1_UN_Whole_C1_X3SC4_L49295_HCLW7DSX5",
#               "COHP_49422_1_UN_Whole_C1_X3SC4_L49296_HCLW7DSX5")

#d10x.data <- sapply(batchids,  function(i){
#  sample_id <- sapply(strsplit(i,split="_") , "[[", 8)
#  d10x_toc <- Read10X(file.path(dataset_loc,i,inputfiles_toc))
#  d10x_tod <- Read10X(file.path(dataset_loc,i,inputfiles_tod))
#  colnames(d10x_toc) <- paste(sample_id, sep="_", colnames(d10x_toc))
#  colnames(d10x_tod) <- paste(sample_id, sep="_", colnames(d10x_tod))
#  d10x_tod
#  d10x_toc
  # Run SoupX
#  sc <- SoupChannel(d10x_tod, d10x_toc, calcSoupProfile = FALSE)
#  sc <- estimateSoup(sc)
#  toc_seu <- CreateSeuratObject(d10x_toc)
#  toc_seu <- SCTransform(toc_seu, vst.flavor="v2")
#  toc_seu <- RunPCA(toc_seu)
#  toc_seu <- RunUMAP(toc_seu, dims = 1:20)
#  toc_seu <- FindNeighbors(toc_seu, dims = 1:20)
#  toc_seu <- FindClusters(toc_seu, resolution = 0.5)
    ## Add meta data to soupX object
#  sc <- setClusters(sc, setNames(toc_seu$seurat_clusters, rownames(toc_seu@meta.data)))
    ## Estimate contamination (automated method) & adjust counts
#  sc <- autoEstCont(sc)
#  out <- adjustCounts(sc)
  
  # Create Seurat object using corrected count matrix
#  d10x_seu <- CreateSeuratObject(out, min.features = 200, names.field = 1, names.delim = "_")
#  d10x_seu <- PercentageFeatureSet(object = d10x_seu, pattern = "^MT-", col.name = "percent.mt")
#  d10x_seu <- subset(d10x_seu, subset = nCount_RNA > 500 & percent.mt < 25)
  # In my experience, FindDoublet doesn't work very well for sctransformed data
  # Using SCTransform, nFeature_SCT in doublets are not doubled than singlets, more like equal number
#  d10x_seu <- NormalizeData(d10x_seu) %>% 
#  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
#    ScaleData () %>%
#    RunPCA(verbose = FALSE) %>%
#    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
#    FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
#    FindClusters(resolution = 1, verbose = FALSE)
  # Run DoubletFinder to identify doublets
#  sweep.list_test <- paramSweep_v3(d10x_seu, PCs = 1:20, sct = F)
#  sweep.stats_test <- summarizeSweep(sweep.list_test, GT = F)
#  bcmvn <- find.pK(sweep.stats_test)
  # Select nExp_poi
#  homotypic.prop <- modelHomotypic(d10x_seu@meta.data$seurat_clusters)           
#  nExp_poi <- round(0.05*nrow(d10x_seu@meta.data))  ## Assuming 5% doublet formation rate - from 10X manual
#  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#  pk <- as.numeric(as.vector(bcmvn$pK)[which.max(bcmvn$BCmetric)])
  ## Run DoubletFinder 
#  d10x_seu <- doubletFinder_v3(d10x_seu, PCs = 1:20, pN = 0.25, pK = pk, 
#                               nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
#  d10x_seu <- doubletFinder_v3(d10x_seu, PCs = 1:20, pN = 0.25, pK = pk, 
#                               nExp = nExp_poi.adj, 
#                               reuse.pANN = colnames(d10x_seu@meta.data[7]), sct = F)
#  names(d10x_seu@meta.data)[9] <- "doublet_finder" ## rename the column so all objects have same column name
#  names(d10x_seu@meta.data)[7] <- "pANN_col"
#  names(d10x_seu@meta.data)[8] <- "DF.classifications_col"
#  d10x_seu
#})

# Merge all objects
# coh.ob <- Reduce(merge, d10x.data)
# saveRDS(coh.ob, file = "Appendiceal_allsamples_soupX_doubletfinder_noSCT.rds")
# rm(d10x.data)

# Write out the matrix file for submission
# writeMM(coh.ob@assays$RNA@counts, file = "matrix.mtx")
# write(x = rownames(coh.ob@assays$RNA@counts), file = "genes.tsv")
# write(x = colnames(coh.ob@assays$RNA@counts), file = "barcodes.tsv")

# Add in meta data 
# meta.data <- read.csv(file = "/labs/banovich/Notebook/Linh/Appendiceal_Seurat/Appendiceal_meta.data.csv",
#                      header = T, stringsAsFactors = F)
# coh.ob@meta.data$Diagnosis <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                              from = meta.data$Library_ID,
#                                              to = as.character(meta.data$Clinical_Diagnosis))
# coh.ob@meta.data$Status <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                           from = meta.data$Library_ID,
#                                           to = as.character(meta.data$Status))
# coh.ob@meta.data$Gender <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                           from = meta.data$Library_ID,
#                                           to = as.character(meta.data$Gender))
# coh.ob@meta.data$Age <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                        from = meta.data$Library_ID,
#                                        to = as.character(meta.data$Age))
# coh.ob@meta.data$Ethnicity <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                              from = meta.data$Library_ID,
#                                              to = as.character(meta.data$Ethnicity))
# coh.ob@meta.data$Pathology <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                              from = meta.data$Library_ID,
#                                              to = as.character(meta.data$Pathology_Reviewed))
# coh.ob@meta.data$Flowcell <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                             from = meta.data$Library_ID,
#                                             to = as.character(meta.data$Flowcell_ID))
# coh.ob@meta.data$Pathology2 <- plyr::mapvalues(x = coh.ob@meta.data$orig.ident,
#                                               from = meta.data$Library_ID,
#                                               to = as.character(meta.data$Pathology_abb))

# Save the meta data file
# write.csv(coh.ob@meta.data, file = "Appendiceal_meta_data_allcells_noannotation.csv")

#===============================================================================
# Seurat integration and cell annotation start from here
#===============================================================================
# Read in the dgCMatrix count matrix file
sc_data <- Read10X(data.dir = ".",
                    gene.column = 1)

# Read in the meta data file
meta_data <- read.csv("Appendiceal_meta_data_allcells_noannotation.csv",
                      header = T)
rownames(meta_data) <- meta_data$X
meta_data$X <- NULL
head(meta_data)

# Create Seurat object
coh.ob <- CreateSeuratObject(sc_data,
                          project = "scRNA appendiceal",
                          meta.data = meta_data)

# Check doublet cells from DoubletFinder
VlnPlot(coh.ob, "nFeature_RNA", group.by = "orig.ident",split.by = "doublet_finder",
        split.plot = T, pt.size = 0)
table(coh.ob$doublet_finder)
coh.ob <- subset(coh.ob, subset = doublet_finder == "Singlet")

# Run integration pipeline to remove batch effects caused by different sequencing batches
seu_list <- SplitObject(coh.ob, split.by = "Flowcell")
seu_list <- lapply(seu_list, function(x) SCTransform(x, vst.flavor="v2",
                                                     return.only.var.genes=FALSE,
                                                     method = "glmGamPoi",
                                                     vars.to.regress = "percent.mt"))
features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = features)
seu_list <- lapply(seu_list, FUN = RunPCA, features = features)
coh.anchors <- FindIntegrationAnchors(object.list = seu_list, 
                                      normalization.method = "SCT",
                                      anchor.features = features, 
                                      reduction = "rpca",
                                      k.anchor = 10) # I tried different k.anchor value and 10 seems to work best

# Add cell type marker genes into the features to integrate for downstream granular annotation
more_features <- c(features, "EPCAM","CD68","CD79A","CD3E","PECAM1","DCN","THY1",
                   "CXCR2","FCGR3B", "CD4","CCR7","IL7R","CD8A","CD8B","IL2RA",
                   "FOXP3", "NKG7", "GNLY","FCGR3A","MKI67","PCNA","IL1B","NLRP3",
                   "CD14","CLEC12A","FCN1","SPP1","MKI67","CD1C","LILRA4","BATF3",
                   "INHBA","LYVE1","C1QC","LAMP3","FCGR3A","KIT","LYZ","MS4A7",
                   "COL1A2","FAP","VIM", "FOXF1", "TAGLN", "LUM","RGS5", "SOX6",
                   "IL18","WT1","MSLN","ACTA2", "FAP","ACTA2","MYL9","HOPX","IL6",
                   "CXCL12","CXCL14","DPT","CD74","HLA-DRA","HLA-DPA1","HLA-DRB1",
                   "IGHG1","IGHA1","IGHD","MS4A1","LRMP","JCHAIN","LGR5","ASCL2",
                   "OLFM4","CEACAM6","MUC2","MUC5B","TFF3","SPINK4","CLCA1","SPDEF",
                   "FCGBP","CA2","SLC26A2","FABP1", "FABP2","CEACAM1","CHGA","TPH1",
                   "NEUROD1","UBE2C","TOP2A","SPDEF","WFDC2","ASCL2","SMOC2","BEST2",
                   "OTOP2","POU2F3","LRMP","TRPM5","HAND2","PHOX2B","HAND2","TUBB2B",
                   "S100B", "PHOX2B", "PLP1", "ACKR1", "VWF","GJA4", "HEY1", "PROX1", 
                   "PDPN")

var_list <- lapply(seu_list, function(x) FindVariableFeatures(x, nfeatures = 10000))
all_genes <- lapply(var_list, VariableFeatures) %>% Reduce(intersect, .)

mt_genes <- grep("^MT", all_genes, value = TRUE)
ribo_genes <- grep("^RB", all_genes, value = TRUE)
all_genes <- all_genes[!all_genes %in% c(mt_genes, ribo_genes)]
coh.combined.sct <- IntegrateData(anchorset = coh.anchors, normalization.method = "SCT",
                                  features.to.integrate = unique(c(all_genes, more_features)))

coh.combined.sct <- RunPCA(coh.combined.sct, verbose = FALSE)

# ----------------------------------
# Supplementary Figure 1
# ----------------------------------
# Select number of PCs 
sdev <- prcomp(t(GetAssayData(coh.combined.sct, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:30]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(7, 30, by = 4), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(7, 30, by = 4), 
                       method = 'perpendicular line', aggregate = "voting")

# Run dimensional reduction and check for batch effects on Flowcell
coh.combined.sct <- RunUMAP(coh.combined.sct, reduction = "pca", dims = 1:15)

DimPlot(coh.combined.sct, group.by = "Flowcell")
rm(coh.anchors, seu_list, coh.ob)

# ==============================================================================
# Level 1 annotation
# ==============================================================================
# Run FindNeighbors and FindClusters using the PCs selected above
coh.combined.sct <- FindNeighbors(coh.combined.sct, dims = 1:15)
coh.combined.sct <- FindClusters(coh.combined.sct, resolution = 0.1)

# Level 1 annotation
DimPlot(coh.combined.sct, label = T) + NoLegend() -> p1
FeaturePlot(coh.combined.sct, c("EPCAM","CD68","CD79A","CD3D","PECAM1","DCN","KIT",
                                "CXCR2","FCGR3B"), min.cutoff = "q9", slot = "scale.data",
            ncol = 3) -> p2

# ----------------------------------
# Supplementary Figure 1
# ----------------------------------
p1 + p2 

DotPlot(coh.combined.sct, features = c("EPCAM","CD68","KIT","TPSB2","CD79A","CD3E","PECAM1",
                                       "THY1","DCN","COL3A1","CXCR2",'FCGR3B'))
# FCGR3B (CD16B) is a surface protein on NK, neutrophil, monocyte and macrophage

onion <- as.character(coh.combined.sct@meta.data$seurat_clusters)
onion[onion == 8] <- "Epithelial"
onion[onion %in% c(4,10)] <- "Myeloid cells"
onion[onion == 9] <- "FCGR3B+ neutrophils"
onion[onion %in% c(2,6,7)] <- "B cells"
onion[onion %in% c(0,3)] <- "T cells"
onion[onion == 1] <- "Mesenchymal cells"
onion[onion == 5] <- "Endothelial"
coh.combined.sct$Celltype1 <- onion

# UMAP plot 
DimPlot(coh.combined.sct, group.by = "Celltype1", label = T, repel = T) + NoLegend()

# DotPlot for marker gene expression 
DotPlot(coh.combined.sct,
        features = c("EPCAM","FABP1","MUC2","CD68","CD14","FCGR3A","CD79A","IGHD",
                     "MS4A1","CD3E","CD3D","IL7R","DCN","THY1","LUM","PECAM1",
                     "ACKR1","VWF","CXCR2","FCGR3B","KIT"),
        group.by = "Celltype1", cols = "RdBu") +
  coord_flip() +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=15))

# Save this object
saveRDS(coh.combined.sct, 
        file = "/scratch/lbui/Appendiceal_data/Appendiceal_integrated_rpca_CT1.rds")

# ==============================================================================
# Further annotate the T cell population
# ==============================================================================
# Subset out T cells
t_cells <- subset(coh.combined.sct,
                  subset = Celltype1 == "T cells")

# Find Variable Features, rescale the integrated assay and rerun PCA to get new PCs
DefaultAssay(t_cells) <- "RNA"
t_cells <- FindVariableFeatures(t_cells, nfeatures = 3000, selection.method = "vst")
DefaultAssay(t_cells) <- "integrated"
t_cells <- ScaleData(t_cells, vars.to.regress = "percent.mt")

# Select number of PCs for dimensional reduction
t_cells <- RunPCA(t_cells)
sdev <- prcomp(t(GetAssayData(t_cells, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:20]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(6, 20, by = 2), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(6, 20, by = 2), 
                       method = 'perpendicular line', aggregate = "voting")

# Unsupervised clustering and CT2 annotation
t_cells <- RunUMAP(t_cells, dims = 1:15)
t_cells <- FindNeighbors(t_cells, dims = 1:15)
t_cells <- FindClusters(t_cells, resolution = 0.5)

# ----------------------------------
# Supplementary Figure 6a
# ----------------------------------
DimPlot(t_cells, label = T, repel = T, raster = TRUE) + NoLegend() -> p1
FeaturePlot(t_cells, c("CD4","IL7R","CCR7","CD8A","CD8B","FOXP3",
                       "NKG7","GNLY","MKI67"),
            min.cutoff = "q7") -> p2
p1+p2

DotPlot(t_cells, features = c("CD4","CCR7","IL7R","CD8A","CD8B","IL2RA","FOXP3",
                              "NKG7", "GNLY","FCGR3A","MKI67","PCNA","CD3E"))

onion <- as.character(t_cells$seurat_clusters)
onion[onion %in% c(0,12)] <- "CD8+ T cells"
onion[onion %in% c(1,2,3,4,6,8,9,10,11,13,15,16,17)] <- "CD4+ T cells"
onion[onion == 14] <- "Proliferating T cells"
onion[onion == 7] <- "Tregs"
onion[onion == 5] <- "NK cells"
t_cells$Celltype2 <- onion

# Check annotation with violin plot
colors_list <- c("#FF7671", "#B2B33A", "#00BF80", "#00B0F3", "#EE6CF0")
Stacked_VlnPlot(t_cells, 
                features = c("CD4","IL7R","CCR7","CD8A","CD8B","FOXP3","NKG7","GNLY",
                             "MKI67","PCNA"), 
                group.by = "Celltype2", assay = "SCT", x_lab_rotate = TRUE,
                colors_use = colors_list) 

# UMAP figures
DimPlot(t_cells, group.by = "Celltype2")
DimPlot(t_cells, group.by = "Flowcell")
DimPlot(t_cells, group.by = "orig.ident")

# Save object
saveRDS(t_cells, file = "Appendiceal_Tcells_CT2.rds")

# ==============================================================================
# Further annotate the Myeloid cell population
# ==============================================================================
# Subset out myeloid cells
myeloid <- subset(coh.combined.sct,
                  subset = Celltype1 == "Myeloid cells")

# Find variable features and rescale the integrated assay
DefaultAssay(myeloid) <- "RNA"
myeloid <- FindVariableFeatures(myeloid, nfeatures = 3000, selection.method = "vst")
DefaultAssay(myeloid) <- "integrated"
myeloid <- ScaleData(myeloid, vars.to.regress = "percent.mt")
myeloid <- RunPCA(myeloid) 

# Select number of PCs for dimensional reduction
sdev <- prcomp(t(GetAssayData(myeloid, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:30]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(6, 20, by = 2), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(6, 20, by = 2), 
                       method = 'perpendicular line', aggregate = "voting")

# Run UMAP and find clusters
myeloid <- RunUMAP(myeloid, dims = 1:14)
myeloid <- FindNeighbors(myeloid, dims = 1:14)
myeloid <- FindClusters(myeloid, resolution = 1.5)

# Check to see if there's any batch effect 
DimPlot(myeloid, group.by = "Flowcell") #looks good
DimPlot(myeloid, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()

# ----------------------------------
# Supplementary Figure 9
# ----------------------------------
# Annotation Level 2
cluster_col <- nord("aurora",24)
names(cluster_col) <- as.character(unique(myeloid$seurat_clusters))
DefaultAssay(myeloid) <- "SCT"
plot_features <- c("CD14","FCGR3A","FCN1","LILRA4","CD1C","BATF3","NLRP3","PLTP",
                   "LYZ","IL1B","C1QC","C1QA","SPP1","LUM")
source("/home/lbui/SC_scripts/MMRF_CRISPR/CART_plot_functions.R") #code from Heini
create_dotplot_heatmap_horizontal(myeloid, 
                                  plot_features = plot_features,
                                  group_var = "seurat_clusters",
                                  group_colors = cluster_col,
                                  column_title = "Myeloid markers")

# Annotate cells
onion <- as.character(myeloid$seurat_clusters)
onion[onion %in% c(0,1,7,19,21)] <- "Macrophages"
onion[onion %in% c(5,6,11)] <- "Monocyte-like"
onion[onion %in% c(10,14,23)] <- "SPP1+ macrophages"
onion[onion %in% c(3,4,8,9,13)] <- "C1Qhi monocytes"
onion[onion %in% c(18)] <- "cDC1"
onion[onion %in% c(2)] <- "cDC2"
onion[onion %in% c(12,15,20)] <- "pDCs"
onion[onion %in% c(16,17,22)] <- "Mesenchymal cells"
myeloid$Celltype2 <- onion

# Remove mesenchymal cells
myeloid2 <- subset(myeloid, subset = Celltype2 != "Mesenchymal cells")

# Save object
saveRDS(myeloid2, 
        file = "Appendiceal_Myeloidcells_CT2.rds")

# Adjust CT1 label on the large object
coh.combined.sct$Celltype1 <- ifelse(rownames(coh.combined.sct@meta.data) %in%
                                       rownames(myeloid@meta.data[myeloid@meta.data$Celltype2 == "Mesenchymal cells",]),
                                     "Mesenchymal cells", coh.combined.sct$Celltype1)

# ==============================================================================
# Further annotate the Mesenchymal cell population
# ==============================================================================
# Subset out mesenchymal cells
mesen <- subset(coh.combined.sct,
                subset = Celltype1 == "Mesenchymal cells")

# Find variable features and rescale the integrated assay
DefaultAssay(mesen) <- "RNA"
mesen <- FindVariableFeatures(mesen, nfeatures = 3000, selection.method = "vst")
DefaultAssay(mesen) <- "integrated"
mesen <- ScaleData(mesen, vars.to.regress = "percent.mt")
mesen <- RunPCA(mesen)

# Select number of PCs for dimensional reduction
sdev <- prcomp(t(GetAssayData(mesen, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:30]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(5, 30, by = 5), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(5, 30, by = 5), 
                       method = 'perpendicular line', aggregate = "voting")

# Run UMAP and find clusters
mesen <- RunUMAP(mesen, dims = 1:10)
mesen <- FindNeighbors(mesen, dims = 1:10)
mesen <- FindClusters(mesen, resolution = 1)

# Check to see if there's still batch effect 
DimPlot(mesen, group.by = "Flowcell") 

# ----------------------------------
# Supplementary Figure 11a-b
# ----------------------------------
DimPlot(mesen, group.by = "seurat_clusters")
cluster_col <- nord("aurora",21)
names(cluster_col) <- as.character(unique(mesen$seurat_clusters))
DefaultAssay(mesen) <- "SCT"
plot_features <- c("THY1", "COL1A2", "VIM", "TAGLN","FAP","RGS5","MSLN",
                   "WT1","ACTA2")
source("/home/lbui/SC_scripts/MMRF_CRISPR/CART_plot_functions.R") #code from Heini
create_dotplot_heatmap_horizontal(mesen, 
                                  plot_features = plot_features,
                                  group_var = "seurat_clusters",
                                  group_colors = cluster_col,
                                  column_title = "Mesen markers")

# Annotation Level 2
DimPlot(mesen, label = T, repel = T) + NoLegend() 
VlnPlot(mesen, c("THY1", "COL1A2", "VIM", "TAGLN","FAP","RGS5","MSLN",
                     "WT1","ACTA2"), pt.size = 0) 

onion <- as.character(mesen$seurat_clusters)
onion[onion %in% c(2,3,6,7,14,15,18)] <- "CAFs"
onion[onion %in% c(4,5,8,10,13,17)] <- "Fibroblasts"
onion[onion %in% c(0,1,11,12,16,20)] <- "Myofibroblasts"
onion[onion == 9] <- "Pericytes"
onion[onion == 19] <- "SMC"
mesen$Celltype2 <- onion

# Subset out CAFs to further split into different subtypes
caf <- subset(mesen, subset = Celltype2 == "CAFs")
## Make plot for FAP expression (Supplementary Figure 11c)
FeaturePlot(caf, "FAP", split.by = "Status") & theme(legend.position = c(1,0.2))
## Adjust labels since none of the CAF cell in control expresses FAP 
mesen$Celltype2 <- ifelse(mesen$Celltype2 == "CAFs" & mesen$Status == "Control",
                          "Fibroblasts", mesen$Celltype2)
## Re-subset out CAFs to further split into different subtypes
caf <- subset(mesen, subset = Celltype2 == "CAFs")
DefaultAssay(caf) <- "integrated"
caf <- RunPCA(caf)
caf <- RunUMAP(caf, dims = 1:6)
caf <- FindNeighbors(caf, dims = 1:6)
caf <- FindClusters(caf, resolution = 0.5)
DimPlot(caf, label = T, repel = T) + NoLegend()

## Annotate CAFs
Idents(caf) <- as.character(caf$seurat_clusters)
caf_markers <- FindAllMarkers(caf, assay = "RNA")
top_genes <- caf_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
gene_list <- c("FAP",unique(top_genes$gene),"UPK3B","MMP2","ACTA2","CDH1","MSLN","WT1", #mesothelial markers
               "IL6","CXCL12","DPT","HOPX","MYL9","CD74","HLA-DRB1","HLA-DRA") #caf subtype markers 
gene_list <- gene_list[gene_list != "MT-ATP8"] #remove MT gene

# ----------------------------------
# Supplementary Figure 11d-e
# ----------------------------------
cluster_col <- nord("aurora",10)
names(cluster_col) <- as.character(unique(caf$seurat_clusters))
DefaultAssay(caf) <- "SCT"
source("/home/lbui/SC_scripts/MMRF_CRISPR/CART_plot_functions.R") #code from Heini
# Supplementary figure (marker dotplots)
create_dotplot_heatmap_horizontal(caf, 
                                  plot_features = unique(gene_list),
                                  group_var = "seurat_clusters",
                                  group_colors = cluster_col,
                                  column_title = "CAF markers")
DimPlot(caf, group.by = "seurat_clusters")

## Assign clusters to cell types
onion <- as.character(caf$seurat_clusters)
onion[onion %in% c(1,4,5)] <- "myCAFs"
onion[onion %in% c(0,3)] <- "apCAFs"
onion[onion %in% c(2,7,8)] <- "iCAFs"
onion[onion %in% c(6)] <- "fiCAFs"
onion[onion == 9] <- "Mesothelial"
caf$Celltype3 <- onion

## Add CAF subtype into the mesen object
mesen$Celltype3 <- plyr::mapvalues(x = rownames(mesen@meta.data),
                                   from = rownames(caf@meta.data),
                                   to = as.character(caf$Celltype3))
mesen$Celltype3 <- ifelse(mesen$Celltype3 %in% as.character(caf$Celltype3),
                          mesen$Celltype3, mesen$Celltype2)

# Make some plots
ct3_col_mesen_all <- c("Myofibroblasts" = "#00BF80", "SMC" = "#EE6CF0", 
                       "Fibroblasts" = "#A4A525", "Pericytes" = "#00B0F3", 
                       "iCAFs" = "#92AA25", "fiCAFs" = "#00BA43","apCAFs" = "#FF7671",
                       "myCAFs"="#00B9E1", "Mesothelial" = "#ff65ae")
DimPlot(mesen, group.by = "Celltype3", cols = ct3_col_mesen_all) 
DimPlot(caf, group.by = "Celltype3")

# VlnPlot for CAF subtypes
Idents(caf) <- as.character(caf$Celltype3)
Stacked_VlnPlot(seurat_object = caf, 
                features = c("FAP","ACTA2","MYL9","HOPX","IL6","CXCL12","CXCL14",
                             "DPT","CD74", "HLA-DRA","HLA-DRB1","KRT8","UPK3B"),
                x_lab_rotate = TRUE, colors_use = ct3_col_mesen_all,
                assay="SCT")

# Readjust CT2 to move Mesothelial out of CAFs
mesen$Celltype2 <- ifelse(mesen$Celltype2 == "CAFs" & mesen$Celltype3 == "Mesothelial",
                          "Mesothelial", mesen$Celltype2)

# Save object
saveRDS(mesen, file = "Appendiceal_Mesenchymal_CT2.rds")

# ==============================================================================
# Further annotate the B cell population
# ==============================================================================
# Subset out B cells
bcell <- subset(coh.combined.sct,
                subset = Celltype1 == "B cells")

# Find variable features and rescale the integrated assay
DefaultAssay(bcell) <- "RNA"
bcell <- FindVariableFeatures(bcell, nfeatures = 3000,selection.method = "vst")
DefaultAssay(bcell) <- "integrated"
bcell <- ScaleData(bcell, vars.to.regress = "percent.mt")
bcell <- RunPCA(bcell)

# Select number of PCs for dimensional reduction
sdev <- prcomp(t(GetAssayData(bcell, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:20]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(6, 20, by = 2), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(6, 20, by = 2), 
                       method = 'perpendicular line', aggregate = "voting")

# Run UMAP and find clusters
bcell <- RunUMAP(bcell, dims = 1:15)
bcell <- FindNeighbors(bcell, dims = 1:15)
bcell <- FindClusters(bcell, resolution = 1)

# Check to see if there's still batch effect 
DimPlot(bcell, group.by = "Flowcell") 

# ----------------------------------
# Supplementary Figure 6c
# ----------------------------------
# Annotation Level 2
DimPlot(bcell, label = T, repel = T) + NoLegend() -> p1
FeaturePlot(bcell, 
            features = c("IGHG1","IGHA1","IGHD","MS4A1","LRMP","JCHAIN"), 
            min.cutoff = "q9", slot = "scale.data") -> p2
p1+p2

DotPlot(bcell, features = c("IGHG1","IGHA1","IGHD","MS4A1","LRMP","JCHAIN",
                                         "CD79A","KIT"))

onion <- as.character(bcell$seurat_clusters)
onion[onion %in% c(15,17,18,21,22)] <- "Plasma B"
onion[onion %in% c(8,11,20)] <- "GALTB"
onion[onion %in% c(0,1,2,3,4,5,6,7,9,10,12,14,19,23,24,25,27)] <- "Follicular B"
onion[onion %in% c(13,16,28)] <- "GCBcell"
onion[onion == 26] <- "Mast cells"
bcell$Celltype2 <- onion

DimPlot(bcell, group.by = "Celltype2")

VlnPlot(bcell, features = c("IGHG1","IGHA1","IGHD","MS4A1","LRMP","KIT"), 
        group.by = "Celltype2", pt.size = 0)

# Move mast cells out of the B cell lineage
bcell$Celltype1 <- ifelse(bcell$Celltype2 == "Mast cells", "Mast cells", "B cells")

# ----------------------------------
# Supplementary Figure 6d
# ----------------------------------
bcell2 <- subset(bcell, subset = Celltype2 != "Mast cells")
VlnPlot(bcell2, features = c("IGHG1","IGHA1","IGHD","MS4A1","LRMP","CD79A"), 
        group.by = "Celltype2", pt.size = 0, assay = "SCT")

# Save object
saveRDS(bcell2, 
        file = "Appendiceal_Bcells_CT2.rds")

# ==============================================================================
# Further annotate the epithelial cells 
# For CNV analysis, see the Epithelial_infercnv.R script
# ==============================================================================
epi <- subset(coh.combined.sct,
               subset = Celltype1 == "Epithelial")

# Find variable features and rescale the intgrated assay
DefaultAssay(epi) <- "RNA"
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 3000)
DefaultAssay(epi) <- "integrated"
epi <- ScaleData(epi, vars.to.regress = "percent.mt")
epi <- RunPCA(epi)

# Select number of PCs
sdev <- prcomp(t(GetAssayData(epi, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:30]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(5, 30, by = 5), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(5, 30, by = 5), 
                       method = 'perpendicular line', aggregate = "voting")

# Dimensional reduction and find clusters
epi <- RunUMAP(epi, dims = 1:15)
epi <- FindNeighbors(epi, dims = 1:15)
epi <- FindClusters(epi, resolution = 0.5)
DimPlot(epi)

# ----------------------------------
# Supplementary Figure 4a-b
# ----------------------------------
# Annotation
DimPlot(epi, label = T, repel = T) + NoLegend()-> p1
FeaturePlot(epi, c("CA2","SPINK4","MUC2","TFF3","MUC5B","CEACAM6"), 
            min.cutoff = "q7") -> p2
p1+p2

DefaultAssay(epi) <- "SCT" #check gene expression in SCT assay
DimPlot(epi, label = T, repel = T) + NoLegend()-> p1
FeaturePlot(epi, c("CA2","SPINK4","MUC2","TFF3","MUC5B","CEACAM6"), 
            min.cutoff = "q7") -> p2
p1+p2

# Assign cell types to clusters
onion <- as.character(epi$seurat_clusters)
onion[onion %in% c(0,4)] <- "MUC5Bhi cells"
onion[onion %in% c(1,2,3,5,6,8,10,11,14)] <- "Goblet-like cells"
onion[onion %in% c(12,13)] <- "SPINK4hi cells"
onion[onion %in% c(7,9)] <- "Enterocytes"
epi$Celltype2 <- onion

# Make dimplot
DimPlot(epi, group.by = "Celltype2")

# Save object
saveRDS(epi, file = "Appendiceal_Epithelial_CT2.rds") 

# ==============================================================================
# Further annotate the Endothelial cells
# ==============================================================================
# Subset out the Endothelial cells
endo <- subset(coh.combined.sct,
              subset = Celltype1 == "Endothelial")

# Find variable features and rescale the intgrated assay
DefaultAssay(endo) <- "RNA"
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = 3000)
DefaultAssay(endo) <- "integrated"
endo <- ScaleData(endo, vars.to.regress = "percent.mt")
endo <- RunPCA(endo)

# Select number of PCs
sdev <- prcomp(t(GetAssayData(endo, 
                              assay = "integrated", 
                              slot = "scale.data")))$sdev[1:30]
pc_res <- findPC(sdev = sdev, method = "all", 
                 number = seq(5, 30, by = 5), 
                 figure = T)
selected_pcs <- findPC(sdev = sdev, number = seq(5, 30, by = 5), 
                       method = 'perpendicular line', aggregate = "voting")

# Dimensional reduction and find clusters
endo <- RunUMAP(endo, dims = 1:20)
endo <- FindNeighbors(endo, dims = 1:20)
endo <- FindClusters(endo, resolution = 0.5)

# Check batch effects
DimPlot(endo, group.by = "Flowcell")

# Cell type 2 annotation
DimPlot(endo, label = T, repel = T) + NoLegend() -> p1
FeaturePlot(endo, c("ACKR1", "VWF","GJA4", "HEY1","PROX1", "PDPN"), 
            min.cutoff = "q8") -> p2
p1+p2

VlnPlot(endo, c("ACKR1", "VWF","GJA4", "HEY1","PROX1", "PDPN"),
        pt.size = 0)

onion <- as.character(endo$seurat_clusters)
onion[onion %in% c(0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,18)] <- "Venous Endothelial cell"
onion[onion %in% c(7,17,19)] <- "Arterial Endothelial cells"
endo@meta.data$Celltype2 <- onion
DimPlot(endo, group.by = "Celltype2", split.by = "Status")

VlnPlot(endo, features = c("ACKR1", "VWF","GJA4", "HEY1"),
        group.by = "Celltype2", pt.size = 0) 

# Save object
saveRDS(endo, 
        file = "Appendiceal_Endothelial_CT2.rds")

# ==============================================================================
# Add all CT2 annotation to Appendiceal all lineages object
# ==============================================================================
ct2_list <- list(epi@meta.data, t_cells@meta.data, myeloid2@meta.data, 
                 mesen@meta.data, bcell@meta.data, endo@meta.data)

onion <- lapply(ct2_list, function(xx){
  as.data.frame(cbind(rownames(xx), xx$Celltype2, xx$Celltype1))
})

celltype2_meta <- do.call(rbind, onion)
colnames(celltype2_meta) <- c("Barcode","Celltype2","Celltype1")

coh.combined.sct$Celltype2 <- plyr::mapvalues(x=rownames(coh.combined.sct@meta.data),
                                              from=celltype2_meta$Barcode,
                                              to=celltype2_meta$Celltype2)
coh.combined.sct$Celltype2 <- ifelse(coh.combined.sct$Celltype1 == "FCGR3B+ neutrophils",
                                    "FCGR3B+ neutrophils", coh.combined.sct$Celltype2)
coh.combined.sct$Celltype1 <- ifelse(coh.combined.sct$Celltype2 == "Mast cells", 
                                     "Mast cells", coh.combined.sct$Celltype1)
coh.combined.sct$Celltype3 <- plyr::mapvalues(x=rownames(coh.combined.sct@meta.data),
                                              from=rownames(mesen@meta.data),
                                              to=mesen$Celltype3)

coh.combined.sct$Celltype3 <- ifelse(coh.combined.sct$Celltype3 %in% mesen$Celltype3,
                                     coh.combined.sct$Celltype3, coh.combined.sct$Celltype2)

saveRDS(coh.combined.sct, 
        file = "/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")

# Save tables with cell numbers (Supplementary Table 3)
write.csv(table(coh.combined.sct$Celltype3, coh.combined.sct$Status),
          file = "Appendiceal_CT3_Status_cellnumber.csv")
write.csv(table(coh.combined.sct$Celltype3, coh.combined.sct$Pathology2),
          file = "Appendiceal_CT3_Pathology_cellnumber.csv")

# ----------------------------------
# Supplementary Figure 2
# ----------------------------------
ct3_cols <- read.csv("/scratch/lbui/RStudio_folder/20230820/CT3_color_codes.csv")
ct3_cols$color <- paste0("#",ct3_cols$color)
DimPlot(coh.combined.sct, group.by = "Celltype3", cols = ct3_cols$color)

features <- c("CD4","IL7R","CCR7","CD8A","CD8B","FOXP3","NKG7","GNLY","FCN1",
              "S100A9","S100A8","FCGR3A","LST1","LILRB2","LILRA4","GZMB","IL3RA",
              "CPA3","CLEC9A","FLT3","IDO1","CD1C","FCER1A","HLA-DQA1",
              "LAMP3","CCR7","FSCN1","INHBA","IL1RN","CCL4",
              "NLRP3","ERG","IL1B","LYVE1","PLTP","SEPP1","C1QC",
              "C1QA","APOE","HLA-DQA1","HLA-DQB1","CD163","CD68",
              "IL1B","CCL4","CXCL2","CXCR4","AREG","EGFR1","SPP1",
              "ISG15","FN1", "COL1A2","DCN","THY1","FAP","VIM", "TAGLN", 
              "LUM","RGS5", "SOX6","IL18","WT1", "MSLN","ACTA2",
              "MYL9","HOPX","IL6","CXCL12","CXCL14", "DPT","CD74",
              "HLA-DRA","HLA-DPA1","HLA-DRB1","CD3E","IGHG1","IGHA1","IGHD",
              "MS4A1","LRMP","JCHAIN","KIT","TPSAB1","LGR5","ASCL2","OLFM4",
              "CEACAM6","MUC2","MUC5B","TFF3","SPINK4","CLCA1","SPDEF","FCGBP",
              "CA2","SLC26A2","FABP1", "POU2F3","LRMP","TRPM5", "MKI67",
              "ACKR1", "VWF","GJA4", "HEY1","PROX1", "PDPN","FCGR3B","CXCR2")

DotPlot(coh.combined.sct, features = unique(features), 
        group.by = "Celltype3", cols = "RdBu", assay = "SCT") + 
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8))

