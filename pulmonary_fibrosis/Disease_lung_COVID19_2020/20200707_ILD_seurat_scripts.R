# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-07-07
# Title: scRNA-seq data analysis for Covid-19 manuscript (revised)
# ====================================

# =====================================
# Enviornment variables
# =====================================
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()
options(future.globals.maxSize = 4096*1024^2 )
set.seed(2811)

# =====================================
# Load libraries
# =====================================
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RCurl)
library(gplots)
library(BSDA)
library(reshape)

# =====================================
# Load all the objects
# =====================================
# Banovich - Kropski dataset
bano <- readRDS("/scratch/lbui/RStudio_folder/20200511/20200511_ILD_batchcorrected.rds")
# Kaminski dataset
load("/scratch/agutierrez/IPF/R/Seurat/Objects/Clean.Soup107.forGEO.Seurat.Robj")
# Pittbursg dataset
pitt <- readRDS("/scratch/agutierrez/IPF/R/Misc/GSE128033/20200402/GSE128033_ILD_sct.rds")
# Northwestern dataset
northwt <- readRDS("/scratch/agutierrez/IPF/R/Misc/GSE122960/20200401/GSE122960_ILD_sct.rds")

# ======================================
# Preprocessing 
# ======================================
# Save Kaminski meta.data as separate file for annonation after merging
kaminski_metadata <- as.data.frame(cbind(intron.exon.107$orig.ident, intron.exon.107$Library_Identity,
                           intron.exon.107$Disease_Identity, intron.exon.107$Subject_Identity))
rownames(kaminski_metadata) <- NULL
kaminski_metadata$V1 <- NULL
kaminski_metadata <- unique(kaminski_metadata)
colnames(kaminski_metadata) <- c("Library_Identity","Diagnosis","Subject_Identity")
onion <- as.character(kaminski_metadata$Diagnosis)
onion[onion == "Control"] <- "Control"
onion[onion %in% c("COPD","IPF")] <- "Disease"
kaminski_metadata$Status <- onion

write.csv(kaminski_metadata, file = "Kaminski_metadata.csv")

# Reorganize Kaminski data set and Bano - Kropski dataset to keep only the counts
kaminski <- intron.exon.107
kaminski@meta.data <- kaminski@meta.data[ ,c("Library_Identity","nCount_RNA",
                                                "nFeature_RNA","percent.mito")]

colnames(kaminski@meta.data) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt")
rm(intron.exon.107)

# ======================================
# Split data set into population and merge individual populations
# ======================================
# Add a dataset identity column to each dataset
bano@meta.data$dataset <- "VUMC/TGen"
kaminski@meta.data$dataset <- "Yale/BWH"
northwt@meta.data$dataset <- "Northwestern"
pitt@meta.data$dataset <- "Pittsburg"

# Splitting Kaminski object
kaminski <- SCTransform(kaminski)
kaminski <- RunPCA(kaminski)
kaminski <- RunUMAP(kaminski, dims = 1:20)
kaminski <- FindNeighbors(kaminski, dims = 1:20)
kaminski <- FindClusters(kaminski, resolution = 0.01)

DotPlot(kaminski, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(kaminski@meta.data$seurat_clusters)
onion[onion %in% c(0,1,3,6)] <- "Immune"
onion[onion == 2] <- "Epithelial"
onion[onion == 5] <- "Endothelial"
onion[onion == 4] <- "Mesenchymal"
kaminski@meta.data$population <- onion

DimPlot(kaminski, group.by = "population")
DotPlot(kaminski, features = c("EPCAM","PTPRC","PECAM1"), group.by = "population")
saveRDS(kaminski, file = "20200521_Kaminski_pop.rds")

# Austin has already ran SCTransform, RunPCA and RunUMAP for the Pitt and Northwestern dataset
# Splitting Pittsburg object
pitt <- FindNeighbors(pitt, dims = 1:20)
pitt <- FindClusters(pitt, resolution = 0.01)

DotPlot(pitt, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(pitt@meta.data$seurat_clusters)
onion[onion %in% c(0,5)] <- "Immune"
onion[onion == 1] <- "Epithelial"
onion[onion %in% c(2,4)] <- "Endothelial"
onion[onion == 3] <- "Mesenchymal"
pitt@meta.data$population <- onion

DimPlot(pitt, group.by = "population")
saveRDS(pitt, file = "20200521_Pittsburg_pop.rds")

# Splitting Northwestern object
northwt <- FindNeighbors(northwt, dims = 1:20)
northwt <- FindClusters(northwt, resolution = 0.01)

DotPlot(northwt, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(northwt@meta.data$seurat_clusters)
onion[onion == 0] <- "Immune"
onion[onion %in% c(1,3)] <- "Epithelial"
onion[onion == 2] <- "Endothelial"
onion[onion == 4] <- "Mesenchymal"
northwt@meta.data$population <- onion

DimPlot(northwt, group.by = "population")
saveRDS(northwt, file = "20200521_Northwestern_pop.rds")

# Merge all dataset
ild <- merge(x=bano, y=c(kaminski, pitt, northwt))
saveRDS(ild, file = "20200526_ILD_pop.rds") # moved this file into scratch/lbui/Covid19_Seurat

# ==============================================================================
# Revised on 20200707
# ==============================================================================
# Load object
ild <- readRDS("/scratch/lbui/Covid19_Seurat/20200526_ILD_pop.rds")
# Add meta data
meta.data <- read.csv("/scratch/lbui/Covid19_Seurat/Metadata_Covid19_edited.csv", header = T)
ild@meta.data$Diagnosis <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                            from = meta.data$Library_ID,
                                            to = as.character(meta.data$Diagnosis))
ild@meta.data$Status <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                         from = meta.data$Library_ID,
                                         to = as.character(meta.data$Status))
ild@meta.data$Source_detail <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                        from = meta.data$Library_ID,
                                        to = as.character(meta.data$Source_detail))

ild@meta.data$Sample_Name <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                               from = meta.data$Library_ID,
                                               to = as.character(meta.data$Sample_Name))

ild@meta.data$Age <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                     from = meta.data$Library_ID,
                                     to = as.character(meta.data$Age))

ild@meta.data$Ethnicity <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                           from = meta.data$Library_ID,
                                           to = as.character(meta.data$Race))

ild@meta.data$Smoking_status <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                                from = meta.data$Library_ID,
                                                to = as.character(meta.data$Smoking_Status))
onion <- ild@meta.data$Diagnosis 

onion[onion == "Control"] <- "Control"
onion[onion == "IPF"] <- "IPF"
onion[onion == "COPD"] <- "COPD"
onion[onion %in% c("cHP","IPAF","CTD-ILD","ILD","NSIP","Sarcoidosis",
                   "Hypersensitivity pneumonitis")] <- "Other-ILD"

ild@meta.data$Diagnosis2 <- onion

# Filter out cells with low number of nfeature_RNA
ild <- PercentageFeatureSet(object = ild, pattern = "^MT-", col.name = "percent.mt")
smoothScatter(ild@meta.data$percent.mt, ild@meta.data$nFeature_RNA)
ild <- subset(ild, subset = nFeature_RNA > 1000 & percent.mt < 25)

# Save the object
saveRDS(ild, file = "/scratch/lbui/Covid19_Seurat/20200707_ILD_pop.rds")

# ==============================
# EPITHELIAL POPULATION
# ==============================
# Subset out the epithelial cells from ILD object
epi <- subset(ild, cells=rownames(ild@meta.data[ild@meta.data$population == "Epithelial",]))

# Seurat pipeline for SCTransform and clustering
epi <- SCTransform(epi, batch_var = "Source_detail") 
epi <- RunPCA(epi)
epi <- RunUMAP(epi, dims = 1:20, verbose = F)
epi <- FindNeighbors(epi, dims = 1:20)
epi <- FindClusters(epi, resolution = 1.5)
DimPlot(epi, group.by = "dataset")

DotPlot(epi, features = c("FOXJ1","AGER","PDPN","CAV1","EMP2", 
                          "SFTPC","ABCA3","LAMP3","KRT17","KRT5"))
DotPlot(epi, features = c("SCGB3A2","SCGB1A1", "MGP","MUC5B","MUC5AC"))
DotPlot(epi, features = c("EPCAM","PTPRC","PECAM1","LUM"))
DotPlot(epi, features = c("SFTPC","SCGB1A1","SCGB3A2","FOXJ1","KRT5","AGER"))
DotPlot(epi, features = c("MKI67","CHGA","CALCA","FOXI1"))

# Annotation time - Level 2
onion <- as.character(epi@meta.data$seurat_clusters)
onion[onion %in% c(0,1,3,8,18,23,38,44,45,46,50,53,56)] <- "Ciliated Cells"
onion[onion %in% c(2,5,6,7,9,11,12,16,22,24,27,31,33,36,37,40,47,49,55)] <- "AT2"
onion[onion %in% c(14,35)] <- "Transitional AT2"
onion[onion %in% c(4,32,39)] <- "Basal"
onion[onion %in% c(10,13,15,17,19,21,25,28,30,48,51)] <- "Secretory Cells"
onion[onion %in% c(34)] <- "KRT5-/KRT17+"
onion[onion %in% c(26)] <- "AT1"
onion[onion %in% c(57)] <- "PNEC/Ionocytes"
onion[onion %in% c(20,21,29,41,42,43,52,54,58)] <- "Doublets"

epi@meta.data$CellTypeSimple <- onion

# Annotation time - Level 3
onion <- as.character(epi@meta.data$seurat_clusters)
onion[onion %in% c(0,1,3,8,18,23,38,44,45,46,50,53)] <- "Ciliated Cells"
onion[onion == 56] <- "Differentiating Ciliated Cells"
onion[onion %in% c(2,5,6,7,9,11,12,16,22,24,27,31,33,36,37,40,47,49,55)] <- "AT2"
onion[onion %in% c(14,35)] <- "Transitional AT2"
onion[onion %in% c(4,32,39)] <- "Basal"
onion[onion %in% c(15,17,28,48,51)] <- "Club Cells"
onion[onion %in% c(10,13,19,21,25,30)] <- "Goblet Cells"
onion[onion %in% c(34)] <- "KRT5-/KRT17+"
onion[onion %in% c(26)] <- "AT1"
onion[onion %in% c(57)] <- "PNEC/Ionocytes"
onion[onion %in% c(20,21,29,41,42,43,52,54,58)] <- "Doublets"

epi@meta.data$CellType1 <- onion

# Annotation time - Level 4
onion <- as.character(epi@meta.data$seurat_clusters)
onion[onion %in% c(0,1,3,8,18,23,38,44,45,46,50,53)] <- "Ciliated Cells"
onion[onion == 56] <- "Differentiating Ciliated Cells"
onion[onion %in% c(2,5,6,7,9,11,12,16,22,24,27,31,33,36,37,40,47,49,55)] <- "AT2"
onion[onion %in% c(14,35)] <- "Transitional AT2"
onion[onion %in% c(4,32,39)] <- "Basal"
onion[onion %in% c(15,17,51)] <- "SCGB3A2+"
onion[onion %in% c(28,48)] <- "SCGB3A2+/SCGB1A1+"
onion[onion %in% c(21)] <- "MUC5AC+"
onion[onion %in% c(10,13,19,25,30)] <- "MUC5B+"
onion[onion %in% c(34)] <- "KRT5-/KRT17+"
onion[onion %in% c(26)] <- "AT1"
onion[onion %in% c(57)] <- "PNEC/Ionocytes"
onion[onion %in% c(20,21,29,30,41,42,43,52,54,58)] <- "Doublets"

epi@meta.data$CellType2 <- onion

# Remove doublets 
epi2 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$CellTypeSimple != "Doublets", ]))
epi2 <-SCTransform(epi2, batch_var = "Source_detail")
epi2 <- RunPCA(epi2)
epi2 <- RunUMAP(epi2, dims = 1:20, verbose = F)

# Make some plots to check the annotation
DimPlot(epi2, group.by = "CellTypeSimple")
DimPlot(epi2, group.by = "CellType")
DotPlot(epi2, features = c("AGER","PDPN","CAV1","EMP2", "SFTPC","ABCA3","LAMP3",
                          "KRT17","KRT5","FOXJ1","SCGB3A2","SCGB1A1", "MUC5B",
                          "MUC5AC","FOXI1"))

# Subset out the PNEC/Ionocyte cluster for annotation
c57 <- subset(epi2, cells=rownames(epi2@meta.data[epi2@meta.data$CellTypeSimple == "PNEC/Ionocytes",]))
c57 <- FindVariableFeatures(c57, nfeatures = 3000)
c57 <- ScaleData(c57)
c57 <- RunPCA(c57)
c57 <- RunUMAP(c57, dims = 1:5, verbose = F)
c57 <- FindNeighbors(c57, dims = 1:5)
c57 <- FindClusters(c57, resolution = 0.1)
DotPlot(c57, features=c("FOXI1", "ASCL3", "TFCP2L1", "CFTR","CHGA" , "CALCA"))

onion <- as.character(c57@meta.data$seurat_clusters)
onion[onion %in% c(0,2)] <- "PNECs"
onion[onion == 1] <- "Ionocytes"
c57@meta.data$CellType2 <- onion

# Combine the PNEC/Ionocyte cells into the Epi object
epi2 <- subset(epi2, cells=rownames(epi2@meta.data[epi2@meta.data$CellTypeSimple != "PNEC/Ionocytes",]))
epi2 <- merge(epi2, c57)
epi2 <- SCTransform(epi2, batch_var = "Source_detail")
epi2 <- RunPCA(epi2)
epi2 <- RunUMAP(epi2, dims = 1:17, verbose = F)

# Save the object
saveRDS(epi2, file = "20200708_Epi_annotated_noDoublets.rds")
rm(epi)

# ==============================
# MESENCHYMAL POPULATION
# ==============================
# Subset out the mesothelial cells from ILD 
meso <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Mesenchymal", ]))

# Seurat pipeline for SCTransform and clustering
meso <- SCTransform(meso, batch_var = "Source_detail") 
meso <- RunPCA(meso)
meso <- RunUMAP(meso, dims = 1:20, verbose = F)
meso <- FindNeighbors(meso, dims = 1:20)
meso <- FindClusters(meso, resolution = 1.5)
DimPlot(meso)

# Annotation time - Level 2
DotPlot(meso, features = c("MSLN", "UPK3B", "HP", "WT1")) #Mesothelial
DotPlot(meso, features = c("ACTA2", "PDGFRB", "MYH11", "TAGLN", "DES", "ACTG2")) #Smooth Muscle Cells
DotPlot(meso, features =  c("ACTA2", "PDGFRB", "RGS5", "HIGD1B", "GJA4")) #Pericytes
DotPlot(meso, features =  c("LUM", "DCN", "PDGFRA")) #Fibroblasts
DotPlot(meso, features = c("PTPRC","PECAM1","EPCAM"))

onion <- as.character(meso@meta.data$seurat_clusters)
onion[onion %in% c(0,1,2,5,6,7,9,10,11,12,13,15,17,18,23,27)] <- "Fibroblasts"
onion[onion %in% c(3,20)] <- "Smooth Muscle Cells"
onion[onion %in% c(4,26)] <- "Mesothelial"
onion[onion %in% c(8,16,19,21,29)] <- "Pericytes"
onion[onion %in% c(14,22,24,25,28,30)] <- "Doublets"

meso@meta.data$CellTypeSimple <- onion
DimPlot(meso, group.by = "CellTypeSimple")

# Annotation time - Level 3 & 4 
DotPlot(meso, features = c("PLIN2","HAS1", "TWIST1"))
DotPlot(meso, features = c("MYLK","ACTA2","COL8A1", "COL1A1","WNT2"))

onion <- as.character(meso@meta.data$seurat_clusters)
onion[onion %in% c(6)] <- "PLIN2+ Fibroblasts"
onion[onion %in% c(18)] <- "HAS1 High Fibroblasts"
onion[onion == 1] <- "Myofibroblasts"
onion[onion %in% c(0,1,2,5,7,9,10,11,12,13,15,17,23,27)] <- "Fibroblasts"
onion[onion %in% c(3,20)] <- "Smooth Muscle Cells"
onion[onion %in% c(4,26)] <- "Mesothelial"
onion[onion %in% c(8,16,19,21,29)] <- "Pericytes"
onion[onion %in% c(14,22,24,25,28,30)] <- "Doublets"

meso@meta.data$CellType1 <- onion
meso@meta.data$CellType2 <- onion
DimPlot(meso, group.by = "CellType1")

# Remove doublets
meso2 <- subset(meso, cells = rownames(meso@meta.data[meso@meta.data$CellTypeSimple != "Doublets", ]))
meso2 <- SCTransform(meso2, batch_var = "Source_detail") 
meso2 <- RunPCA(meso2)
meso2 <- RunUMAP(meso2, dims = 1:20, verbose = F)
DimPlot(meso2, group.by = "CellType1")

# Save the object
saveRDS(meso2, "20200709_Mesenchymal_noDoublets.rds")
rm(meso)

# ==============================
# ENDOTHELIAL POPULATION
# ==============================
# Subset out the Endothelial cells from the ILD object
endo <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Endothelial", ]))

# Seurat pipeline for SCTransform and clustering
endo <- SCTransform(endo, batch_var = "Source_detail") 
endo <- RunPCA(endo)
endo <- RunUMAP(endo, dims = 1:20, verbose = F)
endo <- FindNeighbors(endo, dims = 1:20)
endo <- FindClusters(endo, resolution = 1)

# Annotation time - Level 2 & 3
DotPlot(endo, features = c("PTPRC", "EPCAM", "PECAM1", "VWF","LUM")) 
DotPlot(endo, features = c("VWF", "CCL21","PROX1","PDPN")) 

onion <- as.character(endo@meta.data$seurat_clusters)
onion[onion %in% c(0,1,2,3,4,5,6,7,9,10,11,12,13,14,16,17,21,25,26,27)] <- "Vascular Endothelial Cells"
onion[onion %in% c(8,18,22,28)] <- "Lymphatic Endothelial Cells"
onion[onion %in% c(15,19,20,23,24,29)] <- "Doublets"
endo@meta.data$CellTypeSimple <- onion
endo@meta.data$CellType1 <- onion
DimPlot(endo, group.by = "CellTypeSimple")

# Annotation time - Level 4
DotPlot(endo, features = c("CA4","VIPR1","RGCC","CYB5A")) #Capillary
DotPlot(endo, features = c("DKK2","GJA5","BMX","HEY1")) # Arterial
DotPlot(endo, features = c("PLA1A","CPE","PTGDS","ACKR1")) # Venous
DotPlot(endo, features = c("SPRY1","PLVAP","COL15A1","VWA1","MYC")) # Bronchial Vessels

onion <- as.character(endo@meta.data$seurat_clusters)
onion[onion %in% c(1,7,9,17,21,27)] <- "Capillary"
onion[onion %in% c(0,16)] <- "Arterial"
onion[onion %in% c(4,6,11,25)] <- "Venous"
onion[onion %in% c(2,3,5,6,10,12,13,14,26)] <- "Bronchial Vessels"
onion[onion %in% c(8,18,22,28)] <- "Lymphatic Endothelial Cells"
onion[onion %in% c(15,19,20,23,24,29)] <- "Doublets"
endo@meta.data$CellType2 <- onion
DimPlot(endo, group.by = "CellType2")

# Remove doublets
endo2 <- subset(endo, cells = rownames(endo@meta.data[endo@meta.data$CellTypeSimple != "Doublets", ]))
endo2 <- SCTransform(endo2, batch_var = "Source_detail") 
endo2 <- RunPCA(endo2)
endo2 <- RunUMAP(endo2, dims = 1:20, verbose = F)
DimPlot(endo2, group.by = "CellType2")

saveRDS(endo2, "20200709_Endothelial_noDoublets.rds")
rm(endo)

# ==============================
# IMMUNE CELLS
# ==============================
# Subset out the Immune cells from ILD object
immune <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Immune", ]))

# Split the immune object into smaller batches for SCTransform (performed on the cluster)
batch <- read.csv("batch_for_sct_immune2.csv", header = T)
immune1 <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$orig.ident %in% batch$Batch_1, ]))
immune2 <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$orig.ident %in% batch$Batch_2, ]))
rm(immune)

# Perform SCTransform for each subset object
immune1 <- SCTransform(immune1, batch_var = "Source_detail") 
immune2 <- SCTransform(immune2, batch_var = "Source_detail") 

saveRDS(immune1, file = "20200606_Immune_subset1_sct.rds")
saveRDS(immune2, file = "20200606_Immune_subset2_sct.rds")

# Merge all object
immune <- merge(x=immune1, y=immune2)
rm(immune1, immune2)

# Seurat pipeline
immune <- NormalizeData(immune)
immune <- FindVariableFeatures(immune, nfeatures = 3000)
immune <- ScaleData(immune)
immune <- RunPCA(immune)
immune <- RunUMAP(immune, dims = 1:20, verbose = F)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 1.5)

DimPlot(immune, group.by = "dataset")
DimPlot(immune)
saveRDS(immune, file = "20200708_Immune_notannotated.rds")

# Annotation time - Level 2
DotPlot(immune, features = c("PTPRC", "EPCAM", "VWF", "ACTA2","PECAM1","MKI67","CD3E"))
DotPlot(immune, features = c("MS4A1","CD19","CD79A")) # B cells
DotPlot(immune, features = c("JCHAIN","IGHG1","IGLL5")) # Plasma cells
DotPlot(immune, features = c("CPA3","KIT")) # Mast cells
DotPlot(immune, features = c("S100A12","FCN1","S100A9","LYZ","CD14","ITGAL")) # Monocytes
DotPlot(immune, features = c("FCER1A", "CD1C", "CLEC9A","LAMP3","PLD4")) # cDCs
DotPlot(immune, features = c("LILRA4","CLEC4C","JCHAIN","IRF8","LILRB4","CD3E","CD4")) # pDCs
DotPlot(immune, features = c("LYZ", "MARCO", "FCGR1A", "C1QA","APOC1")) # Macrophages
DotPlot(immune, features = c("NCR1","KLRB1","NKG7","CD8A","GNLY")) # NK Cells
DotPlot(immune, features = c("CD3E")) # T cells

onion <- as.character(immune@meta.data$seurat_clusters)

onion[onion %in% c(22,44,64) ] <- "B Cells"
onion[onion %in% c(50,55) ] <- "Plasma Cells"
onion[onion %in% c(47,51,61)] <- "Mast Cells"
onion[onion %in% c(7,12,31,34,52,56)] <- "Monocytes"
onion[onion %in% c(8,10,35,46,54)] <- "cDCs"
onion[onion %in% c(60,66)] <- "pDCs"
onion[onion %in% c(0,1,2,3,4,5,6,9,11,15,16,19,21,23,24,26,27,28,29,30,32,33,36,
                   37,38,39,40,42,43,45,48,49,58,62)] <- "Macrophages"
onion[onion %in% c(17,25)] <- "NK Cells"
onion[onion %in% c(13,14,18,20,25,39,41,53,57,59,63)] <- "T Cells"
onion[onion %in% c(48,65)] <- "Doublets"

immune@meta.data$CellTypeSimple <- onion
DimPlot(immune, group.by = "CellTypeSimple")

# Annotation time - Level 3 
DotPlot(immune, features = c("MKI67","CD1","LYZ")) # Proliferating Macrophages
DotPlot(immune, features = c("MKI67","CD1","CD3E")) # Proliferating T Cells
DotPlot(immune, features = c("FOXP3")) # TRegs
DotPlot(immune, features = c("CD8A")) # CD8 T Cells
DotPlot(immune, features = c("CD4","IL7R","CD8A")) # CD4 T Cells

onion <- as.character(immune@meta.data$seurat_clusters)
onion[onion %in% c(22,44,64) ] <- "B Cells"
onion[onion %in% c(50,55) ] <- "Plasma Cells"
onion[onion %in% c(47,51,61)] <- "Mast Cells"
onion[onion %in% c(7,12,31,34,52,56)] <- "Monocytes"
onion[onion %in% c(8,10,35,46,54)] <- "cDCs"
onion[onion %in% c(60,66)] <- "pDCs"
onion[onion %in% c(0,1,2,3,4,5,6,9,11,15,16,19,21,23,24,26,27,28,29,30,32,33,36,
                   38,39,40,42,43,45,48,58,62)] <- "Macrophages"
onion[onion %in% c(17,25)] <- "NK Cells"
onion[onion %in% c(13,14,18,20,25,39,41,57,59)] <- "T Cells"
onion[onion %in% c(48,65)] <- "Doublets"
onion[onion %in% c(37,49)] <- "Proliferating Macrophages"
onion[onion %in% c(53,63)  ] <- "Proliferating T Cells"

immune@meta.data$CellType1 <- onion
DimPlot(immune, group.by = "CellType1")

# Annotation time - Level 4 
DotPlot(immune, features = c("MKI67","CD1","LYZ")) # Proliferating Macrophages
DotPlot(immune, features = c("MKI67","CD1","CD3E")) # Proliferating T Cells
DotPlot(immune, features = c("FOXP3","CD3E")) # TRegs
DotPlot(immune, features = c("CD8A","CD3E")) # CD8 T Cells
DotPlot(immune, features = c("CD4","CD3E")) # CD4 T Cells

onion <- as.character(immune@meta.data$seurat_clusters)
onion[onion %in% c(22,44,64) ] <- "B Cells"
onion[onion %in% c(50,55) ] <- "Plasma Cells"
onion[onion %in% c(47,51,61)] <- "Mast Cells"
onion[onion %in% c(7,12,31,34,52,56)] <- "Monocytes"
onion[onion %in% c(8,10,35,46,54)] <- "cDCs"
onion[onion %in% c(60,66)] <- "pDCs"
onion[onion %in% c(0,1,2,3,4,5,6,9,11,15,16,19,21,23,24,26,27,28,29,30,32,33,36,
                   38,39,40,42,43,45,48,58,62)] <- "Macrophages"
onion[onion %in% c(17,25)] <- "NK Cells"
onion[onion %in% c(48,65)] <- "Doublets"
onion[onion %in% c(37,49)] <- "Proliferating Macrophages"
onion[onion %in% c(53,63)  ] <- "Proliferating T Cells"
onion[onion %in% c(20) ] <- "TRegs"
onion[onion %in% c(13,18,25,39,41) ] <- "CD8 T Cells"
onion[onion %in% c(14,57,59)] <- "CD4 T Cells"

immune@meta.data$CellType2 <- onion
DimPlot(immune, group.by = "CellType2")

DotPlot(immune, features = c("MARCO","CD3E","MS4A1","JCHAIN","CPA3","S100A12",
                             "FCER1A","LILRA4","NKG7","EPCAM","ACTA2"),
        group.by = "CellTypeSimple")
DotPlot(immune, features = c("MARCO","CD3E","MS4A1","JCHAIN","CPA3","S100A12",
                             "FCER1A","LILRA4","NKG7","MKI67","EPCAM","ACTA2",
                             "CD4","CD8A","FOXP3","IL7R"),
        group.by = "CellType2")

# Remove Doublet cells
immune2 <- subset(immune, cells = rownames(immune@meta.data[immune@meta.data$CellTypeSimple != "Doublets", ]))
immune2 <- FindVariableFeatures(immune2, nfeatures = 3000)
immune2 <- ScaleData(immune2)
immune2 <- RunPCA(immune2)
immune2 <- RunUMAP(immune2, dims = 1:20, verbose = F)
DimPlot(immune2, group.by = "CellType2")

# Save objects
saveRDS(immune2, "20200709_Immune_noDoublets.rds")
rm(immune)

# ==========================================
# Merge all objects
# ==========================================
# Without doublets
ild <- merge(x=epi2, y=c(meso2, immune2, endo2))
rm(epi2, endo2, meso2, immune2)
ild <- FindVariableFeatures(ild, nfeatures = 3000)
ild <- ScaleData(ild)
ild <- RunPCA(ild)
ild <- RunUMAP(ild, dims = 1:30, verbose = F)
DimPlot(ild, group.by = "CellType", label = T, repel = T) + NoLegend()

# Save object
saveRDS(ild, file = "/scratch/lbui/Covid19_Seurat/20200709_ILD_noDoublets.rds")

# Get some cell counts
write.csv(table(ild@meta.data$dataset), file = "20200710_ILD_dataset_counts.csv")
write.csv(table(ild@meta.data$Diagnosis2), file = "20200710_ILD_diagnosis_counts.csv")
write.csv(table(ild@meta.data$dataset,ild@meta.data$population), 
          file = "20200710_ILD_pop_counts.csv")
write.csv(table(ild@meta.data$dataset, ild@meta.data$Diagnosis2), 
          file = "20200710_ILD_diagnosis_counts2.csv")
write.csv(table(ild@meta.data$CellType2, ild@meta.data$Diagnosis2), 
          file = "20200712_ILD_diagnosis_counts2.csv")






