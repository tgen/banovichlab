# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2021-02-04
# Title: scRNA-seq data analysis for Covid-19 manuscript (rerun with Seurat v4)
# Note: 
# Seurat v4 does not support batch_var in SCTransform anymore, so I recommend:
# 1. Seurat version v4.0.0 (released 01/30/2021)
# 2. sctransform: remotes::install_github("ChristophH/sctransform@4342164")
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
library(sctransform)
library(Seurat)
library(dplyr)

#=====================================
# Load all the objects 
# These objects were originally downloaded and read into Seurat v3 on 2020/05, I reprocessed  with Seurat v4 to keep everything consistent
# Please contact us directly if you want to use these objects. We do, however, provide a merged object with Level 1 annotation in Figshare  
#=====================================
# Banovich - Kropski dataset
bano <- readRDS("ILD_bano.rds") #GSE135893 and unpublished data
# Yale/BMW dataset
load("Clean.Soup107.forGEO.Seurat.Robj") #GSE136831
# Pittbursg dataset
pitt <- readRDS("Pittsburg.rds") #GSE128033
# Northwestern dataset
northwt <- readRDS("Northwestern.rds") #GSE122960

# ======================================
# Rerun SCTransform and split data set into population and merge individual populations
# ======================================
# Adjust Kaminski object to match with others
kaminski_metadata <- as.data.frame(cbind(intron.exon.107$orig.ident, intron.exon.107$Library_Identity,
                                         intron.exon.107$Disease_Identity, intron.exon.107$Subject_Identity))
colnames(kaminski_metadata) <- c("orig.ident","Library_Identity","Disease_Identity","Subject_Identity")
kaminski <- intron.exon.107
kaminski@meta.data <- kaminski@meta.data[ ,c("Library_Identity","nCount_RNA",
                                             "nFeature_RNA","percent.mito")]
colnames(kaminski@meta.data) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt")
rm(intron.exon.107)

# Add a dataset identity column to each dataset
bano@meta.data$dataset <- "VUMC/TGen"
kaminski@meta.data$dataset <- "Yale/BWH"
northwt@meta.data$dataset <- "Northwestern"
pitt@meta.data$dataset <- "Pittsburg"

# Filter out cells with low number of nfeature_RNA & high mitochondria reads
bano <- PercentageFeatureSet(object = bano, pattern = "^MT-", col.name = "percent.mt")
kaminski <- PercentageFeatureSet(object = kaminski, pattern = "^MT-", col.name = "percent.mt")
northwt <- PercentageFeatureSet(object = northwt, pattern = "^MT-", col.name = "percent.mt")
pitt <- PercentageFeatureSet(object = pitt, pattern = "^MT-", col.name = "percent.mt")

pdf("Smoothscaterplot.pdf")
smoothScatter(bano@meta.data$percent.mt, bano@meta.data$nFeature_RNA) 
smoothScatter(kaminski@meta.data$percent.mt, kaminski@meta.data$nFeature_RNA) 
smoothScatter(northwt@meta.data$percent.mt, northwt@meta.data$nFeature_RNA) 
smoothScatter(pitt@meta.data$percent.mt, pitt@meta.data$nFeature_RNA) 
dev.off()

bano <- subset(bano, subset = nFeature_RNA > 1000 & percent.mt < 25)
kaminski <- subset(kaminski, subset = nFeature_RNA > 1000 & percent.mt < 25)
northwt <- subset(northwt, subset = nFeature_RNA > 1000 & percent.mt < 25)
pitt <- subset(pitt, subset = nFeature_RNA > 1000 & percent.mt < 25)

# Splitting TGEN/VUMC object
bano <- SCTransform(bano)
bano <- RunPCA(bano)
bano <- RunUMAP(bano, dims = 1:20)
bano <- FindNeighbors(bano, dims = 1:20)
bano <- FindClusters(bano, resolution = 0.01)

DotPlot(bano, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(bano@meta.data$seurat_clusters)
onion[onion %in% c(0,3,7,8)] <- "Immune"
onion[onion %in% c(1,2,6)] <- "Epithelial"
onion[onion == 4] <- "Endothelial"
onion[onion == 5] <- "Mesenchymal"
bano@meta.data$population <- onion

# Splitting Yale object
kaminski <- SCTransform(kaminski)
kaminski <- RunPCA(kaminski)
kaminski <- RunUMAP(kaminski, dims = 1:20)
kaminski <- FindNeighbors(kaminski, dims = 1:20)
kaminski <- FindClusters(kaminski, resolution = 0.05)

DotPlot(kaminski, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(kaminski@meta.data$seurat_clusters)
onion[onion %in% c(0,1,2,5,8,9,11)] <- "Immune"
onion[onion %in% c(3,4)] <- "Epithelial"
onion[onion %in%  c(7,10)] <- "Endothelial"
onion[onion == 6] <- "Mesenchymal"
kaminski@meta.data$population <- onion

# Splitting Pittsburg object
pitt <- SCTransform(pitt)
pitt <- RunPCA(pitt)
pitt <- RunUMAP(pitt, dims = 1:20)
pitt <- FindNeighbors(pitt, dims = 1:20)
pitt <- FindClusters(pitt, resolution = 0.01)

DotPlot(pitt, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(pitt@meta.data$seurat_clusters)
onion[onion %in% c(0,4,6)] <- "Immune"
onion[onion == 1] <- "Epithelial"
onion[onion %in% c(2,5)] <- "Endothelial"
onion[onion == 3] <- "Mesenchymal"
pitt@meta.data$population <- onion

# Splitting Northwestern object
northwt <- SCTransform(northwt)
northwt <- RunPCA(northwt)
northwt <- RunUMAP(northwt, dims = 1:20)
northwt <- FindNeighbors(northwt, dims = 1:20)
northwt <- FindClusters(northwt, resolution = 0.02)

DotPlot(northwt, features = c("EPCAM","PTPRC","PECAM1"))

onion <- as.character(northwt@meta.data$seurat_clusters)
onion[onion %in% c(0,3)] <- "Immune"
onion[onion %in% c(1,2)] <- "Epithelial"
onion[onion == 4] <- "Endothelial"
onion[onion == 5] <- "Mesenchymal"
northwt@meta.data$population <- onion

# Merge all dataset
ild <- merge(x=bano, y=c(kaminski, pitt, northwt))

# ===================================
# Add meta data and save objects
# ===================================
# Add meta data
meta.data <- read.csv("Sample_metadata_information.csv",
                      header = T )
ild@meta.data$Diagnosis <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                           from = meta.data$Library_ID,
                                           to = as.character(meta.data$Diagnosis))
ild@meta.data$Status <- plyr::mapvalues(x = ild@meta.data$orig.ident,
                                        from = meta.data$Library_ID,
                                        to = as.character(meta.data$Status))
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

# Save the object
saveRDS(bano, file = "BanoVUMC_pop.rds")
saveRDS(kaminski, file = "Kaminski_pop.rds")
saveRDS(pitt, file = "Pittsburg_pop.rds")
saveRDS(northwt, file = "Northwestern_pop.rds")
saveRDS(ild, file = "ILD_alldataset_population_annotated.rds")

# Remove files to free memory
rm(bano,kaminski,pitt,northwt)

# Remove the SCT assay to reduce size of the file for Figshare
ild@assays$SCT <- NULL
saveRDS(ild, file="ILD_alldataset_population_noSCT.rds")

# ==============================================================================
# Running cell type annotation per population (start from here with the submitted object)
# We provided this object through Figshare (see link on Github page)
# ==============================================================================
# Read in the ILD object
ild <- readRDS("ILD_alldataset_population_noSCT.rds")

# ---------------------------------------
# Epithelial population
# ---------------------------------------
# Subset out the epithelial cells from ILD object
epi <- subset(ild, cells=rownames(ild@meta.data[ild@meta.data$population == "Epithelial",]))

# Seurat pipeline for SCTransform and clustering
epi <- SCTransform(epi, batch_var="dataset")
epi <- RunPCA(epi)
epi <- RunUMAP(epi, dims = 1:20, verbose = F)
epi <- FindNeighbors(epi, dims = 1:20)
epi <- FindClusters(epi, resolution = 1.5)
DimPlot(epi, group.by = "dataset")

DotPlot(epi, features = c("FOXJ1","AGER","PDPN","CAV1","EMP2", 
                          "SFTPC","ABCA3","LAMP3","KRT17","KRT5"))
DotPlot(epi, features = c("SCGB3A2","SCGB1A1", "MGP","MUC5B","MUC5AC"))
DotPlot(epi, features = c("EPCAM","PTPRC","PECAM1","LUM","VWF"))
DotPlot(epi, features = c("SFTPC","SCGB1A1","SCGB3A2","FOXJ1","KRT5","AGER"))
DotPlot(epi, features = c("CHGA","CALCA","FOXI1"))

# Annotation time - Level 2
onion <- as.character(epi@meta.data$seurat_clusters)
onion[onion %in% c(2,3,7,11,12,17,20,24,29,35,40,49,51,52,56,66)] <- "Ciliated Cells"
onion[onion %in% c(0,1,4,5,6,8,15,16,19,22,23,25,26,43,45,47,54,59)] <- "AT2"
onion[onion %in% c(21)] <- "Transitional AT2"
onion[onion %in% c(9,33,42)] <- "Basal"
onion[onion %in% c(10,13,18,27,28,30,31,34,36,37,38,44,46,50,53,57,64,65)] <- "Secretory Cells"
onion[onion %in% c(39,41)] <- "KRT5-/KRT17+"
onion[onion %in% c(32,60,61)] <- "AT1"
onion[onion %in% c(55)] <- "PNEC/Ionocytes"
onion[onion %in% c(14,48,58,62,63)] <- "Doublets"

epi@meta.data$CellTypeSimple <- onion

# Annotation time - Level 3
onion <- as.character(epi@meta.data$seurat_clusters)
onion[onion %in% c(2,3,7,11,12,17,20,24,29,35,40,49,51,52,56,66)] <- "Ciliated Cells"
onion[onion %in% c(0,1,4,5,6,8,15,16,19,22,23,25,26,43,45,47,54,59)] <- "AT2"
onion[onion %in% c(21)] <- "Transitional AT2"
onion[onion %in% c(9,33,42)] <- "Basal"
onion[onion %in% c(10,13,23,36,50,18,30,31,37,65)] <- "Club Cells"
onion[onion %in% c(46,64,27,28,34,38,44,48,53,57,64)] <- "Goblet Cells"
onion[onion %in% c(39,41)] <- "KRT5-/KRT17+"
onion[onion %in% c(32,60,61)] <- "AT1"
onion[onion %in% c(55)] <- "PNEC/Ionocytes"
onion[onion %in% c(14,48,58,62,63)] <- "Doublets"

epi@meta.data$CellType1 <- onion

# Annotation time - Level 4
onion <- as.character(epi@meta.data$seurat_clusters)
onion[onion %in% c(2,3,7,11,12,17,20,24,29,35,40,49,51,52,56,66)] <- "Ciliated Cells"
onion[onion %in% c(0,1,4,5,6,8,15,16,19,22,23,25,26,43,45,47,54,59)] <- "AT2"
onion[onion %in% c(21)] <- "Transitional AT2"
onion[onion %in% c(9,33,42)] <- "Basal"
onion[onion %in% c(10,13,23,36,50)] <- "SCGB3A2+"
onion[onion %in% c(18,30,31,37,65)] <- "SCGB3A2+/SCGB1A1+"
onion[onion %in% c(46,64)] <- "MUC5AC+"
onion[onion %in% c(27,28,34,38,44,48,53,57,64)] <- "MUC5B+"
onion[onion %in% c(39,41)] <- "KRT5-/KRT17+"
onion[onion %in% c(32,60,61)] <- "AT1"
onion[onion %in% c(55)] <- "PNEC/Ionocytes"
onion[onion %in% c(14,48,58,62,63)] <- "Doublets"

epi@meta.data$CellType2 <- onion

# Remove doublets 
epi2 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$CellTypeSimple != "Doublets", ]))
epi2 <- RunUMAP(epi2, dims=1:20)

# Save the object
saveRDS(epi2, file="Epithelial_noDoublets.rds")
saveRDS(epi, file="Epithelial_wDoublets.rds")

# ---------------------------------------
# Mesenchymal population
# ---------------------------------------
# Subset out the mesothelial cells from ILD 
meso <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Mesenchymal", ]))

# Seurat pipeline for SCTransform and clustering
meso <- SCTransform(meso, batch_var = "dataset") 
meso <- RunPCA(meso)
meso <- RunUMAP(meso, dims = 1:20, verbose = F)
meso <- FindNeighbors(meso, dims = 1:20)
meso <- FindClusters(meso, resolution = 1.5)
DimPlot(meso)
DimPlot(meso, group.by = "dataset2")

# Annotation time - Level 2
DotPlot(meso, features = c("MSLN", "UPK3B", "HP", "WT1")) #Mesothelial
DotPlot(meso, features = c("ACTA2", "PDGFRB", "MYH11", "TAGLN", "DES", "ACTG2")) #Smooth Muscle Cells
DotPlot(meso, features =  c("ACTA2", "PDGFRB", "RGS5", "HIGD1B", "GJA4")) #Pericytes
DotPlot(meso, features =  c("LUM", "DCN", "PDGFRA")) #Fibroblasts
DotPlot(meso, features = c("PTPRC","PECAM1","EPCAM"))

onion <- as.character(meso@meta.data$seurat_clusters)
onion[onion %in% c(0,1,3,4,6,7,8,9,10,13,14,15,17,18,19,23,27)] <- "Fibroblasts"
onion[onion %in% c(2,20,21,25,28)] <- "Smooth Muscle Cells"
onion[onion %in% c(5,22)] <- "Mesothelial"
onion[onion %in% c(11,16,24)] <- "Pericytes"
onion[onion %in% c(12,26,29)] <- "Doublets"

meso@meta.data$CellTypeSimple <- onion
DimPlot(meso, group.by = "CellTypeSimple")

# Annotation time - Level 3 & 4 
DotPlot(meso, features = c("PLIN2","HAS1", "TWIST1"))
DotPlot(meso, features = c("MYLK","ACTA2","COL8A1", "COL1A1","WNT2"))

onion <- as.character(meso@meta.data$seurat_clusters)
onion[onion %in% c(3)] <- "PLIN2+ Fibroblasts"
onion[onion %in% c(5)] <- "HAS1 High Fibroblasts"
onion[onion %in% c(1,15) ] <- "Myofibroblasts"
onion[onion %in% c(0,4,6,7,8,9,10,13,14,17,18,19,23,27)] <- "Fibroblasts"
onion[onion %in% c(2,20,21,25,28)] <- "Smooth Muscle Cells"
onion[onion %in% c(5,22)] <- "Mesothelial"
onion[onion %in% c(11,16,24)] <- "Pericytes"
onion[onion %in% c(12,26,29)] <- "Doublets"

meso@meta.data$CellType1 <- onion
meso@meta.data$CellType2 <- onion
DimPlot(meso, group.by = "CellType1")
meso@meta.data$CellType2 <- meso@meta.data$CellType1

# Remove doublets
meso2 <- subset(meso, cells = rownames(meso@meta.data[meso@meta.data$CellTypeSimple != "Doublets", ]))
DimPlot(meso2, group.by = "CellType1")

# Save the object
saveRDS(meso2, file="Mesenchymal_noDoublets.rds")
saveRDS(meso, file="Mesenchymal_wDoublets.rds")
rm(meso)

# ---------------------------------------
# Endothelial population
# ---------------------------------------
# Subset out the Endothelial cells from the ILD object
endo <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Endothelial", ]))

# Seurat pipeline for SCTransform and clustering
endo <- SCTransform(endo, batch_var="dataset") 
endo <- RunPCA(endo)
endo <- RunUMAP(endo, dims = 1:20, verbose = F)
endo <- FindNeighbors(endo, dims = 1:20)
endo <- FindClusters(endo, resolution = 1)

# Annotation time - Level 2 & 3
DotPlot(endo, features = c("PTPRC", "EPCAM", "PECAM1", "VWF","LUM")) 
DotPlot(endo, features = c("VWF", "CCL21","PROX1","PDPN")) 

onion <- as.character(endo@meta.data$seurat_clusters)
onion[onion %in% c(0,1,2,3,4,5,6,7,8,9,13,14,16,17,18,19,20,23,24,25,26)] <- "Vascular Endothelial Cells"
onion[onion %in% c(10,11,12,22)] <- "Lymphatic Endothelial Cells"
onion[onion %in% c(15,21,27)] <- "Doublets"
endo@meta.data$CellTypeSimple <- onion
endo@meta.data$CellType1 <- onion
DimPlot(endo, group.by = "CellTypeSimple")

# For Endothelial population, stop here
endo@meta.data$CellType1 <- endo@meta.data$CellTypeSimple
endo@meta.data$CellType2 <- endo@meta.data$CellTypeSimple

# Remove doublets
endo2 <- subset(endo, cells = rownames(endo@meta.data[endo@meta.data$CellTypeSimple != "Doublets", ]))

saveRDS(endo2, file="Endothelial_noDoublets.rds")
saveRDS(endo, file="Endothelial_wDoublets.rds")
rm(endo)

# ---------------------------------------
# Immune population
# ---------------------------------------
# Subset out the Immune cells from ILD object
immune <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Immune", ]))

# Perform SCTransform for each subset object
immune <- SCTransform(immune, batch_var = "dataset") 
immune <- RunPCA(immune)
immune <- RunUMAP(immune, dims = 1:20, verbose = F)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 1.5)
saveRDS(immune, file = "20210827_Immune_SCT.rds")

DimPlot(immune, group.by = "dataset")
DimPlot(immune)

# Annotation time - Level 2
DotPlot(immune, features = c("PTPRC", "EPCAM", "VWF", "ACTA2","PECAM1","LUM"))
DotPlot(immune, features = c("MS4A1","CD19","CD79A")) # B cells
DotPlot(immune, features = c("JCHAIN","IGHG1","IGLL5")) # Plasma cells
DotPlot(immune, features = c("CPA3","KIT")) # Mast cells
DotPlot(immune, features = c("S100A12","FCN1","S100A9","LYZ","CD14","ITGAL")) # Monocytes
DotPlot(immune, features = c("FCER1A", "CD1C", "CLEC9A")) # cDCs
DotPlot(immune, features = c("LILRA4","CLEC4C","JCHAIN")) # pDCs
DotPlot(immune, features = c("LYZ", "MARCO", "FCGR1A", "C1QA","APOC1")) # Macrophages
DotPlot(immune, features = c("CD3E", "KLRB1", "NKG7", "NCR1","GNLY","CD8A")) # NK Cells
DotPlot(immune, features = c("CD3E", "CD8A")) # T cells

onion <- as.character(immune@meta.data$seurat_clusters)

onion[onion %in% c(22,42,56) ] <- "B Cells"
onion[onion %in% c(47,50,53) ] <- "Plasma Cells"
onion[onion %in% c(41,48,69)] <- "Mast Cells"
onion[onion %in% c(52,57,64,66)] <- "Monocytes"
onion[onion %in% c(18,32,39,40)] <- "cDCs"
onion[onion %in% c(45)] <- "pDCs"
onion[onion %in% c(0,1,2,3,7,8,9,10,11,12,13,14,15,16,19,20,21,23,24,25,26,27,28,
                   29,30,31,33,36,37,38,44,48,51,54,55,57,58,59,60,61,62,63,65,67,
                   68)] <- "Macrophages"
onion[onion %in% c(17)] <- "NK Cells"
onion[onion %in% c(4,5,6,34,35,43,46)] <- "T Cells"
onion[onion %in% c(49)] <- "Doublets"

immune@meta.data$CellTypeSimple <- onion
DimPlot(immune, group.by = "CellTypeSimple")

# Annotation time - Level 3 
DotPlot(immune, features = c("MKI67","LYZ","CD3E")) # Proliferating Cells

onion <- as.character(immune@meta.data$seurat_clusters)
onion[onion %in% c(22,42,56) ] <- "B Cells"
onion[onion %in% c(47,50,53) ] <- "Plasma Cells"
onion[onion %in% c(41,48,69)] <- "Mast Cells"
onion[onion %in% c(52,57,64,66)] <- "Monocytes"
onion[onion %in% c(18,32,39,40)] <- "cDCs"
onion[onion %in% c(45)] <- "pDCs"
onion[onion %in% c(0,1,2,3,7,8,9,10,11,12,13,14,15,16,19,20,21,23,24,25,26,27,28,
                   29,30,31,33,36,37,38,48,51,55,57,58,59,60,61,62,63,65,67,
                   68)] <- "Macrophages"
onion[onion %in% c(17)] <- "NK Cells"
onion[onion %in% c(4,5,6,34,35,43,46)] <- "T Cells"
onion[onion %in% c(49)] <- "Doublets"
onion[onion %in% c(44,54)] <- "Proliferating Macrophages"

immune@meta.data$CellType1 <- onion
DimPlot(immune, group.by = "CellType1")

# Annotation time - Level 4 

## Subset out T cells
tcell <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$CellTypeSimple == "T Cells",]))
tcell <- FindVariableFeatures(tcell, nfeatures=3000)
tcell <- ScaleData(tcell)
tcell <- RunPCA(tcell)
tcell <- RunUMAP(tcell,dims=1:14)
tcell <- FindNeighbors(tcell, dims=1:14)
tcell <- FindClusters(tcell, resolution = 0.5)

DotPlot(tcell, features = c("FOXP3","CD3E","CD8A","CD4","MKI67")) 

onion <- as.character(tcell@meta.data$seurat_clusters)

onion[onion %in% c(8,11)] <- "Proliferating T Cells"
onion[onion %in% c(0,1,2,3) ] <- "CD8 T Cells"
onion[onion %in% c(4,5,6,7,9,10)] <- "CD4 T Cells"
tcell$CellType2 <- onion

# Add the Tcells annotation to the immune meta.data
immune@meta.data$CellType2 <- plyr::mapvalues(x = rownames(immune@meta.data),
                                              from = rownames(tcell@meta.data),
                                              to = as.character(tcell@meta.data$CellType2))

onion <- as.character(immune@meta.data$CellType1)
onion1 <- as.character(immune@meta.data$CellType2)

output <- ifelse(onion1 %in% as.character(tcell@meta.data$CellType2), onion1, onion)

immune@meta.data$CellType2 <- output

DotPlot(immune, features = c("MARCO","CD3E","MS4A1","JCHAIN","CPA3","S100A12",
                             "FCER1A","LILRA4","NKG7","EPCAM","ACTA2"),
        group.by = "CellTypeSimple")
DotPlot(immune, features = c("MARCO","CD3E","MS4A1","JCHAIN","CPA3","S100A12",
                             "FCER1A","LILRA4","NKG7","MKI67","EPCAM","ACTA2",
                             "CD4","CD8A","FOXP3","IL7R"),
        group.by = "CellType2")

# Remove Doublet cells
immune2 <- subset(immune, cells = rownames(immune@meta.data[immune@meta.data$CellTypeSimple != "Doublets", ]))
DimPlot(immune2, group.by = "CellType2")

# Save objects
saveRDS(immune2, file= "Immune_noDoublets.rds")
saveRDS(immune, file="Immune_wDoublets.rds")
rm(immune)

# ==========================================
# Merge all objects
# Note:  Seurat v4 does not support batch_var in SCTransform, solution:
# 1. Downgrade to Seurat v3 for merging (I used v3.1.4)
# 2. Install SeuratObject separately for v3 through CRAN
# ==========================================
# Set DefaultAssay to RNA for merging
DefaultAssay(epi2) <- "RNA"
DefaultAssay(endo2) <- "RNA"
DefaultAssay(meso2) <- "RNA"
DefaultAssay(immune2) <- "RNA"

# Without doublets
ild <- merge(x=epi2, y=c(meso2, immune2, endo2))
rm(epi2, endo2, meso2, immune2)
ild <- NormalizeData(ild)
ild <- FindVariableFeatures(ild, nfeatures = 3000, selection.method = "vst")
ild <- ScaleData(ild)
ild <- RunPCA(ild)
ild <- RunUMAP(ild, dims = 1:28, verbose = F)
DimPlot(ild, group.by = "CellType2", label = T, repel = T) + NoLegend()

# Save object
saveRDS(ild, "ILD_noDoublets.rds")

# Get some cell counts
write.csv(table(ild@meta.data$dataset), file = "20200710_ILD_dataset_counts.csv")
write.csv(table(ild@meta.data$Diagnosis2), file = "20200710_ILD_diagnosis_counts.csv")
write.csv(table(ild@meta.data$dataset,ild@meta.data$population), 
          file = "20200710_ILD_pop_counts.csv")
write.csv(table(ild@meta.data$dataset, ild@meta.data$Diagnosis2), 
          file = "20200710_ILD_diagnosis_counts2.csv")
write.csv(table(ild@meta.data$CellType2, ild@meta.data$Diagnosis2), 
          file = "20200712_ILD_diagnosis_counts2.csv")






