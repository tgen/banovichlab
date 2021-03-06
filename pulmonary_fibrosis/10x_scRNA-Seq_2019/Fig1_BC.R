# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
#             Austin J. Gutierrez, agutierrez@tgen.org
#             Arun C. Habermann, arun.c.habermann@vumc.org
#             Stephanie L. Yahn, syahn@tgen.org
# Date: 11/06/2019
# Description: Sample proccessing workflow for IPF project in Seurat V3
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
set.seed(8211)

# ======================================
# Load libraries
# ======================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(ade4)
library(Matrix)

# ======================================
# Define the location of the Cellranger output files to be read into R
# Define the location of barcodes.tsv, genes.tsv, matrix.mtx
# Create a character vector of sample file names
# Read data from each sample folder and combine
# ======================================
# dataset_loc <- "/scratch/agutierrez/10x_fastq/Outs/IPF/"
# inputfiles_loc <- "outs/filtered_feature_bc_matrix"
# batchids <- list.files(path = dataset_loc, pattern = "IPF*",
#                        full.names = FALSE, recursive = FALSE,
#                        ignore.case = FALSE, include.dirs = FALSE)
# 
# d10x.data <- sapply(batchids,  function(i){
#   sample_id <- sapply(strsplit(i,split="_") , "[[", 8)
#   d10x <- Read10X(file.path(dataset_loc,i,inputfiles_loc))
#   colnames(d10x) <- paste(sample_id, sep="_", colnames(d10x))
#   d10x
# })

# ======================================
# do.call is required to create a single S4 object of class dgCMatrix
# (genes x cells) otherwise, cbind will create a 2x1 matrix list of two 
# separate S4 objects of class dgCMatrix (genes x cells)
# ======================================
# experiment.data <- do.call("cbind", d10x.data)

# Save the dgCMatrix
# writeMM(experiment.data, file = "matrix.mtx")
# write(x = rownames(experiment.data), file = "genes.tsv")
# write(x = colnames(experiment.data), file = "barcodes.tsv")

# =============================================================================
# Start from here, point data.dir to the directory with the downloaded raw data
# =============================================================================
# Create Seurat Object from the dgCMatrix
# Read in the dgCMatrix file:
# ======================================
ild_data <- Read10X(data.dir = ".",
                    gene.column = 1)

# ======================================
# Keep all cells with at least 200 detected genes,
# ======================================
ild <- CreateSeuratObject(ild_data,
                          project = "scRNA lung",
                          min.features = 200)

# ======================================
# Add more information into the meta.data
# ======================================
# meta.data <- read.csv("/scratch/lbui/20190623_Final_version/IPF_master.csv", header = T)
# 
# ild@meta.data$Diagnosis <- plyr::mapvalues(x = ild@meta.data$orig.ident,
#                                           from = meta.data$Library_ID,
#                                           to = as.character(meta.data$Diagnosis))
# 
# ild@meta.data$Sample_Name <- plyr::mapvalues(x = ild@meta.data$orig.ident,
#                                             from = meta.data$Library_ID,
#                                             to = as.character(meta.data$Sample_Name))
# 
# ild@meta.data$Sample_Source <- plyr::mapvalues(x = ild@meta.data$orig.ident,
#                                               from = meta.data$Library_ID,
#                                               to = as.character(meta.data$Sample_Source))
# 
# ild@meta.data$Status <- plyr::mapvalues(x = ild@meta.data$orig.ident,
#                                        from = meta.data$Library_ID,
#                                        to = as.character(meta.data$Status))

# ======================================
# QC Seurat object
# ======================================
ild <- PercentageFeatureSet(object = ild, pattern = "^MT-", col.name = "percent.mt")
pdf("QC_plot.pdf")
VlnPlot(ild, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Make density plots
smoothScatter(ild@meta.data$percent.mt, ild@meta.data$nFeature_RNA)
abline(h = 1000, v = 1.5)
text(1.5,700, "nFeature = 1000, percent.mt = 1.5", adj = c(0, -.1))
dev.off()

# Filter out cells with less than 1000 nFeature and more than 25% percent.mt
ild <- subset(ild, subset = nFeature_RNA > 1000 & percent.mt < 25)

# ======================================
# Seurat object processing
# ======================================  
# Run SCTransform
ild <- SCTransform(object = ild, verbose = T)

# Chose a PC for analysis
ild <- FindVariableFeatures(ild, verbose = T, nfeatures = 3000)
ild <- ScaleData(ild, features = row.names(ild@assays$SCT@data))
ild <- RunPCA(ild)

potato <-  matrix(ncol = 17, nrow = ncol(ild))
potato2 <-  matrix(ncol = 17, nrow = ncol(ild))

j <- 0
for(i in 9:25) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(ild, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 9:25
colnames(potato2) <- 9:25
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(ild), 3000), ]
potato2_new <- potato2[sample(1:ncol(ild), 3000), ]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("ild_PCA.pdf")
plot(c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
     mrand_obs[, 1], ylim = c(0.8, 1))
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 1]))
points(c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
       mrand_obs[, 2], ylim = c(0.8, 1), col = "red")
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(ild)
dev.off()

# ======================================
# RunUMAP, FindClusters, FindNeighbors 
# using the PC from the analysis above
# ======================================
ild <- RunUMAP(object = ild, dims = 1:20, verbose = F)
ild <- FindNeighbors(object = ild, dims = 1:20, verbose = F)
ild <- FindClusters(object = ild, resolution = 0.01, verbose = F)

# ======================================
# Annotate clusters using markers specific for populations
# Immune cells: PTPRC
# Epithelial cells: EPCAM
# Endothelial cells: PECAM1+ PTPRC-
# Mesenchymal cells: no expression of those markers
# ======================================

pdf("ild.pdf")
DimPlot(ild)
FeaturePlot(ild, c("PTPRC", "EPCAM", "PECAM1"))
DotPlot(ild, features = c("PTPRC", "EPCAM", "PECAM1"))
dev.off()

Idents(ild, cells = WhichCells(ild, idents = c(0, 3))) <- "Immune"
Idents(ild, cells = WhichCells(ild, idents = c(1, 2))) <- "Epithelial"
Idents(ild, cells = WhichCells(ild, idents = c(4, 6))) <- "Endothelial"
Idents(ild, cells = WhichCells(ild, idents = c(5))) <- "Mesenchymal"
ild@meta.data$population <- Idents(ild)

# ======================================
# Subset into 4 major populations
# ======================================
immune <- subset(ild, idents = "Immune")
epi <- subset(ild, idents = "Epithelial")
endo <- subset(ild, idents = "Endothelial")
meso <- subset(ild, idents = "Mesenchymal")

# Confirm subsetted data
table(ild@meta.data$population)
table(immune@meta.data$population)
table(epi@meta.data$population)
table(endo@meta.data$population)
table(meso@meta.data$population)

# ======================================
# Make plots to determine which PC 
# to use for each subset
# Mesenchymal subset
# ======================================
meso <- FindVariableFeatures(meso, verbose = T, nfeatures = 3000)
meso <- ScaleData(meso,features = row.names(meso@assays$SCT@data))
meso <- RunPCA(meso)
potato <- matrix(ncol = 20, nrow = ncol(meso))
potato2 <- matrix(ncol = 20, nrow = ncol(meso))

j <- 0
for(i in 6:25) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(meso, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 6:25
colnames(potato2) <- 6:25
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(meso), 3000), ]
potato2_new <- potato2[sample(1:ncol(meso), 3000), ]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("meso_PCA.pdf")
plot(c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
     mrand_obs[, 1], ylim = c(0, 1))
lines(spline(x = c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 1]))
points(c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
       mrand_obs[, 2], ylim = c(0, 1), col = "red")
lines(spline(x = c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 2]), col = "red")

dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

# ======================================
# Endothelial subset
# ======================================
endo <- FindVariableFeatures(endo, verbose = T, nfeatures = 3000)
endo <- ScaleData(endo,features = row.names(endo@assays$SCT@data))
endo <- RunPCA(endo)
potato <- matrix(ncol = 12, nrow = ncol(endo))
potato2 <- matrix(ncol = 12, nrow = ncol(endo))

j <- 0
for(i in 9:20) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(endo, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 9:20
colnames(potato2) <- 9:20
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(endo), 3000),]
potato2_new <- potato2[sample(1:ncol(endo), 3000),]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("endo_PCA.pdf")
plot(c(9,10,11,12,13,14,15,16,17,18,19),
     mrand_obs[, 1], ylim = c(.75, 1))
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19),
             y = mrand_obs[, 1]))
points(c(9,10,11,12,13,14,15,16,17,18,19),
       mrand_obs[, 2], ylim = c(.75, 1), col = "red")
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(endo)
dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

# ======================================
# Epithelial subset
# ======================================
epi <- FindVariableFeatures(epi, verbose = T, nfeatures = 3000)
epi <- ScaleData(epi,features = row.names(epi@assays$SCT@data))
epi <- RunPCA(epi)
potato <- matrix(ncol = 15, nrow = ncol(epi))
potato2 <- matrix(ncol = 15, nrow = ncol(epi))

j <- 0
for(i in 7:21) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(epi, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 7:21
colnames(potato2) <- 7:21
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(epi), 3000),]
potato2_new <- potato2[sample(1:ncol(epi), 3000),]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("epi_PCA.pdf")
plot(c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
     mrand_obs[, 1], ylim = c(.75, 1))
lines(spline(x = c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
             y = mrand_obs[, 1]))
points(c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
       mrand_obs[, 2], ylim = c(.75, 1), col = "red")
lines(spline(x = c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(epi, ndims = 30)
dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

# ======================================
# Immune subset
# ======================================
immune <- FindVariableFeatures(immune, verbose = T, nfeatures = 3000)
immune <- ScaleData(immune, features = row.names(immune@assays$SCT@data))
immune <- RunPCA(immune)
potato <- matrix(ncol = 15, nrow = ncol(immune))
potato2 <- matrix(ncol = 15, nrow = ncol(immune))

j <- 0
for(i in 13:27) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(immune, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, ] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 13:27
colnames(potato2) <- 13:27
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(immune), 3000),]
potato2_new <- potato2[sample(1:ncol(immune), 3000),]

for(i in  1:c(dim(potato_new)[2] - 1)){
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("immune_PCA.pdf")
plot(c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
     mrand_obs[, 1], ylim = c(0, 1))
lines(spline(x = c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
             y = mrand_obs[, 1]))
points(c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
       mrand_obs[, 2], ylim = c(0, 1), col = "red")
lines(spline(x = c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(immune, ndims = 30)
dev.off()

# ======================================
# Run UMAP for the Seurat subset
# objects with the PCs determined from above
# ======================================
epi <- RunUMAP(object = epi, dims = 1:11, verbose = F)
endo <- RunUMAP(object = endo, dims = 1:14, verbose = F)
immune <- RunUMAP(object = immune, dims = 1:21, verbose = F)
meso <- RunUMAP(object = meso, dims = 1:14, verbose = F)

# ======================================
# Epithelial
# ======================================
runFeaturePlots = TRUE
epi <- FindNeighbors(epi, dims = 1:20)
epi <- FindClusters(epi, resolution = 2)

if (runFeaturePlots) {FeaturePlot(epi, features = c("SFTPC", "ABCA3", "AGER", "HOPX"))}
Idents(epi, cells = WhichCells(epi, idents = c(19))) <- "AT1"
Idents(epi, cells = WhichCells(epi, idents = c(3, 9, 10, 13, 18, 22, 25, 28, 29))) <- "AT2"
Idents(epi, cells = WhichCells(epi, idents = c(17, 31))) <- "Transitional AT2"

if (runFeaturePlots) {FeaturePlot(epi, features = c("KRT5", "KRT17", "TP63", "COL1A1"))}
Idents(epi, cells = WhichCells(epi, idents = c(6))) <- "Basal"
Idents(epi, cells = WhichCells(epi, idents = c(23))) <- "KRT5-/KRT17+"

if (runFeaturePlots) {FeaturePlot(epi, features = c("SCGB1A1", "SCGB3A2", "MUC5B", "MUC5AC"))}
Idents(epi, cells = WhichCells(epi, idents = c(8, 20))) <- "MUC5B+"
Idents(epi, cells = WhichCells(epi, idents = c(0, 24))) <- "SCGB3A2+"
if (runFeaturePlots) {FeaturePlot(epi, features = c("MGP", "SCGB1A1", "SCGB3A2"))}
Idents(epi, cells = WhichCells(epi, idents = c(12))) <- "SCGB3A2+ SCGB1A1+"

if (runFeaturePlots) {FeaturePlot(epi, features = c("FOXI1", "CHGA"))}
Idents(epi, cells = WhichCells(epi, idents = c(38))) <- "PNECs/Ionocytes"

if (runFeaturePlots) {FeaturePlot(epi, features = c("FOXJ1", "TMEM190", "CAPS", "HYDIN"))}
Idents(epi, cells = WhichCells(epi, idents = c(1, 2, 4, 5, 7, 11, 16, 21, 27, 35))) <- "Ciliated"
Idents(epi, cells = WhichCells(epi, idents = c(14, 34))) <- "Differentiating Ciliated"

if (runFeaturePlots) {FeaturePlot(epi, features = c("EPCAM", "PTPRC", "PECAM1", "LUM"))}
if (runFeaturePlots) {FeaturePlot(epi, features = c("SFTPC", "SCGB1A1", "SCGB3A2", "FOXJ1", "KRT5", "AGER"))}
Idents(epi, cells = WhichCells(epi, idents = c(15, 26, 32, 33, 36, 37))) <- "Doublets"

if (runFeaturePlots) {FeaturePlot(epi, features = c("MKI67"))}
Idents(epi, cells = WhichCells(epi, idents = c(30))) <- "Proliferating Epithelial Cells"

epi$firstCT <- Idents(epi)

pnecionocyte_subset <- subset(epi, idents = c("PNECs/Ionocytes"))
if (runFeaturePlots) {FeaturePlot(pnecionocyte_subset, features = c("FOXI1", "CHGA"))}
Idents(epi, cells = WhichCells(pnecionocyte_subset, expression = CHGA > 0)) <- "PNECs"
Idents(epi, cells = WhichCells(pnecionocyte_subset, expression = FOXI1 > 0)) <- "Ionocytes"
Idents(epi, cells = WhichCells(pnecionocyte_subset, idents = "PNECs/Ionocytes")) <- "Basal"

muc5b_subset <- subset(epi, idents = c("MUC5B+")) %>% FindNeighbors() %>% FindClusters(resolution = 2)
if (runFeaturePlots) {FeaturePlot(muc5b_subset, features = c("MUC5B", "MUC5AC"))}
Idents(epi, cells = WhichCells(muc5b_subset, idents = c(15))) <- "MUC5AC+ High"

epi$celltype <- Idents(epi)
DimPlot(epi)

keep <- c("AT2","Differentiating Ciliated","PNECs","MUC5AC+ High","Transitional AT2", "SCGB3A2+ SCGB1A1+", "AT1", "Ionocytes", "MUC5B+","SCGB3A2+","Basal","Ciliated","KRT5-/KRT17+", "Proliferating Epithelial Cells")
epi <- subset(epi, cells = row.names(epi@meta.data[epi@meta.data$celltype %in% keep, ]))
DimPlot(epi, group.by = "celltype")

potato <- matrix(ncol = 15, nrow = ncol(epi))
potato2 <- matrix(ncol = 15, nrow = ncol(epi))

j <- 0
for(i in 7:21) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(epi, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 7:21
colnames(potato2) <- 7:21
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(epi), 3000),]
potato2_new <- potato2[sample(1:ncol(epi), 3000),]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("epi_PCA.pdf")
plot(c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
     mrand_obs[, 1], ylim = c(.75, 1))
lines(spline(x = c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
             y = mrand_obs[, 1]))
points(c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
       mrand_obs[, 2], ylim = c(.75, 1), col = "red")
lines(spline(x = c(7,8,9,10,11,12,13,14,15,16,17,18,19,20),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(epi, ndims = 30)
dev.off()

epi <- RunUMAP(epi, dims = 1:11)
DimPlot(epi, group.by = "celltype")

# ======================================
# Immune
# ======================================
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 1)

if (runFeaturePlots) {FeaturePlot(immune, features = c("LYZ", "MARCO", "FCGR1A", "C1QB"))}
Idents(immune, cells = WhichCells(immune, idents = c(0, 1, 2, 3, 4, 7, 8, 9, 11, 12, 14, 19, 25, 30))) <- "Macrophages"

if (runFeaturePlots) {FeaturePlot(immune, features = c("CD3E", "FOXP3", "IL7R", "CD8A"))}
Idents(immune, cells = WhichCells(immune, idents = c(5, 15, 23, 27, 29))) <- "T Cells"

if (runFeaturePlots) {FeaturePlot(immune, features = c("CD14", "LYZ", "S100A12", "ITGAL"))}
Idents(immune, cells = WhichCells(immune, idents = c(6, 18, 20, 21))) <- "Monocytes"

if (runFeaturePlots) {FeaturePlot(immune, features = c("MS4A1", "CD19", "JCHAIN", "IGHG1"))}
Idents(immune, cells = WhichCells(immune, idents = c(22))) <- "B Cells"
Idents(immune, cells = WhichCells(immune, idents = c(26))) <- "Plasma Cells"

if (runFeaturePlots) {FeaturePlot(immune, features = c("CPA3", "KIT"))}
Idents(immune, cells = WhichCells(immune, idents = c(24))) <- "Mast Cells"

if (runFeaturePlots) {FeaturePlot(immune, features = c("FCER1A", "CD1C", "CLEC9A", "FSCN1"))}
Idents(immune, cells = WhichCells(immune, idents = c(13))) <- "cDCs"
if (runFeaturePlots) {FeaturePlot(immune, features = c("LILRA4", "CLEC4C", "JCHAIN", "TCF4"))}
Idents(immune, cells = WhichCells(immune, idents = c(31))) <- "pDCs"

if (runFeaturePlots) {FeaturePlot(immune, features = c("CD3E", "KLRB1", "NKG7", "NCR1"))}
Idents(immune, cells = WhichCells(immune, idents = c(10))) <- "NK Cells" 

if (runFeaturePlots) {FeaturePlot(immune, features = c("MKI67", "CD3E", "LYZ"))}
if (runFeaturePlots) {FeaturePlot(immune, features = c("PTPRC", "EPCAM", "VWF", "ACTA2"))}
Idents(immune, cells = WhichCells(immune, idents = c(28))) <- "Proliferating T Cells"
Idents(immune, cells = WhichCells(immune, idents = c(16))) <- "Proliferating Macrophages"
Idents(immune, cells = WhichCells(immune, idents = c(17))) <- "Doublets"

immune$celltype <- Idents(immune)
DimPlot(immune)

DimPlot(immune, group.by = "celltype")
immune_keep <- c("cDCs", "Proliferating Macrophages", "Proliferating T Cells", "NK Cells", "Plasma Cells", "moDCs", "Mast Cells", "B Cells", "Monocytes", "pDCs", "T Cells", "Macrophages")
immune <- subset(immune, cells = row.names(immune@meta.data[immune@meta.data$celltype %in% immune_keep, ]))

potato <- matrix(ncol = 15, nrow = ncol(immune))
potato2 <- matrix(ncol = 15, nrow = ncol(immune))

j <- 0
for(i in 13:27) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(immune, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, ] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 13:27
colnames(potato2) <- 13:27
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(immune), 3000),]
potato2_new <- potato2[sample(1:ncol(immune), 3000),]

for(i in  1:c(dim(potato_new)[2] - 1)){
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("immune_PCA.pdf")
plot(c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
     mrand_obs[, 1], ylim = c(0, 1))
lines(spline(x = c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
             y = mrand_obs[, 1]))
points(c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
       mrand_obs[, 2], ylim = c(0, 1), col = "red")
lines(spline(x = c(13,14,15,16,17,18,19,20,21,22,23,24,25,26),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(immune, ndims = 30)
dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

immune <- RunUMAP(immune, dims = 1:20)
DimPlot(immune, group.by = "celltype")

# ======================================
# Mesenchymal
# ======================================
meso <- FindNeighbors(meso, dims = 1:20)
meso <- FindClusters(meso, resolution = 0.5)

if (runFeaturePlots) {FeaturePlot(meso, features = c("LUM", "ACTA2", "WT1", "PDGFRA"))}
Idents(meso, cells = WhichCells(meso, idents = c(2, 3))) <- "Smooth Muscle Cells"
Idents(meso, cells = WhichCells(meso, idents = c(7))) <- "Mesothelial Cells"

if (runFeaturePlots) {FeaturePlot(meso, features = c("LUM", "HAS1", "ACTA2", "PLIN2"))}
Idents(meso, cells = WhichCells(meso, idents = c(1, 5, 6, 12))) <- "Myofibroblasts"
Idents(meso, cells = WhichCells(meso, idents = c(9))) <- "HAS1 High Fibrboblasts"
Idents(meso, cells = WhichCells(meso, idents = c(4))) <- "Fibroblasts"
Idents(meso, cells = WhichCells(meso, idents = c(0, 13))) <- "PLIN2+ Fibroblasts"

if (runFeaturePlots) {FeaturePlot(meso, features = c("PTPRC", "EPCAM", "PECAM1", "VWF"))}
Idents(meso, cells = WhichCells(meso, idents = c(8, 10, 11, 14))) <- "Doublets"

meso$celltype <- Idents(meso)
DimPlot(meso)

DimPlot(meso, group.by = "celltype")
onion <- as.character(meso@meta.data$celltype)
onion[onion == 'HAS1 High Fibrboblasts'] <- "HAS1 High Fibroblasts"
meso@meta.data$celltype <- onion
meso_keep <- c("PLIN2+ Fibroblasts", "Myofibroblasts", "HAS1 High Fibroblasts", "Mesothelial Cells", "Smooth Muscle Cells", "Fibroblasts")
meso <- subset(meso, cells = row.names(meso@meta.data[meso@meta.data$celltype %in% meso_keep, ]))

potato <- matrix(ncol = 20, nrow = ncol(meso))
potato2 <- matrix(ncol = 20, nrow = ncol(meso))

j <- 0
for(i in 6:25) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(meso, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 6:25
colnames(potato2) <- 6:25
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(meso), 3000), ]
potato2_new <- potato2[sample(1:ncol(meso), 3000), ]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("meso_PCA.pdf")
plot(c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
     mrand_obs[, 1], ylim = c(0, 1))
lines(spline(x = c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 1]))
points(c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
       mrand_obs[, 2], ylim = c(0, 1), col = "red")
lines(spline(x = c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 2]), col = "red")

dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

meso <- RunUMAP(meso, dims = 1:16)
DimPlot(meso, group.by = "celltype")

# ======================================
# Endothelial
# ======================================
endo <- FindNeighbors(endo, dims = 1:10)
endo <- FindClusters(endo)

if (runFeaturePlots) {FeaturePlot(endo, features = c("VWF", "EDNRB", "CCL21", "PECAM1"))}
if (runFeaturePlots) {FeaturePlot(endo, features = c("PTPRC", "EPCAM", "LUM", "ACTA2"))}
Idents(endo, cells = WhichCells(endo, idents = c(0:6, 8, 9, 12:15, 18))) <- "Endothelial Cells"
Idents(endo, cells = WhichCells(endo, idents = c(7, 11))) <- "Lymphatic Endothelial Cells"
Idents(endo, cells = WhichCells(endo, idents = c(10, 16, 17))) <- "Doublets"
endo$celltype <- Idents(endo)

DimPlot(endo, group.by = "celltype")
endo_keep <- c("Lymphatic Endothelial Cells", "Endothelial Cells")
endo <- subset(endo, cells = row.names(endo@meta.data[endo@meta.data$celltype %in% endo_keep, ]))

potato <- matrix(ncol = 12, nrow = ncol(endo))
potato2 <- matrix(ncol = 12, nrow = ncol(endo))

j <- 0
for(i in 9:20) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(endo, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 9:20
colnames(potato2) <- 9:20
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(endo), 3000),]
potato2_new <- potato2[sample(1:ncol(endo), 3000),]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("endo_PCA.pdf")
plot(c(9,10,11,12,13,14,15,16,17,18,19),
     mrand_obs[, 1], ylim = c(.75, 1))
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19),
             y = mrand_obs[, 1]))
points(c(9,10,11,12,13,14,15,16,17,18,19),
       mrand_obs[, 2], ylim = c(.75, 1), col = "red")
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19),
             y = mrand_obs[, 2]), col = "red")

ElbowPlot(endo)
dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

endo <- RunUMAP(endo, dims = 1:9)
DimPlot(endo, group.by = "celltype")

# ======================================
# Merge all objects back into ild
# ======================================
ild <- merge(x = epi, y=c(endo, meso, immune))
ild <- FindVariableFeatures(ild, verbose = T, nfeatures = 3000)
ild <- ScaleData(ild, features = row.names(ild@assays$SCT@data))
ild <- RunPCA(ild)

potato <-  matrix(ncol = 17, nrow = ncol(ild))
potato2 <-  matrix(ncol = 17, nrow = ncol(ild))

j <- 0
for(i in 9:25) {
  j <- j + 1
  print(i)
  onion <- RunUMAP(ild, dims = 1:i, verbose = F)
  potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
  potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
}

colnames(potato) <- 9:25
colnames(potato2) <- 9:25
mrand_obs <- NULL
potato_new <- potato[sample(1:ncol(ild), 3000), ]
potato2_new <- potato2[sample(1:ncol(ild), 3000), ]

for(i in  1:c(dim(potato_new)[2] - 1)) {
  onion <- dist(potato_new[, i])
  onion2 <- dist(potato_new[, i + 1])
  onion3 <- mantel.randtest(onion, onion2)$obs
  onion <- dist(potato2_new[, i])
  onion2 <- dist(potato2_new[, i + 1])
  onion4 <- mantel.randtest(onion, onion2)$obs
  mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
}

pdf("ild_PCA.pdf")
plot(c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
     mrand_obs[, 1], ylim = c(0.8, 1))
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 1]))
points(c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
       mrand_obs[, 2], ylim = c(0.8, 1), col = "red")
lines(spline(x = c(9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
             y = mrand_obs[, 2]), col = "red")
ElbowPlot(ild)
dev.off()
rm(onion,onion2,onion3,onion4,mrand_obs,potato,potato2, potato_new,potato2_new)

ild <- RunUMAP(ild, dims = 1:20)

# ======================================
# Apply metadata
# ======================================
meta.data <- read.csv("GSE135893_IPF_metadata.csv", header = T, row.names = 1)
ild@meta.data <- meta.data

epi_meta <- meta.data[rownames(meta.data) %in% rownames(epi@meta.data), ]
epi@meta.data <- epi_meta

endo_meta <- meta.data[rownames(meta.data) %in% rownames(endo@meta.data), ]
endo@meta.data <- endo_meta


immune_meta <- meta.data[rownames(meta.data) %in% rownames(immune@meta.data), ]
immune@meta.data <- immune_meta


meso_meta <- meta.data[rownames(meta.data) %in% rownames(meso@meta.data), ]
meso@meta.data <- meso_meta

# ======================================
# Figure 1: B
# ======================================
DimPlot(ild, group.by = "celltype")

# ======================================
# Figure 1: C
# ======================================
DimPlot(ild, group.by = "Status")
DimPlot(ild, group.by = "Diagnosis")

# ======================================
# Save final files
# ======================================
saveRDS(ild, file = "ILD.rds")
saveRDS(epi, file = "Epithelial.rds")
saveRDS(endo, file = "Endothelial.rds")
saveRDS(immune, file = "Immune.rds")
saveRDS(meso, file = "Mesochymal.rds")
