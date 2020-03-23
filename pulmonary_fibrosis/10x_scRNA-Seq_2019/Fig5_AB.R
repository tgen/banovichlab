# ======================================
# Environment parameters
# ======================================
setwd("/scratch/agutierrez/IPF/R/Seurat")
main_dir <- getwd()
current_date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, current_date), showWarnings = FALSE)
setwd(file.path(main_dir, current_date))

set.seed(2811)
# ======================================
# Load libraries
# ======================================
library(Seurat)
library(slingshot)

# ======================================
# Read in data
# ======================================
ild <- readRDS("/scratch/agutierrez/IPF/R/Seurat/Reference/2019_Release_IPF.rds")
ild <- ApplyMetaData(ild)
epi <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Epithelial", ]))
rm(ild)

# ==============================================================================
# 
# ==============================================================================
# Subset out only ILD samples from Epi population
epi_ild <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$Status == "Disease",]))
rm(epi)

# five celltypes
five <- subset(epi_ild, cells = rownames(epi_ild@meta.data[epi_ild@meta.data$celltype %in%
                                                            c("KRT5-/KRT17+","Transitional AT2", "SCGB3A2+", "AT1", "AT2"),]))

five <- SCTransform(five, batch_var = "Flowcell_ID")
five <- RunPCA(five)
ElbowPlot(five)
five <- RunUMAP(five, dims = 1:10, verbose = F)
DimPlot(five, group.by = "celltype", cols = epi_col)

# ==============================================================================
# 
# ==============================================================================
SCGB3A2 <- subset(five, cells = rownames(five@meta.data[five@meta.data$celltype %in%
                                                             c("KRT5-/KRT17+","Transitional AT2", "SCGB3A2+", "AT1"),]))

SCGB3A2 <- SCTransform(SCGB3A2, batch_var = "Flowcell_ID")
#SCGB3A2 <- RunPCA(SCGB3A2)
ElbowPlot(SCGB3A2)
SCGB3A2 <- RunUMAP(SCGB3A2, dims = 1:10, verbose = F)
DimPlot(SCGB3A2, group.by = "celltype", cols = epi_col)

# ==============================================================================
# 
# ==============================================================================
AT2 <- subset(five, cells = rownames(five@meta.data[five@meta.data$celltype %in%
                                                          c("KRT5-/KRT17+","Transitional AT2", "AT2", "AT1"),]))

AT2 <- SCTransform(AT2, batch_var = "Flowcell_ID")
#AT2 <- RunPCA(AT2)
ElbowPlot(AT2)
AT2 <- RunUMAP(AT2, dims = 1:10, verbose = F)
DimPlot(AT2, group.by = "celltype", cols = epi_col)

# ==============================================================================
# 
# ==============================================================================
# ======================================
# Convert to SingleCellExperiment
# ======================================
sub_sce <- as.SingleCellExperiment(SCGB3A2)

# Dimentionality reduction
onion <- sub_sce@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce@int_colData$reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) #13 MT genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) #4723 RP genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] # total 21641 genes, 6406 cells

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'UMAP',
                           start.clus="SCGB3A2+", end.clus=c("KRT5-/KRT17+", "AT1"))
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

epi_color <- as.character(SCGB3A2@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(SCGB3A2@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#F1A5C4"
onion[onion == "SCGB3A2+"] <- "#F659DD"
onion[onion == "AT1"] <- "#548BC5"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("SCGB3A2+","Transitional AT2", "KRT5-/KRT17+", "AT1"), 
       col = c("#F659DD","#F1A5C4","#03AF21","#548BC5"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

# ==============================================================================
# 
# ==============================================================================
# ======================================
# Convert to SingleCellExperiment
# ======================================
sub_sce <- as.SingleCellExperiment(AT2)

# Dimentionality reduction
onion <- sub_sce@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce@int_colData$reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) #13 MT genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) #4723 RP genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] # total 21641 genes, 6406 cells

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'UMAP',
                           start.clus="AT2", end.clus=c("KRT5-/KRT17+", "AT1"))
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

epi_color <- as.character(AT2@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(AT2@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#F1A5C4"
onion[onion == "AT2"] <- "#EE7342"
onion[onion == "AT1"] <- "#548BC5"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "KRT5-/KRT17+", "AT1"), 
       col = c("#EE7342","#F1A5C4","#03AF21","#548BC5"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)
lines(slingCurves(sub_slingshot)$curve2, lwd=3)
