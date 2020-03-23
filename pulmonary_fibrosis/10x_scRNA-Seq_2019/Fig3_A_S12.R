# ==============================================================================
# Author(s) : Austin J. Gutierrez, agutierrez@tgen.org
# Date: 13/03/2020
# Description: Figure 3: A
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
setwd("/scratch/agutierrez/IPF/R/Seurat")
main_dir <- getwd()
current_date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, current_date), showWarnings = FALSE)
setwd(file.path(main_dir, current_date))

options(future.globals.maxSize = 2048*1024^2 )
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

withRNA_epi <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$population == "Epithelial", ])) 

epi_col <- c("AT1"="#548BC5","AT2"="#EE7342","Basal"="#D38402","Ciliated"="#B19302",
             "Differentiating Ciliated" = "#9C9966","KRT5-/KRT17+"="#03AF21",
             "MUC5AC+ High" = "#003366","MUC5B+"="#00B2DB","Proliferating Epithelial Cells"=
               "#B874FF", "Transitional AT2"="#F1A5C4", "SCGB3A2+"="#F659DD",
             "SCGB3A2+ SCGB1A1+"="#9900CC")

# ======================================
# 
# ======================================
# Subset out only AT2 samples from Epi population
fig3 <- subset(withRNA_epi, cells = rownames(withRNA_epi@meta.data[withRNA_epi@meta.data$celltype %in%
                                                                     c("Transitional AT2", "AT2", "AT1", "SCGB3A2+"),]))
dim(fig3@meta.data)

fig3 <- SCTransform(fig3, batch_var = "Flowcell_ID")
fig3 <- RunPCA(fig3)
fig3 <- RunUMAP(fig3, dims = 1:5)
DimPlot(fig3, cols = epi_col)


# ======================================
# Convert to SingleCellExperiment
# ======================================
sub_sce <- as.SingleCellExperiment(fig3)

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
                           start.clus="AT2", end.clus="AT1")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
epi_color <- as.character(fig3@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(fig3@meta.data))

onion <- epi_color
onion[onion == "SCGB3A2+"] <- "#F659DD"
onion[onion == "Transitional AT2"] <- "#F1A5C4"
onion[onion == "AT2"] <- "#EE7342"
onion[onion == "AT1"] <- "#548BC5"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "AT1", "SCGB3A2+"), 
       col = c("#EE7342","#F1A5C4","#548BC5", "#F659DD"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

# ======================================
# Save lineage
# ======================================
umap_lineage <- sub_slingshot$slingPseudotime_2
names(umap_lineage) <- colnames(sub_slingshot)

# ======================================
# PCA
# ======================================
sub_sce <- as.SingleCellExperiment(fig3)

# Dimentionality reduction
onion <- sub_sce@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce@int_colData$reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) #13 MT genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) #4723 RP genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] # total 21641 genes, 6406 cells

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'PCA',
                           start.clus="AT2", end.clus="AT1")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
plot(reducedDims(sub_slingshot)$PCA, col = epi_color,
     pch = 16, asp = 1, xlab = "PC_1", ylab = "PC_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "AT1", "SCGB3A2+"), 
       col = c("#EE7342","#F1A5C4","#548BC5", "#F659DD"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

# ======================================
# Save lineage
# ======================================
pca_lineage <- sub_slingshot$slingPseudotime_2
names(pca_lineage) <- colnames(sub_slingshot)

# ======================================
# Save PC_1
# ======================================
pc_1 <- fig3@reductions$pca[[,1]]

# ======================================
# Slap everything together
# ======================================
pca_lineage <- as.data.frame(pca_lineage)
umap_lineage <- as.data.frame(umap_lineage)
pc_1 <- as.data.frame(pc_1[,1])

onion <- cbind(pc_1, pca_lineage, umap_lineage, cell_color)
colnames(onion) <- c("pc_1", "pca_lineage", "umap_lineage", "cell_color")

# Remove KRT5
onion <- onion[onion$cell_color != "#F659DD",]

# ==============================================================================
# 
# ==============================================================================
cor <- cor.test(onion$pca_lineage, onion$pc_1, method = "spearman")
cor <- round(as.numeric(cor$estimate), 3)
plot(as.numeric(onion$pca_lineage),
     as.numeric(onion$pc_1),
     col = onion$cell_color,
     pch = 16,
     xlab = "PCA_lineage",
     ylab = "PC_1",
     main = "AT2_to_AT1",
     sub = paste("spearman", cor, sep = "_"))

cor <- cor.test(onion$umap_lineage, onion$pc_1, method = "spearman")
cor <- round(as.numeric(cor$estimate), 3)
plot(as.numeric(onion$umap_lineage),
     as.numeric(onion$pc_1),
     col = onion$cell_color,
     pch = 16,
     xlab = "UMAP_lineage",
     ylab = "PC_1",
     main = "AT2_to_AT1",
     sub = paste("spearman", cor, sep = "_"))

cor <- cor.test(onion$umap_lineage, onion$pca_lineage, method = "spearman")
cor <- round(as.numeric(cor$estimate), 3)
plot(as.numeric(onion$umap_lineage),
     as.numeric(onion$pca_lineage),
     col = onion$cell_color,
     pch = 16,
     xlab = "UMAP_lineage",
     ylab = "PCA_lineage",
     main = "AT2_to_AT1",
     sub = paste("spearman", cor, sep = "_"))


# ==============================================================================
# ==============================================================================

# ======================================
# Convert to SingleCellExperiment
# ======================================
sub_sce <- as.SingleCellExperiment(fig3)

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
                           start.clus="SCGB3A2+", end.clus="AT1")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
epi_color <- as.character(fig3@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(fig3@meta.data))

onion <- epi_color
onion[onion == "SCGB3A2+"] <- "#F659DD"
onion[onion == "Transitional AT2"] <- "#F1A5C4"
onion[onion == "AT2"] <- "#EE7342"
onion[onion == "AT1"] <- "#548BC5"
epi_color <- onion

cell_color <- as.data.frame(fig3@meta.data)
cell_color <-  cell_color[ -c(1:12, 14:20) ]
cell_color[cell_color == "Transitional AT2"] <- "#F1A5C4"
cell_color[cell_color == "AT2"] <- "#EE7342"
cell_color[cell_color == "AT1"] <- "#548BC5"
cell_color[cell_color == "SCGB3A2+"] <- "#F659DD"


plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "AT1", "SCGB3A2+"), 
       col = c("#EE7342","#F1A5C4","#548BC5", "#F659DD"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

# ======================================
# Save lineage
# ======================================
umap_lineage <- sub_slingshot$slingPseudotime_2
names(umap_lineage) <- colnames(sub_slingshot)

# ======================================
# PCA
# ======================================
sub_sce <- as.SingleCellExperiment(fig3)

# Dimentionality reduction
onion <- sub_sce@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce@int_colData$reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) #13 MT genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) #4723 RP genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] # total 21641 genes, 6406 cells

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'PCA',
                           start.clus="SCGB3A2+", end.clus="AT1")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
plot(reducedDims(sub_slingshot)$PCA, col = epi_color,
     pch = 16, asp = 1, xlab = "PC_1", ylab = "PC_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "AT1", "SCGB3A2+"), 
       col = c("#EE7342","#F1A5C4","#548BC5", "#F659DD"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

# ======================================
# Save lineage
# ======================================
pca_lineage <- sub_slingshot$slingPseudotime_2
names(pca_lineage) <- colnames(sub_slingshot)

# ======================================
# Save PC_1
# ======================================
pc_1 <- fig3@reductions$pca[[,2]]

# ======================================
# Slap everything together
# ======================================
pca_lineage <- as.data.frame(pca_lineage)
umap_lineage <- as.data.frame(umap_lineage)
pc_1 <- as.data.frame(pc_1[,1])

onion <- cbind(pc_1, pca_lineage, umap_lineage, cell_color)
colnames(onion) <- c("pc_1", "pca_lineage", "umap_lineage", "cell_color")

# Remove AT2
onion <- onion[onion$cell_color != "#EE7342",]

# ==============================================================================
# 
# ==============================================================================
cor <- cor.test(onion$pca_lineage, onion$pc_1, method = "spearman")
cor <- round(as.numeric(cor$estimate), 3)
plot(as.numeric(onion$pca_lineage),
     as.numeric(onion$pc_1),
     col = onion$cell_color,
     pch = 16,
     xlab = "PCA_lineage",
     ylab = "PC_2",
     main = "SCGB3A2+_to_AT1",
     sub = paste("spearman", cor, sep = "_"))

cor <- cor.test(onion$umap_lineage, onion$pc_1, method = "spearman")
cor <- round(as.numeric(cor$estimate), 3)
plot(as.numeric(onion$umap_lineage),
     as.numeric(onion$pc_1),
     col = onion$cell_color,
     pch = 16,
     xlab = "UMAP_lineage",
     ylab = "PC_2",
     main = "SCGB3A2+_to_AT1",
     sub = paste("spearman", cor, sep = "_"))

cor <- cor.test(onion$umap_lineage, onion$pca_lineage, method = "spearman")
cor <- round(as.numeric(cor$estimate), 3)
plot(as.numeric(onion$umap_lineage),
     as.numeric(onion$pca_lineage),
     col = onion$cell_color,
     pch = 16,
     xlab = "UMAP_lineage",
     ylab = "PCA_lineage",
     main = "SCGB3A2+_to_AT1",
     sub = paste("spearman", cor, sep = "_"))
