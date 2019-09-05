# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 20/08/2019
# Description: Code for Figure 2: G
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
set.seed(12345)

# ======================================
# Load libraries
# ======================================
library(Seurat, quietly = TRUE, verbose = FALSE)
library(slingshot, quietly = TRUE, verbose = FALSE)
library(scater, quietly = TRUE, verbose = FALSE)
library(RColorBrewer, quietly = TRUE, verbose = FALSE)
library(gam, quietly = TRUE, verbose = FALSE)
library(clusterExperiment, quietly = TRUE, verbose = FALSE)
library(destiny, quietly = TRUE, verbose = FALSE)
library(mclust)
library(rgl)
library(ggplot2)
library(dplyr)
library(gplots)
library(stats)

# ==============================================================================
# FIGURE 2: TRAJECTORY ANALYSIS FOR AT1, TRANS AT2, AT2 AND SCGB3A2+ (both ILD and CONTROL)
# ==============================================================================
# Read in the objects 
epi <- readRDS("Epithelial.rds")

# KRT5-/KRT17+, Transitional AT2, AT2, AT1 and Club cells
krt5_6pop <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$celltype %in%
                                                                  c("KRT5-/KRT17+",
                                                                    "Transitional AT2",
                                                                    "AT2",
                                                                    "AT1",
                                                                    "SCGB3A2+",
                                                                    "SCGB3A2+ SCGB1A1+"),]))
krt5_6pop <- FindVariableFeatures(krt5_6pop, verbose = F, nfeatures = 3000)
krt5_6pop <- ScaleData(krt5_6pop, verbose = F)
krt5_6pop <- RunUMAP(krt5_6pop, dims = 1:10, verbose = F)

sub1 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("AT1","Transitional AT2", "SCGB3A2+"),]))
sub2 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("AT1","Transitional AT2", "AT2"),]))
sub1 <- FindVariableFeatures(sub1, nfeatures = 3000)
sub1 <- ScaleData(sub1)
sub1 <- RunPCA(sub1)
sub1 <- RunUMAP(sub1, dims = 1:10)

sub2 <- FindVariableFeatures(sub2, nfeatures = 3000)
sub2 <- ScaleData(sub2, features = row.names(sub2@assays$SCT@data), verbose = F)
sub2 <- RunPCA(sub2)
sub2 <- RunUMAP(sub2, dims = 1:6)

sub1@meta.data$celltype <- factor(sub1@meta.data$celltype, levels = c("SCGB3A2+","Transitional AT2","AT1"))
sub2@meta.data$celltype <- factor(sub2@meta.data$celltype, levels = c("AT2","Transitional AT2","AT1"))

# Convert to SingleCellExperiment
sub_sce1 <- as.SingleCellExperiment(sub1)
sub_sce2 <- as.SingleCellExperiment(sub2)

# Dimentionality reduction
onion <- sub_sce1@reducedDims@listData$PCA[,1:20]
sub_sce1@reducedDims@listData$PCA <- onion

onion <- sub_sce2@reducedDims@listData$PCA[,1:20]
sub_sce2@reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce1), ignore.case = F, value = T) 
sub_sce1 <- sub_sce1[!rownames(sub_sce1) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce1), ignore.case = F, value = T) 
sub_sce1 <- sub_sce1[!rownames(sub_sce1) %in% temp2,] 

temp <- grep( "^MT-", rownames(sub_sce2), ignore.case = F, value = T) 
sub_sce2 <- sub_sce2[!rownames(sub_sce2) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce2), ignore.case = F, value = T) 
sub_sce2 <- sub_sce2[!rownames(sub_sce2) %in% temp2,] 

# Run Slingshot
sub_slingshot1 <- slingshot(sub_sce1, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="SCGB3A2+", end.clus="AT1")
print(SlingshotDataSet(sub_slingshot1))

sub_slingshot2 <- slingshot(sub_sce2, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="AT1")
print(SlingshotDataSet(sub_slingshot2))

# Plot the trajectory
epi_color <- as.character(sub2@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(sub2@meta.data))

onion <- epi_color
onion[onion == "SCGB3A2+"] <- "#F659DD"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "AT2"] <- "#EE7342"
onion[onion == "AT1"] <- "#F76A62"
epi_color <- onion

# Set the pseudotime variable
t1 <- sub_slingshot1$slingPseudotime_1 
t2 <- sub_slingshot2$slingPseudotime_1 

gene.list1 <- c("AGER","ABCA3","SCGB3A2","SFTPC")
# Prepare data for loess plot for marker genes - SCGB3A2, Trans AT2, AT1
loess_data1 = as.data.frame(sub1@assays$SCT@data[gene.list1,])
loess_data1 = loess_data1[,order(t1)]
temp1 <- loess_data1
temp1 <- t(temp1)
temp1 = as.data.frame(temp1)
temp1$index = 1:nrow(temp1)
temp1$ct = sub1@meta.data$celltype[order(t1)]

# Prepare data for loess plot for marker genes - AT2, Trans AT2, AT1
loess_data2 = as.data.frame(sub2@assays$SCT@data[gene.list1,])
loess_data2 = loess_data2[,order(t2)]
temp2 <- loess_data2
temp2 <- t(temp2)
temp2 = as.data.frame(temp2)
temp2$index = 1:nrow(temp2)
temp2$ct = sub2@meta.data$celltype[order(t2)]

# ======================================
# Figure 2: G
# ======================================
pdf(file = paste(date, "Figure_2_Loess_plots.pdf", sep = "_"))
p1 <- ggplot(temp1, aes(y = AGER, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .2, fill=ct)) + guides(fill=guide_legend())
p2 <- ggplot(temp1, aes(y = ABCA3, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 2)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
p3 <- ggplot(temp1, aes(y = SFTPC, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 7)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .3, fill=ct)) + guides(fill=guide_legend())
p4 <- ggplot(temp1, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4.5)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .2, fill=ct)) + guides(fill=guide_legend())

p5 <- ggplot(temp2, aes(y = AGER, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .2, fill =ct)) + guides(fill=guide_legend())
p6 <- ggplot(temp2, aes(y = ABCA3, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 2)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
p7 <- ggplot(temp2, aes(y = SFTPC, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 7)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .3, fill=ct)) + guides(fill=guide_legend())
p8 <- ggplot(temp2, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4.5)) +
  geom_tile(aes(x = index, y= 0, color = ct, height = .2, fill=ct)) + guides(fill=guide_legend())
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
dev.off()

