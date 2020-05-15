  # ==============================================================================
  # Author(s) : Linh T. Bui, lbui@tgen.org
  # Date: 20/08/2019
  # Description: Code for Figure 3: B
  # ==============================================================================

# Environment parameters
set.seed(12345)

# Load libraries
library(Seurat, quietly = TRUE, verbose = FALSE)
library(slingshot, quietly = TRUE, verbose = FALSE)
library(gridExtra)
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

# Read in the objects 
sub <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200513/4pop_ctrl_and_disease.rds")
sub1 <- subset(sub, cells=rownames(sub@meta.data[sub@meta.data$celltype %in% c("AT2", "Transitional AT2","AT1"),]))
sub2 <- subset(sub, cells=rownames(sub@meta.data[sub@meta.data$celltype %in% c("SCGB3A2+", "Transitional AT2","AT1"),]))

sub_sce1 <- as.SingleCellExperiment(sub1)
sub_sce2 <- as.SingleCellExperiment(sub2)

# Dimentionality reduction
onion <- sub_sce1@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce1@int_colData$reducedDims@listData$PCA <- onion

onion2 <- sub_sce2@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce2@int_colData$reducedDims@listData$PCA <- onion2

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce1), ignore.case = F, value = T) 
sub_sce1 <- sub_sce1[!rownames(sub_sce1) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce1), ignore.case = F, value = T) 
sub_sce1 <- sub_sce1[!rownames(sub_sce1) %in% temp2,] 

temp3 <- grep( "^MT-", rownames(sub_sce2), ignore.case = F, value = T) 
sub_sce2 <- sub_sce2[!rownames(sub_sce2) %in% temp3,]
temp4 <- grep( "^RP", rownames(sub_sce2), ignore.case = F, value = T) 
sub_sce2 <- sub_sce2[!rownames(sub_sce2) %in% temp4,] 

# Run Slingshot
sub_slingshot1 <- slingshot(sub_sce1, clusterLabels = "celltype", reducedDim = 'UMAP',
                           start.clus="AT2", end.clus="AT1")
summary(sub_slingshot1$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot1))

sub_slingshot2 <- slingshot(sub_sce2, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="SCGB3A2+", end.clus="AT1")
summary(sub_slingshot2$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot2))

# Arrange levels
sub1@meta.data$celltype <- factor(sub1@meta.data$celltype, levels = c("AT2","Transitional AT2","AT1"))
sub2@meta.data$celltype <- factor(sub2@meta.data$celltype, levels = c("SCGB3A2+","Transitional AT2","AT1"))

# Set the pseudotime variable
print(SlingshotDataSet(sub_slingshot1))
print(SlingshotDataSet(sub_slingshot2))

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

# Make the Loess plot
p1 <- ggplot(temp1, aes(y = AGER, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 4.2)) 
#  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) 
p2 <- ggplot(temp1, aes(y = ABCA3, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 1.6)) 
#  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
p3 <- ggplot(temp1, aes(y = SFTPC, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 7)) 
#  geom_tile(aes(x = index, y= 0, color = ct, height = .3, fill=ct)) + guides(fill=guide_legend())
p4 <- ggplot(temp1, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 6)) 
 # geom_tile(aes(x = index, y= 0, color = ct, height = .2, fill=ct)) + guides(fill=guide_legend())

p5 <- ggplot(temp2, aes(y = AGER, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 4.2)) 
#  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill =ct)) + guides(fill=guide_legend())
p6 <- ggplot(temp2, aes(y = ABCA3, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 1.6)) 
# geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
p7 <- ggplot(temp2, aes(y = SFTPC, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, 7)) 
# geom_tile(aes(x = index, y= 0, color = ct, height = .3, fill=ct)) + guides(fill=guide_legend())
p8 <- ggplot(temp2, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0,6)) 
# geom_tile(aes(x = index, y= 0, color = ct, height = .2, fill=ct)) + guides(fill=guide_legend())
  
pdf(file = "Fig3B_loess_final.pdf", width = 11, height = 8.5)
p1 + p5
p2 + p6
p3 + p7
p4 + p8
dev.off()

