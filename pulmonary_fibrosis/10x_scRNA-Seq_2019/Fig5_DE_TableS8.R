# ================================
# Author: Linh Bui, lbui@tgen.org
# Date: 2020-03-11
# Science Advances - Revision Figure 5D, 5E and Table S8
# ================================

library(Seurat, quietly = TRUE, verbose = FALSE)
library(RColorBrewer, quietly = TRUE, verbose = FALSE)
library(gam, quietly = TRUE, verbose = FALSE)
library(destiny, quietly = TRUE, verbose = FALSE)
library(mclust)
library(ggplot2)
library(dplyr)
library(slingshot)
library(heatmap3)
library(Matrix)
library(clusterExperiment)
library(rgl)
library(gplots)
library(qvalue)

set.seed(12345)
# Set date and create out folder
getwd()
Sys.Date()
main_dir <- "~/Desktop/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# ==============================================================================
# FIGURE 5D: AT2 - Transitional AT2 - KRT5-/KRT17+
# ==============================================================================

# Read in the Seurat object and subset out AT1 cells
at2 <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200310/Final_AT2.rds")
krt5 <- subset(at2, cells=rownames(at2@meta.data[at2@meta.data$celltype %in% c("AT2","Transitional AT2","KRT5-/KRT17+"),]))

# Convert to SingleCellExperiment
sub_sce <- as.SingleCellExperiment(krt5)

# Dimentionality reduction
onion <- sub_sce@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce@int_colData$reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) 
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) 
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] 

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'UMAP',
                           start.clus="AT2", end.clus="KRT5-/KRT17+")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
epi_color <- as.character(krt5@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(krt5@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "AT2"] <- "#EE7342"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "KRT5-/KRT17+"), 
       col = c("#EE7342","#FE627D","#03AF21"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)

# Set the pseudotime variable
t <- sub_slingshot$slingPseudotime_1 #  lineage 1: KRT5-/KRT17

# Extract the gene expression matrix
Y <- assay(sub_slingshot)

# Fit a GAM with a loess term for pseudotime 
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = "Umap_AT2_gampval_allgenes_noAT1.rds")

# Calculate qvalue of gam.pval
gam.qval <- qvalue(gam.pval)

# Select the top 400 genes that change over pseudotime
topgenes <- names(sort(gam.qval$qvalues, decreasing = FALSE))[1:400] 
heatdata <- assay(sub_slingshot)[rownames(assay(sub_slingshot)) %in% topgenes, 
                                 order(t, na.last = NA)]
allgenes <- factor(rownames(Y))
topgenes2 <- sort(gam.qval$qvalues, decreasing = FALSE)[1:400] 
write.csv(topgenes2, file = "20200311_AT2_KRT5_topgenes_qval.csv")

# Scale the data
heatdata <- t(scale(t(heatdata)))

# Trim z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -3] = -3

# Order genes based on average expression in 25 cells/cell type
mean_matrix = matrix(nrow = 400, ncol = 245)
j=0
potato = for(i in seq(1,6107, 25)){
  j=j+1
  if(i != 6101){
    k=i+25
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
  else{
    k=i+6
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
}
colnames(mean_matrix) = 1:245

mean_matrix = as.data.frame(mean_matrix)

gene_order = as.numeric(apply(mean_matrix, 1, function(xx){names(xx)[xx == max(xx)][1] }))
heatdata_onion = heatdata[order(gene_order),]
color.palette  <- colorRampPalette(c("blue","black","yellow"), space = "Lab")(n=100)

# Make heatmap
# Using heatmap.2 from gplots
# I used heatmap.2 to extract out rowInd value, I tried with plotHeatmap and the index in rownames(ce) is not same index as what used in the heatmap
pdf("20200311_AT2_KRT5_heatmap2.pdf")
h <- heatmap.2(heatdata_onion,
               trace = "none",
               col = color.palette, 
               cexRow = 0.2,
               symbreaks = F,
               scale = "none",
               labCol=F,
               main = "KRT5-/KRT17, AT2, Transitional AT2-UMAP", 
               dendrogram = "row",
               Colv = F) 
dev.off()

# Using plotHeatmap from clusterExperiment (with pseudotime bar)
heatclus <- sub_slingshot$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata_onion, heatclus, transformation = function(x){x})

pdf("20200311_AT2_KRT5_CE2.pdf")  
plotHeatmap(ce, 
            clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', 
            fontsize=15,
            annLegend=T,
            colorScale = color.palette) 
dev.off()

# Heatmap for publication with group break (I manually selected where the breaks are)
# Using heatmap.2 from gplots
pdf("20200311_AT2_KRT5_heatmap2_split.pdf")
heatmap.2(heatdata_onion,
               trace = "none",
               col = color.palette, 
               cexRow = 0.2,
               symbreaks = F,
               scale = "none",
               labCol=F,
               main = "KRT5-/KRT17, AT2, Transitional AT2-UMAP", 
               dendrogram = "row",
               Colv = F,
               rowsep = c(54,112,240)) 
dev.off()

# Using plotHeatmap from ClusterExperiment
allgenes2 <- rev(h$rowInd) # splitting genes in group based on gene index
genes1 <- allgenes2[1:63] 
genes2 <- allgenes2[64:117]
genes3 <- allgenes2[118:153]
genes4 <- allgenes2[154:400]
clusters<- list(genes1,genes2,genes3,genes4)

pdf("20200311_AT2_KRT5_CE_split.pdf")  
plotHeatmap(ce, 
            clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', 
            fontsize=15,
            annLegend=T,
            colorScale = color.palette,
            clusterFeatures = T,
            clusterFeaturesData = clusters,
            nBlankLines = 2) 
dev.off()

# Extract out geneIDs for each bin
allgenes <- rownames(heatdata_onion)[rev(h$rowInd)]
transition1_genes <- allgenes[1:63] 
transition2_genes <- allgenes[64:117]
transition3_genes <- allgenes[118:153]
transition4_genes <- allgenes[154:400]
clustergenes <- list(transition1_genes,transition2_genes,transition3_genes,
                     transition4_genes)
n.obs <- sapply(clustergenes, length) # check length of each vector in the list
seq.max <- seq_len(max(n.obs))
genes <- sapply(clustergenes, "[", i = seq.max) # add NAs into empty slots to create equal length

write.table(genes, 
            file = "AT2_KRT5_Umapheatmap_genes.csv",
            quote = F, sep = ",", row.names = F)

# ==============================================================================
# FIGURE 5E: SCGB3A2+ - Transitional AT2 - KRT5-/KRT17+
# ==============================================================================

# Read in the object and subset out AT1 population
scgb3a2 <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200228/New_SCGB3A2.rds")
krt5 <- subset(scgb3a2, 
               cells=rownames(scgb3a2@meta.data[scgb3a2@meta.data$celltype %in% c("SCGB3A2+","Transitional AT2","KRT5-/KRT17+"),]))

# Convert to SingleCellExperiment
sub_sce <- as.SingleCellExperiment(krt5)

# Dimentionality reduction
onion <- sub_sce@int_colData$reducedDims@listData$PCA[,1:20]
sub_sce@int_colData$reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) #13 MT genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) #2118 RP genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] # total 17183 genes, 4312 cells

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'UMAP',
                           start.clus="SCGB3A2+", end.clus="KRT5-/KRT17+")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
epi_color <- as.character(krt5@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(krt5@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "SCGB3A2+"] <- "#F659DD"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("SCGB3A2+","Transitional AT2", "KRT5-/KRT17+"), 
       col = c("#F659DD","#FE627D","#03AF21"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)

# Set the pseudotime variable
t <- sub_slingshot$slingPseudotime_1 #  lineage 1: KRT5-/KRT17

# Extract the gene expression matrix
Y <- assay(sub_slingshot)

# Fit a GAM with a loess term for pseudotime 
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = "Umap_SCGB3A2_gampval_allgenes.rds")

# Calculate qvalue of gam.pval
gam.qval <- qvalue(gam.pval)

# Select the top 400 genes that change over pseudotime
topgenes <- names(sort(gam.qval$qvalues, decreasing = FALSE))[1:400] 
heatdata <- assay(sub_slingshot)[rownames(assay(sub_slingshot)) %in% topgenes, 
                                 order(t, na.last = NA)]
allgenes <- factor(rownames(Y))
topgenes2 <- sort(gam.qval$qvalues, decreasing = FALSE)[1:400] 
write.csv(topgenes2, file = "20200311_SCGB3A2_KRT5_topgenes_qval.csv")

# Scale the data
heatdata <- t(scale(t(heatdata)))

# Trim z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -3] = -3

# Order genes based on average expression in 25 cells/cell type
mean_matrix = matrix(nrow = 400, ncol = 173)
j=0
potato = for(i in seq(1,4312, 25)){
  j=j+1
  if(i != 4301){
    k=i+25
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
  else{
    k=i+11
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
}
colnames(mean_matrix) = 1:173

mean_matrix = as.data.frame(mean_matrix)

gene_order = as.numeric(apply(mean_matrix, 1, function(xx){names(xx)[xx == max(xx)][1] }))
heatdata_onion = heatdata[order(gene_order),]
color.palette  <- colorRampPalette(c("blue","black","yellow"), space = "Lab")(n=100)

# Save the heatmap object and divide the tree for different transition stages
# Using heatmap.2 from gplots package
pdf("20200311_SCGB3A2_KRT5_UMAP.pdf")
h <- heatmap.2(heatdata_onion,
               trace = "none",
               col = color.palette, 
               cexRow = 0.2,
               symbreaks = F,
               scale = "none",
               labCol=F,
               main = "SCGB3A2, Transitional AT2, KRT5-/KRT17, -UMAP", 
               dendrogram = "row",
               Colv = F) 
dev.off()

# Using plotHeatmap from ClusterExperiment (this one has the pseudotime bar on top)
heatclus <- sub_slingshot$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata_onion, heatclus, transformation = function(x){x})

pdf("20200311_SCGB3A2_KRT5_CE2.pdf")  
p <- plotHeatmap(ce, 
            clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', 
            fontsize=15,
            annLegend=T,
            colorScale = color.palette)
dev.off()

# Extract out geneID for each bin
allgenes <- rownames(heatdata_onion)[rev(h$rowInd)]
transition1_genes <- allgenes[1:49] 
transition2_genes <- allgenes[50:142]
transition3_genes <- allgenes[142:283]
transition4_genes <- allgenes[284:400]
clustergenes <- list(transition1_genes,transition2_genes,transition3_genes,
                     transition4_genes)
n.obs <- sapply(clustergenes, length) # check length of each vector in the list
seq.max <- seq_len(max(n.obs))
at2_heatgenes <- sapply(clustergenes, "[", i = seq.max) # add NAs into empty slots to create equal length

write.table(at2_heatgenes, 
            file = "SCGB3A2_KRT5_Umapheatmap_genes.csv",
            quote = F, sep = ",", row.names = F)

# Make heatmap with 4 different bins for publication
# Using heatmap.2 in gplots
pdf("20200311_SCGB3A2_KRT5_heatmap2.split.pdf")
heatmap.2(heatdata_onion,
          trace = "none",
          col = color.palette, 
          cexRow = 0.2,
          symbreaks = F,
          scale = "none",
          labCol=F,
          main = "SCGB3A2, Transitional AT2, KRT5-/KRT17, -UMAP", 
          dendrogram = "row",
          Colv = F,
          rowsep = c(48,138,274)) 
dev.off()

# Using plotHeatmap in ClusterExperiment (this one has the pseudotime bar on top)
allgenes2 <- rev(h$rowInd)
genes1 <- allgenes2[1:49] 
genes2 <- allgenes2[50:142]
genes3 <- allgenes2[143:283]
genes4 <- allgenes2[284:400]
clusters<- list(genes1,genes2,genes3,genes4)

pdf("20200311_SCGB3A2_KRT5_CE_split.pdf")  
plotHeatmap(ce, 
            clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', 
            fontsize=15,
            annLegend=T,
            colorScale = color.palette,
            clusterFeatures = T,
            clusterFeaturesData = clusters,
            nBlankLines = 2) 
dev.off()

# ==============================================================================
# TABLE S8: Trajectory genes (with q-values) in the heatmap 
# ==============================================================================

# AT2 - KRT5-/KRT17+ trajectory
at2_heatgenes <- read.csv("~/Desktop/RStudio_folder/20200311/AT2_KRT5_Umapheatmap_genes.csv", header = T)
at2_qval <- read.csv("20200311_AT2_KRT5_topgenes_qval.csv", header = T)

colnames(at2_heatgenes) <- c("bin1","bin2","bin3","bin4")
colnames(at2_qval) <- c("GeneID", "qvalues")
at2_heatgenes$bin1.qval <- plyr::mapvalues(x = at2_heatgenes$bin1,
                                           from = at2_qval$GeneID,
                                           to = as.character(at2_qval$qvalues))
at2_heatgenes$bin2.qval <- plyr::mapvalues(x = at2_heatgenes$bin2,
                                           from = at2_qval$GeneID,
                                           to = as.character(at2_qval$qvalues))
at2_heatgenes$bin3.qval <- plyr::mapvalues(x = at2_heatgenes$bin3,
                                           from = at2_qval$GeneID,
                                           to = as.character(at2_qval$qvalues))
at2_heatgenes$bin4.qval <- plyr::mapvalues(x = at2_heatgenes$bin4,
                                           from = at2_qval$GeneID,
                                           to = as.character(at2_qval$qvalues))
at2_heatgenes <- at2_heatgenes[,c("bin1","bin1.qval","bin2","bin2.qval",
                                  "bin3","bin3.qval","bin4","bin4.qval")]
write.csv(at2_heatgenes, file = "20200313_AT2_KRT5_heatmapgenes.csv")

# SCGB3A2 - KRT5-/KRT17+ trajectory
scgb3a2_heatgenes <- read.csv("~/Desktop/RStudio_folder/20200311/SCGB3A2_KRT5_Umapheatmap_genes.csv", header = T)
scgb3a2_qval <- read.csv("20200311_SCGB3A2_KRT5_topgenes_qval.csv", header = T)

colnames(scgb3a2_heatgenes) <- c("bin1","bin2","bin3","bin4")
colnames(scgb3a2_qval) <- c("GeneID", "qvalues")
scgb3a2_heatgenes$bin1.qval <- plyr::mapvalues(x = scgb3a2_heatgenes$bin1,
                                           from = scgb3a2_qval$GeneID,
                                           to = as.character(scgb3a2_qval$qvalues))
scgb3a2_heatgenes$bin2.qval <- plyr::mapvalues(x = scgb3a2_heatgenes$bin2,
                                           from = scgb3a2_qval$GeneID,
                                           to = as.character(scgb3a2_qval$qvalues))
scgb3a2_heatgenes$bin3.qval <- plyr::mapvalues(x = scgb3a2_heatgenes$bin3,
                                           from = scgb3a2_qval$GeneID,
                                           to = as.character(scgb3a2_qval$qvalues))
scgb3a2_heatgenes$bin4.qval <- plyr::mapvalues(x = scgb3a2_heatgenes$bin4,
                                           from = scgb3a2_qval$GeneID,
                                           to = as.character(scgb3a2_qval$qvalues))
scgb3a2_heatgenes <- scgb3a2_heatgenes[,c("bin1","bin1.qval","bin2","bin2.qval",
                                  "bin3","bin3.qval","bin4","bin4.qval")]
write.csv(scgb3a2_heatgenes, file = "20200313_SCGB3A2_KRT5_heatmapgenes.csv")
