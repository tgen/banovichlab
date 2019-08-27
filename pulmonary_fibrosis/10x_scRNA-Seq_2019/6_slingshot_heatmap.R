# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 11/06/2019
# Description: 
# ==============================================================================
# Load libraries
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
library(topGO)
library(org.Hs.eg.db)
library(stats)

set.seed(12345)

getwd()
Sys.Date()
main_dir <- "/Volumes/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# ==============================================================================
# Figure 3G - Pseudotime HEATMAP
# ==============================================================================

# Read in data
krt5 <- readRDS("/Volumes/scratch/lbui/RStudio_folder/Related_files/190807_krt5_ild.rds")

# Convert to SingleCellExperiment
sub_sce <- as.SingleCellExperiment(krt5)

# Dimentionality reduction
onion <- sub_sce@reducedDims@listData$PCA[,1:20]
sub_sce@reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce), ignore.case = F, value = T) #13 MT genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce), ignore.case = F, value = T) #4723 RP genes
sub_sce <- sub_sce[!rownames(sub_sce) %in% temp2,] # total 21641 genes, 6406 cells

# Run Slingshot
sub_slingshot <- slingshot(sub_sce, clusterLabels = "celltype", reducedDim = 'UMAP',
                           start.clus="AT2", end.clus="KRT5-/KRT17+")
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Figure 3G
# Plot the trajectory
epi_color <- as.character(krt5@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(krt5@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "AT2"] <- "#EE7342"
#onion[onion == "AT1"] <- "#F76A62"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "KRT5-/KRT17+", "AT1"), 
       col = c("#EE7342","#FE627D","#03AF21","#F76A62"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

# ----------------------------------------
## AT2, Transitional AT2 and KRT5-/KRT17+
# ----------------------------------------
# Set the pseudotime variable
t <- sub_slingshot$slingPseudotime_1 #  lineage 1: KRT5-/KRT17

# Extract the gene expression matrix
Y <-assay(sub_slingshot)

# Fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = paste(date, "KRT5_AT2_gampval_allgenes.rds", sep = "_"))

# Select the top 400 genes that change over pseudotime
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:400] 
heatdata <- assay(sub_slingshot)[rownames(assay(sub_slingshot)) %in% topgenes, 
                                 order(t, na.last = NA)]
allgenes <- factor(rownames(Y))

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

# Create a Cluster Experiment object for heatmap
heatclus <- sub_slingshot$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata_onion, heatclus, transformation = function(x){x})
​
# Plot the Heatmap
pdf(file = paste(date, "KRT5_3A2_pseudotime_heatmap.CE.pdf", sep = "_"))
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            clusterFeatures = T, colorScale=color.palette, treeheight = 0)
dev.off()
## Manually split cluster for adding a blank line between groups
h <- heatmap.2(heatdata_onion, trace = "none",
               col = color.palette, cexRow = 0.2,
               symbreaks = F, scale = "none", labCol=F, 
               main = "KRT5-/KRT17, AT2, Transitional AT2", 
               Colv = F, dendrogram = "row",rowsep = c(70,123,260))
saveRDS(h, file = paste(date, "KRT5_pseudotime_heatmap.rds", sep = "_"))
​
ind <- rev(h$rowInd)
transition1_ind <- ind[1:70] 
transition2_ind <- ind[71:123]
transition3_ind <- ind[124:260]
transition4_ind <- ind[261:400]
cluster_split <- list(transition1_ind,transition2_ind,transition3_ind,transition4_ind)
​
pdf(file = paste(date, "KRT5_pseudotime_heatmap.split.pdf", sep = "_"))
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            clusterFeatures = T, colorScale=color.palette, treeheight = 0,
            clusterFeaturesData = cluster_split, nBlankLines = 1) 
dev.off()
