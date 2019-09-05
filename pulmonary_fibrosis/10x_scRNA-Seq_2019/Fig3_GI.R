# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 20/08/2019
# Description: Code for Figure 3: G, I, J
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
library(topGO)
library(org.Hs.eg.db)
library(stats)

# ==============================================================================
# Figure 3G - Pseudotime HEATMAP
# ==============================================================================
# Read in data
epi <- readRDS("Epithelial.rds")

# Subset out only ILD samples from Epi population
epi_ild <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$Status == "ILD",]))

# KRT5-/KRT17+ and Transitional AT2 and AT2
krt5 <- subset(epi_ild, cells = rownames(epi_ild@meta.data[epi_ild@meta.data$celltype %in%
                  c("KRT5-/KRT17+","Transitional AT2", "AT2"),]))
krt5 <- FindVariableFeatures(object = krt5, verbose = F, nfeatures = 3000)
krt5 <- ScaleData(object = krt5, verbose = F)
krt5 <- RunUMAP(krt5, dims = 1:10, verbose = F)


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

# Plot the trajectory
epi_color <- as.character(krt5@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(krt5@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "AT2"] <- "#EE7342"
#onion[onion == "AT1"] <- "#F76A62"
epi_color <- onion

# ======================================
# Figure 3: G
# ======================================
plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2", "KRT5-/KRT17+", "AT1"), 
       col = c("#EE7342","#FE627D","#03AF21","#F76A62"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

lin <- getLineages(sub_sce@reducedDims@listData$UMAP, sub_sce$celltype, 
                   start.clus = "AT2", end.clus = c("KRT5-/KRT17+","AT1"))
plot(reducedDims(sub_slingshot)$UMAP, col = epi_color, asp = 1, pch = 16, cex=0.6)
lines(lin, lwd = 5, col = 'black', show.constraints = TRUE)

plot(sub_slingshot$slingPseudotime_1, col=epi_color)
plot(sub_slingshot$slingPseudotime_2, col=epi_color)
plot(sub_slingshot@int_metadata$slingshot@curves$curve1$dist_ind, col=epi_color)
plot(sub_slingshot@int_metadata$slingshot@curves$curve2$dist_ind, col=epi_color)

saveRDS(sub_slingshot,file = paste(date, "KRT5_ILD_slingshot.rds", sep = "_"))

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

# Plot the Heatmap
clab <- sample(c("#EE7342","#FE627D","#03AF21"), length(heatclus), replace = T, prob = NULL)
pdf(file = paste(date, "KRT5_pseudotime_heatmap.pdf", sep = "_"))
heatmap.2(heatdata_onion, trace = "none",
          col = color.palette, cexRow = 0.2,
          symbreaks = F, scale = "none", labCol=F, ColSideColors=epi_color,
          main = "KRT5-/KRT17, AT2, Transitional AT2", 
          Colv = F, dendrogram = "row", split = clusters) #rowsep = c(35,81,138))
dev.off()
hr <- hclust(dist(heatdata_onion), method = "average")
clusters <- dendextend::cutree(hr, k=100)

h <- heatmap.2(heatdata_onion, trace = "none",
               col = color.palette, cexRow = 0.2,
               symbreaks = F, scale = "none", labCol =F,
               main = "KRT5-/KRT17+, AT2 and Transitional AT2 heatmap", 
               Colv = F, dendrogram = "row", rowsep = c(70,123,251))

# Extract out the geneID based on dendrogram order

heatcluster <- list(transition1_genes,transition2_genes,transition3_genes,
                    transition4_genes)
heatclus <- sub_slingshot$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata_onion, heatclus, transformation = function(x){x})

# ======================================
# Figure 3: I
# ======================================
pdf("190805_test2.pdf")  # for publication with pseudotime bar on top
#makeBlankData(ce, groupsOfFeatures = heatcluster, 
 #             groupsOfSamples = features_col,
  #          nBlankFeatures = 4, nBlankSamples = 1)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            #clusterFeaturesData = heatcluster, 
            nBlankLines =2, colorScale = color.palette) 
dev.off()