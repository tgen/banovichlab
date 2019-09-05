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

pdf("190805_test2.pdf")  # for publication with pseudotime bar on top
#makeBlankData(ce, groupsOfFeatures = heatcluster, 
 #             groupsOfSamples = features_col,
  #          nBlankFeatures = 4, nBlankSamples = 1)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            #clusterFeaturesData = heatcluster, 
            nBlankLines =2, colorScale = color.palette) 
dev.off()

# ----------------------------------------
## AT2, Transitional AT2 and AT1
# ----------------------------------------

# Set the pseudotime variable
t <- sub_slingshot$slingPseudotime_2 #  lineage 2: AT1

# Extract the gene expression matrix
Y <-assay(sub_slingshot)

# Fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = "./190730_Figure-4_plots/AT1_AT2_gampval_allgenes.rds")

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

# Make the heatmap plot
pdf(file = paste(date, "KRT5_heatmap.pdf", sep = "_"))
heatmap.2(heatdata_onion, trace = "none",
          col = color.palette, cexRow = 0.2,
          symbreaks = F, scale = "none", labCol = F, 
          main = "AT1, AT2 and Transitional AT2 - pseudotime heatmap", 
          Colv = F, dendrogram = "row", rowsep = c(62,136))
dev.off()

# Use cluster experiment to get the pseudotime 
heatclus <- sub_slingshot$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus, transformation = function(x){x})

# Plot the Heatmap
pdf("./190730_Figure-4_plots/190731_AT1_AT2.clusterexp.pdf") 
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            clusterFeatures = T, col=color.palette) 
dev.off()

h <- heatmap.2(heatdata_onion, trace = "none",
               col = color.palette, cexRow = 0.2,
               symbreaks = F, scale = "none", labCol = F, 
               main = "AT1, AT2 and Transitional AT2 heatmap", 
               Colv = F, dendrogram = "row", rowsep = c(70,123,251))

# ==============================================================================
# FIGURE 3E - GENE ONTOLOGY ENRICHMENT ANALYSIS
# ==============================================================================
# Extract out the geneID based on dendrogram order
allgenes <- rownames(heatdata_onion)[rev(h$rowInd)]
transition1_genes <- allgenes[1:70] 
transition2_genes <- allgenes[71:123]
transition3_genes <- allgenes[124:260]
transition4_genes <- allgenes[261:400]

# Perform TopGO - ILD - KRT5-/KRT17+
geneNames <- rownames(Y)
geneList <- factor(as.integer(geneNames %in% transition3_genes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new('topGOdata', 
              ontology = "BP", 
              allGenes = geneList, # all genes used in the pseudotime GAM pval calculation
              geneSel = transition3_genes, # genes in the 1st section of the heatmap, transition from AT2 to transitional AT2
              annot=annFUN.org, 
              mapping = 'org.Hs.eg.db',  
              ID = 'symbol')
allGO = usedGO(object = GOdata)  
allGO2 = genesInTerm(GOdata) # for geneID

resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
sig_test <- sigGenes(GOdata)
tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(allGO), numChar = 120)

# Calculate p_val_adj 
pval_adj <- p.adjust(tab$raw.p.value, method="BH")
tab$p_val_adj <- pval_adj
tab

# Save output files
write.table(tab, file = paste(date, "KRT5_bin3_TopGO_BP.csv", sep = "_"), sep = ",")

# Make a bar chart 
tab <- tab[tab$p_val_adj<0.1,]
tab <- tab[,c("GO.ID","Term","p_val_adj")]
tab$Term <- gsub(" [a-z]*\\.\\.\\.$", "", tab$Term)
tab$Term <- gsub("\\.\\.\\.$", "", tab$Term)
tab$Term <- paste(tab$GO.ID, tab$Term, sep=", ")
tab$Term <- factor(tab$Term, levels=rev(tab$Term))
tab$p_val_adj <- as.numeric(tab$p_val_adj)

bin_2 <- ggplot(tab, aes(x=Term, y=-log10(p_val_adj))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(tab$p_val_adj)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

# ==============================================================================
# FIGURE 3: VIOLIN PLOT FOR HEATMAP GENES
# ==============================================================================
krt5 <- readRDS("20190806_KRT5_ILD&Control.rds")

my_levels <- c("AT2","Transitional AT2","KRT5-/KRT17+")
krt5@meta.data$celltype <- factor(krt5@meta.data$celltype, levels = my_levels)

# Look for share genes between the pseudotime genes and DEG for AT2 and Trans AT2
at2_deg <- read.csv("~/Desktop/IPF_project/DE_analysis/Disease_vs_Control/AT2 _disease_vs_control .csv")
transat2_deg <- read.csv("~/Desktop/IPF_project/DE_analysis/Disease_vs_Control/Transitional AT2 _disease_vs_control .csv")

intersect(transition1_genes, at2_deg[at2_deg$abs.avg_logFC >= .5,]$X)
#[1] "CEBPD" "DUSP6" "FASN"  "PGC"   "DMBT1"
VlnPlot(krt5, c("DUSP6", "FASN","DMBT1"),split.by = "Status", 
        pt.size = 0, group.by = "celltype")

intersect(transition2_genes, transat2_deg[transat2_deg$abs.logFC. >= .5,]$X)
#[1] "SEPW1"   "MIF"     "DYNLL1"  "TAGLN2"  "CRIP2"   "CRIP1"   "CALM2"   "TMSB4X"  
#"S100A10" "TRAM1"   "CD9"     "CEACAM6" "CLIC5" 

intersect(transition3_genes, transat2_deg[transat2_deg$abs.logFC. >= .5,]$X)
#[1] "CCND1"  "CCND2"  "VIM"    "LGALS1" "KRT8"   "IL32"   "ANXA2" 

# ==============================================================================
# TRAJECTORY ANALYSIS FOR AT1, TRANS AT2 AND AT2 (CONTROL ONLY) - SUPPLEMENTARY
# ==============================================================================

# Read in the objects 
control <- readRDS("/Users/lbui/Desktop/201907_Slingshot_related/Only_Control/190806_control_only.rds")

# Convert to SingleCellExperiment
sub_sce <- as.SingleCellExperiment(control)

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
                           start.clus="AT2", end.clus=c("AT1","KRT5-/KRT17+"))
summary(sub_slingshot$slingPseudotime_1)
print(SlingshotDataSet(sub_slingshot))

# Plot the trajectory
epi_color <- as.character(control@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(control@meta.data))

onion <- epi_color
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "AT2"] <- "#EE7342"
onion[onion == "AT1"] <- "#F76A62"
epi_color <- onion

plot(reducedDims(sub_slingshot)$UMAP, col = epi_color,
     pch = 16, asp = 1, xlab = "UMAP_1", ylab = "UMAP_2", cex=0.6)   
legend(par(xpd = T), x= "topleft", pch = c(20), 
       legend = c("AT2","Transitional AT2",  "AT1", "KRT5-/KRT17+"), 
       col = c("#EE7342","#FE627D","#F76A62","#03AF21"), bty = 'n')
lines(slingCurves(sub_slingshot)$curve1, lwd=3)
lines(slingCurves(sub_slingshot)$curve2, lwd=3)

saveRDS(sub_slingshot,file = paste(date, "AT1_control_slingshot.rds", sep = "_"))

# Set the pseudotime variable
t <- sub_slingshot$slingPseudotime_2 

# Extract the gene expression matrix
Y <-assay(sub_slingshot)

# Fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = paste(date, "AT1_control_gampval_allgenes.rds", sep = "_"))

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

mean_matrix = matrix(nrow = 400, ncol = 164)
j=0
potato = for(i in seq(1,4095, 25)){
  j=j+1
  if(i != 4076){
    k=i+25
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
  else{
    k=i+19
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
}
colnames(mean_matrix) = 1:164

mean_matrix = as.data.frame(mean_matrix)

gene_order = as.numeric(apply(mean_matrix, 1, function(xx){names(xx)[xx == max(xx)][1] }))
heatdata_onion = heatdata[order(gene_order),]
color.palette  <- colorRampPalette(c("blue","black","yellow"), space = "Lab")(n=100)

# Create a Cluster Experiment object for heatmap
heatclus <- sub_slingshot$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata_onion, heatclus, transformation = function(x){x})

# Plot the Heatmap
pdf(file = paste(date, "AT1_pseudotime_heatmap.pdf", sep = "_"))
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            colorScale=color.palette)
dev.off()


# ==============================================================================
# FIGURE S6: HEATMAP TRAJECTORY ANALYSIS FOR AT1, TRANS AT2, AT2 AND SCGB3A2+ (both ILD and CONTROL)
# ==============================================================================

# Read in the slingshot object, created on 20190813 (described in desktop/IPF_related/Scripts/20190813_SCGB3A2_trajectory.R)

sub_slingshot1 <- readRDS("/scratch/lbui/RStudio_folder/20190813/20190813_3A2_AT1_control&ild_slingshot.rds")
sub_slingshot2 <- readRDS("/scratch/lbui/RStudio_folder/20190813/20190813_AT2_AT1_control&ild_slingshot.rds")

# --------------------------------------------------
# SCGB3A2, Trans. AT2 and AT1 heatmap
# -------------------------------------------------

# Set the pseudotime variable
t <- sub_slingshot1$slingPseudotime_1 

# Extract the gene expression matrix
Y <-assay(sub_slingshot1)

# Fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = paste(date, "3A2_AT1_gampval_allgenes.rds", sep = "_"))

# Select the top 400 genes that change over pseudotime
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:400] 
heatdata <- assay(sub_slingshot1)[rownames(assay(sub_slingshot1)) %in% topgenes, 
                                  order(t, na.last = NA)]
allgenes <- factor(rownames(Y))

# Scale the data
heatdata <- t(scale(t(heatdata)))

# Trim z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -3] = -3

# Order genes based on average expression in 25 cells/cell type

mean_matrix = matrix(nrow = 400, ncol = 207)
j=0
potato = for(i in seq(1,5151, 25)){
  j=j+1
  if(i != 5151){
    k=i+25
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
  else{
    k=i+0
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
}
colnames(mean_matrix) = 1:207

mean_matrix = as.data.frame(mean_matrix)

gene_order = as.numeric(apply(mean_matrix, 1, function(xx){names(xx)[xx == max(xx)][1] }))
heatdata_onion = heatdata[order(gene_order),]
color.palette  <- colorRampPalette(c("blue","black","yellow"), space = "Lab")(n=100)

# Create a Cluster Experiment object for heatmap
heatclus <- sub_slingshot1$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus, transformation = function(x){x})

# Plot the Heatmap
pdf(file = paste(date, "AT1_pseudotime_heatmap.pdf", sep = "_"))
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            colorScale=color.palette, treeheight=0)
dev.off()

# ---------------------------------------------------
# AT2, Trans. AT2 and AT1 heatmap
# ---------------------------------------------------

# Set the pseudotime variable
t <- sub_slingshot2$slingPseudotime_1 

# Extract the gene expression matrix
Y <-assay(sub_slingshot2)

# Fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}) 

saveRDS(gam.pval, file = paste(date, "AT2_AT1_gampval_allgenes.rds", sep = "_"))

# Select the top 400 genes that change over pseudotime
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:400] 
heatdata <- assay(sub_slingshot2)[rownames(assay(sub_slingshot2)) %in% topgenes, 
                                  order(t, na.last = NA)]
allgenes <- factor(rownames(Y))

# Scale the data
heatdata <- t(scale(t(heatdata)))

# Trim z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -3] = -3

# Order genes based on average expression in 25 cells/cell type

mean_matrix = matrix(nrow = 400, ncol = 450)
j=0
potato = for(i in seq(1,11242, 25)){
  j=j+1
  if(i != 11226){
    k=i+25
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
  else{
    k=i+16
    tmp = heatdata[,i:k]
    mean_matrix[,j] = rowMeans(tmp)
  }
}
colnames(mean_matrix) = 1:450

mean_matrix = as.data.frame(mean_matrix)

gene_order = as.numeric(apply(mean_matrix, 1, function(xx){names(xx)[xx == max(xx)][1] }))
heatdata_onion = heatdata[order(gene_order),]
color.palette  <- colorRampPalette(c("blue","black","yellow"), space = "Lab")(n=100)

# Create a Cluster Experiment object for heatmap
heatclus <- sub_slingshot2$celltype[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata_onion, heatclus, transformation = function(x){x})

# Plot the Heatmap
pdf(file = paste(date, "AT1_pseudotime_heatmap2.pdf", sep = "_"))
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
            visualizeData = 'transformed', fontsize=15,annLegend=T,
            colorScale=color.palette, treeheight=0)
dev.off()



