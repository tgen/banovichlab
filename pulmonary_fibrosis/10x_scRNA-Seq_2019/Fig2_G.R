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
# FIGURE 2: TRAJECTORY ANALYSIS FOR AT1, TRANS AT2, AT2 AND SCGB3A2+ (both ILD and CONTROL)
# ==============================================================================
# Read in the objects 
krt5_6pop <- readRDS("/Volumes/scratch/lbui/201907_Slingshot_related/190731_krt5_6pop.rds")

# Subset 
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
saveRDS(sub_slingshot1, file = paste(date, "3A2_AT1_control&ild_slingshot.rds", sep = "_"))

sub_slingshot2 <- slingshot(sub_sce2, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="AT1")
print(SlingshotDataSet(sub_slingshot2))
saveRDS(sub_slingshot2, file = paste(date, "AT2_AT1_control&ild_slingshot.rds", sep = "_"))

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


pdf(file = paste(date, "190816_Figure2_Loess_plots", sep = "_"))
p1 <- ggplot(temp1, aes(y = AGER, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4)) 
p2 <- ggplot(temp1, aes(y = ABCA3, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 2)) 
p3 <- ggplot(temp1, aes(y = SFTPC, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 7)) 
p4 <- ggplot(temp1, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4.5)) 

p5 <- ggplot(temp2, aes(y = AGER, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4)) 
p6 <- ggplot(temp2, aes(y = ABCA3, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 2)) 
p7 <- ggplot(temp2, aes(y = SFTPC, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 7)) 
p8 <- ggplot(temp2, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4.5)) 
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
dev.off()

# ==============================================================================
# FIGURE 3: TRAJECTORY ANALYSIS FOR KRT5-/KRT17+, TRANS AT2, AT2 AND SCGB3A2+ 
# ==============================================================================
# AT2, TRANS AT2, KRT5-/KRT17+ Loess plot

# Subset 
sub2_control <- subset(sub2, cells = rownames(sub2@meta.data[sub2@meta.data$Status == "Control",]))

sub2_control <- NormalizeData(sub2_control)
sub2_control <- FindVariableFeatures(sub2_control, nfeatures = 3000)
sub2_control <- ScaleData(sub2_control)
sub2_control <- RunPCA(sub2_control)
sub2_control <- RunUMAP(sub2_control, dims = 1:7)

krt5 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("KRT5-/KRT17+","Transitional AT2", "AT2"),]))
krt5_ild <- subset(krt5, cells = rownames(krt5@meta.data[krt5@meta.data$Status == "ILD",]))

krt5_ild <- FindVariableFeatures(krt5_ild, nfeatures = 3000)
krt5_ild <- ScaleData(krt5_ild)
krt5_ild <- RunPCA(krt5_ild)
krt5_ild <- RunUMAP(krt5_ild, dims = 1:7)

# Convert to SingleCellExperiment
sub_sce3 <- as.SingleCellExperiment(sub2_control)
sub_sce4 <- as.SingleCellExperiment(krt5_ild)

# Dimentionality reduction
onion <- sub_sce3@reducedDims@listData$PCA[,1:20]
sub_sce3@reducedDims@listData$PCA <- onion

onion <- sub_sce4@reducedDims@listData$PCA[,1:20]
sub_sce4@reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce3), ignore.case = F, value = T)  # 13 MT genes
sub_sce3 <- sub_sce3[!rownames(sub_sce3) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce3), ignore.case = F, value = T) # 4723 RB genes
sub_sce3 <- sub_sce3[!rownames(sub_sce3) %in% temp2,] 

temp <- grep( "^MT-", rownames(sub_sce4), ignore.case = F, value = T) 
sub_sce4 <- sub_sce4[!rownames(sub_sce4) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce4), ignore.case = F, value = T) 
sub_sce4 <- sub_sce4[!rownames(sub_sce4) %in% temp2,] 

# Run Slingshot
sub_slingshot3 <- slingshot(sub_sce3, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="AT1")
print(SlingshotDataSet(sub_slingshot3))
saveRDS(sub_slingshot3, file = paste(date, "AT1_AT2_control_slingshot.rds", sep = "_"))

sub_slingshot4 <- slingshot(sub_sce4, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="KRT5-/KRT17+")
print(SlingshotDataSet(sub_slingshot4))
saveRDS(sub_slingshot4, file = paste(date, "AT2_KRT5_ild_slingshot.rds", sep = "_"))

# Set the pseudotime variable
t3 <- sub_slingshot3$slingPseudotime_1 
t4 <- sub_slingshot4$slingPseudotime_1 

# Create a list of interesting genes
gene.list <- c("CLDN18", "ABCA7","COBLL1","C1orf167","ITGB6","ITGB1","ITGA3","CEBPD",
               "CD44","CDKN1A","CD151","GNAS","REL","ATF3","KLF13","TFCP2L1","MYRF",
               "JUN","ETS2","NFKBIZ","FOSB","ETV5","JUNB","BHLHE40","EGR1","ELF3",
               "CREB3L1","NR1D1","IRX3","ETV1","CSRNP1","ZNF385B","BCL6","CEBPA",
               "STAT3","XBP1","NME2","HES4","ELK3","SOX9","SOX4","HES2", "SOX2", "HOPX",
               "FN1","NR1D1","EGF","NFIX","NFIC", "ETS1","SP1","SPI1","ZNF740","MZF1",
               "NFIA")

# Prepare data for loess plot for marker genes - AT2, Trans AT2, AT1 (Control conditions)
heatdata3 <- assay(sub_slingshot3)[rownames(assay(sub_slingshot3)) %in% gene.list, 
                                   order(t3, na.last = NA)]
loess_data3 = as.data.frame(sub2_control@assays$SCT@data[gene.list,])
loess_data3 = loess_data3[,order(t3)]
temp3 <- loess_data3
temp3 <- t(temp3)
temp3 = as.data.frame(temp3)
temp3$index = 1:nrow(temp3)
temp3$ct = sub2_control@meta.data$celltype[order(t3)]

# Prepare data for loess plot for marker genes - AT2, Trans AT2, KRT5-/KRT17+ (ILD conditions)
heatdata4 <- assay(sub_slingshot4)[rownames(assay(sub_slingshot4)) %in% gene.list, 
                                   order(t4, na.last = NA)]
loess_data4 = as.data.frame(krt5_ild@assays$SCT@data[gene.list,])
loess_data4 = loess_data4[,order(t4)]
temp4 <- loess_data4
temp4 <- t(temp4)
temp4 = as.data.frame(temp4)
temp4$index = 1:nrow(temp4)
temp4$ct = krt5_ild@meta.data$celltype[order(t4)]

# Plotting time
p1 <- ggplot(temp3, aes(y = NFIA, x = index)) + geom_smooth(method = loess)  +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
p2 <- ggplot(temp3, aes(y = NFIX, x = index)) + geom_smooth(method = loess)  +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) + geom_point()
p3 <- ggplot(temp3, aes(y = NFIC, x = index)) + geom_smooth(method = loess)    + 
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) + geom_point()
p4 <- ggplot(temp4, aes(y = NFIA, x = index)) + geom_smooth(method = loess)  + 
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) 
p5 <- ggplot(temp4, aes(y = NFIX, x = index)) + geom_smooth(method = loess)    +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) +geom_point()
p6 <- ggplot(temp4, aes(y = NFIC, x = index)) + geom_smooth(method = loess)  +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) + geom_point()
multiplot(p1,p2,p3,p4,p5,p6, cols=2)
+ coord_cartesian(ylim = c(0, 1.7))

### Note on gene selection
intersect(setdiff(topgenes_at1, krt5_heatmap_genes$V1), tf_database1$HGNC.symbol) #(AT1 heatmap genes and KRT5 heatmap genes, TF)
# "FOSB"    "TFCP2L1" "ETV5"    "JUNB"    "BHLHE40" "EGR1"    "CEBPA"   "ELF3"    "CREB3L1" "NR1D1"   "IRX3"   
#"ETV1"    "ZNF385B" "CSRNP1"  "BCL6" 

intersect(setdiff(krt5_heatmap_genes$V1, topgenes_at1), tf_database1$HGNC.symbol)
# "STAT3" "REL"   "XBP1"  "NME2"  "HES4"  "ELK3"  "SOX9"  "SOX4"  "HES2" 

intersect(intersect(topgenes_at1, krt5_heatmap_genes$V1), tf_database1$HGNC.symbol)
# "CEBPD" "MYRF"  "JUN"   "ETS2" # Load libraries
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

set.seed(12345)

getwd()
Sys.Date()
main_dir <- "/Volumes/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# ==============================================================================
# FIGURE 2: TRAJECTORY ANALYSIS FOR AT1, TRANS AT2, AT2 AND SCGB3A2+ (both ILD and CONTROL)
# ==============================================================================
# Read in the objects 
krt5_6pop <- readRDS("/Volumes/scratch/lbui/201907_Slingshot_related/190731_krt5_6pop.rds")

# Subset 
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
saveRDS(sub_slingshot1, file = paste(date, "3A2_AT1_control&ild_slingshot.rds", sep = "_"))

sub_slingshot2 <- slingshot(sub_sce2, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="AT1")
print(SlingshotDataSet(sub_slingshot2))
saveRDS(sub_slingshot2, file = paste(date, "AT2_AT1_control&ild_slingshot.rds", sep = "_"))

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


pdf(file = paste(date, "190816_Figure2_Loess_plots", sep = "_"))
p1 <- ggplot(temp1, aes(y = AGER, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4)) 
p2 <- ggplot(temp1, aes(y = ABCA3, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 2)) 
p3 <- ggplot(temp1, aes(y = SFTPC, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 7)) 
p4 <- ggplot(temp1, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4.5)) 

p5 <- ggplot(temp2, aes(y = AGER, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4)) 
p6 <- ggplot(temp2, aes(y = ABCA3, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 2)) 
p7 <- ggplot(temp2, aes(y = SFTPC, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 7)) 
p8 <- ggplot(temp2, aes(y = SCGB3A2, x = index)) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 4.5)) 
multiplot(p1, p2, p3, p4, p5, p6, p7, p8, cols=2)
dev.off()

# ==============================================================================
# FIGURE 3: TRAJECTORY ANALYSIS FOR KRT5-/KRT17+, TRANS AT2, AT2 AND SCGB3A2+ 
# ==============================================================================
# AT2, TRANS AT2, KRT5-/KRT17+ Loess plot

# Subset 
sub2_control <- subset(sub2, cells = rownames(sub2@meta.data[sub2@meta.data$Status == "Control",]))

sub2_control <- NormalizeData(sub2_control)
sub2_control <- FindVariableFeatures(sub2_control, nfeatures = 3000)
sub2_control <- ScaleData(sub2_control)
sub2_control <- RunPCA(sub2_control)
sub2_control <- RunUMAP(sub2_control, dims = 1:7)

krt5 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("KRT5-/KRT17+","Transitional AT2", "AT2"),]))
krt5_ild <- subset(krt5, cells = rownames(krt5@meta.data[krt5@meta.data$Status == "ILD",]))

krt5_ild <- FindVariableFeatures(krt5_ild, nfeatures = 3000)
krt5_ild <- ScaleData(krt5_ild)
krt5_ild <- RunPCA(krt5_ild)
krt5_ild <- RunUMAP(krt5_ild, dims = 1:7)

# Convert to SingleCellExperiment
sub_sce3 <- as.SingleCellExperiment(sub2_control)
sub_sce4 <- as.SingleCellExperiment(krt5_ild)

# Dimentionality reduction
onion <- sub_sce3@reducedDims@listData$PCA[,1:20]
sub_sce3@reducedDims@listData$PCA <- onion

onion <- sub_sce4@reducedDims@listData$PCA[,1:20]
sub_sce4@reducedDims@listData$PCA <- onion

# Remove mitochondria and ribosomal genes
temp <- grep( "^MT-", rownames(sub_sce3), ignore.case = F, value = T)  # 13 MT genes
sub_sce3 <- sub_sce3[!rownames(sub_sce3) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce3), ignore.case = F, value = T) # 4723 RB genes
sub_sce3 <- sub_sce3[!rownames(sub_sce3) %in% temp2,] 

temp <- grep( "^MT-", rownames(sub_sce4), ignore.case = F, value = T) 
sub_sce4 <- sub_sce4[!rownames(sub_sce4) %in% temp,]
temp2 <- grep( "^RP", rownames(sub_sce4), ignore.case = F, value = T) 
sub_sce4 <- sub_sce4[!rownames(sub_sce4) %in% temp2,] 

# Run Slingshot
sub_slingshot3 <- slingshot(sub_sce3, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="AT1")
print(SlingshotDataSet(sub_slingshot3))
saveRDS(sub_slingshot3, file = paste(date, "AT1_AT2_control_slingshot.rds", sep = "_"))

sub_slingshot4 <- slingshot(sub_sce4, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="KRT5-/KRT17+")
print(SlingshotDataSet(sub_slingshot4))
saveRDS(sub_slingshot4, file = paste(date, "AT2_KRT5_ild_slingshot.rds", sep = "_"))

# Set the pseudotime variable
t3 <- sub_slingshot3$slingPseudotime_1 
t4 <- sub_slingshot4$slingPseudotime_1 

# Create a list of interesting genes
gene.list <- c("CLDN18", "ABCA7","COBLL1","C1orf167","ITGB6","ITGB1","ITGA3","CEBPD",
               "CD44","CDKN1A","CD151","GNAS","REL","ATF3","KLF13","TFCP2L1","MYRF",
               "JUN","ETS2","NFKBIZ","FOSB","ETV5","JUNB","BHLHE40","EGR1","ELF3",
               "CREB3L1","NR1D1","IRX3","ETV1","CSRNP1","ZNF385B","BCL6","CEBPA",
               "STAT3","XBP1","NME2","HES4","ELK3","SOX9","SOX4","HES2", "SOX2", "HOPX",
               "FN1","NR1D1","EGF","NFIX","NFIC", "ETS1","SP1","SPI1","ZNF740","MZF1",
               "NFIA")

# Prepare data for loess plot for marker genes - AT2, Trans AT2, AT1 (Control conditions)
heatdata3 <- assay(sub_slingshot3)[rownames(assay(sub_slingshot3)) %in% gene.list, 
                                   order(t3, na.last = NA)]
loess_data3 = as.data.frame(sub2_control@assays$SCT@data[gene.list,])
loess_data3 = loess_data3[,order(t3)]
temp3 <- loess_data3
temp3 <- t(temp3)
temp3 = as.data.frame(temp3)
temp3$index = 1:nrow(temp3)
temp3$ct = sub2_control@meta.data$celltype[order(t3)]

# Prepare data for loess plot for marker genes - AT2, Trans AT2, KRT5-/KRT17+ (ILD conditions)
heatdata4 <- assay(sub_slingshot4)[rownames(assay(sub_slingshot4)) %in% gene.list, 
                                   order(t4, na.last = NA)]
loess_data4 = as.data.frame(krt5_ild@assays$SCT@data[gene.list,])
loess_data4 = loess_data4[,order(t4)]
temp4 <- loess_data4
temp4 <- t(temp4)
temp4 = as.data.frame(temp4)
temp4$index = 1:nrow(temp4)
temp4$ct = krt5_ild@meta.data$celltype[order(t4)]

# Plotting time
p1 <- ggplot(temp3, aes(y = NFIA, x = index)) + geom_smooth(method = loess)  +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
p2 <- ggplot(temp3, aes(y = NFIX, x = index)) + geom_smooth(method = loess)  +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) + geom_point()
p3 <- ggplot(temp3, aes(y = NFIC, x = index)) + geom_smooth(method = loess)    + 
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) + geom_point()
p4 <- ggplot(temp4, aes(y = NFIA, x = index)) + geom_smooth(method = loess)  + 
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) 
p5 <- ggplot(temp4, aes(y = NFIX, x = index)) + geom_smooth(method = loess)    +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) +geom_point()
p6 <- ggplot(temp4, aes(y = NFIC, x = index)) + geom_smooth(method = loess)  +
  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend()) + geom_point()
multiplot(p1,p2,p3,p4,p5,p6, cols=2)
+ coord_cartesian(ylim = c(0, 1.7))

### Note on gene selection
intersect(setdiff(topgenes_at1, krt5_heatmap_genes$V1), tf_database1$HGNC.symbol) #(AT1 heatmap genes and KRT5 heatmap genes, TF)
# "FOSB"    "TFCP2L1" "ETV5"    "JUNB"    "BHLHE40" "EGR1"    "CEBPA"   "ELF3"    "CREB3L1" "NR1D1"   "IRX3"   
#"ETV1"    "ZNF385B" "CSRNP1"  "BCL6" 

intersect(setdiff(krt5_heatmap_genes$V1, topgenes_at1), tf_database1$HGNC.symbol)
# "STAT3" "REL"   "XBP1"  "NME2"  "HES4"  "ELK3"  "SOX9"  "SOX4"  "HES2" 

intersect(intersect(topgenes_at1, krt5_heatmap_genes$V1), tf_database1$HGNC.symbol)
# "CEBPD" "MYRF"  "JUN"   "ETS2" 