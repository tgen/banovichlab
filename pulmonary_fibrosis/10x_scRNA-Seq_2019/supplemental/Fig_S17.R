# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
#             Austin J. Gutierrez, agutierrez@tgen.org
# Date: 20/08/2019
# Description: Create loess plots using pseudotime ordering
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
set.seed(12345)
main_dir <- "/Users/agutierrez/Documents/projects/single_cell/IPF/loess"
date <- gsub("-", "", Sys.Date())
dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

# =====================================
# Load libraries
# =====================================
library(Seurat)
library(slingshot)
library(dplyr)
library(scater)

# =====================================
# Read in the object(s) 
# =====================================
epi <- readRDS("~/Documents/projects/single_cell/IPF/loess/epi.rds")


# =====================================
# Make subsets 
# =====================================
sub1 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$celltype %in%
                                                     c("KRT5-/KRT17+","Transitional AT2", "SCGB3A2+"), ]))

sub1 <- subset(sub1, cells = rownames(sub1@meta.data[sub1@meta.data$Status %in%
                                                       c("ILD"), ]))

sub2 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$celltype %in%
                                                                 c("KRT5-/KRT17+","Transitional AT2", "AT2"), ]))

sub2 <- subset(sub2, cells = rownames(sub2@meta.data[sub2@meta.data$Status %in%
                                                               c("ILD"), ]))

# =====================================
# Seurat workflow
# =====================================
sub1 <- FindVariableFeatures(sub1, nfeatures = 3000)
sub1 <- ScaleData(sub1)
sub1 <- RunPCA(sub1)
sub1 <- RunUMAP(sub1, dims = 1:10)

sub2 <- FindVariableFeatures(sub2, nfeatures = 3000)
sub2 <- ScaleData(sub2)
sub2 <- RunPCA(sub2)
sub2 <- RunUMAP(sub2, dims = 1:7)

# =====================================
# Convert to SingleCellExperiment
# =====================================
sub_sce1 <- as.SingleCellExperiment(sub1)
sub_sce2 <- as.SingleCellExperiment(sub2)

# =====================================
# Run Slingshot
# =====================================
sub_slingshot1 <- slingshot(sub_sce1, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="SCGB3A2+", end.clus="KRT5-/KRT17+")

sub_slingshot2 <- slingshot(sub_sce2, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="KRT5-/KRT17+")


# =====================================
# Set the pseudotime variable
# =====================================
t1 <- sub_slingshot1$slingPseudotime_1 
t2 <- sub_slingshot2$slingPseudotime_1 

# =====================================
# Genes of interest
# =====================================
gene.list <- c("SFTPC","AGER", "KRT17", "SCGB3A2", "COL1A1")
plot.val1 <- c(3, 0.8, 3, 4.25, 2)
plot.val2 <- c(7, 0.8, 2.75, 2, 2)

# =====================================
# Prepare date for loess plots
# =====================================
loess_data1 <- as.data.frame(sub1@assays$SCT@data[gene.list, ])
loess_data1 <- loess_data1[,order(t1)]
temp1 <- loess_data1
temp1 <- t(temp1)
temp1 <- as.data.frame(temp1)
temp1$index <- 1:nrow(temp1)
temp1$ct <- sub1@meta.data$celltype[order(t1)]

loess_data2 <- as.data.frame(sub2@assays$SCT@data[gene.list, ])
loess_data2 <- loess_data2[, order(t2)]
temp2 <- loess_data2
temp2 <- t(temp2)
temp2 <- as.data.frame(temp2)
temp2$index <- 1:nrow(temp2)
temp2$ct = sub2@meta.data$celltype[order(t2)]

plot.val1 <- c(3, 0.8, 3, 4.25, 2)
plot.val2 <- c(7, 0.8, 2.75, 2, 2)
test <- c(sort(temp1[, 1])[.80*length(temp1[, 1])],
          quantile(temp1[,2], 0.80),
          quantile(temp1[,3], 0.80), 
          quantile(temp1[,4], 0.80), 
          quantile(temp1[,5], 0.80))

max(temp1[,1])
# =====================================
# Its loess time
# =====================================
pdf(file = paste(date, "loess_sup_genes.pdf", sep = "_"), width = 16, height = 8.5)
plot_list <- list()
j = 1
k = 6
for (i in 1:length(gene.list)) {
  print(paste(i, gene.list[i], sep = "_"))
  p1 <- ggplot(temp1, aes_string(y = gene.list[i] , x = temp1[, (ncol(temp1) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, plot.val1[i]))
  #geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  p2 <- ggplot(temp2, aes_string(y = gene.list[i], x = temp2[, (ncol(temp2) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, plot.val2[i]))
  #geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  plot_list[[j]] <- p1
  j = j + 1
  plot_list[[k]] <- p2
  k = k + 1
}
multiplot(plotlist = plot_list, cols = 2)

dev.off()





