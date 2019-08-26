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
krt5_6pop <- readRDS("/Volumes/scratch/lbui/201907_Slingshot_related/190731_krt5_6pop.rds")

# =====================================
# Make subsets 
# =====================================
sub1 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("KRT5-/KRT17+","Transitional AT2", "SCGB3A2+"), ]))

sub1 <- subset(sub1, cells = rownames(sub1@meta.data[sub1@meta.data$Status %in%
                                                                 c("ILD"), ]))

sub2 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("AT1","Transitional AT2", "AT2"), ]))

sub2_control <- subset(sub2, cells = rownames(sub2@meta.data[sub2@meta.data$Status %in%
                                                                 c("Control"), ]))

no_basal <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[!krt5_6pop@meta.data$celltype %in% 
                                                                 c("Basal"), ]))

krt5 <- subset(krt5_6pop, cells = rownames(krt5_6pop@meta.data[krt5_6pop@meta.data$celltype %in%
                                                                 c("KRT5-/KRT17+","Transitional AT2", "AT2"), ]))

krt5_ild <- subset(krt5, cells = rownames(krt5@meta.data[krt5@meta.data$Status %in% 
                                                           c("ILD"), ]))

# =====================================
# Reorder levels for plotting
# =====================================
my_levels <- c("Ciliated","Differentiating Ciliated","SCGB3A2+ SCGB1A1+",
               "SCGB3A2+","MUC5B+","MUC5AC+ High","Proliferating Epithelial Cells",
               "Basal","AT2","Transitional AT2","AT1","KRT5-/KRT17+","Fibroblasts",
               "Myofibroblasts","PLIN2+ Fibroblasts","HAS1 High Fibroblasts")

krt5_6pop@meta.data$celltype <- factor(krt5_6pop@meta.data$celltype, levels = my_levels)

my_levels <- c("AT2","SCGB3A2+","Transitional AT2","AT1","KRT5-/KRT17+")

no_basal@meta.data$celltype <- factor(no_basal@meta.data$celltype, levels = my_levels)

# =====================================
# Seurat workflow
# =====================================
sub1 <- FindVariableFeatures(sub1, nfeatures = 3000)
sub1 <- ScaleData(sub1)
sub1 <- RunPCA(sub1)
sub1 <- RunUMAP(sub1, dims = 1:10)

sub2_control <- FindVariableFeatures(sub2_control, nfeatures = 3000)
sub2_control <- ScaleData(sub2_control)
sub2_control <- RunPCA(sub2_control)
sub2_control <- RunUMAP(sub2_control, dims = 1:7)

krt5_ild <- FindVariableFeatures(krt5_ild, nfeatures = 3000)
krt5_ild <- ScaleData(krt5_ild)
krt5_ild <- RunPCA(krt5_ild)
krt5_ild <- RunUMAP(krt5_ild, dims = 1:7)

# =====================================
# Convert to SingleCellExperiment
# =====================================
sub_sce1 <- as.SingleCellExperiment(sub1)
sub_sce3 <- as.SingleCellExperiment(sub2_control)
sub_sce4 <- as.SingleCellExperiment(krt5_ild)

# =====================================
# Run Slingshot
# =====================================
sub_slingshot1 <- slingshot(sub_sce1, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="SCGB3A2+", end.clus="KRT5-/KRT17+")

sub_slingshot3 <- slingshot(sub_sce3, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="AT1")

sub_slingshot4 <- slingshot(sub_sce4, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="KRT5-/KRT17+")                            

# =====================================
# Set the pseudotime variable
# =====================================
t1 <- sub_slingshot1$slingPseudotime_1 
t3 <- sub_slingshot3$slingPseudotime_1 
t4 <- sub_slingshot4$slingPseudotime_1 

# =====================================
# Genes of interest
# =====================================
gene.list <- c("NR1D1", "SOX4", "SOX9", "AREG", "ZMAT3", "PMEPA1", "TPM1", "MDK")
plot.val <- c(1.25, 2.5, .6, 2, .825, 1.25, 2.25, 2.25)

# =====================================
# Prepare date for loess plots
# =====================================
loess_data1 <- as.data.frame(sub1@assays$SCT@data[gene.list, ])
loess_data1 <- loess_data1[, order(t1)]
temp1 <- loess_data1
temp1 <- t(temp1)
temp1 <- as.data.frame(temp1)
temp1$index <- 1:nrow(temp1)
temp1$ct <- sub1@meta.data$celltype[order(t1)]

loess_data3 <- as.data.frame(sub2_control@assays$SCT@data[gene.list, ])
loess_data3 <- loess_data3[, order(t3)]
temp3 <- loess_data3
temp3 <- t(temp3)
temp3 <- as.data.frame(temp3)
temp3$index <- 1:nrow(temp3)
temp3$ct = sub2_control@meta.data$celltype[order(t3)]

loess_data4 <- as.data.frame(krt5_ild@assays$SCT@data[gene.list, ])
loess_data4 <- loess_data4[, order(t4)]
temp4 <- loess_data4
temp4 <- t(temp4)
temp4 <- as.data.frame(temp4)
temp4$index = 1:nrow(temp4)
temp4$ct = krt5_ild@meta.data$celltype[order(t4)]

# =====================================
# Its loess time
# =====================================
pdf(file = paste(date, "loess_final_genes.pdf", sep = "_"), width = 16, height = 8.5)
for (i in 1:length(gene.list)) {
  print(paste(i, gene.list[i], sep = "_"))
  p1 <- ggplot(temp3, aes_string(y = gene.list[i] , x = temp3[, (ncol(temp3) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, plot.val[i]))
  #geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  p2 <- ggplot(temp4, aes_string(y = gene.list[i], x = temp4[, (ncol(temp4) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, plot.val[i]))
  #geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  p3 <- ggplot(temp1, aes_string(y = gene.list[i], x = temp1[, (ncol(temp1) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, plot.val[i]))
  #geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  multiplot(p1, p2, p3, cols=3)
}
dev.off()

# =====================================
# Vln plots
# =====================================
pdf(file = paste(date, "vln_final_genes.pdf", sep = "_"))
plot.list <- list()
for (i in 1:length(gene.list)) {
  print(paste(i, gene.list[i], sep = "_"))
  plot.list[[i]] <- VlnPlot(no_basal, gene.list[i], split.by = "Status", group.by = "celltype", pt.size = 0)
}
plot.list
dev.off()
