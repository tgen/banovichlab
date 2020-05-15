# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
#             Austin J. Gutierrez, agutierrez@tgen.org
# Date: 03/11/2020
# Description: Create loess plots using pseudotime ordering (Figure 5F and Figure S14)
# ==============================================================================

# ======================================
# Environment parameters
# ======================================
set.seed(12345)

# =====================================
# Load libraries
# =====================================
library(Seurat)
library(slingshot)
library(dplyr)
library(gridExtra)
library(tidyverse)

# Set date and create out folder
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# =====================================
# Read in the object(s) 
# =====================================
at2 <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200310/Final_AT2.rds")
scgb3a2 <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200310/Final_SCGB3A2.rds")
at1 <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200306/AT2_AT1_control.rds")

# =====================================
# Make subsets 
# =====================================
sub1 <- subset(at2, cells = rownames(at2@meta.data[at2@meta.data$celltype %in%
                                                     c("KRT5-/KRT17+","Transitional AT2", "AT2"), ]))

sub2 <- subset(scgb3a2, cells = rownames(scgb3a2@meta.data[scgb3a2@meta.data$celltype %in%
                                                             c("KRT5-/KRT17+","Transitional AT2", "SCGB3A2+"), ]))

# =====================================
# Convert to SingleCellExperiment
# =====================================
sub_sce1 <- as.SingleCellExperiment(sub1)
sub_sce2 <- as.SingleCellExperiment(sub2)

# =====================================
# Run Slingshot
# =====================================
sub_slingshot1 <- slingshot(sub_sce1, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="AT2", end.clus="KRT5-/KRT17+")

sub_slingshot2 <- slingshot(sub_sce2, clusterLabels = "celltype", reducedDim = 'UMAP',
                            start.clus="SCGB3A2+", end.clus="KRT5-/KRT17+")

sub_slingshot3 <- readRDS("/Volumes/scratch/agutierrez/IPF/R/Seurat/20200306/AT2_AT1_control_slingshot.rds")                       

# =====================================
# Set the pseudotime variable
# =====================================
t1 <- sub_slingshot1$slingPseudotime_1 
t2 <- sub_slingshot2$slingPseudotime_1 
t3 <- sub_slingshot3$slingPseudotime_1 

# =====================================
# Genes of interest
# =====================================
gene.list <- c("NR1D1", "SOX4", "SOX9", "ZMAT3", "MDK", "CDKN1A")
plot.val <- c(1.25, 2.5, .6, .825, 2.25,2)

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

loess_data2 <- as.data.frame(sub2@assays$SCT@data[gene.list, ])
loess_data2 <- loess_data2[, order(t2)]
temp2 <- loess_data2
temp2 <- t(temp2)
temp2 <- as.data.frame(temp2)
temp2$index <- 1:nrow(temp2)
temp2$ct = sub2@meta.data$celltype[order(t2)]

loess_data3 <- as.data.frame(at1@assays$SCT@data[gene.list, ])
loess_data3 <- loess_data3[, order(t3)]
temp3 <- loess_data3
temp3 <- t(temp3)
temp3 <- as.data.frame(temp3)
temp3$index = 1:nrow(temp3)
temp3$ct = at1@meta.data$celltype[order(t3)]

# ======================================
# FIGURE 5F
# ======================================
# Loess plot
update_geom_defaults("line", list(colour = 'blue', linetype = 0.5))
theme_set(theme_grey(base_size=6))
pdf(file = "Fig5F_loess_final_genes.pdf", width = 3, height = 1.5)
for (i in 1:length(gene.list)) {
  print(paste(i, gene.list[i], sep = "_"))
  p1 <- ggplot(temp3, aes_string(y = gene.list[i] , x = temp3[, (ncol(temp3) -1)])) + 
    geom_smooth(method=loess, level=1-1e-10) + coord_cartesian(ylim = c(0, plot.val[i])) + scale_linetype()
  #  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  p2 <- ggplot(temp1, aes_string(y = gene.list[i], x = temp1[, (ncol(temp1) -1)])) + 
    geom_smooth(method = loess, level=1-1e-10) + coord_cartesian(ylim = c(0, plot.val[i])) 
  #  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  p3 <- ggplot(temp2, aes_string(y = gene.list[i], x = temp2[, (ncol(temp2) -1)])) + 
    geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, plot.val[i])) 
  #  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  grid.arrange(p1,p2,p3, ncol=3)
}
dev.off()

# Violin plot 
my_levels <- c("AT2","SCGB3A2+","Transitional AT2","AT1","KRT5-/KRT17+")
krt5_5pop@meta.data$celltype <- factor(krt5_5pop@meta.data$celltype, levels = my_levels)

pdf(file = "Fig5F_vln_final_genes.pdf")
plot.list <- list()
for (i in 1:length(gene.list)) {
  print(paste(i, gene.list[i], sep = "_"))
  plot.list[[i]] <- VlnPlot(krt5_5pop, gene.list[i], split.by = "Status", group.by = "celltype", pt.size = 0)
}
plot.list
dev.off()

# ======================================
# FIGURE S14
# ======================================

# Genes of interest
gene.list <- c("SFTPC","AGER","KRT17","SCGB3A2","COL1A1")
plot.val <- c(8.0, 1.0, 3.0, 8.0, 2.0)

loess_data1 <- as.data.frame(sub1@assays$SCT@data[gene.list, ])
loess_data1 <- loess_data1[, order(t1)]
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

# Make the plots
interleave <- function(a, b) {
  shorter <- if (length(a) < length(b)) a else b
  longer  <- if (length(a) >= length(b)) a else b
  slen <- length(shorter)
  llen <- length(longer)
  index.short <- (1:slen) + llen
  names(index.short) <- (1:slen)
  lindex <- (1:llen) + slen
  names(lindex) <- 1:llen
  sindex <- 1:slen
  names(sindex) <- 1:slen
  index <- c(sindex, lindex)
  index <- index[order(names(index))]
  return(c(a, b)[index])
}

pdf(file = "20200514_FigS14_loess_original.pdf", width = 16, height = 8.5)
plot_list1 <- list()
plot_list2 <- list()
for (i in 1:length(gene.list)) {
  print(paste(i, gene.list[i], sep = "_"))
  plot_list1[[i]] <- ggplot(temp1, aes_string(y = gene.list[i], x = temp1[, (ncol(temp1) -1)])) + 
    geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, plot.val[i])) 
  #  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
  plot_list2[[i]] <- ggplot(temp2, aes_string(y = gene.list[i], x = temp2[, (ncol(temp2) -1)])) + 
    geom_smooth(method = loess,level=1-1e-10) + coord_cartesian(ylim = c(0, plot.val[i])) 
  #  geom_tile(aes(x = index, y= 0, color = ct, height = .1, fill=ct)) + guides(fill=guide_legend())
}
plot_list <- interleave(plot_list1, plot_list2)
plot_list <- plyr::compact(plot_list)
grid.arrange(grobs=plot_list, ncol=2)

dev.off()

# Get the pseudotime bar
p1 <- ggplot(temp1, aes_string(y = "SFTPC", x = temp1[, (ncol(temp1) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 8))  + 
  geom_tile(aes(x = index, y= 0, color = ct, height = .5, fill=ct)) + guides(fill=guide_legend())
p2 <- ggplot(temp2, aes_string(y = "SFTPC", x = temp2[, (ncol(temp1) -1)])) + geom_smooth(method = loess) + coord_cartesian(ylim = c(0, 8))  + 
  geom_tile(aes(x = index, y= 0, color = ct, height = .5, fill=ct)) + guides(fill=guide_legend())





