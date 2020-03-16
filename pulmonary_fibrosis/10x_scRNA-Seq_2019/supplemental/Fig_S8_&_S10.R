# ================================
# Author: Linh Bui, lbui@tgen.org
# Date: 2020-03-13
# Science Advances: FIGURE S8, FIGURE S10
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
# Read in the object
# ==============================================================================

epi <- readRDS("/Volumes/scratch/lbui/IPF_seurat_objects/190623_EpiA_sizereduced.rds")
pop5 <- subset(epi, cells=rownames(epi@meta.data[epi@meta.data$celltype %in%
                                                   c("AT2","AT1","KRT5-/KRT17+",
                                                     "SCGB3A2+", "Transitional AT2"),]))

# ==============================================================================
# FIGURE S8
# ==============================================================================

VlnPlot(pop5, "SCGB3A2", split.by = "Status", group.by = "celltype", pt.size = 0)
onion <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Status == "Control",]))
onion2 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Status == "ILD",]))

VlnPlot(onion, "SCGB3A2", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(onion2, "SCGB3A2", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()

# ==============================================================================
# Fig S10
# ==============================================================================

temp1 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Diagnosis == "NSIP",]))
temp2 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Diagnosis == "cHP",]))
temp3 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Diagnosis == "IPF",]))
temp4 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Diagnosis == "Unclassifiable ILD",]))
temp5 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Diagnosis == "sacroidosis",]))
temp6 <- subset(pop5, cells=rownames(pop5@meta.data[pop5@meta.data$Diagnosis == "Control",]))

VlnPlot(temp1, "COL1A1", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp2, "COL1A1", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp3, "COL1A1", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp4, "COL1A1", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp5, "COL1A1", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp6, "COL1A1", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()

VlnPlot(temp1, "KRT17", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp2, "KRT17", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp3, "KRT17", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp4, "KRT17", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp5, "KRT17", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()
VlnPlot(temp6, "KRT17", group.by = "celltype", pt.size = 0, cols = epi_col) + NoLegend()



