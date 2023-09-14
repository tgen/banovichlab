# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/07
# Description: M1/M2 module score for myeloid cells
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
# ==============================================================================
# SET UP THE ENVIRONMENT VARIABLES 
# ==============================================================================
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

options(future.globals.maxSize = 4096*1024^2 )
set.seed(12345)

# ======================================
# Load libraries
# ======================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(ade4)
library(Matrix)
library(ggpubr)
library(RCurl)
library(reshape2)
library(ggrepel)
library(data.table)
library(grid)
library(UpSetR)
library(ComplexHeatmap)
library(nord)
library(circlize)
library(RColorBrewer)
library(scCustomize)
library(tidyverse)

# ==============================================================================
# Read in the Appendiceal object and subset out myeloid cell population
# ==============================================================================
coh.combined.sct <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")
myeloid <- subset(coh.combined.sct, subset = Celltype1 == "Myeloid cells")

# ==============================================================================
# Supplementary Figure 10 - Macrophages M1/M2 scores
# ==============================================================================
# M1/M2 score for macrophages
macro <- subset(myeloid, subset = Celltype2 %in% c("Macrophages","SPP1+ TAMs",
                                                   "C1Qhi monocytes","Monocyte-like"))

M1_genes <- list(c("IL1A","IL1B","IL6","NOS2","TLR2","TLR4","CD80","CD86"))
M2_genes <- list(c("CSF1R","MRC1","PPARG","ARG1","CD163","CLEC10A","CLEC7A",
                   "PDCD1LG2","RETNLB"))

macro <- AddModuleScore(macro, features = M1_genes, name = "M1genes", assay = "SCT")
macro <- AddModuleScore(macro, features = M2_genes, name = "M2genes", assay = "SCT")

m1m2_data <- as.data.frame(cbind(macro@meta.data$orig.ident, macro@meta.data$Pathology2,
                                 macro@meta.data$M1genes1, macro@meta.data$M2genes1,
                                 macro@meta.data$Celltype2))
colnames(m1m2_data) <- c("Ident","Pathology","M1","M2","Celltype")
m1m2_data$M1 <- as.numeric(as.character(m1m2_data$M1))
m1m2_data$M2 <- as.numeric(as.character(m1m2_data$M2))

# Boxplot for M1M2 scores 
m1m2_melt <- melt(m1m2_data)
my_comparisons <- list(c("M2","M1"))
ggboxplot(m1m2_melt, x="variable", y="value", fill = "variable", outlier.shape = 20) +
  ggtitle ("M1/M2 gene signature score per diagnosis") +
  facet_grid(Celltype~Pathology) +
  stat_compare_means(comparisons = my_comparisons, paired = TRUE, hide.ns = TRUE,
                     method = "t.test", label.y = 0.8) +
  NoLegend() + 
  theme_bw() 
