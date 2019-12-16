# ======================================
# Load libraries
# ======================================
library(Seurat)

# ======================================
# Read in Seurat object
# ======================================
epi <- readRDS("Epithelial.rds")

# ======================================
# Figure S: 15
# ======================================
DotPlot(epi, c("SOX9",
               "SOX4",
               "HOPX",
               "NAPSA",
               "NKX2-1",
               "CDH2",
               "VIM",
               "GDF15",
               "CDKN2A",
               "CDKN1A",
               "ITGAV",
               "TGB6",
               "MDK",
               "MMP7",
               "TM4SF1",
               "COL1A1",
               "FN1",
               "KRT17",
               "KRT5"))