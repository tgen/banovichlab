# ======================================
# Load libraries
# ======================================
library(Seurat)

# ======================================
# Read in Seurat object
# ======================================
ild <- readRDS("ILD.rds")
epi <- readRDS("Epithelial.rds")
endo <- readRDS("Endothelial.rds")
immune <- readRDS("Immune.rds")
meso <- readRDS("Mesochymal.rds")

# ======================================
# Figure S: 3A
# ======================================
FeaturePlot(ild, c("PTPRC",
                   "EPCAM",
                   "PECAM1"))

# ======================================
# Figure S: 3C
# ======================================
ElbowPlot(endo)
ElbowPlot(epi)
ElbowPlot(immune)
ElbowPlot(meso)