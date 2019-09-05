# ==============================================================================
# Author(s) : Austin J. Gutierrez, agutierrez@tgen.org
# Date: 20/08/2019
# Description: Figure 2 A, B
# ==============================================================================
# ======================================
# Load libraries
# ======================================
library(Seurat)

# ======================================
# Read in Seurat object
# ======================================
epi <- readRDS("Epithelial.rds")

# ======================================
# Figure 2: A
# ======================================
DimPlot(epi, group.by = "celltype")

# ======================================
# Figure 2: B
# ======================================
gene_list <- c("AGER", "ABCA3", "SFTPC", "SCGB3A2", "SCGB1A1", "MUC5B", "KRT5")

VlnPlot(object = epi,
        ncol = 1,
        features = gene_list,
        pt.size = 0) + NoLegend()