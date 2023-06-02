#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 10/1/2022
# Description: CAR T product scGSVA scoring
#==============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(qdap)
library(googlesheets4)
library(scGSVA)

#==============================================================================
# Helper functions
#==============================================================================

source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/CART/CART_colors_themes.R")
source("/home/hnatri/CART/CART_plot_functions.R")

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)
setwd("/scratch/hnatri/CART/")

#==============================================================================
# Import Seurat objects and gene sets
#==============================================================================

product <- readRDS("/labs/banovich/BCTCSF/Heini/product_projecTILs_modules.rds")

# Canonical T cell markers
gs4_deauth()
canonical_markers  <- gs4_get(url)
sheet_names(canonical_markers)
canonical_markers <- read_sheet(canonical_markers, sheet = "T cells, gene sets")
head(canonical_markers)

canonical_markers_small <- canonical_markers[,which(colnames(canonical_markers) %in% c("RNA", "Bigger_gene_sets"))]
canonical_markers_small <- canonical_markers_small[complete.cases(canonical_markers_small),]
colnames(canonical_markers_small) <- c("GeneID", "Annot")

# Genes for module scoring
memory_markers <- setdiff(canonical_markers[which(canonical_markers$Bigger_gene_sets=="Memory"),]$RNA, c(NA))
dysfunction_markers <- setdiff(canonical_markers[which(canonical_markers$Bigger_gene_sets=="Dysfunction"),]$RNA, c(NA))
cytotoxic_markers <- setdiff(canonical_markers[which(canonical_markers$Bigger_gene_sets=="Cytotoxic"),]$RNA, c(NA))

#==============================================================================
# Seurat AddModuleScore
#==============================================================================

# Adding module scores
product <- AddModuleScore(product,
                          list(memory_markers),
                          nbin = 24,
                          ctrl = 10,
                          k = FALSE,
                          assay = "integrated_sct",
                          name = "memory",
                          seed = 1,
                          search = FALSE)

product <- AddModuleScore(product,
                          list(dysfunction_markers),
                          nbin = 24,
                          ctrl = 10,
                          k = FALSE,
                          assay = "integrated_sct",
                          name = "dysfunction",
                          seed = 1,
                          search = FALSE)

product <- AddModuleScore(product,
                          list(cytotoxic_markers),
                          nbin = 24,
                          ctrl = 10,
                          k = FALSE,
                          assay = "integrated_sct",
                          name = "cytotoxic",
                          seed = 1,
                          search = FALSE)

# Plotting
head(product@meta.data)
p1 <- FeaturePlot(product, features = c("memory1", "dysfunction1", "cytotoxic1"), reduction = "wnn.umap", ncol = 3) &
    my_theme &
    coord_fixed()

FeaturePlot(product, features = c("memory1", "dysfunction1", "cytotoxic1"), reduction = "wnn.umap", ncol = 2, split.by = "Manufacture") &
    my_theme &
    coord_fixed()

p3 <- VlnPlot(product, features = c("memory1", "dysfunction1", "cytotoxic1"), group.by = "c_cluster", pt.size = 0, cols = product_cluster_col_c, ncol = 3) &
    my_theme

p1 / p3

#==============================================================================
# scGSVA analysis
#==============================================================================
  
hsko <- buildAnnot(species="human", keytype="SYMBOL", anntype="GO")
hsko@species
hsko@anntype <- "custom"
hsko@keytype

typeof(hsko@annot)
head(hsko@annot)

hsko@annot <- as.data.frame(canonical_markers_small)

DefaultAssay(product)
res <- scgsva(product, hsko)

res_df <- as.data.frame(res)

rownames(res_df)

identical(rownames(res_df), colnames(product))

product$Cytotoxic <- res_df$Cytotoxic
product$Dysfunction <- res_df$Dysfunction
product$Memory <- res_df$Memory

saveRDS(product, "/labs/banovich/BCTCSF/Heini/product_projecTILs_modules_scGSVA.rds")

# Plotting
FeaturePlot(product, features = c("Memory", "Dysfunction", "Cytotoxic"), reduction = "wnn.umap", ncol = 3) &
    my_theme &
    coord_fixed()

FeaturePlot(product, features = c("Memory", "Dysfunction", "Cytotoxic"), reduction = "wnn.umap", ncol = 2, split.by = "Manufacture") &
    my_theme &
    coord_fixed()

VlnPlot(product, features = c("Memory", "Dysfunction", "Cytotoxic"), group.by = "c_cluster", pt.size = 0, cols = product_cluster_col_c, ncol = 3) &
    my_theme

