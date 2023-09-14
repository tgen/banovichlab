# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/06/14
# Description: Supplementary Figure 8
# ==============================================================================

# ==============================================================================
# SET UP THE ENVIRONMENT VARIABLES 
# ==============================================================================
# Set up working directory
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

options(future.globals.maxSize = 4096*1024^2 )
set.seed(12345)

# Load required libraries
#Sys.unsetenv("GITHUB_PAT")
#devtools::install_github("saeyslab/nichenetr")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(nord)
library(circlize)
library(RColorBrewer)
library(tidyverse)

# ==============================================================================
# Read in the T cell object
# ==============================================================================
# Read in the Seurat object
coh.combined.sct <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")

# Subset out T cells 
t_cells <- subset(coh.combined.sct,
                      subset = Celltype1 %in% c("T cells"))

# Read in the csv file with immune gene list
icg_all <- read.csv("/scratch/lbui/Appendiceal_data/ICGs_HUetal_2021.csv", header = F)

# Heatmap with average expression
# Get cluster average expression
DefaultAssay(t_cells) <- "RNA"
var_features <- c(VariableFeatures(t_cells),icg_all$V1, marker_all$gene)
t_cells <- NormalizeData(t_cells)
t_cells <- ScaleData(t_cells, features = var_features, vars.to.regress = "percent.mt")
tcell_ave <- AverageExpression(t_cells, slot = "scale.data", assay = "RNA",
                               group.by = c("Celltype2","Pathology2"),
                               return.seurat = T)

tcell_ave@meta.data$Pathology <- sapply(strsplit(rownames(tcell_ave@meta.data),"_"), `[`, 2)
tcell_data <- tcell_ave@assays$RNA@scale.data[rownames(tcell_ave@assays$RNA@scale.data) %in% 
                                                unique(icg_all$V1),]

tcell_data <- tcell_data[which(rownames(tcell_data) %in% unique(icg_all$V1)),]

tcell_all <- GetAssayData(t_cells, assay = "RNA", slot = "scale.data")
tcell_data <- tcell_all[rownames(tcell_all) %in% icg_all$V1,]

# Create metadata for heatmap col annotation
tcell_cells_metadata <- data.frame("cell" = rownames(tcell_ave@meta.data),
                                   "Celltype" =  as.character(tcell_ave$orig.ident),
                                   "Pathology" = as.character(tcell_ave$Pathology))
tcell_cells_metadata <- tcell_cells_metadata[order(tcell_cells_metadata$Celltype,
                                                   tcell_cells_metadata$Pathology),]
tcell_data <- tcell_data[,tcell_cells_metadata$cell]

# Make heatmap
col_ha <- HeatmapAnnotation(
  df = data.frame(Celltype = tcell_cells_metadata$Celltype,
                  Pathology = tcell_cells_metadata$Pathology),
  annotation_height = unit(4, "mm"),
  col = col
)

col_fun = colorRamp2(c(-1, 0, 2), c("cadetblue4", "white", "coral2"))
Heatmap(as.matrix(tcell_data), 
        name = "Scaled_Exp", 
        column_split = tcell_cells_metadata$Celltype,
        column_title = "Cells", row_title = "Gene expression",
        col = col_fun,
        use_raster = FALSE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8))) 




