# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/06/14
# Description: Figure 2 - Epithelial cells + Fig S4c
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
library(Seurat)
library(dplyr)
library(ggplot2)
library(ade4)
library(ggpubr)
library(RCurl)
library(reshape2)
library(ggrepel)
library(data.table)
library(grid)
library(ComplexHeatmap)
library(nord)
library(circlize)
library(RColorBrewer)
library(scCustomize)
library(tidyverse)

# ==============================================================================
# Read in the Epithelial with infercnv data
# ==============================================================================
epi <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_Epithelial_infercnv_HMMi3.rds")

# ==============================================================================
# Figure 2a
# ==============================================================================
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
DimPlot(epi, group.by = "Pathology2", cols = pathology_col, label = TRUE, repel = TRUE) +
  NoLegend()

# ==============================================================================
# Figure 2b
# ==============================================================================
# Make Featureplot
cnv_plot <- c("proportion_dupli_chr11","proportion_loss_chr11", "proportion_dupli_chr21",
              "proportion_loss_chr4","proportion_dupli_chr9", "proportion_dupli_chr16",
              "proportion_dupli_chr19","proportion_dupli_chr20", "proportion_dupli_chr7", 
              "proportion_loss_chr8", "proportion_dupli_chr1","proportion_dupli_chr13")
FeaturePlot(epi, features = cnv_plot, ncol = 4) + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

# ==============================================================================
# Figure 2c
# ==============================================================================
epi_cols <- c("Enterocytes" = "#96a925", "Goblet-like cells" = "#00c085",
              "MUC5Bhi cells" = "#00a9fb", "SPINK4hi cells" = "#ff6bcb")
DimPlot(epi, group.by = "Celltype2", cols = epi_cols, label = TRUE, repel = TRUE) +
  NoLegend()

# ==============================================================================
# Figure 2d
# ==============================================================================
Stacked_VlnPlot(epi, features = c("FABP1","EPCAM","CEACAM5","CEACAM6","SPINK4",
                                  "MUC2","REG4","MUC5B","CA2"),
                group.by = "Celltype2", pt.size = 0, 
                colors_use = epi_cols, assay="SCT")

# Supplementary Figure 4c
Stacked_VlnPlot(epi, features = c("FABP1","EPCAM","CEACAM5","CEACAM6","SPINK4",
                                  "MUC2","TFF3","REG4","SPDEF","MUC5B","CA2"),
                group.by = "Pathology2", pt.size = 0, 
                colors_use = pathology_col, assay="SCT")

# ==============================================================================
# Figure 2e
# ==============================================================================
library(speckle)
# Cell proportion
propeller_test <- propeller(clusters = epi$Celltype2, 
                            sample = epi$orig.ident, 
                            group = epi$Pathology)

write.csv(propeller_test, file = "Epithelial_CT2_propellertest.csv")

# Plot cell type proportions
plotCellTypeProps(clusters=epi$Celltype2, sample=epi$Pathology2) 

# ==============================================================================
# Figure 2f - DEG MUC5B+ vs. SPINK4+ and other goblet-like cells
# ==============================================================================
# Subset out the enterocytes
subset_ob <- subset(epi, subset = Celltype2 != "Enterocytes")

# Run FindMarkers to find DEGs MUC5Bhi vs. other cancerous goblet cells
muc5b.1 <- FindMarkers(epi, group.by = "Celltype2",
                     ident.1 = "MUC5B+ goblet cells",
                     ident.2 = "Goblet cells",
                     test.use = "negbinom",
                     assay = "RNA",
                     latent.vars = "Flowcell",
                     logfc.threshold = 0)

write.csv(muc5b.1, file = "Appendiceal_Epithelial_MUC5B_vs_cancerousgoblet_DEGs_RNA.csv")

# Run FindMarkers to find DEGs MUC5Bhi vs. SPINK4hi goblet cells
muc5b.2 <- FindMarkers(epi, group.by = "Celltype2",
                     ident.1 = "MUC5B+ goblet cells",
                     ident.2 = "SPINK4+ goblet cells",
                     test.use = "negbinom",
                     assay = "RNA",
                     latent.vars = "Flowcell",
                     logfc.threshold = 0)

write.csv(muc5b.2, file = "Appendiceal_Epithelial_MUC5B_vs_SPNIK4goblet_DEGs_RNA.csv")

# Run FindMarkers to find DEGs SPINK4hi vs. goblet cells
muc5b.3 <- FindMarkers(epi, group.by = "Celltype2",
                       ident.1 = "SPINK4+ goblet cells",
                       ident.2 = "Goblet cells",
                       test.use = "negbinom",
                       assay = "RNA",
                       latent.vars = "Flowcell",
                       logfc.threshold = 0)

write.csv(muc5b.3, file = "Appendiceal_Epithelial_SPINK4_vs_goblet_DEGs_RNA.csv")

# Heatmap for the top 20 DEGs in each comparison
muc5b.1_genes <- rownames(muc5b.1[order(-muc5b.1$avg_log2FC),][1:20,])
muc5b.2_genes <- rownames(muc5b.2[order(-muc5b.2$avg_log2FC),][1:20,])
muc5b.3_genes <- rownames(muc5b.3[order(-muc5b.3$avg_log2FC),][1:20,])

hm_genes <- unique(c(muc5b.1_genes, muc5b.2_genes, muc5b.3_genes))
hm_genes <- top_genes$gene
ribo <- grep("^RP", hm_genes, value = TRUE)
mito <- grep("^MT", hm_genes, value = TRUE)
hm_genes <- hm_genes[!hm_genes %in% c(ribo,mito)]
hm_sub <- subset(epi, subset = Celltype2 != "Enterocytes")
DefaultAssay(hm_sub) <- "RNA"
hm_sub <- FindVariableFeatures(hm_sub)
hm_sub <- NormalizeData(hm_sub)
hm_sub <- ScaleData(hm_sub, features = rownames(hm_sub))

epi_cols2 <- c("Goblet cells" = "#00c085",
              "MUC5B+ goblet cells" = "#00a9fb", "SPINK4+ goblet cells" = "#ff6bcb")
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
col <- list(Celltype_col = epi_cols2,
            pathology_col = pathology_col)

epi_all_markers <- GetAssayData(hm_sub,assay = "RNA", slot = "scale.data")
epi_hm_markers <- epi_all_markers[which(rownames(epi_all_markers) %in% c(unique(hm_genes),"SPINK4")),]

## Create metadata for heatmap col annotation
epi_cells_metadata <- data.frame("cell" = rownames(hm_sub@meta.data),
                                 "Celltype" =  as.character(hm_sub$Celltype2),
                                 "Pathology" = as.character(hm_sub$Pathology2))
epi_cells_metadata <- epi_cells_metadata[order(epi_cells_metadata$Celltype,
                                               epi_cells_metadata$Pathology),]
epi_hm_markers <- epi_hm_markers[,epi_cells_metadata$cell]

## Set up col and row annotation
col_ha <- HeatmapAnnotation(
  df = data.frame(Celltype = epi_cells_metadata$Celltype,
                  Pathology = epi_cells_metadata$Pathology),
  annotation_height = unit(4, "mm"),
  col = col
)

# Heatmap colors
col_fun = colorRamp2(c(-2, 0, 2), c("cadetblue4", "white", "coral2"))

Heatmap(as.matrix(epi_hm_markers), 
        name = "Scaled_Exp", 
        column_split = epi_cells_metadata$Celltype,
        column_title = "Cells", row_title = "Gene expression",
        col = col_fun,
        use_raster = TRUE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        top_annotation = col_ha,
       # km = 5, gap = unit(2, "mm"),
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8))
) 

# ==============================================================================
# Figure 2g - DEGs LAMN vs LGMA in SPINK4hicells
# ==============================================================================
# Subset out the SPINK4+ cells
spink4_sub <- subset(epi, subset = Celltype2 == "SPINK4hi cells")

# Run FindMarkers
spink4_deg <- FindMarkers(spink4_sub,
                          group.by = "Pathology2",
                          ident.1 = "LAMN",
                          ident.2 = "LGMA",
                          test.use = "negbinom",
                          assay = "RNA",
                          slot = "counts",
                          logfc.threshold = 0)

# Make Volcano plots
spink4_deg$geneID <- rownames(spink4_deg)
ribo <- grep("^RP", spink4_deg$geneID, value = TRUE)

spink4_deg$diffexpressed <- "NO"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 
spink4_deg$diffexpressed[spink4_deg$avg_log2FC > 2 & spink4_deg$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.05, set as "DOWN"
spink4_deg$diffexpressed[spink4_deg$avg_log2FC < -2 & spink4_deg$p_val_adj < 0.05] <- "DOWN"
spink4_deg$diffexpressed[spink4_deg$geneID %in% ribo] <- "NO"
spink4_deg$delabel <- NA
spink4_deg$delabel[spink4_deg$diffexpressed != "NO"] <- spink4_deg$geneID[spink4_deg$diffexpressed != "NO"]

ggplot(spink4_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                    col = diffexpressed, label = delabel)) + 
  geom_point(size = 3) + 
  geom_text_repel(max.overlaps = 20) + 
  theme_bw() + 
  scale_color_manual(values = c("blue", "grey80", "red")) + 
  geom_vline(xintercept=c(-2, 2), col="black", linetype="dashed")  
write.csv(spink4_deg, file = "Epithelial_SPINK4cells_LAMN_LGMA_DEGs.csv")


