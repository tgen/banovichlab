# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/06/14
# Description: Supplementary Figure 5 - Epithelial cells
# ==============================================================================

# =======================================
# SET UP THE ENVIRONMENT VARIABLES 
# =======================================
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
# Supplementary Figure 5
# ==============================================================================
# Figure S5a - heatmap for different pathologies vs control
Idents(epi) <- as.character(epi$Pathology2)
path_deg <- FindAllMarkers(epi, test.use = "negbinom",
                           assay = "RNA")
write.csv(path_deg, file="Epithelial_CT2_allmarkers_RNA.csv")
hm_genes <- path_deg %>% group_by(cluster) %>% top_n(20, avg_log2FC)

patient_col <- nord("aurora", 16)
names(patient_col) <- as.character(unique(epi@meta.data$orig.ident))

pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
col = list(Pathology = pathology_col,
           Patient = patient_col)

epi_all_markers <- GetAssayData(epi,assay = "integrated", slot = "scale.data")
epi_hm_markers <- epi_all_markers[which(rownames(epi_all_markers) %in% unique(hm_genes$gene)),]

## Create metadata for heatmap col annotation
epi_cells_metadata <- data.frame("cell" = rownames(epi@meta.data),
                                 "Patient" =  as.character(epi$orig.ident),
                                 "Pathology" = as.character(epi$Pathology2))
epi_cells_metadata <- epi_cells_metadata[order(epi_cells_metadata$Pathology,
                                               epi_cells_metadata$Patient),]
epi_hm_markers <- epi_hm_markers[,epi_cells_metadata$cell]

## Set up col and row annotation
col_ha <- HeatmapAnnotation(
  df = data.frame(Pathology = epi_cells_metadata$Pathology,
                  Patient = epi_cells_metadata$Patient),
  annotation_height = unit(4, "mm"),
  col = col
)

# Heatmap colors
col_fun = colorRamp2(c(-1, 0, 2), c("cadetblue4", "white", "coral2"))

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
        row_names_gp = gpar(fontsize = 5), # Text size for row names
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8))
) 

# Perform GO enrichment for each pathology (Figure S5b)
epi_genes <- split(path_deg, f=path_deg$cluster)
epi_genes <- lapply(epi_genes, function(x) x[x$p_val_adj <= 0.1 & x$avg_log2FC >= 1,]$gene)
ribo <- lapply(epi_genes, function(x) grep("^RP", x, value = TRUE))
ribo <- unlist(ribo)
epi_genes <- lapply(epi_genes, function(x) x[!x %in% ribo]) #remove ribosomal genes
epi_genes <- lapply(epi_genes, function(x){
  bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
})

epi_genesentrez <- lapply(epi_genes, function(x) x$ENTREZID)
str(epi_genesentrez)

compared_go <- compareCluster(geneCluster = epi_genesentrez, 
                              fun = "enrichGO", OrgDb = org.Hs.eg.db, ont="BP",
                              pAdjustMethod = "BH")
dotplot(compared_go, showCategory=5, includeAll=TRUE, font.size=10,
        color="qvalue") 
compared_go <- setReadable(compared_go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(compared_go)
write.csv(compared_go, file = "Epithelial_pathology_clusterprofiler_BP_out.csv")

compared_pw <- compareCluster(geneCluster = epi_genesentrez, 
                              fun = "enrichPathway", 
                              pAdjustMethod = "BH")
dotplot(compared_pw, showCategory=5,includeAll=TRUE, font.size=8,
        color="qvalue") 

# Figure S5c - DEGs LAMN vs LGMA in goblet-like cells
# Subset out the SPINK4+ cells
gob_sub <- subset(epi, subset = Celltype2 == "Goblet-like cells")

# Run FindMarkers
gob_deg <- FindMarkers(gob_sub,
                       group.by = "Pathology2",
                       ident.1 = "LAMN",
                       ident.2 = "LGMA",
                       test.use = "negbinom",
                       assay = "RNA",
                       slot = "counts",
                       logfc.threshold = 0)

# Make Volcano plots
gob_deg$geneID <- rownames(gob_deg)
ribo <- grep("^RP", gob_deg$geneID, value = TRUE)
gob_deg$diffexpressed <- "NO"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 
gob_deg$diffexpressed[gob_deg$avg_log2FC > 2 & gob_deg$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.05, set as "DOWN"
gob_deg$diffexpressed[gob_deg$avg_log2FC < -1 & gob_deg$p_val_adj < 0.05] <- "DOWN"
gob_deg$diffexpressed[gob_deg$geneID %in% ribo] <- "NO"
gob_deg$delabel <- NA
gob_deg$delabel[gob_deg$diffexpressed != "NO"] <- gob_deg$geneID[gob_deg$diffexpressed != "NO"]

ggplot(gob_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                    col = diffexpressed, label = delabel)) + 
  geom_point(size = 3) + 
  geom_text_repel(max.overlaps = 20) + 
  theme_bw() + 
  scale_color_manual(values = c("blue","grey80", "red")) + 
  geom_vline(xintercept=c(-1, 2), col="black", linetype="dashed")  

write.csv(gob_deg, file = "Epithelial_Goblet_like_LAMN_LGMA_DEGs.csv")
dim(gob_deg[gob_deg$avg_log2FC > 1 & gob_deg$p_val_adj < 0.05,]) # 600 8
dim(gob_deg[gob_deg$avg_log2FC < -1 & gob_deg$p_val_adj < 0.05,]) # 6 8
