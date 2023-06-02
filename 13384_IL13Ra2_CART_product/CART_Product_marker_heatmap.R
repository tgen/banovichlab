#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 02/12/2022
# Description: CAR T product marker heatmap
#==============================================================================

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(qvalue)
library(ComplexHeatmap)
library(nord)
library(circlize)
library(ggmin)
library(gridExtra)
library(ggpubr)
library(DescTools)
library(qdap)
library(plyr)
library(viridis)
library(colourvalues)

#==============================================================================
# Helper functions
#==============================================================================

source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/CART/CART_colors_themes.R")
source("/home/hnatri/CART/CART_plot_functions.R")

product_cluster_col_c <- product_cluster_col_woutsmall
names(product_cluster_col_c) <- paste0("C", names(product_cluster_col_c))

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)
setwd("/scratch/hnatri/CART/")

#==============================================================================
# Import Seurat objects
#==============================================================================

product <- readRDS("/labs/banovich/BCTCSF/Heini/product_projecTILs_modules_scGSVA.rds")

DefaultAssay(product) <- "RNA"

product$Batch <- gsub("Batch", "", product$Batch)
product$Batch <- as.numeric(product$Batch)

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha = HeatmapAnnotation(
  foo = 1:10, 
  bar = sample(letters[1:3], 10, replace = TRUE),
  col = list(foo = col_fun,
             bar = c("a" = "red", "b" = "green", "c" = "blue")
  )
)

#==============================================================================
# Marker gene expression heatmaps
#==============================================================================

# Canonical T cell markers
gs4_deauth()
canonical_markers  <- gs4_get("https://docs.google.com/spreadsheets/d/186PIcU3Jb0Xm3IuZgENEChluykvwsK4JnHF8FbeJASA/edit?usp=sharing")
sheet_names(canonical_markers)
canonical_markers <- read_sheet(canonical_markers, sheet = "T cells, gene sets")
head(canonical_markers)

# Downsampling for plotting (Seurat DoHeatmap uses ggplot, doesn't work with
# >30k cells. downsample = # of cells per identity class. Same problem with
# ComplexHeatmap).

DimPlot(product, group.by = "cluster", reduction = "wnn.umap", label=T) + NoLegend()
VlnPlot(product, group.by = "cluster", features = c("nCount_RNA", "nFeature_RNA"), pt.size=0)
VlnPlot(product, group.by = "cluster", features = c("nCount_Protein", "nFeature_Protein"), pt.size=0)

# 30000
product_downsample <- subset(product, downsample = 30000/length(unique(product@meta.data$cluster)))
table(product_downsample$cluster)

DimPlot(product_downsample, group.by = "cluster", reduction = "wnn.umap", label=T) + NoLegend()

# Rescaling all features
DefaultAssay(product_downsample) <- "Protein"
VariableFeatures(product_downsample) <- rownames(product_downsample) # all protein features
product_downsample <- NormalizeData(product_downsample)
product_downsample <- ScaleData(product_downsample)
                #vars.to.regress = c("percent.mt", "S.Score", "G2M.Score", "ab_panel"))

DefaultAssay(product_downsample) <- "RNA"
VariableFeatures(product_downsample) <- rownames(product_downsample)
product_downsample <- NormalizeData(product_downsample)
product_downsample <- ScaleData(product_downsample)
                #vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

# Top markers for each cluster, using the full object
cluster_markers_rna <- presto::wilcoxauc(product,
                                         group_by = "cluster",
                                         #groups_use = "",
                                         assay = "data",
                                         seurat_assay = "RNA")
cluster_markers_prot <- presto::wilcoxauc(product,
                                          group_by = "cluster",
                                          #groups_use = "",
                                          assay = "data",
                                          seurat_assay = "Protein")

markers <- presto::top_markers(cluster_markers_rna, n = 20, auc_min = 0.5, pct_in_min = 50, pct_out_max = 100)

top_cluster_markers_rna <- markers %>%
    dplyr::select(-rank) %>% 
    unclass() %>% 
    stack() %>%
    pull(values) %>%
    unique() %>%
    .[!is.na(.)]

length(top_cluster_markers_rna)

# For each cluster, top DEGs between manufacture types
product_manufacture_wilcoxauc_rna <- lapply(unique(product$cluster), function(xx){
    message(xx)
    product_subset <- subset(product, subset = cluster == xx)
    if (length(unique(product_subset$Manufacture))<2){
        return(NULL)
    } else {
        markers <- presto::wilcoxauc(product_subset,
                                     group_by = "Manufacture",
                                     groups_use = c("TCM", "Tnmem"), # TCM vs. Tnmem
                                     assay = "data",
                                     seurat_assay = "RNA")
        markers$cluster <- xx
        
        return(markers)
    }
})
names(product_manufacture_wilcoxauc_rna) <- unique(product$cluster)
product_manufacture_wilcoxauc_rna[sapply(product_manufacture_wilcoxauc_rna, is.null)] <- NULL

product_manufacture_wilcoxauc_rna_df <- as.data.frame(do.call(rbind, product_manufacture_wilcoxauc_rna))
product_manufacture_wilcoxauc_rna_df <- product_manufacture_wilcoxauc_rna_df[which(product_manufacture_wilcoxauc_rna_df$group=="Tnmem"),]
product_manufacture_wilcoxauc_rna_df_sig <- product_manufacture_wilcoxauc_rna_df[which(product_manufacture_wilcoxauc_rna_df$pct_in > 50 | product_manufacture_wilcoxauc_rna_df$pct_out > 50),]
product_manufacture_wilcoxauc_rna_df_sig <- product_manufacture_wilcoxauc_rna_df_sig[which(product_manufacture_wilcoxauc_rna_df_sig$padj < 0.01 & product_manufacture_wilcoxauc_rna_df_sig$auc > 0.6),]

# For each cluster, top DEGs between response and non-response
product_response_wilcoxauc_rna <- lapply(unique(product$cluster), function(xx){
    message(xx)
    product_subset <- subset(product, subset = cluster == xx)
    if (length(unique(product_subset$binary_outcome))<2){
        return(NULL)
    } else {
        markers <- presto::wilcoxauc(product_subset,
                                     group_by = "binary_outcome",
                                     groups_use = c("PD", "CR/SD/PR"),
                                     assay = "data",
                                     seurat_assay = "RNA")
        markers$cluster <- xx
        
        return(markers)
    }
})
names(product_response_wilcoxauc_rna) <- unique(product$cluster)
product_response_wilcoxauc_rna[sapply(product_response_wilcoxauc_rna, is.null)] <- NULL

product_response_wilcoxauc_rna_df <- as.data.frame(do.call(rbind, product_response_wilcoxauc_rna))
product_response_wilcoxauc_rna_df <- product_response_wilcoxauc_rna_df[which(product_response_wilcoxauc_rna_df$group=="CR/SD/PR"),]
product_response_wilcoxauc_rna_df_sig <- product_response_wilcoxauc_rna_df[which(product_response_wilcoxauc_rna_df$pct_in > 50 | product_response_wilcoxauc_rna_df$pct_out > 50),]
product_response_wilcoxauc_rna_df_sig <- product_response_wilcoxauc_rna_df_sig[which(product_response_wilcoxauc_rna_df_sig$padj < 0.01 & product_response_wilcoxauc_rna_df_sig$auc > 0.6),]

# Top markers and canonical markers
all_plot_markers <- unique(c(top_cluster_markers_rna, setdiff(canonical_markers$RNA, c(NA))))
length(all_plot_markers)
# Adding DEGs between manufacture types and response groups
all_plot_markers <- unique(c(all_plot_markers, product_manufacture_wilcoxauc_rna_df_sig$feature, product_response_wilcoxauc_rna_df_sig$feature))
# Removing RP and MT (all response DEGs are RPSs)
all_plot_markers <- all_plot_markers[-grep("RP", all_plot_markers)]
all_plot_markers <- all_plot_markers[-grep("MT-", all_plot_markers)]

# A heatmap with annotations
mt_ratio_col <- colorRamp2(c(min(product_downsample$percent.mt),
                             mean(product_downsample$percent.mt),
                             max(product_downsample$percent.mt)),
                           c("darkgreen", "khaki1", "red4"))

# colourvalues::colour_values(number_vector)
memory_col <- colorRamp2(quantile(product_downsample$Memory, seq(0, 1, by = 0.25)), magma(5))
dysfunction_col <- colorRamp2(quantile(product_downsample$Dysfunction, seq(0, 1, by = 0.25)), magma(5))
cytotoxic_col <- colorRamp2(quantile(product_downsample$Cytotoxic, seq(0, 1, by = 0.25)), magma(5))

cont_col <- list(percent.mt = mt_ratio_col,
                 dysfunction = cytotoxic_col,
                 memory = cytotoxic_col,
                 cytotoxic = cytotoxic_col)

# Gene expression data from the RNA assay
cellpop_rna_all_markers <- GetAssayData(product_downsample, slot = "scale.data", assay = "RNA")
cellpop_rna_all_markers <- cellpop_rna_all_markers[which(rownames(cellpop_rna_all_markers) %in% all_plot_markers),]
dim(cellpop_rna_all_markers)

# Order by cluster, then by response
cellpop_cells_metadata <- data.frame("cell" = rownames(product_downsample@meta.data),
                                     "cluster" = product_downsample@meta.data$c_cluster,
                                     "batch" = product_downsample@meta.data$Batch,
                                     "response" = product_downsample@meta.data$response,
                                     "phase" = product_downsample@meta.data$Phase,
                                     #"manu" = cellpop_object@meta.data$Manufacture,
                                     #"ab_panel" = product_downsample@meta.data$ab_panel,
                                     #"CD4_CD8_ratio" = cellpop_object@meta.data$CD4_CD8_ratio,
                                     "nCount_RNA" = product_downsample@meta.data$nCount_RNA,
                                     "nCount_Protein" = product_downsample@meta.data$nCount_Protein,
                                     "nFeature_RNA" = product_downsample@meta.data$nFeature_RNA,
                                     "percent.mt" = product_downsample@meta.data$percent.mt,
                                     "memory" = product_downsample@meta.data$Memory,
                                     "dysfunction" = product_downsample@meta.data$Dysfunction,
                                     "cytotoxic" = product_downsample@meta.data$Cytotoxic)
cellpop_cells_metadata <- cellpop_cells_metadata[order(cellpop_cells_metadata$cluster),]
cellpop_rna_all_markers <- as.data.frame(cellpop_rna_all_markers)[cellpop_cells_metadata$cell]

# RNA heatmap
col$cluster <- product_cluster_col_c
col_ha <- HeatmapAnnotation(
  df = data.frame(cluster = cellpop_cells_metadata$cluster,
                  batch = cellpop_cells_metadata$batch,
                  response = cellpop_cells_metadata$response,
                  cell_cycle_phase = cellpop_cells_metadata$phase,
                  percent.mt = cellpop_cells_metadata$percent.mt,
                  memory = cellpop_cells_metadata$memory,
                  cytotoxic = cellpop_cells_metadata$cytotoxic,
                  dysfunction = cellpop_cells_metadata$dysfunction
                  #manufacture = cellpop_cells_metadata$manu,
                  #ab_panel = cellpop_cells_metadata$ab_panel
                  ),
  annotation_height = unit(4, "mm"),
  col = c(col, cont_col)
)

# Only annotating selected genes
memory_markers <- setdiff(canonical_markers[which(canonical_markers$Bigger_gene_sets=="Memory"),]$RNA, c(NA))
dysfunction_markers <- setdiff(canonical_markers[which(canonical_markers$Bigger_gene_sets=="Dysfunction"),]$RNA, c(NA))
cytotoxic_markers <- setdiff(canonical_markers[which(canonical_markers$Bigger_gene_sets=="Cytotoxic"),]$RNA, c(NA))
labels <- unique(c(memory_markers, dysfunction_markers, cytotoxic_markers))
label_positions <- grep(paste(labels, collapse = "|"), rownames(cellpop_rna_all_markers))

row_ha_labels <- rowAnnotation(link = anno_mark(at = label_positions, 
                               labels = labels, 
                               labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))


# Quantile-normalizing rows for plotting
cellpop_rna_all_markers_qqnorm <- t(apply(cellpop_rna_all_markers, 1, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(cellpop_rna_all_markers_qqnorm) <- colnames(cellpop_rna_all_markers)
rownames(cellpop_rna_all_markers_qqnorm) <- rownames(cellpop_rna_all_markers)

# Heatmap colors
col_fun2 <- colorRamp2(quantile(as.matrix(cellpop_rna_all_markers_qqnorm), seq(0, 1, by = 0.25)), viridis(5))

set.seed(1234)

plot_func <- function(){
  hm_rna  <- Heatmap(as.matrix(cellpop_rna_all_markers_qqnorm), 
                     name = "Gene exp", # Title of legend
                     #col = col_fun2,
                     col = viridis(100),
                     use_raster = T,
                     #column_title = "Cells", row_title = "Gene expression",
                     column_title = NULL, row_title = NULL,
                     column_split = cellpop_cells_metadata$cluster,
                     #row_split = all_markers_rna_sct$Marker_type,
                     #use_raster = TRUE,
                     show_column_names = FALSE,
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     row_km = 10,
                     cluster_columns = FALSE,
                     top_annotation = col_ha,
                     right_annotation = row_ha_labels,
                     height = nrow(cellpop_rna_all_markers)*unit(1, "mm"))
                     #bottom_annotation = bottom_ha,
                     #row_names_gp = gpar(fontsize = 7),  # Text size for row names
                     #row_names_side = "right")
      
  heatmap <- draw(hm_rna)
}

#for(i in 1:length(text_list)) {
#  decorate_annotation("marker_type_text", slice = i, {
#    #grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
#    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left", gp=gpar(fontsize=10))
#  })
#}

p <- plot_func()

# Saving to a file
filename <- "/home/hnatri/CART/product_rna_heatmap.pdf"
pdf(file = filename,
    width = 10,
    height = 12)
p
dev.off()

png(file="/home/hnatri/CART/product_rna_heatmap.png",
    width=1000, height=1200)
p
dev.off()

