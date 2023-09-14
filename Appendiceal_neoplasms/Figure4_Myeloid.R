# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/03
# Description: Myeloid cell figure generation
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

# Load libraries
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

# Rerun PCA and UMAP to recluster myeloid cells (same as the subset + annotation CT2)
# this is only for visualization purposes
myeloid$Celltype2 <- as.character(myeloid$Celltype2)
DefaultAssay(myeloid) <- "RNA"
myeloid <- FindVariableFeatures(myeloid, nfeatures = 3000, selection.method = "vst")
DefaultAssay(myeloid) <- "integrated"
myeloid <- ScaleData(myeloid, vars.to.regress = "percent.mt")
myeloid <- RunPCA(myeloid) 
myeloid <- RunUMAP(myeloid, dims = 1:14)

# ==============================================================================
# Figure 4a
# ==============================================================================
myeloid_cols <- c("C1Qhi monocytes" = "#ef862f","cDC1" = "#c19d25", "cDC2" = "#ecab2c",
              "Macrophages" = "#00bfc3", "Monocyte-like" = "#00b1f1", 
              "pDCs" = "#cb7dfb", "SPP1+ macrophages" = "#ff64b7")
DimPlot(myeloid, group.by = "Celltype2", label = T, repel = T, cols = myeloid_cols) + 
  NoLegend()

# ==============================================================================
# Figure 4b
# ==============================================================================
# Dotplot for cell marker gene expression
source("/home/lbui/SC_scripts/MMRF_CRISPR/CART_plot_functions.R") #code from Heini
DefaultAssay(myeloid) <- "SCT"
create_dotplot_heatmap_horizontal(myeloid, 
                                  plot_features = c("CD14","FCGR3A","FCN1","S100A8",
                                                    "S100A9","LILRA4","IL3RA","CD1C",
                                                    "FCER1A","CLEC9A","BATF3","NLRP3",
                                                    "APOE","LYVE1","PLTP","IL1B",
                                                    "C1QA","C1QB","C1QC","CD163",
                                                    "CCL4","SPP1","CD68"),
                                  group_var = "Celltype2",
                                  group_colors = myeloid_cols,
                                  column_title = "Myeloid markers")

# ==============================================================================
# Figure 4c
# ==============================================================================
library(speckle)

propeller_test <- propeller(clusters = myeloid$Celltype2, 
                            sample = myeloid$orig.ident, 
                            group = myeloid$Pathology)

write.csv(propeller_test, file = "Myeloid_CT2_propellertest.csv")

# Plot cell type proportions
myeloid$Celltype2 <- as.character(myeloid$Celltype2)
plotCellTypeProps(clusters=myeloid$Celltype2, sample=myeloid$Pathology2) 

# ==============================================================================
# Figure 4d-4f - C1Qhi monocytes
# ==============================================================================
# Subset out C1Qhi monocytes
c1q_mono <- subset(myeloid, subset = Celltype2 == "C1Qhi monocytes")
c1q_mono <- PrepSCTFindMarkers(c1q_mono)

# Make vlnplot (Fig 4d)
Stacked_VlnPlot(c1q_mono, features = c("HLA-DRA","HLA-DPA1","HLA-DRB1","SPP1",
                                       "FN1","NFKB1"),
                group.by = "Pathology2", colors_use = pathology_col, assay="SCT")

# Run FindMarkers
pathology_test <- c("LGMA","MHNA","LAMN","GCA")
c1q_deg <- list()
j=0
for(i in unique(pathology_test)){
  j=j+1
  print(i)
  c1q_deg[[j]] <- FindMarkers(c1q_mono,
                              ident.1 = i,
                              ident.2 = "Normal",
                              group.by = "Pathology2",
                              test.use = "negbinom",
                              assay = "SCT")
}
names(c1q_deg) <- pathology_test

# Save DE files as csv
for(i in 1:length(c1q_deg)){
  write.table(c1q_deg[[i]], 
              paste(gsub("/", "", names(c1q_deg[i])), 
                    "_vs_control_DEGs_C1Qhimono", ".csv"), sep =",", quote = F)
}

# Make upset plot
onion <- lapply(c1q_deg, function(xx){row.names(xx[xx$p_val_adj <= .1 & abs(xx$avg_log2FC) >= 1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(c1q_deg, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1 & abs(xx$avg_log2FC) >= 1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(c1q_deg)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

## Figure 4e
upset(upset_d_vs_c, nsets = 4, text.scale = 2, show.numbers = T, 
      keep.order = T,
      sets = c("LAMN", "LGMA","MHNA","GCA")) 

# Perform GO enrichment for each CT in the lineage (Fig 4f)
library(clusterProfiler)
library(org.Hs.eg.db)
c1q_genes <- lapply(c1q_deg, function(x) rownames(x[x$p_val_adj <= 0.1 & x$avg_log2FC >= 1,]))
c1q_ribo <- lapply(c1q_genes, function(x) grep("^RP", x, value = TRUE))
c1q_ribo <- unlist(c1q_ribo)
c1q_genes <- lapply(c1q_genes, function(x) x[!x %in% c1q_ribo]) #remove ribosomal genes
c1q_genes <- lapply(c1q_genes, function(x){
  bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
})

c1q_genesentrez <- lapply(c1q_genes, function(x) x$ENTREZID)
str(c1q_genesentrez)

compared_go <- compareCluster(geneCluster = c1q_genesentrez, 
                              fun = "enrichGO", OrgDb = org.Hs.eg.db, ont="BP",
                              pAdjustMethod = "BH")
dotplot(compared_go, showCategory=5, includeAll=FALSE, font.size=10,
        color="qvalue") 

compared_pw <- compareCluster(geneCluster = c1q_genesentrez, 
                              fun = "enrichPathway", 
                              pAdjustMethod = "BH")
dotplot(compared_pw, showCategory=4,includeAll=FALSE, font.size=10,
        color="qvalue") 

# ==============================================================================
# Figure 4g - Vlnplot macrophages collagen genes
# ==============================================================================
macro <- subset(myeloid, subset = Celltype2 %in% c("Macrophages"))
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
Stacked_VlnPlot(macro, features = c("COL1A1","COL1A2","COL3A1"), 
                group.by = "Pathology2", assay = "SCT", colors_use = pathology_col)

# ==============================================================================
# Figure 4h - Heatmap for other CT2 different pathologies
# ==============================================================================
# Only run for CT2 with >40 cells per pathologies
# No normal cells in these comparisons
subset_ob <- subset(myeloid, subset = Celltype2 %in% 
                      c("cDC2","Macrophages","Monocyte-like",'pDCs','SPP1+ TAMs'))
subset_ob <- subset(subset_ob, subset = Pathology2 %in% c("LAMN","LGMA","MHNA"))

app_list = list()
j=0
for(i in unique(subset_ob@meta.data$Celltype2)){
  j=j+1
  app_list[[j]] <- subset(subset_ob, 
                          subset = Celltype2 == i)
}
for(i in 1:length(app_list)){
  names(app_list) <- lapply(app_list, function(xx){paste(unique(xx@meta.data$Celltype2))})
}

# Perform PrepSCTFindMarkers to ensure all fixed values are set properly
app_list <- lapply(app_list, function(x) PrepSCTFindMarkers(x))

# Run DEG using the SCT pearson residuals
subset_deg <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype2))
  Idents(xx) <- as.character(xx$Pathology2)
  FindAllMarkers(xx, 
                test.use = "negbinom",
                logfc.threshold = 0,
                assay = "SCT",
                latent.vars = "Flowcell")
})

# Write the DEGs into csv file
subset_deg <- mapply(cbind, subset_deg, "celltype"=names(subset_deg), SIMPLIFY=F)
subset_deg_unlist <- do.call(rbind, subset_deg)
write.csv(subset_deg_unlist, file = "Myeloid_CT2_LAMN_LGMA_MHNA_DEG_all_SCT.csv",
          row.names = F)

# Make a heatmap with complexheatmap
onion <- lapply(subset_deg, function(xx) xx %>% group_by(cluster) %>% top_n(10, avg_log2FC))
onion <- lapply(onion, function(xx) xx$gene)
onion <- unlist(onion)
ribo <- grep("^RP", onion, value = TRUE)
mito <- grep("^MT", onion, value = TRUE)
onion <- onion[!onion %in% c(ribo,mito)]

pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
col = list(Celltype = myeloid_cols,
           Pathology = pathology_col)

subset_ave <- AverageExpression(subset_ob, group.by = c("Celltype2","Pathology2"), 
                                return.seurat = TRUE, assays = "integrated")
subset_ave@meta.data$Pathology <- sapply(strsplit(rownames(subset_ave@meta.data),"_"), `[`, 2)
subset_ave@meta.data$Celltype <- subset_ave@meta.data$orig.ident

subset_data <- subset_ave@assays$integrated@scale.data[rownames(subset_ave@assays$integrated@scale.data) %in% 
                                          unique(onion),]

## Create metadata for heatmap col annotation
subset_cells_metadata <- data.frame("cell" = rownames(subset_ave@meta.data),
                                   "Celltype" =  as.character(subset_ave$Celltype),
                                   "Pathology" = as.character(subset_ave$Pathology))
subset_cells_metadata <- subset_cells_metadata[order(subset_cells_metadata$Celltype,
                                                     subset_cells_metadata$Pathology),]
subset_data <- subset_data[,subset_cells_metadata$cell]

## Set up col and row annotation
col_ha <- HeatmapAnnotation(
  df = data.frame(Celltype = subset_cells_metadata$Celltype,
                  Pathology = subset_cells_metadata$Pathology),
  annotation_height = unit(4, "mm"),
  col = col
)

# Heatmap colors
col_fun = colorRamp2(c(-2, 0, 2), c("cadetblue4", "white", "coral2"))

#subset_data2 <- t(scale(t(subset_data)))
Heatmap(as.matrix(subset_data), 
        name = "Scaled_Exp", 
        column_split = subset_cells_metadata$Celltype,
        column_title = "Cells", row_title = "Gene expression",
        col = col_fun,
        use_raster = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
       # km = 5,
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8))
) 

