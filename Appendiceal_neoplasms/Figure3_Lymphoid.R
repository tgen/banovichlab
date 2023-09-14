# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/06/14
# Description: Contain codes used to generate:
# Figure 3 - Lymphoid cells
# Supplementary Figure 6b
# Supplementary Figure 7
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
# Prepare the Seurat object for data analyses and plotting
# ==============================================================================
# Read in the Seurat object
coh.combined.sct <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")

# Subset out T cells and B cells
lymphoid_ob <- subset(coh.combined.sct,
                      subset = Celltype1 %in% c("B cells","T cells"))

# Adjust the SCT assays for FindMarkers function
lymphoid_ob <- PrepSCTFindMarkers(lymphoid_ob)

# Save this object (PrepSCTFindMarkers takes quite some time so it's better to save the object)
saveRDS(lymphoid_ob, file = "Appendiceal_Lymphoid_SCTprep.rds")

# ==============================================================================
# Figure 3a - Umap plot CT2
# ==============================================================================
# Rerun PCA and UMAP to group cells
lymphoid_ob <- RunPCA(lymphoid_ob)
lymphoid_ob <- RunUMAP(lymphoid_ob, dims = 1:15)

# Dimplot
lymphoid_cols <- c("CD4+ T cells" = "#FF8F8A", "CD8+ T cells" = "#CE9C48",
                   "Tregs" = "#F4B990", "NK cells" = "#57B3DE","GALTB" = "#56AF57",
                   "Follicular B" = "#ABB966", "GCBcell" = "#116b4f", 
                   "Plasma B" = "#93abd5", "Proliferating T cells" = "#BD77B2")
DimPlot(lymphoid_ob, group.by = "Celltype2", label = TRUE, repel = TRUE,
        cols = lymphoid_cols) + 
  NoLegend()

# ==============================================================================
# Figure 3b: Lymphoid cell proportion
# ==============================================================================
library(speckle)

propeller_test <- propeller(clusters = lymphoid_ob$Celltype2, 
                            sample = lymphoid_ob$orig.ident, 
                            group = lymphoid_ob$Pathology)

write.csv(propeller_test, file = "Lymphoid_CT2_propellertest.csv")

# Plot cell type proportions
plotCellTypeProps(clusters=lymphoid_ob$Celltype2, sample=lymphoid_ob$Pathology2) 
  coord_flip()

plotCellTypeProps(clusters=lymphoid_ob$Pathology2, sample=lymphoid_ob$Celltype2) +
  coord_flip()

# ==============================================================================
# Figure 3c-d: T cell DEGs cancerous vs. control
# ==============================================================================
# Subset out T-cells
t_cells <- subset(lymphoid_ob, subset = Celltype1 == "T cells")
t_cells <- PrepSCTFindMarkers(t_cells)

# DEG between cancerous vs. control for each CT
app_list = list()
j=0
for(i in unique(t_cells@meta.data$Celltype2)){
  j=j+1
  app_list[[j]] <- subset(t_cells, 
                          subset = Celltype2 == i)
}
for(i in 1:length(app_list)){
  names(app_list) <- lapply(app_list, function(xx){paste(unique(xx@meta.data$Celltype2))})
}

# Perform PrepSCTFindMarkers to ensure all fixed values are set properly
app_list <- lapply(app_list, function(x) PrepSCTFindMarkers(x))

# Run DEG using the SCT pearson residuals
cancerous_vs_control <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype2))
  if(length(unique(xx@meta.data$Status)) > 1) {
    FindMarkers(xx, 
                group.by = "Status", 
                ident.1 = "Cancerous", 
                ident.2 = "Control",
                test.use = "negbinom",
                logfc.threshold = 0,
                assay = "SCT",
                latent.vars = "Flowcell")
  } 
  else{
    return(NULL)
  } 
})
saveRDS(cancerous_vs_control, file = "Appendiceal_Tcell_cancerous_vs_control_DEGs.rds")
# Write csv files
for(i in 1:length(cancerous_vs_control)){
  write.table(cancerous_vs_control[[i]], 
              paste(gsub("/", "", unique(app_list[[i]]@meta.data$Celltype2)), 
                    "_cancerous_vs_control", ".csv"), sep =",", quote = F)
}

# Figure 3c: DEG cancerous vs. control upset plot
onion <- lapply(cancerous_vs_control, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(cancerous_vs_control, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(cancerous_vs_control)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

upset(upset_d_vs_c, nsets = 5, text.scale = 2, show.numbers = T, 
      keep.order = T,
      sets = c("CD4+ T cells", "Tregs","CD8+ T cells", "NK cells","Proliferating T cells")) 

# Figure 3d: violin plots for genes upregulated in the cancerous samples
heatmap_genes <- rownames(upset_d_vs_c[rowSums(upset_d_vs_c) > 3,]) #genes present in at least 3 comparisons
ribo <- grep("^RP", heatmap_genes, value = TRUE)
mito <- grep("^MT", heatmap_genes, value = TRUE)
heatmap_genes <- heatmap_genes[!heatmap_genes %in% c(ribo,mito)]

onion <- lapply(cancerous_vs_control, function(xx){ row.names(xx[order(-xx$avg_log2FC),][1:20,])})
onion <- unique(unlist(onion))

pdf("Tcell_Cancerous_control_vlnplot.pdf")
plot=list()
k=0
for(i in intersect(onion, heatmap_genes)){ # vlnplot for genes in at least 3 comparison, high logFC
  k=k+1
  plot[[k]] <- VlnPlot(t_cells, features = i, group.by = "Celltype2", 
          split.by = "Status", pt.size = 0, 
          split.plot = T, assay = "SCT")
  }
plot
dev.off()

# ==============================================================================
# Figure 3e: T cell DEGs pathology vs. control
# ==============================================================================
app_list = list()
j=0
for(i in unique(t_cells@meta.data$Celltype2)){
  j=j+1
  app_list[[j]] <- subset(t_cells, 
                          subset = Celltype2 == i)
}
for(i in 1:length(app_list)){
  names(app_list) <- lapply(app_list, function(xx){paste(unique(xx@meta.data$Celltype2))})
}

# Perform PrepSCTFindMarkers to ensure all fixed values are set properly
app_list <- lapply(app_list, function(x) PrepSCTFindMarkers(x))

# Run DEG using the SCT pearson residuals
lamn_vs_control <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype2))
  if(length(unique(xx@meta.data$Pathology2)) > 1) {
    FindMarkers(xx, 
                group.by = "Pathology2", 
                ident.1 = "LAMN", 
                ident.2 = "Normal",
                test.use = "negbinom",
                logfc.threshold = 0,
                assay = "SCT",
                latent.vars = "Flowcell",
                min.cells.group = 1)
  } 
  else{
    return(NULL)
  } 
})
saveRDS(lamn_vs_control, file = "Appendiceal_Tcell_LAMN_Control_DEG.rds")

for(i in 1:length(lamn_vs_control)){
  write.table(lamn_vs_control[[i]], 
              paste(gsub("/", "", unique(app_list[[i]]@meta.data$Celltype2)), 
                    "_LAMN_vs_control", ".csv"), sep =",", quote = F)
}

mhna_vs_control <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype2))
  if(length(unique(xx@meta.data$Pathology2)) > 1) {
    FindMarkers(xx, 
                group.by = "Pathology2", 
                ident.1 = "MHNA", 
                ident.2 = "Normal",
                test.use = "negbinom",
                logfc.threshold = 0,
                assay = "SCT",
                latent.vars = "Flowcell")
  } 
  else{
    return(NULL)
  } 
})
saveRDS(mhna_vs_control, file = "Appendiceal_Tcell_MHNA_Control_DEG.rds")

for(i in 1:length(mhna_vs_control)){
  write.table(mhna_vs_control[[i]], 
              paste(gsub("/", "", unique(app_list[[i]]@meta.data$Celltype2)), 
                    "_MHNA_vs_control", ".csv"), sep =",", quote = F)
}

gca_vs_control <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype2))
  if(length(unique(xx@meta.data$Pathology2)) > 1) {
    FindMarkers(xx, 
                group.by = "Pathology2", 
                ident.1 = "GCA", 
                ident.2 = "Normal",
                test.use = "negbinom",
                logfc.threshold = 0,
                assay = "SCT",
                latent.vars = "Flowcell")
  } 
  else{
    return(NULL)
  } 
})
saveRDS(gca_vs_control, file = "Appendiceal_Tcell_GCA_Control_DEG.rds")

for(i in 1:length(gca_vs_control)){
  write.table(gca_vs_control[[i]], 
              paste(gsub("/", "", unique(app_list[[i]]@meta.data$Celltype2)), 
                    "_GCA_vs_control", ".csv"), sep =",", quote = F)
}

lgma_vs_control <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype2))
  if(length(unique(xx@meta.data$Pathology2)) > 1) {
    FindMarkers(xx, 
                group.by = "Pathology2", 
                ident.1 = "LGMA", 
                ident.2 = "Normal",
                test.use = "negbinom",
                logfc.threshold = 0,
                assay = "SCT",
                latent.vars = "Flowcell")
  } 
  else{
    return(NULL)
  } 
})
saveRDS(lgma_vs_control, file = "Appendiceal_Tcell_LGMA_Control_DEG.rds")

for(i in 1:length(lgma_vs_control)){
  write.table(lgma_vs_control[[i]], 
              paste(gsub("/", "", unique(app_list[[i]]@meta.data$Celltype2)), 
                    "_LGMA_vs_control", ".csv"), sep =",", quote = F)
}

# Fig3e: Make a bar plot for number of DEGs in all CT
DE_table1 <- NULL
for (i in 1:length(lamn_vs_control)){
  DE_table1 <- rbind(DE_table1, nrow(lamn_vs_control[[i]][lamn_vs_control[[i]]$p_val_adj <= .1,]))
}
row.names(DE_table1) <- names(lamn_vs_control)

DE_table2 <- NULL
for (i in 1:length(mhna_vs_control)){
  DE_table2 <- rbind(DE_table2, nrow(mhna_vs_control[[i]][mhna_vs_control[[i]]$p_val_adj <= .1,]))
}
row.names(DE_table2) <- names(mhna_vs_control)

DE_table3 <- NULL
for (i in 1:length(gca_vs_control)){
  DE_table3 <- rbind(DE_table3, nrow(gca_vs_control[[i]][gca_vs_control[[i]]$p_val_adj <= .1,]))
}
row.names(DE_table3) <- names(gca_vs_control)

DE_table4 <- NULL
for (i in 1:length(lgma_vs_control)){
  DE_table4 <- rbind(DE_table4, nrow(lgma_vs_control[[i]][lgma_vs_control[[i]]$p_val_adj <= .1,]))
}
row.names(DE_table4) <- names(lgma_vs_control)

DE_table <- cbind(DE_table1, DE_table2, DE_table3, DE_table4)
colnames(DE_table) <- c("LAMN", "MHNA","GCA","LGMA")
DE_table <- melt(DE_table)
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
ggplot(DE_table, aes(x=Var1, y=value, fill=Var2)) +
  geom_col(position="dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x="Cell Type",
       y = "Number of differentially expressed genes")  +
  scale_fill_manual(values=pathology_col)

# ---------------------------------------
# Supplementary Fig 7: Make an Upset plot for shared of those 4 comparisons
# ---------------------------------------
deg_list <- list(lamn_vs_control[[5]], mhna_vs_control[[5]], gca_vs_control[[5]], 
                 lgma_vs_control[[5]]) #change object # in the list for different CT
names(deg_list) <- c("LAMN", "MHNA","GCA","LGMA")
onion <- lapply(deg_list, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(deg_list, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(deg_list)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

upset(upset_d_vs_c, nsets = 4, text.scale = 2, show.numbers = T, 
      keep.order = T,
      sets = c("LAMN", "MHNA","GCA","LGMA")) 

# ---------------------------------------
# Supplementary Figure 6b
# ---------------------------------------
deg_cd4 <- list(lamn_vs_control[[1]], mhna_vs_control[[1]], gca_vs_control[[1]], 
                lgma_vs_control[[1]])
deg_treg <- list(lamn_vs_control[[2]], mhna_vs_control[[2]], gca_vs_control[[2]], 
                 lgma_vs_control[[2]])
deg_cd8 <- list(lamn_vs_control[[3]], mhna_vs_control[[3]], gca_vs_control[[3]], 
                lgma_vs_control[[3]])
deg_nk <- list(lamn_vs_control[[4]], mhna_vs_control[[4]], gca_vs_control[[4]], 
               lgma_vs_control[[4]])
onion1 <- lapply(deg_cd4, function(xx){ row.names(xx[order(-xx$avg_log2FC),][1:10,])})
onion2 <- lapply(deg_treg, function(xx){ row.names(xx[order(-xx$avg_log2FC),][1:10,])})
onion3 <- lapply(deg_cd8, function(xx){ row.names(xx[order(-xx$avg_log2FC),][1:10,])})
onion4 <- lapply(deg_nk, function(xx){ row.names(xx[order(-xx$avg_log2FC),][1:10,])})

onion <- unique(c(unlist(onion1), unlist(onion2), unlist(onion3), unlist(onion4)),
                fromLast = TRUE)
ribo <- grep("^RP", onion, value = TRUE)
mito <- grep("^MT", onion, value = TRUE)
onion <- onion[!onion %in% c(ribo,mito)]

celltype_col <- c("CD4+ T cells" = "#FF8F8A", "CD8+ T cells" = "#CE9C48",
                  "Tregs" = "#F4B990", "NK cells" = "#57B3DE", 
                  "Proliferating T cells" = "#BD77B2")
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")

col = list(Celltype = celltype_col,
           Pathology = pathology_col)

tcell_ave <- AverageExpression(t_cells, group.by = c("Celltype2","Pathology2"), return.seurat = TRUE,
                               assays = "SCT")
tcell_ave@meta.data$Pathology <- sapply(strsplit(rownames(tcell_ave@meta.data),"_"), `[`, 2)
tcell_ave@meta.data$Celltype <- tcell_ave@meta.data$orig.ident

tcell_data <- tcell_ave@assays$SCT@data[rownames(tcell_ave@assays$SCT@data) %in% 
                                          onion,]

## Create metadata for heatmap col annotation
tcell_cells_metadata <- data.frame("cell" = rownames(tcell_ave@meta.data),
                                   "Celltype" =  as.character(tcell_ave$Celltype),
                                   "Pathology" = as.character(tcell_ave$Pathology))
tcell_cells_metadata <- tcell_cells_metadata[order(tcell_cells_metadata$Celltype,
                                                   tcell_cells_metadata$Pathology),]
tcell_data <- tcell_data[,tcell_cells_metadata$cell]

## Set up col and row annotation
col_ha <- HeatmapAnnotation(
  df = data.frame(Celltype = tcell_cells_metadata$Celltype,
                  Pathology = tcell_cells_metadata$Pathology),
  annotation_height = unit(4, "mm"),
  col = col
)

# Heatmap colors
col_fun = colorRamp2(c(-2, 0, 2), c("cadetblue4", "white", "coral2"))

tcell_data2 <- t(scale(t(tcell_data)))
tcell_data2 <- na.omit(tcell_data2)
Heatmap(as.matrix(tcell_data2), 
        name = "Scaled_Exp", 
        column_split = tcell_cells_metadata$Celltype,
        column_title = "Cells", row_title = "Gene expression",
        col = col_fun,
        use_raster = FALSE,
        show_column_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 5), # Text size for row names
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8))
) 
# ==============================================================================
# Figure 3f - cytotoxic/exhaustion marker heatmap
# ==============================================================================
# Read in the marker file
marker_all <- read.csv("/scratch/lbui/Appendiceal_data/tcell_markers.csv",
                       header = T)

# Split to have T cell markers as a separate heatmap
tcell_markers <- marker_all[marker_all$type == "T cells",]
marker_all <- marker_all[marker_all$type != "T cells",]

# Define colors for each level of categorical variables
patient_col <- nord("aurora", 16)
names(patient_col) <- as.character(unique(t_cells@meta.data$orig.ident))

celltype_col <- c("CD4+ T cells" = "#F3766E", "CD8+ T cells" = "#D79128",
                  "NK cells" = "#38ACE2", "Tregs" = "#F4B990",
                  "Proliferating T cells" = "#BD77B2")
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")

marker_col <- c("Activation" = nord("afternoon_prarie", 8)[1],
                "Cytotoxic" = nord("afternoon_prarie", 8)[2],
                "Exhaustion" = nord("afternoon_prarie", 8)[3],
                "Memory" = nord("afternoon_prarie", 8)[4],
                "Memory/activation" = nord("afternoon_prarie", 8)[5],
                "Proliferation" = nord("afternoon_prarie", 8)[6],
                "NK" = nord("afternoon_prarie", 8)[7],
                "Regulatory" = nord("afternoon_prarie", 8)[8])

col = list(Celltype = celltype_col,
           Patient = patient_col,
           Pathology = pathology_col,
           Marker = marker_col)

# Heatmap colors
col_fun = colorRamp2(c(-1, 0, 2), c("cadetblue4", "white", "coral2"))

# Heatmap with average expression
# Get cluster average expression
tcell_ave <- AverageExpression(t_cells, slot = "scale.data", assay = "RNA",
                               group.by = c("Celltype2","orig.ident","Pathology2"),
                               return.seurat = T)

tcell_ave@meta.data$Patient <- sapply(strsplit(rownames(tcell_ave@meta.data),"_"), `[`, 2)
tcell_ave@meta.data$Pathology <- sapply(strsplit(rownames(tcell_ave@meta.data),"_"), `[`, 3)
tcell_data <- tcell_ave@assays$RNA@scale.data[rownames(tcell_ave@assays$RNA@scale.data) %in% 
                                                unique(marker_all$gene),]
tcell_data <- tcell_data[which(rownames(tcell_data) %in% unique(marker_all$gene)),]
all_markers_sct <- marker_all[which(marker_all$gene %in% rownames(tcell_data)),]
all_markers_sct <- all_markers_sct[!duplicated(all_markers_sct$gene),]

tcell_data <- tcell_data[order(match(rownames(tcell_data), all_markers_sct$gene)), ,drop = FALSE]

# Create metadata for heatmap col annotation
tcell_cells_metadata <- data.frame("cell" = rownames(tcell_ave@meta.data),
                                   "Celltype" =  as.character(tcell_ave$orig.ident),
                                   "Patient" = as.character(tcell_ave$Patient),
                                   "Pathology" = as.character(tcell_ave$Pathology))
tcell_cells_metadata <- tcell_cells_metadata[order(tcell_cells_metadata$Celltype,
                                                   tcell_cells_metadata$Pathology,
                                                   tcell_cells_metadata$Patient),]
tcell_data <- tcell_data[,tcell_cells_metadata$cell]

# Make heatmap
col_ha <- HeatmapAnnotation(
  df = data.frame(Celltype = tcell_cells_metadata$Celltype,
                  Pathology = tcell_cells_metadata$Pathology,
                  Patient = tcell_cells_metadata$Patient),
  annotation_height = unit(4, "mm"),
  col = col
)
row_ha <- rowAnnotation(
  df = data.frame(marker_type = all_markers_sct$type),
  annotation_height = unit(4, "mm"),
  col = col
)

Heatmap(as.matrix(tcell_data), 
        name = "Scaled_Exp", 
        column_split = tcell_cells_metadata$Celltype,
        row_split = all_markers_sct$type,
        column_title = "Cells", row_title = "Gene expression",
        col = col_fun,
        use_raster = FALSE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = col_ha,
        right_annotation = row_ha,
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8))
) 

# ==============================================================================
# Figure 3g: Follicular B DEGs
# ==============================================================================
folB_ob <- subset(lymphoid_ob, subset = Celltype2 == "Follicular B")
Idents(folB_ob) <- as.character(folB_ob$Pathology2)
folB_ob <- PrepSCTFindMarkers(folB_ob)
pathology_group <- c("LAMN","LGMA", "MHNA", "GCA")
folB_deg = list()
j=0
for(i in pathology_group){
  j=j+1
  folB_deg[[j]] <- FindMarkers(folB_ob,
                               group.by = "Pathology2",
                               ident.1 = i,
                               ident.2 = "Normal",
                               test.use = "negbinom",
                               assay = "SCT",
                               latent.vars = "Flowcell")
}

names(folB_deg) <- pathology_group
saveRDS(folB_deg, file = "Appendiceal_FollicularB_pathologyvscontrol.rds")

# Save DE files as csv
for(i in 1:length(folB_deg)){
  write.table(folB_deg[[i]], 
              paste(gsub("/", "", names(folB_deg[i])), 
                    "_vs_control_DEGs_FollicularB", ".csv"), sep =",", quote = F)
}

# Make an Upset plot for shared of those 4 comparisons (Supplementary Figure 7)
onion <- lapply(folB_deg, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(folB_deg, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(folB_deg)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

upset(upset_d_vs_c, nsets = 4, text.scale = 2, show.numbers = T, 
      keep.order = T,
      sets = c("LAMN","LGMA", "MHNA", "GCA")) 

# Select the top 20 DEGs per pathology for heatmap
onion <- lapply(folB_deg, function(xx){ row.names(xx[order(-xx$avg_log2FC),][1:20,])})
onion <- unique(unlist(onion))
ribo <- grep("^RP", onion, value = TRUE)
mito <- grep("^MT", onion, value = TRUE)
heatmap_genes <- onion[!onion %in% c(ribo,mito)]

DoHeatmap(folB_ob, features = unique(heatmap_genes), group.by = "Pathology2", 
          draw.lines = TRUE, raster = TRUE) 

# Pick some genes for vlnplot
Stacked_VlnPlot(folB_ob, features = c("KMT2E","MARCH1","CD69","IGLC2","CD83","CDK14"),
                group.by = "Pathology2", assay = "SCT", pt.size = 0,
                colors_use = pathology_col)

