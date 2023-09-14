# ==============================================================================
# Date: 2023/02/16
# Author: Linh T. Bui (lbui@tgen.org)
# Appendiceal cancer manuscript - Figure 1 codes
# ==============================================================================
# ==============================================================================
# SET UP THE ENVIRONMENT VARIABLES 
# ==============================================================================
# Set the working directory
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
library(scCustomize)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(RCurl)
library(reshape2)
library(ggrepel)
library(ComplexHeatmap)

# ==============================================================================
# Read in the Seurat object
# ==============================================================================
coh.combined.sct <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")

# ==============================================================================
# Figure 1B
# ==============================================================================
# Bar plot for number of cells per status/pathology
status_data <- as.data.frame(table(coh.combined.sct$Status))
colnames(status_data) <- c("Group","Counts")
ggplot(status_data, aes(x=Group, y= Counts, fill = Group)) +
  geom_bar(stat="identity", width=0.5) + theme_bw() +
  ylab("Number of cells") + NoLegend()

pathology_data <- as.data.frame(table(coh.combined.sct$Pathology2))
colnames(pathology_data) <- c("Group","Counts")
ggplot(pathology_data, aes(x=Group, y= Counts, fill = Group)) +
  geom_bar(stat="identity", width=0.8) + theme_bw() +
  ylab("Number of cells") + xlab("Pathology group") + NoLegend()

# ==============================================================================
# Figure 1C: heatmap for lineage specific genes
# ==============================================================================
# Run FindAllMarkers to find conserved markers per lineage
Idents(coh.combined.sct) <- coh.combined.sct$Celltype1
coh.combined.sct <- PrepSCTFindMarkers(coh.combined.sct)
all_markers <- FindAllMarkers(coh.combined.sct, assay = "SCT")
write.csv(all_markers, file = "Appendiceal_CT1_allmarkers_SCT.csv")

# Define top 5 genes for heatmap
topgenes <- all_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
heatmap_features <- c("CD3E","CD3D","IL7R",topgenes[topgenes$cluster == "T cells",]$gene,
                      "MS4A1","IGHD","CD79A",topgenes[topgenes$cluster == "B cells",]$gene,
                      "KIT",topgenes[topgenes$cluster == "Mast cells",]$gene,
                      "MUC2","FABP1","EPCAM",topgenes[topgenes$cluster == "Epithelial",]$gene,
                      "CD14","CD68","FCGR3A",topgenes[topgenes$cluster == "Myeloid cells",]$gene,
                      "FCGR3B","CXCR2",topgenes[topgenes$cluster == "FCGR3B+ neutrophils",]$gene,
                      "VWF","ACKR1","PECAM1",topgenes[topgenes$cluster == "Endothelial",]$gene,
                      "LUM","THY1","DCN",topgenes[topgenes$cluster == "Mesenchymal cells",]$gene)

# Order the cell lineage on heatmap
coh.combined.sct$Celltype1 <- factor(coh.combined.sct$Celltype1,
                                     levels = c("T cells","B cells","Mast cells","Epithelial",
                                                "Myeloid cells","FCGR3B+ neutrophils",
                                                "Endothelial","Mesenchymal cells"))

# Assign colors for cell lineages
hm_color <- c("B cells"="#d2331f","Endothelial"="#0072ae","Epithelial"="#43B425",
              "FCGR3B+ neutrophils"="#5f1a1a","Mesenchymal cells"="#00B6E8",
              "Myeloid cells"="#A58BFB","T cells"="#FF62D5", "Mast cells" = "#00785f")

# Generate the heatmap with all cells or subset
Idents(coh.combined.sct) <- "Celltype1"
coh_subset <- subset(coh.combined.sct,downsample=500)
pdf("COH_CT1_heatmap_downsample.pdf")
DoHeatmap(coh_subset, features = unique(heatmap_features), 
          group.by = "Celltype1", group.colors = hm_color) +
  scale_fill_gradientn(colors = c("cadetblue4", "white", "coral2"))
dev.off()

pdf("COH_CT1_heatmap_full.pdf")
DoHeatmap(coh.combined.sct, features = unique(heatmap_features), 
          group.by = "Celltype1",  group.colors = hm_color) +
  scale_fill_gradientn(colors = c("cadetblue4", "white", "coral2"))
dev.off()

# ==============================================================================
# Figure 1D - Umap plot for Status and Flowcell
# ==============================================================================
DimPlot(coh.combined.sct, group.by = "Status", raster = FALSE) 
DimPlot(coh.combined.sct, group.by = "orig.ident", raster = FALSE)

# ==============================================================================
# Figure 1E - Umap plot for CT1 lineage
# ==============================================================================
DimPlot(coh.combined.sct, group.by = "Celltype1", label = T,
        cols = hm_color) + NoLegend()

# ==============================================================================
# Figure 1F - Cell proportion and mutation profile plot
# ==============================================================================
# Cell proportion data prep
coh.combined.sct$Celltype1 <- as.character(coh.combined.sct$Celltype1)
sub <- list()
j=0
for(i in unique(coh.combined.sct$orig.ident)){
  j=j+1
  sub[[j]] <- subset(coh.combined.sct,
                     cells = row.names(coh.combined.sct@meta.data[coh.combined.sct@meta.data$orig.ident == i, ]))
}
for(i in 1:length(sub)){
  names(sub) <- lapply(sub, function(xx){paste(unique(xx@meta.data$orig.ident))})
}

data_list <- lapply(sub, function(xx){as.data.frame(table(xx$Celltype1))})
data_list <- lapply(data_list, function(xx){cbind(xx,prop.table(xx$Freq)*100)})
data_table <- melt(data_list)
data_table <- data_table[!data_table$variable == "Freq",]
colnames(data_table) <- c("Celltype","Variable","Value","Patient_ID")
data_table.df <- data_table[,-2] # remove variable column
data_table.wide <- reshape(data_table.df, idvar = "Patient_ID", #convert to wide format
                           timevar = "Celltype", direction = "wide")
colnames(data_table.wide) <- c("Patient_ID","B cells", "Endothelial","Epithelial",
                               "FCGR3B+ neutrophils","Mast cells","Mesenchymal cells",
                               "Myeloid cells","T cells")

# Adjust patient ID to group into different pathology groups
level_order <- c("L35148","L35150","L35144","L35151","L35146","L35147","K28668",
                 "L49296","L35145","L42572","L49295","L35149","L35152","K28667",
                 "L49293","L49294")
data_table.wide <- data_table.wide[match(level_order, data_table.wide$Patient_ID),] 

# Count data
count_data <- as.data.frame(table(coh.combined.sct$orig.ident))
colnames(count_data) <- c("Patient_ID","Count")
count_data <- count_data[match(level_order, count_data$Patient_ID),]

# Pathology group color annotation
pathology_col <- c("LAMN" = "#FF8D68", "LGMA" = "#8AA0C9", "MHNA" = "#EE8BC2",
                   "GCA" = "#57C2A6", "Normal" = "#B0722C")
col = list(Pathology = pathology_col)

# Mutation data
mutation_data <- read.csv("/scratch/lbui/Appendiceal_data/appendix_mutation_data.csv")

# Adjust pathology to abbreviation
onion <- as.character(mutation_data$Pathology.Reviewed)
onion[onion == "Normal Appendix"] <- "Normal"
onion[onion == "Low-grade appendiceal mucinous neoplasm"] <- "LAMN"
onion[onion == "Low-grade mucinous adenocarcinoma"] <- "LGMA"
onion[onion == "Mod/High-grade Non-mucinous adenocarcinoma"] <- "MHNA"
onion[onion == "Goblet-cell adenocarcinoma"] <- "GCA"
mutation_data$Pathology_abb <- onion

# Combine cell proportion, mutation data and cell count into 1 matrix
mutation_data <- mutation_data[match(level_order, mutation_data$Library_ID),]

data_matrix <- cbind(mutation_data, data_table.wide[2:9],count_data[,2])
colnames(data_matrix) <- c("Patient_ID", "Sample", "Pathology_reviewed","Diagnosis",
                           "Status","Gender","Age","KRAS","GNAS","APC","TP53","Pathology",
                           "B cells","Endothelial","Epithelial","FCGR3B+ neutrophils",
                           "Mast cells","Mesenchymal cells","Myeloid cells","T cells",
                           "Total cell number")
data_matrix[is.na(data_matrix)] <- 0

# Generate the plot
color.palette  <- colorRampPalette(c("light grey","#3A3B3C"))(n=1000)
white.line <- function(j, i, x, y, w, h, fill) { 
  grid.lines(x = c(x - w/2, x + w / 2), y = c(y + h / 2, y + h / 2), 
             gp = gpar(col = 'white', lwd = 1)) }
fa <- colnames(data_matrix[,8:11])

lgd1 = Legend(
  title = gt_render("Mutation"), 
  labels = gt_render(c("Mutant", "Wild-type")),
  legend_gp = gpar(fill = 1:2)
)

ht1 <- Heatmap(as.matrix(data_matrix[,8:11]),
               col = color.palette,
               cluster_rows = F,
               cluster_columns = F,
               na_col = "white",
               column_title = "Mutations",
               column_names_side = "top",
               column_split = fa,
               cell_fun = white.line,
               show_heatmap_legend = F,
               width = ncol(data_matrix[,8:11])*unit(5, "mm")) 

# Cell proportion bar plot
ct1_color <- c("#d2331f","#0072ae","#43B425","#5f1a1a","#00785f","#00B6E8","#A58BFB","#FF62D5")
ha1 = rowAnnotation("Total cell proportion (%)"= anno_barplot(data_matrix[,13:20], 
                                                              gp = gpar(fill = ct1_color), border =F, bar_width = 1,
                                                              width = unit(12, "cm"), cell_fun = white.line))

# Total cell number plot (didn't use this in the figure)
#ha2 = rowAnnotation("Total cell number" = anno_points(data_matrix[,19],
#                                                      border = F,
#                                                      gp = gpar(fill = "grey")))

# Sample pathology label
ht2 <- Heatmap(data_matrix[,12], name = "Pathology", 
               col = c("#F3766E","#55BB82","#A4A433","#4AAFF0","#D775EC"),
               width = unit(5, "mm"),
               cell_fun = white.line,
               cluster_rows = F,
               show_row_names = F,
               row_km = 5,
               row_split = data_matrix$Pathology,
               row_title = NULL,
               gap = unit(2,"mm"))

# Add all plots together
ht_list = ht2 + #rowAnnotation(rn = anno_text(data_matrix$Patient_ID, 
                #                             location = unit(0, "npc"), 
                #                             just = "left")) + 
  ha1 + ht1 

lgd2 = Legend(
  title = gt_render("Cell type"), 
  labels = gt_render(c("B cells", "Endothelial", "Epithelial", "FCGR3B+ neutrophils",
                       "Mast cells","Mesenchymal cells", "Myeloid cells", "T cells")),
  legend_gp = gpar(fill = ct1_color)
)
lgd_list = list(lgd1, lgd2)
draw(ht_list, ht_gap = unit(2, "mm"), annotation_legend_list = lgd_list)

