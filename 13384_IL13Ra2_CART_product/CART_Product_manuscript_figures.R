#==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 10/1/2022
# Description: CAR T product manuscript figures
#==============================================================================

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)
library(stringr)
library(UpSetR)
library(RColorBrewer)
library(qvalue)
library(ComplexHeatmap)
library(nord)
library(circlize)
library(gridExtra)
library(ggpubr)
library(DescTools)
library(scProportionTest)
library(qdap)
library(ggrepel)
library(survival)
library(survminer)
library(presto)
library(googlesheets4)
library(ggbeeswarm)
library(ggforce)

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
# Import data
#==============================================================================

product <- readRDS("/labs/banovich/BCTCSF/Heini/product_projecTILs_modules_scGSVA.rds")

# Previously removed cluster 12, 13, 14
product <- subset(product, subset = cluster %in% c("12", "13", "14"), invert = T)
product$cluster <- factor(product$cluster, levels = c("0", as.character(seq(1, 11))))

product$response <- factor(product$response, levels = c("Complete Response (CR)", "Stable Disease (SD)", "Partial Response (PR)", "Progression Disease (PD)"))
product$cloneType <- factor(product$cloneType, levels = c("Single (0 < X <= 1)", "Small (1 < X <= 5)", "Medium (5 < X <= 20)", "Large (20 < X <= 100)"))
product$cluster <- factor(product$cluster, levels = c("0", as.character(seq(1, 11))))
product$c_cluster <- paste0("C", product$cluster)
product$c_cluster <- factor(product$c_cluster, levels = c(paste0("C", c("0", as.character(seq(1, 13))))))

# Adding a new binary response variable
binary_outcome <- product$response
binary_outcome <- ifelse(binary_outcome == "Progression Disease (PD)", "PD",
                         ifelse(binary_outcome %in% c("Stable Disease (SD)", "Partial Response (PR)", "Complete Response (CR)"), "CR/SD/PR", "NA"))
product$binary_outcome <- binary_outcome

# Canonical T cell markers
gs4_deauth()
canonical_markers  <- gs4_get(url)
sheet_names(canonical_markers)
canonical_markers <- read_sheet(canonical_markers, sheet = "T cells, gene sets")

# Samples over and below median overall survival
upn_surv <- distinct(product@meta.data[,which(colnames(product@meta.data) %in% c("survival_time_months_from_surgery", "UPN", "Diagnosis.Histology", "response"))])
upn_surv <- upn_surv[complete.cases(upn_surv),]
upn_surv$median_surv <- ifelse(upn_surv$UPN %in% upn_surv[which(upn_surv$survival_time_months_from_surgery>median(as.numeric(upn_surv$survival_time_months_from_surgery))),]$UPN, "above", "below")

ggplot(upn_surv, aes(response, after_stat(count))) + geom_bar(aes(fill = median_surv), position = "dodge")

# Adding below/above to the Seurat object
median_surv <- ifelse(product$UPN %in% upn_surv[which(upn_surv$survival_time_months_from_surgery>median(as.numeric(upn_surv$survival_time_months_from_surgery))),]$UPN, "above",
                      ifelse(product$UPN %in% upn_surv[which(upn_surv$survival_time_months_from_surgery<median(as.numeric(upn_surv$survival_time_months_from_surgery))),]$UPN, "below", "NA"))
product@meta.data$median_surv <- median_surv

# GBM only
# Samples over and below median overall survival
product_gbm <- subset(product, subset = Diagnosis.Histology == "Glioblastoma, NOS")
upn_surv_gbm <- distinct(product_gbm@meta.data[,which(colnames(product_gbm@meta.data) %in% c("survival_time_months_from_surgery", "UPN", "Diagnosis.Histology", "response"))])
upn_surv_gbm <- upn_surv_gbm[complete.cases(upn_surv_gbm),]

upn_surv_gbm$median_surv <- ifelse(upn_surv_gbm$UPN %in% upn_surv_gbm[which(upn_surv_gbm$survival_time_months_from_surgery>median(as.numeric(upn_surv_gbm$survival_time_months_from_surgery))),]$UPN, "above", "below")

ggplot(upn_surv_gbm, aes(response, after_stat(count))) + geom_bar(aes(fill = median_surv), position = "dodge")

# Adding below/above to the Seurat object
median_surv_gbm <- ifelse(product_gbm$UPN %in% product_gbm[which(product_gbm$survival_time_months_from_surgery>median(as.numeric(product_gbm$survival_time_months_from_surgery))),]$UPN, "above",
                      ifelse(product_gbm$UPN %in% product_gbm[which(product_gbm$survival_time_months_from_surgery<median(as.numeric(product_gbm$survival_time_months_from_surgery))),]$UPN, "below", "NA"))
product_gbm@meta.data$median_surv <- median_surv_gbm

#==============================================================================
# TCM vs. Tn/mem
#==============================================================================

# Scatter plot of cell type proportions
prop_table <- as.data.frame(table(product@meta.data[,"cluster"], as.character(product@meta.data[,"Manufacture"])))
colnames(prop_table) <- c("Cluster", "Manufacture", "Freq")
prop_table <- spread(prop_table, Cluster, Freq)
# Converting to percetange
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, Cluster, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)
prop_table <- spread(prop_table, Manufacture, Freq)

prop_table$Cluster <- paste0("C", prop_table$Cluster)

celltype_prop_scatter <- ggplot(prop_table, aes(x = TCM, y = Tnmem, color = Cluster, label = Cluster)) +
    geom_point() + 
    scale_color_manual(name = "Cluster", values = product_cluster_col_c) + 
    theme_bw() +
    manuscript_theme +
    ylab("Tn/mem") + 
    xlab("TCM") +
    xlim(0, 25) + 
    ylim(0, 25) +
    geom_text_repel(size = 2) +
    theme(legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
    theme(panel.background = element_rect(colour = "black")) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

celltype_prop_scatter

celltype_prop_scatter <- celltype_prop_scatter + coord_fixed(ratio = 1)

# Saving as a pdf
filename <- "/home/hnatri/CART/Product_manuscript/product_manufacture_celltypeprop_scatter.pdf"
pdf(file = filename,
    #res = 100,
    width = 3,
    height = 3)
celltype_prop_scatter
dev.off()

# Significance testing with propeller
library(speckle)

# Run propeller testing for cell type proportion differences between the two 
# groups
prop_res <- propeller(clusters = product$cluster, sample = product$UPN,
                      group = product$Manufacture)

#==============================================================================
# DimPlots
#==============================================================================

# Cluster
dimplot_product_cluster <- DimPlot(product, reduction = "wnn.umap", group.by = "c_cluster", cols = product_cluster_col_c, label = T, label.box = F, repel = T, raster = T, label.size = 2) +
    ggtitle("") +
    xlab("WNN UMAP 1") +
    ylab("WNN UMAP 2") +
    theme_classic() +
    NoLegend() +
    manuscript_theme +
    coord_fixed(1) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

# Manufacture
manufacture_col_pretty <- manufacture_col
names(manufacture_col_pretty) <- gsub("Tnmem", "Tn/mem", names(manufacture_col_pretty))
product$Manufacture_pretty <- gsub("Tnmem", "Tn/mem", product$Manufacture)
dimplot_product_manufacture <- DimPlot(product, reduction = "wnn.umap", group.by = "Manufacture_pretty", cols = manufacture_col_pretty, raster = T) +
    ggtitle("") +
    xlab("WNN UMAP 1") +
    ylab("WNN UMAP 2") +
    theme_classic() +
    #NoLegend() +
    my_theme + manuscript_theme +
    theme(legend.position = c(0.8,0.9)) +
    theme(legend.text = element_text(size=6)) + 
    theme(legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) +
    coord_fixed(1) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))


# CD8
featureplot_product_cd8 <- FeaturePlot(product, features = c("CD8A"), reduction="wnn.umap", raster=T) +
    ggtitle("") +
    xlab("WNN UMAP 1") +
    ylab("WNN UMAP 2") +
    theme_classic() +
    #NoLegend() +
    theme_classic() +
    my_theme + manuscript_theme +
    theme(legend.position = c(0.8,0.9)) +
    scale_color_gradientn(colors = c("gray89", "hotpink3", "deeppink3")) +
    theme(legend.key.height = unit(0.2, 'cm'), 
          legend.key.width = unit(0.4, 'cm'), 
          legend.text = element_text(size=6)) + +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) +
    coord_fixed(1) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))


DefaultAssay(product) <- "Protein"
featureplot_product_cd8_cite <- FeaturePlot(product, features = c("CD8"), reduction="wnn.umap", raster=T) +
    ggtitle("") +
    xlab("WNN UMAP 1") +
    ylab("WNN UMAP 2") +
    theme_classic() +
    #NoLegend() +
    theme_classic() +
    my_theme + manuscript_theme +
    theme(legend.position = c(0.8,0.9)) +
    scale_color_gradientn(colors = c("gray89", "hotpink3", "deeppink3")) +
    theme(legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(0.4, 'cm'),
          legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) +
    coord_fixed(1) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))


DefaultAssay(product) <- "RNA"

# Proportion of mitochondrial reads
featureplot_product_mt <- FeaturePlot(product, features = c("percent.mt"), reduction="wnn.umap", raster=T) +
    ggtitle("") +
    xlab("WNN UMAP 1") +
    ylab("WNN UMAP 2") +
    theme_classic() +
    my_theme + manuscript_theme +
    theme(legend.position = c(0.8,0.9)) +
    scale_color_gradientn(colors = c("darkgreen", "khaki1", "red4")) +
    theme(legend.key.height = unit(0.3, 'cm'), 
          legend.key.width = unit(0.4, 'cm'),
          legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) +
    coord_fixed(1) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

dimplot_product_phase <- DimPlot(product, reduction = "wnn.umap", group.by = "Phase", cols = cell_cycle_col, raster = T) +
    ggtitle("") +
    xlab("WNN UMAP 1") +
    ylab("WNN UMAP 2") +
    theme_classic() +
    my_theme + manuscript_theme +
    theme(legend.position = c(0.8,0.9)) +
    theme(legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) +
    coord_fixed(1) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

dimplot_product_phase

# Patching all plots together
(dimplot_product_cluster | featureplot_product_cd8) /
 (dimplot_product_manufacture | featureplot_product_mt) /
    (dimplot_product_phase | plot_spacer())

(dimplot_product_cluster | dimplot_product_manufacture) /
    (featureplot_product_cd8 | featureplot_product_cd8_cite)

filename <- "/home/hnatri/CART/Product_manuscript/product_4_patched_dimplots.pdf"
pdf(file = filename,
    width = 4, # inches
    height = 4) # inches
(dimplot_product_cluster | dimplot_product_manufacture) /
    (featureplot_product_cd8 | featureplot_product_cd8_cite)
dev.off()

FeaturePlot(product, features = c("Memory", "Dysfunction", "Cytotoxic"), reduction = "wnn.umap", ncol = 2, split.by = "Manufacture", min.cutoff = 0) &
    ggtitle("") &
    xlab("WNN UMAP 1") &
    ylab("WNN UMAP 2") &
    theme_classic() &
    #NoLegend() +
    my_theme & manuscript_theme &
    theme(legend.position = c(0.2,0.9)) &
    scale_color_gradientn(colors = c("white", "hotpink3", "deeppink3")) &
    theme(#legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        #legend.title = element_text(size=14),
        legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) &
    coord_fixed(1) &
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

FeaturePlot(product, features = c("Memory", "Dysfunction", "Cytotoxic"), reduction = "wnn.umap", ncol = 3, min.cutoff = 0) &
    #ggtitle("") &
    xlab("WNN UMAP 1") &
    ylab("WNN UMAP 2") &
    theme_classic() &
    #NoLegend() +
    my_theme & manuscript_theme &
    theme(legend.position = c(0.2,0.9)) &
    #scale_color_gradientn(colors = c("white", "hotpink3", "deeppink3")) &
    scale_color_viridis(option="magma") &
    theme(#legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        #legend.title = element_text(size=14),
        legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) &
    coord_fixed(1) &
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

memory <- FeaturePlot(product, features = "Memory", reduction = "wnn.umap", min.cutoff = 0, raster = T) &
    #ggtitle("") &
    xlab("WNN UMAP 1") &
    ylab("WNN UMAP 2") &
    theme_classic() &
    #NoLegend() +
    my_theme & manuscript_theme &
    theme(legend.position = c(0.8,0.9)) &
    #scale_color_gradientn(colors = c("white", "hotpink3", "deeppink3")) &
    scale_color_viridis(option="magma", limits = c(0, 0.8)) &
    theme(#legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        #legend.title = element_text(size=14),
        legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) &
    coord_fixed(1) &
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

cytotoxic <- FeaturePlot(product, features = "Cytotoxic", reduction = "wnn.umap", min.cutoff = 0, raster = T) &
    #ggtitle("") &
    xlab("WNN UMAP 1") &
    ylab("WNN UMAP 2") &
    theme_classic() &
    #NoLegend() +
    my_theme & manuscript_theme &
    theme(legend.position = c(0.8,0.9)) &
    #scale_color_gradientn(colors = c("white", "hotpink3", "deeppink3")) &
    scale_color_viridis(option="magma", limits = c(0, 0.8)) &
    theme(#legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        #legend.title = element_text(size=14),
        legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) &
    coord_fixed(1) &
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

dysfunction <- FeaturePlot(product, features = "Dysfunction", reduction = "wnn.umap", min.cutoff = 0, raster = T) &
    #ggtitle("") &
    xlab("WNN UMAP 1") &
    ylab("WNN UMAP 2") &
    theme_classic() &
    #NoLegend() +
    my_theme & manuscript_theme &
    theme(legend.position = c(0.8,0.9)) &
    #scale_color_gradientn(colors = c("white", "hotpink3", "deeppink3")) &
    scale_color_viridis(option="magma", limits = c(0, 0.8)) &
    theme(#legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        #legend.title = element_text(size=14),
        legend.text = element_text(size=6)) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(-1,5,1,1))) &
    coord_fixed(1) &
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

filename <- "/home/hnatri/CART/Product_manuscript/signature_featureplot.pdf"
pdf(file = filename,
    width = 6, # inches
    height = 2) # inches
(memory | cytotoxic | dysfunction)
dev.off()

#==============================================================================
# Heatmap dotplot
#==============================================================================

# log-normalizing and scaling all features in the RNA assay. Scaling so that
# all features can be visualized using the same color scale
DefaultAssay(product) <- "RNA"
product <- NormalizeData(product)
VariableFeatures(product) <- rownames(product)
product <- ScaleData(product)

DefaultAssay(product) <- "Protein"
product <- NormalizeData(product)
VariableFeatures(product) <- rownames(product)
product <- ScaleData(product)

# Finding cluster markers for a heatmap
product_plot <- product
DefaultAssay(product_plot) <- "RNA"

# Removing RP and MT genes
keep_features <- rownames(product_plot)
keep_features <- keep_features[!grepl("RP", keep_features)]
keep_features <- keep_features[!grepl("MT-", keep_features)]

product_plot <- subset(product_plot, features = keep_features)

markers <- presto::wilcoxauc(product_plot, group_by = "cluster", assay = "data", seurat_assay = "Protein")
#markers <- top_markers(markers, n = 3, auc_min = 0.5, pct_in_min = 20, pct_out_max = 100)
markers <- presto::top_markers(markers, n = 8, auc_min = 0.5, pct_in_min = 30, pct_out_max = 100)
#markers <- top_markers(markers, n = 5, auc_min = 0.5)

markers

all_markers <- markers %>%
    dplyr::select(-rank) %>% 
    unclass() %>% 
    stack() %>%
    pull(values) %>%
    unique() %>%
    .[!is.na(.)]

length(all_markers)

product_plot$c_cluster <- paste0("C", product_plot$cluster)

p <- DotPlot(object = product_plot, features = all_markers, group.by = "c_cluster")
p

dotplot_heatmap <- create_dotplot_heatmap(seurat_object = product_plot,
                                          plot_features = all_markers,
                                          group_var = "c_cluster",
                                          group_colors = product_cluster_col_c,
                                          column_title = "Cluster")

dotplot_heatmap <- create_dotplot_heatmap_horizontal(seurat_object = product_plot,
                                                     plot_features = all_markers,
                                                     group_var = "c_cluster",
                                                     group_colors = product_cluster_col_c,
                                                     column_title = "")

dotplot_heatmap

filename <- "/home/hnatri/CART/Product_manuscript/product_heatmap.pdf"
pdf(file = filename,
    width = 6, # inches
    height = 6) # inches
dotplot_heatmap
dev.off()

# Canonical markers
DefaultAssay(product) <- "Protein"
marker_dotplot_heatmap <- create_dotplot_heatmap_horizontal(seurat_object = product,
                                                            plot_features = unique(canonical_markers$Protein),
                                                            group_var = "c_cluster",
                                                            group_colors = product_cluster_col_c,
                                                            column_title = "")

#==============================================================================
# Violin plots
#==============================================================================

DefaultAssay(product) <- "RNA"
VlnPlot(product, group.by = "Manufacture", features = c("PTPRC", "CD4", "CD8A", "SELL", "CD27", "CCR7"), ncol = 3, pt.size = 0, cols = col$manufacture) &
    xlab("") &
    ylab("Gene expression") &
    my_theme

DefaultAssay(product) <- "Protein"
VlnPlot(product, group.by = "Manufacture", features = c("CD45RA", "CD4", "CD8", "CD62L", "CD27", "CD197-CCR7"), ncol = 3, pt.size = 0, cols = col$manufacture) &
    xlab("") &
    ylab("Protein expression") &
    my_theme

#==============================================================================
# Outcome comparison
#==============================================================================

# Scatter plot of cell type proportions
prop_table <- as.data.frame(table(product@meta.data[,"cluster"], as.character(product@meta.data[,"binary_outcome"])))
colnames(prop_table) <- c("Cluster", "Outcome", "Freq")
prop_table <- spread(prop_table, Cluster, Freq)
# Converting to percetange
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, Cluster, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)
prop_table <- spread(prop_table, Outcome, Freq)

prop_table$Cluster <- paste0("C", prop_table$Cluster)
colnames(prop_table) <- c("Cluster", "CR_SD_PR", "PD")

celltype_prop_scatter_1 <- ggplot(prop_table, aes_string(x = "PD", y = "CR_SD_PR", color = "Cluster", label = "Cluster")) +
    geom_point() + 
    scale_color_manual(name = "Cluster", values = product_cluster_col_c) + 
    theme_bw() +
    manuscript_theme +
    ylab("CR\\/SD\\/PR") + 
    xlab("PD") +
    xlim(0, 25) + 
    ylim(0, 25) +
    geom_text_repel(size = 2) +
    theme(legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
    theme(panel.background = element_rect(colour = "black")) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

celltype_prop_scatter_1

celltype_prop_scatter_1 <- celltype_prop_scatter_1 + coord_fixed(ratio = 1)

# Above or below median survival
product_tcm <- subset(product_gbm, subset = Manufacture == "TCM")
product_tnmem <- subset(product_gbm, subset = Manufacture == "Tnmem")

prop_table <- as.data.frame(table(product@meta.data[,"cluster"], as.character(product@meta.data[,"median_surv"])))
colnames(prop_table) <- c("Cluster", "Median_surv", "Freq")
prop_table <- spread(prop_table, Cluster, Freq)
# Converting to percetange
prop_table[,2:length(prop_table)] <- (prop_table[,2:length(prop_table)]/rowSums(prop_table[,2:length(prop_table)]))*100
prop_table <- gather(prop_table, Cluster, Freq, names(prop_table)[2:length(names(prop_table))], factor_key=TRUE)
prop_table <- spread(prop_table, Median_surv, Freq)

prop_table$Cluster <- paste0("C", prop_table$Cluster)
colnames(prop_table) <- c("Cluster", "Above", "Below")

celltype_prop_scatter_11 <- ggplot(prop_table, aes_string(x = "Below", y = "Above", color = "Cluster", label = "Cluster")) +
    geom_point() + 
    scale_color_manual(name = "Cluster", values = product_cluster_col_c) + 
    theme_bw() +
    manuscript_theme +
    ylab("Above") + 
    xlab("Below") +
    xlim(0, 35) + 
    ylim(0, 35) +
    geom_text_repel(size = 2) +
    theme(legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
    theme(panel.background = element_rect(colour = "black")) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

celltype_prop_scatter_11

celltype_prop_scatter_11 <- celltype_prop_scatter_11 + coord_fixed(ratio = 1)

# Saving as a pdf
filename <- "/home/hnatri/CART/Product_manuscript/product_mediansurv_celltypeprop_scatter_tnmem_gbm.pdf"
pdf(file = filename,
    #res = 100,
    width = 2,
    height = 2)
celltype_prop_scatter_11
dev.off()

# Tnmem
product_tnmem <- subset(product, subset = Manufacture == "Tnmem")
prop_table2 <- as.data.frame(table(product_tnmem@meta.data[,"cluster"], as.character(product_tnmem@meta.data[,"binary_outcome"])))
colnames(prop_table2) <- c("Cluster", "Outcome", "Freq")
prop_table2 <- spread(prop_table2, Cluster, Freq)
# Converting to percetange
prop_table2[,2:length(prop_table2)] <- (prop_table2[,2:length(prop_table2)]/rowSums(prop_table2[,2:length(prop_table2)]))*100
prop_table2 <- gather(prop_table2, Cluster, Freq, names(prop_table2)[2:length(names(prop_table2))], factor_key=TRUE)
prop_table2 <- spread(prop_table2, Outcome, Freq)

prop_table2$Cluster <- paste0("C", prop_table2$Cluster)
colnames(prop_table2) <- c("Cluster", "CR_SD_PR", "PD")

celltype_prop_scatter_2 <- ggplot(prop_table2, aes_string(x = "PD", y = "CR_SD_PR", color = "Cluster", label = "Cluster")) +
    geom_point() + 
    scale_color_manual(name = "Cluster", values = product_cluster_col_c) + 
    theme_bw() +
    manuscript_theme +
    ylab("CR\\/SD\\/PR") + 
    xlab("PD") +
    xlim(0, 35) + 
    ylim(0, 35) +
    geom_text_repel(size = 2) +
    theme(legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
    theme(panel.background = element_rect(colour = "black")) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

celltype_prop_scatter_2

celltype_prop_scatter_2 <- celltype_prop_scatter_2 + coord_fixed(ratio = 1)

# TCM
product_tcm <- subset(product, subset = Manufacture == "TCM")
prop_table3 <- as.data.frame(table(product_tcm@meta.data[,"cluster"], as.character(product_tcm@meta.data[,"binary_outcome"])))
colnames(prop_table3) <- c("Cluster", "Outcome", "Freq")
prop_table3 <- spread(prop_table3, Cluster, Freq)
# Converting to percetange
prop_table3[,2:length(prop_table3)] <- (prop_table3[,2:length(prop_table3)]/rowSums(prop_table3[,2:length(prop_table3)]))*100
prop_table3 <- gather(prop_table3, Cluster, Freq, names(prop_table3)[2:length(names(prop_table3))], factor_key=TRUE)
prop_table3 <- spread(prop_table3, Outcome, Freq)

prop_table3$Cluster <- paste0("C", prop_table3$Cluster)
colnames(prop_table3) <- c("Cluster", "CR_SD_PR", "PD")

celltype_prop_scatter_3 <- ggplot(prop_table3, aes_string(x = "PD", y = "CR_SD_PR", color = "Cluster", label = "Cluster")) +
    geom_point() + 
    scale_color_manual(name = "Cluster", values = product_cluster_col_c) + 
    theme_bw() +
    manuscript_theme +
    ylab("CR\\/SD\\/PR") + 
    xlab("PD") +
    xlim(0, 25) + 
    ylim(0, 25) +
    geom_text_repel(size = 2) +
    theme(legend.position = "none") + 
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
    theme(panel.background = element_rect(colour = "black")) +
    theme(plot.margin = margin(0.3, 0.3, 0.3,0.3, "cm"))

celltype_prop_scatter_3

celltype_prop_scatter_3 <- celltype_prop_scatter_3 + coord_fixed(ratio = 1)

(celltype_prop_scatter_1 | celltype_prop_scatter_2 | celltype_prop_scatter_3)

# Saving as a pdf
filename <- "/home/hnatri/CART/Product_manuscript/product_outcome_celltypeprop_scatter.pdf"
pdf(file = filename,
    #res = 100,
    width = 6,
    height = 2)
(celltype_prop_scatter_1 | celltype_prop_scatter_2 | celltype_prop_scatter_3)
dev.off()


#==============================================================================
# Survival kaplan-meier
#==============================================================================

metadata <- product@meta.data
#metadata_survival <- merge(metadata, survival)
metadata <- metadata[,which(colnames(metadata) %in% c("UPN", "Diagnosis.Histology", "Manufacture", "survival_time_months_from_surgery", "death_status"))]
metadata <- distinct(metadata)

metadata$survival_time_months_from_surgery <- as.numeric(metadata$survival_time_months_from_surgery)
metadata$death_status <- as.numeric(metadata$death_status)

metadata_gbm <- metadata[which(metadata$Diagnosis.Histology=="Glioblastoma, NOS"),]
metadata_gbm <- metadata_gbm[,which(colnames(metadata_gbm) %in% c("UPN", "Manufacture", "survival_time_months_from_surgery", "death_status"))]
metadata_gbm <- distinct(metadata_gbm)

metadata$Manufacture <- gsub("Tnmem", "Tn/mem", metadata$Manufacture)

# Segregate survival metrics, Kaplan Meier on Tnmem and TCM
fit <- survfit(Surv(survival_time_months_from_surgery, death_status) ~ Manufacture,
               data = metadata, type="kaplan-meier")
surv_pvalue(fit)

names(fit$strata) <- gsub("Manufacture=", "", names(fit$strata))

surv_plot_manufacture <- ggsurvplot(fit, conf.int = 0.95,
                                    censor= F,
                                    #pval = T,
                                    cols = "strata",
                                    palette = as.character(manufacture_col),
                                    ggtheme = theme_bw(base_size=6),
                                    legend.title="") +
    labs(x = "Months from surgery", y = "Survival")

surv_plot_manufacture$plot

surv_plot <- surv_plot_manufacture$plot

surv_plot <- surv_plot +
    manuscript_theme +
    theme_bw() +
    theme(panel.background = element_rect(colour = "black")) +
    theme(legend.position = c(0.8, 0.9)) +
    theme(legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    theme(legend.title=element_blank(),
          legend.margin=margin(c(1,4,4,4)))

surv_plot
    
filename <- "/home/hnatri/CART/Product_manuscript/surv_plot_manufacture.pdf"
pdf(file = filename,
    width = 3, # inches
    height = 3) # inches
surv_plot
dev.off()
