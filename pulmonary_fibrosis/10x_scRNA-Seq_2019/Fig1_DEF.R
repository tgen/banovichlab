# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
#             Austin J. Gutierrez, agutierrez@tgen.org
# Date: 11/06/2019
# Description: Figure 1: D, E, F
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
set.seed(2811)

# ======================================
# Load libraries
# ======================================
library(Seurat)
library(dplyr)
library(ggplot2)

# ===========================================
# Read in Seurat object
# ===========================================
ild <- readRDS("ILD.rds")

# ===========================================
# Running DE analysis using negative binominal test
# ===========================================
ild_list = list()
j <- 0
for(i in unique(ild@meta.data$celltype)){
  j <- j + 1
  ild_list[[j]] <- SubsetData(ild,
                              cells = row.names(ild@meta.data[ild@meta.data$celltype == i, ]))
}
for(i in 1:length(ild_list)){
  names(ild_list) <- lapply(ild_list, function(xx){paste(unique(xx@meta.data$celltype))})
}

disease_vs_control <- lapply(ild_list, function(xx){
  print(unique(xx@meta.data$celltype))
  if(length(unique(xx@meta.data$Status)) > 1) {
    FindMarkers(xx,
                group.by = "Status",
                ident.1 = "ILD",
                ident.2 = "Control",
                test.use = "negbinom")
  } 
  else{
    return(NULL)
  } 
})

DE_table <- NULL
for (i in 1:length(disease_vs_control)){
  DE_table <- rbind(DE_table, length(disease_vs_control[[i]]$p_val_adj <= .1))
}

row.names(DE_table) <- names(disease_vs_control)
DE_table

# ======================================
# Remove cell types with less than 50 cells/condition
# ======================================
DE_table <- DE_table[-c(1,9,18,26,28,30), ]
DE_table <- as.data.frame(DE_table)

# ======================================
# Add a column for population
# ======================================
onion <- as.character(row.names(DE_table))

onion[onion %in% c("B Cells",
                   "cDCs",
                   "NK Cells",
                   "Mast Cells",
                   "Plasma Cells", 
                   "Monocytes",
                   "Macrophages",
                   "Proliferating Macrophages",
                   "T Cells", 
                   "Proliferating T Cells")] <- "Immune"

onion[onion %in% c("Basal",
                   "Proliferating Epithelial Cells",
                   "Differentiating Ciliated", 
                   "Ciliated", "SCGB3A2+ SCGB1A1+",
                   "SCGB3A2+",
                   "MUC5B+", 
                   "Transitional AT2",
                   "AT2",
                   "AT1")] <- "Epithelial"

onion[onion %in% c("Lymphatic Endothelial Cells",
                   "Endothelial Cells")] <- "Endothelial"

onion[onion %in% c("Fibroblasts",
                   "Myofibroblasts",
                   "Smooth Muscle Cells")] <- "Mesenchymal"

DE_table$population <- onion
colnames(DE_table) <- c("counts", "population")

# ======================================
# Figure 1: D
# ======================================
# Make a bar graph for DE genes for each cell type in Disease vs. Control
ggplot(DE_table, 
       aes( x = reorder(row.names(DE_table),-counts),
            y = counts,
            color = population,
            fill = population)) + 
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(~population,
             scales = "free",
             space = "free") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x="Cell Type",
       y = "Number of differentially expressed genes") +
  NoLegend() 

# ======================================
epi_mes <- subset(ild,cells=rownames(ild@meta.data[ild@meta.data$population %in%
                                                      c("Epithelial", "Mesenchymal"),]))
epi_fibro_CT <- c("HAS1 High Fibroblasts",
                  "PLIN2+ Fibroblasts","Myofibroblasts",
                  "Fibroblasts",
                  "KRT5-/KRT17+",
                  "Transitional AT2",
                  "AT2",
                  "AT1",
                  "SCGB3A2+ SCGB1A1+",
                  "SCGB3A2+",
                  "Proliferating Epithelial Cells",
                  "MUC5B+",
                  "MUC5AC+ High",
                  "Differentiating Ciliated",
                  "Ciliated",
                  "Basal")

epi_fibro <- subset(epi_mes, cells = rownames(epi_mes@meta.data[epi_mes@meta.data$celltype 
                                                                %in% epi_fibro_CT,]))
ecm_genes <- c("ACTA2",
               "COL12A1",
               "COL14A1",
               "COL18A1",
               "COL1A1",
               "COL3A1",
               "COL4A1",
               "COL4A1",
               "COL4A2",
               "COL6A1",
               "COL6A2",
               "COL8A1",
               "DMBT1",
               "DPT",
               "EFEMP1",
               "EFEMP2",
               "EPCAM",
               "FBLN1",
               "FBLN2",
               "FBLN5",
               "FBN1",
               "FN1",
               "HSPG2",
               "KRT5",
               "KRT17",
               "LAMA2",
               "LAMA3",
               "LAMA5",
               "LAMB2",
               "LAMB3",
               "LAMC1",
               "LAMC2",
               "LTBP1",
               "MFAP2",
               "MFAP5",
               "MGP",
               "NPNT",
               "PDGFRA",
               "POSTN",
               "SERPINA1",
               "SERPINA3",
               "TGFBI",
               "TINAGL1",
               "TNC",
               "VCAN")

my_levels <- c("Ciliated",
               "Differentiating Ciliated",
               "SCGB3A2+ SCGB1A1+",
               "SCGB3A2+",
               "MUC5B+",
               "MUC5AC+ High",
               "Proliferating Epithelial Cells",
               "Basal",
               "Transitional AT2",
               "AT2",
               "AT1",
               "KRT5-/KRT17+",
               "Fibroblasts",
               "Myofibroblasts",
               "PLIN2+ Fibroblasts",
               "HAS1 High Fibroblasts")

epi_fibro@meta.data$celltype <- factor(epi_fibro@meta.data$celltype, levels = my_levels)
# ======================================
# Figure 1: E
# ======================================
DoHeatmap(subset(epi_fibro, downsample = 100),
          features=ecm_genes,
          group.by = "celltype", 
          angle = 45) + NoLegend()

# ======================================
epi_fibro@meta.data$Status <- factor(epi_fibro@meta.data$Status, 
                                     levels=c("ILD", "Control"))

# ======================================
# Figure 1: F
# ======================================
VlnPlot(epi_fibro, c("COL1A1","CDKN2A","MMP7","MUC5B","FN1","SMAD3"), 
        group.by = "celltype",
        split.by = "Status",
        pt.size = 0,
        adjust = 2,
        combine = FALSE) 




