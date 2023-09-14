# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/07
# Description: Mesenchymal cell figures
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
#Sys.unsetenv("GITHUB_PAT")
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
library(speckle)
library(CellChat)

# ==============================================================================
# Read in the Mesenchymal object
# ==============================================================================
mesen <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_Mesenchymal_CT2.rds")

# ==============================================================================
# Figure 5a
# ==============================================================================
ct2_color <- c("Myofibroblasts" = "#75cdc1", "CAFs" = "#dd7237", "Mesothelial" = "#0078b1",
               "Fibroblasts" = "#12a037", "Pericytes" = "#f12c89", "SMC" = "#ff91e4")
DimPlot(mesen, group.by = "Celltype2", cols = ct2_color)

# ==============================================================================
# Figure 5b
# ==============================================================================
propeller_test <- propeller(clusters = mesen$Celltype2, 
                            sample = mesen$orig.ident, 
                            group = mesen$Pathology)

write.csv(propeller_test, file = "Mesenchymal_CT2_propellertest.csv")

# Plot cell type proportions
plotCellTypeProps(clusters=mesen$Celltype2, sample=mesen$Pathology2) 

# ==============================================================================
# Figure 5c
# ==============================================================================
# Subset out CAFs
caf <- subset(mesen, subset = Celltype2 == "CAFs")
DefaultAssay(caf) <- "integrated"
caf <- RunPCA(caf)
caf <- RunUMAP(caf, dims = 1:6)

caf_col <- nord("aurora",4)
DimPlot(caf, group.by = "Celltype3", cols = caf_col)

# ==============================================================================
# Figure 5d
# ==============================================================================
# DEG between different CAF subtypes
Idents(caf) <- caf$Celltype3
caf_degs <- FindAllMarkers(caf,
                           test.use = "negbinom",
                           assay = "RNA",
                           latent.vars = "Flowcell")

write.csv(caf_degs, file = "CAF_subtype_allmarkers_RNA.csv")
top_genes <- caf_degs %>% group_by(cluster) %>% top_n(10, avg_log2FC)
ribo <- grep("^RP", top_genes$gene, value = TRUE)

DoHeatmap(subset(caf, downsample=700), 
          features = unique(top_genes[!top_genes$gene %in% ribo,]$gene), 
          group.by = "Celltype3",
          draw.lines = F, group.colors=caf_col) +
  scale_fill_gradientn(colors = c("cadetblue4", "white", "coral2"))

# ==============================================================================
# Figure 5e
# ==============================================================================
# Cell proportion significant test
propeller_test <- propeller(clusters = caf$Celltype3, 
                            sample = caf$orig.ident, 
                            group = caf$Pathology)

write.csv(propeller_test, file = "Mesenchymal_CAF_CT3_propellertest.csv")

# Plot cell type proportions
plotCellTypeProps(clusters=caf$Celltype3, sample=caf$Pathology2) 

# ==============================================================================
# Figure 5f
# ==============================================================================
# Perform GO enrichment for each CT in the lineage
library(clusterProfiler)
library(org.Hs.eg.db)
caf_deg_list <- split(caf_degs, f = caf_degs$cluster)
caf_genes <- lapply(caf_deg_list, function(x) (x[x$p_val_adj <= 0.1,]$gene))
caf_ribo <- lapply(caf_genes, function(x) grep("^RP", x, value = TRUE))
caf_mito <- lapply(caf_genes, function(x) grep("^MT", x, value = TRUE))
caf_ribo <- unlist(caf_ribo)
caf_mito <- unlist(caf_mito)
caf_genes <- lapply(caf_genes, function(x) x[!x %in% c(caf_ribo,caf_mito)]) #remove ribosomal genes
caf_genes <- lapply(caf_genes, function(x){
  bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
})

caf_genesentrez <- lapply(caf_genes, function(x) x$ENTREZID)
str(caf_genesentrez)

compared_go <- compareCluster(geneCluster = caf_genesentrez, 
                              fun = "enrichGO", OrgDb = org.Hs.eg.db, ont="BP",
                              pAdjustMethod = "BH")
dotplot(compared_go, showCategory=4, includeAll=FALSE, font.size=10,
        color="qvalue") +
  theme( 
    legend.text=element_text(size=5),
    legend.position = "top",
    legend.box = "vertical",
    text = element_text(size=6))

compared_pw <- compareCluster(geneCluster = caf_genesentrez, 
                              fun = "enrichPathway", 
                              pAdjustMethod = "BH")
dotplot(compared_pw, showCategory=5,includeAll=FALSE, font.size=10,
        color="qvalue") 

# ==============================================================================
# Figure 5g - Cellchat for SPP1+ TAMs, C1Qhi monocytes, CAFs
# ==============================================================================
# Read in the whole object
coh.combined.sct <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")

# Subset out the cells of interest
subset_ct <- c("SPP1+ TAMs", "C1Qhi monocytes", "iCAFs", "Fibroblasts","Monocyte-like",
               "Macrophages", "fiCAFs","myCAFs","apCAFs")
subset_ob <- subset(coh.combined.sct, subset = Celltype3 %in% subset_ct)
subset_ob <- subset(subset_ob, subset = Pathology2 != "Normal")

# Set the ligand-receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# Set database to use
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

# Create a Cellchat object
subset_ob <- PrepSCTFindMarkers(subset_ob)
cellchat_ob <- createCellChat(object = subset_ob, 
                              meta = subset_ob@meta.data, 
                              group.by = "Celltype3",
                              assay = "SCT")

# Add meta data into the Cellchat object
cellchat_ob <- addMeta(cellchat_ob, meta = cellchat_ob@meta)
cellchat_ob <- setIdent(cellchat_ob, ident.use = "Celltype3") # set "labels" as default cell identity
levels(cellchat_ob@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_ob@idents)) # number of cell

cellchat_ob@DB <- CellChatDB.use

# Subset data to only keep genes related to the analysis to save memory
cellchat_ob <- subsetData(cellchat_ob)
future::plan("multisession", workers = 4) # do parallel

cellchat_ob <- identifyOverExpressedGenes(cellchat_ob)
cellchat_ob <- identifyOverExpressedInteractions(cellchat_ob) 

# Compute the communication probability and infer cellular communication network
cellchat_ob <- computeCommunProb(cellchat_ob)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_ob <- filterCommunication(cellchat_ob, min.cells = 40)

# Infer the cell-cell communication at a signaling pathway level
cellchat_ob <- computeCommunProbPathway(cellchat_ob)

# Calculate the aggregated cell-cell communication network
cellchat_ob <- aggregateNet(cellchat_ob)

# Compute the network centrality scores
cellchat_ob <- netAnalysis_computeCentrality(cellchat_ob, slot.name = "netP")

cellchat_ob

pdf("Myeloid_CAFs_cellchat_SCT.pdf", width = 20, height =16)
cellchat_cols <- c("C1Qhi monocytes" = "#ef862f", "SPP1+ TAMs" = "#9488C1",
                   "Fibroblasts" = "#12a037", "apCAFs" = "#bf616a", 
                   "myCAFs" = "#b48ead", "iCAFs" = "#bbc28b", "fiCAFs"="#d99d79",
                   "Monocyte-like" = "#69cdea", "Macrophages" = "#00bfc3")
netVisual_chord_gene(cellchat_ob, sources.use = c(2,7), 
                     targets.use = c(1,3,4,5,8), color.use = cellchat_cols, #slot.name = "netP",
                     lab.cex = 0.5,legend.pos.x = 30)
netVisual_chord_gene(cellchat_ob, sources.use = c(6,9), 
                     targets.use = c(1,3,4,5,8), color.use = cellchat_cols, #slot.name = "netP",
                     lab.cex = 0.5,legend.pos.x = 30)
dev.off()

