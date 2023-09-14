# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/07
# Description: Mesenchymal cell granular annotation + figures
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

# ======================================
# Load libraries
# ======================================
#Sys.unsetenv("GITHUB_PAT")
library(Seurat)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Read in the Mesenchymal object
# ==============================================================================
mesen <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_Mesenchymal_CT2.rds")
# Subset out CAFs
caf <- subset(mesen, subset = Celltype2 == "CAFs")

# ==============================================================================
# DEG pathology in each CAF subtype
# ==============================================================================
app_list = list()
j=0
for(i in unique(caf@meta.data$Celltype3)){
  j=j+1
  app_list[[j]] <- subset(caf, 
                          subset = Celltype3 == i)
}
for(i in 1:length(app_list)){
  names(app_list) <- lapply(app_list, function(xx){paste(unique(xx@meta.data$Celltype3))})
  
}

# Run DEG using the RNA assay (PrepSCTfindmarkers gave TRUE/FALSE error, maybe too low cell numbers)
caf_pathology_deg <- lapply(app_list, function(xx){
  print(unique(xx@meta.data$Celltype3))
  Idents(xx) <- as.character(xx$Pathology2)   
  FindAllMarkers(xx, 
                 test.use = "negbinom",
                 logfc.threshold = 0,
                 assay = "RNA",
                 latent.vars = "Flowcell")
} )

saveRDS(caf_pathology_deg, file = "CAF_pathology_RNA_DEG.rds")
for(i in 1:length(caf_pathology_deg)){
  write.table(caf_pathology_deg[[i]], 
              paste(gsub("/", "", unique(app_list[[i]]@meta.data$Celltype3)), 
                    "_pathology_DEGs_findallmarkers", ".csv"), sep =",", quote = F)
}
# Select top 10 DEGs for heatmap
top_genes <- lapply(caf_pathology_deg, function(xx) 
  xx %>% group_by(cluster) %>% top_n(10, avg_log2FC))
top_genes.df <- lapply(top_genes, function(xx) xx$gene)
hm_genes <- unique(unlist(top_genes.df))
ribo <- grep("^RP", hm_genes, value = TRUE)
hm_genes <- hm_genes[!hm_genes %in% ribo]

# iCAFs
icaf_genes <- caf_pathology_deg[[1]] %>% group_by(cluster) %>% top_n(20, avg_log2FC)
icaf <- subset(caf, subset = Celltype3 == "iCAFs")
DoHeatmap(icaf, features = unique(icaf_genes$gene), group.by = "Pathology2", 
          assay = "SCT" , group.colors = pathology_col)

# fiCAFs
ficaf_genes <- caf_pathology_deg[[2]] %>% group_by(cluster) %>% top_n(25, avg_log2FC)
ficaf <- subset(caf, subset = Celltype3 == "fiCAFs")
DoHeatmap(ficaf, features = unique(ficaf_genes$gene), group.by = "Pathology2", 
          assay = "SCT", group.colors = pathology_col)

# myCAFs
mycaf_genes <- caf_pathology_deg[[3]] %>% group_by(cluster) %>% top_n(25, avg_log2FC)
mycaf <- subset(caf, subset = Celltype3 == "myCAFs")
DoHeatmap(mycaf, features = unique(mycaf_genes$gene), group.by = "Pathology2", 
          assay = "SCT", group.colors = pathology_col)

# apCAFs
apcaf_genes <- caf_pathology_deg[[4]] %>% group_by(cluster) %>% top_n(25, avg_log2FC)
ribo <- grep("^RP", apcaf_genes$gene, value = TRUE)
apcaf_genes <- apcaf_genes[!apcaf_genes$gene %in% ribo,]
apcaf <- subset(caf, subset = Celltype3 == "apCAFs")
DoHeatmap(apcaf, features = unique(apcaf_genes$gene),group.by = "Pathology2", 
          assay = "SCT", group.colors = pathology_col)

