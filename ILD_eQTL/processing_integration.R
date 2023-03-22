# ==============================================================================
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: Dec 2021
# Description: Sample proccessing for the ILD project in Seurat V4
# ==============================================================================

# ======================================
# Load libraries
# ======================================

library(data.table)
library(googlesheets4)
library(RCurl)
library(qdap)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ade4)
library(Matrix)
library(DoubletFinder)
library(tidyverse)
library(harmony)
library(findPC)

# ======================================
# Environment parameters
# ======================================

getwd()
#setwd("/scratch/hnatri/ILD/Seurat_objects")

options(future.globals.maxSize = 4096*1024^2 )
set.seed(2811)

# ======================================
# Helper functions
# ======================================

source("utilities.R")
source("lung_celltype_markers.R")
source("lung_celltype_colors.R")

# ======================================
# Metadata
# ======================================

# Preventing trying to getcredentials
gs4_deauth()

metadata <- gs4_get("https://docs.google.com/spreadsheets/d/1rCrE8KAelXHX_MHxi48119ZLmJhs29QZb1eGcc11mg0/edit?usp=sharing")
sheet_names(metadata)

ipf_metadata <- read_sheet(metadata, sheet = "IPF") 
hashed_samples_metadata <- read_sheet(metadata, sheet = "More_Less_hashing") 

# Fix Library_IDs
ipf_metadata$Library_ID <- as.character(ipf_metadata$Library_ID)
ipf_metadata$Library_ID <- multigsub(c("-less", "-more", "-negative"),
                                     c("-1", "-2", "-0"),
                                     ipf_metadata$Library_ID)
# Fix sites
ipf_metadata$Processing_site <- as.character(ipf_metadata$Processing_site)
ipf_metadata$Processing_site[ipf_metadata$Processing_site %in% c("LB", "LP","MC","JCR", "LP/RK", "MC/RK", "RK")] <- "TGen"
ipf_metadata$Processing_site[ipf_metadata$Processing_site %in% c("NW","CC","CC/CH","CH","VU","CT")] <- "Vanderbilt"

# COPD samples to be removed
unique(ipf_metadata$Clinical_Diagnosis)
copd_samples <- ipf_metadata[which(ipf_metadata$Clinical_Diagnosis=="COPD"),]$Library_ID
ipf_metadata <- ipf_metadata[-which(ipf_metadata$Clinical_Diagnosis=="COPD"),]

# ======================================
# Define the location of the Cellranger output files to be read into R
# Define the location of barcodes.tsv, genes.tsv, matrix.mtx
# Create a character vector of sample file names
# Read data from each sample folder and combine
# ======================================

dataset_loc <- "/labs/banovich/SingleCell/CellRanger/3_1_0/Ensemble_93/PipelineData/Projects/IPF/CellRangerOuts/GeneExpression/"
inputfiles_loc <- "outs/filtered_feature_bc_matrix"
batchids <- list.files(path = dataset_loc, pattern = "IPF*",
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE)

batchids <- batchids[-grep(".csv", batchids)]

d10x_data <- sapply(batchids, function(i){
  message(i)
  sample_id <- sapply(strsplit(i,split="_") , "[[", 8)
  #message(sample_id)
  d10x <- Read10X(file.path(dataset_loc,i,inputfiles_loc))
  colnames(d10x) <- paste(sample_id, sep="_", colnames(d10x))

  d10x_seu <- CreateSeuratObject(d10x, min.features = 0, min.cells = 0, names.field = 1, names.delim = "_")
  d10x_seu <- PercentageFeatureSet(object = d10x_seu, pattern = "^MT-", col.name = "percent.mt")
  
  # Add sample metadata
  d10x_seu$Library_ID <- d10x_seu$orig.ident
  sample_metadata <- ipf_metadata[which(ipf_metadata$Library_ID==i),]
  d10x_seu@meta.data$Sample_Name = plyr::mapvalues(x = d10x_seu@meta.data$Library_ID,
                                                     from = ipf_metadata$Library_ID,
                                                     to = ipf_metadata$Sample_Name)
  d10x_seu@meta.data$Sample_Type = plyr::mapvalues(x = d10x_seu@meta.data$Library_ID,
                                                     from = ipf_metadata$Library_ID,
                                                     to = ipf_metadata$Sample_Type)
  d10x_seu@meta.data$Flowcell_ID = plyr::mapvalues(x = d10x_seu@meta.data$Library_ID,
                                                    from = ipf_metadata$Library_ID,
                                                    to = ipf_metadata$Flowcell_ID)
  
  # Changing cell names
  d10x_seu <- RenameCells(d10x_seu, paste0(d10x_seu@meta.data$Sample_Name, "_", d10x_seu@meta.data$Sample_Type, "_", colnames(d10x_seu)))
  
  d10x_seu
})

# ======================================
# Processing hashed less fibrotic and more fibrotic samples
# ======================================

hashed_dataset_loc <- "/labs/banovich/SingleCell/CellRanger/3_1_0/Ensemble_93/PipelineData/Projects/IPF/CellRangerOuts/CITE"
inputfiles_loc <- "outs/filtered_feature_bc_matrix"
hashed_batchids <- list.files(path = hashed_dataset_loc, pattern = "*AB*",
                              full.names = FALSE, recursive = FALSE,
                              ignore.case = FALSE, include.dirs = TRUE)

hashed_10x_data <- sapply(hashed_batchids,  function(i){
    message(i)
    # i <- "F01724-GEX_F01725-AB"
    # Reading 10X data
    hashed_ild <- Read10X(file.path(hashed_dataset_loc,i,inputfiles_loc),
                          gene.column = 2,
                          cell.column = 1,
                          unique.features = TRUE,
                          strip.suffix = FALSE)
    
    GEX <- hashed_ild[[1]]
    AB <- hashed_ild[[2]]
    
    # Identifying hashing antibodies
    antibodies <- names(rowSums(AB)[order(-rowSums(AB))][1:2])
    
    # Remove unused antibodies
    AB <- AB[rownames(AB) %in% antibodies,]
    rowSums(AB)
    
    # Removing cells that have no HTO counts
    drop_cells <- colnames(AB[,which(colSums(AB)==0)])
    AB <- AB[,!colnames(AB) %in% drop_cells]
    # Find shared barcodes
    shared_barcodes <- intersect(colnames(GEX), colnames(AB))
    # Subset to only chared barcodes
    GEX <- GEX[, shared_barcodes]
    AB <- as.matrix(AB[, shared_barcodes])
    
    # Creating the object with the expression matrix
    hashed_ild <- CreateSeuratObject(counts = GEX)
    hashed_ild <- PercentageFeatureSet(hashed_ild, pattern = "^MT-", col.name = "percent.mt")
    
    # Creating the hash assay
    hashed_ild[["Hash"]] <- CreateAssayObject(counts = AB)
    
    # Normalizing
    hashed_ild <- NormalizeData(hashed_ild, assay = "Hash", normalization.method = "CLR")
    hashed_ild <- ScaleData(hashed_ild, assay = "Hash")
    
    # Assign single cells back to their sample origins
    # Demultiplexing based on the hashing antibodies
    hashed_ild <- HTODemux(hashed_ild, assay = "Hash", positive.quantile = 0.99, verbose = F)
    head(hashed_ild@meta.data)
    
    # Global hash classification
    #table(hash_test$Hash_classification.global)
    
    # Subset to remove doublets (need to use backtics because of the special characters)
    hashed_ild <- subset(hashed_ild, subset = `Hash_classification.global` == "Singlet")
    
    # Add sample metadata
    hashed_sample_metadata <- hashed_samples_metadata[which(hashed_samples_metadata$Library_ID==i),]
    hashed_ild@meta.data$Sample_Name = plyr::mapvalues(x = hashed_ild@meta.data$hash.ID,
                                                       from = hashed_sample_metadata$Hash_AB_name,
                                                       to = hashed_sample_metadata$Sample_Name)
    hashed_ild@meta.data$Sample_Type = plyr::mapvalues(x = hashed_ild@meta.data$hash.ID,
                                                      from = hashed_sample_metadata$Hash_AB_name,
                                                      to = hashed_sample_metadata$Sample_Type)
    hashed_ild@meta.data$Hash_Library_IDD = plyr::mapvalues(x = hashed_ild@meta.data$hash.ID,
                                                      from = hashed_sample_metadata$Hash_AB_name,
                                                      to = hashed_sample_metadata$Library_ID)
    
    hashed_ild@meta.data$Library_ID <- sapply(strsplit(hashed_ild@meta.data$Hash_Library_ID, split="-") , "[[", 1)
    
    # Changing cell names
    hashed_ild <- RenameCells(hashed_ild, paste0(hashed_ild@meta.data$Sample_Name, "_", hashed_ild@meta.data$Sample_Type, "_", colnames(hashed_ild)))

    hashed_ild
})

# Combining with the GEX libraries
hashed_10x_data <- lapply(hashed_10x_data, function(xx){
    #xx <- hashed_10x_data[[1]]
    xx@meta.data$Hash_Library_ID <- as.character(xx@meta.data$Library_ID)
    #sapply(strsplit(ild_test@meta.data$Library_ID, split="-") , "[[", 1)
    xx@meta.data$Library_ID <- sapply(strsplit(xx@meta.data$Hash_Library_ID, split="-") , "[[", 1)
    
    xx
})

gex_hashed_10x_data <- c(d10x_data, hashed_10x_data)
names(gex_hashed_10x_data)

saveRDS(gex_hashed_10x_data, "/scratch/hnatri/ILD/Seurat_objects/gex_hashed_10x_data.rds")
saveRDS(hashed_10x_data, "/scratch/hnatri/ILD/Seurat_objects/hashed_10x_data.rds")
saveRDS(d10x_data, "/scratch/hnatri/ILD/Seurat_objects/gex_10x_data.rds")
#gex_hashed_10x_data <- readRDS("gex_hashed_10x_data.rds")

# ======================================
# QC
# ======================================

# Merge all Seurat objects
ild <- Reduce(merge, gex_hashed_10x_data)

# Removing COPD samples
ild <- subset(ild, subset = orig.ident %in% setdiff(ild@meta.data$orig.ident, copd_samples))
ild <- PercentageFeatureSet(object = ild, pattern = "^MT-", col.name = "percent.mt")

plot_ild <- ild
Idents(plot_ild) <- "ident"

VlnPlot(plot_ild, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Make density plots
smoothScatter(ild@meta.data$percent.mt, ild@meta.data$nFeature_RNA)
abline(h = 1000, v = 25)
text(1.5,700, "nFeature = 1000, percent.mt = 25", adj = c(0, -.1))

filename <- "/home/hnatri/ILD_processing_annotation/smoothscatter.pdf"
pdf(file = filename,
    width = 4, # The width of the plot in inches
    height = 3) # The height of the plot in inches
smoothScatter(ild@meta.data$percent.mt, ild@meta.data$nFeature_RNA)
abline(h = 1000, v = 25)
text(1.5,700, "nFeature = 1000, percent.mt = 25", adj = c(0, -.1))
dev.off()

# Filter out low read count cells
ild2 <- subset(ild, subset = nFeature_RNA > 1000 & percent.mt < 25)
length(unique(ild2@meta.data$orig.ident))

# Count number of cells before and after filter
table(ild@meta.data$orig.ident)
write.csv(table(ild@meta.data$orig.ident), file = "cellcount_prefiltered.csv")
write.csv(table(ild2@meta.data$orig.ident), file = "cellcount_postfiltered.csv")

# Save the object
saveRDS(ild, file = "/scratch/hnatri/ILD/Seurat_objects/prefiltered.rds")
saveRDS(ild2, file = "/scratch/hnatri/ILD/Seurat_objects/postfiltered.rds")

# Which samples were filtered out?
setdiff(ild@meta.data$orig.ident, ild2@meta.data$orig.ident)
# "F01510" "F01869" "F01481" "F01480" "F01635"

ild_all <- ild2

#======================================
# Add metadata
# ======================================

ild_all@meta.data$orig.ident <- gsub("VU", "", ild_all@meta.data$orig.ident)
ild_all@meta.data$orig.ident <- gsub("T", "", ild_all@meta.data$orig.ident)

ipf_metadata <- ipf_metadata[which(ipf_metadata$Library_ID %in% ild_all@meta.data$orig.ident),]

# Count Sample_Names that have multiple samples
ipf_metadata[ipf_metadata$Sample_Name %in% unique(ipf_metadata$Sample_Name[duplicated(ipf_metadata$Sample_Name)]),]
length(unique(ipf_metadata[ipf_metadata$Sample_Name %in% unique(ipf_metadata$Sample_Name[duplicated(ipf_metadata$Sample_Name)]),]$Sample_Name))

# How many patients with IPF?
unique(ipf_metadata$Clinical_Diagnosis)
length(unique(ipf_metadata[which(ipf_metadata$Clinical_Diagnosis=="IPF"),]$Sample_Name))

ctrl_inds <- c("F00409",
               "F01157",
               "F01174",
               "F01365",
               "F01366",
               "F01367",
               "F01394",
               "F01403",
               "F01483",
               "F01494",
               "F01508",
               "F01511",
               "F01513",
               "F01607",
               "F01639",
               "F01641",
               "F01851",
               "F01853",
               "F01874",
               "F02509",
               "F02522",
               "F02524",
               "F02526",
               "F02528",
               "F02607",
               "F02609",
               "F02611",
               "F02613",
               "VUHD103",
               "VUHD65",
               "VUHD66",
               "VUHD67",
               "VUHD68",
               "VUHD70",
               "VUILD59-1",
               "VUILD59-2",
               "VUHD74",
               "VUHD76",
               "VUHD77",
               "VUHD87",
               "VUHD101",
               "THD0010",
               "THD014",
               "THD016")

ctrl_inds <- gsub("VU", "", ctrl_inds)
ctrl_inds <- gsub("T", "", ctrl_inds)

# Match metadata to seurat object
ild_all@meta.data$Diagnosis <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                  from = ipf_metadata$Library_ID,
                                                  to = as.character(ipf_metadata$Clinical_Diagnosis))
ild_all@meta.data$Sample_Name <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                    from = ipf_metadata$Library_ID,
                                                    to = as.character(ipf_metadata$Sample_Name))
ild_all@meta.data$Sample_Source <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                      from = ipf_metadata$Library_ID,
                                                      to = as.character(ipf_metadata$Sample_Source))
ild_all@meta.data$Sample_Type <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                    from = ipf_metadata$Library_ID,
                                                    to = as.character(ipf_metadata$Sample_Type))
ild_all@meta.data$Tobacco <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                from = ipf_metadata$Library_ID,
                                                to = as.character(ipf_metadata$Tobacco))
ild_all@meta.data$Status <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                               from = ipf_metadata$Library_ID,
                                               to = as.character(ipf_metadata$Status))
ild_all@meta.data$Flowcell_ID <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                    from = ipf_metadata$Library_ID,
                                                    to = as.character(ipf_metadata$Flowcell_ID))
ild_all@meta.data$Processing_site <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                        from = ipf_metadata$Library_ID,
                                                        to = as.character(ipf_metadata$Processing_site))
ild_all@meta.data$Gender <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                               from = ipf_metadata$Library_ID,
                                               to = as.character(ipf_metadata$Gender))
ild_all@meta.data$Age <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                            from = ipf_metadata$Library_ID,
                                            to = as.character(ipf_metadata$Age))
ild_all@meta.data$Ethnicity <- plyr::mapvalues(x = ild_all@meta.data$Library_ID,
                                                  from = ipf_metadata$Library_ID,
                                                  to = as.character(ipf_metadata$Ethnicity))

ild_all@meta.data$Sample_Type <- gsub("Control\n", "Control", ild_all@meta.data$Sample_Type)

saveRDS(ild_all, file = "/scratch/hnatri/ILD/Seurat_objects/ILD_all_metadata.rds")

# Extract metadata
md <- ild_all@meta.data %>% as.data.table

# Count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md[, .N, by = c("orig.ident")]
hist(md[, .N, by = c("orig.ident")]$N, breaks = 50)
head(sort(md[, .N, by = c("orig.ident")]$N))

head(ild_all@meta.data)

unique(sapply(strsplit(batchids,split="_") , "[[", 8))
setdiff(unique(sapply(strsplit(batchids,split="_") , "[[", 8)), unique(ild_all@meta.data$orig.ident))

# Calculating cell cycle scores
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ild_all.cc <- CellCycleScoring(ild_all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ild_all <- ild_all.cc

# Visualize the distribution of cell cycle markers across
RidgePlot(ild_all.cc, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

saveRDS(ild_all, file = "/scratch/hnatri/ILD/Seurat_objects/ILD_all_metadata_ccscore.rds")

# =====================================
# Batch integration using rPCA
# =====================================

ild_all <- readRDS("/scratch/hnatri/ILD/Seurat_objects/ILD_merged_DoubletFinder_allLibs.rds")
ild_all.list <- SplitObject(object = ild_all, split.by = "Flowcell_ID")
counter <- 0
ild_all.list <- lapply(X = ild_all.list, FUN = function(x) {
  counter <<- counter + 1
  message(counter)
  x <- NormalizeData(x, verbose = T)
  x <- FindVariableFeatures(x, verbose = T)
})

features <- SelectIntegrationFeatures(object.list = ild_all.list)
ild_all.list <- lapply(X = ild_all.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})

# Identify anchors and integrate datasets using rPCA
anchors <- FindIntegrationAnchors(object.list = ild_all.list,
                                  reference = c(6, 12, 18, 24),
                                  reduction = "rpca",
                                  anchor.features = 3000,
                                  dims = 1:30,
                                  #k.anchor = 20,
                                  verbose = T)

ild_all.integrated <- IntegrateData(anchorset = anchors,
                                    dims = 1:30,
                                    #features.to.integrate = all_features,
                                    verbose = T)

# TODO: 3k variable features
ild_all.integrated <- FindVariableFeatures(ild_all.integrated,
                                           nfeatures = 3000)
ild_all.integrated.features <- VariableFeatures(ild_all.integrated)
ild_all.integrated <- ScaleData(ild_all.integrated,
                                verbose = T)
ild_all.integrated <- RunPCA(ild_all.integrated,
                             reduction.name = "integrated_pca",
                             features = ild_all.integrated.features,
                             verbose = T)

# Finding the numbers of PCs for UMAP
best_pcs_ild_all.integrated <- get_pcs(ild_all.integrated, reduction_name = "integrated_pca")
# Constructing UMAP
# n.neighbors = 50, min.dist = 0.5, spread = 1, 
ild_all.integrated <- RunUMAP(ild_all.integrated,
                              dims = 1:17,
                              reduction = "integrated_pca",
                              reduction.name = "integrated_umap",
                              verbose = T)
ild_all.integrated <- FindNeighbors(ild_all.integrated,
                                    dims = 1:17,
                                    reduction = "integrated_pca")
ild_all.integrated.clusters <- FindClusters(object = ild_all.integrated,
                                            resolution = c(0.01, 0.025, 0.05,
                                                           0.10, 0.40, 0.80),
                                            graph.name = "integrated_snn")

ild_all <- ild_all.integrated.clusters
ild_all@meta.data$seurat_clusters <- ild_all@meta.data$integrated_snn_res.0.025
Idents(ild_all) <- "integrated_snn_res.0.025"

# Annotation for the four cell populations:
# PTPRC+ (immune cells), EPCAM+ (epithelial cells), PECAM1+/PTPRC− 
# (endothelial cells), and PTPRC−/EPCAM−/PECAM1− (mesenchymal cells)
DimPlot(ild_all, group.by = "seurat_clusters")
DefaultAssay(ild_all) <- "RNA"
DotPlot(ild_all, features = c("EPCAM","PECAM1","PTPRC","LUM","DCN"))
FeaturePlot(ild_all, c("EPCAM","PECAM1","PTPRC","LUM","DCN"))

# Cluster 6 looks like doublets. Proportions of doublets in each cluster:
head(ild_all@meta.data)
doublet_summary <- ild_all@meta.data %>% group_by(seurat_clusters, doublet_finder) %>%
  summarize(n = n())

doublet_summary <- pivot_wider(doublet_summary, names_from = doublet_finder, values_from = n) 
doublet_summary$doublet_prop <- doublet_summary$Doublet/(doublet_summary$Doublet+doublet_summary$Singlet)
doublet_summary <- arrange(doublet_summary, desc(doublet_prop))

cellpops <- as.character(ild_all@meta.data$seurat_clusters)
cellpops[cellpops %in% c(0,4,8)] <- "Immune"
cellpops[cellpops %in% c(3,7)] <- "Endothelial"
cellpops[cellpops %in% c(1,2)] <- "Epithelial"
cellpops[cellpops %in% c(5)] <- "Mesenchymal"

ild_all@meta.data$population <- cellpops

doublets <- as.character(ild_all@meta.data$seurat_clusters)
doublets <- ifelse(doublets==6, "doublets", "singlets")

ild_all@meta.data$doublets <- doublets
DimPlot(ild_all, group.by = "doublets")

# Adding the index type to the metadata
doubleindex_libs <- readLines("/home/hnatri/ILD_processing_annotation/2021_7_13_12_48_52_lpeter_Sample_List.tsv")

libindex <- ild_all@meta.data$orig.ident
libindex <- ifelse(libindex %in% doubleindex_libs, "double", "single")
ild_all@meta.data$lib_index_type <- libindex

# Checking overlap with the genotype data
gt_samples <- readLines("/home/hnatri/ILD_eQTL/Data/gencove_sample_n132.tsv")
length(unique(ild_all@meta.data$Sample_Name))
scrnaseq_samples <- ild_all@meta.data$Sample_Name
scrnaseq_samples <- gsub("VU", "", scrnaseq_samples)
scrnaseq_samples <- gsub("TI", "", scrnaseq_samples)

scrnaseq_samples <- gsub("HD101 \\(We wrote HD103 on day of processing\\)", "HD101", scrnaseq_samples)
scrnaseq_samples <- gsub("ILD92 More" , "ILD92", scrnaseq_samples)
scrnaseq_samples <- gsub("HD84", "HD084", scrnaseq_samples)
scrnaseq_samples <- gsub("HD85", "HD085", scrnaseq_samples)
scrnaseq_samples <- gsub("HD92", "HD092", scrnaseq_samples)
scrnaseq_samples <- gsub("HD94", "HD094", scrnaseq_samples)
scrnaseq_samples <- gsub("HD95", "HD095", scrnaseq_samples)
scrnaseq_samples <- gsub("HD98", "HD098", scrnaseq_samples)
scrnaseq_samples <- gsub("HD69", "HD069", scrnaseq_samples)
scrnaseq_samples <- gsub("HD65", "HD065", scrnaseq_samples)
scrnaseq_samples <- gsub("HD66", "HD066", scrnaseq_samples)
scrnaseq_samples <- gsub("HD67", "HD067", scrnaseq_samples)
scrnaseq_samples <- gsub("HD68", "HD068", scrnaseq_samples)
scrnaseq_samples <- gsub("HD70", "HD070", scrnaseq_samples)
scrnaseq_samples <- gsub("HD87", "HD087", scrnaseq_samples)

ild_all@meta.data$clean_Sample_Name <- scrnaseq_samples
scrnaseq_samples <- unique(scrnaseq_samples)

gt_samples <- gsub("VU", "", gt_samples)
gt_samples <- gsub("TI", "", gt_samples)
length(intersect(scrnaseq_samples, gt_samples))
setdiff(scrnaseq_samples, gt_samples)
setdiff(gt_samples, scrnaseq_samples)

# Adding disease status to TILD103
ild_all@meta.data[which(ild_all@meta.data$Sample_Name=="TILD103"),]$Status <- "Disease"

# Saving the object
saveRDS(ild_all, "/scratch/hnatri/ILD/Seurat_objects/ild_all.6_12_18_24_integrated.clusters.4pop.210704.rds")
