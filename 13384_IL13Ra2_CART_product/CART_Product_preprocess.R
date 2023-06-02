#==============================================================================#
# Original author(s) : Stephanie L. Yahn,
#                      Heini M. Natri hnatri@tgen.org
# Date: 2021/11/30
# Description: Preprocessing of leuk PBMC and Product samples for the CAR T
# project
#==============================================================================#

#==============================================================================#
# Load libraries
#==============================================================================#

#remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)
library(stringr)
library(readxl)
library(glmGamPoi)
library(gridExtra)

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)

#==============================================================================#
# Batch 1 preprocessing
# Processed separately because it relies on Demuxlet for demultiplexing
#==============================================================================#

# Read in output from CellRanger (2 data types: Gene Expression and Antibody Capture)
# Austin is working on moving the data to /labs/, currently batches still split
# in three locations within /labs/.
# This reads the data in as a list, because it's multimodal
beet1 = Read10X(data.dir = "/labs/banovich/BCTCSF/Outs/TGen_Ref_F02323-GEX_F02325-AB/outs/filtered_feature_bc_matrix/")

# Simplify antibody names
# 189 antibodies in this pool
# 3 different antibody panels, some instances where the number differs, e.g. multiple
# CD4s with different antibody name
rownames(beet1$`Antibody Capture`)
rownames(beet1$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(beet1$`Antibody Capture`)), ' ')
rownames(beet1$`Antibody Capture`) = str_remove(rownames(beet1$`Antibody Capture`), "mouse_")
rownames(beet1$`Antibody Capture`) = str_remove(rownames(beet1$`Antibody Capture`), "rat_")
rownames(beet1$`Antibody Capture`) = str_remove(rownames(beet1$`Antibody Capture`), "human_")

# Set up the Seurat object
batch1 = CreateSeuratObject(counts = beet1[["Gene Expression"]])
# Calculating the % of MT reads
batch1 = PercentageFeatureSet(batch1, pattern = "^MT-", col.name = "percent.mt")
# Adding a protein assay based on the antibody matrix. Underscores get replaced
# with dashes.
batch1[["Protein"]] = CreateAssayObject(beet1[["Antibody Capture"]][, colnames(x = batch1)])

# Add demuxlet metadata to assign patient identities to each cell.
# read in the output from demuxlet
demuxlet = fread(file ="/labs/banovich/BCTCSF/Stephanie/Demuxlet/Batch1/1_TGen_Ref_F02323-GEX_F02325-AB_4xexomes.best")

# Obtain the patient ID and the assignment from the BEST column. IDs are based
# on our exome sequencing data.

# Assignment: doublet or singlet
assignment = sapply(demuxlet$BEST, function(x) {
  strsplit(x,"-")[[1]][[1]]
})

# Patient ID
patient_id = sapply(demuxlet$BEST, function(x) {
  strsplit(x,"-")[[1]][[2]]
})

# Adding the new metadata to the object
CellsMeta = batch1@meta.data
CellsMeta["Exome_Sample_Name"] = patient_id
CellsMeta["Demultiplex_Assignment"] = assignment
batch1 = AddMetaData(batch1, CellsMeta)
table(batch1$Demultiplex_Assignment)

# Add sample metadata
# UPN = unique patient number
# Cycle = which infusion
# Day = how many days after the infusion
# Using the exome sample names to map over the new metadata features
b1.meta.data = read_excel("/labs/banovich/BCTCSF/Stephanie/Batches_metadata_forR.xlsx", sheet = 1)
b1.meta.data = b1.meta.data %>% drop_na(UPN)

batch1@meta.data$UPN = plyr::mapvalues(x = batch1@meta.data$Exome_Sample_Name,
                                       from = b1.meta.data$Exome_Sample_Name,
                                       to = as.character(b1.meta.data$UPN))

batch1@meta.data$Sample_Type = plyr::mapvalues(x = batch1@meta.data$UPN,
                                               from = b1.meta.data$UPN,
                                               to = as.character(b1.meta.data$Sample_Type))

batch1@meta.data$Cycle = plyr::mapvalues(x = batch1@meta.data$UPN,
                                         from = b1.meta.data$UPN,
                                         to = as.character(b1.meta.data$Cycle))

batch1@meta.data$Day = plyr::mapvalues(x = batch1@meta.data$UPN,
                                       from = b1.meta.data$UPN,
                                       to = as.character(b1.meta.data$Day))

batch1@meta.data$Manufacture = plyr::mapvalues(x = batch1@meta.data$UPN,
                                               from = b1.meta.data$UPN,
                                               to = as.character(b1.meta.data$Manufacture))

# Used later to keep the leukapheresis PBMCs
batch1$Cycle_Day = paste("Cycle", batch1$Cycle, "_Day", batch1$Day, sep = "")

# Adding batch
batch1@meta.data["Batch"] = "Batch1"

# keep only singlets
Idents(batch1)
Idents(batch1) = batch1$Demultiplex_Assignment
batch1_filtered = subset(batch1, idents = "SNG")
length(colnames(batch1_filtered))
length(colnames(batch1))

# Remove unwanted cells/keep wanted cells
# Keeping high quality, singlets, sample types that we want
# subset = nFeature_RNA > 500 & percent.mt < 25 & nCount_Protein < 10000
VlnPlot(batch1_filtered, features = c("nFeature_RNA", "nCount_RNA", "nFeature_Protein", "nCount_Protein", "percent.mt"), ncol = 3)
summary(batch1_filtered$nCount_Protein)
hist(batch1_filtered$nCount_Protein)
plot(ecdf(batch1_filtered$nCount_RNA))

batch1_filtered = subset(batch1_filtered, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10 & nCount_Protein < 10000)
length(colnames(batch1_filtered))
length(colnames(batch1))

# remove CSF and Tumor samples
batch1_filtered = subset(batch1_filtered, 
                         cells = row.names(batch1_filtered@meta.data[!(batch1_filtered@meta.data$Sample_Type %in% c("Tumor", "CSF")), ]))
table(batch1_filtered$Sample_Type)
# Currently has leukapheresis PBMCs and PBMCs matched to CSF
table(batch1_filtered$Manufacture)
# keep only leuk_PBMC and Product samples, remove PBMCs matched to CSF
Idents(batch1_filtered) = batch1_filtered$Cycle_Day
batch1_filtered = subset(batch1_filtered, idents = "CycleNA_DayNA")
table(batch1_filtered$Manufacture)

# Normalize and scale protein data
# due to the unspecific binding background signal, log-normalization doesn't work well for CITEseq protein data
# instead, Seurat recommends centered log-ratio (CLR) normalization computed independently for each feature
batch1_filtered = NormalizeData(batch1_filtered, assay = "Protein", normalization.method = "CLR")
batch1_filtered = ScaleData(batch1_filtered, assay = "Protein")

rm(CellsMeta, demuxlet, b1.meta.data)

#==============================================================================#
# Batches 2-19 preprocessing
#==============================================================================#

# Using the hashtag oligos instead of demuxlet

# Read in output from CellRanger (2 data types: Gene Expression and Antibody Capture)
beet2 = Read10X("/labs/banovich/BCTCSF/Outs/F02397-GEX_F02398-AB/outs/filtered_feature_bc_matrix/")
beet3 = Read10X("/labs/banovich/BCTCSF/Outs/F02574-GEX_F02575-AB/outs/filtered_feature_bc_matrix/")
beet4 = Read10X("/labs/banovich/BCTCSF/Outs/F02588-GEX_F02589-AB/outs/filtered_feature_bc_matrix/")

beet5 = Read10X("/labs/banovich/BCTCSF/Outs/IL13OP/IL13OP_F02660-GEX_F02661-AB/outs/filtered_feature_bc_matrix/")
beet6 = Read10X("/labs/banovich/BCTCSF/Outs/IL13OP/IL13OP_F02677-GEX_F02678-AB/outs/filtered_feature_bc_matrix/")

data_loc = "/labs/banovich/SingleCell/CellRanger/5_0_0/Ensemble_98/PipelineData/Projects/BCTCSF/CellRangerOuts/CITE/"
beet7 = Read10X(paste0(data_loc, "IL13OP_F02736-GEX_F02755-FB/outs/filtered_feature_bc_matrix/"))
beet8 = Read10X(paste0(data_loc, "IL13OP_F02981-GEX_F02982-FB/outs/filtered_feature_bc_matrix/"))
beet9 = Read10X(paste0(data_loc, "IL13OP_F02985-GEX_F02986-FB/outs/filtered_feature_bc_matrix/"))
beet10 = Read10X(paste0(data_loc, "IL13OP_F03246-GEX_F03247-FB/outs/filtered_feature_bc_matrix/"))
beet11 = Read10X(paste0(data_loc, "IL13OP_F03254-GEX_F03256-FB/outs/filtered_feature_bc_matrix/"))
beet12 = Read10X(paste0(data_loc, "IL13OP_F03253-GEX_F03255-FB/outs/filtered_feature_bc_matrix/"))
beet13 = Read10X(paste0(data_loc, "IL13OP_F03263-GEX_F03265-FB/outs/filtered_feature_bc_matrix/"))
beet14 = Read10X(paste0(data_loc, "IL13OP_F03264-GEX_F03266-FB/F03264-GEX_F03266-FB/outs/filtered_feature_bc_matrix/"))
beet15 = Read10X(paste0(data_loc, "IL13OP_F03285-GEX_F03286-FB/outs/filtered_feature_bc_matrix/"))
beet16 = Read10X(paste0(data_loc, "IL13OP_F03287-GEX_F03288-FB/outs/filtered_feature_bc_matrix/"))
beet17 = Read10X(paste0(data_loc, "IL13OP_F03289-GEX_F03290-FB/F03289-GEX_F03290-FB/outs/filtered_feature_bc_matrix/"))
beet18 = Read10X(paste0(data_loc, "IL13OP_F03291-GEX_F03292-FB/F03291-GEX_F03292-FB/outs/filtered_feature_bc_matrix/"))
beet19 = Read10X(paste0(data_loc, "IL13OP_F03293-GEX_F03294-FB/outs/filtered_feature_bc_matrix/"))

# Simplify antibody names
# Creating a list of objects
bt.list = ls(pattern="beet")
bt.list = str_sort(bt.list, numeric = TRUE)
bt.list = do.call("list", mget(bt.list))

# Looping through all the objects
for (i in 1:length(bt.list)) {
  rownames(bt.list[[i]]$`Antibody Capture`) = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', rownames(bt.list[[i]]$`Antibody Capture`)), ' ')
  rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "mouse_")
  rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "rat_")
  rownames(bt.list[[i]]$`Antibody Capture`) = str_remove(rownames(bt.list[[i]]$`Antibody Capture`), "human_")
}

# Within the antibody capture, there are CITE-seq antibodies and the hashing
# antibodies, need to separate
# Cell hashing antibody names
hash_antibodies = c("TotalSeqC0251_Hashtag1", 
                    "TotalSeqC0252_Hashtag2",  
                    "TotalSeqC0253_Hashtag3", 
                    "TotalSeqC0254_Hashtag4",
                    "TotalSeqC0256_Hashtag6",
                    "TotalSeqC0257_Hashtag7", 
                    "TotalSeqC0258_Hashtag8",
                    "TotalSeqC0259_Hashtag9")

# Not all CITEseq antibodies were used for CITEseq batches 14, 17, and 18, 
# but cellranger was run as if they were, so need to filter the antibodies for these batches
# to remove false positive counts
# all_abs sheet (from Google Drive > feature_barcoding)
all_abs = read.csv(file = "/home/hnatri/CART/8plex_feature_reference - BCTCSF_FR.csv")
all_abs$ab_name = strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\2', all_abs$name), ' ')
all_abs$ab_name = as.character(all_abs$ab_name)
all_abs$ab_name = str_remove(all_abs$ab_name, "mouse_")
all_abs$ab_name = str_remove(all_abs$ab_name, "rat_")
all_abs$ab_name = str_remove(all_abs$ab_name, "human_")

# panel_137 sheet (from Google Drive > feature_barcoding)
panel_137 = read.csv(file = "/home/hnatri/CART/TS-C human panel_137_Antibodies (#399905).xlsx - BarcodeList.csv")
# all_abs sheet has ab names in the correct format to match matrix, panel_137 does not
panel_137 = subset(all_abs, id %in% panel_137$DNA_ID)

# filter antibodies in batches 14, 17, and 18
keep_prots = c(panel_137$ab_name, hash_antibodies)
bt.list$beet14$`Antibody Capture` = bt.list$beet14$`Antibody Capture`[rownames(bt.list$beet14$`Antibody Capture`) %in% keep_prots,]
bt.list$beet17$`Antibody Capture` = bt.list$beet17$`Antibody Capture`[rownames(bt.list$beet17$`Antibody Capture`) %in% keep_prots,]
bt.list$beet18$`Antibody Capture` = bt.list$beet18$`Antibody Capture`[rownames(bt.list$beet18$`Antibody Capture`) %in% keep_prots,]

# Iterating through the beet list to create a batch list and a filtered object
# list
batch.list = list()
batch.list_filtered = list()
for (i in 2:length(bt.list)) { #excluding Batch 1
  #i <- 2
  message(names(bt.list)[i])
  # Set up the Seurat object
  # Splitting out the gene expression and the antibody matrices
  GEX = bt.list[[i]][[1]]
  AB = bt.list[[i]][[2]]
  # Creating the object with the expression matrix
  batch.list[[i]] = CreateSeuratObject(counts = GEX)
  batch.list[[i]] = PercentageFeatureSet(batch.list[[i]], pattern = "^MT-", col.name = "percent.mt")
  
  hash = AB[rownames(AB) %in% hash_antibodies,]
  citeseq = AB[!(rownames(AB) %in% hash_antibodies),]
  
  # Creating the protein and hash assays
  batch.list[[i]][["Protein"]] = CreateAssayObject(counts = citeseq)
  batch.list[[i]][["Hash"]] = CreateAssayObject(counts = hash)
}

# Smooth-scatter plot of MT reads and RNA counts
batch.list2 <- batch.list
batch.list2[[1]] = batch1
bt_merge <- merge(x = batch.list2[[1]], y = batch.list2[2:length(batch.list2)])
#bt_merge <- PercentageFeatureSet(object = bt_merge, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(bt_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Make density plots
# subset = nFeature_RNA > 500 & nCount_RNA > 5000 & percent.mt < 10 & nCount_Protein < 10000)
pdf("/home/hnatri/CART/QC_plot_nFeature.pdf")
smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nFeature_RNA)
abline(h = 500, v = 10)
text(1.5,700, "nFeature_RNA = 500, percent.mt = 10", adj = c(0, -.1))
dev.off()

pdf("/home/hnatri/CART/QC_plot_nCount.pdf")
smoothScatter(bt_merge@meta.data$percent.mt, bt_merge@meta.data$nCount_RNA)
abline(h = 1000, v = 10)
text(1.5,1200, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))
dev.off()

smoothScatter(bt_merge@meta.data$percent.mt, log(bt_merge@meta.data$nCount_RNA))
abline(h = log(1000), v = 10)
text(1.5,700, "nCount_RNA = 1000, percent.mt = 10", adj = c(0, -.1))
#dev.off()

# Normalizing and filtering
for (i in 2:length(batch.list)) { 
  # Normalizing and scaling
  batch.list[[i]] = NormalizeData(batch.list[[i]], assay = "Protein", normalization.method = "CLR")
  batch.list[[i]] = ScaleData(batch.list[[i]], assay = "Protein")
  batch.list[[i]] = NormalizeData(batch.list[[i]], assay = "Hash", normalization.method = "CLR")
  batch.list[[i]] = ScaleData(batch.list[[i]], assay = "Hash")
  
  # Assign single cells back to their sample origins
  # Demultiplexing based on the hashing antibodies
  # Singlets kept based on "hash classification global"
  batch.list[[i]] = HTODemux(batch.list[[i]], assay = "Hash", positive.quantile = 0.99, verbose = F)
  
  # Add sample metadata
  meta.data = read_excel("/labs/banovich/BCTCSF/Stephanie/Batches_metadata_forR.xlsx", sheet = i)
  meta.data = meta.data %>% drop_na(UPN)
  batch.list[[i]]@meta.data$UPN = plyr::mapvalues(x = batch.list[[i]]@meta.data$hash.ID,
                                                  from = meta.data$CellHashing_Ab,
                                                  to = as.character(meta.data$UPN))
  
  batch.list[[i]]@meta.data$Sample_Type = plyr::mapvalues(x = batch.list[[i]]@meta.data$UPN,
                                                          from = meta.data$UPN,
                                                          to = as.character(meta.data$Sample_Type))
  
  batch.list[[i]]@meta.data$Cycle = plyr::mapvalues(x = batch.list[[i]]@meta.data$UPN,
                                                    from = meta.data$UPN,
                                                    to = as.character(meta.data$Cycle))
  
  batch.list[[i]]@meta.data$Day = plyr::mapvalues(x = batch.list[[i]]@meta.data$UPN,
                                                  from = meta.data$UPN,
                                                  to = as.character(meta.data$Day))
  batch.list[[i]]@meta.data$Manufacture = plyr::mapvalues(x = batch.list[[i]]@meta.data$UPN,
                                                          from = meta.data$UPN,
                                                          to = as.character(meta.data$Manufacture))
  
  batch.list[[i]]$Cycle_Day = paste("Cycle", batch.list[[i]]$Cycle, "_Day", batch.list[[i]]$Day, sep = "")
  
  batch.list[[i]]@meta.data["Batch"] = paste("Batch", i, sep = "")
  
  # Keeping singlets only
  Idents(batch.list[[i]]) = batch.list[[i]]$Hash_classification.global
  # VlnPlot(batch.list[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_Protein", "nFeature_Protein", "nFeature_Hash"), pt.size = 0)
  batch.list_filtered[[i]] = subset(batch.list[[i]], idents = "Singlet")
}

# Filtering for features, counts, and proportion of MT reads
# Some outliers in nCounts_Protein >10,000
for (i in 2:length(batch.list_filtered)) { 
  batch.list_filtered[[i]] = subset(batch.list_filtered[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10 & nCount_Protein < 10000)
}

for (i in 2:length(batch.list_filtered)) { 
  # remove CSF and Tumor samples
  batch.list_filtered[[i]] = subset(batch.list_filtered[[i]], 
                                    cells = row.names(batch.list_filtered[[i]]@meta.data[!(batch.list_filtered[[i]]@meta.data$Sample_Type %in% c("Tumor", "CSF")), ]))
  
  # keep only leuk_PBMC and Product samples
  Idents(batch.list_filtered[[i]]) = batch.list_filtered[[i]]$Cycle_Day
  batch.list_filtered[[i]] = subset(batch.list_filtered[[i]], idents = "CycleNA_DayNA")
  
  # Scale protein/antibody data (must re-scale after subsetting)
  batch.list_filtered[[i]] = ScaleData(batch.list_filtered[[i]], assay = "Protein")
  batch.list_filtered[[i]] = ScaleData(batch.list_filtered[[i]], assay = "Hash")
}

# Remove cells erroneously identified as “TotalSeqC0259-Hashtag9” in batch11
# and batch15 as these batches did not include this hashing antibody
batch.list_filtered[[11]] = subset(batch.list_filtered[[11]], cells = row.names(batch.list_filtered[[11]]@meta.data[!(batch.list_filtered[[11]]@meta.data$UPN %in% "TotalSeqC0259-Hashtag9"), ]))
batch.list_filtered[[15]] = subset(batch.list_filtered[[15]], cells = row.names(batch.list_filtered[[15]]@meta.data[!(batch.list_filtered[[15]]@meta.data$UPN %in% "TotalSeqC0259-Hashtag9"), ]))
# Scale protein/antibody data (must re-scale after subsetting)
batch.list_filtered[[11]] = ScaleData(batch.list_filtered[[11]], assay = "Protein")
batch.list_filtered[[11]] = ScaleData(batch.list_filtered[[11]], assay = "Hash")
batch.list_filtered[[15]] = ScaleData(batch.list_filtered[[15]], assay = "Protein")
batch.list_filtered[[15]] = ScaleData(batch.list_filtered[[15]], assay = "Hash")

# Add batch 1 to lists
batch.list[[1]] = batch1
batch.list_filtered[[1]] = batch1_filtered

# Name list elements
batch_names = paste("Batch", seq(1,19), sep = "")
names(batch.list) = batch_names
batch_filt_names = paste(batch_names, "_filtered", sep = "")
names(batch.list_filtered) = batch_filt_names

# # cells in each batch
lapply(names(batch.list), function(xx){
  message(xx, " ", nrow(batch.list[[xx]]@meta.data))
})

lapply(names(batch.list_filtered), function(xx){
  message(xx, " ", nrow(batch.list_filtered[[xx]]@meta.data))
})

# Remove Protein assay from non-CITEseq batches
nonCITE = c("Batch7_filtered", "Batch8_filtered", "Batch9_filtered", "Batch10_filtered", "Batch11_filtered", 
            "Batch12_filtered", "Batch13_filtered", "Batch15_filtered", "Batch16_filtered", "Batch19_filtered")
for (i in nonCITE) {
  batch.list_filtered[[i]]$Protein = NULL
}

# SCTransform normalization of RNA counts (replaces NormalizeData(),
# ScaleData(), and FindVariableFeatures()). SCTransform also supports using
# glmGamPoi package which substantially improves the speed of the learning
# procedure
for (i in 1:length(batch.list_filtered)) {
  DefaultAssay(batch.list_filtered[[i]]) = "RNA"
  batch.list_filtered[[i]] = SCTransform(batch.list_filtered[[i]],
                                         method = "glmGamPoi",
                                         vars.to.regress = c("percent.mt"),
                                         verbose = T)
}

# Add cell cycle score
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#batch.list_filtered[["Batch19_filtered"]] <- NormalizeData(batch.list_filtered[["Batch19_filtered"]])

for(i in 1:length(batch.list_filtered)){
  message(names(batch.list_filtered)[i])
  DefaultAssay(batch.list_filtered[[i]]) <- "SCT"
  batch.list_filtered[[i]] <- CellCycleScoring(batch.list_filtered[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = F) 
}

# Renormalizing with cell cycle scores
for (i in 1:length(batch.list_filtered)) {
  DefaultAssay(batch.list_filtered[[i]]) = "RNA"
  batch.list_filtered[[i]] = SCTransform(batch.list_filtered[[i]],
                                         method = "glmGamPoi",
                                         vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                         verbose = T)
}

saveRDS(batch.list_filtered, file = "/scratch/hnatri/CART/batch_list_filtered_211207.rds")

