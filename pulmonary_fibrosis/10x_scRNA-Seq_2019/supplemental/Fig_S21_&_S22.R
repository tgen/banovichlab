# ==============================================================================
# GSE128033 Processing
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
setwd("/scratch/agutierrez/IPF/R/Misc/GSE128033")
main_dir <- getwd()
current_date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, current_date), showWarnings = FALSE)
setwd(file.path(main_dir, current_date))

getwd()

set.seed(2811)
# ======================================
# Load libraries
# ======================================
library(Seurat)
library(ggplot2)

# ======================================
# Define the location of the Cellranger output files to be read into R
# Define the location of barcodes.tsv, genes.tsv, matrix.mtx
# Create a character vector of sample file names
# Read data from each sample folder and combine
# ======================================
dataset_loc <- "/scratch/agutierrez/IPF/R/Misc/GSE128033/Data/"
batchids <- list.files(path = dataset_loc, pattern = "GSM*",
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE)

d10x.data <- sapply(batchids,  function(i){
  sample_id <- sapply(strsplit(i,split="_"), "[[", 1)
  d10x <- Read10X(file.path(dataset_loc,i))
  colnames(d10x) <- paste(sample_id, sep="_", colnames(d10x))
  d10x
})

# ======================================
# do.call is required to create a single S4 object of class dgCMatrix
# (genes x cells) otherwise, cbind will create a 2x1 matrix list of two 
# separate S4 objects of class dgCMatrix (genes x cells)
# ======================================
#experiment.data <- do.call("cbind", d10x.data)
n <- max(sapply(d10x.data, nrow)) 
experiment.data <- do.call(cbind, lapply(d10x.data, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 

# Save the dgCMatrix
writeMM(experiment.data, file = "matrix.mtx")
write(x = rownames(experiment.data), file = "genes.tsv")
write(x = colnames(experiment.data), file = "barcodes.tsv")

# ======================================
# Create Seurat Object from the dgCMatrix
# Read in the dgCMatrix file:
# ======================================
ild_data <- Read10X(data.dir = ".", gene.column=1)

# ======================================
# Keep all cells with at least 200 detected genes,
# Use unique sample name from the 8th field of initial
# input identity (from batchids)
# ======================================
ild <- CreateSeuratObject(
  ild_data,
  project = "scRNA lung",
  min.features = 200,
  names.field = 1,
  names.delim = "_")

# Free up some memory
rm(d10x.data)
rm(experiment.data)
rm(ild_data)

# save file
saveRDS(ild, file = "GSE128033_raw.rds")

# ======================================
# QC Seurat object
# ======================================
ild <- PercentageFeatureSet(object = ild, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(ild, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Make density plots
smoothScatter(ild@meta.data$percent.mt, ild@meta.data$nFeature_RNA)
abline(h = 1000, v = 1.5)
text(1.5,700, "nFeature = 1000, percent.mt = 1.5", adj = c(0, -.1))

# Filter out cells with less than 1000 nFeature and more than 25% percent.mt
ild <- subset(ild, subset = nFeature_RNA > 1000 & percent.mt < 25)

# ======================================
# Seurat object processing
# ======================================  
# Run SCTransform
ild <- SCTransform(object = ild, verbose = T, batch_var = "orig.ident")

# Chose a PC for analysis
ild <- FindVariableFeatures(ild, verbose = T, nfeatures = 3000)
ild <- ScaleData(ild, features = row.names(ild@assays$SCT@data))
ild <- RunPCA(ild)
ElbowPlot(ild)
ild <- RunUMAP(object = ild, dims = 1:10, verbose = F)
DimPlot(ild)
saveRDS(ild, file = "GSE128033_ILD_sct.rds")

# ======================================
# Transfer Anchors 
# ======================================
reference <- readRDS("/scratch/agutierrez/IPF/R/Seurat/Reference/2019_Release_IPF.rds")

anchors <- FindTransferAnchors(reference = reference,
                               query = ild,
                               dims = 1:10,
                               npcs = NULL)

predictions <- TransferData(anchorset = anchors,
                            refdata = reference@meta.data$celltype,
                            weight.reduction = "pca",
                            dims = 1:10,
                            verbose = T)

ild <- AddMetaData(ild, metadata = predictions)
pred <- prop.table(table(ild@meta.data$predicted.id))
ref <- prop.table(table(reference@meta.data$celltype))
ref <- ref[rownames(ref) %in% rownames(pred)]
ild@meta.data[8:length(colnames(ild@meta.data))] <- NULL
ild@meta.data$Study <- "GSE128033"
GSE128033 <- ild
rm(ild)

comparison <- cbind(pred, ref)
comparison <- (comparison * 100)

write.csv(comparison, "comparison.csv", quote = F)
saveRDS(GSE128033, "GSE128033.rds")

# ==============================================================================
# GSE122960 Processing
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
setwd("/scratch/agutierrez/IPF/R/Misc/GSE122960")
main_dir <- getwd()
current_date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, current_date), showWarnings = FALSE)
setwd(file.path(main_dir, current_date))

getwd()

set.seed(2811)
# ======================================
# Load libraries
# ======================================
library(Seurat)
library(RCurl)
library(DoubletFinder)
library(Matrix)
library(ade4)

# ======================================
#
# ======================================
# ======================================
# Define the location of the Cellranger output files to be read into R
# Define the location of barcodes.tsv, genes.tsv, matrix.mtx
# Create a character vector of sample file names
# Read data from each sample folder and combine
# ======================================
dataset_loc <- "/scratch/agutierrez/IPF/R/Misc/GSE122960/Data/"
#inputfiles_loc <- "outs/filtered_feature_bc_matrix"
batchids <- list.files(path = dataset_loc, pattern = "GSM*",
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE,
                       all.files = TRUE)

d10x.data <- sapply(batchids,  function(i){
  sample_id <- sapply(strsplit(i,split="_"), "[[", 1)
  d10x <- Read10X_h5(file.path(dataset_loc, i, "/filtered_gene_bc_matrices_h5.h5"))
  colnames(d10x) <- paste(sample_id, sep="_", colnames(d10x))
  d10x
})

# ======================================
# do.call is required to create a single S4 object of class dgCMatrix
# (genes x cells) otherwise, cbind will create a 2x1 matrix list of two 
# separate S4 objects of class dgCMatrix (genes x cells)
# ======================================
#experiment.data <- do.call("cbind", d10x.data)

n <- max(sapply(d10x.data, nrow)) 
experiment.data <- do.call(cbind, lapply(d10x.data, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 

# Save the dgCMatrix
writeMM(experiment.data, file = "matrix.mtx")
write(x = rownames(experiment.data), file = "genes.tsv")
write(x = colnames(experiment.data), file = "barcodes.tsv")
# ======================================
# Create Seurat Object from the dgCMatrix
# Read in the dgCMatrix file:
# ======================================
ild_data <- Read10X(data.dir = ".", gene.column=1)

# ======================================
# Keep all cells with at least 200 detected genes,
# Use unique sample name from the 8th field of initial
# input identity (from batchids)
# ======================================
ild <- CreateSeuratObject(
  ild_data,
  project = "scRNA lung",
  min.features = 200,
  names.field = 1,
  names.delim = "_")

# Free up some memory
rm(d10x.data)
rm(experiment.data)
rm(ild_data)

# save file
saveRDS(ild, file = "GSE122960_raw.rds")

# ======================================
# QC Seurat object
# ======================================
ild <- PercentageFeatureSet(object = ild, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(ild, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Make density plots
smoothScatter(ild@meta.data$percent.mt, ild@meta.data$nFeature_RNA)
abline(h = 1000, v = 1.5)
text(1.5,700, "nFeature = 1000, percent.mt = 1.5", adj = c(0, -.1))

# Filter out cells with less than 1000 nFeature and more than 25% percent.mt
ild <- subset(ild, subset = nFeature_RNA > 1000 & percent.mt < 25)

# ======================================
# Seurat object processing
# ======================================  
# Run SCTransform
ild <- SCTransform(object = ild, verbose = T, batch_var = "orig.ident")

# Chose a PC for analysis
ild <- FindVariableFeatures(ild, verbose = T, nfeatures = 3000)
ild <- ScaleData(ild, features = row.names(ild@assays$SCT@data))
ild <- RunPCA(ild)
ElbowPlot(ild)
ild <- RunUMAP(object = ild, dims = 1:15, verbose = F)
DimPlot(ild)
saveRDS(ild, file = "GSE122960_ILD_sct.rds")

# ======================================
# Transfer Anchors 
# ======================================
reference <- readRDS("/scratch/agutierrez/IPF/R/Seurat/Reference/2019_Release_IPF.rds")
anchors <- FindTransferAnchors(reference = reference,
                               query = ild,
                               dims = 1:15,
                               npcs = NULL)

predictions <- TransferData(anchorset = anchors,
                            refdata = reference@meta.data$celltype,
                            weight.reduction = "pca",
                            dims = 1:15,
                            verbose = T)

ild <- AddMetaData(ild, metadata = predictions)
pred <- prop.table(table(ild@meta.data$predicted.id))
ref <- prop.table(table(reference@meta.data$celltype))
ref <- ref[rownames(ref) %in% rownames(pred)]
ild@meta.data[8:length(colnames(ild@meta.data))] <- NULL
ild@meta.data$Study <- "GSE122960"
GSE122960 <- ild
rm(ild)

comparison <- cbind(pred, ref)
comparison <- (comparison * 100)

write.csv(comparison, "comparison.csv", quote = F)
saveRDS(GSE122960, "GSE122960.rds")

# ==============================================================================
# Prepare data for plotting
# ==============================================================================
GSE135893 <- reference
GSE135893@meta.data$Study <- "GSE135893"
rm(reference)

GSE122960@meta.data$celltype <- GSE122960@meta.data$predicted.id
GSE122960@meta.data$predicted.id <- NULL
GSE128033@meta.data$celltype <- GSE128033@meta.data$predicted.id
GSE128033@meta.data$predicted.id <- NULL

head(GSE128033@meta.data)
head(GSE122960@meta.data)
head(GSE135893@meta.data)

combined <- merge(x = GSE128033, y=c(GSE122960, GSE135893))
keep <-  c("KRT5-/KRT17+",
           "SCGB3A2+",
           "Transitional AT2",
           "PLIN2+ Fibroblasts")

small <- subset(combined, cells = rownames(combined@meta.data[combined@meta.data$celltype %in% keep, ]))
small <- FindVariableFeatures(small, verbose = T, nfeatures = 3000)
small <- ScaleData(small, features = row.names(small@assays$SCT@data))
small <- RunPCA(small)

Authors <- as.character(small@meta.data$Study)

onion <- Authors
onion[onion == "GSE128033"] <- "Morse et al."
onion[onion == "GSE122960"] <- "Reyfman et al."
onion[onion == "GSE135893"] <- "Habermann et al."
small@meta.data$Authors <- onion

# ======================================
# Figure S: 22
# ======================================
DotPlot(small,
        features = c("AGER",
                     "SFTPC",
                     "SCGB3A2",
                     "SCGB1A1",
                     "PLIN2",
                     "COL1A1",
                     "KRT17",
                     "KRT5"),
        group.by = "celltype",
        split.by = "Authors",
        cols = c("red",
                 "green",
                 "blue"))

# ======================================
# Markers for each celltype
# ======================================
transitional_at2_markers <- FindMarkers(small,
                                        ident.1 = "Transitional AT2",
                                        ident.2 = c("SCGB3A2+",
                                                    "PLIN2+ Fibroblasts",
                                                    "KRT5-/KRT17+"),
                                        min.pct = 0.25)

transitional_at2_markers
transitional_at2_markers <- transitional_at2_markers[order(transitional_at2_markers$avg_logFC), ]
transitional_at2_markers <- rownames(transitional_at2_markers[1:50, ])

SCGB3A2_markers <- FindMarkers(small,
                               ident.1 = "SCGB3A2+",
                               ident.2 = c("Transitional AT2",
                                           "PLIN2+ Fibroblasts",
                                           "KRT5-/KRT17+"),
                               min.pct = 0.25)

SCGB3A2_markers
SCGB3A2_markers <- SCGB3A2_markers[order(SCGB3A2_markers$avg_logFC), ]
SCGB3A2_markers <- rownames(SCGB3A2_markers[1:50, ])

PLIN2_markers <- FindMarkers(small,
                             ident.1 = "PLIN2+ Fibroblasts",
                             ident.2 = c("SCGB3A2+",
                                         "Transitional AT2",
                                         "KRT5-/KRT17+"),
                             min.pct = 0.25)

PLIN2_markers
PLIN2_markers <- PLIN2_markers[order(PLIN2_markers$avg_logFC), ]
PLIN2_markers <- rownames(PLIN2_markers[1:50, ])

KRT5_markers <- FindMarkers(small,
                            ident.1 = "KRT5-/KRT17+",
                            ident.2 = c("SCGB3A2+",
                                        "PLIN2+ Fibroblasts",
                                        "Transitional AT2"),
                            min.pct = 0.25)

KRT5_markers
KRT5_markers <- KRT5_markers[order(KRT5_markers$avg_logFC), ]
KRT5_markers <- rownames(KRT5_markers[1:50, ])

gene_list <- c(transitional_at2_markers, SCGB3A2_markers, PLIN2_markers, KRT5_markers)

# Manually reorder gene list after RP & MT, removal
gene_list<-c("KRT5", "KRT7","KRT17", "KRT19","HOPX","MMP7","SLC34A2","CEACAM6","NAPSA","PLIN2", "THBS1","CCL2","VIM","COL1A1","COL6A2",
             "COL1A2","COL3A1", "IGFBP4","IGFBP6","IGFBP7","DCN","CFD","MGP","FBLN1","LUM","SERPINF1","CCDC80","TIMP1","SPARCL1",
             "PTGDS","PLA2G2A","SFRP2","C7","CTSL","CLU","C3", "AGER","ADIRF","SCGB3A2","SCGB3A1","SCGB1A1","SFTPB","SLPI",
             "WFDC2","SFTPC","SFTPA2","SFTPA1","CD74","CYB5A","CXCL17","SFTA2","BPIFB1",
             "FABP5","PIGR")

genes[!genes %in% gene_list]

keep <- rownames(small@meta.data[small@meta.data$Study %in% "GSE135893", ])
GSE135893 <- subset(small, cells = keep)
GSE135893 <- FindVariableFeatures(GSE135893, verbose = T, nfeatures = 3000)
GSE135893 <- ScaleData(GSE135893, features = row.names(GSE135893@assays$SCT@data))
GSE135893 <- RunPCA(GSE135893)
Idents(GSE135893) <- "celltype"

keep <-  rownames(small@meta.data[small@meta.data$Study == "GSE122960", ])
GSE122960 <- subset(small, cells = keep)
GSE122960 <- FindVariableFeatures(GSE122960, verbose = T, nfeatures = 3000)
GSE122960 <- ScaleData(GSE122960, features = row.names(GSE122960@assays$SCT@data))
GSE122960 <- RunPCA(GSE122960)
Idents(GSE122960) <- "celltype"

keep <-  rownames(small@meta.data[small@meta.data$Study == "GSE128033", ])
GSE128033 <- subset(small, cells = keep)
GSE128033 <- FindVariableFeatures(GSE128033, verbose = T, nfeatures = 3000)
GSE128033 <- ScaleData(GSE128033, features = row.names(GSE128033@assays$SCT@data))
GSE128033 <- RunPCA(GSE128033)
Idents(GSE128033) <- "celltype"

# ======================================
# Figure S: 23
# ======================================
DoHeatmap(subset(GSE135893, downsample = 100), features = gene_list, group.by = "celltype") + ggtitle("Habermann et al.")
DoHeatmap(subset(GSE122960, downsample = 100), features = gene_list, group.by = "celltype") + ggtitle("Reyfman et al.")
DoHeatmap(subset(GSE128033, downsample = 100), features = gene_list, group.by = "celltype") + ggtitle("Morse et al.")

table(GSE135893@meta.data$celltype)
table(GSE122960@meta.data$celltype)
table(GSE128033@meta.data$celltype)


