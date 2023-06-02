#==============================================================================#
# Author(s) : Stephanie L. Yahn,
#             Heini M Natri, hnatri@tgen.org
# Date: 2021/12/06
# Description: Cell type annotation of CITEseq and non-CITEseq products and 
# leukPBMCs
#==============================================================================#

#==============================================================================#
# Load libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/utilities.R")

#==============================================================================#
# CITEseq and nonCITEseq products
#==============================================================================#

products = readRDS("/scratch/hnatri/CART/CITEproducts_integrated_wnn_211212.rds")

DimPlot(products, reduction = "wnn.umap", label = T, group.by = "wsnn_res.1.2")

# Adjusting dimensions
DefaultAssay(products) <- "integrated_protein"
products = RunUMAP(products,
                   dims = 1:10,
                   reduction = "integrated_protein_pca",
                   reduction.key = "integrated_protein_umap_",
                   reduction.name = "integrated_protein_umap",
                   return.model = TRUE)
products = FindNeighbors(products,
                         dims = 1:10,
                         reduction = "integrated_protein_pca",
                         graph.name = c("integrated_protein_nn",
                                        "integrated_protein_snn"))
# resolution 0.5
products = FindClusters(products,
                        resolution = c(0.1,0.3,0.5,0.8,1),
                        graph.name = "integrated_protein_snn")

DimPlot(products, reduction = "integrated_protein_umap", label = T, group.by = "integrated_protein_snn_res.0.5")
FeaturePlot(products, reduction = "integrated_protein_umap", features = c("CD4", "CD8"))

Idents(products) <- "wsnn_res.1.2"
products@meta.data$seurat_clusters_2 <- products@meta.data$seurat_clusters
products@meta.data$seurat_clusters <- products@meta.data$wsnn_res.1.2

DefaultAssay(products) <- "integrated_protein"
FeaturePlot(products, features = c("CD14", "CD16", "CD3", "CD19", "CD11b", "CD11c"), 
            min.cutoff = "q01", max.cutoff = "q99", reduction="wnn.umap")
FeaturePlot(products, features = c("CD8", "CD4", "CD19"), reduction="wnn.umap")

DefaultAssay(products) <- "integrated_sct"
FeaturePlot(products, features = c("CD8A", "CD4"), reduction="wnn.umap")
FeaturePlot(products, features = c("percent.mt", "nCount_RNA", "nFeature_RNA", "nCount_Protein"), reduction="wnn.umap")
VlnPlot(products, features = c("percent.mt", "nCount_RNA", "nFeature_RNA", "nCount_Protein"), group.by = "wsnn_res.1.2", pt.size = 0, ncol = 2)

# Control antibodies
DefaultAssay(products) <- "Protein"
FeaturePlot(products, features = c("IgG2a-K-isotype-Ctrl", "IgG1-K-isotype-Ctrl"), reduction = "wnn.umap")

#==============================================================================#
# Transfer CITEseq product annotations to non-CITEseq products
#==============================================================================#

# Integrated non-CITEseq products
nonCITE.products <- readRDS(file = "/scratch/hnatri/CART/nonCITEproducts_integrated_211212.rds")

DefaultAssay(nonCITE.products) <- "SCT"
FeaturePlot(nonCITE.products, features = c("CD4", "CD8A"), reduction = "integrated_sct_umap")
FeaturePlot(nonCITE.products, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), reduction = "integrated_sct_umap")
VlnPlot(nonCITE.products, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), group.by = "integrated_sct_snn_res.0.8", pt.size = 0, ncol = 2)

DefaultAssay(products) <- "integrated_sct"
DefaultAssay(nonCITE.products) <- "integrated_sct"

# integrated_protein_snn_res.0.7
anchors = FindTransferAnchors(
  reference = products,
  query = nonCITE.products,
  normalization.method = "SCT",
  reference.reduction = "integrated_sct_pca",
  reference.assay = "integrated_sct",
  query.assay = "integrated_sct",
  dims = 1:30,
  recompute.residuals = FALSE # If using SCT as a normalization method, compute query Pearson residuals using the reference SCT model parameters.
)

nonCITE.products = MapQuery(
  anchorset = anchors,
  query = nonCITE.products,
  reference = products,
  refdata = list(
    integrated_protein_snn_res.0.5 = "integrated_protein_snn_res.0.5",
    integrated_protein_snn_res.0.8 = "integrated_protein_snn_res.0.8",
    wsnn_res.0.5 = "wsnn_res.0.5",
    wsnn_res.0.8 = "wsnn_res.0.8",
    wsnn_res.1 = "wsnn_res.1",
    predicted_ADT = "Protein",
    predicted_integrated_protein = "integrated_protein"
  ),
  reference.reduction = "integrated_sct_pca", 
  reduction.model = "wnn.umap" # TODO: integrated_sct_umap, integrated_protein_umap
)

DimPlot(nonCITE.products, group.by= "Batch", reduction = "ref.umap")
DimPlot(products, group.by= "Batch", reduction = "wnn.umap")

p1 <- DimPlot(products, reduction = "wnn.umap",
              group.by = "integrated_protein_snn_res.0.5",
              label = TRUE, label.size = 4, repel = TRUE)

p2 <- DimPlot(nonCITE.products, reduction = "ref.umap",
              group.by = "predicted.integrated_protein_snn_res.0.5",
              label = TRUE, label.size = 4, repel = TRUE)

p1 + p2

p3 <- DimPlot(products, reduction = "wnn.umap",
              group.by = "wsnn_res.1",
              label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

p4 <- DimPlot(nonCITE.products, reduction = "ref.umap",
              group.by = "predicted.wsnn_res.1",
              label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

p3 + p4

DefaultAssay(nonCITE.products) <- "integrated_sct"
FeaturePlot(nonCITE.products, features = c("CD4", "CD8A"), reduction = "ref.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

DefaultAssay(nonCITE.products) <- "SCT"
FeaturePlot(nonCITE.products, features = c("CD4", "CD8A"), reduction = "ref.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

DefaultAssay(products) <- "SCT"
FeaturePlot(products, features = c("CD4", "CD8A"), reduction = "wnn.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

DefaultAssay(nonCITE.products) <- "predicted_ADT"
FeaturePlot(nonCITE.products, features = c("CD4", "CD8"), reduction = "ref.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")
FeaturePlot(nonCITE.products, features = c("predictedADT_CD4", "predictedADT_CD8"), reduction = "ref.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")


DefaultAssay(nonCITE.products) <- "predicted_integrated_protein"
FeaturePlot(nonCITE.products, features = c("CD4", "CD8"), reduction = "ref.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")
FeaturePlot(nonCITE.products, features = c("predictedADT_CD4", "predictedADT_CD8"), reduction = "ref.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

DefaultAssay(products) <- "integrated_sct"
FeaturePlot(products, features = c("CD4", "CD8A"), reduction = "wnn.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

DefaultAssay(products) <- "SCT"
FeaturePlot(products, features = c("CD4", "CD8A"), reduction = "wnn.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

DefaultAssay(products) <- "Protein"
FeaturePlot(products, features = c("CD4", "CD8"), reduction = "wnn.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

# Saving the objects
saveRDS(nonCITE.products, "/scratch/hnatri/CART/noncite_products_211213.rds")
saveRDS(products, "/scratch/hnatri/CART/cite_products_211213.rds")

# Merging the CITEseq and non-CITEseq objects
# Renaming the assays
nonCITE.products@assays
nonCITE.products <- RenameAssays(nonCITE.products, "predicted_integrated_protein" = "integrated_protein")
nonCITE.products <- RenameAssays(nonCITE.products, "predicted_ADT" = "Protein")

names(nonCITE.products@reductions) <- gsub("ref.umap", "wnn.umap", names(nonCITE.products@reductions))

DefaultAssay(nonCITE.products) <- "SCT"
FeaturePlot(nonCITE.products, features = c("CD4", "CD8A"), reduction = "wnn.umap", 
            min.cutoff = "q01", ncol = 1, keep.scale = "all")

cite_noncite_product <- merge(x=nonCITE.products, y=products, merge.dr = c("wnn.umap"))

# Adding the AB panel to the metadata
cite_noncite_product@meta.data$ab_panel <- ifelse(cite_noncite_product@meta.data$Batch %in% c("Batch1"), "TS-C_99328",
                                                  ifelse(cite_noncite_product@meta.data$Batch %in% c("Batch2", "Batch3", "Batch4", "Batch5", "Batch6"), "TS-C_99814",
                                                         ifelse(cite_noncite_product@meta.data$Batch %in% c("Batch14", "Batch17", "Batch18"), "TS-C_399905",
                                                                ifelse(cite_noncite_product@meta.data$Batch %in% c("Batch20"), "TS-A_399907", "nonCITEseq"))))

DimPlot(cite_noncite_product, group.by = "wsnn_res.0.8", reduction = "wnn.umap", label = T)
DimPlot(cite_noncite_product, group.by = "wsnn_res.1", reduction = "wnn.umap", label = T)
DimPlot(cite_noncite_product, group.by = "Batch", reduction = "wnn.umap")
DimPlot(cite_noncite_product, group.by = "ab_panel", reduction = "wnn.umap")
DimPlot(cite_noncite_product, group.by = "Manufacture", reduction = "wnn.umap")

# Adding a new column with cluster IDs for CITEseq and nonCITEseq batches
cite_noncite_product@meta.data$cluster <- cite_noncite_product@meta.data$wsnn_res.1
cite_noncite_product@meta.data$cluster[is.na(cite_noncite_product@meta.data$wsnn_res.1)] <- cite_noncite_product@meta.data$predicted.wsnn_res.1[is.na(cite_noncite_product@meta.data$cluster)]

dim(cite_noncite_product@meta.data)
table(cite_noncite_product@meta.data$cluster)
DimPlot(cite_noncite_product, group.by = "cluster", reduction = "wnn.umap", label = T) + NoLegend()

# Setting and reordering Idents
Idents(cite_noncite_product) <- "cluster"
Idents(cite_noncite_product) <- factor(x = Idents(cite_noncite_product), levels = sort(as.numeric(levels(cite_noncite_product))))
DimPlot(cite_noncite_product, group.by = "cluster", split.by = "ab_panel", reduction = "wnn.umap")

#Add outcome to metadata
outcomes = read.csv(file = "/labs/banovich/BCTCSF/Stephanie/20200720_Combined_data.csv")
cite_noncite_product$response = plyr::mapvalues(cite_noncite_product$UPN, from = outcomes$UPN, to = outcomes$Best.Response.for..Brain.Disease...RANO.Criteria)
cite_noncite_product$response[cite_noncite_product$response %in% c("208", "NA", "Partial Response (PR)")] <- "NA"
cite_noncite_product$response[is.na(cite_noncite_product$response)] <- "NA"

DimPlot(cite_noncite_product, group.by = "response", reduction = "wnn.umap")
DimPlot(cite_noncite_product, group.by = "cluster", reduction = "wnn.umap", label = T) + NoLegend()
DimPlot(cite_noncite_product, group.by = "Batch", reduction = "wnn.umap")
DimPlot(cite_noncite_product, group.by = "ab_panel", reduction = "wnn.umap")
DimPlot(cite_noncite_product, group.by = "Manufacture", reduction = "wnn.umap")
DimPlot(cite_noncite_product, group.by = "Phase", reduction = "wnn.umap")

DefaultAssay(cite_noncite_product) <- "Protein"
FeaturePlot(cite_noncite_product, features = c("CD4", "CD8"), reduction = "wnn.umap")
VlnPlot(cite_noncite_product, features = c("CD4", "CD8"), pt.size=0)

DefaultAssay(cite_noncite_product) <- "SCT"
FeaturePlot(cite_noncite_product, features = c("CD4", "CD8A"), reduction = "wnn.umap")
VlnPlot(cite_noncite_product, features = c("CD4", "CD8A"), pt.size=0)

saveRDS(cite_noncite_product, "/scratch/hnatri/CART/cite_noncite_product_211213.rds")

# Adding the SCT assay from the CITEseq+nonCITEseq object to the object
# with the imputed protein levels for nonCITEseq

# CITEseq and nonCITEseq products with WNN
product <- readRDS("/scratch/hnatri/CART/cite_noncite_product_211213.rds")
# SCT integration across all batches was done separately
product_sct <- readRDS("/scratch/hnatri/CART/ALLbatches_product_integrated_rPCA_SCT_211211.rds")

DefaultAssay(product) <- "integrated_sct"

# Keeping cells present in both objects
product_sct <- subset(product_sct,
                      cells = colnames(product))

# Reclustering after subsetting
product_sct <- ScaleData(product_sct,
                         vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                         verbose = F)
product_sct <- RunPCA(product_sct,
                      reduction.name = "integrated_sct_pca",
                      verbose = F)
get_pcs(product_sct, reduction_name = "integrated_sct_pca")
product_sct = RunUMAP(product_sct,
                      reduction = "integrated_sct_pca",
                      reduction.name = "integrated_sct_umap",
                      dims = 1:10,
                      return.model = TRUE)
product_sct = FindNeighbors(product_sct,
                            reduction = "integrated_sct_pca",
                            graph.name = c("integrated_sct_nn", "integrated_sct_snn"),
                            dims = 1:10)
# resolution 0.2
product_sct = FindClusters(product_sct,
                           resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                           graph.name = "integrated_sct_snn")

DimPlot(product_sct, group.by = "integrated_sct_snn_res.0.5", reduction = "integrated_sct_umap")
DimPlot(product_sct, group.by = "Batch", reduction = "integrated_sct_umap")

# Assigning the new assay and dimensionality reductions
product@assays$integrated_sct <- product_sct@assays$integrated_sct
product@reductions$integrated_sct_pca <- product_sct@reductions$integrated_sct_pca
product@reductions$integrated_sct_umap <- product_sct@reductions$integrated_sct_umap

# Adding a new column with cluster IDs for CITEseq and nonCITEseq batches
product@meta.data$cluster <- product@meta.data$wsnn_res.1
product@meta.data$cluster[is.na(product@meta.data$wsnn_res.1)] <- product@meta.data$predicted.wsnn_res.1[is.na(product@meta.data$cluster)]

DefaultAssay(product) <- "integrated_protein"
product$cluster <- as.numeric(as.character(product$cluster))
DimPlot(product, reduction = "wnn.umap", group.by = "cluster", label = T)

# Rescaling the protein data with ab_panel (including nonCITEseq) as a covariate
product <- ScaleData(product, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score", "ab_panel"))

# Saving the merged object
saveRDS(product, "/scratch/hnatri/CART/cite_noncite_product_integrated_merged_211213.rds")

