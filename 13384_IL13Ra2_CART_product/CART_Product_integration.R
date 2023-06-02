#==============================================================================#
# Author(s) : Stephanie L. Yahn,
#             Heini M Natri, hnatri@tgen.org
# Date: 2021/11/30
# Description: Dataset integration for the CAR T project
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyr)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/Utilities/utilities.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Integrate RNA for all batches, Product
#==============================================================================#

# Processing was done according to CART_Product_leukPBMC_preprocess.R

# The batch_list_filtered_filtered object is a list containing all 19 batches
# The batches have been QC filtered, subset to contain only leuk PBMC and Product, 
# RNA has been normalized and scaled using SCTransform
# Protein has been normalized and scaled using NormalizeData(normalization.method = "CLR") and ScaleData()

batch_list_filtered = readRDS(file = "/labs/banovich/BCTCSF/Heini/batch_list_filtered_211207.rds")

# 213 (in the clinical paper) only has 15 cells after filtering
names(batch_list_filtered)
unique(batch_list_filtered[["Batch16_filtered"]]$UPN)
"213" %in% batch_list_filtered[["Batch16_filtered"]]$UPN

# Renaming cells
batch_list_filtered <- lapply(batch_list_filtered, function(x){
  renamed.assay <- RenameCells(x,
                               new.names = paste0(unique(x$Batch), "_", colnames(x)))
}
)

# Product only, all batches
products_list = list()
for (i in 1:length(batch_list_filtered)) {
  Idents(batch_list_filtered[[i]]) = batch_list_filtered[[i]]$Sample_Type
  products_list[[i]] = subset(batch_list_filtered[[i]], idents = "Product")
  # rerun SCTransform because we're working with a subset
  DefaultAssay(products_list[[i]]) = "RNA"
  products_list[[i]] = SCTransform(products_list[[i]],
                                   method = "glmGamPoi",
                                   vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                   verbose = F)
}

# Set nfeatures to 1000 to avoid error at IntegrateData
products.features = SelectIntegrationFeatures(object.list = products_list,
                                              nfeatures = 1000)

# Running PCA for rPCA integration
for (i in 1:length(batch_list_filtered)) {
  products_list[[i]] = RunPCA(products_list[[i]],
                              features = products.features,
                              verbose = F)
}

products_list = PrepSCTIntegration(object.list = products_list,
                                   anchor.features = products.features, 
                                   verbose = F)

# Identify anchors and integrate the datasets based on RNA
# Using k=20 neighbors to find anchors (default = 5)?
# Using rPCA to avoid "problem too large" error with CCA
products.anchors = FindIntegrationAnchors(object.list = products_list,
                                          normalization.method = "SCT", 
                                          anchor.features = products.features,
                                          #reference = c(1,17,18,3),
                                          k.anchor = 20,
                                          dims = 1:30,
                                          reduction = "rpca")
# Default value for k.weight causes an error here (not enough cells for some 
# samples?)
products.integrated = IntegrateData(anchorset = products.anchors,
                                    normalization.method = "SCT", 
                                    new.assay.name = "integrated_sct",
                                    #k.weight = 80,
                                    verbose = T)
unique(products.integrated$Batch)
# No need to run ScaleData if you've used SCT integration
products.integrated = RunPCA(products.integrated,
                             reduction.name = "integrated_sct_pca",
                             verbose = F)
get_pcs(products.integrated, reduction_name = "integrated_sct_pca")
products.integrated = RunUMAP(products.integrated,
                              reduction = "integrated_sct_pca",
                              reduction.name = "integrated_sct_umap",
                              dims = 1:11,
                              return.model = TRUE)
products.integrated = FindNeighbors(products.integrated,
                                    reduction = "integrated_sct_pca",
                                    dims = 1:11,
                                    graph.name = c("integrated_sct_nn",
                                                   "integrated_sct_snn"))
# resolution 0.2
products.integrated = FindClusters(products.integrated,
                                   resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                   graph.name = "integrated_sct_snn")

products.integrated@meta.data$ab_panel <- ifelse(products.integrated@meta.data$Batch %in% c("Batch1"), "TS-C_99328",
                                                 ifelse(products.integrated@meta.data$Batch %in% c("Batch2", "Batch3", "Batch4", "Batch5", "Batch6"), "TS-C_99814",
                                                        ifelse(products.integrated@meta.data$Batch %in% c("Batch14", "Batch17", "Batch18"), "TS-C_399905",
                                                               ifelse(products.integrated@meta.data$Batch %in% c("Batch20"), "TS-A_399907", "nonCITE"))))

DimPlot(products.integrated, group.by = "ab_panel", split.by = "integrated_sct_snn_res.0.2", ncol=4, reduction = "integrated_sct_umap")
DimPlot(products.integrated, group.by = "Batch", reduction = "integrated_sct_umap")
DimPlot(products.integrated, group.by = "Manufacture", reduction = "integrated_sct_umap")
DimPlot(products.integrated, group.by = "UPN", reduction = "integrated_sct_umap")
DimPlot(products.integrated, group.by = "integrated_sct_snn_res.0.2", reduction = "integrated_sct_umap")
FeaturePlot(products.integrated, features = "CD8A", reduction = "integrated_sct_umap")
FeaturePlot(products.integrated, features = c("nCount_RNA", "nFeature_RNA","nCount_Protein", "nFeature_Protein", "percent.mt"), reduction = "integrated_sct_umap")

saveRDS(products.integrated,  "/scratch/hnatri/CART/ALLbatches_product_integrated_rPCA_SCT_211211.rds")

#==============================================================================#
# Integrate CITEseq Products (SCT and Protein)
#==============================================================================#

batch_list_filtered = readRDS(file = "/scratch/hnatri/CART/batch_list_filtered_211207.rds")

CITE = c("Batch1_filtered", "Batch2_filtered", "Batch3_filtered", "Batch4_filtered", "Batch5_filtered", 
         "Batch6_filtered", "Batch14_filtered", "Batch17_filtered", "Batch18_filtered")

# Keeping only CITEseq Products
products_list = list()
for (i in CITE) {
  Idents(batch_list_filtered[[i]]) = batch_list_filtered[[i]]$Sample_Type
  products_list[[i]] = subset(batch_list_filtered[[i]], idents = "Product")
  # rerun SCTransform because we're working with a subset
  DefaultAssay(products_list[[i]]) = "RNA"
  products_list[[i]] = SCTransform(products_list[[i]],
                                   method = "glmGamPoi",
                                   vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                   verbose = T)
}

options(future.globals.maxSize = 30000 * 1024^2)

# Set nfeatures to 1000 to avoid error at IntegrateData
products.features = SelectIntegrationFeatures(object.list = products_list,
                                              nfeatures = 1000)

# Running PCA for rPCA integration
for (i in 1:length(products_list)) {
  products_list[[i]] = RunPCA(products_list[[i]],
                              features = products.features,
                              verbose = F)
}

products_list = PrepSCTIntegration(object.list = products_list,
                                   anchor.features = products.features, 
                                   verbose = F)


# Identify anchors and integrate the datasets based on RNA using CCA
products.anchors = FindIntegrationAnchors(object.list = products_list,
                                          normalization.method = "SCT", 
                                          anchor.features = products.features,
                                          reduction = "rpca",
                                          dims = 1:30,
                                          k.anchor = 20)
CITE.products.integrated = IntegrateData(anchorset = products.anchors,
                                         normalization.method = "SCT", 
                                         new.assay.name = "integrated_sct",
                                         verbose = F)

# No need to run ScaleData if you've used SCT integration pipeline
CITE.products.integrated = RunPCA(CITE.products.integrated,
                                  reduction.name = "integrated_sct_pca",
                                  verbose = F)
get_pcs(CITE.products.integrated, reduction_name = "integrated_sct_pca")
CITE.products.integrated = RunUMAP(CITE.products.integrated,
                                   reduction = "integrated_sct_pca",
                                   reduction.name = "integrated_sct_umap",
                                   dims = 1:11,
                                   return.model = TRUE)
CITE.products.integrated = FindNeighbors(CITE.products.integrated,
                                         reduction = "integrated_sct_pca",
                                         dims = 1:11,
                                         graph.name = c("integrated_sct_nn",
                                                        "integrated_sct_snn"))
# resolution 0.2
CITE.products.integrated = FindClusters(CITE.products.integrated,
                                        resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                        graph.name = "integrated_sct_snn")

DimPlot(CITE.products.integrated, group.by = "Batch", reduction = "integrated_sct_umap")
DimPlot(CITE.products.integrated, group.by = "Manufacture", reduction = "integrated_sct_umap")
DimPlot(CITE.products.integrated, group.by = "UPN", reduction = "integrated_sct_umap")
DimPlot(CITE.products.integrated, group.by = "integrated_sct_snn_res.0.2", reduction = "integrated_sct_umap")
FeaturePlot(CITE.products.integrated, features = "CD8A", reduction = "integrated_sct_umap")

# Integrate protein
for (i in CITE) {
  DefaultAssay(products_list[[i]]) = "Protein"
  # use all Protein features for dimensional reduction
  VariableFeatures(products_list[[i]]) = rownames(products_list[[i]][["Protein"]])
}

protein.features = SelectIntegrationFeatures(object.list = products_list)

for (i in CITE) {
  products_list[[i]] = ScaleData(products_list[[i]],
                                 assay = "Protein",
                                 features = protein.features,
                                 verbose = F,
                                 vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  products_list[[i]] = RunPCA(products_list[[i]],
                              assay = "Protein",
                              features = protein.features,
                              verbose = F)
}

protein.anchors = FindIntegrationAnchors(object.list = products_list,
                                         reduction = "rpca", 
                                         dims = 1:30,
                                         k.anchor = 20,
                                         verbose = F)
protein.integrated = IntegrateData(anchorset = protein.anchors,
                                   dims = 1:30,
                                   new.assay.name = "integrated_protein")

protein.integrated@meta.data$ab_panel <- ifelse(protein.integrated@meta.data$Batch %in% c("Batch1"), "TS-C_99328",
                                                ifelse(protein.integrated@meta.data$Batch %in% c("Batch2", "Batch3", "Batch4", "Batch5", "Batch6"), "TS-C_99814",
                                                       ifelse(protein.integrated@meta.data$Batch %in% c("Batch14", "Batch17", "Batch18"), "TS-C_399905",
                                                              ifelse(protein.integrated@meta.data$Batch %in% c("Batch20"), "TS-A_399907", "nonCITE"))))

unique(protein.integrated@meta.data$ab_panel)
DefaultAssay(protein.integrated) = "integrated_protein"
protein.integrated = ScaleData(protein.integrated,
                               assay ="integrated_protein",
                               verbose = F,
                               vars.to.regress = c("percent.mt", "S.Score", "G2M.Score","ab_panel"))
protein.integrated = RunPCA(protein.integrated,
                            assay ="integrated_protein",
                            reduction.name = "integrated_protein_pca",
                            verbose = F)
get_pcs(protein.integrated, reduction = "integrated_protein_pca")
protein.integrated = RunUMAP(protein.integrated,
                             dims = 1:11,
                             reduction = "integrated_protein_pca",
                             reduction.name = "integrated_protein_umap",
                             return.model = TRUE)
protein.integrated = FindNeighbors(protein.integrated,
                                   dims = 1:11,
                                   reduction = "integrated_protein_pca",
                                   graph.name = c("integrated_protein_nn",
                                                  "integrated_protein_snn"))
# resolution 0.5
protein.integrated = FindClusters(protein.integrated,
                                  resolution = c(0.1,0.3,0.5,0.8,1),
                                  graph.name = "integrated_protein_snn")
DimPlot(protein.integrated, group.by = "Batch", reduction = "integrated_protein_umap")
DimPlot(protein.integrated, group.by = "Manufacture", reduction = "integrated_protein_umap")
DimPlot(protein.integrated, group.by = "integrated_protein_snn_res.0.5", reduction = "integrated_protein_umap", label = T)
DimPlot(protein.integrated, group.by = "UPN", reduction = "integrated_protein_umap")
FeaturePlot(protein.integrated, features = "CD8", reduction = "integrated_protein_umap")

# Control antibodies
DefaultAssay(protein.integrated)<- "Protein"
FeaturePlot(protein.integrated, features = c("IgG2a-K-isotype-Ctrl", "IgG1-K-isotype-Ctrl"), reduction = "integrated_protein_umap",
            min.cutoff = "q01", max.cutoff = "q99")
VlnPlot(protein.integrated, features = c("IgG2a-K-isotype-Ctrl", "IgG1-K-isotype-Ctrl"),
        group.by = "integrated_protein_snn_res.0.5", pt.size=0)
FeaturePlot(protein.integrated, features = c("CD8"))

# [remove clusters that stains for most proteins including controls, then re-cluster]
# Removing clusters that stain for control antibodies
protein.integrated = subset(protein.integrated,
                            cells = row.names(protein.integrated@meta.data[!(protein.integrated@meta.data$integrated_protein_snn_res.0.5 %in% c("6", "10")), ]))
# Rescaling
protein.integrated <- ScaleData(protein.integrated,
                                vars.to.regress = c("percent.mt", "S.Score", "G2M.Score", "ab_panel"))
# Reclustering protein
DefaultAssay(protein.integrated) = "integrated_protein"

get_pcs(protein.integrated, reduction = "integrated_protein_pca")
protein.integrated = RunUMAP(protein.integrated,
                             dims = 1:11,
                             reduction = "integrated_protein_pca",
                             reduction.name = "integrated_protein_umap",
                             return.model = TRUE)
FeaturePlot(protein.integrated, features = c("CD4", "CD8"), reduction = "integrated_protein_umap")

protein.integrated = FindNeighbors(protein.integrated,
                                   dims = 1:11,
                                   reduction = "integrated_protein_pca",
                                   graph.name = c("integrated_protein_nn",
                                                  "integrated_protein_snn"))
# resolution 0.5
protein.integrated = FindClusters(protein.integrated,
                                  resolution = c(0.1,0.3,0.5,0.8,1),
                                  graph.name = "integrated_protein_snn")
DimPlot(protein.integrated, group.by = "Batch")
DimPlot(protein.integrated, group.by = "Manufacture")
DimPlot(protein.integrated, group.by = "UPN")
DimPlot(protein.integrated, group.by = "integrated_protein_snn_res.0.5")

# Control antibodies
DefaultAssay(protein.integrated)<- "Protein"
FeaturePlot(protein.integrated, features = c("IgG2a-K-isotype-Ctrl", "IgG1-K-isotype-Ctrl"), reduction = "integrated_protein_umap",
            min.cutoff = "q01", max.cutoff = "q99")
VlnPlot(protein.integrated, features = c("IgG2a-K-isotype-Ctrl", "IgG1-K-isotype-Ctrl"),
        group.by = "integrated_protein_snn_res.0.5", pt.size=0)
FeaturePlot(protein.integrated, features = c("CD8"))

# Removing the same cells from CITE.batch.integrated
DefaultAssay(CITE.products.integrated) <- "integrated_sct"
CITE.products.integrated <- subset(CITE.products.integrated,
                                   cells = colnames(protein.integrated))
CITE.products.integrated <- ScaleData(CITE.products.integrated,
                                      vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
# Reclustering
get_pcs(CITE.products.integrated, reduction_name = "integrated_sct_pca")
CITE.products.integrated = RunUMAP(CITE.products.integrated,
                                   reduction = "integrated_sct_pca",
                                   reduction.name = "integrated_sct_umap",
                                   dims = 1:11,
                                   return.model = TRUE)
CITE.products.integrated = FindNeighbors(CITE.products.integrated,
                                         reduction = "integrated_sct_pca",
                                         graph.name = c("integrated_sct_nn", "integrated_sct_snn"),
                                         dims = 1:11)
# resolution 0.2
CITE.products.integrated = FindClusters(CITE.products.integrated,
                                        resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                        graph.name = "integrated_sct_snn")

DimPlot(CITE.products.integrated, group.by = "integrated_sct_snn_res.0.5", reduction = "integrated_sct_umap")
DimPlot(CITE.products.integrated, group.by = "Batch", reduction = "integrated_sct_umap")
DimPlot(CITE.products.integrated, group.by = "Manufacture", reduction = "integrated_sct_umap")
DimPlot(CITE.products.integrated, group.by = "UPN", reduction = "integrated_sct_umap")
#length(unique(CITE.batch.integrated@meta.data$UPN))
FeaturePlot(CITE.products.integrated, features = c("CD8A"), reduction = "integrated_sct_umap")

# Adding the integrated protein assay and reductions to CITE.products.integrated
# object
CITE.products.integrated[["integrated_protein"]] = protein.integrated[["integrated_protein"]]
CITE.products.integrated[["integrated_protein_pca"]] = protein.integrated[["integrated_protein_pca"]]
CITE.products.integrated[["integrated_protein_umap"]] = protein.integrated[["integrated_protein_umap"]]

saveRDS(CITE.products.integrated, file = "/scratch/hnatri/CART/CITEproducts_integrated_211212.rds")

#q(save="no")

#==============================================================================#
# Weighted Nearest Neighbor (WNN) of CITEseq Products (SCT and Protein)
#==============================================================================#

# Already ran RunPCA() for the integrated_sct and integrated_protein assays, 
# so can move directly to FindMultiModalNeighbors()
# CITEseq batches only
# Product only
CITE.products.integrated = FindMultiModalNeighbors(CITE.products.integrated,
                                                   reduction.list = list("integrated_sct_pca", "integrated_protein_pca"),
                                                   dims.list = list(1:10, 1:20))

CITE.products.integrated = RunUMAP(CITE.products.integrated,
                                   nn.name = "weighted.nn",
                                   reduction.name = "wnn.umap",
                                   reduction.key = "wnnUMAP_",
                                   return.model = TRUE)
# resolution 0.3
CITE.products.integrated = FindClusters(CITE.products.integrated,
                                        graph.name = "wsnn",
                                        algorithm = 3,
                                        resolution = c(0.1,0.3,0.5,0.8,1,1.2),
                                        verbose = FALSE)

# Visualizations
DimPlot(CITE.products.integrated, group.by = "wsnn_res.1.2", label = T, repel = T, reduction = 'wnn.umap') + NoLegend()
DimPlot(CITE.products.integrated, group.by = "Manufacture", reduction = 'wnn.umap')

DefaultAssay(CITE.products.integrated) <- "Protein"
FeaturePlot(CITE.products.integrated, features = c("CD4", "CD8"), reduction = 'wnn.umap')
VlnPlot(CITE.products.integrated, features = c("CD4", "CD8"), group.by = "wsnn_res.1", pt.size = 0) + NoLegend()

table(CITE.products.integrated@meta.data$wsnn_res.1.2)

DefaultAssay(CITE.products.integrated) <- "RNA"
FeaturePlot(CITE.products.integrated, features = c("CD4", "CD8A"), reduction = 'wnn.umap')
VlnPlot(CITE.products.integrated, features = c("CD4", "CD8A"), group.by = "wsnn_res.0.5", pt.size = 0) + NoLegend()

saveRDS(CITE.products.integrated, file = "/scratch/hnatri/CART/CITEproducts_integrated_wnn_211212.rds")

p1 = DimPlot(CITE.products.integrated, reduction = 'wnn.umap', group.by = "Batch")
p2 = DimPlot(CITE.products.integrated, reduction = 'wnn.umap', group.by = 'Manufacture')
p1 + p2

p3 = DimPlot(CITE.products.integrated, reduction = 'integrated_sct_umap', 
             group.by = "integrated_sct_snn_res.0.2", label = F) + NoLegend()
p4 = DimPlot(CITE.products.integrated, reduction = 'integrated_protein_umap',
             group.by = "integrated_protein_snn_res.0.5", label = F) + NoLegend()
p5 = DimPlot(CITE.products.integrated, reduction = 'wnn.umap', 
             group.by = "wsnn_res.0.3", label = F) + NoLegend()

p3 + p4 + p5

DimPlot(CITE.products.integrated, reduction = 'wnn.umap', label = T, split.by = "Sample_Type") + NoLegend()
FeaturePlot(CITE.products.integrated, reduction = "wnn.umap", features = feature_markers,
            min.cutoff = "q02", max.cutoff = "q99", ncol = 3)

DimPlot(CITE.products.integrated, reduction = "wnn.umap", label = F, split.by = "UPN", ncol = 9) + NoLegend()

#==============================================================================#
# Integrate nonCITEseq Products (SCT)
#==============================================================================#

nonCITE = c("Batch7_filtered", "Batch8_filtered", "Batch9_filtered", "Batch10_filtered", "Batch11_filtered", 
            "Batch12_filtered", "Batch13_filtered", "Batch15_filtered", "Batch16_filtered", "Batch19_filtered")

# keep only nonCITEseq Products
products_list = list()
for (i in nonCITE) {
  Idents(batch_list_filtered[[i]]) = batch_list_filtered[[i]]$Sample_Type
  products_list[[i]] = subset(batch_list_filtered[[i]], idents = "Product")
  # rerun SCTransform because we're working with a subset
  DefaultAssay(products_list[[i]]) = "RNA"
  products_list[[i]] = SCTransform(products_list[[i]],
                                   method = "glmGamPoi",
                                   vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                                   verbose = F)
}

options(future.globals.maxSize = 30000 * 1024^2)

# Set nfeatures to 1000 to avoid error
products.features = SelectIntegrationFeatures(object.list = products_list,
                                              nfeatures = 1000)

# Running PCA for rPCA integration
for (i in 1:length(products_list)) {
  products_list[[i]] = RunPCA(products_list[[i]],
                              features = products.features,
                              verbose = F)
}

products_list = PrepSCTIntegration(object.list = products_list,
                                   anchor.features = products.features, 
                                   verbose = F)


# Identify anchors and integrate the datasets based on RNA using CCA
products.anchors = FindIntegrationAnchors(object.list = products_list,
                                          normalization.method = "SCT", 
                                          anchor.features = products.features,
                                          reduction = "rpca",
                                          dims = 1:30,
                                          k.anchor = 20)
nonCITE.products.integrated = IntegrateData(anchorset = products.anchors,
                                            normalization.method = "SCT", 
                                            new.assay.name = "integrated_sct",
                                            #k.weight = 80,
                                            verbose = T)

# No need to run ScaleData after SCT integration
nonCITE.products.integrated = RunPCA(nonCITE.products.integrated,
                                     reduction.name = "integrated_sct_pca",
                                     verbose = F)
get_pcs(nonCITE.products.integrated, reduction_name = "integrated_sct_pca")
nonCITE.products.integrated = RunUMAP(nonCITE.products.integrated,
                                      reduction = "integrated_sct_pca",
                                      reduction.name = "integrated_sct_umap",
                                      dims = 1:11,
                                      return.model = TRUE)
nonCITE.products.integrated = FindNeighbors(nonCITE.products.integrated,
                                            reduction = "integrated_sct_pca",
                                            dims = 1:11,
                                            graph.name = c("integrated_sct_nn",
                                                           "integrated_sct_snn"))
# resolution 0.2
nonCITE.products.integrated = FindClusters(nonCITE.products.integrated,
                                           resolution = c(0.1,0.2,0.3,0.5,0.8,1),
                                           graph.name = "integrated_sct_snn")

DimPlot(nonCITE.products.integrated, group.by = "Batch", reduction = "integrated_sct_umap")
DimPlot(nonCITE.products.integrated, group.by = "Manufacture", reduction = "integrated_sct_umap")
DimPlot(nonCITE.products.integrated, group.by = "UPN", reduction = "integrated_sct_umap")
DimPlot(nonCITE.products.integrated, group.by = "integrated_sct_snn_res.0.2", reduction = "integrated_sct_umap")
DefaultAssay(nonCITE.products.integrated) <- "SCT"
FeaturePlot(nonCITE.products.integrated, features = c("CD4", "CD8A"), reduction = "integrated_sct_umap")

saveRDS(nonCITE.products.integrated, file = "/scratch/hnatri/CART/nonCITEproducts_integrated_211212.rds")
