# ==============================================================================
# Author: Linh T. Bui (lbui@tgen.org)
# Date: 2023/03
# Purpose: running infercnv analysis for epithelial cells
# ==============================================================================
# ==============================================================================
# Environment variables
# ==============================================================================
set.seed(12345)

getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# Load libraries
library(infercnv)
library(Seurat)

# ==============================================================================
# Prepare files for infercnv on the Epithelial object 
# ==============================================================================
# Read in the current annotated file
epi <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_Epithelial_CT2.rds")

# Create a cellAnnotations file
# I used the Pathology for this, added the orig.ident on the seurat cluster numbers
cellAnnotations <- as.data.frame(cbind(rownames(epi@meta.data),
                                       epi@meta.data$Pathology))

colnames(cellAnnotations) <- NULL
cellAnnotations[,1] <- gsub("-", "_", cellAnnotations[,1])
rownames(cellAnnotations) <- cellAnnotations[,1]
cellAnnotations[,1] <- NULL
write.table(cellAnnotations, file="cellAnnotations.txt",sep= "\t",
            quote = F)

# Extract out the count matrix
singleCell.counts.matrix <- GetAssayData(epi, slot="counts", assay = "RNA")

# Change hyphen to underscore (infercnv will give erorrs if there's hyphen in cellbarcode)
colnames(singleCell.counts.matrix) <- gsub("-", "_", colnames(singleCell.counts.matrix))

write.table(round(singleCell.counts.matrix, digits=3), 
            file='singleCell.counts.matrix', quote=F, sep="\t")

# ==============================================================================
# Run Infercnv
# ==============================================================================
# Create the infercnv object
# Note: need to save all the related files in the working directory
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix="singleCell.counts.matrix",
                                    annotations_file="cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="/home/lbui/GRCh38-2020-A/GRCh38_ensemble98_gene_ordering_file.txt", 
                                    #this file has geneID, chr, start, end
                                    ref_group_names=c("Normal Appendix"),
                                    min_max_counts_per_cell=c(100,Inf))

# Perform infercnv operations to reveal CNV signal
infercnv_obj <- infercnv::run(infercnv_obj,
                             min_cells_per_gene = 10, #default is 3
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="HMMi3_output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=FALSE,   # group cells by CNV profile, not by sample
                             denoise=T,
                             num_threads=8,
                             analysis_mode='subclusters', # to identify tumor clusters
                             no_plot=TRUE,
                             HMM=T,
                             HMM_type="i3", # 3 states (neutral, depletion, duplication)
                             hclust_method='ward.D2',
                             sd_amplifier=3,  # sets midpoint for logistic
                             noise_logistic=TRUE, # turns gradient filtering on
                             num_ref_groups=2, #1st run plot looks like there are 2 groups in normal cells
                             tumor_subcluster_partition_method='random_trees',
                             k_obs_groups = 5) # split observation groups into 5 to find malignant + normal cells

saveRDS(infercnv_obj, file = "COH_epi_infercnv_subcluster_HMMi3_splitnormalcancerous.rds")

# Make a plot
plot_cnv(infercnv_obj,
         out_dir=".",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         cluster_by_groups=FALSE,
         x.center=1,
         x.range="auto",
         hclust_method='ward.D2',
        # color_safe_pal=TRUE, #using a color blindness safe palette
         output_filename="COH_epi_infercnv_subcluster_HMMi3_k5",
         output_format="pdf",
         k_obs_groups = 5,
         dynamic_resize=0)

# Visualize CNV values in UMAP space
## Adjust cell barcode in seurat object to match with inferncv (inferncv uses "_", not "-")
rownames(epi@meta.data) <- gsub("-", "_", rownames(epi@meta.data))
colnames(epi@assays$RNA@counts) <- rownames(epi@meta.data)
colnames(epi@assays$RNA@data) <- rownames(epi@meta.data)

## Add HMM info into Seurat object 
epi <- infercnv::add_to_seurat(infercnv_output_path="/scratch/lbui/RStudio_folder/20230727/HMMi3_output_dir/",
                              seurat_obj=epi, 
                              top_n=10)

## Adjust cell barcodes to match with other assays in the Seurat object
substring(rownames(epi@meta.data), 24, 24) <- "-"
colnames(epi@assays$RNA@counts) <- rownames(epi@meta.data)
colnames(epi@assays$RNA@data) <- rownames(epi@meta.data)
saveRDS(epi, file = "20230727_Epi_infercnv_HMMi3.rds")

## Make Featureplot
cnv_plot <- c("proportion_dupli_chr11","proportion_loss_chr11", "proportion_dupli_chr21",
              "proportion_loss_chr4","proportion_dupli_chr9","proportion_dupli_chr12",
              "proportion_dupli_chr19","proportion_dupli_chr20", "proportion_dupli_chr7", 
              "proportion_loss_chr8", "proportion_dupli_chr1","proportion_dupli_chr13")
FeaturePlot(epi, features = cnv_plot, ncol = 4) + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

## Featureplot with scaled data
cnv_plot <- c("proportion_scaled_dupli_chr1", "proportion_scaled_dupli_chr7", 
              "proportion_scaled_dupli_chr9", "proportion_scaled_dupli_chr11", 
              "proportion_scaled_dupli_chr19","proportion_scaled_dupli_chr20",
              "proportion_scaled_dupli_chr21", "proportion_scaled_loss_chr4")
FeaturePlot(epi, features = cnv_plot, ncol = 4) + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

