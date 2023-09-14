# ==============================================================================
# Author(s) : Linh T. Bui, lbui@tgen.org
# Date: 2023/08/02
# Description: Codes used to generate plots for Figure 6
# ==============================================================================

# ==============================================================================
# SET UP THE ENVIRONMENT VARIABLES 
# ==============================================================================
# Create a working directory to save files
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
#devtools::install_github("saeyslab/multinichenetr")
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)
library(RColorBrewer)

# ==============================================================================
# Read in the COH Appendix Seurat object
# ==============================================================================
coh.combined.sct <- readRDS("/scratch/lbui/Appendiceal_data/Appendiceal_integratedrpca_alllineages_final.rds")

# ====================================
# Prepare the data for MultiNichenet analysis
# ====================================
# subset out the 4 pathology groups
subset_ob <- subset(coh.combined.sct,
                    subset = Pathology2 %in% c("LAMN","LGMA","MHNA","GCA"))
subset_ob <- PrepSCTFindMarkers(subset_ob)

# Adjust CT2 label to group cell types with very few cells together
tumor <- c("Goblet-like cells","MUC5Bhi cells","SPINK4hi cells")
endo <- c("Venous Endothelial cell","Arterial Endothelial cells")
DCs <- c("pDCs","cDC1","cDC2")
bcells <- c("Follicular B","GALTB","GCBcell","Plasma B")
onion <- as.character(subset_ob$Celltype3)
onion <- ifelse(onion %in% epi, "Epithelial", 
                ifelse(onion %in% endo, "Endothelial", 
                       ifelse(onion %in% DCs, "Dendritic cells", 
                              ifelse(onion %in% bcells, "B cells", onion))))
subset_ob$Celltype3_new <- onion

# Convert to singlecellexperiment object
sce <- as.SingleCellExperiment(subset_ob, assay = "SCT") #run this for the pathology comparison

# Load the Nichenet network and matrices
# I had trouble reading the file directly so I downloaded to isilon and read from my scratch
#organism = "human"
#if(organism == "human"){
#  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
#    distinct(ligand, receptor) %>% 
#    mutate(ligand = make.names(ligand), receptor = make.names(receptor))
#  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
#  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
#  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
#} else if(organism == "mouse"){
#  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
#    distinct(ligand, receptor) %>% 
#    mutate(ligand = make.names(ligand), receptor = make.names(receptor))
#  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
#  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
#}

lr_network <- readRDS("/scratch/lbui/Appendiceal_data/lr_network_human_21122021.rds")
lr_network <- lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor) %>% 
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))
ligand_target_matrix <- readRDS("/scratch/lbui/Appendiceal_data/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>% make.names()

# Update gene symbols in scRNAseq to match with Nichenet network data
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

# ====================================
# Extract cell abundance and expression information
# ====================================
# Define group, sample ID and Celltypes
sample_id = "orig.ident"
group_id = "Pathology2" #compare among pathology groups
celltype_id = "Celltype3_new"
covariates = "Flowcell"
batches = NA

# Define sender and receiver (use ALl vs. All in this comparison)
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

# Extract cell type abundance and expression information
min_cells = 5 #min cells per group per sample (recommended in the vignette is 10)
abundance_expression_info <- get_abundance_expression_info(sce = sce, 
                                                           sample_id = sample_id, 
                                                           group_id = group_id, 
                                                           celltype_id = celltype_id, 
                                                           min_cells = min_cells, 
                                                           senders_oi = senders_oi, 
                                                           receivers_oi = receivers_oi, 
                                                           lr_network = lr_network, 
                                                           batches = batches)

# Plot cell type abundance result
abundance_expression_info$abund_plot_sample

# ====================================
# Perform genome-wide differential expression analysis of receiver and sender cell types
# ====================================
# Define the contrasts and covariates of interest for the DE analysis
contrasts_oi = c("'LAMN-(LGMA+MHNA+GCA)/3','LGMA-(LAMN+MHNA+GCA)/3','MHNA-(LGMA+LAMN+GCA)/3','GCA-(LGMA+LAMN+MHNA)/3'") #make sure there's no space 
contrast_tbl = tibble(contrast = 
                        c("LAMN-(LGMA+MHNA+GCA)/3", "LGMA-(LAMN+MHNA+GCA)/3",
                          "MHNA-(LGMA+LAMN+GCA)/3", "GCA-(LGMA+LAMN+MHNA)/3"), 
                      group = c("LAMN","LGMA","MHNA","GCA"))

# Define parameters for the analysis
logFC_threshold = 0.25
p_val_threshold = 0.05
fraction_cutoff = 0.05
p_val_adj = FALSE
empirical_pval = FALSE
top_n_target = 250
cores_system = 8
n.cores = min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type

# Define the weights of the prioritization of both expression, differential expression and NicheNet activity information
# I used Multinichenet default parameters
prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)
prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

# All group, sample, cell type, batch and covariate names should be syntactically valid 
SummarizedExperiment::colData(sce)$Celltype3_new <- SummarizedExperiment::colData(sce)$Celltype3_new %>% 
  make.names()
SummarizedExperiment::colData(sce)$Pathology2 <- SummarizedExperiment::colData(sce)$Pathology2 %>% make.names()
SummarizedExperiment::colData(sce)$Status <- SummarizedExperiment::colData(sce)$Status %>% make.names()
SummarizedExperiment::colData(sce)$Flowcell <- SummarizedExperiment::colData(sce)$Flowcell %>% make.names()

# Perform multinichenet analysis
multinichenet_output <- multi_nichenet_analysis(sce = sce, 
                                                celltype_id = celltype_id, 
                                                sample_id = sample_id, 
                                                group_id = group_id,
                                                lr_network = lr_network, 
                                                ligand_target_matrix = ligand_target_matrix, 
                                                contrasts_oi = contrasts_oi, 
                                                contrast_tbl = contrast_tbl, 
                                                batches = batches, 
                                                covariates = covariates,
                                                prioritizing_weights = prioritizing_weights, 
                                                min_cells = 5, 
                                                logFC_threshold = logFC_threshold, 
                                                p_val_threshold = p_val_threshold,  
                                                fraction_cutoff = fraction_cutoff, 
                                                p_val_adj = p_val_adj, 
                                                empirical_pval = empirical_pval, 
                                                top_n_target = top_n_target, 
                                                n.cores = n.cores, 
                                                sender_receiver_separate = FALSE, 
                                                verbose = TRUE)

# Save the output file and the prioritization table
saveRDS(multinichenet_output,
        file = "Appendiceal_Multinichenet_pathology_SCT_mincell5_allpathologies.rds")
write.csv(multinichenet_output$prioritization_tables$group_prioritization_tbl,
          file = "Appendiceal_Multinichenet_SCT_group_prioritization_tbl.csv")

# Check the output
multinichenet_output$celltype_info$avg_df %>% head()
multinichenet_output$celltype_info$frq_df %>% head()
multinichenet_output$celltype_info$avg_df_group %>% head()
multinichenet_output$celltype_info$frq_df_group %>% head()
multinichenet_output$celltype_info$rel_abundance_df %>% head()
multinichenet_output$celltype_de %>% head()
multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% head()
multinichenet_output$prioritization_tables$group_prioritization_tbl %>% head()

# Circos plot of top-prioritized links
prioritized_tbl_oi_top20 <- get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 
                                             20, rank_per_group = TRUE) #select top 20 per group
prioritized_tbl_oi <- multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_top20$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_top20)

prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers <- union(prioritized_tbl_oi$sender %>% 
                             unique(), prioritized_tbl_oi$receiver %>% 
                             unique()) %>% 
  sort()

# Set up colors (same color scheme as dimplot)
ct3_cols <- c("#d2331f","#ef862f","#ff8f8a","#ce9c48","#0072ae","#00c085",
              "#12a037","#bbc28b","#00bfc3","#57b3de","#ff64b7","#f4b990")
colors_sender = ct3_cols %>%
  magrittr::set_names(senders_receivers)
colors_receiver = ct3_cols %>% 
  magrittr::set_names(senders_receivers)

# Make circos plots
pdf("Multinichenet_LAMN_LGMA_MHNA_top20_circos_SCTassay_allpathologies.pdf")
circos_list <- make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
dev.off()

# Visualization of scaled ligand-receptor pseudobulk products and ligand activity
plot_oi <- make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, 
                                             prioritized_tbl_oi_top20,
                                             widths = c(8,3,3))
plot_oi

# Plot ligand activities
ligands_oi <- multinichenet_output$prioritization_tables$ligand_activities_target_de_tbl %>% 
  inner_join(contrast_tbl) %>% 
  group_by(group, receiver) %>% 
  distinct(ligand, receiver, group, activity) %>% 
  top_n(5, activity) %>% 
  pull(ligand) %>% 
  unique()

plot_oi <- make_ligand_activity_plots(multinichenet_output$prioritization_tables, 
                                      ligands_oi, contrast_tbl, widths = NULL)
plot_oi

