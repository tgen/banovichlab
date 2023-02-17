# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 12/15/2022
# Description: ILD cell type proportions by genotype
# ==============================================================================

# ======================================
# Import libraries
# ======================================

library(Seurat)
library(plyr)
library(dplyr)
library(qdap)
library(data.table)
library(tidyverse)
library(speckle)
library(vcfR)
library(googlesheets4)
library(ggplot2)

# ======================================
# Helper functions
# ======================================

source("/home/hnatri/Utilities/utilities.R")

# ======================================
# Environment variables
# ======================================

set.seed(1234)

# Colors
source("/home/hnatri/ILD_processing_annotation/lung_celltype_colors.R")
gs4_deauth()
paper_celltype_colors  <- gs4_get("https://docs.google.com/spreadsheets/d/1NnLlftpsqq_JEaINBMOe2T4qB1_XPihrcl1mNDim8B4/edit?usp=sharing")
sheet_names(paper_celltype_colors)
paper_celltype_colors <- read_sheet(paper_celltype_colors, sheet = "Sheet1")
head(paper_celltype_colors)

paper_celltype_colors[which(paper_celltype_colors[,c("Pretty cell type name")]=="cDC2"),]

lineage_cols <- paper_celltype_colors$lineage_color
names(lineage_cols) <- paper_celltype_colors$Lineage
celltype_cols <- paper_celltype_colors$celltype_color
names(celltype_cols) <- paper_celltype_colors$`Pretty cell type name`

# ======================================
# Import data
# ======================================

# Seurat objects
object_list <- readRDS("/labs/banovich/IPF/eQTL/object_list.rds")
object_names <- names(object_list)

# All GWAS risk variants
gwas_table <- gs4_get("https://docs.google.com/spreadsheets/d/1ZqMRkxeNFE6Ln10IQfdkXlgjG3jeTFo8L5QSlr9I9yU/edit?usp=sharing")
gwas_table <- read_sheet(gwas_table, sheet = "Sheet1")

gwas_table$rsid <- sapply(strsplit(gwas_table$`Variant and risk allele`, split='-', fixed=TRUE), `[`, 1)
gwas_table$risk_allele <- sapply(strsplit(gwas_table$`Variant and risk allele`, split='-', fixed=TRUE), `[`, 2)
gwas_table$chr_pos <- paste0(gwas_table$chr, "_", gwas_table$pos)

# mashr top SNPs
mashr_snp_info <- fread("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/filtered_MAF-HWE-INDPW_snp-info.txt", header=T)
mashr_sighits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-significant-eQTL.rds")
mashr_tophits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-top-eQTL.rds")

mashr_tophits$rowData$chr_pos <- paste0("chr", mashr_tophits$rowData$snp_chr, "_", mashr_tophits$rowData$snp_loc)

# Genotype data
vcf <- read.vcfR("/labs/banovich/IPF/eQTL/qtl_mapping_vcf/filtered_MAF-HWE-INDPW.vcf", verbose = FALSE)
gt_data <- extract.gt(vcf,
                      element = "GT",
                      mask = FALSE,
                      as.numeric = F,
                      return.alleles = FALSE,
                      IDtoRowNames = TRUE,
                      extract = TRUE,
                      convertNA = TRUE)
var_info <- vcf@fix
rm(vcf)

var_info <- as.data.frame(var_info)
gt_data <- as.data.frame(gt_data)

gt_data <- t(gt_data)

# Fixing sample names
rownames(gt_data) <- gsub("VU", "", rownames(gt_data))
rownames(gt_data) <- gsub("TI", "", rownames(gt_data))
rownames(gt_data) <- gsub("HD65", "HD065", rownames(gt_data))
rownames(gt_data) <- gsub("HD66", "HD066", rownames(gt_data))
rownames(gt_data) <- gsub("HD67", "HD067", rownames(gt_data))
rownames(gt_data) <- gsub("HD68", "HD068", rownames(gt_data))
rownames(gt_data) <- gsub("HD69", "HD069", rownames(gt_data))
rownames(gt_data) <- gsub("HD70", "HD070", rownames(gt_data))
rownames(gt_data) <- gsub("HD84", "HD084", rownames(gt_data))
rownames(gt_data) <- gsub("HD85", "HD085", rownames(gt_data))
rownames(gt_data) <- gsub("HD92", "HD092", rownames(gt_data))
rownames(gt_data) <- gsub("HD94", "HD094", rownames(gt_data))
rownames(gt_data) <- gsub("HD95", "HD095", rownames(gt_data))
rownames(gt_data) <- gsub("HD98", "HD098", rownames(gt_data))

object_list <- lapply(object_list, function(xx){
    xx$clean_Sample_Name <- gsub("VU", "", xx$clean_Sample_Name)
    xx$clean_Sample_Name <- gsub("TI", "", xx$clean_Sample_Name)
    
    subset(xx, subset = clean_Sample_Name %in% rownames(gt_data))
})
names(object_list) <- object_names

# ======================================
# Running propeller
# ======================================

# Running the analysis for immune cells and other cells
names(object_list)
non_immune <- merge(x = object_list[[1]], y = object_list[3:4])

# Running separately for ILD and CTRL
unique(non_immune$Status)
non_immune_ILD <- subset(non_immune, subset = Status == "Disease")
non_immune_CTRL <- subset(non_immune, subset = Status == "Control")
immune_ILD <- subset(object_list[["immune"]], subset = Status == "Disease")
immune_CTRL <- subset(object_list[["immune"]], subset = Status == "Control")

propeller_object_list <- list("non_immune" = non_immune,
                              "immune" = object_list[["immune"]],
                              "non_immune_ILD" = non_immune_ILD,
                              "non_immune_CTRL" = non_immune_CTRL,
                              "immune_ILD" = immune_ILD,
                              "immune_CTRL" = immune_CTRL)
propeller_object_names <- names(propeller_object_list)

# IPF GWAS meta-analysis nominal stats
ipf_gwas <- fread("/labs/banovich/IPF/GWAS/allen_gwas_harmonized.tsv.gz")
dim(ipf_gwas)
# Selecting nominally significant hits with a relaxed threshold
ipf_gwas_sig <- ipf_gwas[which(ipf_gwas$p_value < 1e-6),]
ipf_gwas_sig <- ipf_gwas_sig[which(ipf_gwas_sig$rsid %in% colnames(gt_data)),]

# Only keeping SNPs we want to test for
gt_data <- gt_data[,which(colnames(gt_data) %in% unique(mashr_sighits$rowData$variant_id))]

# Looping through all Seurat objects
propeller_res_list <- lapply(names(propeller_object_list), function(xx){
    # xx <- names(propeller_object_list)[1]
    transformed_props <- speckle::getTransformedProps(clusters = propeller_object_list[[xx]]$manual_annotation_1,
                                                      sample = as.character(propeller_object_list[[xx]]$clean_Sample_Name),
                                                      transform = "asin")
    
    # Looping through variants, adding genotype info to the Seurat object
    # Running genome-wide
    propeller_var_res <- lapply(unique(mashr_sighits$rowData$variant_id), function(variant){
        message(xx, variant)
        # variant <- "rs6067410"
        propeller_object_list[[xx]]@meta.data$variant <- NA
        propeller_object_list[[xx]]@meta.data$variant <- plyr::mapvalues(x = propeller_object_list[[xx]]@meta.data$clean_Sample_Name,
                                                                     from = rownames(gt_data),
                                                                     to = gt_data[,variant])
        propeller_object_list[[xx]]@meta.data$variant[-which(propeller_object_list[[xx]]@meta.data$clean_Sample_Name %in% rownames(gt_data))] <- NA
        
        # Creating the design matrix
        if(length(unique(propeller_object_list[[xx]]@meta.data$variant))<2){
            return(NULL)
        }
        sample_group <- distinct(propeller_object_list[[xx]]@meta.data[,which(colnames(propeller_object_list[[xx]]@meta.data) %in% c("clean_Sample_Name", "variant"))])
        design <- model.matrix(~0+sample_group$variant)
        colnames(design) <- gsub("sample_group\\$variant", "", colnames(design))
        
        not_na <- sample_group[complete.cases(sample_group),]$clean_Sample_Name
        if(length(not_na)<2){
            return(NULL)
        }
        
        transformed_props[[1]] <- transformed_props[[1]][,which(colnames(transformed_props[[1]]) %in% not_na)]
        transformed_props[[2]] <- transformed_props[[2]][,which(colnames(transformed_props[[2]]) %in% not_na)]
        transformed_props[[3]] <- transformed_props[[3]][,which(colnames(transformed_props[[3]]) %in% not_na)]
        
        # Running propeller
        res <- speckle::propeller.anova(prop.list = transformed_props,
                                        design = design,
                                        coef = 1:ncol(design),
                                        robust = TRUE,
                                        trend = FALSE,
                                        sort = TRUE)
        
        res$var <- variant
        res$celltype <- rownames(res)

        res
    })
    # Results to a single dataframe
    propeller_var_res <- propeller_var_res[!sapply(propeller_var_res,is.null)]
    propeller_var_res_df <- rbindlist(propeller_var_res, use.names=TRUE, fill=TRUE)
    
    write.table(propeller_var_res_df, paste0("/scratch/hnatri/ILD/", xx, "_sigeqtl_propeller_res.tsv"), sep = "\t", quote = F)
    
    propeller_var_res_df
})
names(propeller_res_list) <- propeller_object_names


