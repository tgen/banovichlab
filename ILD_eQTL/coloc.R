# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 08/21/2021
# Description: eQTL-GWAS colocalization analysis using coloc v4
# ==============================================================================

# ======================================
# Import libraries
# ======================================

#library(VariantAnnotation)
#library(LDlinkR)
#library(snpStats)
#library(susieR)
library(coloc)
library(biomaRt)
#library(patchwork)
library(googlesheets4)
library(optparse)
library(dplyr)
library(org.Hs.eg.db)
library(readJDX)
library(data.table)

# ======================================
# Helper functions
# ======================================

source("utilities.R")
source("lung_celltype_markers.R")
source("lung_celltype_colors.R")

# ======================================
# Environment variables
# ======================================

set.seed(1234)
setwd("/home/hnatri/ILD_eQTL/")

# ======================================
# Parsing command line arguments
# ======================================

option.list <- list(
    make_option("--celltype", type="character", action="store", default="epithelial_AT2", help="Celltype"),
    make_option("--gwas", type="character", action="store", default="allen_gwas", help="GWAS"),
    make_option("--eqtl_type", type="character", action="store", default="int", help="sc or int")
)

opt <- parse_args(OptionParser(option_list=option.list))
sinkall(filename = paste0("/scratch/hnatri/ILD/coloc_", opt$celltype, ":", opt$gwas, ".Rout"))

message(opt$celltype)
message(opt$gwas)
celltype <- c(opt$celltype)
gwas <- c(opt$gwas)

# ======================================
# Import eQTL and GWAS summary statistics, genotype data, and MAF info
# ======================================

# GWAS summary statistics
gwas_list <- list("ukbb_gwas" = "/labs/banovich/IPF/GWAS/ukbb_gwas_harmonized.tsv.gz",
                  "allen_gwas" = "/labs/banovich/IPF/GWAS/allen_gwas_harmonized.tsv.gz",
                  "gtex_lung" = "/labs/banovich/IPF/GWAS/gtex_lung_harmonized.tsv.gz",
                  "gtex_blood" = "/labs/banovich/IPF/GWAS/gtex_blood_harmonized.tsv.gz",
                  "gtex_brain" = "/labs/banovich/IPF/GWAS/gtex_brain_harmonized.tsv.gz",
                  "east_asian_gwas" = "/labs/banovich/IPF/GWAS/east_asian_gwas_harmonized.tsv.gz",
                  "asthma_adult_gwas" = "/labs/banovich/IPF/GWAS/asthma_adult_gwas_harmonized.tsv.gz",
                  "asthma_child_gwas" = "/labs/banovich/IPF/GWAS/asthma_child_gwas_harmonized.tsv.gz")

# Gene lists
new_gene_list <- lapply(seq(1, 7), function(i){
    genes <- readLines(paste0("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/2022-09-21_eGenesMulticelltypeClusterLists/2022-09-21_eQTL-heatmap-allBetas-nonUniqueclusterk", i, ".txt"))
    
    genes
})
new_genes <- unlist(new_gene_list)
length(new_genes)
new_genes <- sapply(strsplit(new_genes,"\\|"), `[`, 1)
new_genes <- unique(new_genes)
genes <- new_genes

# Ensembl IDs
genes <- unique(mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENSEMBL"))
genes <- as.character(na.omit(genes))

length(unique(genes))

# Significant mashR hits
mashr_sighits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-significant-eQTL.rds")

# Minor allele frequencies
# ILD MAF
ild_maf <- fread("/labs/banovich/IPF/eQTL/snpid_rsid.tsv")
colnames(ild_maf) <- c("chr", "pos", "rsid_snpid", "ref", "alt", "drop")
ild_maf$rsid <- sapply(strsplit(ild_maf$rsid_snpid,"\\;"), `[`, 1)
ild_maf$snp_id <- paste0("chr", ild_maf$chr, "_", ild_maf$pos, "_", ild_maf$ref, "_", ild_maf$alt)
ild_maf$rsid_alt <- paste0(ild_maf$rsid, "_", ild_maf$alt)
ild_maf$rsid_ref <- paste0(ild_maf$rsid, "_", ild_maf$ref)
mashr_snp_info <- fread("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/filtered_MAF-HWE-INDPW_snp-info.txt", header=T)

# Nominal stats
if (celltype %in% c("gtex_lung", "gtex_blood", "gtex_brain")){
    eqtl_nom_stats <- fread(gwas_list[[celltype]][1])
} else {
    path <- ifelse(opt$eqtl_type == "int", "/labs/banovich/IPF/eQTL/int-eQTL/", "/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/")
    eqtl_nom_stats <- fread(paste0(path, celltype, ".tsv"), header=T, sep="\t")
    eqtl_nom_stats <- distinct(eqtl_nom_stats)
    colnames(eqtl_nom_stats) <- gsub("snp_rsid", "rsid", colnames(eqtl_nom_stats))
    colnames(eqtl_nom_stats) <- gsub("gene_id", "feature_id", colnames(eqtl_nom_stats))
    
    limix_res <- fread(paste0("/labs/banovich/IPF/eQTL/limix_summaryStats_20220701_allSNPS/", celltype, "/qtl_results_all.txt"), header=T, sep="\t")
    colnames(limix_res) <- gsub("snp_id", "rsid", colnames(limix_res))
    limix_res <- merge(limix_res, ild_maf, by = "rsid") # MAF comes from the LIMIX output
    limix_maf <- distinct(limix_res[,c("snp_id", "rsid", "maf", "snp_chromosome", "snp_position", "ref", "alt")])
    
    # Adding MAF info from LIMIX outputs
    eqtl_nom_stats <- merge(eqtl_nom_stats, limix_maf, by = "rsid")
    
    # Sample size
    eqtl_sample_n <- mashr_sighits$colData[celltype, "n_samples"]
    
    # SE from SD
    eqtl_nom_stats$se <- eqtl_nom_stats$sds/sqrt(eqtl_sample_n)
    # Varbeta from SE
    eqtl_nom_stats$varbeta <- (eqtl_nom_stats$se*eqtl_nom_stats$se)
    
    # Assessed allele is the REF allele. Changing to match GWAS (REF = non-effect allele)
    # Adding allele info
    eqtl_nom_stats$snp_id <- paste0("chr", eqtl_nom_stats$snp_chromosome, "_", eqtl_nom_stats$snp_position, "_", eqtl_nom_stats$alt, "_", eqtl_nom_stats$ref)
    
    colnames(eqtl_nom_stats) <- gsub("posterior_means", "beta", colnames(eqtl_nom_stats))
    colnames(eqtl_nom_stats) <- gsub("lfsr", "p_value", colnames(eqtl_nom_stats))
    eqtl_nom_stats$n_samples <- eqtl_sample_n
    
    eqtl_nom_stats$rsid_ref <- paste0(eqtl_nom_stats$rsid, "_", eqtl_nom_stats$ref)
}


# All GWAS risk variants
# Preventing trying to get credentials
gs4_deauth()

gwas_table <- gs4_get("https://docs.google.com/spreadsheets/d/1ZqMRkxeNFE6Ln10IQfdkXlgjG3jeTFo8L5QSlr9I9yU/edit?usp=sharing")
gwas_table <- read_sheet(gwas_table, sheet = "Sheet1")

all_gwas_annotated_genes <- na.omit(as.character(sapply(strsplit(gwas_table$ENSEMBL, split=', ', fixed=TRUE), `[`, 1:2)))
length(all_gwas_annotated_genes)

# 91 nominally significant GWAS variant nearest genes
gwas_nearest_genes <- read.table("/labs/banovich/IPF/eQTL/GWAS_var_p1e-12_nearest_genes.tsv", header = T, sep = "\t")

# Making sure we're including genes annotated for GWAS risk variants
all_sig_egenes <- unique(c(genes, all_gwas_annotated_genes, unique(gwas_nearest_genes$ensembl)))

# Subsetting the eQTL dataset to only include selected genes
eqtl_nom_stats <- eqtl_nom_stats[which(eqtl_nom_stats$feature_id %in% all_sig_egenes),]

# All tested SNPs for significant eGenes
eqtl_snps <- unique(eqtl_nom_stats$snp_id)

# Sample sizes
eqtl_celltype_gwas_ninds <- rep(NA, 8)
names(eqtl_celltype_gwas_ninds) <- c("ukbb_gwas", "allen_gwas", "gtex_lung", "gtex_blood", "gtex_brain",
                                     "asthma_adult_gwas", "asthma_child_gwas", "copd_gwas")
eqtl_celltype_gwas_ninds["ukbb_gwas"] <- 449404
eqtl_celltype_gwas_ninds["allen_gwas"] <- 11259
eqtl_celltype_gwas_ninds["asthma_adult_gwas"] <- 327253
eqtl_celltype_gwas_ninds["asthma_child_gwas"] <- 314633
eqtl_celltype_gwas_ninds["east_asian_gwas"] <- 178020
eqtl_celltype_gwas_ninds["gtex_lung"] <- 515
eqtl_celltype_gwas_ninds["gtex_blood"] <- 670
eqtl_celltype_gwas_ninds["gtex_brain"] <- 205
eqtl_celltype_gwas_ninds

# n_samples for GTEx eQTLs
if(celltype %in% c("gtex_lung", "gtex_blood", "gtex_brain")){
    eqtl_nom_stats$n_samples <- eqtl_celltype_gwas_ninds[[celltype]]
}

# ======================================
# Running coloc
# ======================================

message("Running coloc")
run_coloc <- function(celltype, gwas, gwas_data, gene){
  # Data and sample sizes
  if (!(celltype %in% c("gtex_lung", "gtex_blood", "gtex_brain")) & length(unlist(strsplit(celltype, "_")))>1){
      ct <- unlist(strsplit(celltype, "_"))[2]
  } else {
      ct <- celltype
  }

  # Subsetting eQTL summary statistics
  gene_eqtl_stats <- eqtl_nom_stats[which(eqtl_nom_stats$feature_id == gene),]
  
  eqtl_N <- unique(gene_eqtl_stats$n_samples)
  gwas_N <- as.numeric(eqtl_celltype_gwas_ninds[gwas])
  
  # Subsetting GWAS summary statistics
  gene_gwas_stats <- gwas_data
  if (gwas %in% c("gtex_lung", "gtex_blood", "gtex_brain")){
      gene_gwas_stats <- gene_gwas_stats[which(gene_gwas_stats$feature_id == gene),]
  }
  
  if (nrow(gene_gwas_stats)<100){
      message("Not enough variants")
      return(NULL)
  }
  
  # Subsetting genotype data to retain the tested SNPs
  gene_eqtl_stats <- gene_eqtl_stats[which(gene_eqtl_stats$snp_id %in% gene_gwas_stats$snp_id),]
  gene_gwas_stats <- gene_gwas_stats[which(gene_gwas_stats$snp_id %in% gene_eqtl_stats$snp_id),]
  
  # Some rows are duplicated and not removed by distict() over the whole df?
  gene_eqtl_stats <- distinct(gene_eqtl_stats, snp_id, .keep_all = TRUE)
  
  # Matching SNP order
  gene_eqtl_stats <- gene_eqtl_stats[match(gene_gwas_stats$snp_id, gene_eqtl_stats$snp_id),]
  
  n_occur <- data.frame(table(gene_gwas_stats$snp_id))
  n_occur[n_occur$Freq > 1,]
  
  if (nrow(gene_eqtl_stats)<100){
    message("Not enough variants")
    return(NULL)
  }
  
  if (length(n_occur[n_occur$Freq > 1,]$Var1)>0){
    gene_eqtl_stats <- gene_eqtl_stats[-(which(gene_eqtl_stats$snp_id %in% n_occur[n_occur$Freq > 1,]$Var1)),]
    gene_gwas_stats <- gene_gwas_stats[-(which(gene_gwas_stats$snp_id %in% n_occur[n_occur$Freq > 1,]$Var1)),]
  }

  shared_snps <- gene_eqtl_stats$snp_id

  if (nrow(gene_eqtl_stats)<100){
    message("Not enough variants")
    return(NULL)
  }
  
  # Retaining variants where MAF is >0 but <1
  gene_eqtl_stats <- gene_eqtl_stats[which(gene_eqtl_stats$maf<1 & gene_eqtl_stats$maf>0),]
  gene_gwas_stats <- gene_gwas_stats[which(gene_gwas_stats$snp_id %in% gene_eqtl_stats$snp_id),]
  
  gene_gwas_stats <- gene_gwas_stats[which(gene_gwas_stats$maf<1 & gene_gwas_stats$maf>0),]
  gene_eqtl_stats <- gene_eqtl_stats[which(gene_eqtl_stats$snp_id %in% gene_gwas_stats$snp_id),]
  
  gene_eqtl_stats <- gene_eqtl_stats[match(gene_gwas_stats$snp_id, gene_eqtl_stats$snp_id),]
  
  if (nrow(gene_eqtl_stats)<100){
      message("Not enough variants")
      return(NULL)
  }
  
  # Constructing the data object for coloc
  datalist_names <- c("snp","position","MAF","type","N","pvalues")
  eqtl_data_list <- vector("list", length(datalist_names))
  names(eqtl_data_list) <- datalist_names
  eqtl_data_list$snp <- gene_eqtl_stats$snp_id
  eqtl_data_list$position <- gene_eqtl_stats$snp_position
  eqtl_data_list$MAF <- gene_eqtl_stats$maf
  eqtl_data_list$pvalues <- gene_eqtl_stats$p_value
  eqtl_data_list$type <- "quant"
  eqtl_data_list$N <- eqtl_N
  
  gwas_data_list <- vector("list", length(datalist_names))
  names(gwas_data_list) <- datalist_names
  gwas_data_list$snp <- gene_gwas_stats$snp_id
  gwas_data_list$position <- gene_gwas_stats$snp_position
  gwas_data_list$MAF <- gene_gwas_stats$maf
  gwas_data_list$pvalues <- gene_gwas_stats$p_value
  gwas_data_list$type <- "quant"
  gwas_data_list$N <- gwas_N
  
  # Basic coloc
  abf_res <- coloc.abf(dataset1=eqtl_data_list,
                       dataset2=gwas_data_list)
  
  return(abf_res)
}

# An empty list for the results
reslist_names <- paste0(celltype, ":", gwas)
res_list <- vector("list", length(reslist_names))
names(res_list) <- reslist_names

# For each cell type, running coloc with GWAS for significant eGenes
for (ct in celltype){
  sig_genes <- all_sig_egenes
  
  gene_res_list <- vector("list", length(sig_genes))
  names(gene_res_list) <- sig_genes
  
  for (gwas in gwas){
      if (ct==gwas){
          return(NULL)
      }
      # Importing GWAS data
      gwas_data <- fread(gwas_list[[gwas]])
      res_list[[paste0(ct, ":", gwas)]] <- gene_res_list
      for (gene in sig_genes){
        message(ct, " ", gwas, " ", gene)
        coloc_res <- run_coloc(ct, gwas, gwas_data, gene)
        
        res_list[[paste0(ct, ":", gwas)]][[gene]] <- coloc_res
    }
  }
}

names(res_list)

res_file <- paste0("/scratch/hnatri/ILD/mashr_coloc_abf_res_list_", celltype, ":", gwas, "_", opt$eqtl_type, ".rds")
saveRDS(res_list, res_file)

sinkall()
