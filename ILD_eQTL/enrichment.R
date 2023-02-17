# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 10/14/2022
# Description: eQTL-GWAS enrichment testing
# ==============================================================================

# ======================================
# Import libraries
# ======================================

library(biomaRt)
library(org.Hs.eg.db)
library(patchwork)
library(googlesheets4)
library(dplyr)
library(data.table)
library(nullranges)

# ======================================
# Helper functions
# ======================================

source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/ILD_processing_annotation/lung_celltype_markers.R")
source("/home/hnatri/ILD_processing_annotation/lung_celltype_colors.R")

# ======================================
# Environment variables
# ======================================

set.seed(1234)

# ======================================
# Import data, generate null distribution
# ======================================

# New gene lists
new_gene_list <- lapply(seq(1, 7), function(i){
    genes <- readLines(paste0("/scratch/hnatri/ILD/2022-09-21_eGenesMulticelltypeClusterLists/2022-09-21_eQTL-heatmap-allBetas-nonUniqueclusterk", i, ".txt"))
    genes <- as.data.frame(genes)
    genes$cluster <- i
    
    genes
})
new_genes <- do.call(rbind, new_gene_list)
new_genes$gene_symbol <- sapply(strsplit(new_genes$genes,"\\|"), `[`, 1)
new_genes$rsid <- sapply(strsplit(new_genes$genes,"\\|"), `[`, 2)

# Ensembl IDs
new_genes$ensembl <- mapIds(org.Hs.eg.db, keys = new_genes$gene_symbol, keytype = "SYMBOL", column = "ENSEMBL")

# SNP info
ild_maf <- fread("/labs/banovich/IPF/eQTL/snpid_rsid.tsv")
colnames(ild_maf) <- c("chr", "pos", "rsid_snpid", "ref", "alt", "drop")
ild_maf$rsid <- sapply(strsplit(ild_maf$rsid_snpid,"\\;"), `[`, 1)
ild_maf$snp_id <- paste0("chr", ild_maf$chr, "_", ild_maf$pos, "_", ild_maf$ref, "_", ild_maf$alt)
ild_maf$rsid_alt <- paste0(ild_maf$rsid, "_", ild_maf$alt)
ild_maf$rsid_ref <- paste0(ild_maf$rsid, "_", ild_maf$ref)
mashr_snp_info <- fread("/scratch/hnatri/ILD/2022-08-10_38celltypes-mashr/filtered_MAF-HWE-INDPW_snp-info.txt", header=T)
mashr_sighits <- readRDS("/scratch/hnatri/ILD/2022-08-10_38celltypes-mashr/mashr_summary_stats-significant-eQTL.rds")
mashr_tophits <- readRDS("/scratch/hnatri/ILD/2022-08-10_38celltypes-mashr/mashr_summary_stats-top-eQTL.rds")
mashr_sighits$rowData$snp_id <- paste0("chr", mashr_sighits$rowData$snp_id)

new_genes <- merge(new_genes, mashr_snp_info, by = "rsid")

# TSS info
genes_gtf <- read.table("/labs/banovich/SingleCell/CellRanger/3_1_0/Ensemble_93/PipelineData/Projects/IPF/References/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf", sep = "\t")
gene_ids <- sapply(strsplit(genes_gtf$V9," "), `[`, 2)
gene_ids <- gsub(";", "", gene_ids)
genes_gtf$feature_id <- gene_ids
genes_gtf <- dplyr::filter(genes_gtf, V3 == "gene")
genes_gtf <- dplyr::select(genes_gtf, V1, V4, V5, feature_id)
colnames(genes_gtf) <- c("chr", "start", "stop", "ensembl")

new_genes <- merge(new_genes, genes_gtf, by = "ensembl")
new_genes$distance_to_TSS <- new_genes$start-new_genes$snp_loc

# mashr stats, keeping non-significant only
celltypes <- readLines("/home/hnatri/ILD_eQTL/coloc_mashr_celltypes.tsv")
celltypes <- c("endothelial_Lymphatic", "epithelial_AT2", "epithelial_AT1", "epithelial_Basal", "immune_Alveolarmacrophage", "immune_Bcells", "immune_Monocyte", "mesenchymal_MatrixFB")
mashr_res_list <- lapply(celltypes, function(celltype){
    eqtl_nom_stats <-  fread(paste0("/scratch/hnatri/ILD/2022-08-10_38celltypes-mashr/", celltype, ".tsv"), header=T, sep="\t")
    eqtl_nom_stats$snp_id <- paste0("chr", eqtl_nom_stats$snp_id)
    colnames(eqtl_nom_stats) <- gsub("feature_id", "ensembl", colnames(eqtl_nom_stats))
    
    # eQTLs not significant for any celltype
    eqtl_nom_stats$ensembl_rsid <- paste0(eqtl_nom_stats$ensembl, "|", eqtl_nom_stats$snp_rsid)
    eqtl_nom_stats <- eqtl_nom_stats[-which(eqtl_nom_stats$ensembl_rsid %in% rownames(mashr_sighits$rowData)),]
    #length(intersect(eqtl_nom_stats$snp_id, mashr_snp_info$snp_id)    # TSS info
    eqtl_nom_stats <- merge(eqtl_nom_stats, genes_gtf, by = "ensembl")
    eqtl_nom_stats <- merge(eqtl_nom_stats, mashr_snp_info, by = "snp_id")
    eqtl_nom_stats$distance_to_TSS <- eqtl_nom_stats$start-eqtl_nom_stats$snp_loc
    eqtl_nom_stats <- eqtl_nom_stats[,c("ensembl", "rsid", "snp_id", "chr", "snp_loc", "start", "distance_to_TSS")]
    
    eqtl_nom_stats
})

mashr_res_nonsig <- do.call(rbind, mashr_res_list)
mashr_res_nonsig <- dplyr::distinct(mashr_res_nonsig)

# Null set of eQTLs
colnames(new_genes) <- gsub("start", "gene_start", colnames(new_genes))
colnames(mashr_res_nonsig) <- gsub("start", "gene_start", colnames(mashr_res_nonsig))

new_genes$group <- "sig"
mashr_res_nonsig$group <- "nonsig"

mashr_sighits$rowData$group <- "sig"
colnames(mashr_sighits$rowData) <- gsub("snp_chr", "chr", colnames(mashr_sighits$rowData))
colnames(mashr_sighits$rowData) <- gsub("feature_id", "ensembl", colnames(mashr_sighits$rowData))

mashr_sighits_df <- as.data.frame(mashr_sighits$rowData)
mashr_sighits_df <- merge(mashr_sighits_df, genes_gtf, by = "ensembl")
colnames(mashr_sighits_df) <- gsub("start", "gene_start", colnames(mashr_sighits_df))
colnames(mashr_sighits_df) <- gsub("chr.x", "chr", colnames(mashr_sighits_df))

# Taking a smaller subset
#mashr_sighits_subset <- mashr_sighits_df[sample(rownames(mashr_sighits_df), size = 10000),]

all_eqtls <- plyr::rbind.fill(mashr_sighits_df, mashr_res_nonsig)

# IPF GWAS meta-analysis nominal stats
ipf_gwas <- fread("/labs/banovich/IPF/GWAS/allen_gwas_harmonized.tsv.gz")
ipf_gwas <- ipf_gwas[which(ipf_gwas$snp_id %in% mashr_snp_info$snp_id),]
# Selecting nominally significant hits with a relaxed threshold
ipf_gwas_sig <- ipf_gwas[which(ipf_gwas$p_value < 1e-6),]

all_eqtls <- all_eqtls[which(all_eqtls$snp_id %in% ipf_gwas$snp_id),]
all_eqtls$abs_distance_to_TSS <- abs(all_eqtls$gene_start-all_eqtls$snp_loc)
all_eqtls$abs_distance_to_TSS[is.na(all_eqtls$abs_distance_to_TSS)]
all_eqtls$chr <- as.numeric(all_eqtls$chr)

# Need to remove rownames from the dataframe, otherwise the resulting mcols
# object will be invalid
rownames(all_eqtls) <- NULL

all_eqtls_granges <- makeGRangesFromDataFrame(all_eqtls,
                                              keep.extra.columns=T,
                                              ignore.strand=T,
                                              seqinfo=NULL,
                                              seqnames.field="chr",
                                              start.field="snp_loc",
                                              end.field="snp_loc")

null_set <- matchRanges(focal = all_eqtls_granges[all_eqtls_granges$group=="sig"],
                        pool = all_eqtls_granges[all_eqtls_granges$group=="nonsig"],
                        covar = ~abs_distance_to_TSS)
                        #method = "stratified",
                        #replace = TRUE)

plotPropensity(null_set)

# ======================================
# Enrichment
# ======================================

eqtl_gwas_snp_overlap <- as.data.frame(matrix(nrow=length(tested_shared_snps), ncol=15))
rownames(eqtl_gwas_snp_overlap) <- tested_shared_snps
colnames(eqtl_gwas_snp_overlap) <- c("multistate", "global", "unique", "nonsig", "anysig", "immune", "epithelial", "k1", "k2", "k3", "k4", "k5", "k6", "k7", "gwas")

eqtl_gwas_snp_overlap$multistate <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% mashr_sighits$rowData[which(mashr_sighits$rowData$type=="multi-state"),]$rsid, 1, 0)
eqtl_gwas_snp_overlap$global <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% mashr_sighits$rowData[which(mashr_sighits$rowData$type=="global"),]$rsid, 1, 0)
eqtl_gwas_snp_overlap$unique <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% mashr_sighits$rowData[which(mashr_sighits$rowData$type=="unique"),]$rsid, 1, 0)
eqtl_gwas_snp_overlap$anysig <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% mashr_sighits$rowData$rsid, 1, 0)
eqtl_gwas_snp_overlap$nonsig <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% null_set$rsid, 1, 0)
eqtl_gwas_snp_overlap$gwas <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% ipf_gwas_sig$rsid, 1, 0)

eqtl_gwas_snp_overlap$immune <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% new_genes[which(new_genes$cluster %in% c(4, 6, 7)),]$rsid, 1, 0)
eqtl_gwas_snp_overlap$epithelial <- ifelse(rownames(eqtl_gwas_snp_overlap) %in% new_genes[which(new_genes$cluster %in% c(3, 5)),]$rsid, 1, 0)

# Adding cluster eQTLs
for(i in seq(1, 7)){
    snps <- new_genes[which(new_genes$cluster==i),]$rsid
    eqtl_gwas_snp_overlap[which(rownames(eqtl_gwas_snp_overlap) %in% snps), c(paste0("k", i))] <- 1
    eqtl_gwas_snp_overlap[-which(rownames(eqtl_gwas_snp_overlap) %in% snps), c(paste0("k", i))] <- 0
}

head(eqtl_gwas_snp_overlap)

# Are lineage or cell type specific eQTLs more likely to to be significant in GWAS?

# A matrix with numbers/proportions of GWAS and non-GWAS SNPs for each class
prop_matrix <- matrix(NA, nrow = 14, ncol = 4)
colnames(prop_matrix) <- c("gwas", "non_gwas", "prop_gwas", "prop_non_gwas")
rownames(prop_matrix) <- colnames(eqtl_gwas_snp_overlap[1:14])

for(i in rownames(prop_matrix)){
    prop_matrix[i, "gwas"] <- table(eqtl_gwas_snp_overlap[,i], eqtl_gwas_snp_overlap$gwas)[2,2]
    prop_matrix[i, "non_gwas"] <- table(eqtl_gwas_snp_overlap[,i], eqtl_gwas_snp_overlap$gwas)[2,1]
}

prop_matrix
prop_matrix <- as.data.frame(prop_matrix)
prop_matrix$prop_gwas <- prop_matrix$gwas/(prop_matrix$gwas+prop_matrix$non_gwas)*100
prop_matrix$prop_non_gwas <- prop_matrix$non_gwas/(prop_matrix$gwas+prop_matrix$non_gwas)*100

res <- fisher.test(table(eqtl_gwas_snp_overlap[,"global"], eqtl_gwas_snp_overlap$gwas))
res$estimate

# Adding Fisher's test results
for(i in rownames(prop_matrix)){
    prop_matrix[i, "fishers_p"] <- fisher.test(table(eqtl_gwas_snp_overlap[,i], eqtl_gwas_snp_overlap$gwas))$p.value
    prop_matrix[i, "fishers_or"] <- fisher.test(table(eqtl_gwas_snp_overlap[,i], eqtl_gwas_snp_overlap$gwas))$estimate
}

