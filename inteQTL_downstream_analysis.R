# ==============================================================================
# Author(s) : Heini M. Natri, hnatri@tgen.org
# Date: 12/15/2022
# Description: ILD int-eQTL downstream analysis
# ==============================================================================

# ======================================
# Import libraries
# ======================================

suppressMessages({library(plyr)
                  library(dplyr)
                  library(qdap)
                  library(data.table)
                  library(tidyverse)
                  library(Seurat)
                  library(scCustomize)
                  library(googlesheets4)
                  library(QTLExperiment)
                  library(ComplexHeatmap)})

# ======================================
# Helper functions
# ======================================

source("/home/hnatri/Utilities/utilities.R")

# ======================================
# Environment variables
# ======================================

set.seed(1234)

# Colors
gs4_deauth()
paper_celltype_colors  <- gs4_get("https://docs.google.com/spreadsheets/d/1NnLlftpsqq_JEaINBMOe2T4qB1_XPihrcl1mNDim8B4/edit?usp=sharing")
sheet_names(paper_celltype_colors)
paper_celltype_colors <- read_sheet(paper_celltype_colors, sheet = "Sheet1")

paper_celltype_colors[which(paper_celltype_colors[,c("Pretty cell type name")]=="cDC2"),]

paper_celltype_colors$`Cell type level` <- gsub(" ", "", paper_celltype_colors$`Cell type level`)

lineage_cols <- paper_celltype_colors$lineage_color
names(lineage_cols) <- paper_celltype_colors$Lineage
celltype_cols <- paper_celltype_colors$celltype_color
names(celltype_cols) <- paper_celltype_colors$`Pretty cell type name`
level2_color <- paper_celltype_colors$level2_color
names(level2_color) <- paper_celltype_colors$Annotation_level_2

# ======================================
# Import data
# ======================================

# Seurat objects
object_list <- readRDS("/labs/banovich/IPF/eQTL/object_list.rds")
object_names <- names(object_list)

unique(object_list[[1]]$Status)
object_list_disease <- lapply(object_list, function(xx){
    subset(xx, subset = Status == "Disease")
})
names(object_list_disease) <- object_names

object_list_ctrl <- lapply(object_list, function(xx){
    subset(xx, subset = Status == "Control")
})
names(object_list_ctrl) <- object_names

# DEGs between ILD and CTRL with a relaxed threshold
deg_list <- lapply(object_list, function(object){
    object_deg_list <- lapply(unique(object$manual_annotation_1), function(xx){
    message(xx)
    object_subset <- subset(object, subset = manual_annotation_1 == xx)
    if (length(unique(object_subset$Status))<2){
        return(NULL)
    } else {
        markers <- presto::wilcoxauc(object_subset,
                                     group_by = "Status",
                                     groups_use = c("Disease", "Control"),
                                     assay = "data",
                                     seurat_assay = "RNA")
        markers$celltype <- xx
        
        return(markers)
    }
})
    names(object_deg_list) <- unique(object$manual_annotation_1)
    object_deg_list[sapply(object_deg_list, is.null)] <- NULL
    
    object_deg_df <- as.data.frame(do.call(rbind, object_deg_list))
    object_deg_df <- object_deg_df[which(object_deg_df$group=="Disease"),]
    
    object_deg_df
})

deg_df <- as.data.frame(do.call(rbind, deg_list))

#write.table(deg_df, "/scratch/hnatri/ILD/deg_df.tsv", sep = "\t", quote = F, row.names = F)
#deg_df <- read.table("/scratch/hnatri/ILD/deg_df.tsv", sep = "\t", header = T)
deg_df$gene_celltype <- paste0(deg_df$feature, "_", deg_df$celltype)
degs_sig_df <- deg_df[which(deg_df$padj<0.1),]
#degs_sig_df <- degs_sig_df[which(degs_sig_df$pct_in>0.5 | degs_sig_df$pct_out>0.5),]

# int-eQTL summary stats
inteqtl_all <- readRDS("/labs/banovich/IPF/eQTL/20220830_interaction-eQTL-mashr/mashr_applied_significant_ieQTL_all.rds")
inteqtl_top <- readRDS("/labs/banovich/IPF/eQTL/20220830_interaction-eQTL-mashr/mashr_applied_significant_ieQTL_top.rds")

inteqtl_top_sig <- as.data.frame(assay(inteqtl_top, "significant"))
inteqtl_top_sig$gene_snp <- rownames(inteqtl_top_sig)
inteqtl_top_sig <- pivot_longer(inteqtl_top_sig, cols = colnames(assay(inteqtl_top, "significant")),
                                names_to = "celltype", values_to = "significant")
inteqtl_top_sig$gene <- sapply(strsplit(inteqtl_top_sig$gene_snp, split='|', fixed=TRUE), `[`, 1)
inteqtl_top_sig$snp <- sapply(strsplit(inteqtl_top_sig$gene_snp, split='|', fixed=TRUE), `[`, 2)

setdiff(inteqtl_top_sig$celltype, paper_celltype_colors$`Pretty cell type name (old)`)
inteqtl_top_sig$pretty_celltype <- plyr::mapvalues(x = inteqtl_top_sig$celltype,
                                                   from = paper_celltype_colors$`Pretty cell type name (old)`,
                                                   to = paper_celltype_colors$`Pretty cell type name`)

inteqtl_top_sig$gene_celltype <- paste0(inteqtl_top_sig$gene, "_", inteqtl_top_sig$celltype)
inteqtl_top_sig <- inteqtl_top_sig[which(inteqtl_top_sig$significant==TRUE),]

# MAFs for each cell type
maf_res_df <- read.table("/scratch/hnatri/ILD/maf_res_df.tsv", sep = "\t", header = T)
maf_res_df_wide <- pivot_wider(maf_res_df, id_cols = c("snp", "celltype"), values_from = "maf", names_from = "group")
maf_res_df_wide <- maf_res_df_wide[which(maf_res_df_wide$ILD>0.05 & maf_res_df_wide$CTRL>0.05),]

maf_res_df_wide$pretty_celltype <- plyr::mapvalues(x = maf_res_df_wide$celltype,
                from = paper_celltype_colors$`mashR result name`,
                to = paper_celltype_colors$`Pretty cell type name`)
maf_res_df_wide$snp_celltype <- paste0(maf_res_df_wide$snp, "_", maf_res_df_wide$pretty_celltype)

# Dropping SNPs where MAF for ILD or CTRL is <0.05
inteqtl_top_sig$snp_celltype <- paste0(inteqtl_top_sig$snp, "_", inteqtl_top_sig$celltype)
inteqtl_top_sig <- inteqtl_top_sig[which(inteqtl_top_sig$snp_celltype %in% maf_res_df_wide$snp_celltype),]

# Saving to a file
#write.table(inteqtl_top_sig, "/scratch/hnatri/ILD/inteqtl_sig_maf_filtered.tsv", sep = "\t", quote = F, row.names = F)
#inteqtl_top_sig <- read.table("/scratch/hnatri/ILD/inteqtl_sig_maf_filtered.tsv", sep = "\t", header = T)
#inteqtl_top_sig <- distinct(inteqtl_top_sig[,2:4])

# For each eQTL, finding the top SNP for that cell type (the top object includes
# all significant SNPs)
lfsr_df <- as.data.frame(assay(inteqtl_top, "lfsr"))
lfsr_df$gene_snp <- rownames(lfsr_df)
lfsr_df <- pivot_longer(lfsr_df, cols = colnames(assay(inteqtl_top, "lfsr")),
                        names_to = "celltype", values_to = "lfsr")
lfsr_df$gene <- sapply(strsplit(lfsr_df$gene_snp, split='|', fixed=TRUE), `[`, 1)
lfsr_df$snp <- sapply(strsplit(lfsr_df$gene_snp, split='|', fixed=TRUE), `[`, 2)

# Merging
inteqtl_top_sig$gene_snp_celltype <- paste0(inteqtl_top_sig$gene_snp, "_", inteqtl_top_sig$celltype)
lfsr_df$gene_snp_celltype <- paste0(lfsr_df$gene_snp, "_", lfsr_df$celltype)
inteqtl_top_sig <- merge(inteqtl_top_sig, lfsr_df, by = "gene_snp_celltype")

# Filtering to retain top SNPs
inteqtl_top_sig <- as.data.frame(inteqtl_top_sig)
inteqtl_topsnp_sig <- inteqtl_top_sig %>% group_by(gene_celltype) %>%
    dplyr::slice(which.min(lfsr))

inteqtl_topsnp_sig <- inteqtl_topsnp_sig[,c("gene_snp_celltype", "gene_snp.x", "celltype.x", "gene.x", "snp.x", "pretty_celltype", "gene_celltype", "lfsr")]
colnames(inteqtl_topsnp_sig) <- c("gene_snp_celltype", "gene_snp", "celltype", "gene", "snp", "pretty_celltype", "gene_celltype", "lfsr")

#write.table(inteqtl_topsnp_sig, "/labs/banovich/IPF/eQTL/inteQTL_topSNPs.tsv", quote = F, sep = "\t", row.names = F)
#inteqtl_topsnp_sig <- read.table("/labs/banovich/IPF/eQTL/inteQTL_topSNPs.tsv", sep = "\t", header = T)

#lfsr_df$gene <- sapply(strsplit(lfsr_df$gene_snp,"\\|"), `[`, 1)
#egene_plot_data_005 <- lfsr_df[which(lfsr_df$lfsr<0.05),]
#
## If eQTL is significant in any celltype with lfsr<0.05, it will be considered
## significant in other cell types with lfsr<0.1
#egene_plot_data_010 <- lfsr_df[which(lfsr_df$lfsr<0.10),]
#egene_plot_data_010 <- egene_plot_data_010[which(egene_plot_data_010$gene_snp %in% egene_plot_data_005$gene_snp),]

# Barplots of top int-eQTL and eGenes
eqtl_plot_data <- as.data.frame(table(inteqtl_topsnp_sig$pretty_celltype))
#eqtl_plot_data <- as.data.frame(table(inteqtl_top_sig$pretty_celltype))
eqtl_barplot <- ggplot(eqtl_plot_data, aes(x=reorder(Var1, -Freq), y=Freq, fill = Var1)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name = "Cell type", values = celltype_cols) +
    theme_bw() +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("") +
    ylab("# top int-eQTL") +
    #theme(legend.text = element_text(size=6)) +
    theme(plot.margin = margin(10,10,25,30, "pt"))

filename <- "/home/hnatri/ILD_eQTL/top_inteqtl_barplot.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 8,
    height = 5)

eqtl_barplot

dev.off()

egene_plot_data <- distinct(inteqtl_topsnp_sig[,c("pretty_celltype", "gene.x")])
egene_plot_data <- as.data.frame(table(egene_plot_data$pretty_celltype))

egene_barplot <- ggplot(egene_plot_data, aes(x=reorder(Var1, -Freq), y=Freq, fill = Var1)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name = "Cell type", values = celltype_cols) +
    theme_bw() +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("") +
    ylab("# int-eGenes") +
    #theme(legend.text = element_text(size=6)) +
    theme(plot.margin = margin(10,10,25,30, "pt"))

filename <- "/home/hnatri/ILD_eQTL/int-egene_barplot.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 8,
    height = 5)

egene_barplot

dev.off()

# Out of the significant eGene-cell type pairs, which genes were DE between ILD
# and CTRL in that cell type?
setdiff(degs_sig_df$celltype, inteqtl_topsnp_sig$celltype)
setdiff(inteqtl_topsnp_sig$celltype, degs_sig_df$celltype)
#degs_sig_df$celltype <- gsub("Differentiating Ciliated", "Differentiating ciliated", degs_sig_df$celltype)

length(intersect(degs_sig_df$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# What proportion were expressed by >30% cells in both groups?
length(intersect(deg_df[which(deg_df$pct_in>30 & deg_df$pct_out>30),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# ...and had logFC below threshold?
length(intersect(deg_df[which(deg_df$pct_in>30 & deg_df$pct_out>30 & abs(deg_df$logFC)<0.2),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# Which cell types were these?
exp_both_small_logfc <- inteqtl_topsnp_sig[which(inteqtl_topsnp_sig$gene_celltype %in% deg_df[which(deg_df$pct_in>30 & deg_df$pct_out>30 & abs(deg_df$logFC)<0.2),]$gene_celltype),]

#write.table(exp_both_small_logfc, "/scratch/hnatri/ILD/exp_both_small_logfc.tsv", sep = "\t", quote = F, row.names = F)
exp_both_small_logfc <- read.table("/scratch/hnatri/ILD/exp_both_small_logfc.tsv", sep = "\t", header = T)
exp_both_small_logfc[which(exp_both_small_logfc$gene.x=="MUC5B"),]

exp_both_small_logfc[which(exp_both_small_logfc$gene.x=="DSP"),]$snp.x
exp_both_small_logfc[which(exp_both_small_logfc$gene.x=="DSP" & exp_both_small_logfc$snp.x=="rs1358909"),]$celltype.x

exp_both_small_logfc_plot <- as.data.frame(table(exp_both_small_logfc$celltype.x))

exp_both_small_logfc_barplot <- ggplot(exp_both_small_logfc_plot, aes(x=reorder(Var1, -Freq), y=Freq, fill = Var1)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name = "Cell type", values = celltype_cols) +
    theme_bw() +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("") +
    ylab("# int-eGenes") +
    #theme(legend.text = element_text(size=6)) +
    theme(plot.margin = margin(10,10,25,30, "pt"))

filename <- "/home/hnatri/ILD_eQTL/int-egene_exp_both_small_logfc_barplot.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 8,
    height = 5)

exp_both_small_logfc_barplot

dev.off()

# ...and had logFC above threshold?
length(intersect(deg_df[which(deg_df$pct_in>30 & deg_df$pct_out>30 & abs(deg_df$logFC)>1),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# What proportion were expressed by >30% in either group and had logFC above threshold?
length(intersect(deg_df[which(deg_df$pct_in>30 | deg_df$pct_out>30 & abs(deg_df$logFC)>2),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# Less than 30% in one or both groups?
length(intersect(deg_df[which(deg_df$pct_in<30 | deg_df$pct_out<30),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# Expressed by >30% in ILD or CTRL only?
length(intersect(deg_df[which(deg_df$pct_in<30 & deg_df$pct_out>30),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))
length(intersect(deg_df[which(deg_df$pct_in>30 & deg_df$pct_out<30),]$gene_celltype, inteqtl_topsnp_sig$gene_celltype))/length(unique(inteqtl_topsnp_sig$gene_celltype))

# Relationship between int-eQTL beta and logFC sign
inteqtl_topsnp_sig_de <- merge(inteqtl_topsnp_sig, deg_df, by = "gene_celltype")
inteqtl_topsnp_sig_de_beta <- merge(inteqtl_topsnp_sig_de, beta_df_long, by = "gene_snp_celltype")
plot(x=inteqtl_topsnp_sig_de_beta$logFC, y=inteqtl_topsnp_sig_de_beta$value,
     xlab = "logFC", ylab = "beta")

inteqtl_top_sig_small <- distinct(inteqtl_topsnp_sig[,c("gene_snp.x", "celltype.x")])
inteqtl_topsnp_sig_de <- inteqtl_topsnp_sig[which(inteqtl_topsnp_sig$gene_celltype %in% degs_sig_df$gene_celltype),]

# Plotting basic stats
plot(inteqtl_top@colData$n_samples, inteqtl_top@colData$nSignificant)

top_df <- as.data.frame(inteqtl_top@colData)

# All significant int-eQTLs
ggplot(data = top_df, aes(x = top_df$n_samples, y = top_df$nSignificant, color = top_df$Annotation_level_2)) +
    geom_point() +
    scale_color_manual(name = "Annotation", values = level2_color) +
    theme_bw() + 
    ylab("# significant int-eQTL") +
    xlab("# samples")

lm(top_df$nSignificant ~ top_df$n_samples)    
summary(lm(top_df$nSignificant ~ top_df$n_samples))

# Finding the # int-eQTL target genes for each celltype
lsfr_df <- as.data.frame(assay(inteqtl_top, "lfsr"))
lsfr_df$gene_snp <- rownames(lsfr_df)
lsfr_df_long <- pivot_longer(lsfr_df,
                             cols = colnames(assay(inteqtl_top, "lfsr")),
                             names_to = "celltype")
lsfr_df_long$gene <- sapply(strsplit(lsfr_df_long$gene_snp, split='|', fixed=TRUE), `[`, 1)
lsfr_df_long$snp <- sapply(strsplit(lsfr_df_long$gene_snp, split='|', fixed=TRUE), `[`, 2)
#lsfr_df_long$celltype <- gsub("Differentiating ciliated", "Differentiating Ciliated", lsfr_df_long$celltype)
lsfr_df_long$gene_celltype <- paste0(lsfr_df_long$gene, "_", lsfr_df_long$celltype)
lsfr_df_long$snp_celltype <- paste0(lsfr_df_long$snp, "_", lsfr_df_long$celltype)
#lsfr_df_long_sig <- lsfr_df_long[which(lsfr_df_long$value<0.01),]
lsfr_df_long <- lsfr_df_long[which(lsfr_df_long$snp_celltype %in% maf_res_df_wide$snp_celltype),]
n_sig_egenes <- sapply(unique(lsfr_df_long_sig$celltype), function(xx){
    length(unique(lsfr_df_long_sig[which(lsfr_df_long_sig$celltype==xx),]$gene))
})

n_sig_egenes <- as.data.frame(n_sig_egenes)
n_sig_egenes$Pretty.cell.type.name <- rownames(n_sig_egenes)

# What proportion of int-eGenes were expressed more in ILD?
head(inteqtl_topsnp_sig_de)
inteqtl_topsnp_sig_de_sig <- inteqtl_topsnp_sig_de[which(inteqtl_topsnp_sig_de$pval<0.1),]
dim(inteqtl_topsnp_sig_de_sig)
dim(inteqtl_topsnp_sig_de_sig[which(inteqtl_topsnp_sig_de_sig$logFC>0),])

top_df <- merge(top_df, n_sig_egenes, by = "Pretty.cell.type.name")

ggplot(data = top_df, aes(x = top_df$n_samples, y = top_df$n_sig_egenes.x, color = top_df$Annotation_level_2)) +
    geom_point() +
    scale_color_manual(name = "Annotation", values = level2_color) +
    theme_bw() + 
    ylab("# significant int-eGenes") +
    xlab("# samples")

lm(top_df$n_sig_egenes.x ~ top_df$n_samples)    
summary(lm(top_df$n_sig_egenes.x ~ top_df$n_samples))

# Distance from TSS and effect sizes
beta_df <- as.data.frame(assay(inteqtl_top, "betas"))
beta_df$gene_snp <- rownames(beta_df)
beta_df_long <- pivot_longer(beta_df,
                             cols = colnames(assay(inteqtl_top, "betas")),
                             names_to = "celltype")
beta_df_long$gene <- sapply(strsplit(beta_df_long$gene_snp, split='|', fixed=TRUE), `[`, 1)
beta_df_long$snp <- sapply(strsplit(beta_df_long$gene_snp, split='|', fixed=TRUE), `[`, 2)
#beta_df_long$celltype <- gsub("Differentiating ciliated", "Differentiating Ciliated", beta_df_long$celltype)
beta_df_long$snp_celltype <- paste0(beta_df_long$snp, "_", beta_df_long$celltype)
beta_df_long$gene_celltype <- paste0(beta_df_long$gene, "_", beta_df_long$celltype)
beta_df_long$gene_snp_celltype <- paste0(beta_df_long$gene_snp, "_", beta_df_long$celltype)

# MAF filtered and top SNPs only
beta_df_long_sig <- beta_df_long[which(beta_df_long$gene_snp_celltype %in% inteqtl_topsnp_sig$gene_snp_celltype),]

inteqtl_topsnp_sig$beta <- plyr::mapvalues(x = inteqtl_topsnp_sig$gene_snp_celltype,
                                        from = beta_df_long_sig$gene_snp_celltype,
                                        to = beta_df_long_sig$value)

# Adding distance to TSS
# TSS info
genes_gtf <- read.table("/labs/banovich/SingleCell/CellRanger/3_1_0/Ensemble_93/PipelineData/Projects/IPF/References/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf", sep = "\t")
head(genes_gtf)
gene_ids <- sapply(strsplit(genes_gtf$V9," "), `[`, 2)
gene_ids <- gsub(";", "", gene_ids)
gene_names <- sapply(strsplit(genes_gtf$V9," "), `[`, 6)
gene_names <- gsub(";", "", gene_names)
genes_gtf$feature_id <- gene_ids
genes_gtf$feature_name <- gene_names
genes_gtf <- dplyr::filter(genes_gtf, V3 == "gene")
genes_gtf <- dplyr::select(genes_gtf, V1, V4, V5, feature_id, feature_name)
colnames(genes_gtf) <- c("chr", "start", "stop", "ensembl", "gene")

# SNP info
ild_maf <- fread("/labs/banovich/IPF/eQTL/snpid_rsid.tsv")
colnames(ild_maf) <- c("chr", "pos", "rsid_snpid", "ref", "alt", "drop")
ild_maf$rsid <- sapply(strsplit(ild_maf$rsid_snpid,"\\;"), `[`, 1)
ild_maf$snp_id <- paste0("chr", ild_maf$chr, "_", ild_maf$pos, "_", ild_maf$ref, "_", ild_maf$alt)
ild_maf$rsid_alt <- paste0(ild_maf$rsid, "_", ild_maf$alt)
ild_maf$rsid_ref <- paste0(ild_maf$rsid, "_", ild_maf$ref)
mashr_snp_info <- fread("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/filtered_MAF-HWE-INDPW_snp-info.txt", header=T)

# Adding # of significant cell types, SNP info, gene info
sig_gene_vars <- as.data.frame(table(inteqtl_topsnp_sig$gene_snp))
colnames(sig_gene_vars) <- c("gene_snp", "n_celltypes")

inteqtl_topsnp_sig <- merge(inteqtl_topsnp_sig, sig_gene_vars, by = "gene_snp")
inteqtl_topsnp_sig$snp_celltype <- paste0(inteqtl_topsnp_sig$snp, "_", inteqtl_topsnp_sig$celltype)
#beta_df_long_sig <- beta_df_long_sig[which(beta_df_long_sig$snp_celltype %in% maf_res_df_wide$snp_celltype),]
hist(inteqtl_topsnp_sig$n_celltypes)
table(inteqtl_topsnp_sig$n_celltypes)
inteqtl_topsnp_sig$type <- ifelse(inteqtl_topsnp_sig$n_celltypes==1, "unique", "multi-state")

colnames(inteqtl_topsnp_sig) <- gsub("snp", "rsid", colnames(inteqtl_topsnp_sig))
inteqtl_topsnp_sig <- merge(inteqtl_topsnp_sig, mashr_snp_info, by = "rsid")
inteqtl_topsnp_sig <- merge(inteqtl_topsnp_sig, genes_gtf, by = "gene")

inteqtl_topsnp_sig$abs_distance_TSS <- abs(inteqtl_topsnp_sig$start - inteqtl_topsnp_sig$snp_loc)

#write.table(inteqtl_topsnp_sig, "/scratch/hnatri/ILD/inteQTL_topSNPs.tsv", sep = "\t", quote = F, row.names = F)
inteqtl_topsnp_sig <- read.table("/scratch/hnatri/ILD/inteQTL_topSNPs.tsv", sep = "\t", header = T)

colnames(inteqtl_topsnp_sig) <- gsub("gene_rsid_celltype", "gene_snp_celltype", colnames(inteqtl_topsnp_sig))

library(ggpubr)
comparisons <- list(#c("global", "multi-state"),
                    #c("global", "unique"),
                    c("multi-state", "unique"))

t.test(abs(beta_df_long_sig_top[which(beta_df_long_sig_top$type=="multi-state"),]$abs_distance_TSS),
       abs(beta_df_long_sig_top[which(beta_df_long_sig_top$type=="unique"),]$abs_distance_TSS))

t.test(abs(beta_df_long_sig_top[which(beta_df_long_sig_top$type=="multi-state"),]$beta),
       abs(beta_df_long_sig_top[which(beta_df_long_sig_top$type=="unique"),]$beta))

beta_df_long_sig_top <- as.data.frame(beta_df_long_sig_top)
beta_df_long_sig_top$type <- factor(beta_df_long_sig_top$type, levels = c("unique", "multi-state"))

p1 <- ggplot(beta_df_long_sig_top, aes(x=type, y=abs(abs_distance_TSS)/1e+06)) +
    geom_violin(scale="width") + 
    geom_boxplot(width=0.1, outlier.shape = NA) +
    #scale_color_manual(name = "Manufacture", values = manufacture_col) + 
    theme_classic() +
    theme(text = element_text(size = 6)) +
    xlab("") + 
    ylab("|distance to TSS (Mb)|") +
    ggpubr::stat_compare_means(method = "t.test", comparisons = comparisons, label="p.format") +
    # stat_compare_means(comparison=my_comparisons,label="p.format",method="wilcox.test"
    theme(legend.position = "none")

# Beta
unique(beta_df_long_sig_top$type)
#beta_df$type <- gsub("multistate", "multi-state", beta_df$type)
p2 <- ggplot(beta_df_long_sig_top, aes(x=type, y=abs(beta))) +
    geom_violin(scale="width") + 
    geom_boxplot(width=0.1, outlier.shape = NA) +
    #scale_color_manual(name = "Manufacture", values = manufacture_col) + 
    theme_classic() +
    theme(text = element_text(size = 6)) +
    #manuscript_theme +
    xlab("") + 
    ylab("|effect size|") +
    ggpubr::stat_compare_means(method = "t.test", comparisons = comparisons, label="p.format") +
    # stat_compare_means(comparison=my_comparisons,label="p.format",method="wilcox.test"
    theme(legend.position = "none")

p1 | p2

filename <- "/home/hnatri/ILD_eQTL/inteQTL_distTSS_beta.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 5,
    height = 2.5)

p1 | p2

dev.off()

# Comparison with ct-eQTLs
# MashR results
mashr_sighits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-significant-eQTL.rds")
mashr_tophits <- readRDS("/labs/banovich/IPF/eQTL/2022-08-10_38celltypes-mashr/mashr_summary_stats-top-eQTL.rds")
mashr_tophits$rowData$sc_int <- "sc-eQTL"

# Histograms of # of significant cell types
# Counting the # of significant cell types for each int-eQTL
nsig_inteqtls <- as.data.frame(table(inteqtl_topsnp_sig$gene_rsid))
colnames(nsig_inteqtls) <- c("gene_snp", "nSignificant")
nsig_inteqtls$sc_int <- "int-eQTL"

histogram_df <- rbind(mashr_tophits$rowData[,c("sc_int", "nSignificant")],
                      nsig_inteqtls[,c("sc_int", "nSignificant")])

length(unique(inteqtl_topsnp_sig$gene_rsid))

histogram <- ggplot(as.data.frame(histogram_df), aes(x = nSignificant, fill = sc_int)) +
    geom_histogram() + 
    scale_fill_manual(name = "eQTL type", values = c("azure3", "darkblue")) + 
    theme_minimal() +
    ylab("# of top-eQTL") +
    xlab("Significant in # cell types")

histogram

filename <- "/home/hnatri/ILD_eQTL/inteQTL_sceQTL_nSighistogram.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 4,
    height = 2)

histogram

dev.off()

# % of cells expressing each feature in ILD vs CTRL
#p <- DotPlot(ild, features = unique(beta_df_long_sig_top$gene.x), group.by = "Status")
#p_data <- p$data
percent_exp_ild <- lapply(names(object_list_disease), function(xx){
    #percent <- Percent_Expressing(seurat_object = object_list_disease[[xx]], features = unique(beta_df_long_sig_top$gene.x), group_by = "manual_annotation_1")
    #percent$gene_name <- rownames(percent)
#    
    #percent <- pivot_longer(percent, cols = setdiff(colnames(percent), c("gene_name")))
    p <- DotPlot(object_list_disease[[xx]], features = unique(beta_df_long_sig_top$gene.x), group.by = "manual_annotation_1")
    p_data <- p$data
    p_data$status <- "ild"
    
    p_data
    })

percent_exp_ctrl <- lapply(names(object_list_ctrl), function(xx){
    p <- DotPlot(object_list_ctrl[[xx]], features = unique(beta_df_long_sig_top$gene.x), group.by = "manual_annotation_1")
    p_data <- p$data
    p_data$status <- "ctrl"
    
    p_data
    })

# Into a single dataframe
percent_exp_ild_df <- do.call("rbind", percent_exp_ild)
percent_exp_ctrl_df <- do.call("rbind", percent_exp_ctrl)

colnames(percent_exp_ild_df) <- c("ild_avg_exp", "ild_pct_exp", "gene", "celltype", "ild_avg_exp_scaled", "status")
colnames(percent_exp_ctrl_df) <- c("ctrl_avg_exp", "ctrl_pct_exp", "gene", "celltype", "ctrl_avg_exp_scaled", "status")

percent_exp_ild_df$gene_celltype <- paste0(percent_exp_ild_df$gene, "_", percent_exp_ild_df$celltype)
percent_exp_ctrl_df$gene_celltype <- paste0(percent_exp_ctrl_df$gene, "_", percent_exp_ctrl_df$celltype)

setdiff(beta_df_long_sig_top$gene_celltype.x, percent_exp_ild_df$gene_celltype)

# Combining with the int-eQTL df
head(beta_df_long_sig_top)
colnames(beta_df_long_sig_top) <- gsub("gene.x", "gene", colnames(beta_df_long_sig_top))
colnames(beta_df_long_sig_top) <- gsub("gene_celltype.x", "gene_celltype", colnames(beta_df_long_sig_top))
sig_top_df <- merge(beta_df_long_sig_top, percent_exp_ild_df, by = "gene_celltype")
sig_top_df <- merge(sig_top_df, percent_exp_ctrl_df, by = "gene_celltype")
sig_top_df$expdiff <- sig_top_df$ctrl_avg_exp - sig_top_df$ild_avg_exp
sig_top_df$pctdiff <- sig_top_df$ctrl_pct_exp - sig_top_df$ild_pct_exp
head(sig_top_df)

# Adding DEG info
head(degs_df)
sig_top_df <- merge(sig_top_df, degs_df, by = "gene_celltype")

# Plotting some metrics
sig_top_df_small <- sig_top_df[,setdiff(colnames(sig_top_df), NA)]
write.table(sig_top_df_small, "/scratch/hnatri/ILD/sig_top_df.tsv", sep = "\t", quote = F, row.names = F)
sig_top_df_small <- read.table("/scratch/hnatri/ILD/sig_top_df.tsv", header = T, sep = "\t")

hist(sig_top_df_small$p_value, main = "int-eQTL mashr lfsr", xlab = "int-eQTL mashr lfsr")
sig_top_df_small <- sig_top_df_small[which(sig_top_df_small$p_value<0.001),]
sig_top_df_small$exp_in_both <- ifelse(sig_top_df_small$gene_celltype %in% sig_top_df_small[which(sig_top_df_small$ctrl_pct_exp>30 & sig_top_df_small$ild_pct_exp>30),]$gene_celltype, "yes", "no")
sig_top_df_small$exp_in_both <- factor(sig_top_df_small$exp_in_both, levels = c("yes", "no"))

# Which genes had int-eQTLs and were expressed in >30% of cells in both groups
# in that celltype?
select_genes <- sig_top_df_small[which(sig_top_df_small$ctrl_pct_exp>30 & sig_top_df_small$ild_pct_exp>30),]

head(sig_top_df_small)
sig_top_df_small$snp_celltype <- paste0(sig_top_df_small$rsid, "_", sig_top_df_small$celltype.x)
setdiff(maf_res_df_wide$pretty_celltype, sig_top_df_small$celltype.x)
setdiff(sig_top_df_small$celltype.x, maf_res_df_wide$pretty_celltype)
setdiff(sig_top_df_small$snp_celltype, maf_res_df_wide$snp_celltype)
sig_top_df_small <- sig_top_df_small[which(sig_top_df_small$snp_celltype %in% maf_res_df_wide$snp_celltype),]

p1 <- ggplot(sig_top_df_small, aes(x=abs(logFC.x), y = abs(beta))) +
    #geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    xlab("abs(logFC)") +
    theme(legend.position = "none")

p2 <- ggplot(sig_top_df_small, aes(x=abs(log(expdiff+1)), y = abs(beta))) +
    #geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    theme(legend.position = "none")

p3 <- ggplot(sig_top_df_small, aes(x=abs(pctdiff), y = abs(beta))) +
    #geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    theme(legend.position = "none")

p4 <- ggplot(sig_top_df_small, aes(x=log(ild_avg_exp+1), y = abs(beta))) +
    #geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    theme(legend.position = "none")

p5 <- ggplot(sig_top_df_small, aes(x=log(ctrl_avg_exp+1), y = abs(beta))) +
    #geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    theme(legend.position = "none")

p6 <- ggplot(sig_top_df_small, aes(x=ild_pct_exp, y = abs(beta))) +
    #geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    theme(legend.position = "none")

p7 <- ggplot(sig_top_df_small, aes(x=ctrl_pct_exp, y = abs(beta))) +
    #$geom_pointdensity() +
    #scale_colour_gradientn(colors = terrain.colors(n=8), limits=c(1, 40000)) +
    geom_point(aes(colour = exp_in_both)) +
    scale_color_manual(values = c("brown2", "azure3")) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    theme(legend.position = "none")

p1 <- patchwork::wrap_elements(ggExtra::ggMarginal(p1, groupColour = T, groupFill = T))
p2 <- patchwork::wrap_elements(ggExtra::ggMarginal(p2, groupColour = T, groupFill = T))
p3 <- patchwork::wrap_elements(ggExtra::ggMarginal(p3, groupColour = T, groupFill = T))
p4 <- patchwork::wrap_elements(ggExtra::ggMarginal(p4, groupColour = T, groupFill = T))
p5 <- patchwork::wrap_elements(ggExtra::ggMarginal(p5, groupColour = T, groupFill = T))
p6 <- patchwork::wrap_elements(ggExtra::ggMarginal(p6, groupColour = T, groupFill = T))
p7 <- patchwork::wrap_elements(ggExtra::ggMarginal(p7, groupColour = T, groupFill = T))

patched <- (p1 | p2 | p3 ) / (p4 | p5 | p6) /  (p7 | plot_spacer() | plot_spacer())
patched #+ plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ======================================
# Heatmap of # int-eQTL and DEGs
# ======================================

# Matrix of the # of int-eGenes for each celltype
#sig_res_mx <- table(beta_df_long_sig_top$gene.x, beta_df_long_sig_top$celltype.x)
beta_df_long_sig_top_n <- beta_df_long_sig_top %>% group_by(celltype.x) %>% summarise(n_egenes = n_distinct(gene))
#sig_res_mx <- beta_df_long_sig_top_n
colnames(beta_df_long_sig_top_n) <- gsub("celltype.x", "celltype", colnames(beta_df_long_sig_top_n))

# Adding # of DEGs
head(degs_sig_df)
setdiff(beta_df_long_sig_top$celltype.x, paper_celltype_colors$`Pretty cell type name`)
setdiff(degs_sig_df$celltype, paper_celltype_colors$`Pretty cell type name`)
length(unique(beta_df_long_sig_top_n$celltype))
length(unique(degs_sig_df$celltype))
#degs_sig_df$celltype <- gsub("Differentiating Ciliated", "Differentiating ciliated", degs_sig_df$celltype)
setdiff(beta_df_long_sig_top_n$celltype, degs_sig_df$celltype)
setdiff(degs_sig_df$celltype, beta_df_long_sig_top_n$celltype)
degs_sig_df_n <- degs_sig_df %>% group_by(celltype) %>% summarise(n_degs = n_distinct(feature))
sig_res_mx <- merge(beta_df_long_sig_top_n, degs_sig_df_n, by = "celltype")

setdiff(sig_res_mx$celltype, paper_celltype_colors$`Pretty cell type name`) 
setdiff(paper_celltype_colors$`Pretty cell type name`, sig_res_mx$celltype) 

rownames(sig_res_mx) <- sig_res_mx$celltype
dim(sig_res_mx)

# Adding pairs with 0 significant genes
setdiff(paper_celltype_colors$`Pretty cell type name`, rownames(sig_res_mx))
setdiff(rownames(sig_res_mx), paper_celltype_colors$`Pretty cell type name`)

sig_res_mx_add <- matrix(0, nrow=length(setdiff(paper_celltype_colors$`Pretty cell type name`, rownames(sig_res_mx))), ncol=ncol(sig_res_mx))
rownames(sig_res_mx_add) <- setdiff(paper_celltype_colors$`Pretty cell type name`, rownames(sig_res_mx))
colnames(sig_res_mx_add) <- colnames(sig_res_mx)

#sig_res_mx_wzeros <- rbind(sig_res_mx, sig_res_mx_add)
#sig_res_mx_wzeros <- rbind(sig_res_mx_wzeros, rep(NA, ncol(sig_res_mx_wzeros)))

# Adding summary rows for the four cell lineages and all celltypes
#sig_res_mx_wzeros[c("All cell types"),] <- colSums(sig_res_mx_wzeros[setdiff(rownames(sig_res_mx_wzeros), "All cell types"),])
colnames(sig_res_mx)

# Adding numbers of eGenes that overlap with DEGs
sig_res_mx$n_overlap <- NA
sig_res_mx$prop_overlap <- NA
head(sig_res_mx)
for(celltype in sig_res_mx$celltype){
    egenes <- beta_df_long_sig_top[which(beta_df_long_sig_top$celltype.x==celltype),]$gene 
    degs <- degs_sig_df[which(degs_sig_df$celltype==celltype),]$feature
    
    n_overlapping <- length(intersect(egenes, degs))
    prop_overlapping <- n_overlapping/length(unique(egenes, degs))
    
    sig_res_mx[celltype, "n_overlap"] <- n_overlapping
    sig_res_mx[celltype, "prop_overlap"] <- prop_overlapping
}

setdiff(rownames(sig_res_mx), c("AT1",
          "Transitional AT2",
          "AT2",
          "Proliferating - Epi",
          "Secretory - SCGB1A1+/MUC5B+",
          "Secretory - SCGB1A1+/SCGB3A2+",
          "Secretory - SCGB3A2+",
          "KRT5-/KRT17+",
          "Ciliated",
          "Differentiating ciliated",
          "Basal",
          "NK",
          "CD4",
          "CD8/NKT",
          "Plasma",
          "B cells",
          "cDC1",
          "moDC",
          "Macrophage - SPP1+",
          "Alveolar macrophage",
          "Monocyte",
          "Inflammatory monocyte",
          "Monocyte-derived macrophage",
          "Mast",
          "pDC",
          "cDC2",
          "Proliferating - Imm",
          "Lymphatic",
          "Arteriole",
          "Peribronchiolar",
          "Venule",
          "CA4+ capillary",
          "Capillary",
          "Matrix FB",
          "SMC",
          "WNT2+ FB",
          "Mesothelial",
          "Pericyte"))

setdiff(names(celltype_cols), rownames(sig_res_mx))
setdiff(rownames(sig_res_mx), names(celltype_cols))

# Reordering rows
#rownames(sig_res_mx) <- gsub("Lung", "GTEx Lung", rownames(sig_res_mx))
#rownames(sig_res_mx) <- gsub("Blood", "GTEx Whole Blood", rownames(sig_res_mx))
#rownames(sig_res_mx) <- gsub("Brain", "GTEx Brain Cortex", rownames(sig_res_mx))
setdiff(c("AT1",
          "Transitional AT2",
          "AT2",
          "Proliferating - Epi",
          "Secretory - SCGB1A1+/MUC5B+",
          "Secretory - SCGB1A1+/SCGB3A2+",
          "Secretory - SCGB3A2+",
          #"KRT5-/KRT17+",
          "Ciliated",
          "Differentiating ciliated",
          "Basal",
          "NK",
          "CD4",
          "CD8/NKT",
          "Plasma",
          "B cells",
          "cDC1",
          "moDC",
          #"Macrophage - SPP1+",
          "Alveolar macrophage",
          "Monocyte",
          "Inflammatory monocyte",
          "Monocyte-derived macrophage",
          "Mast",
          #"pDC",
          #"cDC2",
          "Proliferating - Imm",
          "Lymphatic",
          "Arteriole",
          #"Peribronchiolar",
          "Venule",
          #"CA4+ capillary",
          #"Capillary",
          #"Matrix FB",
          "SMC",
          #"WNT2+ FB",
          #"Mesothelial",
          "Pericyte"), rownames(sig_res_mx))
sig_res_mx <- sig_res_mx[c("AT1",
                           "Transitional AT2",
                           "AT2",
                           "Proliferating - Epi",
                           "Secretory - SCGB1A1+/MUC5B+",
                           "Secretory - SCGB1A1+/SCGB3A2+",
                           "Secretory - SCGB3A2+",
                           #"KRT5-/KRT17+",
                           "Ciliated",
                           "Differentiating ciliated",
                           "Basal",
                           "NK",
                           "CD4",
                           "CD8/NKT",
                           "Plasma",
                           "B cells",
                           "cDC1",
                           "moDC",
                           #"Macrophage - SPP1+",
                           "Alveolar macrophage",
                           "Monocyte",
                           "Inflammatory monocyte",
                           "Monocyte-derived macrophage",
                           "Mast",
                           #"pDC",
                           #"cDC2",
                           "Proliferating - Imm",
                           "Lymphatic",
                           "Arteriole",
                           #"Peribronchiolar",
                           "Venule",
                           #"CA4+ capillary",
                           #"Capillary",
                           #"Matrix FB",
                           "SMC",
                           #"WNT2+ FB",
                           #"Mesothelial",
                           "Pericyte"),]

# Row annotations
lineage <- data.frame("lineage" = c(rep("Epithelial", 10),
                                    rep("Immune", 13),
                                    rep("Endothelial", 3),
                                    rep("Mesenchymal", 2)))
                                    #rep("GTEx", 3)))
lineage_cols <- c(lineage_cols, "GTEx" = "azure3", "All cell types" = "white")

row_annot <- rowAnnotation(df = lineage,
                           col = list("lineage" = lineage_cols),
                           width = unit(0.8, "cm"))

# Constructing the heatmap
sig_res_mx$perc_overlap <- sig_res_mx$prop_overlap*100
sig_res_mx_egenes <- as.matrix(sig_res_mx[,2])
sig_res_mx_degs <- as.matrix(sig_res_mx[,3])
# prop_overlap
sig_res_mx_overlap <- as.matrix(sig_res_mx[,6])
rownames(sig_res_mx_egenes) <- rownames(sig_res_mx)
rownames(sig_res_mx_degs) <- rownames(sig_res_mx)
rownames(sig_res_mx_overlap) <- rownames(sig_res_mx)
cn <- c("perc_overlap")

max(sig_res_mx$n_egenes)
max(sig_res_mx$n_degs)
max(sig_res_mx$perc_overlap)

heatmap <- Heatmap(sig_res_mx_overlap,
                   name = "int-eGene-DEG\noverlap (%)",
                   #name = "Significant\nDEGs",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_column_names = FALSE,
                   row_names_gp = grid::gpar(fontsize = 7),
                   #column_names_gp = grid::gpar(fontsize = 7),
                   heatmap_legend_param = list(legend_gp = gpar(fontsize = 7),
                                               title_gp = gpar(fontsize = 7)),
                   bottom_annotation = HeatmapAnnotation(
                       text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right", gp = gpar(fontsize = 6)),
                       annotation_height = max_text_width(cn)),
                   row_names_side = "left",
                   col = colorRamp2(c(0, 100), c("white", "chartreuse4")), # chartreuse4
                   width = 1*unit(4, "mm"), 
                   height = nrow(sig_res_mx_egenes)*unit(4, "mm"),
                   #cell_fun = function(j, i, x, y, width, height, fill){
                   #    grid.text(sprintf("%.0f", sig_res_mx_egenes[i, j]), x, y, gp = gpar(fontsize = 6))
                   #},
                   #top_annotation = col_annot,
                   left_annotation = row_annot)

heatmap

draw(heatmap)

filename <- "/home/hnatri/ILD_eQTL/heatmap_integene_deg_percent_overlap.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 4,
    height = 7)

heatmap

dev.off()

cn <- c("n_egenes")
heatmap <- Heatmap(sig_res_mx_egenes,
                   name = "int-eGenes",
                   #name = "Significant\nDEGs",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_column_names = FALSE,
                   row_names_gp = grid::gpar(fontsize = 7),
                   #column_names_gp = grid::gpar(fontsize = 7),
                   heatmap_legend_param = list(legend_gp = gpar(fontsize = 7),
                                               title_gp = gpar(fontsize = 7)),
                   bottom_annotation = HeatmapAnnotation(
                       text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right", gp = gpar(fontsize = 6)),
                       annotation_height = max_text_width(cn)),
                   row_names_side = "left",
                   col = colorRamp2(c(0, 6027), c("white", "brown3")), # chartreuse4
                   width = 1*unit(4, "mm"), 
                   height = nrow(sig_res_mx_egenes)*unit(4, "mm"),
                   #cell_fun = function(j, i, x, y, width, height, fill){
                   #    grid.text(sprintf("%.0f", sig_res_mx_egenes[i, j]), x, y, gp = gpar(fontsize = 6))
                   #},
                   #top_annotation = col_annot,
                   left_annotation = row_annot)

heatmap

draw(heatmap)

filename <- "/home/hnatri/ILD_eQTL/heatmap_integene.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 4,
    height = 7)

heatmap

dev.off()

cn <- c("n_egenes")
heatmap <- Heatmap(sig_res_mx_degs,
                   name = "DEGs",
                   #name = "Significant\nDEGs",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_column_names = FALSE,
                   row_names_gp = grid::gpar(fontsize = 7),
                   #column_names_gp = grid::gpar(fontsize = 7),
                   heatmap_legend_param = list(legend_gp = gpar(fontsize = 7),
                                               title_gp = gpar(fontsize = 7)),
                   bottom_annotation = HeatmapAnnotation(
                       text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right", gp = gpar(fontsize = 6)),
                       annotation_height = max_text_width(cn)),
                   row_names_side = "left",
                   col = colorRamp2(c(0, 12896), c("white", "darkblue")), # chartreuse4
                   width = 1*unit(4, "mm"), 
                   height = nrow(sig_res_mx_egenes)*unit(4, "mm"),
                   #cell_fun = function(j, i, x, y, width, height, fill){
                   #    grid.text(sprintf("%.0f", sig_res_mx_egenes[i, j]), x, y, gp = gpar(fontsize = 6))
                   #},
                   #top_annotation = col_annot,
                   left_annotation = row_annot)

heatmap

draw(heatmap)

filename <- "/home/hnatri/ILD_eQTL/heatmap_degs.pdf"
pdf(file = filename,
    #units="in",
    #res = 100,
    width = 4,
    height = 7)

heatmap

dev.off()
