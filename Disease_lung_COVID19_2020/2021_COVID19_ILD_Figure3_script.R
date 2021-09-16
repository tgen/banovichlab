# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-04-12
# Title: scRNA-seq data analysis for Covid-19 manuscript - Figure 3
# ====================================

# ==============================================================================
# SET UP THE ENVIRONMENT VARIABLES 
# ==============================================================================
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()
options(future.globals.maxSize = 4096*1024^2 )
set.seed(2811)

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RCurl)
library(gplots)
library(BSDA)
library(reshape2)
library(UpSetR)
library(calibrate)
library(dplyr)
library(ggrepel)
library(broom)
library(rstatix)
library(sctransform)
library(topGO)
library(org.Hs.eg.db)

# ==============================================================================
# Read in the Epithelial object (containing all cell types)
# ==============================================================================
epi <- readRDS("/scratch/lbui/Covid19_ILD_objects/20210204_Epithelial_noDoublets.rds")

# ==============================================================================
# Figure 3A: Gene expression correlation
# ==============================================================================
# Included in the scripts for Fig S10 (2021_COVID19_ILD_Supplementary_figures.R)

# ==============================================================================
# Figure 3B Boxplot for some candidates AT2 cells
# ==============================================================================
epi_list = list()
j <- 0
for(i in unique(epi@meta.data$CellType2)){
  j <- j + 1
  epi_list[[j]] <- subset(epi,
                          cells = row.names(epi@meta.data[epi@meta.data$CellType2 == i, ]))
}
for(i in 1:length(epi_list)){
  names(epi_list) <- lapply(epi_list, function(xx){paste(unique(xx@meta.data$CellType2))})
}
summary(epi_list)

genelist2 <- c("ACE2","TMPRSS2","NPR1","CTSL","FURIN","BSG","ADAM17","HSPA5","ITGB6","NT5DC1",
               "SOCS1","SOCS2","CSF3","FAM46C","UBD","REC8","ELF1","CLEC4D","LY6E")

# Calculate p_adj_value using FindMarkers neg binom test
copd_vs_control <- lapply(epi_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "COPD",
                ident.2 = "Control",
                test.use = "negbinom", # using this test: the slot is set to "counts"
                features = genelist2,
                latent.vars = c("dataset","Age","Ethinicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

ipf_vs_control <- lapply(epi_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "IPF",
                ident.2 = "Control",
                test.use = "negbinom",
                features = genelist2,
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

other_vs_control <- lapply(epi_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "Other-ILD",
                ident.2 = "Control",
                test.use = "negbinom", 
                features = genelist2,
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

# Save the DEG files
for(i in 1:length(epi_list)){
  write.table(copd_vs_control[[i]], paste(gsub("/", "", unique(epi_list[[i]]@meta.data$CellType2)), "_copd_vs_control", ".csv"), sep =",", quote = F)
}
for(i in 1:length(epi_list)){
  write.table(ipf_vs_control[[i]], paste(gsub("/", "", unique(epi_list[[i]]@meta.data$CellType2)), "_ipf_vs_control", ".csv"), sep =",", quote = F)
}
for(i in 1:length(epi_list)){
  write.table(other_vs_control[[i]], paste(gsub("/", "", unique(epi_list[[i]]@meta.data$CellType2)), "_other_vs_control", ".csv"), sep =",", quote = F)
}

# Extract out only the p_adj_values from the DEG lists
p_val1=list()
for(i in 1:length(copd_vs_control)){
  colnames(copd_vs_control[[i]]) <- c("p_val","avg_logFC","pct.1","pct.2",names(copd_vs_control[i]))
  p_val1[[i]] <- copd_vs_control[[i]][,4:5]
  p_val1[[i]]$GeneID <- rownames(p_val1[[i]])
}
copd_pval <- Reduce(
  function(x, y) merge(x, y, by = "GeneID", all = T),
  lapply(p_val1, function(x) { x$id <- rownames(x); x }))
copd_pval$pct.2.x <- NULL
copd_pval$pct.2.y <- NULL
copd_pval$pct.2 <- NULL
copd_pval$id.x <- NULL
copd_pval$id <- NULL
copd_pval$id.y <- NULL
rownames(copd_pval) <- copd_pval[,1]

p_val2=list()
for(i in 1:length(ipf_vs_control)){
  colnames(ipf_vs_control[[i]]) <- c("p_val","avg_logFC","pct.1","pct.2",names(ipf_vs_control[i]))
  p_val2[[i]] <- ipf_vs_control[[i]][,4:5]
  p_val2[[i]]$GeneID <- rownames(p_val2[[i]])
}
ipf_pval <- Reduce(
  function(x, y) merge(x, y, by = "GeneID", all = T),
  lapply(p_val2, function(x) { x$id <- rownames(x); x }))
ipf_pval$pct.2.x <- NULL
ipf_pval$pct.2.y <- NULL
ipf_pval$pct.2 <- NULL
ipf_pval$id.x <- NULL
ipf_pval$id <- NULL
ipf_pval$id.y <- NULL
rownames(ipf_pval) <- ipf_pval[,1]

p_val3=list()
for(i in 1:length(other_vs_control)){
  colnames(other_vs_control[[i]]) <- c("p_val","avg_logFC","pct.1","pct.2",names(other_vs_control[i]))
  p_val3[[i]] <- other_vs_control[[i]][,4:5]
  p_val3[[i]]$GeneID <- rownames(p_val3[[i]])
}
other_pval <- Reduce(
  function(x, y) merge(x, y, by = "GeneID", all = T),
  lapply(p_val3, function(x) { x$id <- rownames(x); x }))
other_pval$pct.2.x <- NULL
other_pval$pct.2.y <- NULL
other_pval$pct.2 <- NULL
other_pval$id.x <- NULL
other_pval$id <- NULL
other_pval$id.y <- NULL
rownames(other_pval) <- other_pval[,1]

copd_pval$group1 <- "Control"
ipf_pval$group1 <- "Control"
other_pval$group1 <- "Control"
copd_pval$group2 <- "COPD"
ipf_pval$group2 <- "IPF"
other_pval$group2 <- "Other-ILD"

# Combine all files to generate a file with all p_adj_values 
pval <- rbind(melt(copd_pval), melt(ipf_pval), melt(other_pval))
colnames(pval) <- c("GeneID","group1","group2","CellType","p_adj")
write.csv(pval, file = "20210207_Epi_Boxplot_p_adj_val.csv")

# Add significant into the p_val file
onion <- pval$p_adj
onion[onion <= 0.05] <- "**"
onion[onion > 0.05 & onion <= 0.1] <- "*"
onion[onion > 0.1] <- NA
pval$significant <- onion 

# Extract out gene counts
assayData <- GetAssayData(epi, slot = "counts", assay = "SCT")
gene1 <- data.frame(assayData[rownames(assayData) == "ITGB6",])
gene2 <- data.frame(assayData[rownames(assayData) == "SOCS2",])
gene3 <- data.frame(assayData[rownames(assayData) == "LY6E",])
gene4 <- data.frame(assayData[rownames(assayData) == "CSF3",])

## Adding metadata
colnames(gene1) <- c("counts")
gene1$CellType <- epi@meta.data$CellType2
gene1$Diagnosis <- epi@meta.data$Diagnosis2
gene1$GeneID <- "ITGB6"

colnames(gene2) <- c("counts")
gene2$CellType <- epi@meta.data$CellType2
gene2$Diagnosis <- epi@meta.data$Diagnosis2
gene2$GeneID <- "SOCS2"

colnames(gene3) <- c("counts")
gene3$CellType <- epi@meta.data$CellType2
gene3$Diagnosis <- epi@meta.data$Diagnosis2
gene3$GeneID <- "LY6E"

colnames(gene4) <- c("counts")
gene4$CellType <- epi@meta.data$CellType2
gene4$Diagnosis <- epi@meta.data$Diagnosis2
gene4$GeneID <- "CSF3"

# Combine all gene count values for all genes
plot_data <- rbind(gene1,gene2,gene3,gene4)

# Only plot AT2 cells
plot_data.df <- plot_data[plot_data$CellType == "AT2",]
plot_data.df <- plot_data.df[plot_data.df$counts > 0,]

# Extract out p_adj_val
stat.test <- pval[pval$GeneID %in% c("ITGB6","SOCS2","LY6E","CSF3"),]
stat.test <- stat.test[stat.test$CellType == "AT2",]
stat.test <- stat.test[order(stat.test$GeneID),]
stat.test$y.position <- c(40,45,50,22,24,26,30,32,34,16,18,20)

# Plotting
## Boxplot
plot_data.df2 <- plot_data.df[plot_data.df$GeneID == "SOCS2",]
stat.test2 <- stat.test[stat.test$GeneID == "SOCS2",]
ggplot(plot_data.df2, aes(x=Diagnosis, y=counts, color = Diagnosis)) + 
  geom_boxplot(outlier.size = -1,width = 0.75, na.rm = T) +
  geom_jitter(position = position_jitter(0.15)) + 
  theme_bw() +
  stat_pvalue_manual(stat.test2, size=4, tip.length = 0.01, step.group.by = "GeneID",
                     label = "significant") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=0)) +
  ylim(0,20) + 
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=10))

# ==============================================================================
# Figure 3C: DEG and upset plot
# ==============================================================================
# Subset out AT2 cells
at2 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$CellTypeSimple == "AT2",]))

# Subset out ACE2+ AT2 cells
ace2 <- subset(at2,  ACE2 > 0, slot = "counts")
# Count number of ACE2+ cells
table(ace2@meta.data$Diagnosis2)
#Control      COPD       IPF Other-ILD 
#427        73       138       141 

# Subset out control vs. cld cells
at2.control <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$Diagnosis2 == "Control",]))
at2.cld <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$Diagnosis2 != "Control",]))

# Get the cell BC for ACE2+ cells
cellBC <- rownames(ace2@meta.data)

# Add a column into the meta data for ACE2+ cells
onion <- rownames(at2.control@meta.data)
output <- ifelse(onion %in% cellBC, "ACE2+", "ACE2-")
at2.control@meta.data$ACE2 <- output

onion <- rownames(at2.cld@meta.data)
output <- ifelse(onion %in% cellBC, "ACE2+", "ACE2-")
at2.cld@meta.data$ACE2 <- output

# Run a DEG analysis
at2.control_deg <- FindMarkers(at2.control,
                               group.by = "ACE2",
                               ident.1 = "ACE2+",
                               ident.2 = "ACE2-",
                               test.use = "negbinom",
                               latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                               assay ="SCT",
                               logfc.threshold = 0)
write.csv(at2.control_deg, file = "20210220_AT2_control_ACE2_deg.csv")

at2.cld_deg <- FindMarkers(at2.cld,
                           group.by = "ACE2",
                           ident.1 = "ACE2+",
                           ident.2 = "ACE2-",
                           test.use = "negbinom",
                           latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                           assay ="SCT",
                           logfc.threshold = 0)
write.csv(at2.cld_deg, file = "20210220_AT2_cld_ACE2_deg.csv")

# DEG for ACE2+ in Control vs. CLD samples
ace2_CLD_con <- FindMarkers(ace2,
                            group.by = "Status",
                            ident.1 = "Disease",
                            ident.2 = "Control",
                            test.use = "negbinom",
                            latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                            assay ="SCT",
                            logfc.threshold = 0)
write.csv(ace2_CLD_con, file = "20210222_AT2_ACE2_CLD_Control_deg.csv")

# DEG for AT2 in Control vs. CLD samples
at2_CLD_con <- FindMarkers(at2,
                           group.by = "Status",
                           ident.1 = "Disease",
                           ident.2 = "Control",
                           test.use = "negbinom",
                           latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                           assay ="SCT",
                           logfc.threshold = 0)
write.csv(ace2_CLD_con, file = "20210220_AT2_CLD_Control_deg.csv")

# Create a list
deg_list <- list(at2.control_deg, at2.cld_deg, ace2_CLD_con, at2_CLD_con)
names(deg_list) <- c("ACE2+/-_Control","ACE2+/-_CLD","ACE2_CLD_Control", "AT2_CLD_Control")

# Make an Upset plot for shared of those 3 comparisons
onion <- lapply(deg_list, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(deg_list, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(deg_list)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

upset(upset_d_vs_c, nsets = 4, text.scale = 2, show.numbers = T)

# ==============================================================================
# Figure 3E: ACE2+ correlated gene analysis
# ==============================================================================
# Subset out AT2 ACE2+ cells
at2 <- subset(epi, cells=rownames(epi@meta.data[epi@meta.data$CellTypeSimple == "AT2",]))
at2.df <- subset(at2, ACE2 > 0)

ace2p_list = list()
j <- 0
for(i in unique(at2.df@meta.data$Diagnosis2)){
  j <- j + 1
  ace2p_list[[j]] <- subset(at2.df,
                            cells = row.names(at2.df@meta.data[at2.df@meta.data$Diagnosis2 == i, ]))
}
for(i in 1:length(ace2p_list)){
  names(ace2p_list) <- lapply(ace2p_list, function(xx){paste(unique(xx@meta.data$Diagnosis2))})
}
summary(ace2p_list)

# Perform Spearman correlation analysis
cor_list1 = list()
matrix.1 = list()
matrix_mod = list()
gene = list()
j <- 0
for(i in 1:length(ace2p_list)){
  j = j+1
  matrix.1[[j]] <- ace2p_list[[i]]@assays$SCT@data
  matrix_mod[[j]] <- as.matrix(matrix.1[[j]])
  gene[[j]] <- as.numeric(matrix_mod[[j]]["ACE2",])
  cor_list1[[j]] <- apply(matrix_mod[[j]],1,
                          function(x){cor.test(gene[[j]],x,method = "spearman",
                                               adjust = "BH", exact = F,
                                               na.action = "na.exclude")})
}

onion <- lapply(cor_list1, function(x) as.data.frame(sapply(x, function(xx) cbind(xx$p.value,xx$estimate))))
onion <- lapply(onion, function(x) as.data.frame(t(x)))
onion <- lapply(onion, function(x) x[!is.na(x[,2]),])
colnames <-  c("p_value","rho")
onion <- lapply(onion, setNames, colnames)

for(i in 1:length(onion)){
  onion[[i]]$geneID <- rownames(onion[[i]])
  onion[[i]]$Set <- names(ace2p_list[i])
}

# Prepare data for plotting
plot_data <- do.call(rbind, onion)
#plot_data <- plot_data[!duplicated(plot_data$geneID),]
plot_data <- plot_data[plot_data$geneID != "ACE2",]
temp1 <- grep( "^MT-", rownames(plot_data), ignore.case = F, value = T)
temp2 <- grep( "^RP", rownames(plot_data), ignore.case = F, value = T)
plot_data <- plot_data[!plot_data$geneID %in% temp1,]
plot_data <- plot_data[!plot_data$geneID %in% temp2,]

copd.gene <- c("OAS1","CLEC3B","WNT11","IGSF6","ACKR4","MST1R","LST1",
               "OASL","IFI6","IFIT1","IFIT2","IGFBP2","FOXA3","TNFAIP8L2",
               "IGLV3-10","SOX9")
ipf.gene <- c("NXF3","RMI1","ITPKB","IGFBP7","CLDN10","CES3","CD93",
              "IGSF22","FKBP10","MCM4","CDH1","SYT16","SP4","ITPKB")
control.gene <- c("SCN7A","AC245052.4")
other.gene <- c("IL10ORB-DT","ITLN1","NOD2","IGFBP2","PRKX","MUC4","KRT13",
                "CLDN4","MUC16","ITGB8","TNFRSF11B","NOD2")

plot_data$genelabels <- ""
plot_data$genelabels <- ifelse(plot_data$geneID %in% copd.gene & plot_data$Set == "COPD" 
                               |plot_data$geneID %in% ipf.gene & plot_data$Set == "IPF"
                               |plot_data$geneID %in% other.gene & plot_data$Set == "Other-ILD"
                               |plot_data$geneID %in% control.gene & plot_data$Set == "Control"
                               , TRUE,FALSE)

# Plot correlation 
pdf("20210207_ACE2_correlation_SCT_data_clean.pdf")
ggplot(plot_data, aes(x=-log10(p_value), y=rho, color = Set)) +
  geom_point() +
  xlab(expression("-Log"[10]*"(p_value)")) +
  ylab(expression("Spearman rho values")) +
  theme_bw() +
  geom_hline(aes(yintercept = quantile(rho, 0.01, na.rm = TRUE)), 
             size = 0.5,color="black", linetype="dashed") +
  geom_hline(aes(yintercept = quantile(rho, 0.99, na.rm = TRUE)), 
             size=0.5,color="black", linetype="dashed") +
  geom_text_repel(label = ifelse(plot_data$genelabels == TRUE, 
                                 plot_data$geneID,""),
                  max.overlaps=1000) +
  #geom_text_repel(data=subset(plot_data,abs(rho) >= 0.5 & -log10(p_value) >= 3),
  #                            aes(-log10(p_value), rho, label = geneID),
  #                                 size = 3,max.overlaps=500) +
  ggtitle("ACE2 correlated genes")
dev.off()

write.csv(plot_data, file = "20210207_ACE2_correlation_SCT_data.csv")

# GO analysis with TopGO
# Prepare the TopGO genelist for all genes in the correlation analysis
cor_list <- plot_data
cor_genes <- as.vector(cor_list[,1])
names(cor_genes) <- rownames(cor_list)

# Filter out significant ACE2 correlated genes (p_value <= 0.03, 99th quantile rho)
quantile(cor_list$rho,0.99)
genelist <- cor_list[cor_list$rho >= quantile(cor_list$rho,0.99),]

# Count number of significant genes
dim(genelist[genelist$Set == "COPD",])
# 330   4
dim(genelist[genelist$Set == "Control",])
# 2 4
dim(genelist[genelist$Set == "IPF",])
# 108   4
dim(genelist[genelist$Set == "Other-ILD",])
# 268   4

# Extract out the significant genelist in the CLD samples
cld_genes <- as.vector(genelist[genelist$Diagnosis != "Control",][,2])
names(cld_genes) <- rownames(genelist[genelist$Diagnosis != "Control",])

# Run TopGO for the significant genes among all the correlated genes in the result table
GOdata <- new('topGOdata', 
              ontology = "BP", 
              allGenes = cor_genes, 
              geneSelectionFun = function(x){return(x %in% cld_genes)},
              annot=annFUN.org, 
              mapping = 'org.Hs.eg.db', 
              ID = 'symbol')

# Run statistic tests
allGO = usedGO(object = GOdata) 
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher
tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(allGO), numChar = 120)
tab

#performing BH correction on our p values
p.adj <- round(p.adjust(tab$raw.p.value,method="BH"),digits = 4)
tab <- cbind(tab, p.adj)

# Save results
write.csv(tab, file = "20210208_ACE2_CLD_correlated_TopGO_Fisher.csv")
