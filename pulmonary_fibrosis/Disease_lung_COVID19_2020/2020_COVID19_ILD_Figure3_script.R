# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-07-03
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
library(reshape)
library(UpSetR)
library(calibrate)
library(dplyr)
library(ggrepel)
library(broom)
library(rstatix)

# ==============================================================================
# Read in the Immune object 
# ==============================================================================
immune <- readRDS("/scratch/lbui/Covid19_Seurat/20200709_Immune_noDoublets.rds")

# ==============================================================================
# Figure 3A: Immune cell proportion
# ==============================================================================
# Subset to 4 objects
control <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$Diagnosis2 == "Control",]))
copd <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$Diagnosis2 == "COPD",]))
ipf <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$Diagnosis2 == "IPF",]))
other <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$Diagnosis2 == "Other-ILD",]))

# Extract out total number of cells per sample per Diagnosis
onion1 <- as.data.frame(table(control$Sample_Name,control$CellType2))
onion2 <- as.data.frame(table(copd$Sample_Name,copd$CellType2))
onion3 <- as.data.frame(table(ipf$Sample_Name,ipf$CellType2))
onion4 <- as.data.frame(table(other$Sample_Name,other$CellType2))

onion1$Diagnosis <- "Control"
onion2$Diagnosis <- "COPD"
onion3$Diagnosis <- "IPF"
onion4$Diagnosis <- "Other-ILD"

# Name the columns
colnames(onion1) <- c("Sample","CellType","Counts","Diagnosis")
colnames(onion2) <- c("Sample","CellType","Counts","Diagnosis")
colnames(onion3) <- c("Sample","CellType","Counts","Diagnosis")
colnames(onion4) <- c("Sample","CellType","Counts","Diagnosis")

# Calculate percentage of each CT in all Immune cells
onion1$Percent <- onion1$Counts/length(rownames(immune@meta.data)) * 100
onion2$Percent <- onion2$Counts/length(rownames(immune@meta.data)) * 100
onion3$Percent <- onion3$Counts/length(rownames(immune@meta.data)) * 100
onion4$Percent <- onion4$Counts/length(rownames(immune@meta.data)) * 100

# Combine all data into 1 data frame
data_table <- rbind(onion1, onion2, onion3, onion4)

# Plotting
data_table$Diagnosis <- factor(data_table$Diagnosis,
                               levels = c("Other-ILD","IPF","COPD","Control"))
ggplot(data_table, aes(x=CellType, y=Percent,fill = Diagnosis)) + 
  coord_flip() +
  geom_boxplot(outlier.size = -1) + 
  ylim(0,1) +
  theme_bw() +
  ylab("Proportion Immune cells")

# Tukey post-hoc test for multitesting
data_table$Diagnosis <- as.factor(data_table$Diagnosis)
levels(data_table$Diagnosis)

res.aov <- data_table %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis, p.adjust.method = "bonferroni") 

write.csv(as.data.frame(res.aov), file = "20200730_Immune_cell_proportion_Tukey.csv")

# ==============================================================================
# Figure 3B: Immune cell heatmap
# ==============================================================================
immune_list = list()
j <- 0
for(i in unique(immune@meta.data$CellType2)){
  j <- j + 1
  immune_list[[j]] <- subset(immune,
                             cells = row.names(immune@meta.data[immune@meta.data$CellType2== i, ]))
}
for(i in 1:length(immune_list)){
  names(immune_list) <- lapply(immune_list, function(xx){paste(unique(xx@meta.data$CellType2))})
}
summary(immune_list)

immune_list <- immune_list[names(immune_list) != "TRegs"]

# Run DEG for Covid-19 candidate genes    
genelist <- c("ACE2","TMPRSS2","CTSL","CTSB","FURIN","PCSK5","PCSK7","ADAM17",
              "BSG","PIKFYVE","TPCN2","ACE","ISG15", "IFI44", "IFI27", "CXCL10", 
              "RSAD2", "IFIT1", "IFI44L", "CCL8", "XAF1", "GBP1", "IRF7", 
              "CEACAM1","IFNB1","IFNG","IFNAR1","IFNAR2","IFNGR1","IFNGR2","TRIM28",
              "TLR7","MX1","STAT1","TBK1","CCR2","CXCL10","IFI6","LY6E",
              "TRIM27","TNF","IL1B","IL6","IL6R","IL6ST","TGFB1","NFKB1",
              "NFKB2","CEBPB","AREG","FCGR3A","IL10","IFITM1","IFITM3","ISG20",
              "CD163","IL1R2","MRC1","HAVCR2","LGALS9","S100A8","S100A9",
              "HLA-DRA", "HLA-DQA1", "HLA-DQA2", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", 
              "HLA-DQB2", "HLA-DRB5", "HLA-DQB1", "HLA-DMA", "HLA-DMB")

heatmap_deg <- lapply(immune_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  FindMarkers(xx,
              group.by = "Status",
              ident.1 = "Disease",
              ident.2 = "Control",
              features = unique(genelist),
              test.use = "negbinom",
              latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
              logfc.threshold = 0)
})
saveRDS(heatmap_deg, file = "20200721_Immune_disease_control_deglist.rds")

# Save the DEG files
for(i in 1:length(immune_list)){
  write.table(heatmap_deg[[i]], paste(gsub("/", "", unique(immune_list[[i]]@meta.data$CellType2)), "_disease_vs_control", ".csv"), sep =",", quote = F)
}

# Name the list with CT and create a new list with geneID and CT for logFC
temp=list()
for(i in 1:length(heatmap_deg)){
  names(heatmap_deg) <- lapply(immune_list, function(xx){paste(unique(xx@meta.data$CellType2))})
  colnames(heatmap_deg[[i]]) <- c("p_val",names(heatmap_deg[i]),"pct.1","pct.2","p_val_adj")
  temp[[i]] <- heatmap_deg[[i]][,1:2]
  temp[[i]]$GeneID <- rownames(temp[[i]])
}
summary(heatmap_deg) 
summary(temp)

# Sort the dataframe based on geneID
sorted_temp <- lapply(temp, function(df){
  df[order(df$GeneID),]
})

## Prepare data for the heatmap
heatmap_data <- Reduce(
  function(x, y) merge(x, y, by = "GeneID", all = T),
  lapply(sorted_temp, function(x) { x$id <- rownames(x); x }))
heatmap_data$p_val.x <- NULL
heatmap_data$id.x <- NULL
heatmap_data$p_val.y <- NULL
heatmap_data$id.y <- NULL
heatmap_data$p_val <- NULL
heatmap_data$id <- NULL
rownames(heatmap_data) <- heatmap_data[,1]

heatmap_data[,1] <- NULL

## Arrange row names according to the genelist order
heatmap_data$id <- rownames(heatmap_data)
keyDF <- data.frame(key=genelist,weight=1:length(genelist))
merged <- merge(heatmap_data,keyDF,by.x='id',by.y='key',all.x=T,all.y=F)
res <- merged[order(merged$weight),]
res$weight <- NULL
rownames(res) <- res$id
res$id <- NULL

## Binary heatmap
color.palette  <- colorRampPalette(c("ivory2","orange1"))(n=1000)
res[res < 0] <- 0
res[res > 0] <- 1
res$NA_counts <- rowSums(is.na(res)) # count number of NA per gene
res <- res[res$NA_counts <= 6,] # remove genes with NA in more than 50% of cell types

heatmap.2(as.matrix(res[,-13]),
          trace = "none",
          symbreaks = F,
          symm = F,
          symkey = F,
          scale = "none",
          cexCol = 0.8,
          cexRow = 0.6,
          col = color.palette,
          main = "Immune logFC Disease vs. control per CT",
          Colv = F,
          Rowv = F,
          margins = c(10, 9))

binomtest <- capture.output(binom.test(length(res[res > 0]), length(res[res >= 0]), .5))
writeLines(binomtest, con = file("20200626_Immune_heatmap_binomtest_ipf_up.txt"))

# ==============================================================================
# Figure 3C: Boxplot with p-adj-value per diagnosis per CT
# ==============================================================================
# Split the ild object
immune_list = list()
j <- 0
for(i in unique(immune@meta.data$CellType2)){
  j <- j + 1
  immune_list[[j]] <- subset(immune,
                          cells = row.names(immune@meta.data[immune@meta.data$CellType2 == i, ]))
}
for(i in 1:length(immune_list)){
  names(immune_list) <- lapply(immune_list, function(xx){paste(unique(xx@meta.data$CellType2))})
}
summary(immune_list)

# Remove cell types with very few cells
immune_list <- immune_list[names(immune_list) != "TRegs"] # not enough cells for neg binom model
immune_list <- immune_list[names(immune_list) != "pDCs"] # not enough cells for neg binom model

# Create a vector for genes of interest
genelist <- c("IL6R","IL6ST","SOCS1","SOCS2","CCL2","CCL3","CXCL8","AREG","FCGR3A",
              "S100A8","S100A9","GZMB","LAG3")

# Calculate p_adj_value using FindMarkers neg binom test
copd_vs_control <- lapply(immune_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "COPD",
                ident.2 = "Control",
                test.use = "negbinom", # using this test: the slot is set to "counts"
                features = genelist,
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

ipf_vs_control <- lapply(immune_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "IPF",
                ident.2 = "Control",
                test.use = "negbinom",
                features = genelist,
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

other_vs_control <- lapply(immune_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "Other-ILD",
                ident.2 = "Control",
                test.use = "negbinom", 
                features = genelist,
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                min.cells.group = 1,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

# Save the DEG files
for(i in 1:length(immune_list)){
  write.table(copd_vs_control[[i]], paste(gsub("/", "", unique(immune_list[[i]]@meta.data$CellType2)), "_copd_vs_control", ".csv"), sep =",", quote = F)
}
for(i in 1:length(immune_list)){
  write.table(ipf_vs_control[[i]], paste(gsub("/", "", unique(immune_list[[i]]@meta.data$CellType2)), "_ipf_vs_control", ".csv"), sep =",", quote = F)
}
for(i in 1:length(immune_list)){
  write.table(other_vs_control[[i]], paste(gsub("/", "", unique(immune_list[[i]]@meta.data$CellType2)), "_other_vs_control", ".csv"), sep =",", quote = F)
}

saveRDS(copd_vs_control, file = "20200715_Immune_COPD_vs_Control.rds")
saveRDS(ipf_vs_control, file = "20200715_Immune_IPF_vs_Control.rds")
saveRDS(other_vs_control, file = "20200715_Immune_Other_vs_Control.rds")

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
write.csv(pval, file = "20200722_Immune_Boxplot_p_adj_val.csv")

# Add significant into the p_val file
onion <- pval$p
onion[onion <= 0.05] <- "**"
onion[onion > 0.05 & onion <= 0.1] <- "*"
onion[onion > 0.1] <- NA
pval$significant <- onion 

## Extract out the count data for box plot
assayData <- GetAssayData(immune, slot = "counts")
gene1 <- data.frame(assayData[rownames(assayData) == "AREG",])
gene2 <- data.frame(assayData[rownames(assayData) == "CXCL8",])
gene3 <- data.frame(assayData[rownames(assayData) == "CCL3",])
gene4 <- data.frame(assayData[rownames(assayData) == "SOCS1",])
gene5 <- data.frame(assayData[rownames(assayData) == "IL6ST",])

## Adding metadata
colnames(gene1) <- c("counts")
gene1$CellType <- immune@meta.data$CellType2
gene1$Diagnosis <- immune@meta.data$Diagnosis2
gene1$GeneID <- "AREG"

colnames(gene2) <- c("counts")
gene2$CellType <- immune@meta.data$CellType2
gene2$Diagnosis <- immune@meta.data$Diagnosis2
gene2$GeneID <- "CXCL8"

colnames(gene3) <- c("counts")
gene3$CellType <- immune@meta.data$CellType2
gene3$Diagnosis <- immune@meta.data$Diagnosis2
gene3$GeneID <- "CCL3"

colnames(gene4) <- c("counts")
gene4$CellType <- immune@meta.data$CellType2
gene4$Diagnosis <- immune@meta.data$Diagnosis2
gene4$GeneID <- "SOCS1"

colnames(gene5) <- c("counts")
gene5$CellType <- immune@meta.data$CellType2
gene5$Diagnosis <- immune@meta.data$Diagnosis2
gene5$GeneID <- "IL6ST"

# Combine all counts into 1 data.frame
plot_data <- rbind(gene1, gene2, gene3, gene4, gene5)

# Extract out p_adj_val
stat.test <- pval[pval$GeneID %in% unique(plot_data$GeneID),]
write.csv(stat.test, file = "20200813_Immune_Fig3C_negbinom.csv")

# Plotting the cDCs
plot_data.m1 <- plot_data[plot_data$CellType == "cDCs",]
plot_data.m1.1 <- plot_data.m1[plot_data.m1$GeneID == "IL6ST",]

ggplot(plot_data.m1.1, aes(x=Diagnosis, y=counts, color = Diagnosis)) + 
  geom_boxplot(outlier.size = 0) +
  geom_jitter(position = position_jitter(0.2)) + 
  theme_bw() +
  theme(axis.title.x=element_text(size=0), 
        axis.text.x=element_text(size=0,angle = 45, hjust=1)) +
  theme(axis.title.y=element_text(size=0),
        axis.text.y=element_text(size=12)) +
  NoLegend()

# Plotting the macrophages
plot_data.m2 <- plot_data[plot_data$CellType == "Macrophages",]
plot_data.m2.1 <- plot_data.m1[plot_data.m1$GeneID == "IL6ST",]

ggplot(plot_data.m2.1, aes(x=Diagnosis, y=counts, color = Diagnosis)) + 
  geom_boxplot(outlier.size = 0) +
  geom_jitter(position = position_jitter(0.2)) + 
  theme_bw() +
  theme(axis.title.x=element_text(size=0), 
        axis.text.x=element_text(size=0,angle = 45, hjust=1)) +
  theme(axis.title.y=element_text(size=0),
        axis.text.y=element_text(size=12)) +
  NoLegend()

# Plotting the monocytes
plot_data.m3 <- plot_data[plot_data$CellType == "Monocytes",]
plot_data.m3.1 <- plot_data.m1[plot_data.m1$GeneID == "IL6ST",]

ggplot(plot_data.m3.1, aes(x=Diagnosis, y=counts, color = Diagnosis)) + 
  geom_boxplot(outlier.size = 0) +
  geom_jitter(position = position_jitter(0.2)) + 
  theme_bw() +
  theme(axis.title.x=element_text(size=0), 
        axis.text.x=element_text(size=0,angle = 45, hjust=1)) +
  theme(axis.title.y=element_text(size=0),
        axis.text.y=element_text(size=12)) +
  NoLegend()

# ==============================================================================
# Figure 3D,E,F and Figure S8: Immune cell AddModuleScore plots
# ==============================================================================
# ---------------------------
# IFNa genes (Figure S8A)
# ---------------------------
ifn.gene <- list(c("ISG15","IFI44","IFI27","CXCL10","RSAD2","IFIT1","IFI44L",
                   "CCL8","XAF1","GBP1","IRF7","CEACAM1")) # Unterman et al. 2020 (preprint)
immune <- AddModuleScore(object = immune, features = ifn.gene, name = 'IFNscores',
                         assay = "SCT")
ifn_data <- as.data.frame(cbind(immune@meta.data$orig.ident, immune@meta.data$CellType2,
                                immune@meta.data$Diagnosis2,immune@meta.data$IFNscores1))
colnames(ifn_data) <- c("Ident","CellType","Diagnosis","Scores")
ifn_data$Diagnosis <- as.factor(ifn_data$Diagnosis)
levels(ifn_data$Diagnosis)
ifn_data$Scores <- as.numeric(as.character(ifn_data$Scores))

ifn_tukey <- ifn_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(ifn_tukey), file = "20200812_Immune_IFNa_Tukey.csv")

ggplot(ifn_data,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = 0) + #remove outliers, to include outliers change size to 0
  ggtitle ("IFNa scores per diagnosis") +
  theme_bw() +
  # ylim(-0.5,0.5) + # remove this option to include all outliers
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# ---------------------------
# IL6 genes (Figure S8B)
# ---------------------------
il6.gene <- list(c("ARID5A","SOCS3","PIM1","BCL3","BATF","MYC")) # Unterman et al. 2020 (preprint)
immune <- AddModuleScore(object = immune, features = il6.gene, name = 'IL6scores',
                         assay = "SCT")
il6_data <- as.data.frame(cbind(immune@meta.data$orig.ident, immune@meta.data$CellType2,
                                immune@meta.data$Diagnosis2,immune@meta.data$IL6scores1))
colnames(il6_data) <- c("Ident","CellType","Diagnosis","Scores")
il6_data$Diagnosis <- as.factor(il6_data$Diagnosis)
levels(il6_data$Diagnosis)
il6_data$Scores <- as.numeric(as.character(il6_data$Scores))

il6_tukey <- il6_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(il6_tukey), file = "20200812_Immune_IL6_Tukey.csv")

ggplot(il6_data,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = -1) + 
  ggtitle ("IL6 scores per diagnosis") +
  theme_bw() +
  ylim(-0.5,1) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# ---------------------------
# HLA type II genes (Figure 3D, S8C)
# ---------------------------
hla.gene <- list(c("HLA-DRA", "HLA-DQA1","HLA-DQA2","HLA-DPA1","HLA-DRB1","HLA-DPB1",
                   "HLA-DQB2","HLA-DRB5","HLA-DQB1","HLA-DMA","HLA-DMB")) # Unterman et al. 2020 (preprint)
immune <- AddModuleScore(object = immune, features = hla.gene, name = 'HLAscores',
                         assay = "SCT")
hla_data <- as.data.frame(cbind(immune@meta.data$orig.ident, immune@meta.data$CellType2,
                                immune@meta.data$Diagnosis2,immune@meta.data$HLAscores1))
colnames(hla_data) <- c("Ident","CellType","Diagnosis","Scores")
hla_data$Diagnosis <- as.factor(hla_data$Diagnosis)
levels(hla_data$Diagnosis)
hla_data$Scores <- as.numeric(as.character(hla_data$Scores))

hla_tukey <- hla_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(hla_tukey), file = "20200812_Immune_HLA_Tukey.csv")

ggplot(hla_data,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = -1) + 
  ggtitle ("HLA Type II scores per diagnosis") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# ---------------------------
# Cytotoxic scores (Figure 3E & Figure S8D)
# ---------------------------
# Unterman et al. 2020 (preprint): Cytotoxic T cell markers (NKG7, KLRG1, PRF1, GZMH)
# Unterman et al. 2020 (preprint): Cytotoxic/pro-inflammatory cytokines (PRF1, GZMH, IFNG)
# Szabo et al., 2019 (Nat Comm.): Cytotoxicity associated genes (GNLY, GZMB, GZMK)
cyto.gene <- list(c("PRF1", "GZMH", "IFNG","NKG7", "KLRG1", "PRF1","GNLY","GZMB",
                    "GZMK")) 
immune <- AddModuleScore(object = immune, features = cyto.gene, name = 'Cytotoxicityscores',
                         assay = "SCT")
cyto_data <- as.data.frame(cbind(immune@meta.data$orig.ident, immune@meta.data$CellType2,
                                 immune@meta.data$Diagnosis2,immune@meta.data$Cytotoxicityscores1))
colnames(cyto_data) <- c("Ident","CellType","Diagnosis","Scores")
cyto_data <- cyto_data[cyto_data$CellType %in% c("CD4 T Cells", "CD8 T Cells",
                                                 "NK Cells","Proliferating T Cells",
                                                 "TRegs"),]
cyto_data$Diagnosis <- as.factor(cyto_data$Diagnosis)
levels(cyto_data$Diagnosis)
cyto_data$Scores <- as.numeric(as.character(cyto_data$Scores))

cyto_tukey <- cyto_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(cyto_tukey), file = "20200812_Immune_Cytotoxic_Tukey_plot2.csv")

ggplot(cyto_data,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = -1) + 
  ggtitle ("Cytotoxicity scores per diagnosis (PRF1,GZMH,IFNG,NKG7,KLRG1,PRF1,GNLY,
           GZMB,GZMK)") +
  theme_bw() +
  ylim(-0.5,2.5) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# ---------------------------
# T Cells exhaustion scores (Figure 3F & Figure S8E)
# ---------------------------
exhaustion.gene <- list(c("LAG3","TIGIT","PDCD1","CTLA4","HAVCR2","TOX", # Zhanng et al. 2020 (preprint)
                          "PRDM1","MAF"))  # Unterman et al., 2020 (preprint)
immune <- AddModuleScore(object = immune, features = exhaustion.gene, name = 'Exhaustionscores',
                         assay = "SCT")
exha_data <- as.data.frame(cbind(immune@meta.data$orig.ident, immune@meta.data$CellType2,
                                 immune@meta.data$Diagnosis2,immune@meta.data$Exhaustionscores1))
colnames(exha_data) <- c("Ident","CellType","Diagnosis","Scores")
exha_data <- exha_data[exha_data$CellType %in% c("CD4 T Cells","CD8 T Cells","TRegs",
                                                 "NK Cells"),]
exha_data$Diagnosis <- as.factor(exha_data$Diagnosis)
levels(exha_data$Diagnosis)
exha_data$Scores <- as.numeric(as.character(exha_data$Scores))

exha_tukey <- exha_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(exha_tukey), file = "20200812_Immune_Exhaustion_Tukey.csv")

ggplot(exha_data,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = -1) + 
  ggtitle ("Exhaustion scores per diagnosis") +
  theme_bw() +
  ylim(-0.25,0.6)+
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 