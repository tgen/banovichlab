# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-07-03
# Title: scRNA-seq data analysis for Covid-19 manuscript - Figure 2
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

# ==============================================================================
# Read in the Epithelial object (containing all cell types)
# ==============================================================================
epi <- readRDS("/scratch/lbui/Covid19_Seurat/20200708_Epi_annotated_noDoublets.rds")

# ==============================================================================
# Figure 2A: Binary heatmap
# ==============================================================================
# Split the Epithelial object into a list of different CT objects
epi_list = list()
j <- 0
for(i in unique(epi@meta.data$CellType1)){
  j <- j + 1
  epi_list[[j]] <- subset(epi,
                          cells = row.names(epi@meta.data[epi@meta.data$CellType1 == i, ]))
}
for(i in 1:length(epi_list)){
  names(epi_list) <- lapply(epi_list, function(xx){paste(unique(xx@meta.data$CellType1))})
}
summary(epi_list)

# Run DEG for Covid-19 candidate genes      
genelist <- c("ACE2","BSG","TMPRSS2","CTSL","CTSB","FURIN","PCSK5","PCSK7",
              "ADAM17","PIKFYVE","TPCN2","AGT","ACE","ITGB6","C1QA","C1QB",
              "C1QC","C2","C3","C4B","IFNAR1","IFNAR2","IFNGR1",
              "IFNGR2","PTPN11","EIF2AK2","EIF2AK3","CXCL1","TRIM27","TRIM28","NFKB1",
              "RNF41","JUN","SOCS1","SOCS2","CSF2","CSF3","ICAM1","CD47",
              "CD44","CCL2","CCL3","FGA","FGG","ATG5","ATG7","BECN1","SQSTM1",
              "MAP1LC3A","MAP1LC3B","ATF6","ERN1","MUC5B") #remove NT5DC1 since it's in non-human primate (Ziegler et al., 2020. CELL)

heatmap_deg <- lapply(epi_list, function(xx){
  print(unique(xx@meta.data$CellType1))
  FindMarkers(xx,
              group.by = "Status",
              ident.1 = "Disease",
              ident.2 = "Control",
              features = genelist,
              test.use = "negbinom",
              latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
              logfc.threshold = 0,
              min.pct = 0)
})

# Save the DEG files
for(i in 1:length(epi_list)){
   write.table(heatmap_deg[[i]], paste(gsub("/", "", unique(epi_list[[i]]@meta.data$CellType1)), 
                                       "_disease_vs_control", ".csv"), sep =",", quote = F)}

  # Name the list with CT and create a new list with geneID and CT for logFC
temp=list()
for(i in 1:length(heatmap_deg)){
  names(heatmap_deg) <- lapply(epi_list, function(xx){paste(unique(xx@meta.data$CellType1))})
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

# Prepare data for the heatmap
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

# Arrange row names according to the genelist order
heatmap_data$id <- rownames(heatmap_data)
keyDF <- data.frame(key=genelist,weight=1:length(genelist))
merged <- merge(heatmap_data,keyDF,by.x='id',by.y='key',all.x=T,all.y=F)
res <- merged[order(merged$weight),]
res$weight <- NULL
rownames(res) <- res$id
res$id <- NULL

## Binary heatmap
color.palette  <- colorRampPalette(c("ivory2","orange1"))(n=1000)
test <- res
test[test < 0] <- 0
test[test > 0] <- 1

# Rearrange cell type orders
test <- test[c("AT1","AT2","Transitional AT2","KRT5-/KRT17+","Goblet Cells",
               "Club Cells","PNEC/Ionocytes","Ciliated Cells",
               "Differentiating Ciliated Cells","Basal")]

heatmap.2(as.matrix(test),
          trace = "none",
          symbreaks = F,
          symm = F,
          symkey = F,
          scale = "none",
          labCol = "",
          cexRow = 0.6,
          col = color.palette,
          main = "Epithelial logFC Disease vs. Control per CT",
          Colv = F,
          Rowv = F,
          add.expr = text(x = seq_along(colnames(test)), y = -2, srt = 45,
                          labels = colnames(test), xpd = TRUE))

binomtest <- capture.output(binom.test(length(test[test > 0]), length(test[test>=0]), .5))
writeLines(binomtest, con = file("20200903_Epi_heatmap_binomtest_disease_control.txt"))

# ==============================================================================
# Figure 2B: AddModuleScore for SARS-CoV-2 entry genes
# ==============================================================================
# AddModuleScore for SARS-CoV-2 genes
entry.gene <- list(c("ACE2","BSG","HSPA5","TMPRSS2","CTSL","FURIN","ADAM17"))
entry.gene2 <- list(c("ACE2","TMPRSS2","CTSL","FURIN","ADAM17")) # Figure S5B

epi <- AddModuleScore(object = epi, features = entry.gene, 
                      name = 'Entryscores', assay = "SCT")
epi <- AddModuleScore(object = epi, features = entry.gene2, 
                      name = 'Entry2scores', assay = "SCT")
# Make plots
VlnPlot(object = epi, features = 'Entryscores1', group.by = "CellType2", 
        pt.size = 0, split.by = "Status") + NoLegend()

# Perform statistic test for Entry scores - 1
entry_data <- as.data.frame(cbind(epi@meta.data$orig.ident, epi@meta.data$CellType1,
                                  epi@meta.data$Diagnosis2,epi@meta.data$Entryscores1))
colnames(entry_data) <- c("Ident","CellType","Diagnosis","Scores")

entry_data$Diagnosis <- as.factor(entry_data$Diagnosis)
levels(entry_data$Diagnosis)
entry_data$Scores <- as.numeric(as.character(entry_data$Scores))

entry_tukey <- entry_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(entry_tukey), file = "20200908_Epi_entryscore1_Tukey.csv")

## Boxplot
celltype.plot <-  c("AT1","AT2","Transitional AT2","KRT5-/KRT17+","Goblet Cells",
                    "Club Cells", "PNECs/Ionocytes","Ciliated Cells",
                    "Differentiating Ciliated Cells","Basal")
entry_data.df <- entry_data[entry_data$CellType %in% celltype.plot,]
entry_data.df$CellType <- factor(entry_data.df$CellType, # reorganize x-axis order
                                 levels=celltype.plot)

ggplot(entry_data.df,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = -1) + 
  ylim(-0.5,0.75) +
  ggtitle ("SARS-CoV-2 entry gene scores, with BSG, HSPA5") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# Perform statistic test for Entry scores - 2
entry_data2 <- as.data.frame(cbind(epi@meta.data$orig.ident, epi@meta.data$CellType1,
                                  epi@meta.data$Diagnosis2,epi@meta.data$Entry2scores1))
colnames(entry_data2) <- c("Ident","CellType","Diagnosis","Scores")

entry_data2$Diagnosis <- as.factor(entry_data2$Diagnosis)
levels(entry_data2$Diagnosis)
entry_data2$Scores <- as.numeric(as.character(entry_data2$Scores))

entry_tukey2 <- entry_data2 %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(entry_tukey2), file = "20200908_Epi_entryscore2_Tukey.csv")

## Boxplot
entry_data2.df <- entry_data2[entry_data2$CellType %in% celltype.plot,]
entry_data2.df$CellType <- factor(entry_data2.df$CellType, # reorganize x-axis order
                                 levels=celltype.plot)

ggplot(entry_data2.df,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = 0) + 
  ggtitle ("SARS-CoV-2 entry gene scores-No BSG, HSPA5") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# ==============================================================================
# Figure 2C: Gene expression correlation
# ==============================================================================
# Included in the scripts for Fig S4

# ==============================================================================
# Figure 2D: Boxplot for some candidates AT2 cells
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
epi_list <- epi_list[names(epi_list) != "Ionocytes"]
epi_list <- epi_list[names(epi_list) != "PNECs"]

genelist2 <- c("ACE2","TMPRSS2","CTSL","FURIN","BSG","ADAM17","HSPA5","ITGB6","NT5DC1",
               "SOCS1","SOCS2","CSF3")

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
write.csv(pval, file = "20200722_Epi_Boxplot_p_adj_val.csv")

# Add significant into the p_val file
onion <- pval$p
onion[onion <= 0.05] <- "**"
onion[onion > 0.05 & onion <= 0.1] <- "*"
onion[onion > 0.1] <- NA
pval$significant <- onion 

# Extract out gene counts
assayData <- GetAssayData(epi, slot = "counts")
gene1 <- data.frame(assayData[rownames(assayData) == "ITGB6",])
gene2 <- data.frame(assayData[rownames(assayData) == "SOCS1",])
gene3 <- data.frame(assayData[rownames(assayData) == "SOCS2",])
gene4 <- data.frame(assayData[rownames(assayData) == "CSF3",])

## Adding metadata
colnames(gene1) <- c("counts")
gene1$CellType <- epi@meta.data$CellType2
gene1$Diagnosis <- epi@meta.data$Diagnosis2
gene1$GeneID <- "ITGB6"

colnames(gene2) <- c("counts")
gene2$CellType <- epi@meta.data$CellType2
gene2$Diagnosis <- epi@meta.data$Diagnosis2
gene2$GeneID <- "SOCS1"

colnames(gene3) <- c("counts")
gene3$CellType <- epi@meta.data$CellType2
gene3$Diagnosis <- epi@meta.data$Diagnosis2
gene3$GeneID <- "SOCS2"

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
stat.test <- pval[pval$GeneID %in% c("ITGB6","SOCS1","SOCS2","CSF3"),]
stat.test <- stat.test[stat.test$CellType == "AT2",]
stat.test <- stat.test[order(stat.test$GeneID),]
stat.test$y.position <- c(40,45,50,22,24,26,12,14,16,16,18,20)

# Plotting
plot_data.df2 <- plot_data.df[plot_data.df$GeneID == "CSF3",]
stat.test2 <- stat.test[stat.test$GeneID == "CSF3",]
ggplot(plot_data.df2, aes(x=Diagnosis, y=counts, color = Diagnosis)) + 
  geom_boxplot(outlier.size = -1,width = 0.75, na.rm = T) +
  geom_jitter(position = position_jitter(0.15)) + 
  theme_bw() +
  stat_pvalue_manual(stat.test2, size=4, tip.length = 0.01, step.group.by = "GeneID",
                 label = "significant") +
  theme(axis.title.x=element_text(size=15), 
        axis.text.x=element_text(size=0)) +
  ylim(0,50) + 
  theme(axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=10))

# ==============================================================================
# Figure 2E: ACE2+ correlated gene analysis
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
plot_data <- plot_data[!duplicated(plot_data$geneID),]
plot_data <- plot_data[plot_data$geneID != "ACE2",]
temp1 <- grep( "^MT-", rownames(plot_data), ignore.case = F, value = T)
temp2 <- grep( "^RP", rownames(plot_data), ignore.case = F, value = T)
plot_data <- plot_data[!plot_data$geneID %in% temp1,]
plot_data <- plot_data[!plot_data$geneID %in% temp2,]

copd.gene <- c("OAS1","CLEC3B","ITGB7","WNT11","IGSF6","ACKR4","MST1R","LST1",
               "OASL","ILLN2","IFI6","IFI27","IFIT1","IFIT2","IFIT3")
ipf.gene <- c("IL31RA","SOX2","DDX17","NXF3","RMI1","GBP6","KIF19","COL1A1",
              "MAPK9","ITPKB","IGFBP2","IGFBP7","GLRA3","CLDN10")
control.gene <- c("SCN7A","AC245052.4","CTA-204B4.2")
other.gene <- c("IL10","ITLN1","CD83","IGFBP2","NRTN","MUC4","KRT13",
                "SFTPB")

plot_data$genelabels <- ""
plot_data$genelabels <- ifelse(plot_data$geneID %in% copd.gene & plot_data$Set == "COPD" 
                               |plot_data$geneID %in% ipf.gene & plot_data$Set == "IPF"
                               |plot_data$geneID %in% other.gene & plot_data$Set == "Other-ILD"
                               |plot_data$geneID %in% control.gene & plot_data$Set == "Control"
                               , TRUE,FALSE)

# Plot correlation 
pdf("20200904_ACE2_correlation_SCT_data.pdf")
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
                                 plot_data$geneID,"")) +
  # geom_text_repel(data=subset(plot_data,abs(rho) >= 0.5 & -log10(p_value) >= 2),
  #                  aes(-log10(p_value), rho, label = geneID),
  #                 size = 3) +
  ggtitle("ACE2 correlated genes")
dev.off()

write.csv(plot_data, file = "20200908_ACE2_correlation_SCT_data.csv")


