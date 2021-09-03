# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-07-03
# Title: scRNA-seq data analysis for Covid-19 manuscript - Figure 2
# Update on 20210124 to update Figure 2d, GO analysis for 2e and add ModuleScore for responsive factors
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
# Figure 2A: Binary heatmap
# ==============================================================================
# Split the Epithelial object into a list of different CT objects
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

# Run DEG for Covid-19 candidate genes
# Update gene list on 20210125 to include ISGs in Martin-Sancho (preprint, 2020)
genelist <- c("ACE2","BSG","TMPRSS2","NPR1","CTSL","CTSB","FURIN","PCSK5","PCSK7",
              "ADAM17","PIKFYVE","TPCN2","AGT","ACE","ITGB6","IFNAR1","IFNAR2",
              "EIF2AK2","EIF2AK3","CD44","IFNGR1","IFNGR2","FAM46C","UBD","REC8",
              "ELF1","CLEC4D","LY6E","SPATS2L","ZBP1","DNAJC6","IFIT3","RGS22",
              "B4GALT5","ISG20","GNB4","SPATA13","NRN1","ERLIN1","APOL2","RAB27A",
              "FZD5","C9orf91","TAGAP","HSPA8","CNP","ETV6","MSR1","BST2","CXCL1",
              "CSF2","CSF3","ICAM1","CD47","CCL2","CCL3","TRIM27","TRIM28","RNF41",
              "NFKB1","JUN","SOCS1","SOCS2","C1QA","C1QB","C1QC","C2","C3","C4B",
              "PTPN11","FGA","FGG","ATG5","ATG7","BECN1","SQSTM1",
              "MAP1LC3A","MAP1LC3B","ATF6","ERN1","MUC5B") 

heatmap_deg <- lapply(epi_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  FindMarkers(xx,
              group.by = "Status",
              ident.1 = "Disease",
              ident.2 = "Control",
              features = genelist,
              test.use = "negbinom",
              latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
              logfc.threshold = 0,
              min.pct = 0,
              assay ="SCT")
})

# Save the DEG files
for(i in 1:length(epi_list)){
   write.table(heatmap_deg[[i]], paste(gsub("/", "", unique(epi_list[[i]]@meta.data$CellType2)), 
                                       "_disease_vs_control", ".csv"), sep =",", quote = F)}
saveRDS(heatmap_deg, file = "20210206_Epi_Fig2A_deg.rds")

  # Name the list with CT and create a new list with geneID and CT for logFC
temp=list()
for(i in 1:length(heatmap_deg)){
  names(heatmap_deg) <- lapply(epi_list, function(xx){paste(unique(xx@meta.data$CellType2))})
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
color.palette  <- colorRampPalette(c("ivory2","orange1","light grey"))(n=1000)
res[res < 0] <- 0
res[res > 0] <- 1
res[is.na(res)] <- 3 # assign NA value to 3 so can plot with grey color

# Rearrange cell type orders
res <- res[c("AT1","AT2","Transitional AT2","KRT5-/KRT17+","MUC5B+","MUC5AC+",
              "SCGB3A2+/SCGB1A1+","SCGB3A2+","PNEC/Ionocytes","Ciliated Cells",
              "Basal")]

## Updated on 04/12/2021 to test only genes with significant p_adj_value (<=0.1)
# Name the list with CT and create a new list with geneID and CT for p_adj_values
heatmap_deg2 = heatmap_deg
temp2=list()
for(i in 1:length(heatmap_deg2)){
  names(heatmap_deg2) <- lapply(epi_list, function(xx){paste(unique(xx@meta.data$CellType2))})
  colnames(heatmap_deg2[[i]]) <- c("p_val","avg_log2FC","pct.1","pct.2",names(heatmap_deg2[i]))
  temp2[[i]] <- heatmap_deg2[[i]][,4:5]
  temp2[[i]]$GeneID <- rownames(temp2[[i]])
}
summary(heatmap_deg2) 
summary(temp2)

# Sort the dataframe based on geneID
sorted_temp2 <- lapply(temp2, function(df){
  df[order(df$GeneID),]
})

# Prepare data for the heatmap
heatmap_data2 <- Reduce(
  function(x, y) merge(x, y, by = "GeneID", all = T),
  lapply(sorted_temp2, function(x) { x$id <- rownames(x); x }))
heatmap_data2$pct.2.x <- NULL
heatmap_data2$id.x <- NULL
heatmap_data2$pct.2.y <- NULL
heatmap_data2$id.y <- NULL
heatmap_data2$pct.2 <- NULL
heatmap_data2$id <- NULL
rownames(heatmap_data2) <- heatmap_data2[,1]
heatmap_data2[,1] <- NULL

# Only keep the genes with p_adj_val <= 0.1 for the binomial test
test1 <- ifelse(heatmap_data2 <= 0.1, TRUE, FALSE)
test2 <- data.frame(matrix(nrow = 81, ncol = 11))
i=0
for(k in(1:length(rownames(test1)))){
  i=i+1
  test2[i,] <- ifelse(test1[i,]==TRUE, heatmap_data[i,],-1)
}
colnames(test2) <- colnames(test1)
rownames(test2) <- rownames(test1)

test2[is.na(test2)] <- -2 # asign NA to -2
binomtest <- capture.output(binom.test(length(test2[test2 > 0]), length(test2[test2 > -1]), .5,
                                       conf.level = 0.95, alternative = "greater"))
writeLines(binomtest, con = file("20210413_Epi_heatmap_binomtest_disease_control.txt"))

# Binomial test for AT2 only
test2.df <- as.data.frame(cbind(rownames(test2),test2[,2]))
rownames(test2.df) <- test2.df$V1
test2.df$V1 <- NULL
test2.df$V2 <- as.numeric(test2.df$V2)
binomtest <- capture.output(binom.test(length(test2.df[test2.df > 0]), length(test2.df[test2.df > -1]), .5,
                                       conf.level = 0.95, alternative = "greater"))
writeLines(binomtest, con = file("20210416_Epi_AT2_binomtest_disease_control.txt"))

# Add border for significant p_val_adj
test1 <- as.data.frame(ifelse(heatmap_data2 <= 0.1, 1, 0))
test1 <- test1 %>%
  slice(match(genelist, rownames(test1)))
test1 <- test1[c("AT1","AT2","Transitional AT2","KRT5-/KRT17+","MUC5B+","MUC5AC+",
             "SCGB3A2+/SCGB1A1+","SCGB3A2+","PNEC/Ionocytes","Ciliated Cells",
             "Basal")]
test1 <- apply(test1, 2, rev)
test1 <- t(test1)
nx=11
ny=81
makeRects <- function(tfMat,border){
  cAbove = expand.grid(1:nx,1:ny)[tfMat,]
  xl=cAbove[,1]-0.49
  yb=cAbove[,2]-0.49
  xr=cAbove[,1]+0.49
  yt=cAbove[,2]+0.49
  rect(xl,yb,xr,yt,border=border,lwd=1)
}            #this is the function to make the rectangles/borders

pdf("plot.pdf")
heatmap.2(as.matrix(res),
          trace = "none",
          symbreaks = F,
          symm = F,
          symkey = F,
          scale = "none",
          labCol = "",
          cexRow = 0.3,
          col = color.palette,
          main = "Epithelial logFC Disease vs. Control per CT",
          Colv = F,
          Rowv = F,
          add.expr = text(x = seq_along(colnames(res)), y = -2, srt = 45,
                          labels = colnames(res), xpd = TRUE,
                          {makeRects(as.matrix(test1),"black")}))
dev.off()

# ==============================================================================
# Figure 2C: AddModuleScore for SARS-CoV-2 entry genes
# ==============================================================================
# AddModuleScore for SARS-CoV-2 genes
entry.gene <- list(c("ACE2","BSG","HSPA5","TMPRSS2","CTSL","FURIN","ADAM17","NRP1"))
entry.gene2 <- list(c("ACE2","TMPRSS2","CTSL","FURIN","ADAM17")) # Figure S5B
res.gene <- list("LY6E","CLEC4D","UBD","ELF1","FAM46C","REC8") # viral entry restriction genes - Martin-Sancho 2020
rep.gene <- list("SPATS2L","ZBP1","DNAJC6","IFIT3","RGS22","IFIT1","IFIT5","B4GALT5") # RNA replication inhibition genes - Martin-Sancho 2020

epi <- AddModuleScore(object = epi, features = entry.gene, 
                      name = 'Entryscores', assay = "SCT")
epi <- AddModuleScore(object = epi, features = entry.gene2, 
                      name = 'Entry2scores', assay = "SCT")
epi <- AddModuleScore(object = epi, features = res.gene, 
                      name = 'Restrictionscores', assay = "SCT")
epi <- AddModuleScore(object = epi, features = rep.gene, 
                      name = 'Replicationscores', assay = "SCT")

# Make plots
VlnPlot(object = epi, features = 'Entryscores1', group.by = "CellType2", 
        pt.size = 0, split.by = "Status") + NoLegend()

# Perform statistic test for Entry scores - 1
entry_data <- as.data.frame(cbind(epi@meta.data$orig.ident, epi@meta.data$CellType2,
                                  epi@meta.data$Diagnosis2,epi@meta.data$Entryscores1))
colnames(entry_data) <- c("Ident","CellType","Diagnosis","Scores")

entry_data$Diagnosis <- as.factor(entry_data$Diagnosis)
levels(entry_data$Diagnosis)
entry_data$Scores <- as.numeric(as.character(entry_data$Scores))

entry_tukey <- entry_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(entry_tukey), file = "20210208_Epi_entryscore1_Tukey.csv")

## Boxplot
celltype.plot <-  c("AT1","AT2","Transitional AT2","KRT5-/KRT17+","MUC5B+",
                    "MUC5AC+","SCGB3A2+/SCGB1A1+", "SCGB3A2+","PNEC/Ionocytes",
                    "Ciliated Cells","Basal")
entry_data.df <- entry_data[entry_data$CellType %in% celltype.plot,]
entry_data.df$CellType <- factor(entry_data.df$CellType, # reorganize x-axis order
                                 levels=celltype.plot)

ggplot(entry_data.df,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = -1) + 
  ylim(-0.5,0.75) +
  ggtitle ("SARS-CoV-2 entry gene scores, with BSG, HSPA5, NRP1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# Perform statistic test for Entry scores - 2
entry_data2 <- as.data.frame(cbind(epi@meta.data$orig.ident, epi@meta.data$CellType2,
                                  epi@meta.data$Diagnosis2,epi@meta.data$Entry2scores1))
colnames(entry_data2) <- c("Ident","CellType","Diagnosis","Scores")

entry_data2$Diagnosis <- as.factor(entry_data2$Diagnosis)
levels(entry_data2$Diagnosis)
entry_data2$Scores <- as.numeric(as.character(entry_data2$Scores))

entry_tukey2 <- entry_data2 %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(entry_tukey2), file = "20210208_Epi_entryscore2_Tukey.csv")

## Boxplot
entry_data2.df <- entry_data2[entry_data2$CellType %in% celltype.plot,]
entry_data2.df$CellType <- factor(entry_data2.df$CellType, # reorganize x-axis order
                                 levels=celltype.plot)

ggplot(entry_data2.df,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = 0) + 
  ggtitle ("SARS-CoV-2 entry gene scores-No BSG, HSPA5") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# Perform statistic test for restriction factor score
res_data <- as.data.frame(cbind(epi@meta.data$orig.ident, epi@meta.data$CellType2,
                                  epi@meta.data$Diagnosis2,epi@meta.data$Restrictionscores1))
colnames(res_data) <- c("Ident","CellType","Diagnosis","Scores")

res_data$Diagnosis <- as.factor(res_data$Diagnosis)
levels(res_data$Diagnosis)
res_data$Scores <- as.numeric(as.character(res_data$Scores))

res_tukey <- res_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(res_tukey), file = "20210208_Epi_restrictionscore1_Tukey.csv")

## Boxplot
res_data.df <- res_data[res_data$CellType %in% celltype.plot,]
res_data.df$CellType <- factor(res_data.df$CellType, # reorganize x-axis order
                                 levels=celltype.plot)

ggplot(res_data.df,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = 0) + 
 # ylim(-0.5,0.75) +
  ggtitle ("SARS-CoV-2 entry restriction gene scores") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

# Perform statistic test for replication factor score
rep_data <- as.data.frame(cbind(epi@meta.data$orig.ident, epi@meta.data$CellType2,
                                epi@meta.data$Diagnosis2,epi@meta.data$Replicationscores1))
colnames(rep_data) <- c("Ident","CellType","Diagnosis","Scores")

rep_data$Diagnosis <- as.factor(rep_data$Diagnosis)
levels(rep_data$Diagnosis)
rep_data$Scores <- as.numeric(as.character(rep_data$Scores))

rep_tukey <- rep_data %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Scores ~ Diagnosis)

write.csv(as.data.frame(rep_tukey), file = "20210208_Epi_replicationscore1_Tukey.csv")

## Boxplot
rep_data.df <- rep_data[rep_data$CellType %in% celltype.plot,]
rep_data.df$CellType <- factor(rep_data.df$CellType, # reorganize x-axis order
                               levels=celltype.plot)

ggplot(rep_data.df,aes(x=CellType, y=Scores, fill=Diagnosis)) +
  geom_boxplot(outlier.size = 0) + 
  # ylim(-0.5,0.75) +
  ggtitle ("SARS-CoV-2 entry replication inhibition gene scores") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 


