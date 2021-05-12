# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-07-18
# Title: scRNA-seq data analysis for Covid-19 manuscript - Supplementary Figures
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
library(RColorBrewer)
library(VennDiagram)
library(rstatix)

# ==============================================================================
# Supplementary Figure 1 -  DotPlot for cell type markers
# ==============================================================================
# Read in the ILD object (containing all cell types)
ild <- readRDS("/scratch/lbui/Covid19_ILD_objects/20210211_ILD_noDoublets.rds")

# Order cell types along the x axis according to population 
cell_order <- c("Lymphatic Endothelial Cells","Vascular Endothelial Cells","AT1",
                "AT2","Transitional AT2","KRT5-/KRT17+","Basal","Ciliated Cells",
                "SCGB3A2+","SCGB3A2+/SCGB1A1+","MUC5AC+","MUC5B+","PNEC/Ionocytes",
                "Club Cells", "Goblet Cells","B Cells","cDCs","Macrophages",
                "Proliferating Macrophages","Mast Cells","Monocytes","NK Cells",
                "pDCs","Plasma Cells","CD4 T Cells","CD8 T Cells" ,"Proliferating T Cells",
                "Fibroblasts","Myofibroblasts","HAS1 High Fibroblasts",
                "PLIN2+ Fibroblasts","Mesothelial","Pericytes","Smooth Muscle Cells")

ild$CellType2 <- factor(ild$CellType2, levels=cell_order)

# Make dotplot
DotPlot(ild, features = c("TAGLN","ACTG2","GJA4","RGS5","HP","WT1","LUM","PDGFRA",
                          "CD3E","JCHAIN","IGHG1","LILRA4","CLEC4C","NKG7","NCR1",
                          "S100A12","FCN1","CPA3","KIT","MARCO","APOC1","CLEC9A", 
                          "CD1C","MS4A1","CD19","MUC5B","SCGB3A2","SCGB1A1",
                          "FOXI1","FOXJ1","KRT5","KRT17","SFTPC","ABCA3","AGER",
                          "VWF","CCL21"), 
        group.by = "CellType2") + 
  theme(axis.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# ==============================================================================
# Supplementary Figure 3 - BSG+, CTSL+ and HSPA5+ cell proportion
# ==============================================================================
# Extract out total number of cells per ident
table1 <- as.data.frame(table(ild$orig.ident, ild$CellType1, ild$Diagnosis2))

# BSG+ cells
sub1 <- subset(ild, BSG > 0, slot = "counts") 
tableS1 <- as.data.frame(table(sub1$orig.ident,sub1$CellType1, sub1$Diagnosis2))
data_tableS1 <- merge(table1, tableS1, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS1) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS1$Percent <- data_tableS1$Count/data_tableS1$Total * 100
data_tableS1$geneid <- "BSG"

# CTSL+ cells
sub2 <- subset(ild, CTSL > 0, slot = "counts")
tableS2 <- as.data.frame(table(sub2$orig.ident,sub2$CellType1, sub2$Diagnosis2))
data_tableS2 <- merge(table1, tableS2, by = c("Var1", "Var2","Var3"))
colnames(data_tableS2) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS2$Percent <- data_tableS2$Count/data_tableS2$Total * 100
data_tableS2$geneid <- "CTSL"

# HSPA5+ cells
sub3 <- subset(ild, HSPA5 > 0, slot = "counts")
tableS3 <- as.data.frame(table(sub3$orig.ident,sub3$CellType1, sub3$Diagnosis2))
data_tableS3 <- merge(table1, tableS3, by = c("Var1", "Var2","Var3"))
colnames(data_tableS3) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS3$Percent <- data_tableS3$Count/data_tableS3$Total * 100
data_tableS3$geneid <- "HSPA5"

# FURIN+ cells
sub4 <- subset(ild, FURIN > 0, slot = "counts")
tableS4 <- as.data.frame(table(sub4$orig.ident,sub4$CellType1, sub4$Diagnosis2))
data_tableS4 <- merge(table1, tableS4, by = c("Var1", "Var2","Var3"))
colnames(data_tableS4) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS4$Percent <- data_tableS4$Count/data_tableS4$Total * 100
data_tableS4$geneid <- "FURIN"

# NRP1+ cells
sub5 <- subset(ild, NRP1 > 0, slot = "counts")
tableS5 <- as.data.frame(table(sub5$orig.ident,sub5$CellType1, sub5$Diagnosis2))
data_tableS5 <- merge(table1, tableS5, by = c("Var1", "Var2","Var3"))
colnames(data_tableS5) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS5$Percent <- data_tableS5$Count/data_tableS5$Total * 100
data_tableS5$geneid <- "NRP1"

# Combine all tables and perform dplyr mean calculations
onion <- rbind(data_tableS1, data_tableS5, data_tableS2, data_tableS3, data_tableS4)
onion_means = onion %>% group_by(geneid, CellType, Diagnosis) %>% dplyr::summarise(Mean = mean(Percent, na.rm = T),
                                                                                   n=n(),
                                                                                   sd = sd(Percent, na.rm = T),
                                                                                   se = sd/sqrt(n))
onion_means = as.data.frame(onion_means)
onion_means = onion_means[onion_means$Mean > 0,]

# Add a column for total of cells expressing each gene per CT
bsg_total <- as.data.frame(table(sub1@meta.data$CellType1, sub1$Diagnosis2))
colnames(bsg_total) <- c("CellType","Diagnosis","Total")
bsg_total$geneid <- "BSG"

ctsl_total <- as.data.frame(table(sub2@meta.data$CellType1, sub2$Diagnosis2))
colnames(ctsl_total) <- c("CellType","Diagnosis","Total")
ctsl_total$geneid <- "CTSL"

hspa5_total <- as.data.frame(table(sub3@meta.data$CellType1, sub3$Diagnosis2))
colnames(hspa5_total) <- c("CellType","Diagnosis","Total")
hspa5_total$geneid <- "HSPA5"

furin_total <- as.data.frame(table(sub4@meta.data$CellType1, sub4$Diagnosis2))
colnames(furin_total) <- c("CellType","Diagnosis","Total")
furin_total$geneid <- "FURIN"

nrp1_total <- as.data.frame(table(sub5@meta.data$CellType1, sub5$Diagnosis2))
colnames(nrp1_total) <- c("CellType","Diagnosis","Total")
nrp1_total$geneid <- "NRP1"

total <- rbind(bsg_total, nrp1_total, ctsl_total, hspa5_total, furin_total)

onion_means <- merge(onion_means, total, by = c("CellType","geneid", "Diagnosis"))

# Tukey test
onion$Diagnosis <- as.factor(onion$Diagnosis)
levels(onion$Diagnosis)

sp_tukey <- onion %>% group_by(CellType,geneid) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis)

# Save the table for supplemental info
write.csv (onion, file = "2020216_FigS2_percentage.allCT.csv")
write.csv (onion_means, file = "20210216_FigS2_percentage_means.allCT.csv")
write.csv (as.data.frame(sp_tukey), file = "20210216_FigS2_tukey.allCT.csv")

# Order CT on the y axis
cell_order <- c("Lymphatic Endothelial Cells","Vascular Endothelial Cells","AT1",
                "AT2","Transitional AT2","KRT5-/KRT17+","Basal","Ciliated Cells",
                "PNEC/Ionocytes","Club Cells", "Goblet Cells","B Cells","cDCs","Macrophages",
                "Proliferating Macrophages","Mast Cells","Monocytes","NK Cells",
                "pDCs","Plasma Cells","T Cells", "Proliferating T Cells",
                "Fibroblasts","Myofibroblasts","HAS1 High Fibroblasts",
                "PLIN2+ Fibroblasts","Mesothelial","Pericytes","Smooth Muscle Cells")
onion_means$CellType <- factor(onion_means$CellType, 
                               levels=cell_order)

gene_order <- c("BSG","NRP1","HSPA5","CTSL","FURIN")
onion_means$geneid <- factor(onion_means$geneid, levels = gene_order)

# Plotting 
ggplot(onion_means, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  facet_grid(~geneid, scales = "free") +  
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  # geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5,
  #           position = position_dodge(width = 1),inherit.aes = TRUE) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + NoLegend() +
  ylab("Fraction of cells") 

# ==============================================================================
# Supplementary Figure 4A - Double positive cells percentage (FURIN)
# ==============================================================================
# Extract out total number of cells per ident
table1 <- as.data.frame(table(ild$orig.ident, ild$CellType1, ild$Diagnosis2))

# ACE2+ FURIN+ cells
sub1 <- subset(ild, ACE2 > 0 & FURIN > 0, slot = "counts") 
tableS1 <- as.data.frame(table(sub1$orig.ident,sub1$CellType1, sub1$Diagnosis2))
data_tableS1 <- merge(table1, tableS1, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS1) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS1$Percent <- data_tableS1$Count/data_tableS1$Total * 100
data_tableS1$geneid <- "ACE2+/FURIN+"

# BSG+ FURIN+ cells
sub2 <- subset(ild, BSG > 0 & FURIN > 0, slot = "counts") 
tableS2 <- as.data.frame(table(sub2$orig.ident,sub2$CellType1, sub2$Diagnosis2))
data_tableS2 <- merge(table1, tableS2, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS2) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS2$Percent <- data_tableS2$Count/data_tableS2$Total * 100
data_tableS2$geneid <- "BSG+/FURIN+"

# HSPA5+ FURIN+ cells
sub3 <- subset(ild, HSPA5 > 0 & FURIN > 0, slot = "counts") 
tableS3 <- as.data.frame(table(sub3$orig.ident,sub3$CellType1, sub3$Diagnosis2))
data_tableS3 <- merge(table1, tableS3, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS3) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS3$Percent <- data_tableS3$Count/data_tableS3$Total * 100
data_tableS3$geneid <- "HSPA5+/FURIN+"

# NRP1+ FURIN+ cells
sub4 <- subset(ild, NRP1 > 0 & FURIN > 0, slot = "counts") 
tableS4 <- as.data.frame(table(sub4$orig.ident, sub4$CellType1, sub4$Diagnosis2))
data_tableS4 <- merge(table1, tableS4, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS4) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS4$Percent <- data_tableS4$Count/data_tableS4$Total * 100
data_tableS4$geneid <- "NRP1+/FURIN+"

# Combine all tables and perform dplyr mean calculations
onion <- rbind(data_tableS1, data_tableS2, data_tableS3, data_tableS4)
onion_means = onion %>% group_by(geneid, CellType, Diagnosis) %>% 
  dplyr::summarise(Mean = mean(Percent, na.rm = T),
                   n=n(),
                   sd = sd(Percent, na.rm = T),
                   se = sd/sqrt(n))

onion_means = as.data.frame(onion_means)
onion_means = onion_means[onion_means$Mean > 0,]

# Add a column for total of cells expressing each gene per CT
ace2_total <- as.data.frame(table(sub1@meta.data$CellType1, sub1$Diagnosis2))
colnames(ace2_total) <- c("CellType","Diagnosis","Total")
ace2_total$geneid <- "ACE2+/FURIN+"

bsg_total <- as.data.frame(table(sub2@meta.data$CellType1, sub2$Diagnosis2))
colnames(bsg_total) <- c("CellType","Diagnosis","Total")
bsg_total$geneid <- "BSG+/FURIN+"

hspa5_total <- as.data.frame(table(sub3@meta.data$CellType1, sub3$Diagnosis2))
colnames(hspa5_total) <- c("CellType","Diagnosis","Total")
hspa5_total$geneid <- "HSPA5+/FURIN+"

nrp1_total <- as.data.frame(table(sub4@meta.data$CellType1, sub4$Diagnosis2))
colnames(nrp1_total) <- c("CellType","Diagnosis","Total")
nrp1_total$geneid <- "NRP1+/FURIN+"

total <- rbind(ace2_total,bsg_total,hspa5_total, nrp1_total)

onion_means <- merge(onion_means, total, by = c("CellType","geneid", "Diagnosis"))

# Tukey test
onion$Diagnosis <- as.factor(onion$Diagnosis)
levels(onion$Diagnosis)

dp_tukey <- onion %>% group_by(CellType,geneid) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis)

# Save the table for supplemental info
write.csv (onion, file = "20210216_FigS4_DP_FURIN_percentage.allCT.csv")
write.csv (onion_means, file = "20210216_FigS4_DP_FURIN_percentage_means.allCT.csv")
write.csv(as.data.frame(dp_tukey), file = "20210216_DP_Epi_FURIN_tukey.csv")

# Make a Venn Diagram for all CT
set1 <- rownames(sub1@meta.data)
set2 <- rownames(sub2@meta.data)
set3 <- rownames(sub3@meta.data)
set4 <- rownames(sub4@meta.data)

myCol <- brewer.pal(4, "Pastel2")

temp <- venn.diagram(
  x = list(set1, set2, set3, set4),
  category.names = c("ACE2+ FURIN+" , "BSG+ FURIN+ " , "HSPA5+ FURIN+","NRP1+ FURIN+"),
  filename = NULL,
  resolution = 500,
  compression = "lzw",
  lwd = 2,
  fill = myCol,
  col = myCol,
  fontfamily = "sans",
  cat.color = myCol,
  cex = 2,
  cat.fontface = 4,
  lty =2,
  height = 1000, 
  width = 500
)

grid.draw(temp)

# Barplot for percentage 
onion_means$CellType <- factor(onion_means$CellType, 
                               levels=cell_order)
onion_means$Diagnosis <- factor(onion_means$Diagnosis,
                                levels=c("Other-ILD","IPF","COPD","Control"))
ggplot(onion_means, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  facet_grid(~geneid, scales = "free") +  
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  #geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5,
  #          position = position_dodge(width = 1),inherit.aes = TRUE) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + 
  ylab("Percentage of cells") 

# ==============================================================================
# Supplementary Figure 4B - Double positive cells percentage (CTSL)
# ==============================================================================
# Extract out total number of cells per ident
table1 <- as.data.frame(table(ild$orig.ident, ild$CellType1, ild$Diagnosis2))

# ACE2+ CTSL+ cells
sub1 <- subset(ild, ACE2 > 0 & CTSL > 0, slot = "counts") 
tableS1 <- as.data.frame(table(sub1$orig.ident,sub1$CellType1, sub1$Diagnosis2))
data_tableS1 <- merge(table1, tableS1, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS1) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS1$Percent <- data_tableS1$Count/data_tableS1$Total * 100
data_tableS1$geneid <- "ACE2+/CTSL+"

# BSG+ CTSL+ cells
sub2 <- subset(ild, BSG > 0 & CTSL > 0, slot = "counts") 
tableS2 <- as.data.frame(table(sub2$orig.ident,sub2$CellType1, sub2$Diagnosis2))
data_tableS2 <- merge(table1, tableS2, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS2) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS2$Percent <- data_tableS2$Count/data_tableS2$Total * 100
data_tableS2$geneid <- "BSG+/CTSL+"

# HSPA5+ CTSL+ cells
sub3 <- subset(ild, HSPA5 > 0 & CTSL > 0, slot = "counts") 
tableS3 <- as.data.frame(table(sub3$orig.ident,sub3$CellType1, sub3$Diagnosis2))
data_tableS3 <- merge(table1, tableS3, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS3) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS3$Percent <- data_tableS3$Count/data_tableS3$Total * 100
data_tableS3$geneid <- "HSPA5+/CTSL+"

# NRP1+ CTSL+ cells
sub4 <- subset(ild, NRP1 > 0 & CTSL > 0, slot = "counts") 
tableS4 <- as.data.frame(table(sub4$orig.ident, sub4$CellType1, sub4$Diagnosis2))
data_tableS4 <- merge(table1, tableS4, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS4) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS4$Percent <- data_tableS4$Count/data_tableS4$Total * 100
data_tableS4$geneid <- "NRP1+/CTSL+"

# Combine all tables and perform dplyr mean calculations
onion <- rbind(data_tableS1, data_tableS2, data_tableS3, data_tableS4)
onion_means = onion %>% group_by(geneid, CellType, Diagnosis) %>% dplyr::summarise(Mean = mean(Percent, na.rm = T),
                                                                                   n=n(),
                                                                                   sd = sd(Percent, na.rm = T),
                                                                                   se = sd/sqrt(n))
onion_means = as.data.frame(onion_means)
onion_means = onion_means[onion_means$Mean > 0,]

# Add a column for total of cells expressing each gene per CT
ace2_total <- as.data.frame(table(sub1@meta.data$CellType1, sub1$Diagnosis2))
colnames(ace2_total) <- c("CellType","Diagnosis","Total")
ace2_total$geneid <- "ACE2+/CTSL+"

bsg_total <- as.data.frame(table(sub2@meta.data$CellType1, sub2$Diagnosis2))
colnames(bsg_total) <- c("CellType","Diagnosis","Total")
bsg_total$geneid <- "BSG+/CTSL+"

hspa5_total <- as.data.frame(table(sub3@meta.data$CellType1, sub3$Diagnosis2))
colnames(hspa5_total) <- c("CellType","Diagnosis","Total")
hspa5_total$geneid <- "HSPA5+/CTSL+"

nrp1_total <- as.data.frame(table(sub4@meta.data$CellType1, sub4$Diagnosis2))
colnames(nrp1_total) <- c("CellType","Diagnosis","Total")
nrp1_total$geneid <- "NRP1+/CTSL+"

total <- rbind(ace2_total,bsg_total,hspa5_total, nrp1_total)

onion_means <- merge(onion_means, total, by = c("CellType","geneid", "Diagnosis"))

# Tukey test
onion$Diagnosis <- as.factor(onion$Diagnosis)
levels(onion$Diagnosis)

dp_tukey <- onion %>% group_by(CellType,geneid) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis)

# Save the table for supplemental info
write.csv (onion, file = "20210216_FigS4_DP_CTSL_percentage.allCT.csv")
write.csv (onion_means, file = "20210216_FigS4_DP_CTSL_percentage_means.allCT.csv")
write.csv(as.data.frame(dp_tukey), file = "20210216_DP_Epi_CTSL_tukey.csv")

# Make a Venn Diagram for all CT
set1 <- rownames(sub1@meta.data)
set2 <- rownames(sub2@meta.data)
set3 <- rownames(sub3@meta.data)
set4 <- rownames(sub4@meta.data)

myCol <- brewer.pal(4, "Pastel2")

temp <- venn.diagram(
  x = list(set1, set2, set3,set4),
  category.names = c("ACE2+ CTSL+" , "BSG+ CTSL+ " , "HSPA5+ CTSL+", "NRP1+ CTSL+"),
  filename = NULL,
  resolution = 500,
  compression = "lzw",
  lwd = 2,
  fill = myCol,
  col = myCol,
  fontfamily = "sans",
  cat.color = myCol,
  cex = 2,
  cat.fontface = 4,
  lty =2,
  height = 1000, 
  width = 500
)

grid.draw(temp)

# Plotting 
onion_means$CellType <- factor(onion_means$CellType, 
                               levels=cell_order)
onion_means$Diagnosis <- factor(onion_means$Diagnosis,
                                 levels=c("Other-ILD","IPF","COPD","Control"))
ggplot(onion_means, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  facet_grid(~geneid, scales = "free") +  
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  #geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5,
  #          position = position_dodge(width = 1),inherit.aes = TRUE) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + 
  ylab("Percentage of cells") +
  NoLegend()

# ==============================================================================
# Figure 3A and Table S10: SARS-CoV-2 gene correlation analysis
# ==============================================================================
epi <- readRDS("/scratch/lbui/Covid19_ILD_objects/20210204_Epithelial_noDoublets.rds")

library(gridExtra)
library(ggplot2)

interleave <- function(a, b) { 
  
  shorter <- if (length(a) < length(b)) a else b
  longer  <- if (length(a) >= length(b)) a else b
  
  slen <- length(shorter)
  llen <- length(longer)
  
  index.short <- (1:slen) + llen
  names(index.short) <- (1:slen)
  
  lindex <- (1:llen) + slen
  names(lindex) <- 1:llen
  
  sindex <- 1:slen
  names(sindex) <- 1:slen
  
  index <- c(sindex, lindex)
  index <- index[order(names(index))]
  
  return(c(a, b)[index])
}

# Separate by disease groups
at2 <- subset(epi, cells = rownames(epi@meta.data[epi@meta.data$CellTypeSimple == "AT2",]))
sub1 <- subset(at2, cells=rownames(at2@meta.data[at2@meta.data$Diagnosis2 == "Control",]))
sub2 <- subset(at2, cells=rownames(at2@meta.data[at2@meta.data$Diagnosis2 == "COPD",]))
sub3 <- subset(at2, cells=rownames(at2@meta.data[at2@meta.data$Diagnosis2 == "IPF",]))
sub4 <- subset(at2, cells=rownames(at2@meta.data[at2@meta.data$Diagnosis2 == "Other-ILD",]))

gene1 <- c("ACE2") #change name correspondingly
gene2 <- c("TMPRSS2","CTSL","FURIN","ADAM17","ITGB6")

plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()
plot_list4 <- list()
plot_index <- 1
for (j in 1:length(gene1)) {
  
  for (i in 1:length(gene2)) {
    onion1 = sub1@assays$SCT@data[gene1[j], ]
    onion2 = sub1@assays$SCT@data[gene2[i], ]
    
    onion1 = as.data.frame(onion1)
    onion2 = as.data.frame(onion2)
    
    onion1$ind = sub1@meta.data[rownames(onion1), ]$orig.ident
    onion2$ind = sub1@meta.data[rownames(onion2), ]$orig.ident
    
    onion1_means = onion1 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion1, na.rm = T))
    onion2_means = onion2 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
    
    onion1_means = as.data.frame(onion1_means)
    onion1_means = onion1_means[onion1_means$Mean > 0,]
    
    onion2_means = as.data.frame(onion2_means)
    onion2_means = onion2_means[onion2_means$Mean > 0,]
    
    means1 = merge(onion1_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
    means1 = means1[,-4]
    
    fit1 <- lm(means1$Mean.y ~ means1$Mean.x)
    pVal1 <- round(anova(fit1)$'Pr(>F)'[1], 4)
    
    if (pVal1 < 1) {
      plot_list1[[plot_index]] <- ggplot(means1, aes(x = Mean.x, y = Mean.y)) +
        geom_point() + 
        ggtitle(paste("AT2_Control", pVal1, sep = "")) + 
        xlab(gene1[j]) + 
        ylab(gene2[i]) + 
        theme(plot.title = element_text(size = 5),
              legend.title = element_blank(),
              axis.title=element_text(size=5), 
              axis.text = element_text(size = 5),
              legend.text = element_text(size = 3)) +
        geom_smooth(method = lm, se = TRUE) + 
        NoLegend()
      
      plot_index <- plot_index + 1
    } 
    
    onion3 = sub2@assays$SCT@data[gene1[j], ]
    onion4 = sub2@assays$SCT@data[gene2[i], ]
    
    onion3 = as.data.frame(onion3)
    onion4 = as.data.frame(onion4)
    
    onion3$ind = sub2@meta.data[rownames(onion3), ]$orig.ident
    onion4$ind = sub2@meta.data[rownames(onion4), ]$orig.ident
    
    onion3_means = onion3 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion3, na.rm = T))
    onion4_means = onion4 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion4, na.rm = T))
    
    onion3_means = as.data.frame(onion3_means)
    onion3_means = onion3_means[onion3_means$Mean > 0,]
    
    onion4_means = as.data.frame(onion4_means)
    onion4_means = onion4_means[onion4_means$Mean > 0,]
    
    means2 = merge(onion3_means, onion4_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
    means2 = means2[,-4]
    
    fit2 <- lm(means2$Mean.y ~ means2$Mean.x)
    pVal2 <- round(anova(fit2)$'Pr(>F)'[1], 4)
    
    if (pVal2 < 1) {
      plot_list2[[plot_index]] <- ggplot(means2, aes(x = Mean.x, y = Mean.y)) +
        geom_point() + 
        ggtitle(paste("AT2_COPD", pVal2, sep = "")) + 
        xlab(gene1[j]) + 
        ylab(gene2[i]) + 
        theme(plot.title = element_text(size = 5),
              legend.title = element_blank(),
              axis.title=element_text(size=5), 
              axis.text = element_text(size = 5),
              legend.text = element_text(size = 3)) +
        geom_smooth(method = lm, se = TRUE) + 
        NoLegend()
      
      plot_index <- plot_index + 1
    } 
    
    onion5 = sub3@assays$SCT@data[gene1[j], ]
    onion6 = sub3@assays$SCT@data[gene2[i], ]
    
    onion5 = as.data.frame(onion5)
    onion6 = as.data.frame(onion6)
    
    onion5$ind = sub3@meta.data[rownames(onion5), ]$orig.ident
    onion6$ind = sub3@meta.data[rownames(onion6), ]$orig.ident
    
    onion5_means = onion5 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion5, na.rm = T))
    onion6_means = onion6 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion6, na.rm = T))
    
    onion5_means = as.data.frame(onion5_means)
    onion5_means = onion5_means[onion5_means$Mean > 0,]
    
    onion6_means = as.data.frame(onion6_means)
    onion6_means = onion6_means[onion6_means$Mean > 0,]
    
    means3 = merge(onion5_means, onion6_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
    means3 = means3[,-4]
    
    fit3 <- lm(means3$Mean.y ~ means3$Mean.x)
    pVal3 <- round(anova(fit3)$'Pr(>F)'[1], 4)
    
    if (pVal3 < 1) {
      plot_list3[[plot_index]] <- ggplot(means3, aes(x = Mean.x, y = Mean.y)) +
        geom_point() + 
        ggtitle(paste("AT2_IPF", pVal3, sep = "")) + 
        xlab(gene1[j]) + 
        ylab(gene2[i]) + 
        theme(plot.title = element_text(size = 5),
              legend.title = element_blank(),
              axis.title=element_text(size=5), 
              axis.text = element_text(size = 5),
              legend.text = element_text(size = 3)) +
        geom_smooth(method = lm, se = TRUE) + 
        NoLegend()
      
      plot_index <- plot_index + 1
    } 
    
    onion7 = sub4@assays$SCT@data[gene1[j], ]
    onion8 = sub4@assays$SCT@data[gene2[i], ]
    
    onion7 = as.data.frame(onion7)
    onion8 = as.data.frame(onion8)
    
    onion7$ind = sub4@meta.data[rownames(onion7), ]$orig.ident
    onion8$ind = sub4@meta.data[rownames(onion8), ]$orig.ident
    
    onion7_means = onion7 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion7, na.rm = T))
    onion8_means = onion8 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion8, na.rm = T))
    
    onion7_means = as.data.frame(onion7_means)
    onion7_means = onion7_means[onion7_means$Mean > 0,]
    
    onion8_means = as.data.frame(onion8_means)
    onion8_means = onion8_means[onion8_means$Mean > 0,]
    
    means4 = merge(onion7_means, onion8_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
    means4 = means4[,-4]
    
    fit4 <- lm(means4$Mean.y ~ means4$Mean.x)
    pVal4 <- round(anova(fit4)$'Pr(>F)'[1], 4)
    
    if (pVal4 < 1) {
      plot_list4[[plot_index]] <- ggplot(means4, aes(x = Mean.x, y = Mean.y)) +
        geom_point() + 
        ggtitle(paste("AT2_OtherILD", pVal4, sep = "")) + 
        xlab(gene1[j]) + 
        ylab(gene2[i]) + 
        theme(plot.title = element_text(size = 5),
              legend.title = element_blank(),
              axis.title=element_text(size=5), 
              axis.text = element_text(size = 5),
              legend.text = element_text(size = 3)) +
        geom_smooth(method = lm, se = TRUE) + 
        NoLegend()
      
      plot_index <- plot_index + 1
    } 
  }
}

plot_lista <- interleave(plot_list1, plot_list2)
plot_listb <- interleave(plot_list3, plot_list4)
plot_list <- interleave(plot_lista, plot_listb)

plot_list <- plyr::compact(plot_list)

pdf(file = paste("20210208_Correlation_plot_SCT_AT2_diseasegroups1.pdf", 
                 sep = "_"), width = 8.5, height = 11)
grid.arrange(grobs = plot_list, ncol = 4)
dev.off()

rm(fit1,fit2,fit3,fit4,gene1,gene2,i,j,means1,means2,means3,means4,onion1,onion1_means,onion2,onion2_means,onion3,onion3_means,onion4,onion4_means,onion5,onion5_means,onion6,onion6_means,onion7,onion7_means,onion8,onion8_means,plot_index,plot_list,plot_list1,plot_list2,plot_list3,plot_list4,plot_lista,plot_listb,pVal1,pVal2,pVal3,pVal4)

# ==============================================================================
# Supplementary Dataset 1: Count Matrix for unpublished dataset
# ==============================================================================
# Load the object
ild <- readRDS("/scratch/lbui/Covid19_ILD_objects/20210211_ILD_noDoublets.rds")

# Extract out gene count matrix for genes used in the study
genelist <- c("ACE2","BSG","TMPRSS2","NPR1","CTSL","CTSB","FURIN","PCSK5","PCSK7",
              "ADAM17","PIKFYVE","TPCN2","AGT","ACE","ITGB6","IFNAR1","IFNAR2",
              "EIF2AK2","EIF2AK3","CD44","IFNGR1","IFNGR2","FAM46C","UBD","REC8",
              "ELF1","CLEC4D","LY6E","SPATS2L","ZBP1","DNAJC6","IFIT3","RGS22",
              "B4GALT5","ISG20","GNB4","SPATA13","NRN1","ERLIN1","APOL2","RAB27A",
              "FZD5","C9orf91","TAGAP","HSPA8","CNP","ETV6","MSR1","BST2","CXCL1",
              "CSF2","CSF3","ICAM1","CD47","CCL2","CCL3","TRIM27","TRIM28","RNF41",
              "NFKB1","JUN","SOCS1","SOCS2","C1QA","C1QB","C1QC","C2","C3","C4B",
              "PTPN11","FGA","FGG","ATG5","ATG7","BECN1","SQSTM1","CXCL10", 
              "MAP1LC3A","MAP1LC3B","ATF6","ERN1","MUC5B","ISG15", "IFI44", "IFI27", 
              "RSAD2", "IFIT1", "IFI44L", "CCL8", "XAF1", "GBP1", "IRF7", 
              "CEACAM1","IFNB1","IFNG","IFNAR1","IFNAR2","IFNGR1","IFNGR2","TRIM28",
              "TLR7","MX1","STAT1","TBK1","CCR2","CXCL10","IFI6","LY6E",
              "TRIM27","TNF","IL1B","IL6","IL6R","IL6ST","TGFB1","NFKB1",
              "NFKB2","CEBPB","AREG","FCGR3A","IL10","IFITM1","IFITM3","ISG20",
              "CD163","IL1R2","MRC1","HAVCR2","LGALS9","S100A8","S100A9",
              "HLA-DRA", "HLA-DQA1", "HLA-DQA2", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", 
              "HLA-DQB2", "HLA-DRB5", "HLA-DQB1", "HLA-DMA", "HLA-DMB")
genelist <- unique(genelist)

# Subset out only the unpublished samples
samples <- c("TILD041" ,"TILD051" ,"VUILD79" ,"VUILD77", "VUHD84" , "TILD049",
             "VUILD74","VUILD75","VUILD76","VUILD78","VUILD81", "THD0006","TILD013",
             "TILD039","TILD037","VUHD072","VUILD67","VUILD68","VUHD073","VUILD69",
             "VUHD075","VUHD078","VUHD080","THD0007","THD0008","THD0009","TILD050",
             "TILD055","TILD062","VUILD73","TILD063","TILD059","VUHD85","THD0012",
             "VUILD80" ,"VU_COPD_29","VU_COPD_30","VU_COPD_34","VU_COPD_38")

# Subset out only those genes
ild.df <- subset(ild, features = genelist, 
                 cells = rownames(ild@meta.data[ild@meta.data$Sample_Name %in% samples, ]))

saveRDS(ild.df, file = "20210216_ILD_COVID_unpublished39samples.rds")

# Export out the count matrix
write.table(as.matrix(GetAssayData(object = ild.df, slot = "counts")), 
            file = "TableS4_ILD_COVID_unpublished39samples_counts.csv", 
            sep = ',', row.names = T, col.names = T, quote = F)

# ==============================================================================
# Supplementary Figure 5: SARS-CoV-2 entry gene expression
# ==============================================================================
epi <- readRDS("/scratch/lbui/Covid19_saved/20210204_Epithelial_noDoublets.rds")
vlncol <- c("Control"="#F8766D","COPD"="#7CAE00","IPF"="#00BFC4",
            "Other-ILD"="#C77CFF")
VlnPlot(epi, c("NRP1"), group.by = "CellType2",split.by = "Diagnosis2", 
        pt.size = 0, assay="SCT", cols=vlncol) + NoLegend()
DotPlot(epi, features = "ACE2", group.by = "CellType2", split.by = "Status")
VlnPlot(epi, features = "ACE2", group.by = "CellType2", slot = "counts")

# ==============================================================================
# Supplementary Figure 6: NRP1
# ==============================================================================
VlnPlot(ild, "NRP1", group.by = "CellType2", pt.size = 0) + NoLegend()
celltype.use <- c("Vascular Endothelial Cells","Macrophages",
                  "Proliferating Macrophages","pDCs","Fibroblasts","Myofibroblasts",
                  "HAS1 High Fibroblasts","PLIN2+ Fibroblasts")
immune <- readRDS("/scratch/lbui/Covid19_saved/20210204_Immune_noDoublets.rds")
endo <- readRDS("/scratch/lbui/Covid19_saved/20210204_Endothelial_noDoublets.rds")
mesen <- readRDS("/scratch/lbui/Covid19_saved/20210204_Mesenchymal_noDoublets.rds")

nrp1_list = list()
j <- 0
for(i in unique(mesen@meta.data$CellType2)){
  j <- j + 1
  nrp1_list[[j]] <- subset(mesen,
                          cells = row.names(mesen@meta.data[mesen@meta.data$CellType2 == i, ]))
}
for(i in 1:length(nrp1_list)){
  names(nrp1_list) <- lapply(nrp1_list, function(xx){paste(unique(xx@meta.data$CellType2))})
}
summary(nrp1_list)
macrophage <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$CellType2 == "Macrophages",]))
pro.macro <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$CellType2 == "Proliferating Macrophages",]))
pdc <- subset(immune, cells=rownames(immune@meta.data[immune@meta.data$CellType2 == "pDCs",]))
vas.endo <- subset(endo, cells=rownames(endo@meta.data[endo@meta.data$CellType2 == "Vascular Endothelial Cells",]))

nrp1_list <- c(nrp1_list, macrophage, pro.macro, pdc, vas.endo)
names(nrp1_list) <- c("Mesothelial","Fibroblasts","Smooth Muscle Cells","Pericytes",
                      "Myofibroblasts","PLIN2+ Fibroblasts","HAS1 High Fibroblasts",
                      "Macrophages","Proliferating Macrophages", "pDCs",
                      "Vascular Endothelial Cells")
nrp1_list <- nrp1_list[names(nrp1_list) %in% celltype.use]

# Calculate p_adj_value using FindMarkers neg binom test
copd_vs_control <- lapply(nrp1_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "COPD",
                ident.2 = "Control",
                test.use = "negbinom", # using this test: the slot is set to "counts"
                features = "NRP1",
                latent.vars = c("dataset","Age","Ethinicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

ipf_vs_control <- lapply(nrp1_list, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "IPF",
                ident.2 = "Control",
                test.use = "negbinom",
                features = "NRP1",
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})

nrp1_list2 <- nrp1_list[-4] # error with HAS since there's no DEG found
other_vs_control <- lapply(nrp1_list2, function(xx){
  print(unique(xx@meta.data$CellType2))
  if(length(unique(xx@meta.data$Diagnosis2)) > 1) {
    FindMarkers(xx,
                group.by = "Diagnosis2",
                ident.1 = "Other-ILD",
                ident.2 = "Control",
                test.use = "negbinom", 
                features = "NRP1",
                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                logfc.threshold = 0,
                assay="SCT")
  } 
  else{
    return(NULL)
  } 
})
for(i in 1:length(nrp1_list)){
  write.table(copd_vs_control[[i]], paste(gsub("/", "", unique(nrp1_list[[i]]@meta.data$CellType2)), "_copd_vs_control", ".csv"), sep =",", quote = F)
}
for(i in 1:length(nrp1_list)){
  write.table(ipf_vs_control[[i]], paste(gsub("/", "", unique(nrp1_list[[i]]@meta.data$CellType2)), "_ipf_vs_control", ".csv"), sep =",", quote = F)
}
for(i in 1:length(nrp1_list2)){
  write.table(other_vs_control[[i]], paste(gsub("/", "", unique(nrp1_list2[[i]]@meta.data$CellType2)), "_other_vs_control", ".csv"), sep =",", quote = F)
}

# ==============================================================================
# Supplementary Figure 8E: Violin plot for Epithelial cells
# ==============================================================================
# Read in Epithelial cells
epi <- readRDS("/scratch/lbui/Covid19_ILD_objects/20210204_Epithelial_noDoublets.rds")

# Violin plots
vlncol <- c("Control"="#F8766D","COPD"="#7CAE00","IPF"="#00BFC4",
            "Other-ILD"="#C77CFF")
genes <- c("CSF3","ITGB6","SOCS1","SOCS2","LY6E")
plot_list <- list()
j=0
for(i in(1:length(genes))){
  j=j+1
  plot_list[[j]] <- VlnPlot(epi, features=genes[i], group.by = "CellType2", split.by = "Diagnosis2",
          split.plot = T, pt.size = 0, cols = vlncol) + NoLegend()
}
pdf("20210217_Fig8_vlnplot.pdf")
plot_list
dev.off()

# ==============================================================================
# Supplementary Figure 8F: Volcano plot for DEG in AT2 CLD vs. Control
# ==============================================================================
# Subset out AT2 cells
at2 <- subset(epi, cells=rownames(epi@meta.data[epi@meta.data$CellTypeSimple == "AT2",]))

# Run DEG
copd_control_at2 <- FindMarkers(at2, 
                                ident.1 = "COPD", 
                                ident.2 = "Control",
                                group.by = "Diagnosis2",
                                test.use = "negbinom",
                                latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                                logfc.threshold = 0,
                                assay="SCT")

ipf_control_at2 <- FindMarkers(at2, 
                               ident.1 = "IPF", 
                               ident.2 = "Control",
                               group.by = "Diagnosis2",
                               test.use = "negbinom",
                               latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                               logfc.threshold = 0,
                               assay="SCT")

other_control_at2 <- FindMarkers(at2, 
                                 ident.1 = "Other-ILD", 
                                 ident.2 = "Control",
                                 group.by = "Diagnosis2",
                                 test.use = "negbinom",
                                 latent.vars = c("dataset","Age","Ethnicity","Smoking_status"),
                                 logfc.threshold = 0,
                                 assay="SCT")

# Save objects
write.csv(copd_control_at2, file = "20210215_AT2_COPD_Control_DEG.csv")
write.csv(ipf_control_at2, file = "20210215_AT2_IPF_Control_DEG.csv")
write.csv(other_control_at2, file = "20210215_AT2_Other_Control_DEG.csv")

# Make a volcano plot
copd_control_at2$genename <- rownames(copd_control_at2)
ipf_control_at2$genename <- rownames(ipf_control_at2)
other_control_at2$genename <- rownames(other_control_at2)

copd_control_at2$Test <- "COPD"
ipf_control_at2$Test <- "IPF"
other_control_at2$Test <- "Other-ILD"

test <- rbind(copd_control_at2, ipf_control_at2, other_control_at2)
temp1 <- grep( "^MT-", rownames(test), ignore.case = F, value = T) # pull out mitochondria genes
temp2 <- grep( "^RP", rownames(test), ignore.case = F, value = T) # pull out ribosomal genes
test <- test[!test$genename %in% temp2,] # remove ribosomal genes 
test <- test[!test$genename %in% temp1,] # remove mitochondria genes 
copd_control_at2 <- copd_control_at2[!copd_control_at2$genename %in% c(temp1,temp2),]
ipf_control_at2 <- ipf_control_at2[!ipf_control_at2$genename %in% c(temp1,temp2),]
other_control_at2 <- other_control_at2[!other_control_at2$genename %in% c(temp1,temp2),]
cols <- c("all" = "light grey","IPF" = "#00BFC4", "COPD" = "#7CAE00", 
          "Other-ILD" = "#C77CFF")

# Draw volcano plot
ggplot(test, aes(x = avg_log2FC, y = -log10(p_val_adj), color = "all")) + 
  geom_point() +
  geom_point(data = test[abs(test$avg_log2FC) >= 1 & test$Test == "IPF" & 
                           test$p_val_adj <= 0.01,], color = "#00BFC4") +
  geom_point(data = test[abs(test$avg_log2FC) >= 1 & test$Test == "COPD" & 
                           test$p_val_adj <= 0.01,], color = "#7CAE00") +
  geom_point(data = test[abs(test$avg_log2FC) >= 1 & test$Test == "Other-ILD" & 
                           test$p_val_adj <= 0.01,], color = "#C77CFF") +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("P value, Log"[10]*"")) +
  theme_bw() +
  geom_text_repel(data=subset(ipf_control_at2,abs(avg_log2FC) >= 1 & p_val_adj <= 0.01),
                  aes(avg_log2FC, -log10(p_val), label = genename, colour="IPF"),
                  size = 3, nudge_x = -0.05, nudge_y = 0, max.overlaps = 100) +
  geom_text_repel(data=subset(copd_control_at2,abs(avg_log2FC) >= 1 & p_val_adj <= 0.01),
                  aes(avg_log2FC, -log10(p_val), label = genename, colour="COPD"),
                  size = 3, nudge_x = 0, nudge_y = 0, max.overlaps = 100) +
  geom_text_repel(data=subset(other_control_at2,abs(avg_log2FC) >= 1 & p_val_adj <= 0.01),
                  aes(avg_log2FC, -log10(p_val), label = genename, colour="Other-ILD"),
                  size = 3, nudge_x = 0.2, nudge_y = -1, max.overlaps = 100) +
  scale_colour_manual(name="Diagnosis",values=cols)  +
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  ggtitle("AT2 cells CLD vs. Control")

# ==============================================================================
# Supplementary Figure 12: Violin plot for all Immune cell types (same as Fig4C)
# ==============================================================================
immune <- readRDS("/scratch/lbui/Covid19_ILD_objects/20210204_Immune_noDoublets.rds")

vlncol <- c("Control"="#F8766D","COPD"="#7CAE00","IPF"="#00BFC4",
            "Other-ILD"="#C77CFF")
VlnPlot(immune, c("S100A9"), group.by = "CellType2",split.by = "Diagnosis2", 
        pt.size = 0, assay="SCT", cols=vlncol) + NoLegend()
