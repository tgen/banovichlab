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
ild <- readRDS("/scratch/lbui/Covid19_Seurat/20200709_ILD_noDoublets.rds")

# Order cell types along the x axis according to population 
cell_order <- c("Lymphatic Endothelial Cells","Vascular Endothelial Cells","AT1",
                "AT2","Basal","Ciliated Cells","Ionocytes","KRT5-/KRT17+","PNEC/Ionocytes",
                "Secretory Cells", "Transitional AT2","B Cells","cDCs","Macrophages",
                "Mast Cells","Monocytes","NK Cells","pDCs","Plasma Cells","T Cells",
                "Fibroblasts","Mesothelial","Pericytes","Smooth Muscle Cells")

ild$CellTypeSimple <- factor(ild$CellTypeSimple, levels=cell_order)

# Make dotplot
DotPlot(ild, features = c("TAGLN","ACTG2","GJA4","RGS5","HP","WT1","LUM","PDGFRA",
                          "CD3E","JCHAIN","IGHG1","LILRA4","CLEC4C","NKG7","NCR1",
                          "S100A12","FCN1","CPA3","KIT","MARCO","APOC1","CLEC9A", 
                          "CD1C","MS4A1","CD19","MUC5B","SCGB3A2","SCGB1A1",
                          "FOXI1","FOXJ1","KRT5","KRT17","SFTPC","ABCA3","AGER",
                          "VWF","CCL21"), 
        group.by = "CellTypeSimple") + 
  theme(axis.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# ==============================================================================
# Supplementary Figure 2 - BSG+, CTSL+ and HSPA5+ cell proportion
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
write.csv (onion, file = "20201123_FigS2_percentage.allCT.csv")
write.csv (onion_means, file = "20201123_FigS2_percentage_means.allCT.csv")
write.csv (as.data.frame(sp_tukey), file = "20201123_FigS2_tukey.allCT.csv")

# Order CT on the y axis
cell_order <- c("Lymphatic Endothelial Cells","Vascular Endothelial Cells","AT1",
                "AT2","Transitional AT2","KRT5-/KRT17+","Basal","Ciliated Cells",
                "Differentiating Ciliated Cells","PNEC/Ionocytes",
                "Club Cells", "Goblet Cells","B Cells","cDCs","Macrophages",
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
# Supplementary Figure 3A - Double positive cells percentage (FURIN)
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
write.csv (onion, file = "20201124_FigS3_DP_FURIN_percentage.allCT.csv")
write.csv (onion_means, file = "20201124_FigS3_DP_FURIN_percentage_means.allCT.csv")
write.csv(as.data.frame(dp_tukey), file = "20201124_DP_Epi_FURIN_tukey.csv")

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
# Supplementary Figure 3B - Double positive cells percentage (CTSL)
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
write.csv (onion, file = "20201124_FigS3_DP_CTSL_percentage.allCT.csv")
write.csv (onion_means, file = "20201124_FigS3_DP_CTSL_percentage_means.allCT.csv")
write.csv(as.data.frame(dp_tukey), file = "20201124_DP_Epi_CTSL_tukey.csv")

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
# Figure 2C and Figure S4: SARS-CoV-2 gene correlation analysis
# ==============================================================================
epi <- readRDS("/scratch/lbui/Covid19_Seurat/20200708_Epi_annotated_noDoublets.rds")

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

gene1 <- c("NRP1") #change name correspondingly
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

pdf(file = paste("20201123_Correlation_plot_SCT_AT2_diseasegroups3.pdf", 
                 sep = "_"), width = 8.5, height = 11)
grid.arrange(grobs = plot_list, ncol = 4)
dev.off()

# ==============================================================================
# Supplementary Table 4: Count Matrix for unpublished dataset
# ==============================================================================
# Load the object
ild <- readRDS("/scratch/lbui/Covid19_Seurat/20200709_ILD_noDoublets.rds")

# Extract out gene count matrix for genes used in the study
genelist <- c("ACE2","BSG","TMPRSS2","CTSL","CTSB","FURIN","PCSK5","PCSK7",
              "ADAM17","PIKFYVE","TPCN2","AGT","ACE","ITGB6","C1QA","C1QB",
              "C1QC","C2","C3","C4B","IFNAR1","IFNAR2","IFNGR1","NRP1","HSPA5",
              "IFNGR2","PTPN11","EIF2AK2","EIF2AK3","CXCL1","TRIM27","TRIM28","NFKB1",
              "RNF41","JUN","SOCS1","SOCS2","CSF2","CSF3","ICAM1","CD47",
              "CD44","CCL2","CCL3","FGA","FGG","ATG5","ATG7","BECN1","SQSTM1",
              "MAP1LC3A","MAP1LC3B","ATF6","ERN1","MUC5B","ISG15", "IFI44", "IFI27", 
              "CXCL10", "RSAD2", "IFIT1", "IFI44L", "CCL8", "XAF1", "GBP1", "IRF7", 
              "CEACAM1","IFNB1","IFNG","IFNAR1","IFNAR2","IFNGR1","IFNGR2","TRIM28",
              "TLR7","MX1","STAT1","TBK1","CCR2","CXCL10","IFI6","LY6E",
              "TRIM27","TNF","IL1B","IL6","IL6R","IL6ST","TGFB1","NFKB1",
              "NFKB2","CEBPB","AREG","FCGR3A","IL10","IFITM1","IFITM3","ISG20",
              "CD163","IL1R2","MRC1","HAVCR2","LGALS9","S100A8","S100A9",
              "HLA-DRA", "HLA-DQA1", "HLA-DQA2", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", 
              "HLA-DQB2", "HLA-DRB5", "HLA-DQB1", "HLA-DMA", "HLA-DMB","IL6R","IL6ST",
              "SOCS1","SOCS2","CCL2","CCL3","CXCL8","AREG","FCGR3A",
              "S100A8","S100A9","GZMB","LAG3","ISG15","IFI44","IFI27","CXCL10",
              "RSAD2","IFIT1","IFI44L","CCL8","XAF1","GBP1","IRF7","CEACAM1",
              "ARID5A","SOCS3","PIM1","BCL3","BATF","MYC","PRF1", "GZMH", "IFNG",
              "NKG7", "KLRG1", "PRF1","GNLY","GZMB","GZMK","LAG3","TIGIT","PDCD1",
              "CTLA4","HAVCR2","TOX","PRDM1","MAF")
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

saveRDS(ild.df, file = "20201201_ILD_COVID_unpublished39samples.rds")

# Export out the count matrix
write.table(as.matrix(GetAssayData(object = ild.df, slot = "counts")), 
            file = "TableS4_ILD_COVID_unpublished39samples_counts.csv", 
            sep = ',', row.names = T, col.names = T, quote = F)

