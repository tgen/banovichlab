# ====================================
# Author: Linh T. Bui, lbui@tgen.org
# Date: 2020-06-29
# Title: scRNA-seq data analysis for Covid-19 manuscript - Figure 1
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
library(grid)
library(rstatix)
library(VennDiagram)
library(RColorBrewer)

# ==============================================================================
# Read in the ILD object (containing all cell types)
# ==============================================================================
ild <- readRDS("/scratch/lbui/Covid19_Seurat/20200709_ILD_noDoublets.rds")

# ==============================================================================
# FIGURE 1A: ACE2+ and TMPRSS2+ cells
# ==============================================================================
# Extract out total number of cells per ident
table1 <- as.data.frame(table(ild$orig.ident,ild$CellType1))

# ACE2+ cells
sub1 <- subset(ild, ACE2 > 0, slot = "counts") 
tableA1 <- as.data.frame(table(sub1$orig.ident,sub1$CellType1))
data_tableA1 <- merge(table1, tableA1, by = c("Var1", "Var2"))
colnames(data_tableA1) <- c("Ident","CellType","Total","Count")
data_tableA1$Percent <- data_tableA1$Count/data_tableA1$Total * 100
data_tableA1$geneid <- "ACE2"

# TMPRSS2+ cells
sub2 <- subset(ild, TMPRSS2 > 0, slot = "counts")
tableA2 <- as.data.frame(table(sub2$orig.ident,sub2$CellType1))
data_tableA2 <- merge(table1, tableA2, by = c("Var1", "Var2"))
colnames(data_tableA2) <- c("Ident","CellType","Total","Count")
data_tableA2$Percent <- data_tableA2$Count/data_tableA2$Total * 100
data_tableA2$geneid <- "TMPRSS2"

# Combine all tables and perform dplyr mean calculations
onion <- rbind(data_tableA1, data_tableA2)
onion_means = onion %>% group_by(CellType, geneid) %>% 
  dplyr::summarise(Mean = mean(Percent, na.rm = T),
                   n=n(),
                   sd = sd(Percent, na.rm = T),
                   se = sd/sqrt(n))
onion_means = as.data.frame(onion_means)
onion_means = onion_means[onion_means$Mean > 0,]

# Add a column for total of cells expressing each gene per CT
ace2_total <- as.data.frame(table(sub1@meta.data$CellType1))
colnames(ace2_total) <- c("CellType","Total")
ace2_total$geneid <- "ACE2"

tmprss2_total <- as.data.frame(table(sub2@meta.data$CellType1))
colnames(tmprss2_total) <- c("CellType","Total")
tmprss2_total$geneid <- "TMPRSS2"

total <- rbind(ace2_total, tmprss2_total)

onion_means <- merge(onion_means, total, by = c("CellType","geneid"))

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

# Plotting
ggplot(onion_means, aes(x=CellType, y= Mean, fill = CellType)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  facet_grid(~geneid, scales = "free") +  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(label = paste(Total,"(",round(Mean,2),")",sep = "")), hjust = -1, size = 3) +
  theme(strip.text.x = element_text(size = 8, colour = "red", angle = 0)) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + NoLegend() +
  ylab("Percentage of cells") 

# Save all files
write.csv(onion, file = "20200713_Fig1A_AllCT_percentage.csv")
write.csv(onion_means, file = "20200713_Fig1A_AllCT_percentage_means.csv")

# ==============================================================================
# FIGURE 1B-1C: ACE2+ AND TMPRSS2+ CELLS
# ==============================================================================
# -------------------------------------
# Figure 1B: ACE2 plot
# -------------------------------------
# Select cell types with > 1% cells expressing ACE2 for plotting
ace2_CT <- onion_means[onion_means$geneid == "ACE2" & onion_means$Mean >= 0.5,]$CellType
sub1 <- subset(sub1, cells = rownames(sub1@meta.data[sub1@meta.data$CellType1 %in% ace2_CT,])) #ACE2+ cells
sub <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$CellType1 %in% ace2_CT,])) # All cells

# Extract out the ACE2+ counts per diagnosis
tableB1 <- as.data.frame(table(sub$orig.ident,sub$CellType2,sub$Diagnosis2))
tableB2 <- as.data.frame(table(sub1$orig.ident,sub1$CellType2,sub1$Diagnosis2))
data_tableB1 <- merge(data.frame(tableB1, row.names=NULL), data.frame(tableB2, row.names=NULL), 
                      by = c("Var1", "Var2", "Var3"), all = TRUE)
colnames(data_tableB1) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableB1$Percent <- data_tableB1$Count/data_tableB1$Total * 100

# Calculate means per CT group
onion_meansB1 = data_tableB1 %>% group_by(CellType, Diagnosis) %>% dplyr::summarise(Mean = mean(Percent, na.rm = T),
                                                                           n=n(),
                                                                           sd = sd(Percent, na.rm = T),
                                                                           se = sd/sqrt(n))
onion_meansB1 = as.data.frame(onion_meansB1)
onion_meansB1 = onion_meansB1[onion_meansB1$Mean > 0,]

# Add total counts of ACE2+ cells per diagnosis into the plot data
ace2_totalB <- as.data.frame(table(sub1@meta.data$CellType2, sub1@meta.data$Diagnosis2))
colnames(ace2_totalB) <- c("CellType","Diagnosis","Total")
ace2_totalB$geneid <- "ACE2"
onion_meansB1 <- merge(onion_meansB1, ace2_totalB, by = c("CellType","Diagnosis"))

# Tukey test (updated since this is better with multitesting)
data_tableB1$Diagnosis <- as.factor(data_tableB1$Diagnosis)
levels(data_tableB1$Diagnosis)

ace2_tukey <- data_tableB1 %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis)

# Save the table for supplemental info
write.csv (data_tableB1, file = "20201013_Fig1B_ACE2_percentage.allCT.csv")
write.csv (onion_meansB1, file = "20201013_Fig1B_ACE2_percentage_means.allCT.csv")
write.csv(as.data.frame(ace2_tukey), file = "20201013_ACE2_Epi_tukey.csv")

# Plotting 
onion_meansB1$Diagnosis <- factor(onion_meansB1$Diagnosis,
                                  levels=c("Other-ILD","IPF","COPD","Control"))
ggplot(onion_meansB1, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5,
            position = position_dodge(width = 1),inherit.aes = TRUE) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + #NoLegend() +
  ylab("Percentage of ACE2+ cells") 

# -------------------------------------
# Figure 1C: TMPRSS2 plot 
# -------------------------------------
# Select cell types with > 10% cells expressing ACE2 for plotting
tmprss2_CT <- onion_means[onion_means$geneid == "TMPRSS2" & round(onion_means$Mean,1) >= 10,]$CellType
sub2 <- subset(sub2, cells = rownames(sub2@meta.data[sub2@meta.data$CellType1 %in% tmprss2_CT,])) #TMPRSS2+ cells
sub <- subset(ild, cells = rownames(ild@meta.data[ild@meta.data$CellType1 %in% tmprss2_CT,])) # All cells

# Extract out the TMPRSS2+ counts per diagnosis
tableB1 <- as.data.frame(table(sub$orig.ident,sub$CellType2,sub$Diagnosis2))
tableB2 <- as.data.frame(table(sub2$orig.ident,sub2$CellType2,sub2$Diagnosis2))
data_tableB2 <-merge(data.frame(tableB1, row.names=NULL), data.frame(tableB2, row.names=NULL), 
                     by = c("Var1", "Var2", "Var3"), all = TRUE)
colnames(data_tableB2) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableB2$Percent <- data_tableB2$Count/data_tableB2$Total * 100

# Calculate means per CT group
onion_meansB2 = data_tableB2 %>% group_by(CellType, Diagnosis) %>% dplyr::summarise(Mean = mean(Percent, na.rm = T),
                                                                                    n=n(),
                                                                                    sd = sd(Percent, na.rm = T),
                                                                                    se = sd/sqrt(n))
onion_meansB2 = as.data.frame(onion_meansB2)
onion_meansB2 = onion_meansB2[onion_meansB2$Mean > 0,]

# Add total counts of TMPRSS2+ cells per diagnosis into the plot data
TMPRSS2_totalB <- as.data.frame(table(sub2@meta.data$CellType2, sub2@meta.data$Diagnosis2))
colnames(TMPRSS2_totalB) <- c("CellType","Diagnosis","Total")
TMPRSS2_totalB$geneid <- "TMPRSS2"
onion_meansB2 <- merge(onion_meansB2, TMPRSS2_totalB, by = c("CellType","Diagnosis"))

# Tukey statistic test
data_tableB2$Diagnosis <- as.factor(data_tableB2$Diagnosis)
levels(data_tableB2$Diagnosis)

tmprss2_tukey <- data_tableB2 %>% group_by(CellType) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis)

# Save the table for supplemental info
write.csv (data_tableB2, file = "20201013_Fig1C_TMPRSS2_percentage.allCT.csv")
write.csv (onion_meansB2, file = "20201013_Fig1C_TMPRSS2_percentage_means.allCT.csv")
write.csv(as.data.frame(tmprss2_tukey), file = "20201013_TMPRSS2_Epi_tukey.csv")

# Plotting 
onion_meansB2$Diagnosis <- factor(onion_meansB2$Diagnosis,
                                  levels=c("Other-ILD","IPF","COPD","Control"))
ggplot(onion_meansB2, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  # scale_fill_manual(values=CT_col) + 
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  # facet_grid(~Diagnosis, scales = "free") +  
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5,
            position = position_dodge(width = 1),inherit.aes = TRUE) +
  # theme(strip.text.x = element_text(size = 8, colour = "red", angle = 0)) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + NoLegend() +
  ylab("Percentage of TMPRSS2+ cells") 

# ==============================================================================
# FIGURE 1D, 1E, 1F - DOUBLE POSITIVE CELL PERCENTAGE
# ==============================================================================
# Extract out total number of cells per ident
table1 <- as.data.frame(table(ild$orig.ident, ild$CellType2, ild$Diagnosis2))

# ACE2+ TMPRSS2+ cells
sub1 <- subset(ild, ACE2 > 0 & TMPRSS2 > 0, slot = "counts") 
tableS1 <- as.data.frame(table(sub1$orig.ident,sub1$CellType2, sub1$Diagnosis2))
data_tableS1 <- merge(table1, tableS1, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS1) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS1$Percent <- data_tableS1$Count/data_tableS1$Total * 100
data_tableS1$geneid <- "ACE2+/TMPRSS2+"

# BSG+ TMPRSS2+ cells
sub2 <- subset(ild, BSG > 0 & TMPRSS2 > 0, slot = "counts") 
tableS2 <- as.data.frame(table(sub2$orig.ident,sub2$CellType2, sub2$Diagnosis2))
data_tableS2 <- merge(table1, tableS2, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS2) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS2$Percent <- data_tableS2$Count/data_tableS2$Total * 100
data_tableS2$geneid <- "BSG+/TMPRSS2+"

# HSPA5+ TMPRSS2+ cells
sub3 <- subset(ild, HSPA5 > 0 & TMPRSS2 > 0, slot = "counts") 
tableS3 <- as.data.frame(table(sub3$orig.ident,sub3$CellType2, sub3$Diagnosis2))
data_tableS3 <- merge(table1, tableS3, by = c("Var1", "Var2", "Var3"))
colnames(data_tableS3) <- c("Ident","CellType","Diagnosis","Total","Count")
data_tableS3$Percent <- data_tableS3$Count/data_tableS3$Total * 100
data_tableS3$geneid <- "HSPA5+/TMPRSS2+"

# Combine all tables and perform dplyr mean calculations
onion <- rbind(data_tableS1, data_tableS2, data_tableS3)
onion_means = onion %>% group_by(geneid, CellType, Diagnosis) %>% dplyr::summarise(Mean = mean(Percent, na.rm = T),
                                                                                   n=n(),
                                                                                   sd = sd(Percent, na.rm = T),
                                                                                   se = sd/sqrt(n))
onion_means = as.data.frame(onion_means)
onion_means = onion_means[onion_means$Mean > 0,]

# Add a column for total of cells expressing each gene per CT
ace2_total <- as.data.frame(table(sub1@meta.data$CellType2, sub1$Diagnosis2))
colnames(ace2_total) <- c("CellType","Diagnosis","Total")
ace2_total$geneid <- "ACE2+/TMPRSS2+"

bsg_total <- as.data.frame(table(sub2@meta.data$CellType2, sub2$Diagnosis2))
colnames(bsg_total) <- c("CellType","Diagnosis","Total")
bsg_total$geneid <- "BSG+/TMPRSS2+"

hspa5_total <- as.data.frame(table(sub3@meta.data$CellType2, sub3$Diagnosis2))
colnames(hspa5_total) <- c("CellType","Diagnosis","Total")
hspa5_total$geneid <- "HSPA5+/TMPRSS2+"

total <- rbind(ace2_total,bsg_total,hspa5_total)

onion_means <- merge(onion_means, total, by = c("CellType","geneid", "Diagnosis"))

# -------------------------------
# Tukey test for significant differences in DP cell percentage per diagnosis group
# -------------------------------
onion$Diagnosis <- as.factor(onion$Diagnosis)
levels(onion$Diagnosis)

dp_tukey <- onion %>% group_by(CellType,geneid) %>%
  na.omit()%>%
  tukey_hsd(Percent ~ Diagnosis)
dp_tukey <- as.data.frame(dp_tukey[dp_tukey$CellType %in% cells.used,])

# Save all the files for supplemental info
write.csv (onion, file = "20201013_FigS3_DP_percentage.allCT.csv")
write.csv (onion_means, file = "20201013_FigS3_DP_percentage_means.allCT.csv")
write.csv(dp_tukey, file = "20201013_DP_Epi_tukey.csv")

# -----------------------------------
# Figure 1D: Venn Diagram for all CT
# -----------------------------------
# Extract out the cell ID from the ACE2+ TMPRSS2+, BSG+ TMPRSS2+ and HSPA5+ TMPRSS2+ cells
set1 <- rownames(sub1@meta.data)
set2 <- rownames(sub2@meta.data)
set3 <- rownames(sub3@meta.data)

# Make venn diagram using VennDiagram package
myCol <- brewer.pal(3, "Pastel2")

temp <- venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("ACE2+ TMPRSS2+" , "BSG+ TMPRSS2+ " , "HSPA5+ TMPRSS2+"),
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

# ------------------------------------
# Figure 1E: ACE2+ TMPRSS2+ cell fraction
# ------------------------------------
# Only plot Epithelial cells
cells.used <- c("AT1","AT2","Transitional AT2","KRT5-/KRT17+","Basal",
                "Ciliated Cells","Differentiating Ciliated Cells","PNEC","Ionocytes",
                "SCGB3A2+/SCGB1A1+", "SCGB3A2+","MUC5AC+","MUC5B+")

onion_means <- onion_means[onion_means$CellType %in% cells.used,]
onion_means1 <- onion_means[onion_means$geneid == "ACE2+/TMPRSS2+",]

# Plotting 
onion_means1$Diagnosis <- factor(onion_means1$Diagnosis,
                                  levels=c("Other-ILD","IPF","COPD","Control"))
ggplot(onion_means1, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5,
            position = position_dodge(width = 1),inherit.aes = TRUE) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  coord_flip() + 
  ylab("Percentage of ACE2+ TMPRSS2+ cells") 

# ------------------------------------
# Figure 1F: BSG+ TMPRSS2+ and HSPA5+ TMPRSS2+ cell fraction
# ------------------------------------
onion_means2 <- onion_means[onion_means$geneid != "ACE2+/TMPRSS2+",]

# Plotting 
ggplot(onion_means2, aes(x=CellType, y= Mean, fill = Diagnosis)) +
  geom_bar(stat="identity",position='dodge', width = 0.8) +
  facet_wrap(~geneid, scales = "free", ncol = 1, nrow = 2) +  
  theme_bw() +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2,
                position=position_dodge(.8)) +
  geom_text(aes(label = paste(Total)), size = 3, hjust = -0.5, angle = 90,
            position = position_dodge(width = 1),inherit.aes = TRUE) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) + 
  scale_y_continuous()+
  ylab("Percentage of cells") 



