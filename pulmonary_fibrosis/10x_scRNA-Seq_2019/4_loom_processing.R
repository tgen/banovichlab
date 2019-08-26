# ==============================================================================
# Author(s) : Austin Gutierrez, agutierrez@tgen.org
# Date : 02/07/19
# Description: Script for processing velocyto loom files
# ==============================================================================
# ======================================
# Environment parameters
# ======================================
set.seed(12345)

# ======================================
# Load libraries
# ======================================
library(loomR)
library(BiocGenerics)
library(velocyto.R)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ade4)
library(Matrix)

# ======================================
# Read in loom files
# ======================================
loom_1 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_THD0001_1_LU_Whole_C1_X5SCR_F01157_H75C7DRXX/velocyto/IPF_THD0001_1_LU_Whole_C1_X5SCR_F01157_H75C7DRXX.loom")
loom_2 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD028_1_MF_Whole_C12_X5SCR_F01380_H75C7DRXX/velocyto/IPF_TILD028_1_MF_Whole_C12_X5SCR_F01380_H75C7DRXX.loom")
loom_3 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD019_1_LF_Whole_C3_X5SCR_F01375_H75C7DRXX/velocyto/IPF_TILD019_1_LF_Whole_C3_X5SCR_F01375_H75C7DRXX.loom")
loom_4 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD68_1_LU_Whole_C1_X5SCR_HD68_XXXXXXXXX/velocyto/IPF_VUHD68_1_LU_Whole_C1_X5SCR_HD68_XXXXXXXXX.loom")
loom_5 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD65_1_MF_Whole_C3_X5SCR_F01392_H75C7DRXX/velocyto/IPF_VUILD64_1_MF_Whole_C3_X5SCR_F01392_H75C7DRXX.loom")
loom_6 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD65_1_LU_Whole_C1_X5SCR_HD65_XXXXXXXXX/velocyto/IPF_VUHD65_1_LU_Whole_C1_X5SCR_HD65_XXXXXXXXX.loom")
loom_7 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD001_1_LF_Whole_C1_X5SCR_F00431_HWYTFBBXX/velocyto/IPF_TILD001_1_LF_Whole_C1_X5SCR_F00431_HWYTFBBXX.loom")
loom_8 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD54_1_LU_Whole_C1_X5SCR_F00207_HMWLCBGX7/velocyto/IPF_VUILD54_1_LU_Whole_C1_X5SCR_F00207_HMWLCBGX7.loom")
loom_9 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD57_1_LU_Whole_C3_X5SCR_ILD57_XXXXXXXXX/velocyto/IPF_VUILD57_1_LU_Whole_C3_X5SCR_ILD57_XXXXXXXXX.loom")
loom_10 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD006_2_MF_Whole_C2_X5SCR_F01303_H75C7DRXX/velocyto/IPF_TILD006_2_MF_Whole_C2_X5SCR_F01303_H75C7DRXX.loom")
loom_11 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD61_1_LU_Whole_C1_X5SCR_ILD61-2_XXXXXXXXX/velocyto/IPF_VUILD61_1_LU_Whole_C1_X5SCR_VUILD61-more_XXXXXXXXX.loom")
loom_12 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD60_1_LU_Whole_C1_X5SCR_ILD60-2_XXXXXXXXX/velocyto/IPF_VUILD60_1_LU_Whole_C1_X5SCR_VUILD60-more_XXXXXXXXX.loom")
loom_13 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD62_1_LU_Whole_C1_X5SCR_ILD62-2_XXXXXXXXX/velocyto/IPF_VUILD62_1_LU_Whole_C1_X5SCR_VUILD62-more_XXXXXXXXX.loom")
loom_14 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD64_1_LF_Whole_C1_X5SCR_F01390_H75C7DRXX/velocyto/IPF_VUILD64_1_LF_Whole_C1_X5SCR_F01390_H75C7DRXX.loom")
loom_15 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_THD0005_1_LU_Whole_C6_X5SCR_F01366_H75C7DRXX/velocyto/IPF_THD0005_1_LU_Whole_C6_X5SCR_F01366_H75C7DRXX.loom")
loom_16 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD71_1_LU_Whole_C1_X5SCR_F01394_H75C7DRXX/velocyto/IPF_VUILD94_1_LU_Whole_C1_X5SCR_F01394_H75C7DRXX.loom")
loom_17 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD006_2_LF_Whole_C1_X5SCR_F01302_H75C7DRXX/velocyto/IPF_TILD006_2_LF_Whole_C1_X5SCR_F01302_H75C7DRXX.loom")
loom_18 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD62_1_LU_Whole_C1_X5SCR_ILD62-1_XXXXXXXXX/velocyto/IPF_VUILD62_1_LU_Whole_C1_X5SCR_VUILD62-less_XXXXXXXXX.loom")
loom_19 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD015_1_MF_Whole_C4_X5SCR_F01173_H75C7DRXX/velocyto/IPF_TILD015_1_MF_Whole_C4_X5SCR_F01173_H75C7DRXX.loom")
loom_20 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_THD0005_1_LU_Whole_C7_X5SCR_F01367_H75C7DRXX/velocyto/IPF_THD0005_1_LU_Whole_C7_X5SCR_F01367_H75C7DRXX.loom")
loom_21 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD019_1_MF_Whole_C4_X5SCR_F01376_H75C7DRXX/velocyto/IPF_TILD019_1_MF_Whole_C4_X5SCR_F01376_H75C7DRXX.loom")
loom_22 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD66_1_LU_Whole_C2_X5SCR_HD66_XXXXXXXXX/velocyto/IPF_VUHD66_1_LU_Whole_C2_X5SCR_HD66_XXXXXXXXX.loom")
loom_23 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD54_1_LU_Whole_C2_X5SCR_F00208_HMWLCBGX7/velocyto/IPF_VUILD54_1_LU_Whole_C2_X5SCR_F00208_HMWLCBGX7.loom")
loom_24 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD60_1_LU_Whole_C1_X5SCR_ILD60-1_XXXXXXXXX/velocyto/IPF_VUILD60_1_LU_Whole_C1_X5SCR_VUILD60-less_XXXXXXXXX.loom")
loom_25 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD028_1_LF_Whole_C11_X5SCR_F01379_H75C7DRXX/velocyto/IPF_TILD028_1_LF_Whole_C11_X5SCR_F01379_H75C7DRXX.loom")
loom_26 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD59_1_LU_Whole_C1_X5SCR_ILD59-1_XXXXXXXXX/velocyto/IPF_VUILD59_1_LU_Whole_C1_X5SCR_ILD59-less_XXXXXXXXX.loom")
loom_27 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD015_1_LF_Whole_C3_X5SCR_F01172_H75C7DRXX/velocyto/IPF_TILD015_1_LF_Whole_C3_X5SCR_F01172_H75C7DRXX.loom")
loom_28 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_THD0002_1_LU_Whole_C3_X5SCR_F01174_H75C7DRXX/velocyto/IPF_THD0002_1_LU_Whole_C3_X5SCR_F01174_H75C7DRXX.loom")
loom_29 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD010_2_MF_MNC_C1_X5SCR_F01214_H75C7DRXX/velocyto/IPF_TILD010_2_MF_MNC_C1_X5SCR_F01214_H75C7DRXX.loom")
loom_30 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_THD0005_1_LU_Whole_C5_X5SCR_F01365_H75C7DRXX/velocyto/IPF_THD0005_1_LU_Whole_C5_X5SCR_F01365_H75C7DRXX.loom")
loom_31 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD58_1_LU_Whole_C1_X5SCR_ILD58_XXXXXXXXX/velocyto/IPF_VUILD58_1_LU_Whole_C1_X5SCR_ILD58_XXXXXXXXX.loom")
loom_32 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD48_1_LU_SCS_C1_X5SCR_F00202_HMWLCBGX7/velocyto/IPF_VUILD48_1_LU_SCS_C1_X5SCR_F00202_HMWLCBGX7.loom")
loom_33 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD67_1_LU_Whole_C1_X5SCR_HD67_XXXXXXXXX/velocyto/IPF_VUHD67_1_LU_Whole_C1_X5SCR_HD67_XXXXXXXXX.loom")
loom_34 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD030_1_MF_Whole_C4_X5SCR_F01386_H75C7DRXX/velocyto/IPF_TILD030_1_MF_Whole_C4_X5SCR_F01386_H75C7DRXX.loom")
loom_35 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD64_1_LF_Whole_C2_X5SCR_F01391_H75C7DRXX/velocyto/IPF_VUILD64_1_LF_Whole_C2_X5SCR_F01391_H75C7DRXX.loom")
loom_36 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_TILD030_1_MF_Whole_C4_X5SCR_F01385_H75C7DRXX/velocyto/IPF_TILD030_1_MF_Whole_C4_X5SCR_F01385_H75C7DRXX.loom")
loom_37 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD59_1_LU_Whole_C1_X5SCR_ILD59-2_XXXXXXXXX/velocyto/IPF_VUILD59_1_LU_Whole_C1_X5SCR_ILD59-more_XXXXXXXXX.loom")
loom_38 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD69_1_LU_Whole_C1_X5SCR_F00409_HWYTFBBXX/velocyto/IPF_VUHD069_1_LU_Whole_C1_X5SCR_F00409_HWYTFBBXX.loom")
loom_39 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD61_1_LU_Whole_C1_X5SCR_ILD61-1_XXXXXXXXX/velocyto/IPF_VUILD61_1_LU_Whole_C1_X5SCR_VUILD61-less_XXXXXXXXX.loom")
loom_40 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD53_1_LU_Whole_C1_X5SCR_ILD53_XXXXXXXXX/velocyto/IPF_VUILD53_1_LU_Whole_C1_X5SCR_ILD53_XXXXXXXXX.loom")
loom_41 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD55_1_LU_Whole_C1_X5SCR_ILD55_XXXXXXXXX/velocyto/IPF_VUILD55_1_LU_Whole_C1_X5SCR_ILD55_XXXXXXXXX.loom")
loom_42 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUILD63_1_LU_Whole_C1_X5SCR_ILD63_XXXXXXXXX/velocyto/IPF_VUILD63_1_LU_Whole_C1_X5SCR_ILD63_XXXXXXXXX.loom")
loom_43 <- read.loom.matrices("/scratch/agutierrez/10x_fastq/Outs/IPF/IPF_VUHD70_1_LU_Whole_C1_X5SCR_HD70_XXXXXXXXX/velocyto/IPF_VUHD70_1_LU_Whole_C1_X5SCR_HD70_XXXXXXXXX.loom")

boris <- c(loom_1, loom_2, loom_3, loom_4, loom_5,
           loom_6, loom_7, loom_8, loom_9, loom_10,
           loom_11, loom_12, loom_13, loom_14, loom_15,
           loom_16, loom_17, loom_18, loom_19, loom_20,
           loom_21, loom_22, loom_23, loom_24, loom_25,
           loom_26, loom_27, loom_28, loom_29, loom_30,
           loom_31, loom_32, loom_33, loom_34, loom_35,
           loom_36, loom_37, loom_38, loom_39, loom_40,
           loom_41, loom_42, loom_43)

ldat <- combine(boris)

ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("_unique.bam","",gsub(":","_",colnames(x)))
  colnames(x) <-  gsub("_unique.bam","",gsub("x","",colnames(x)))
  # names <- do.call(rbind, strsplit(colnames(x), "_"))
  # names <- as.data.frame(names)
  # names <- paste(names[8], names[10],sep = "_")
  # colnames(x) <- names
  x
})

saveRDS(ldat, "new_ldat.rds")

rm(boris)
rm (loom_1, loom_2, loom_3, loom_4, loom_5,
    loom_6, loom_7, loom_8, loom_9, loom_10,
    loom_11, loom_12, loom_13, loom_14, loom_15,
    loom_16, loom_17, loom_18, loom_19, loom_20,
    loom_21, loom_22, loom_23, loom_24, loom_25,
    loom_26, loom_27, loom_28, loom_29, loom_30,
    loom_31, loom_32, loom_33, loom_34, loom_35,
    loom_36, loom_37, loom_38, loom_39, loom_40,
    loom_41, loom_42, loom_43)

emat<- cbind(ldat[[1]],ldat[[4]],ldat[[7]],ldat[[10]],
             ldat[[13]],ldat[[16]],ldat[[19]],ldat[[22]],
             ldat[[25]],ldat[[28]],ldat[[31]],ldat[[34]],
             ldat[[37]],ldat[[40]],ldat[[43]],ldat[[46]],
             ldat[[49]],ldat[[52]],ldat[[55]],ldat[[58]],
             ldat[[61]],ldat[[64]],ldat[[67]],ldat[[70]],
             ldat[[73]],ldat[[76]],ldat[[79]],ldat[[82]],
             ldat[[85]],ldat[[88]],ldat[[91]],ldat[[94]],
             ldat[[97]],ldat[[100]],ldat[[103]],ldat[[106]],
             ldat[[109]],ldat[[112]],ldat[[115]],ldat[[118]],
             ldat[[121]],ldat[[124]],ldat[[127]])

nmat<- cbind(ldat[[2]],ldat[[5]],ldat[[8]],ldat[[11]],
             ldat[[14]],ldat[[17]],ldat[[20]],ldat[[23]],
             ldat[[26]],ldat[[29]],ldat[[32]],ldat[[35]],
             ldat[[38]],ldat[[41]],ldat[[44]],ldat[[47]],
             ldat[[50]],ldat[[53]],ldat[[56]],ldat[[59]],
             ldat[[62]],ldat[[65]],ldat[[68]],ldat[[71]],
             ldat[[74]],ldat[[77]],ldat[[80]],ldat[[83]],
             ldat[[86]],ldat[[89]],ldat[[92]],ldat[[95]],
             ldat[[98]],ldat[[101]],ldat[[104]],ldat[[107]],
             ldat[[110]],ldat[[113]],ldat[[116]],ldat[[119]],
             ldat[[122]],ldat[[125]],ldat[[128]])

smat<- cbind(ldat[[3]],ldat[[6]],ldat[[9]],ldat[[12]],
             ldat[[15]],ldat[[18]],ldat[[21]],ldat[[24]],
             ldat[[27]],ldat[[30]],ldat[[33]],ldat[[36]],
             ldat[[39]],ldat[[42]],ldat[[45]],ldat[[48]],
             ldat[[51]],ldat[[54]],ldat[[57]],ldat[[60]],
             ldat[[63]],ldat[[66]],ldat[[69]],ldat[[72]],
             ldat[[75]],ldat[[78]],ldat[[81]],ldat[[84]],
             ldat[[87]],ldat[[90]],ldat[[93]],ldat[[96]],
             ldat[[99]],ldat[[102]],ldat[[105]],ldat[[108]],
             ldat[[111]],ldat[[114]],ldat[[117]],ldat[[120]],
             ldat[[123]],ldat[[126]],ldat[[129]])


# emat = cbind(ldat[[seq(1,115,3)]])
# nmat = cbind(ldat[seq(2,116,3)])
# smat = cbind(ldat[seq(3,117,3)])

hist(log10(Matrix::rowSums(epi_emat)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene (emat)')
hist(log10(Matrix::rowSums(epi_nmat)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene (nmat)')
hist(log10(Matrix::rowSums(epi_smat)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene (smat)')

# ======================================
# Subset epithelial cells
# ======================================
epi <- readRDS("/scratch/lbui/20190623_Final_version/190623_EpiA_sizereduced.rds")
onion <- as.character(epi@meta.data$Diagnosis)
onion[onion == "sacroidosis"] <- "Sarcoidosis"
epi@meta.data$Diagnosis <- onion


DimPlot(epi, group.by = "celltypeDimplotOrder")

#epi_subset <- ldat[["spliced"]][colnames(ldat[[1]]) %in% rownames(epi@meta.data)]
#epi_subset <- subset(emat, colnames(emat) %in% rownames(epi@meta.data))
epi_emat <- emat[,colnames(emat) %in% rownames(epi@meta.data)]
rm(emat)

saveRDS(epi_emat, "epi_emat.rds")
epi_emat <- readRDS("epi_emat.rds")

epi_nmat <- nmat[,colnames(nmat) %in% rownames(epi@meta.data)]
rm(nmat)

saveRDS(epi_nmat, "epi_nmat.rds")
epi_nmat <- readRDS("epi_nmat.rds")

epi_smat <- smat[,colnames(smat) %in% rownames(epi@meta.data)]
rm(smat)

saveRDS(epi_smat, "epi_smat.rds")
epi_smat <- readRDS("epi_smat.rds")

epi_emb <- epi@reductions$umap@cell.embeddings
pca_emb <- epi@reductions$pca@cell.embeddings
epi_color <- as.character(epi@meta.data$celltype)
epi_color <- setNames(epi_color, rownames(epi@meta.data))

onion <- epi_color
onion[onion == "MUC5AC+ High"] <- "#00B6C8"
onion[onion == "Basal"] <- "#D38402"
onion[onion == "Proliferating Epithelial Cells"] <- "#B874FF"
onion[onion == "Differentiating Ciliated"] <- "#9D9A01"
onion[onion == "Ciliated"] <- "#B19302"
onion[onion == "SCGB3A2+ SCGB1A1+"] <- "#FD56C7"
onion[onion == "Differentiating Ciliated"] <- "#9D9A01"
onion[onion == "SCGB3A2+"] <- "#F659DD"
onion[onion == "MUC5B+"] <- "#00B2DB"
onion[onion == "KRT5-/KRT17+"] <- "#03AF21"
onion[onion == "Transitional AT2"] <- "#FE627D"
onion[onion == "AT2"] <- "#EE7342"
onion[onion == "AT1"] <- "#F76A62"
epi_color <- onion

saveRDS(epi_color, "epi_color.rds")

