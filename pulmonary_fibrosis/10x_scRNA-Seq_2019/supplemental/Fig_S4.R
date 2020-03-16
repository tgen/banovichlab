# ====================================
# Author : Linh T. Bui, lbui@tgen.org
# Date: 2020-03-11
# SCIENECES ADVANCES - IPF: FIGURE S4
# ===================================

# Set up working environment
library(Seurat)
library(dplyr)
library(UpSetR)
library(rlist)

set.seed(12345)

# Set date and create out folder
getwd()
Sys.Date()
main_dir <- "/scratch/lbui/RStudio_folder/"
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()

# Read in Seurat object
ild <- readRDS("/scratch/lbui/IPF_seurat_objects/190623_ILD_annotated_sizereduced.rds")

# Running DE analysis using negative binominal test
ild_list = list()
j=0
for(i in unique(ild@meta.data$celltype)){
  j=j+1
  ild_list[[j]] <- SubsetData(ild, cells = row.names(ild@meta.data[ild@meta.data$celltype == i,]))
}
for(i in 1:length(ild_list)){
  names(ild_list) <- lapply(ild_list, function(xx){paste(unique(xx@meta.data$celltype))})
}

disease_vs_control <- lapply(ild_list, function(xx){
  print(unique(xx@meta.data$celltype))
  if(length(unique(xx@meta.data$Status)) > 1) {
    FindMarkers(xx, group.by = "Status", ident.1 = "Control", ident.2 = "ILD", test.use = "negbinom")
  } 
  else{
    return(NULL)
  } 
})

for(i in 1:length(ild_list)){
  write.table(disease_vs_control[[i]], paste(gsub("/", "", unique(ild_list[[i]]@meta.data$celltype)), "_disease_vs_control", ".csv"), sep =",", quote = F)
}

# Prepare upset plots
## Remove cell types with less than 50 cells before making the upset plot
temp <- c("HAS1 High Fibroblasts","KRT5-/KRT17+","MUC5AC+ High", "pDCs", "PLIN2+ Fibroblasts", "Mesothelial Cells")
temp1 <- list.remove(disease_vs_control, c("HAS1 High Fibroblasts","KRT5-/KRT17+","MUC5AC+ High", "pDCs", "PLIN2+ Fibroblasts", "Mesothelial Cells"))
names <- unique(ild@meta.data$celltype)
names <- names[! names %in% temp]

## Make upset plots (for cell types) with the top 60 intersects
onion <- lapply(temp1, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(temp1, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]], onion2[[6]],
                      onion2[[7]], onion2[[8]], onion2[[9]], onion2[[10]], onion2[[11]], onion2[[12]],
                      onion2[[13]], onion2[[14]], onion2[[15]], onion2[[16]], onion2[[17]], onion2[[18]],
                      onion2[[19]], onion2[[20]], onion2[[21]], onion2[[22]], onion2[[23]], onion2[[24]],
                      onion2[[25]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

pdf("20200311_All_upset.pdf")
upset(upset_d_vs_c, nsets = 25, text.scale = 1, show.numbers = F, nintersects = 60)
dev.off()

## Make upset plots (for cell types per population)
# Epithelial cells
epi <- list.subset (temp1, c("Basal","SCGB3A2+ SCGB1A1+","AT2","Proliferating Epithelial Cells",
                       "SCGB3A2+", "AT1", "Differentiating Ciliated", "MUC5B+",
                       "Ciliated", "Transitional AT2"))

onion <- lapply(epi, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(epi, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]], onion2[[6]],
                      onion2[[7]], onion2[[8]], onion2[[9]], onion2[[10]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(epi)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

pdf("20200311_Epi_upset.pdf")
upset(upset_d_vs_c, nsets = 25, text.scale = 1, show.numbers = F, nintersects = 60)
dev.off()

# Endothelial and Mesenchymal cells
endo <- list.subset (temp1, c("Lymphatic Endothelial Cells", "Endothelial Cells", "Smooth Muscle Cells","Fibroblasts","Myofibroblasts"))

onion <- lapply(endo, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(endo, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]], onion2[[6]],
                      onion2[[7]], onion2[[8]], onion2[[9]], onion2[[10]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(endo)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

pdf("20200311_Endo_mesen_upset.pdf")
upset(upset_d_vs_c, nsets = 25, text.scale = 1, show.numbers = F, nintersects = 60)
dev.off()

# Immune cells
immune <- list.subset (temp1, c("Proliferating Macrophages","Mast Cells", "T Cells", "Plasma Cells",
                                "Macrophages","Proliferating T Cells","NK Cells", "B Cells", "cDCs", "Monocytes"))

onion <- lapply(immune, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(immune, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_d_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_d_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], onion2[[5]], onion2[[6]],
                      onion2[[7]], onion2[[8]], onion2[[9]], onion2[[10]])
upset_d_vs_c <- as.data.frame(upset_d_vs_c)

row.names(upset_d_vs_c) <- onion
colnames(upset_d_vs_c) <- names(immune)

upset_d_vs_c[upset_d_vs_c == T] <- 1
upset_d_vs_c[upset_d_vs_c == F] <- 0

pdf("20200311_Immune_upset.pdf")
upset(upset_d_vs_c, nsets = 25, text.scale = 1, show.numbers = F, nintersects = 60)
dev.of()

