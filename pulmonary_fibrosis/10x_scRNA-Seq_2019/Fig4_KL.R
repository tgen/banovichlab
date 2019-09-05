#==============================================================================#
# Author(s) : Stephanie L. Yahn, syahn@tgen.org
#             https://github.com/Coolgenome/iTALK
# Date: 01/08/2019
# Description: Ligand Receptor Analysis for IPF Project; 
#               Cytoscape Network Plots
#==============================================================================#

#==============================================================================#
#### Load libraries ####
#==============================================================================#
setwd("/scratch/syahn/")
library(Seurat)
library(dplyr)

#==============================================================================#
#### Load iTALK functions ####
#==============================================================================#
## iTALK function to find top x% mean expressed genes, edited by Stephanie
rawParse<-function(data,top_genes=20){
  res=NULL
  cell_group<-unique(data$celltype)
  pb <- progress::progress_bar$new(total = length(cell_group))
  pb$tick(0)
  for(i in cell_group){
    sub_data<-data[data$celltype==i,]
    counts<-t(subset(sub_data,select=-celltype))
    counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
    temp<-data.frame(rowMeans(counts),i,stringsAsFactors = FALSE)
    temp<-temp[order(temp[,1],decreasing=TRUE),]
    temp<-temp[1:ceiling(nrow(temp)*top_genes/100),]
    temp<-temp %>% tibble::rownames_to_column()
    res<-rbind(res,temp)
    pb$tick()
  }
  colnames(res)<-c('gene','exprs','celltype')
  return(res)
}
## iTALK function to find LR pairs, edited by Stephanie
FindLR<-function(data_1){
  load("/scratch/syahn/Seurat/LR_database.rda")
  gene_list_1<-data_1
  gene_list_2<-gene_list_1
  ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_1$gene)
  receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_2$gene)
  ind<-intersect(ligand_ind,receptor_ind)
  FilterTable_1<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
    left_join(gene_list_1[,c('gene','exprs','celltype')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
    dplyr::rename(cell_from_mean_exprs=exprs,cell_from=celltype) %>%
    left_join(gene_list_2[,c('gene','exprs','celltype')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
    dplyr::rename(cell_to_mean_exprs=exprs,cell_to=celltype)
  ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_2$gene)
  receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_1$gene)
  ind<-intersect(ligand_ind,receptor_ind)
  FilterTable_2<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
    left_join(gene_list_2[,c('gene','exprs','celltype')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
    dplyr::rename(cell_from_mean_exprs=exprs,cell_from=celltype) %>%
    left_join(gene_list_1[,c('gene','exprs','celltype')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
    dplyr::rename(cell_to_mean_exprs=exprs,cell_to=celltype)
  FilterTable<-rbind(FilterTable_1,FilterTable_2)
  FilterTable<-FilterTable[!duplicated(FilterTable),]
  res<-as.data.frame(FilterTable) %>% dplyr::rename(ligand=Ligand.ApprovedSymbol,receptor=Receptor.ApprovedSymbol)
  return(res)
}

#==============================================================================#
#### Read in gene expression data ####
#==============================================================================#
epi_mesen <- readRDS(file = "/scratch/lbui/20190623_Final_version/190702_epi_mesenchymal.rds")
epi_mesen <- FindVariableFeatures(epi_mesen, verbose = T, nfeatures = 3000)
epi_mesen <- ScaleData(epi_mesen, features = row.names(epi_mesen@assays$SCT@data))
data <- epi_mesen@assays$SCT@scale.data
data <- as.data.frame(t(data))
data$celltype <- (epi_mesen@meta.data$celltype)
data$Status <- (epi_mesen@meta.data$Status)
## sanity check
setequal(rownames(epi_mesen@meta.data[epi_mesen@meta.data$Status == "ILD",]), 
         rownames(data[data$Status == "ILD",]))

#==============================================================================#
#### find top 20% expressed genes (in each celltype) for Disease cells ####
#==============================================================================#
dis_data <- dplyr:: filter(data, Status == "ILD")
## this step takes several minutes
dis_20pct <- rawParse(dis_data,top_genes=20)

#==============================================================================#
#### Find LR pairs ####
#==============================================================================#
ILD_LRpairs<-NULL
load("/scratch/syahn/Seurat/LR_database.rda")
ILD_LRpairs <- FindLR(dis_20pct)
ILD_LRpairs <- ILD_LRpairs[order(ILD_LRpairs$cell_from_mean_exprs*ILD_LRpairs$cell_to_mean_exprs,decreasing=T),]

#==============================================================================#
####From Mesen to Epi (Disease) ####
#==============================================================================#
mesen_to_epi <- dplyr::filter(ILD_LRpairs, grepl("^Fibro|^PLIN|^Myofibro", cell_from)) %>%
  dplyr::filter(!grepl("ibro|Mesothelial Cells|Smooth Muscle Cells", cell_to))
## sanity check
unique(mesen_to_epi$cell_from)
unique(mesen_to_epi$cell_to)

mesen_cells <- c(as.character(unique(mesen_to_epi$cell_from)))
epi_cells <- c(as.character(unique(mesen_to_epi$cell_to)))

fibro_to_epi <- list()
for (i in 1:length(epi_cells)) {
  fibro_to_epi[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Fibroblasts") %>%
    dplyr::filter(cell_to == epi_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
fibro_to_epi <- bind_rows(fibro_to_epi)

myofibro_to_epi <- list()
for (i in 1:length(epi_cells)) {
  myofibro_to_epi[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Myofibroblasts") %>%
    dplyr::filter(cell_to == epi_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
myofibro_to_epi <- bind_rows(myofibro_to_epi)

plin2_to_epi <- list()
for (i in 1:length(epi_cells)) {
  plin2_to_epi[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "PLIN2+ Fibroblasts") %>%
    dplyr::filter(cell_to == epi_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
plin2_to_epi <- bind_rows(plin2_to_epi)

top5_mesen_to_epi <- rbind(fibro_to_epi,myofibro_to_epi,plin2_to_epi)
top5_mesen_to_epi %>% mutate_if(is.factor, as.character) -> top5_mesen_to_epi
write.table(top5_mesen_to_epi, file = "top5_mesen_to_epi.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)

#==============================================================================#
#### From Epi to Mesen (Disease) ####
#==============================================================================#
epi_to_mesen <- dplyr::filter(ILD_LRpairs, grepl("^Fibro|^PLIN|^Myofibro", cell_to)) %>%
  dplyr::filter(!grepl("ibro|Mesothelial Cells|Smooth Muscle Cells", cell_from))
## sanity check
unique(epi_to_mesen$cell_from)
unique(epi_to_mesen$cell_to)

mesen_cells <- c(as.character(unique(epi_to_mesen$cell_to)))
epi_cells <- c(as.character(unique(epi_to_mesen$cell_from)))

AT1_to_mesen<- list()
for (i in 1:length(mesen_cells)) {
  AT1_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "AT1") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
AT1_to_mesen <- bind_rows(AT1_to_mesen)

AT2_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  AT2_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "AT2") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
AT2_to_mesen <- bind_rows(AT2_to_mesen)

basal_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  basal_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Basal") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
basal_to_mesen <- bind_rows(basal_to_mesen)

ciliated_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  ciliated_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Ciliated") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
ciliated_to_mesen <- bind_rows(ciliated_to_mesen)

diffcil_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  diffcil_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Differentiating Ciliated") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
diffcil_to_mesen <- bind_rows(diffcil_to_mesen)

krt5neg_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  krt5neg_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "KRT5-/KRT17+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
krt5neg_to_mesen <- bind_rows(krt5neg_to_mesen)

muc5ac_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  muc5ac_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "MUC5AC+ High") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
muc5ac_to_mesen <- bind_rows(muc5ac_to_mesen)

muc5b_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  muc5b_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "MUC5B+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
muc5b_to_mesen <- bind_rows(muc5b_to_mesen)

prolifepi_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  prolifepi_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Proliferating Epithelial Cells") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
prolifepi_to_mesen <- bind_rows(prolifepi_to_mesen)

scgb3a2_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  scgb3a2_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "SCGB3A2+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
scgb3a2_to_mesen <- bind_rows(scgb3a2_to_mesen)

scgb3a2_1a1_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  scgb3a2_1a1_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "SCGB3A2+ SCGB1A1+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
scgb3a2_1a1_to_mesen <- bind_rows(scgb3a2_1a1_to_mesen)

transAT2_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  transAT2_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Transitional AT2") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_n(5, cell_from_mean_exprs*cell_to_mean_exprs)
}
transAT2_to_mesen <- bind_rows(transAT2_to_mesen)


top5_epi_to_mesen <- rbind(AT1_to_mesen,AT2_to_mesen,basal_to_mesen,ciliated_to_mesen,
                           diffcil_to_mesen,krt5neg_to_mesen,muc5ac_to_mesen,muc5b_to_mesen,
                           prolifepi_to_mesen,scgb3a2_to_mesen,scgb3a2_1a1_to_mesen,
                           transAT2_to_mesen)
top5_epi_to_mesen %>% mutate_if(is.factor, as.character) -> top5_epi_to_mesen

write.table(top5_epi_to_mesen, file = "top5_epi_to_mesen.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)


#==============================================================================#
#### Read in lists of DE genes for each cell type ####
#==============================================================================#
celltypes <- as.character(unique(epi_mesen$celltype))
celltypes <- celltypes[!celltypes %in% c("Mesothelial Cells", "Smooth Muscle Cells", 
                                         "HAS1 High Fibroblasts")]

setwd("/scratch/lbui/20190731_IPF_Seurat/DE_analysis/Disease_vs_Control/")
list.filenames<-list.files(pattern=".csv$")
list.filenames
list.filenames <- list.filenames[c(1,2,4,6,7,9,11,17,18,19,23,24,27,28,31)]
list.degs<-list()
for (i in 1:length(list.filenames)){
  list.degs[[i]]<-read.csv(list.filenames[i], row.names = NULL)
  colnames(list.degs[[i]]) <- c("gene", "p_val", "logFC", "pct.1", "pct.2", "q.value")
}
names(list.degs)<-celltypes

for (i in 1:length(celltypes)) {
  list.degs[[celltypes[i]]]$celltype <- celltypes[i]
}

setwd("/scratch/syahn/")

#==============================================================================#
#### Cytoscape network plot ####
#==============================================================================#
## mesen DE genes
mesen_deg <- bind_rows(list.degs[mesen_cells])
mesen_L_deg <- list()
for (i in 1:length(mesen_cells)){
  mesen_L_deg[[i]] <- mesen_deg[mesen_deg$celltype==mesen_cells[i],] %>% dplyr::filter(q.value <= 0.1) %>%
    dplyr::filter(gene %in% top5_mesen_to_epi$ligand[grep(paste("^", mesen_cells[i], sep = ""), 
                                                          top5_mesen_to_epi$cell_from)])
}
mesen_L_deg <- bind_rows(mesen_L_deg)
mesen_L_deg <- na.omit(mesen_L_deg)
mesen_L_deg <- mesen_L_deg[,c(1,3,6,7)]
colnames(mesen_L_deg) <- c("ligand", "ligand_logFC", "ligand_qvalue", "cell_from")

mesen_R_deg <- list()
for (i in 1:length(mesen_cells)){
  mesen_R_deg[[i]] <- mesen_deg[mesen_deg$celltype==mesen_cells[i],] %>% dplyr::filter(q.value <= 0.1) %>%
    dplyr::filter(gene %in% top5_epi_to_mesen$receptor[grep(paste("^",mesen_cells[i], sep = ""), 
                                                            top5_epi_to_mesen$cell_to)])
}
mesen_R_deg <- bind_rows(mesen_R_deg)
mesen_R_deg <- na.omit(mesen_R_deg)
mesen_R_deg <- mesen_R_deg[,c(1,3,6,7)]
colnames(mesen_R_deg) <- c("receptor", "receptor_logFC", "receptor_qvalue", "cell_to")


## epi DE genes
epi_deg <- bind_rows(list.degs[epi_cells])
epi_R_deg <- list()
for (i in 1:length(epi_cells)){
  epi_R_deg[[i]] <- epi_deg[epi_deg$celltype==epi_cells[i],] %>% dplyr::filter(q.value <= 0.1) %>%
    dplyr::filter(gene %in% top5_mesen_to_epi$receptor[grep(paste("^",epi_cells[i], sep = ""), 
                                                            top5_mesen_to_epi$cell_to)])
}
epi_R_deg <- bind_rows(epi_R_deg)
epi_R_deg <- na.omit(epi_R_deg)
epi_R_deg <- epi_R_deg[,c(1,3,6,7)]
colnames(epi_R_deg) <- c("receptor", "receptor_logFC", "receptor_qvalue", "cell_to")

epi_L_deg <- list()
for (i in 1:length(epi_cells)){
  epi_L_deg[[i]] <- epi_deg[epi_deg$celltype==epi_cells[i],] %>% dplyr::filter(q.value <= 0.1) %>%
    dplyr::filter(gene %in% top5_epi_to_mesen$ligand[grep(paste("^",epi_cells[i], sep = ""), 
                                                          top5_epi_to_mesen$cell_from)])
}
epi_L_deg <- bind_rows(epi_L_deg)
epi_L_deg <- na.omit(epi_L_deg)
epi_L_deg <- epi_L_deg[,c(1,3,6,7)]
colnames(epi_L_deg) <- c("ligand", "ligand_logFC", "ligand_qvalue", "cell_from")

## add DE gene info to LR dataframes
top5_epi_to_mesen_deg <- merge(top5_epi_to_mesen, epi_L_deg, all.x = T) %>%
  merge(mesen_R_deg, all.x = T)

write.table(top5_epi_to_mesen_deg, file = "top5_epi_to_mesen_deg.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)


top5_mesen_to_epi_deg <- merge(top5_mesen_to_epi, mesen_L_deg, all.x = T) %>%
  merge(epi_R_deg, all.x = T)

write.table(top5_mesen_to_epi_deg, file = "top5_mesen_to_epi_deg.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)