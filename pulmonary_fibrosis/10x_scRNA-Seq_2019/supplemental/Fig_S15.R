# ======================================
# Environment parameters
# ======================================
set.seed(8342)
getwd()
Sys.Date()
main_dir <- getwd()
date <- gsub("-", "", Sys.Date())
dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))
getwd()

# ======================================
# Load libraries
# ======================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

# ======================================
# 
# ======================================
setwd("/scratch/syahn/iTALK")

#==============================================================================#
#### Load iTALK functions ####
#==============================================================================#
## iTALK function to find top x% expressed genes, edited by Stephanie
rawParse<-function(data,top_genes=20,stats='mean'){
  res=NULL
  cell_group<-unique(data$celltype)
  pb <- progress::progress_bar$new(total = length(cell_group))
  pb$tick(0)
  for(i in cell_group){
    sub_data<-data[data$celltype==i,]
    counts<-t(subset(sub_data,select=-celltype))
    counts<-apply(counts,2,function(x) {storage.mode(x) <- 'numeric'; x})
    if(stats=='mean'){
      temp<-data.frame(rowMeans(counts),i,stringsAsFactors = FALSE)
    }else if(stats=='median'){
      temp<-data.frame(apply(counts, 1, FUN = median),i,stringsAsFactors = FALSE)
    }else{
      print('error stats option')
    }
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
FindLR<-function(data_1,data_2=NULL,datatype,database){
  if(is.null(database)){
    load("/scratch/syahn/Seurat/LR_database.rda")
  }
  if(datatype=='mean count'){
    gene_list_1<-data_1
    if(is.null(data_2)){
      gene_list_2<-gene_list_1
    }else{
      gene_list_2<-data_2
    }
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
  }else if(datatype=='DEG'){
    gene_list_1<-data_1
    if(is.null(data_2)){
      gene_list_2<-gene_list_1
    }else{
      gene_list_2<-data_2
    }
    ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_1$gene)
    receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_2$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_1<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(gene_list_1[,c('gene','logFC','q.value','celltype')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_from_logFC=logFC,cell_from_q.value=q.value,cell_from=celltype) %>%
      left_join(gene_list_2[,c('gene','logFC','q.value','celltype')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_to_logFC=logFC,cell_to_q.value=q.value,cell_to=celltype)
    ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_2$gene)
    receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_1$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_2<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(gene_list_2[,c('gene','logFC','q.value','celltype')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_from_logFC=logFC,cell_from_q.value=q.value,cell_from=celltype) %>%
      left_join(gene_list_1[,c('gene','logFC','q.value','celltype')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_to_logFC=logFC,cell_to_q.value=q.value,cell_to=celltype)
    FilterTable<-rbind(FilterTable_1,FilterTable_2)
  }else{
    stop('Error: invalid data type')
  }
  
  FilterTable<-FilterTable[!duplicated(FilterTable),]
  res<-as.data.frame(FilterTable) %>% dplyr::rename(ligand=Ligand.ApprovedSymbol,receptor=Receptor.ApprovedSymbol)
  if(datatype=='DEG'){
    res<-res[!(res$cell_from_logFC==0.0001 & res$cell_to_logFC==0.0001),]
  }
  return(res)
}

#==============================================================================#
#### Read in gene expression data ####
#==============================================================================#
epi_mesen <- readRDS(file = "/scratch/lbui/20190623_Final_version/190702_epi_mesenchymal.rds")
data <- epi_mesen@assays$SCT@scale.data
data <- as.data.frame(t(data))
data$celltype <- (epi_mesen@meta.data$celltype)
data$Status<- (epi_mesen@meta.data$Status)
## test if transposition and addition of metadata label aligns
setequal(rownames(epi_mesen@meta.data[epi_mesen@meta.data$Status == "ILD",]), 
         rownames(data[data$Status == "ILD",]))


#==============================================================================#
#### find top 20% expressed genes (in each celltype) for Disease cells ####
#==============================================================================#
dis_data <- dplyr:: filter(data, Status == "ILD")
dis_20percenthighly_exprs_genes<-rawParse(dis_data,top_genes=20,stats='mean')
unique(dis_20percenthighly_exprs_genes$celltype)
write.table(dis_20percenthighly_exprs_genes, 
            file = "dis_20percenthighly_exprs_genes_allcelltypes.txt", 
            sep = "\t", quote = F)

dis_20percenthighly <- read.delim(file = "dis_20percenthighly_exprs_genes_allcelltypes.txt", 
                                  sep = "\t", header = T)


#==============================================================================#
#### Find LR pairs ####
#==============================================================================#
res<-NULL
load("/scratch/syahn/Seurat/LR_database.rda")
res<-FindLR(dis_20percenthighly,datatype='mean count',database=database)
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),]
write.table(res, file = "epi_mesen_FindLR_ILD_top20percentexpr_allcells.txt", 
            sep = "\t", row.names = F)

ILD_LRpairs <- read.delim(file = "epi_mesen_FindLR_ILD_top20percentexpr_allcells.txt", 
                          sep = "\t", header = T)


#==============================================================================#
#### From Epi to Mesen ####
#==============================================================================#
epi_to_mesen <- dplyr::filter(ILD_LRpairs, grepl("^Fibro|^PLIN|^Myofibro", cell_to)) %>%
  dplyr::filter(!grepl("ibro|Mesothelial Cells|Smooth Muscle Cells", cell_from))
unique(epi_to_mesen$cell_from)
unique(epi_to_mesen$cell_to)

mesen_cells <- c(as.character(unique(epi_to_mesen$cell_to)))
epi_cells <- c(as.character(unique(epi_to_mesen$cell_from)))

## Select the top 20% most highly co-expressed LR pairs
AT1_to_mesen<- list()
for (i in 1:length(mesen_cells)) {
  AT1_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "AT1") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
AT1_to_mesen <- bind_rows(AT1_to_mesen)

AT2_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  AT2_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "AT2") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
AT2_to_mesen <- bind_rows(AT2_to_mesen)

basal_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  basal_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Basal") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
basal_to_mesen <- bind_rows(basal_to_mesen)

ciliated_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  ciliated_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Ciliated") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
ciliated_to_mesen <- bind_rows(ciliated_to_mesen)

diffcil_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  diffcil_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Differentiating Ciliated") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
diffcil_to_mesen <- bind_rows(diffcil_to_mesen)

krt5neg_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  krt5neg_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "KRT5-/KRT17+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
krt5neg_to_mesen <- bind_rows(krt5neg_to_mesen)

muc5ac_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  muc5ac_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "MUC5AC+ High") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
muc5ac_to_mesen <- bind_rows(muc5ac_to_mesen)

muc5b_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  muc5b_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "MUC5B+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
muc5b_to_mesen <- bind_rows(muc5b_to_mesen)

prolifepi_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  prolifepi_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Proliferating Epithelial Cells") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
prolifepi_to_mesen <- bind_rows(prolifepi_to_mesen)

scgb3a2_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  scgb3a2_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "SCGB3A2+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
scgb3a2_to_mesen <- bind_rows(scgb3a2_to_mesen)

scgb3a2_1a1_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  scgb3a2_1a1_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "SCGB3A2+ SCGB1A1+") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
scgb3a2_1a1_to_mesen <- bind_rows(scgb3a2_1a1_to_mesen)

transAT2_to_mesen <- list()
for (i in 1:length(mesen_cells)) {
  transAT2_to_mesen[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Transitional AT2") %>%
    dplyr::filter(cell_to == mesen_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
transAT2_to_mesen <- bind_rows(transAT2_to_mesen)

top20percent_epi_to_mesen <- rbind(AT1_to_mesen,AT2_to_mesen,basal_to_mesen,ciliated_to_mesen,
                                   diffcil_to_mesen,krt5neg_to_mesen,muc5b_to_mesen,muc5ac_to_mesen,
                                   prolifepi_to_mesen,scgb3a2_to_mesen,scgb3a2_1a1_to_mesen,
                                   transAT2_to_mesen)
top20percent_epi_to_mesen %>% mutate_if(is.factor, as.character) -> top20percent_epi_to_mesen


## Find associations between ligand and receptor expression within individuals
LR_data = top20percent_epi_to_mesen
LR_data$ligand_receptor <- paste(LR_data$ligand, LR_data$receptor, sep = "_")

fit_matrix = matrix(ncol = 2, nrow = dim(LR_data)[1])
rownames(fit_matrix) <- paste( LR_data$ligand_receptor, LR_data$cell_from, LR_data$cell_to, sep = "_")

for(i in 1:dim(LR_data)[1]){
  
  onion = epi_mesen@assays$SCT@data[LR_data$ligand[i],]
  onion2 = epi_mesen@assays$SCT@data[LR_data$receptor[i],]
  
  onion = onion[epi_mesen@meta.data$celltype == as.character(LR_data$cell_from[i])]
  onion2 = onion2[epi_mesen@meta.data$celltype == as.character(LR_data$cell_to[i])]
  
  onion = as.data.frame(onion)
  onion2 = as.data.frame(onion2)
  
  onion$ind = do.call(rbind, strsplit(rownames(onion), split = "_"))[,1]
  onion2$ind = do.call(rbind, strsplit(rownames(onion2), split = "_"))[,1]
  
  onion_means = onion %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion, na.rm = T))
  onion2_means = onion2 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
  
  onion_means = as.data.frame(onion_means)
  onion_means = onion_means[onion_means$Mean > 0,]
  
  onion2_means = as.data.frame(onion2_means)
  onion2_means = onion2_means[onion2_means$Mean > 0,]
  
  means = merge(onion_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
  
  lmsum = summary(lm(means[,3] ~ means[,2]))
  fit_matrix[i,1] = ifelse(lmsum$residuals[1] == 0, "NaN", lmsum$coefficients[2,4])
  fit_matrix[i,2] = ifelse(lmsum$residuals[1] == 0, "NaN", lmsum$adj.r.squared)
  
}

LR_data_sig = LR_data[fit_matrix[,1] <= 0.01 & fit_matrix[,1] != "NaN",]
plot_list = list()
for(i in 1:dim(LR_data_sig)[1]){
  
  onion = epi_mesen@assays$SCT@data[LR_data_sig$ligand[i],]
  onion2 = epi_mesen@assays$SCT@data[LR_data_sig$receptor[i],]
  
  onion = onion[epi_mesen@meta.data$celltype == as.character(LR_data_sig$cell_from[i])]
  onion2 = onion2[epi_mesen@meta.data$celltype == as.character(LR_data_sig$cell_to[i])]
  
  onion = as.data.frame(onion)
  onion2 = as.data.frame(onion2)
  
  onion$ind = do.call(rbind, strsplit(rownames(onion), split = "_"))[,1]
  onion2$ind = do.call(rbind, strsplit(rownames(onion2), split = "_"))[,1]
  
  onion$Status = epi_mesen@meta.data[rownames(onion),]$Status
  onion2$Status = epi_mesen@meta.data[rownames(onion2),]$Status
  
  onion_means = onion %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion, na.rm = T))
  onion2_means = onion2 %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
  
  onion_means = as.data.frame(onion_means)
  onion_means = onion_means[onion_means$Mean > 0,]
  
  onion2_means = as.data.frame(onion2_means)
  onion2_means = onion2_means[onion2_means$Mean > 0,]
  
  means = merge(onion_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
  means = means[,-4]
  
  plot_list[[i]] <- ggplot(means, aes(x = Mean.x, y = Mean.y)) + geom_point(aes(color = Status.x), size = 0.5) + 
    ggtitle(paste(LR_data_sig[i,]$cell_from, "to", LR_data_sig[i,]$cell_to)) + 
    xlab(LR_data_sig[i,]$ligand) + ylab(LR_data_sig[i,]$receptor) + 
    theme(plot.title = element_text(size = 7), legend.title = element_blank(), axis.title=element_text(size=8), 
          axis.text = element_text(size = 5), legend.text = element_text(size = 5))
  
}

grid.arrange(grobs = plot_list, ncol = 5)


#==============================================================================#
#### From Mesen to Epi ####
#==============================================================================#
mesen_to_epi <- dplyr::filter(ILD_LRpairs, grepl("^Fibro|^PLIN|^Myofibro", cell_from)) %>%
  dplyr::filter(!grepl("ibro|Mesothelial Cells|Smooth Muscle Cells", cell_to))
unique(mesen_to_epi$cell_from)
unique(mesen_to_epi$cell_to)

mesen_cells <- c(as.character(unique(mesen_to_epi$cell_from)))
epi_cells <- c(as.character(unique(mesen_to_epi$cell_to)))

## Select the top 20% most highly co-expressed LR pairs
fibro_to_epi <- list()
for (i in 1:length(epi_cells)) {
  fibro_to_epi[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Fibroblasts") %>%
    dplyr::filter(cell_to == epi_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
fibro_to_epi <- bind_rows(fibro_to_epi)

myofibro_to_epi <- list()
for (i in 1:length(epi_cells)) {
  myofibro_to_epi[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "Myofibroblasts") %>%
    dplyr::filter(cell_to == epi_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
myofibro_to_epi <- bind_rows(myofibro_to_epi)

plin2_to_epi <- list()
for (i in 1:length(epi_cells)) {
  plin2_to_epi[[i]] <- dplyr::filter(ILD_LRpairs, cell_from == "PLIN2+ Fibroblasts") %>%
    dplyr::filter(cell_to == epi_cells[i]) %>% top_frac(0.2, cell_from_mean_exprs*cell_to_mean_exprs)
}
plin2_to_epi <- bind_rows(plin2_to_epi)

top20percent_mesen_to_epi <- rbind(fibro_to_epi,myofibro_to_epi,plin2_to_epi)
top20percent_mesen_to_epi %>% mutate_if(is.factor, as.character) -> top20percent_mesen_to_epi


## Find associations between ligand and receptor expression within individuals
LR_data2 = top20percent_mesen_to_epi
LR_data2$ligand_receptor <- paste(LR_data2$ligand, LR_data2$receptor, sep = "_")

fit_matrix2 = matrix(ncol = 2, nrow = dim(LR_data2)[1])
rownames(fit_matrix2) <- paste( LR_data2$ligand_receptor, LR_data2$cell_from, LR_data2$cell_to, sep = "_")

for(i in 1:dim(LR_data2)[1]){
  
  onion = epi_mesen@assays$SCT@data[LR_data2$ligand[i],]
  onion2 = epi_mesen@assays$SCT@data[LR_data2$receptor[i],]
  
  onion = onion[epi_mesen@meta.data$celltype == as.character(LR_data2$cell_from[i])]
  onion2 = onion2[epi_mesen@meta.data$celltype == as.character(LR_data2$cell_to[i])]
  
  onion = as.data.frame(onion)
  onion2 = as.data.frame(onion2)
  
  onion$ind = do.call(rbind, strsplit(rownames(onion), split = "_"))[,1]
  onion2$ind = do.call(rbind, strsplit(rownames(onion2), split = "_"))[,1]
  
  onion_means = onion %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion, na.rm = T))
  onion2_means = onion2 %>% group_by(ind) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
  
  
  onion_means = as.data.frame(onion_means)
  onion_means = onion_means[onion_means$Mean > 0,]
  
  onion2_means = as.data.frame(onion2_means)
  onion2_means = onion2_means[onion2_means$Mean > 0,]
  
  
  means = merge(onion_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
  
  lmsum = summary(lm(means[,3] ~ means[,2]))
  fit_matrix2[i,1] = ifelse(lmsum$residuals[1] == 0, "NaN", lmsum$coefficients[2,4])
  fit_matrix2[i,2] = ifelse(lmsum$residuals[1] == 0, "NaN", lmsum$adj.r.squared)
  
}


LR_data_sig2 = LR_data2[fit_matrix2[,1] <= 0.01 & fit_matrix2[,1] != "NaN",]
plot_list_2 = list()


for(i in 1:dim(LR_data_sig2)[1]){
  
  onion = epi_mesen@assays$SCT@data[LR_data_sig2$ligand[i],]
  onion2 = epi_mesen@assays$SCT@data[LR_data_sig2$receptor[i],]
  
  onion = onion[epi_mesen@meta.data$celltype == as.character(LR_data_sig2$cell_from[i])]
  onion2 = onion2[epi_mesen@meta.data$celltype == as.character(LR_data_sig2$cell_to[i])]
  
  onion = as.data.frame(onion)
  onion2 = as.data.frame(onion2)
  
  onion$ind = do.call(rbind, strsplit(rownames(onion), split = "_"))[,1]
  onion2$ind = do.call(rbind, strsplit(rownames(onion2), split = "_"))[,1]
  
  onion$Status = epi_mesen@meta.data[rownames(onion),]$Status
  onion2$Status = epi_mesen@meta.data[rownames(onion2),]$Status
  
  onion_means = onion %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion, na.rm = T))
  onion2_means = onion2 %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
  
  onion_means = as.data.frame(onion_means)
  onion_means = onion_means[onion_means$Mean > 0,]
  
  onion2_means = as.data.frame(onion2_means)
  onion2_means = onion2_means[onion2_means$Mean > 0,]
  
  means = merge(onion_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
  means = means[,-4]
  
  plot_list_2[[i]] <- ggplot(means, aes(x = Mean.x, y = Mean.y)) + geom_point(aes(color = Status.x), size = 0.5) + 
    ggtitle(paste(LR_data_sig2[i,]$cell_from, "to", LR_data_sig2[i,]$cell_to)) + 
    xlab(LR_data_sig2[i,]$ligand) + ylab(LR_data_sig2[i,]$receptor) + 
    theme(plot.title = element_text(size = 7), legend.title = element_blank(), axis.title=element_text(size=8), 
          axis.text = element_text(size = 5), legend.text = element_text(size = 5))
  
}

grid.arrange(grobs = plot_list_2, ncol = 5)

