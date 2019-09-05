library(Seurat)
library(dplyr)
library(UpSetR)
library(calibrate)
library(MAST)
library(DESeq2)

set.seed(2811)


## DEG for control vs IPF, all cell types and make Volcanoplot

ild <- readRDS("ILD.rds")
all_deg <- FindMarkers(ild, group.by = "Diagnosis", ident.1 = "IPF", ident.2 = "Control",
                      test.use = "negbinom", logfc.threshold = 0, min.pct = 0.2)
test <- FindMarkers(ild, group.by = "Diagnosis", ident.1 = "IPF", ident.2 = "Control",
                        logfc.threshold = 0)

onion <- as.num
eric(all_deg$p_val_adj)
onion[onion <= "1e-320"] <- "8.764725e-320"
#onion <- as.numeric(onion)
all_deg$p_val_adj2 <- as.numeric(onion)

plot(all_deg$avg_logFC,-log10(all_deg$p_val_adj2), col="light grey", xlab = "log2(Fold change)"
     , ylab = "-log10(p value)", xlim = c(-1.5,1.5), pch=20, cex=.7)
sel <- which(all_deg$p_val_adj <= .1)
points(all_deg[sel, "avg_logFC"], -log10(all_deg[sel, "p_val_adj"]), col="#EE7342", pch=20,cex=.7)
with(subset(all_deg, p_val_adj<=0.1 & abs(avg_logFC) >= 0.5), textxy(avg_logFC, -log10(p_val_adj), 
                            labs = row.names(all_deg[abs(all_deg$avg_logFC) >= 0.5, ]), cex=0.7, offset =1, 
                                                    pos = 3, col="#EE7342"))
write.table(all_deg, file = "190722_All_IPF_vs_control.csv", sep = ",")

###

epi_mes <- readRDS("190708_Epi_mes.rds")

epi_mes.ipf <- subset(epi_mes, cells = row.names(epi_mes@meta.data
                                                 [epi_mes@meta.data$Diagnosis == c("IPF"),]))
epi_fibro_CT <- c("HAS1 High Fibroblasts","PLIN2+ Fibroblasts","Myofibroblasts",
                  "Fibroblasts","KRT5-/KRT17+","Transitional AT2","AT2","AT1",
                  "SCGB3A2+ SCGB1A1+","SCGB3A2+","Proliferating Epithelial Cells",
                  "MUC5B+","MUC5AC+ High","Differentiating Ciliated","Ciliated","Basal")

epi_fibro <- subset(epi_mes, cells = rownames(epi_mes@meta.data[epi_mes@meta.data$celltype 
                                                                %in% epi_fibro_CT,]))
krt5_basal <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Basal", test.use = "negbinom")
krt5_at2 <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "AT2", test.use = "negbinom")
krt5_at1 <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "AT1", test.use = "negbinom")
krt5_transat2 <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Transitional AT2", test.use = "negbinom")
krt5_ciliated <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Ciliated", test.use = "negbinom")
krt5_diff_ciliated <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Differentiating Ciliated", test.use = "negbinom")
krt5_proli_epi <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Proliferating Epithelial Cells", test.use = "negbinom")
krt5_muc5b <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "MUC5B+", test.use = "negbinom")
krt5_muc5ac<- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "MUC5AC+ High", test.use = "negbinom")
krt5_scgb3a2 <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "SCGB3A2+", test.use = "negbinom")
krt5_scgb3a2_1a1<- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "SCGB3A2+ SCGB1A1+", test.use = "negbinom")
krt5_fibro <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Fibroblasts", test.use = "negbinom")
krt5_myo <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "Myofibroblasts", test.use = "negbinom")
krt5_plin2 <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "PLIN2+ Fibroblasts", test.use = "negbinom")
krt5_has1 <- FindMarkers(epi_mes.ipf, group.by = "celltype", ident.1 = "KRT5-/KRT17+", ident.2 = "HAS1 High Fibroblasts", test.use = "negbinom")

epi_fibro_deg <- list(krt5_basal,krt5_at2,krt5_at1,krt5_transat2,krt5_ciliated,krt5_diff_ciliated,krt5_proli_epi,krt5_muc5ac,krt5_muc5b,
                      krt5_scgb3a2, krt5_scgb3a2_1a1, krt5_fibro,krt5_myo,krt5_plin2,krt5_has1)
names(epi_fibro_deg) <- c("Basal","AT2","AT1","Transitional AT2","Ciliated","Differentiating Ciliated","Proliferating Epithelial",
                          "MUC5B+","MUC5AC+ High", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "Fibroblasts", "Myofibroblasts","PLIN2+ Fibroblasts",
                          "HAS1 High Fibroblasts")

onion <- lapply(epi_fibro_deg, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(epi_fibro_deg, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_krt5_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_krt5_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], 
                      onion2[[5]], onion2[[6]], onion2[[7]], onion2[[8]],
                      onion2[[9]], onion2[[10]], onion2[[11]], onion2[[12]],
                      onion2[[13]], onion2[[14]], onion2[[15]])
upset_krt5_vs_c <- as.data.frame(upset_krt5_vs_c)

row.names(upset_krt5_vs_c) <- onion
colnames(upset_krt5_vs_c) <- c("Basal","AT2","AT1","Transitional AT2","Ciliated","Differentiating Ciliated","Proliferating Epithelial",
                               "MUC5B+","MUC5AC+ High", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "Fibroblasts", "Myofibroblasts","PLIN2+ Fibroblasts",
                               "HAS1 High Fibroblasts")

upset_krt5_vs_c[upset_krt5_vs_c == T] <- 1
upset_krt5_vs_c[upset_krt5_vs_c == F] <- 0

upset(upset_krt5_vs_c, nsets = 15, text.scale = 2, show.numbers = F, keep.order = T, sets = c("Basal",
                                "AT2","AT1","Transitional AT2","Ciliated","Differentiating Ciliated","Proliferating Epithelial",
                                "MUC5B+","MUC5AC+ High", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "Fibroblasts", "Myofibroblasts","PLIN2+ Fibroblasts",
                                "HAS1 High Fibroblasts"), nintersects = 200)

# ======================================
# Figure S: 3
# ======================================
epi_fibro_deg2 <- list(krt5_basal,krt5_at2,krt5_at1,krt5_transat2,krt5_proli_epi,
                      krt5_fibro,krt5_myo,krt5_plin2,krt5_has1)
names(epi_fibro_deg2) <- c("Basal","AT2","AT1","Transitional AT2","Proliferating Epithelial",
                          "Fibroblasts", "Myofibroblasts","PLIN2+ Fibroblasts",
                          "HAS1 High Fibroblasts")

onion <- lapply(epi_fibro_deg2, function(xx){ row.names(xx[xx$p_val_adj <= .1,])})
onion <- unique(unlist(onion))
onion2 <- lapply(epi_fibro_deg2, function(xx) {onion %in% row.names(xx[xx$p_val_adj <= .1,])})
upset_krt5_vs_c <- as.data.frame(onion2, col.names = 1:length(onion2) )
upset_krt5_vs_c <- cbind(onion2[[1]], onion2[[2]], onion2[[3]], onion2[[4]], 
                         onion2[[5]], onion2[[6]], onion2[[7]], onion2[[8]],
                         onion2[[9]])
upset_krt5_vs_c <- as.data.frame(upset_krt5_vs_c)

row.names(upset_krt5_vs_c) <- onion
colnames(upset_krt5_vs_c) <- c("Basal","AT2","AT1","Transitional AT2","Proliferating Epithelial",
                               "Fibroblasts", "Myofibroblasts","PLIN2+ Fibroblasts",
                               "HAS1 High Fibroblasts")

upset_krt5_vs_c[upset_krt5_vs_c == T] <- 1
upset_krt5_vs_c[upset_krt5_vs_c == F] <- 0

upset(upset_krt5_vs_c, nsets = 9, text.scale = 2, show.numbers = F, keep.order = T, 
      sets = c("HAS1 High Fibroblasts","PLIN2+ Fibroblasts","Myofibroblasts","Fibroblasts",
               "Proliferating Epithelial","Transitional AT2","AT1","AT2","Basal"), nintersects = 90)

