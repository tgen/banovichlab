# =====================================
# Function for collapsing list of plots
# =====================================
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

SOX4 <- c("AREG", "ZMAT3", "PMEPA1", "TPM1")
SOX9 <- c("AREG", "ZMAT3", "PMEPA1")
NR1D1 <- c("AREG", "ZMAT3", "PMEPA1", "CRIP1", "TMP1")

motifs <- c("SOX4", "SOX9")
gene_list <- c("CDKN1A",
               "CRIP1",
               "AREG",
               "AHNAK2",
               "CLDN4",
               "ZMAT3",
               "CMTM3",
               "IL32",
               "MARCKS",
               "TMEM59L",
               "TPM1",
               "PMEPA1",
               "HOMER3",
               "TMEM132A",
               "CST6",
               "FSTL1",
               "COL1A1",
               "CDKN2A",
               "PHACTR3",
               "MMP14",
               "STMN1",
               "ABRACL",
               "CMTM7",
               "MYL9",
               "RHOD",
               "TUBB",
               "GAS6",
               "CDH1",
               "MACC1",
               "CALM2",
               "KRT19",
               "TAX1BP3",
               "TRAM1",
               "MDK")

# You will need to subset the Epithelial and Mesenchymal populations from the ILD object.
epi_mesen <- readRDS(file = "epi_mesenchymal.rds")
onion <- union(rownames(newest@meta.data[newest@meta.data$population == "Epithelial", ]),
               rownames(newest@meta.data[newest@meta.data$population == "Mesenchymal",  ]))

sub_1 <- subset(newest, cells = onion)

KRT5_sub <- subset(sub_1, cells = rownames(sub_1@meta.data[sub_1@meta.data$celltype == "KRT5-/KRT17+", ]))
Trans_AT2_sub <- subset(sub_1, cells = rownames(sub_1@meta.data[sub_1@meta.data$celltype == "Transitional AT2", ]))
rm(sub_1, newest)

plot_list1 <- list()
plot_list2 <- list()
plot_index <- 1
for (j in 1:length(motifs)) {
  
  for (i in 1:length(gene_list)) {
    onion = KRT5_sub@assays$SCT@data[motifs[j] , ]
    onion2 = KRT5_sub@assays$SCT@data[gene_list[i], ]
    
    onion = as.data.frame(onion)
    onion2 = as.data.frame(onion2)
    
    onion$ind = do.call(rbind, strsplit(rownames(onion), split = "_"))[,1]
    onion2$ind = do.call(rbind, strsplit(rownames(onion2), split = "_"))[,1]
    
    onion$Status = sub_1@meta.data[rownames(onion),]$Status
    onion2$Status = sub_1@meta.data[rownames(onion2),]$Status
    
    onion_means = onion %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion, na.rm = T))
    onion2_means = onion2 %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
    
    onion_means = as.data.frame(onion_means)
    onion_means = onion_means[onion_means$Mean > 0,]
    
    onion2_means = as.data.frame(onion2_means)
    onion2_means = onion2_means[onion2_means$Mean > 0,]
    
    means = merge(onion_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
    means = means[,-4]
    
    fit <- lm(means$Mean.y ~ means$Mean.x)
    pVal <- round(anova(fit)$'Pr(>F)'[1], 4)
    
    if (pVal < .05) {
      plot_list1[[plot_index]] <- ggplot(means, aes(x = Mean.x, y = Mean.y)) +
      geom_point(aes(color = Status.x), size = .5) + 
      ggtitle(paste("KRT5-/KRT17+ Cells ", pVal, sep = "")) + 
      xlab(motifs[j]) + 
      ylab(gene_list[i]) + 
      theme(plot.title = element_text(size = 5),
      legend.title = element_blank(),
      axis.title=element_text(size=5), 
      axis.text = element_text(size = 5),
      legend.text = element_text(size = 3)) +
      geom_smooth(method = lm, se = TRUE) + 
      NoLegend()
      
      plot_index <- plot_index + 1
    }
  
    onion = Trans_AT2_sub@assays$SCT@data[motifs[j] , ]
    onion2 = Trans_AT2_sub@assays$SCT@data[gene_list[i], ]
    
    onion = as.data.frame(onion)
    onion2 = as.data.frame(onion2)
    
    onion$ind = do.call(rbind, strsplit(rownames(onion), split = "_"))[,1]
    onion2$ind = do.call(rbind, strsplit(rownames(onion2), split = "_"))[,1]
    
    onion$Status = sub_1@meta.data[rownames(onion),]$Status
    onion2$Status = sub_1@meta.data[rownames(onion2),]$Status
    
    onion_means = onion %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion, na.rm = T))
    onion2_means = onion2 %>% group_by(ind, Status) %>% dplyr::summarise(Mean = mean(onion2, na.rm = T))
    
    onion_means = as.data.frame(onion_means)
    onion_means = onion_means[onion_means$Mean > 0,]
    
    onion2_means = as.data.frame(onion2_means)
    onion2_means = onion2_means[onion2_means$Mean > 0,]
    
    means = merge(onion_means, onion2_means, by.x = "ind", by.y = "ind", all.x = F, all.y = F)
    means = means[,-4]
    
    fit <- lm(means$Mean.y ~ means$Mean.x)
    pVal <- round(anova(fit)$'Pr(>F)'[1], 4)
    
    if (pVal < .05) {
      plot_list2[[plot_index]] <- ggplot(means, aes(x = Mean.x, y = Mean.y)) + 
      geom_point(aes(color = Status.x), size = .5) + 
      ggtitle(paste("Transitional AT2 Cells ", pVal, sep = "")) + 
      xlab(motifs[j]) + 
      ylab(gene_list[i]) + 
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
plot_list <- interleave(plot_list1, plot_list2)
plot_list <- plyr::compact(plot_list)
pdf(file = paste("all_sig_genes.pdf", sep = "_"), width = 8.5, height = 11)
#grid.arrange(grobs = plot_list, ncol = 2)

grid.arrange(grobs = plot_list, ncol = 5)
dev.off()






