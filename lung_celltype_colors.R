library(RColorBrewer)
library(ggplot2)
library(nord)
library(circlize)

# Theme for plotting
my_theme <- theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  legend.position = "none",
                  plot.title = element_text(size = 12))

# Colors
epi_col <- c("AT1" = "#548BC5",
             "AT2" = "#EE7342",
             "Proliferating" = "#B874FF",
             "Transitional AT2" = "#F1A5C4",
             "AT2 - low quality" = "#2FC895",
             "Basal" = "#D38402",
             "Ciliated" = "#B19302",
             "PNEC" = "#FFA19D",
             "Differentiating Ciliated" = "#9C9966",
             "KRT5-/KRT17+" = "#03AF21",
             "Secretory - MUC5AC+" = "#003366",
             "Secretory - SCGB1A1+/MUC5B+" = "#00B2DB",
             "Secretory - SCGB3A1+/MUC5B+" = "#F659DD",
             "Secretory - SCGB1A1+/SCGB3A2+" = "#9900CC",
             "Secretory - SCGB3A2+" = "#FF5BC8")

immune_col <- c("Proliferating" = "#009FFA",
                "Monocyte-derived macrophage" = "#FF6B65",
                "moDC" = "#A39922",
                "Inflammatory monocyte" = "#CC8B22",
                "Alveolar macrophage" = "#EB7B27",
                "cDC2" = "#FF57C3",
                "cDC1" = "#00B7BB",
                "Macrophage - SPP1+" = "#00AE21",
                "Monocyte" = "#00AFE0",
                "CD4" = "#6CA421",
                "NK" = "#00B990",
                "Mast" = "#C471FA",
                "pDC" = "#FF5D97",
                "CD8/NKT" = "#00B560",
                "Plasma" = "#F35EE6",
                "B cells" = "#748AFA")

mesen_col <- c("MyoFB" = "#00B13B",
               "HAS1 High Activated FB" = "#FF6B65",
               "SMC" = "#DC68F6",
               "MyoFB - Activated" = "#00B995",
               "Matrix FB" ="#D38721",
               "PLIN2+ FB" = "#4A92FA",
               "Pericyte" = "#00B0DC",
               "Mesothelial" = "#87A022",
               "WNT2+ FB" = "#FF57B9")

endo_col <- c("Endothelial - capillary" = "#FF6B65",
              "Endothelial - venule" = "#00B7BB",
              "Endothelial - arteriole" = "#00ADE5",
              "Lymphatic" = "#009AFA", 
              "Endothelial - inflamed" = "#9B7FFA",
              "Endothelial - CA4+ capillary" = "#FF57CF", 
              "Endothelial - peribronchiolar" = "#FF5C9E")

simple_col <- c("Vascular endothelial" = "#FF6B65",
                "Lymphatic endothelial" = "#009AFA",
                "Secretory" = "#00B2DB")

lung_celltype_cols <- c(epi_col, immune_col, mesen_col, endo_col, simple_col)

# Heatmap colors
col_fun <- colorRamp2(c(-2, 0, 2), c("cadetblue4", "white", "coral2"))

# Flowcell_ID
flowcells <- c("H75C7DRXX",
               "HN3KNDMXX",
               "HNM3GDRXX",
               "HM7W7DRXX",
               "HLKH5DSXY",
               "HV27JDSXY",
               "H3HLJDSX2-H5W7FDSX2",
               "HWYTFBBXX",
               "HC57FDRXX",
               "HM7F5DRXX",
               "HFJKFDSXY",
               "H5LLNDSXX",
               "H5LLFDSXX",
               "H5LLJDSXX",
               "H5LM5DSXX",
               "HCKWNDSXX",
               "HMWLCBGX7",
               "H5LLGDSXX",
               "H5LM3DSXX",
               "HC73FDSXX",
               "HCCM2DSXX",
               "HH7HWDSXX",
               "H5NTYDSXY",
               "H3G2VDSX2")
flowcell_col <- colorRampPalette(brewer.pal(10, "Spectral"))(nb.cols <- length(flowcells))
names(flowcell_col) <- flowcells
colScale_flowcell <- scale_fill_manual(name = "Flowcell_ID", values = flowcell_col)

# Cell cycle
cell_cycle_col <- c("G1" = nord("lumina", 5)[1],
                    "G2M" = nord("lumina", 5)[3],
                    "S" = nord("lumina", 5)[5])

# Cluster
cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- 16)
names(cluster_col) <- c(0, seq(1, 15))
colScale_cluster <- scale_fill_manual(name = "Cluster", values = cluster_col)

ild_col <- list("cluster_col" = cluster_col,
                "flowcell_col" = flowcell_col,
                "cell_cycle_col" = cell_cycle_col)


