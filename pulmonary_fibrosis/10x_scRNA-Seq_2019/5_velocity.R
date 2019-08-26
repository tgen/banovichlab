# ==============================================================================
# Author(s) : Austin Gutierrez, agutierrez@tgen.org
# Date : 02/07/19
# Description: Script for RNA velocity analysis from the veloyto-team
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
# Full epi set
# ======================================
epi_emat <- filter.genes.by.cluster.expression(epi_emat,epi_color,min.max.cluster.average = .1)
epi_nmat <- filter.genes.by.cluster.expression(epi_nmat,epi_color,min.max.cluster.average = .05)
epi_smat <- filter.genes.by.cluster.expression(epi_smat,epi_color,min.max.cluster.average = 0.001)

length(intersect(rownames(epi_emat),rownames(epi_nmat)))
length(intersect(intersect(rownames(epi_emat),rownames(epi_nmat)), rownames(epi_smat)))

arrow.scale = 3
cell.alpha = 0.4
cell.cex = 1
fig.height = 4
fig.width = 4.5
fit.quantile <- 0.01

rvel <- gene.relative.velocity.estimates(emat = epi_emat,
              nmat = epi_nmat, smat = epi_smat, kCells = 50,
              fit.quantile = fit.quantile, diagonal.quantiles = TRUE,
              n.cores = 40)

#saveRDS(rvel, "rvel.rds")

rvel <- readRDS("rvel.rds")

show.velocity.on.embedding.cor(epi_emb, rvel, n = 50, 
                               scale = 'sqrt', cell.colors = epi_color,
                               cex = cell.cex, arrow.scale
                               = arrow.scale, show.grid.flow = T, 
                               min.grid.cell.mass = 1, grid.n = 40,
                               arrow.lwd = 2, main = c("Epi subset_", 40),
                               n.cores = 40)

# ======================================
# 
# ======================================
sub_1 <- readRDS("/scratch/lbui/20190623_Final_version/190719_krt5_trans.rds")

sub_1_emb <- sub_1@reductions$umap@cell.embeddings
sub_1_emat <- epi_emat[,colnames(epi_emat) %in% rownames(sub_1@meta.data)]
sub_1_nmat <- epi_nmat[,colnames(epi_nmat) %in% rownames(sub_1@meta.data)]
sub_1_smat <- epi_smat[,colnames(epi_smat) %in% rownames(sub_1@meta.data)]

sub_1_emat <- filter.genes.by.cluster.expression(sub_1_emat, epi_color, min.max.cluster.average = .1)
sub_1_nmat <- filter.genes.by.cluster.expression(sub_1_nmat, epi_color, min.max.cluster.average = .05)
sub_1_smat <- filter.genes.by.cluster.expression(sub_1_smat, epi_color, min.max.cluster.average = 0.001)

length(intersect(rownames(sub_1_emat),rownames(sub_1_nmat)))
length(intersect(intersect(rownames(sub_1_emat),rownames(sub_1_nmat)), rownames(sub_1_smat)))

arrow.scale <- 3
cell.alpha <- 0.4
cell.cex <- 1
fig.height <- 4
fig.width <- 4.5
fit.quantile <- 0.01

sub_1_rvel <- gene.relative.velocity.estimates(emat = sub_1_emat,
                                         nmat = sub_1_nmat,
                                         smat = sub_1_smat, kCells = 50,
                                         fit.quantile = fit.quantile, diagonal.quantiles = TRUE,
                                         n.cores = 40)

#saveRDS(sub_1_rvel, "sub_1_rvel.rds")

show.velocity.on.embedding.cor(sub_1_emb, sub_1_rvel, n = 100, 
                               scale = 'sqrt', cell.colors = epi_color,
                               cex = cell.cex, arrow.scale
                               = 1, show.grid.flow = T, 
                               min.grid.cell.mass = 1, grid.n = 50,
                               arrow.lwd = 2, main = c("5 pop velocity"),
                               n.cores = 40)

# Trajectory modeling
test <- show.velocity.on.embedding.eu(sub_1_emb, sub_1_rvel, n = 50, scale = 'sqrt', 
                              cell.colors = epi_color, cex = cell.cex,
                              nPcs = 30, sigma = 2.5, show.trajectories = TRUE,
                              diffusion.steps = 400, n.trajectory.clusters = 15
                              , ntop.trajectories = 1, embedding.knn = T,
                              control.for.neighborhood.density = TRUE,
                              n.cores = 40) 
