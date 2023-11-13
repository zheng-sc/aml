#Sketch integration 
#Seurat V5 
#zheng
#Fri Jul 28 11:07:15 2023
#15 sub-libs 778,438 cells
#aml_all_cells_rpca_mrd007

# packages and functions ----------------------------------------------------------------
library(Seurat) #v4.9.9.9058
library(SeuratObject)
library(BPCells)
library(dplyr)

library(cowplot) # ggsave2
library(ggrastr) # geom_point_rast
library(ggplot2)
library(ggrepel)
library(patchwork)

library(stringr) # string manipulation
library(magrittr) # %>% 

library(future)
# functions
source('/home/big/zheng/rscripts/funcs_umapfull.r')
source('/home/big/zheng/rscripts/themes.R')

file.edit('/home/big/zheng/rscripts/funcs_umapfull.r')

# set this option when analyzing large datasets
options(future.globals.maxSize = 256*1024^3)
options(Seurat.object.assay.version = "v5")

future::plan(strategy = "multicore", workers = 32)

setwd("/home/big/zheng_song/aml")

utils::methods(class = 'Seurat')
# annotation MRD007--------------------------------------------------------------
m7_all_cells <- readRDS(paste0(wd_path, '/comb_s15/aml_all_cells_rpca_mrd007/aml_all_cells_rpca_mrd007.Rds'))

DefaultAssay(m7_all_cells) <- 'RNA'
umap <- Feature_rast5(m7_all_cells, colorset = 'gg', 'cluster_full')
umap

p1 <- DimPlot(m7_all_cells, reduction = "umap.full", group.by = "donor", alpha = 0.1)
p2 <- DimPlot(m7_all_cells, reduction = "umap.full", group.by = "transplant", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 2)


umap_transplant7 <- Feature_rast5(m7_all_cells, colorset = 'gg', 'transplant', facets = 'transplant', sz = 0.1)
umap_transplant7

m8_all_cells <- readRDS(paste0(wd_path, '/comb_s15/aml_all_cells_rpca_mrd008/aml_all_cells_rpca_mrd008.Rds'))
umap8 <- Feature_rast5(m8_all_cells, colorset = 'gg', 'cluster_full')
umap8

p1 <- DimPlot(m8_all_cells, reduction = "umap.full", group.by = "donor", alpha = 0.1)
p2 <- DimPlot(m8_all_cells, reduction = "umap.full", group.by = "transplant", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 2)

umap_transplant8 <- Feature_rast5(m8_all_cells, colorset = 'gg', 'transplant', sz = 0.1)
umap_transplant8

umap_donor <- Feature_rast5(m8_all_cells, colorset = 'gg', 'donor', sz = 0.1)
umap_donor

m8_all_cells$
