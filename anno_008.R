#Annotation 
#Seurat V5 
#zheng
#Tue Aug 15 16:55:02 2023
#15 sub-libs 778,438 cells

# packages and functions ----------------------------------------------------------------
library(Seurat) #v4.9.9.9058
library(SeuratObject)
library(BPCells)
library(dplyr)

#devtools::install_github("satijalab/seurat", "seurat5")
#devtools::install_github("satijalab/seurat-data", "seurat5")
#devtools::install_github("satijalab/azimuth", "seurat5")

library(Azimuth)

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

# set this option when analyzing large datasets
options(future.globals.maxSize = 256*1024^3)
options(Seurat.object.assay.version = "v5")

future::plan(strategy = "multicore", workers = 64)

setwd("/home/big/zheng_song/aml")

wd_path <- getwd()

options(warn = 1)

# integrated Seurat object ------------------------------------------------

all_cells <- readRDS(paste0(wd_path, "/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008/aml_all_cells_cca_rpca_harmony_scvi_mrd008.Rds"))

DefaultAssay(all_cells) <- "RNA"

# split transplant ---------------------------------------------------------
#backup Peripheral Blood Stem Cell(pbsc) on disk
pbsc <- subset(all_cells, subset = transplant == "Transplant")

saveRDS(
  object = pbsc,
  file = "aml_transplant_cca_rpca_harmony_scvi_mrd008.Rds",
  destdir = paste0(wd_path, '/comb_s15/aml_transplant_cca_rpca_harmony_scvi_mrd008')
)

pbsc <- readRDS(paste0(wd_path, "/comb_s15/aml_transplant_cca_rpca_harmony_scvi_mrd008/aml_transplant_cca_rpca_harmony_scvi_mrd008.Rds"))

#backup Peripheral Blood Mononuclear Cell(pbmc) on disk
pbmc <- subset(all_cells, subset = transplant == "Non-transplant")

saveRDS(
  object = pbmc,
  file = "aml_non_transplant_cca_rpca_harmony_scvi_mrd008.Rds",
  destdir = paste0(wd_path, '/comb_s15/aml_non_transplant_cca_rpca_harmony_scvi_mrd008')
)

pbmc <- readRDS(paste0(wd_path, "/comb_s15/aml_non_transplant_cca_rpca_harmony_scvi_mrd008/aml_non_transplant_cca_rpca_harmony_scvi_mrd008.Rds"))

rm(all_cells)
rm(pbmc)
rm(pbsc)
rm(merge_cells)

# run azimuth -------------------------------------------------------------
# anno on harmnoy.full reduction ------------------------------------------
#pbmc non-transplant
pbmc[["sketch"]] <- JoinLayers(pbmc[["sketch"]])
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])

# harmnoy.full reduction
#pbmc <- RunUMAP(pbmc, reduction = "harmony.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_", assay = "RNA")

pbmc <- RunAzimuth(query = pbmc, reference = "pbmcref", verbose = T, assay = "RNA")

Feature_rast(pbmc, "predicted.celltype.l3", do.label = T, colorset = "gg")

Feature_rast(pbmc, "predicted.celltype.l2", do.label = T, colorset = "gg") %T>% figsave("/azimuth_anno/pbmc_harmony_anno.pdf", w = 300, h = 150)
Feature_rast(pbmc, "donor", colorset = "gg") %T>% figsave("/azimuth_anno/pbmc_harmony_donor.pdf", w = 300, h = 150)
Feature_rast(pbmc, "sample", colorset = "gg") %T>% figsave("/azimuth_anno/pbmc_harmony_sample.pdf", w = 300, h = 150)

#pbsc transplant
pbsc[["sketch"]] <- JoinLayers(pbsc[["sketch"]])
pbsc[["RNA"]] <- JoinLayers(pbsc[["RNA"]])

# harmnoy.full reduction
# pbsc <- RunUMAP(pbsc, reduction = "harmony.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_", assay = "RNA")

pbsc <- RunAzimuth(query = pbsc, reference = "bonemarrowref", verbose = T, assay = "RNA")

Feature_rast(pbsc, "predicted.celltype.l2", do.label = T, colorset = "gg") %T>% figsave("/azimuth_anno/pbsc_harmony_anno.pdf", w = 300, h = 150)
Feature_rast(pbsc, "donor", colorset = "gg") %T>% figsave("/azimuth_anno/pbsc_harmony_donor.pdf", w = 300, h = 150)
Feature_rast(pbsc, "sample", colorset = "gg") %T>% figsave("/azimuth_anno/pbsc_harmony_sample.pdf", w = 300, h = 150)

# anno on cca.full reduction ------------------------------------------
#pbmc non-transplant
# pbmc[["sketch"]] <- JoinLayers(pbmc[["sketch"]])
# pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])

# integrated.cca.full reduction
pbmc <- RunUMAP(pbmc, reduction = "integrated.cca.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_", assay = "RNA")

pbmc <- RunAzimuth(query = pbmc, reference = "pbmcref", verbose = T, assay = "RNA")

Feature_rast(pbmc, "predicted.celltype.l2", do.label = T, colorset = "gg") %T>% figsave("/azimuth_anno/pbmc_cca_anno.pdf", w = 300, h = 150)
Feature_rast(pbmc, "predicted.celltype.l2", do.label = T, colorset = "gg", facets = "tp") %T>% figsave("/azimuth_anno/pbmc_cca_anno_tp.pdf", w = 600, h = 150)
Feature_rast(pbmc, "donor", colorset = "gg") %T>% figsave("/azimuth_anno/pbmc_cca_donor.pdf", w = 300, h = 150)
Feature_rast(pbmc, "sample", colorset = "gg") %T>% figsave("/azimuth_anno/pbmc_cca_sample.pdf", w = 300, h = 150)

#pbsc transplant
# pbsc[["sketch"]] <- JoinLayers(pbsc[["sketch"]])
# pbsc[["RNA"]] <- JoinLayers(pbsc[["RNA"]])

# integrated.cca.full reduction
pbsc <- RunUMAP(pbsc, reduction = "integrated.cca.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_", assay = "RNA")

pbsc <- RunAzimuth(query = pbsc, reference = "bonemarrowref", verbose = T, assay = "RNA")

Feature_rast(pbsc, "predicted.celltype.l2", do.label = T, colorset = "gg") %T>% figsave("/azimuth_anno/pbsc_cca_anno.pdf", w = 300, h = 150)
Feature_rast(pbsc, "donor", colorset = "gg") %T>% figsave("/azimuth_anno/pbsc_cca_donor.pdf", w = 300, h = 150)
Feature_rast(pbsc, "sample", colorset = "gg") %T>% figsave("/azimuth_anno/pbsc_cca_sample.pdf", w = 300, h = 150)

# collect cell type label ------------------------------------------------
pbmc_type <- pbmc@meta.data %>% select(predicted.celltype.l2, predicted.celltype.l2.score)

pbsc_type <- pbsc@meta.data %>% select(predicted.celltype.l2, predicted.celltype.l2.score)

all_cells_type <- rbind(pbmc_type, pbsc_type)

all_cells_type$bc <- rownames(all_cells_type)

# introduce cell type into all cells
all_cells <- readRDS(paste0(wd_path, "/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008/aml_all_cells_cca_rpca_harmony_scvi_mrd008.Rds"))

all_cells@meta.data %>% colnames()

all_cells$cluster_full <- NULL
all_cells$clusters_full <- NULL
all_cells$cluster_full.score <- NULL
all_cells$clusters_full.score <- NULL

all_cells$bc <- rownames(all_cells@meta.data)

# all_cells_type from pbmc and pbsc which are subset from all cells by transplant
# all_cells@meta.data %<>% left_join(all_cells_type, by = "bc")
# 
# rownames(all_cells@meta.data) <- all_cells@meta.data$bc
# 
# all_cells@meta.data %>% colnames()
# 
# DefaultAssay(all_cells) <- "RNA"

#umap.full on harmony
# Feature_rast5(all_cells, "harmony_clusters_full", colorset = "gg")
# 
# Feature_rast5(all_cells, "predicted.celltype.l2", colorset = "gg")
# 
#umap.full on cca
# all_cells <- RunUMAP(all_cells, reduction = "integrated.cca.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_")
# 
# Feature_rast5(all_cells, "cca_clusters_full", colorset = "gg")
# 
# Feature_rast5(all_cells, "predicted.celltype.l2", colorset = "gg")

# runazimuth on all cells
# reference on bonemarrowref
all_cells[["RNA"]] <- JoinLayers(all_cells[["RNA"]])

all_cells <- RunAzimuth(query = all_cells, reference = "pbmcref", verbose = T, assay = "RNA")

(Feature_rast5(all_cells, "predicted.celltype.l2", do.label = T, colorset = "gg")) %T>% figsave("/azimuth_anno/all_cells_bonemarrowref_anno_l2.pdf", w = 300, h = 150)
(Feature_rast5(all_cells, "predicted.celltype.l1", do.label = T, colorset = "gg")) %T>% figsave("/azimuth_anno/all_cells_bonemarrowref_anno_l1.pdf", w = 300, h = 150)

Feature_rast5(all_cells, "donor", colorset = "gg") %T>% figsave("/azimuth_anno/all_cells_bonemarrowref_donor.pdf", w = 300, h = 150)
Feature_rast5(all_cells, "sample", colorset = "gg") %T>% figsave("/azimuth_anno/all_cells_bonemarrowref_sample.pdf", w = 300, h = 150)
Feature_rast5(all_cells, "TRDV2", colorset = "gg", color_grd = "grd") %T>% figsave("/azimuth_anno/all_cells_bonemarrowref_trdv2.pdf", w = 300, h = 150)

saveRDS(
  object = all_cells,
  file = "aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno.Rds",
  destdir = paste0(wd_path, '/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno')
)

all_cells <- readRDS(paste0(wd_path, "/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno.Rds"))

cca_umap_full <- readRDS(paste0(wd_path, "/aml_all_cells_cca_mrd008_umapfull_reduction.Rds"))

all_cells@reductions$umap.full <- cca_umap_full

Feature_rast5(all_cells, "predicted.celltype.l2", do.label = T, colorset = "gg")

Feature_rast5(all_cells, "cca_clusters_full", do.label = T, colorset = "gg")

Feature_rast5(all_cells, do.label = T, colorset = "gg")

all_cells <- BuildClusterTree(all_cells, assay = "RNA", dims = 1:14, reduction = "integrated.cca.full", reorder = T, verbose = T)

PlotClusterTree(object = all_cells)

hspc <- subset(all_cells, predicted.celltype.l1 == "HSPC")

Feature_rast5(hspc, "predicted.celltype.l2", do.label = T, colorset = "gg")

all_cells$bm_ref_celltype.l2 <- all_cells$predicted.celltype.l2

all_cells$bm_ref_celltype.l1 <- all_cells$predicted.celltype.l1

all_cells <- RunAzimuth(query = all_cells, reference = "pbmcref", verbose = T, assay = "RNA", min.dist = 0.8)

Feature_rast5(all_cells, "predicted.celltype.l1", do.label = T, colorset = "gg", sz = 0.1)

all_cells <- RunUMAP(all_cells, reduction = "integrated_dr", min.dist = 0.8, dims = 1:50, reduction.name = "refumap",
                     reduction.key = "UMAP_")

Feature_rast(all_cells, "predicted.celltype.l2", do.label = T, colorset = "gg", sz = 0.1)

Feature_rast5(all_cells, "predicted.celltype.l2", do.label = T, colorset = "gg", sz = 0.1, facets = "transplant")

Feature_rast5(all_cells, "bm_ref_celltype.l2", do.label = T, colorset = "gg", sz = 0.1, facets = "transplant", color_grd = "grd")

Feature_rast5(all_cells, "CD36", do.label = T, colorset = "gg", sz = 0.1, facets = "transplant", color_grd = "grd")

Feature_rast5(all_cells, "VCAM1", do.label = T, colorset = "gg", sz = 0.1, facets = "transplant", color_grd = "grd")

Feature_rast5(all_cells, "CDH5", do.label = T, colorset = "gg", sz = 0.1, facets = "tp", color_grd = "grd")

Feature_rast5(all_cells, "bm_ref_celltype.l2", do.label = T, colorset = "gg", sz = 0.1, facets = "tp", color_grd = "grd")

hspc_bmref <- subset(all_cells, bm_ref_celltype.l2 %in% c("CLP", "Early Eryth", "EMP", "GMP", "HSC", "Late Eryth", "pre B", "pre-pDC", "pro B", "Prog Mk", "Stromal"))

hspc <- RunAzimuth(query = hspc, reference = "bonemarrowref", verbose = T, assay = "RNA")

DimPlot(all_cells, group.by = "bm_ref_celltype.l1", repel = T, alpha = 0.1, reduction = "umap.full", cols = c())

DimPlot(all_cells, group.by = "bm_ref_celltype.l1", repel = T, alpha = 0.1, reduction = "umap.full", cols = c())

FeaturePlot(all_cells, features = "percent.mt", alpha = 0.1, reduction = "umap.full")

FeaturePlot(all_cells, features = "CSF1R", alpha = 0.1, reduction = "umap.full")

DimPlot(subset(all_cells, sublib == "s1"), group.by = "predicted.celltype.l2", repel = T, alpha = 0.05, reduction = "umap.full", raster=FALSE, label = F, )

DimPlot(all_cells, group.by = "predicted.celltype.l2", repel = T, alpha = 0.05, reduction = "umap.full", raster=FALSE, label = F, )

DimPlot(subset(all_cells, predicted.celltype.l2 == "gdT"), group.by = "predicted.celltype.l2", repel = T, alpha = 1, reduction = "umap.full", raster=FALSE, label = F, )

all_cells

markers <- FindMarkers(all_cells, ident.1 = "35", ident.2 = "38", group.by = "cca_clusters_full", assay = "RNA")

rm(hspc)
rm(hspc_bmref)

all_cells$pbmc_ref_celltype.l3 <- all_cells$predicted.celltype.l3

all_cells$pbmc_ref_celltype.l2 <- all_cells$predicted.celltype.l2

all_cells$pbmc_ref_celltype.l1 <- all_cells$predicted.celltype.l1

all_cells <- RunAzimuth(query = all_cells, reference = "bonemarrowref", verbose = T, assay = "RNA")

DimPlot(subset(all_cells, sublib == "s1"), group.by = "predicted.celltype.l2", repel = T, alpha = 0.05, reduction = "umap.full", raster=FALSE, label = F, )

DimPlot(all_cells, group.by = "predicted.celltype.l2", repel = T, alpha = 0.05, reduction = "umap.full", raster=FALSE, label = F, )

# %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10


DoHeatmap(object = all_cells, )


ClusterCompare(all_cells, id1 = "35", id2 = "38", do.plot = T)

# 
 # '#456c35', '#778f51', '#acc38b', 
 # '#494b64', '#2f1c94','#036fa0', '#80a3db', '#bfb5f5', '#e3e0f3',
 # '#b58824', '#f9ae17', '#ffd16c', '#cfc94e', '#fffa90', '#fdffcc', 
 # '#751a2c', '#b31e3c','#b33a3a', '#ab3a18', '#f2744e', '#ed896f'
# 

pbmc_small <- BuildClusterTree(object = pbmc_small, )
PlotClusterTree(object = pbmc_small)
