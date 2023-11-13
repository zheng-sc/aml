#Sketch integration 
#Seurat V5 
#zheng
#Wed Jul 26 11:05:04 2023
#15 sub-libs 778,438 cells


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
source('/home/big/zheng/rscripts/funcs.r')
source('/home/big/zheng/rscripts/themes.R')

# set this option when analyzing large datasets
options(future.globals.maxSize = 256*1024^3)
options(Seurat.object.assay.version = "v5")

future::plan(strategy = "multicore", workers = 32)

setwd("/home/big/zheng_song/aml")
wd_path <- "/home/big/zheng_song/aml"
mat_path <- "/home/big/zheng/parse_bio/aml/novaseq/analysis/comb_s15/all-well/DGE_filtered"

# read matrix data --------------------------------------------------------
all_mat_h5ad <- BPCells::open_matrix_anndata_hdf5(paste0(mat_path, '/comb_s15_raw.h5ad'))

# Write the matrix to a directory
BPCells::write_matrix_dir(
  mat = all_mat_h5ad,
  dir = paste0(mat_path, '/bpcell_mat')
)

#read mat which is the output of BPCell
all_mat <- BPCells::open_matrix_dir(dir =  paste0(mat_path, '/bpcell_mat'))

# read in cell metadata
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

cell_meta$donor <- str_extract(cell_meta$sample, "MRD\\d\\d\\d")
cell_meta %<>% dplyr::relocate(donor, .before = sample) 

cell_meta$transplant <- str_extract(cell_meta$sample, "Transplant")
cell_meta %<>% dplyr::relocate(transplant, .before = species) 

cell_meta$tp <- str_extract(cell_meta$sample, "Transplant|D30|D100|Before_SCT")
cell_meta %<>% dplyr::relocate(tp, .before = species) 

# create seurat object
all_cells <- Seurat::CreateSeuratObject(counts = all_mat, meta.data = cell_meta)

#setting the initial cell class into a single type to avoid bias
all_cells@meta.data$orig.ident <- factor(rep("pbmc", nrow(all_cells@meta.data)))
Idents(all_cells) <- all_cells@meta.data$orig.ident


# QC ----------------------------------------------------------------------
# cell quality control
# befroe QC 778438  cells
all_cells@meta.data$percent_mt <- PercentageFeatureSet(all_cells, pattern = "^MT-")

VlnPlot(all_cells, pt.size = 0.02, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3) %>% 
  figsave("comb_s15/before_qc.pdf", w = 300, h = 100)

Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "nFeature_RNA", g = "percent_mt", color_grd = "grd", noaxis = F, axis.number = T) %>% 
  figsave("comb_s15/before_qc_scatter.pdf", w = 120, h = 100)
Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "percent_mt", g = "nFeature_RNA", color_grd = "grd", noaxis = F, axis.number = T) %>% 
  figsave("comb_s15/before_qc_scatter_2.pdf", w = 120, h = 100)

all_cells %<>% subset(nFeature_RNA < 10000 & nCount_RNA < 100000 & percent_mt < 50)
# first filter 778371 cells
VlnPlot(all_cells, pt.size = 0.02, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3) %>% 
  figsave("comb_s15/1st_before_qc.pdf", w = 300, h = 100)
Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "nFeature_RNA", g = "percent_mt", color_grd = "grd", noaxis = F, axis.number = T) %>% 
  figsave("comb_s15/1st_before_qc_scatter.pdf", w = 120, h = 100)
Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "percent_mt", g = "nFeature_RNA", color_grd = "grd", noaxis = F, axis.number = T) %>% 
  figsave("comb_s15/1st_before_qc_scatter_2.pdf", w = 120, h = 100)

all_cells %<>% subset(nFeature_RNA < 8500 & nCount_RNA < 60000 & percent_mt < 25)
# second filter 778071 cells 
VlnPlot(all_cells, pt.size = 0.02, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3) %>% 
  figsave("comb_s15/2nd_after_qc.pdf", w = 300, h = 100)
Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "nFeature_RNA", g = "percent_mt", color_grd = "grd", noaxis = F, axis.number = T) %>% 
  figsave("comb_s15/2nd_after_qc_scatter.pdf", w = 120, h = 100)
Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "percent_mt", g = "nFeature_RNA", color_grd = "grd", noaxis = F, axis.number = T) %>% 
  figsave("comb_s15/2nd_after_qc_scatter_2.pdf", w = 120, h = 100)


#normalize
all_cells <- NormalizeData(all_cells)

# split assay into 24 layers
all_cells[["RNA"]] <- split(all_cells[["RNA"]], f = all_cells$donor)

all_cells <- FindVariableFeatures(all_cells, verbose = TRUE)

# switch to analyzing the full dataset (on-disk)
DefaultAssay(all_cells) <- "RNA"
# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(all_cells) <- "sketch"

#Error in gzfile(file, mode) : cannot open the connection
saveRDS(all_cells, '/comb_s15/aml_all_no_int.rds')

# Sample representative cells from each dataset and integration sketched data---------
#sketch 10k cells from total dataset
all_cells <- SketchData(all_cells, assay = "RNA", ncells = 10000, sketched.assay = "sketch", method = "LeverageScore")

DefaultAssay(all_cells) <- "sketch"
all_cells %<>%  FindVariableFeatures(verbose = T) %>% ScaleData(verbose = T) %>% RunPCA(verbose = T)

# integrate the datasets
all_cells <- IntegrateLayers(all_cells, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                          dims = 1:30, k.anchor = 20, reference = which(Layers(all_cells, search = "data") %in% c("data.MRD007")),
                          verbose = F)
# cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "integrated.rpca", dims = 1:30)
all_cells <- FindClusters(all_cells, resolution = 2)

all_cells <- RunUMAP(all_cells, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = T)

#rejoin the layers in the sketched assay for DEG
all_cells[["sketch"]] <- JoinLayers(all_cells[["sketch"]])
c10_markers <- FindMarkers(object = all_cells, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
head(c10_markers)

plot.s1 <- DimPlot(all_cells, group.by = "donor", reduction = "umap")

# integrate the full datasets ---------------------------------------------
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells

all_cells[["sketch"]] <- split(all_cells[["sketch"]], f = all_cells$donor)

all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")

#all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
#                      full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "celltype.manual"))

all_cells <- RunUMAP(all_cells, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
                  reduction.key = "UMAP_full_")

p1 <- DimPlot(all_cells, reduction = "umap.full", group.by = "donor", alpha = 0.1)
p2 <- DimPlot(all_cells, reduction = "umap.full", group.by = "transplant", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 1)

DefaultAssay(all_cells) <- "sketch"
DefaultAssay(all_cells) <- "RNA"

DefaultAssay(all_cells) <- "RNA"
FeaturePlot(all_cells, c("CD14", 'CSF1R'))

VlnPlot(all_cells, c("CD14", 'CSF1R'))

all_cells_007 <- readRDS("/home/big/zheng/parse_bio/aml/novaseq/analysis/comb_s15/all-well/DGE_filtered/aml_all_cells_rpca_mrd007.Rds")

saveRDS(
  object = all_cells_007,
  file = "aml_all_cells_rpca_mrd007.Rds",
  destdir = paste0(wd_path, '/comb_s15/aml_all_cells_rpca_mrd007')
)

all_cells_007 <- readRDS(paste0(wd_path, '/comb_s15/aml_all_cells_rpca_mrd007/aml_all_cells_rpca_mrd007.Rds'))
