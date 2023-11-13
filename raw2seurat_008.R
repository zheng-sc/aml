#Sketch integration 
#Seurat V5 
#zheng
#Fri Jul 28 15:38:03 2023
#15 sub-libs 778,438 cells

# packages and functions ----------------------------------------------------------------
library(Seurat) #v4.9.9.9058
library(SeuratObject)
library(BPCells)
library(dplyr)
library(purrr)

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
wd_path <- getwd()
mat_path <- "/home/big/zheng/parse_bio/aml/novaseq/analysis/comb_s15/all-well/DGE_filtered"

# read matrix data --------------------------------------------------------
all_mat_h5ad <- BPCells::open_matrix_anndata_hdf5(paste0(wd_path, '/comb_s15/comb_s15_raw.h5ad'))

# Write the matrix to a directory
BPCells::write_matrix_dir(
  mat = all_mat_h5ad,
  dir = paste0(wd_path, '/comb_s15/bpcell_mat')
)

#read mat which is the output of BPCell
all_mat <- BPCells::open_matrix_dir(dir =  paste0(wd_path, '/comb_s15/bpcell_mat'))

# read in cell metadata
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

cell_meta$donor <- str_extract(cell_meta$sample, "MRD\\d\\d\\d")
cell_meta %<>% dplyr::relocate(donor, .before = sample) 

cell_meta$transplant <- str_extract(cell_meta$sample, "Transplant")

cell_meta$transplant %<>% tidyr::replace_na("Non-transplant")

cell_meta$transplant %<>% factor(levels = c("Transplant", "Non-transplant"))

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

all_cells <- FindVariableFeatures(all_cells, verbose = TRUE, selection.method = 'vst')

# switch to analyzing the full dataset (on-disk)
DefaultAssay(all_cells) <- "RNA"
# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(all_cells) <- "sketch"


saveRDS(object = all_cells, file = 'aml_all_cells_no_int.rds', destdir = paste0(wd_path, '/comb_s15/aml_all_cells_no_int'))

all_cells <- readRDS("/home/big/zheng_song/aml/comb_s15/aml_all_cells_no_int/aml_all_cells_no_int.rds")

#Sample representative cells from each dataset and integration sketched data---------
#sketch 10k cells from total dataset
all_cells <- SketchData(all_cells, assay = "RNA", ncells = 10000, sketched.assay = "sketch", method = "LeverageScore")

DefaultAssay(all_cells) <- "sketch"
all_cells %<>%  FindVariableFeatures(verbose = T, selection.method = 'vst', nfeatures = 2000)

top10 <- head(VariableFeatures(all_cells), 10)
plot1 <- VariableFeaturePlot(all_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
(plot1 + plot2) %>% figsave('/comb_s15/variableFeatures.pdf', w = 400, h = 100, dpi = 300)

all_cells %<>% ScaleData(verbose = T) %>% RunPCA(verbose = T)

VizDimLoadings(all_cells, dims = 1:10, reduction = "pca")

Seurat::ElbowPlot(all_cells, ndims = 30, reduction = 'pca') %T>% figsave('/comb_s15/elbowplot_pca.pdf', w = 100, h = 100, dpi = 300)

all_cells %<>% JackStraw(num.replicate = 100, dims = 30) %>% ScoreJackStraw(dims = 1:30)

JackStrawPlot(all_cells, dims = 1:30) %T>% figsave('/comb_s15/jackstrawplot_pca.pdf', w = 100, h = 100, dpi = 300)

# CCA integration ---------------------------------------------------------
all_cells <- IntegrateLayers(all_cells, method = CCAIntegration, orig = "pca", new.reduction = "integrated.cca",
                             dims = 1:14, k.anchor = 20, verbose = T)
# cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "integrated.cca", dims = 1:14)
all_cells <- FindClusters(all_cells, resolution = 2, cluster.name = 'cca_clusters')
all_cells <- RunUMAP(all_cells, reduction = "integrated.cca", dims = 1:14, verbose = T)

#rejoin the layers in the sketched assay for DEG
# all_cells[["sketch"]] <- JoinLayers(all_cells[["sketch"]])
# c10_markers <- FindMarkers(object = all_cells, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
# head(c10_markers)

DimPlot(all_cells, group.by = "donor", reduction = "umap", pt.size = 0.1) %T>% figsave('/comb_s15/m8_sketch_cca_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap", pt.size = 0.1) %T>% figsave('/comb_s15/m8_sketch_cca_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "cca_clusters", reduction = "umap", pt.size = 0.1) %T>% figsave('/comb_s15/m8_sketch_cca_umap_cluster.pdf', w = 150, h = 100, dpi = 300)

# integrate the full datasets 
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
#all_cells[["sketch"]] <- split(all_cells[["sketch"]], f = all_cells$donor)

all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.cca")

all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.cca.full",
                         full.reduction = "integrated.cca.full", dims = 1:14, refdata = list(cca_clusters_full = "cca_clusters"))

#old plot m8_full_cca_umap_donor.pdf/ m8_full_cca_umap_cluster.pdf/ m8_full_cca_umap_transplant.pdf
#all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.cca.full",
#full.reduction = "integrated.cca.full", dims = 1:14, refdata = list(cluster_full = "seurat_clusters"))

all_cells <- RunUMAP(all_cells, reduction = "integrated.cca.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_")

DimPlot(all_cells, group.by = "donor", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_cca_umap_donor_new.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "cluster_full", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_cca_umap_cluster_new.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_cca_umap_transplant_new.pdf', w = 150, h = 100, dpi = 300)

saveRDS(all_cells@reductions$umap.full, file = "aml_all_cells_cca_mrd008_umapfull_reduction.Rds")

# RPCA integrate the datasets ---------------------------------------------
all_cells <- IntegrateLayers(all_cells, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
                             dims = 1:14, k.anchor = 20, reference = which(Layers(all_cells, search = "data") %in% c("data.MRD008")),
                             verbose = T)
# cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "integrated.rpca", dims = 1:14)
all_cells <- FindClusters(all_cells, resolution = 2, cluster.name = "rpca_clusters")

all_cells <- RunUMAP(all_cells, reduction = "integrated.rpca", dims = 1:14, verbose = T)

all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")

all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                         full.reduction = "integrated.rpca.full", dims = 1:14, refdata = list(rpca_clusters_full = "rpca_clusters"))

all_cells <- RunUMAP(all_cells, reduction = "integrated.rpca.full", dims = 1:14, reduction.name = "umap.full",
                     reduction.key = "UMAPfull_")

DimPlot(all_cells, group.by = "donor", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_rpca_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_rpca_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "cluster_full", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_rpca_umap_cluster.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "sample", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_rpca_umap_sample.pdf', w = 150, h = 100, dpi = 300)

saveRDS(all_cells@reductions$umap.full, file = "aml_all_cells_rpca_mrd008_umapfull_reduction.Rds")

# scVI integration --------------------------------------------------------
# integrate the datasets
# library(reticulate)
# 
# DefaultAssay(all_cells) <- "sketch"
all_cells[["sketch"]] <- split(all_cells[["sketch"]], f = all_cells$donor)

# set up scvi in conda env(aml_env)
#Sys.setenv(RETICULATE_PYTHON = "/home/zheng_song/.conda/envs/aml_env/bin/python")
# library(reticulate)
# reticulate::use_condaenv("aml_env")
# reticulate::use_python("/home/zheng_song/.conda/envs/aml_env/bin/python3")

#library(reticulate)
#/home/zheng_song/.conda/envs/aml_env/bin/python
DefaultAssay(all_cells) <- "sketch"
all_cells <- IntegrateLayers(object = all_cells, method = scVIIntegration, new.reduction = 'integrated.scvi',
                             dims = 1:14, verbose = T, conda_env = "/home/zheng_song/.conda/envs/aml_env")

# cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "integrated.scvi", dims = 1:14)
all_cells <- FindClusters(all_cells, resolution = 2, cluster.name = "scvi_clusters")

all_cells <- RunUMAP(all_cells, reduction = "integrated.scvi", dims = 1:14, verbose = T)

all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.scvi")

all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.scvi.full",
                         full.reduction = "integrated.scvi.full", dims = 1:14, refdata = list(scvi_clusters_full = "scvi_clusters"))

all_cells <- RunUMAP(all_cells, reduction = "integrated.scvi.full", dims = 1:14, reduction.name = "umap.full",
                     reduction.key = "UMAPfull_")

DimPlot(all_cells, group.by = "donor", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_scvi_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "cluster_full", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_scvi_umap_cluster.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_scvi_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "sample", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_scvi_umap_sample.pdf', w = 150, h = 100, dpi = 300)

saveRDS(all_cells@reductions$umap.full, file = "aml_all_cells_scvi_mrd008_umapfull_reduction.Rds")

# Harmony integration -----------------------------------------------------
all_cells <- IntegrateLayers(all_cells, method = HarmonyIntegration, orig = "pca", new.reduction = "harmony",
                             dims = 1:14, verbose = T)
# cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "harmony", dims = 1:14)
all_cells <- FindClusters(all_cells, resolution = 2, cluster.name = "harmony_clusters")
all_cells <- RunUMAP(all_cells, reduction = "harmony", dims = 1:14, verbose = T)

#rejoin the layers in the sketched assay for DEG
# all_cells[["sketch"]] <- JoinLayers(all_cells[["sketch"]])
# c10_markers <- FindMarkers(object = all_cells, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
# head(c10_markers)

DimPlot(all_cells, group.by = "donor", reduction = "umap", pt.size = 0.1) %T>% figsave('/comb_s15/m8_sketch_harmony_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap", pt.size = 0.1) %T>% figsave('/comb_s15/m8_sketch_harmony_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "harmony_clusters", reduction = "umap", pt.size = 0.1) %T>% figsave('/comb_s15/m8_sketch_harmony_umap_cluster.pdf', w = 150, h = 100, dpi = 300)

# integrate the full datasets 
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
#all_cells[["sketch"]] <- split(all_cells[["sketch"]], f = all_cells$donor)

all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")

all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                         full.reduction = "harmony.full", dims = 1:14, refdata = list(harmony_clusters_full = "harmony_clusters"))

#old plot m8_full_cca_umap_donor.pdf/ m8_full_cca_umap_cluster.pdf/ m8_full_cca_umap_transplant.pdf
#all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.cca.full",
#full.reduction = "integrated.cca.full", dims = 1:14, refdata = list(cluster_full = "seurat_clusters"))

all_cells <- RunUMAP(all_cells, reduction = "harmony.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_")

DimPlot(all_cells, group.by = "donor", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_harmony_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "cluster_full", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_harmony_umap_cluster.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_harmony_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "sample", reduction = "umap.full", alpha = 0.1) %T>% figsave('/comb_s15/m8_full_harmony_umap_sample.pdf', w = 150, h = 100, dpi = 300)

saveRDS(all_cells@reductions$umap.full, file = "aml_all_cells_harmony_mrd008_umapfull_reduction.Rds")

# saveRDS(
#   object = all_cells,
#   file = "aml_all_cells_rpca_mrd008.Rds",
#   destdir = paste0(wd_path, '/comb_s15/aml_all_cells_multi_inte_mrd008')
# )

saveRDS(
  object = all_cells,
  file = "aml_all_cells_cca_rpca_harmony_scvi_mrd008.Rds",
  destdir = paste0(wd_path, '/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008')
)

all_cells <- readRDS(paste0(wd_path, "/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno.Rds"))

# introduce TCR from TRUST4 -----------------------------------------------
mega_bcs <- read.csv(paste0(wd_path, '/mega_barcode.csv'), sep = ',')

mega_bcs$bc1_well <- str_pad(mega_bcs$bc1_well , 2, pad = "0")
mega_bcs$bc2_well <- str_pad(mega_bcs$bc2_well , 2, pad = "0")
mega_bcs$bc3_well <- str_pad(mega_bcs$bc3_well , 2, pad = "0")

mega_bcs

tcr_raw <- c(1:15) %>% map(~read.csv(paste0('~/big/zheng/parse_bio/aml/novaseq/analysis/trust4/S',.x,'_barcode_report.tsv'),sep = '\t') %>% 
                             distinct(X.barcode, .keep_all = T) %>% 
                             mutate(bc1 = str_sub(X.barcode, 17, 24), bc2 = str_sub(X.barcode, 9, 16), bc3 = str_sub(X.barcode, 1, 8)) %>%
                             left_join(select(mega_bcs, 1,2), by = 'bc1') %>%
                             left_join(select(mega_bcs, 3,4), by = 'bc2') %>% 
                             left_join(select(mega_bcs, 5,6), by = 'bc3')  %>% 
                             mutate(bc_wells = paste0(bc1_well, '_', bc2_well,'_', bc3_well,'__s', .x)) %>% 
                             filter(!grepl('NA', bc_wells)))

tcr_raw <- Reduce(rbind, tcr_raw)

tcr_raw %<>% distinct(bc_wells, .keep_all = T)

colnames(tcr_raw)

#IGH/TRB/TRD_information
chain1_names <- paste0(c('v', 'd', 'j', 'c', 'cdr3_nt', 'cdr3_aa',
                         'read_cnt','consensus_id','cdr3_germline_similarity','consensus_full_length'), '_chain1')

#IGK/IGL/TRA/TRG_information
chain2_names <- paste0(c('v', 'd', 'j', 'c', 'cdr3_nt', 'cdr3_aa', 
                         'read_cnt','consensus_id','cdr3_germline_similarity','consensus_full_length'), '_chain2')

tcr_bc <- select(tcr_raw, 2,3,4,13) 
colnames(tcr_bc)

tcr_bc[,chain1_names] <- tcr_bc$chain1 %>% str_split_fixed( ',', 10)

tcr_bc[,chain2_names] <- tcr_bc$chain2 %>% str_split_fixed( ',', 10)

tcr_bc

tcr_bc[,c('chain1', 'chain2')] <- NULL

tcr_bc %>% colnames()

mutate(tcr_bc,v_chain1 = replace(v_chain1, v_chain1 == '*', NA), v_chain2 = replace(v_chain2, v_chain2 == '*', NA))

freq_chain1 <- tcr_bc %>% filter(!is.na(cdr3_aa_chain1) & cdr3_aa_chain1 != '') %>%  group_by(cell_type) %>%  count(cdr3_aa_chain1, name = 'cdr3_chain1_freq' ) %>%  arrange(desc(cdr3_chain1_freq)) %>%  
  mutate(cdr3_chain1_perc = cdr3_chain1_freq/sum(cdr3_chain1_freq)*100) %T>% print()

freq_chain2 <- tcr_bc %>% filter(!is.na(cdr3_aa_chain2) & cdr3_aa_chain2 != '') %>%  group_by(cell_type) %>%  count(cdr3_aa_chain2, name = 'cdr3_chain2_freq' ) %>%  arrange(desc(cdr3_chain2_freq)) %>%  
  mutate(cdr3_chain2_perc = cdr3_chain2_freq/sum(cdr3_chain2_freq)*100) %T>% print()

class(tcr_bc)

tcr_bc$bc_wells  %>%  table() 
tcr_bc$cell_type %>%  table() 

# process gdTCR
gd_tcr <- subset(tcr_bc, cell_type == "gdT")
 
# gd_tcr %<>% left_join(freq_chain1, by = c ('cell_type', 'cdr3_aa_chain1'),  suffix = c('', '')) %>% 
#            left_join(freq_chain2, by = c ('cell_type', 'cdr3_aa_chain2'),  suffix = c('', '')) 

gd_tcr %>% colnames()

colnames(gd_tcr) <- gsub("chain1", replacement = "TRD", colnames(gd_tcr))

colnames(gd_tcr) <- gsub("chain2", replacement = "TRG", colnames(gd_tcr))

colnames(gd_tcr)

#remove pseudo trg gene 

trg_gene <- gd_tcr$v_TRG %>% table() %>% as_data_frame()

colnames(trg_gene) <- c('gene', 'freq')

trg_gene$gene_s <- gsub(pattern = "[*]\\d\\d$", replacement = "", trg_gene$gene)

# Funtional: TRGV 2 3 4 5 8 9 
# non-Funtioal: TRGV 1 10 11
# Pseudogene: TRGV A B 5P 6 7 

trg_gene$func <- case_when(trg_gene$gene_s %in% c("TRGV2", "TRGV3", "TRGV4", "TRGV5", "TRGV8", "TRGV9") ~ 'functional',
                           trg_gene$gene_s %in% c("TRGV1", "TRGV10", "TRGV11") ~ 'non-functional',
                           T ~ 'pseudogene')

colnames(trg_gene) <- c("v_TRG", "freq", "gene_s", "func_TRG")

gd_tcr %<>% left_join(trg_gene[, c('v_TRG', 'func_TRG')], by = 'v_TRG')

gd_tcr$v_TRD %>% table()

gd_tcr %<>% subset(v_TRD %in% c("TRDV1*01", "TRDV2*01", "TRDV2*02", "TRDV2*03", "TRDV3*01", NA))

freq_trd <- gd_tcr %>% filter(!is.na(cdr3_aa_TRD) & cdr3_aa_TRD!= '') %>% count(cdr3_aa_TRD, name = 'cdr3_freq_TRD' ) %>%  arrange(desc(cdr3_freq_TRD)) %>%  
  mutate(cdr3_perc_TRD = cdr3_freq_TRD/sum(cdr3_freq_TRD)*100) %T>% print()

freq_trg <- gd_tcr %>% filter(!is.na(cdr3_aa_TRG) & cdr3_aa_TRG!= '') %>% group_by(func_TRG) %>% count(cdr3_aa_TRG, name = 'cdr3_freq_TRG' ) %>%  arrange(desc(cdr3_freq_TRG)) %>%  
  mutate(cdr3_perc_TRG = cdr3_freq_TRG/sum(cdr3_freq_TRG)*100) %T>% print()

gd_tcr %<>% left_join(freq_trd, by = c ('cdr3_aa_TRD'),  suffix = c('', '')) %>% 
            left_join(freq_trg, by = c ('func_TRG', 'cdr3_aa_TRG'),  suffix = c('', '')) 

gd_tcr %<>% mutate(cdr3_perc_func_TRG = case_when(gd_tcr$func_TRG == 'functional' ~ cdr3_perc_TRG))

colnames(gd_tcr)

write.csv(gd_tcr, file = 'all_clean_gdTCR.csv', sep = ',')

# process abTCR
ab_tcr <- subset(tcr_bc, cell_type == "abT")

ab_tcr %>% colnames()

colnames(ab_tcr) <- gsub("chain1", replacement = "TRB", colnames(ab_tcr))

colnames(ab_tcr) <- gsub("chain2", replacement = "TRA", colnames(ab_tcr))

colnames(ab_tcr)

#remove TRDV abT cells
#TRDV1*01 TRDD3*01 TRAJ52*01 TRAC
ab_tcr %<>% filter(!v_TRA %in% c("TRDV1*01", "TRDV2*01", "TRDV3*01"))

#remove pseudo tra gene
tra_gene <- ab_tcr$v_TRA %>% table() %>% as_data_frame()

colnames(tra_gene) <- c('gene', 'freq')

tra_gene$gene_s <- gsub(pattern = "[*]\\d\\d$", replacement = "", tra_gene$gene)

tra_gene$gene_s 

# Funtional:
# non-Funtioal: TRAV8-7
# Pseudogene: TRAV8-5, TRAV11, TRAV15, TRAV28, TRAV31, TRAV32, TRAV33, and TRAV37 

tra_gene$func <- case_when(tra_gene$gene_s %in% c("TRAV8-5", "TRAV8-7", "TRAV11", "TRAV15", "TRAV28", "TRAV29/DV5", "TRAV31", "TRAV32", "TRAV33", "TRAV37") ~ 'pseudogene',
                           T ~ 'functional')

colnames(tra_gene) <- c("v_TRA", "freq", "gene_s", "func_TRA")

ab_tcr %<>% left_join(tra_gene[, c('v_TRA', 'func_TRA')], by = 'v_TRA')

###########
#remove pseudo tra gene
trb_gene <- ab_tcr$v_TRB %>% table() %>% as_data_frame()

colnames(trb_gene) <- c('gene', 'freq')

trb_gene$gene_s <- gsub(pattern = "[*]\\d\\d$", replacement = "", trb_gene$gene)

trb_gene$gene_s

# Funtional:
# non-Funtioal: TRAV8-7
# Pseudogene: TRAV8-5, TRAV11, TRAV15, TRAV28, TRAV31, TRAV32, TRAV33, and TRAV37 
hgnc_trb <- read.delim("hgnc_trbv.txt")

hgnc_trb$func_TRB <- str_extract(hgnc_trb$Name, "(?<=\\().+?(?=\\))")

colnames(hgnc_trb)[2] <- "gene_s"

trb_gene %<>% left_join(hgnc_trb[, c('gene_s', 'func_TRB')], by = 'gene_s')

trb_gene %<>% mutate(func_TRB = case_when(is.na(trb_gene$func_TRB) ~ "functional",
                    T ~ trb_gene$func_TRB))

colnames(trb_gene) <- c("v_TRB", "freq", "gene_s", "func_TRB")

ab_tcr %<>% left_join(trb_gene[, c('v_TRB', 'func_TRB')], by = 'v_TRB')


# freq / pct
freq_trb <- ab_tcr %>% filter(!is.na(cdr3_aa_TRB) & cdr3_aa_TRB!= '') %>% group_by(func_TRB) %>% count(cdr3_aa_TRB, name = 'cdr3_freq_TRB' ) %>%  arrange(desc(cdr3_freq_TRB)) %>%  
  mutate(cdr3_perc_TRB = cdr3_freq_TRB/sum(cdr3_freq_TRB)*100) %T>% print()

freq_tra <- ab_tcr %>% filter(!is.na(cdr3_aa_TRA) & cdr3_aa_TRA!= '') %>% group_by(func_TRA) %>% count(cdr3_aa_TRA, name = 'cdr3_freq_TRA' ) %>%  arrange(desc(cdr3_freq_TRA)) %>%  
  mutate(cdr3_perc_TRA = cdr3_freq_TRA/sum(cdr3_freq_TRA)*100) %T>% print()

ab_tcr %<>% left_join(freq_trb, by = c ('func_TRB', 'cdr3_aa_TRB'),  suffix = c('', '')) %>% 
            left_join(freq_tra, by = c ('func_TRA', 'cdr3_aa_TRA'),  suffix = c('', '')) 

ab_tcr %<>% mutate(cdr3_perc_func_TRB = case_when(ab_tcr$func_TRB == 'functional' ~ cdr3_perc_TRB))

ab_tcr %<>% mutate(cdr3_perc_func_TRA = case_when(ab_tcr$func_TRA == 'functional' ~ cdr3_perc_TRA))

colnames(ab_tcr)

write.csv(ab_tcr, "all_clean_abTCR.csv")

is.na(ab_tcr$cdr3_freq_TRA) %>% table()

# introduce TCR from rerun TRUST4 -----------------------------------------------
mega_bcs <- read.csv(paste0(wd_path, '/mega_barcode.csv'), sep = ',')

mega_bcs$bc1_well <- str_pad(mega_bcs$bc1_well , 2, pad = "0")
mega_bcs$bc2_well <- str_pad(mega_bcs$bc2_well , 2, pad = "0")
mega_bcs$bc3_well <- str_pad(mega_bcs$bc3_well , 2, pad = "0")

mega_bcs

all_cells$raw_bcs <- paste0(str_pad(all_cells$bc1_wind , 2, pad = "0"), "_", str_pad(all_cells$bc2_wind , 2, pad = "0"), "_", str_pad(all_cells$bc3_wind , 2, pad = "0"))

tcr_raw <- c(1:15) %>% map(~read.csv(paste0('~/big/zheng/parse_bio/aml/novaseq/analysis/trust4_rerun/S',.x,'_barcode_report.tsv'),sep = '\t') %>% 
                             distinct(X.barcode, .keep_all = T) %>% 
                             mutate(bc1 = str_sub(X.barcode, 69, 76 ), 
                                    bc2 = str_sub(X.barcode, 39, 46 ),
                                    bc3 = str_sub(X.barcode, 1, 8 )) %>%
                             left_join(select(mega_bcs, 1,2), by = 'bc1') %>%
                             left_join(select(mega_bcs, 3,4), by = 'bc2') %>% 
                             left_join(select(mega_bcs, 5,6), by = 'bc3')  #%>% 
                             #mutate(bc_wells = paste0(bc1_well, '_', bc2_well,'_', bc3_well,'__s', .x)) %>% 
                             #filter(!grepl('NA', bc_wells))
                             )

# sublib reporder after running 
tcr_raw[[1]]$bc_wells <- paste0(tcr_raw[[1]]$bc1_well, '_', tcr_raw[[1]]$bc2_well,'_', tcr_raw[[1]]$bc3_well,'__s1')
tcr_raw[[2]]$bc_wells <- paste0(tcr_raw[[2]]$bc1_well, '_', tcr_raw[[2]]$bc2_well,'_', tcr_raw[[2]]$bc3_well,'__s2')
tcr_raw[[3]]$bc_wells <- paste0(tcr_raw[[3]]$bc1_well, '_', tcr_raw[[3]]$bc2_well,'_', tcr_raw[[3]]$bc3_well,'__s9')
tcr_raw[[4]]$bc_wells <- paste0(tcr_raw[[4]]$bc1_well, '_', tcr_raw[[4]]$bc2_well,'_', tcr_raw[[4]]$bc3_well,'__s10')
tcr_raw[[5]]$bc_wells <- paste0(tcr_raw[[5]]$bc1_well, '_', tcr_raw[[5]]$bc2_well,'_', tcr_raw[[5]]$bc3_well,'__s11')

tcr_raw[[6]]$bc_wells <- paste0(tcr_raw[[6]]$bc1_well, '_', tcr_raw[[6]]$bc2_well,'_', tcr_raw[[6]]$bc3_well,'__s12')
tcr_raw[[7]]$bc_wells <- paste0(tcr_raw[[7]]$bc1_well, '_', tcr_raw[[7]]$bc2_well,'_', tcr_raw[[7]]$bc3_well,'__s13')
tcr_raw[[8]]$bc_wells <- paste0(tcr_raw[[8]]$bc1_well, '_', tcr_raw[[8]]$bc2_well,'_', tcr_raw[[8]]$bc3_well,'__s14')
tcr_raw[[9]]$bc_wells <- paste0(tcr_raw[[9]]$bc1_well, '_', tcr_raw[[9]]$bc2_well,'_', tcr_raw[[9]]$bc3_well,'__s15')
tcr_raw[[10]]$bc_wells <- paste0(tcr_raw[[10]]$bc1_well, '_', tcr_raw[[10]]$bc2_well,'_', tcr_raw[[10]]$bc3_well,'__s4')

tcr_raw[[11]]$bc_wells <- paste0(tcr_raw[[11]]$bc1_well, '_', tcr_raw[[11]]$bc2_well,'_', tcr_raw[[11]]$bc3_well,'__s5')
tcr_raw[[12]]$bc_wells <- paste0(tcr_raw[[12]]$bc1_well, '_', tcr_raw[[12]]$bc2_well,'_', tcr_raw[[12]]$bc3_well,'__s6')
tcr_raw[[13]]$bc_wells <- paste0(tcr_raw[[13]]$bc1_well, '_', tcr_raw[[13]]$bc2_well,'_', tcr_raw[[13]]$bc3_well,'__s7')
tcr_raw[[14]]$bc_wells <- paste0(tcr_raw[[14]]$bc1_well, '_', tcr_raw[[14]]$bc2_well,'_', tcr_raw[[14]]$bc3_well,'__s8')
tcr_raw[[15]]$bc_wells <- paste0(tcr_raw[[15]]$bc1_well, '_', tcr_raw[[15]]$bc2_well,'_', tcr_raw[[15]]$bc3_well,'__s3')

tcr_raw <- reduce(tcr_raw, rbind) # 43068

tcr_raw %<>% filter(!grepl('NA', bc_wells)) #32490

tcr_raw %<>% distinct(bc_wells, .keep_all = T) #20974

colnames(tcr_raw)

#IGH/TRB/TRD_information
chain1_names <- paste0(c('v', 'd', 'j', 'c', 'cdr3_nt', 'cdr3_aa',
                         'read_cnt','consensus_id','cdr3_germline_similarity','consensus_full_length'), '_chain1')

#IGK/IGL/TRA/TRG_information
chain2_names <- paste0(c('v', 'd', 'j', 'c', 'cdr3_nt', 'cdr3_aa', 
                         'read_cnt','consensus_id','cdr3_germline_similarity','consensus_full_length'), '_chain2')

tcr_bc <- select(tcr_raw, 2,3,4,13) 
colnames(tcr_bc)

tcr_bc[,chain1_names] <- tcr_bc$chain1 %>% str_split_fixed( ',', 10)

tcr_bc[,chain2_names] <- tcr_bc$chain2 %>% str_split_fixed( ',', 10)

tcr_bc

tcr_bc[,c('chain1', 'chain2')] <- NULL

tcr_bc %>% colnames()

tcr_bc <- mutate(tcr_bc,v_chain1 = replace(v_chain1, v_chain1 == '*', NA), v_chain2 = replace(v_chain2, v_chain2 == '*', NA))

# freq_chain1 <- tcr_bc %>% filter(!is.na(cdr3_aa_chain1) & cdr3_aa_chain1 != '') %>%  group_by(cell_type) %>%  count(cdr3_aa_chain1, name = 'cdr3_chain1_freq' ) %>%  arrange(desc(cdr3_chain1_freq)) %>%  
#   mutate(cdr3_chain1_perc = cdr3_chain1_freq/sum(cdr3_chain1_freq)*100) %T>% print()
# 
# freq_chain2 <- tcr_bc %>% filter(!is.na(cdr3_aa_chain2) & cdr3_aa_chain2 != '') %>%  group_by(cell_type) %>%  count(cdr3_aa_chain2, name = 'cdr3_chain2_freq' ) %>%  arrange(desc(cdr3_chain2_freq)) %>%  
#   mutate(cdr3_chain2_perc = cdr3_chain2_freq/sum(cdr3_chain2_freq)*100) %T>% print()

class(tcr_bc)

tcr_bc$bc_wells  %>%  table() 

tcr_bc$cell_type %>%  table() 

# process gdTCR
gd_tcr <- subset(tcr_bc, cell_type == "gdT")

# gd_tcr %<>% left_join(freq_chain1, by = c ('cell_type', 'cdr3_aa_chain1'),  suffix = c('', '')) %>% 
#            left_join(freq_chain2, by = c ('cell_type', 'cdr3_aa_chain2'),  suffix = c('', '')) 

gd_tcr %>% colnames()

colnames(gd_tcr) <- gsub("chain1", replacement = "TRD", colnames(gd_tcr))

colnames(gd_tcr) <- gsub("chain2", replacement = "TRG", colnames(gd_tcr))

colnames(gd_tcr)

#remove pseudo trg gene 

trg_gene <- gd_tcr$v_TRG %>% table() %>% as_data_frame()

colnames(trg_gene) <- c('gene', 'freq')

trg_gene$gene_s <- gsub(pattern = "[*]\\d\\d$", replacement = "", trg_gene$gene)

# Funtional: TRGV 2 3 4 5 8 9 
# non-Funtioal: TRGV 1 10 11
# Pseudogene: TRGV A B 5P 6 7 

trg_gene$func <- case_when(trg_gene$gene_s %in% c("TRGV2", "TRGV3", "TRGV4", "TRGV5", "TRGV8", "TRGV9") ~ 'functional',
                           trg_gene$gene_s %in% c("TRGV1", "TRGV10", "TRGV11") ~ 'non-functional',
                           T ~ 'pseudogene')

colnames(trg_gene) <- c("v_TRG", "freq", "gene_s", "func_TRG")

gd_tcr %<>% left_join(trg_gene[, c('v_TRG', 'func_TRG')], by = 'v_TRG')

gd_tcr$v_TRD %>% table()

gd_tcr %<>% subset(v_TRD %in% c("TRDV1*01", "TRDV2*01", "TRDV2*02", "TRDV2*03", "TRDV3*01", "TRAV14/DV4*02", "TRAV29/DV5*01", "TRAV38-2/DV8*01", NA))

freq_trd <- gd_tcr %>% filter(!is.na(cdr3_aa_TRD) & cdr3_aa_TRD!= '') %>% count(cdr3_aa_TRD, name = 'cdr3_freq_TRD' ) %>%  arrange(desc(cdr3_freq_TRD)) %>%  
  mutate(cdr3_perc_TRD = cdr3_freq_TRD/sum(cdr3_freq_TRD)*100) %T>% print()

freq_trg <- gd_tcr %>% filter(!is.na(cdr3_aa_TRG) & cdr3_aa_TRG!= '') %>% group_by(func_TRG) %>% count(cdr3_aa_TRG, name = 'cdr3_freq_TRG' ) %>%  arrange(desc(cdr3_freq_TRG)) %>%  
  mutate(cdr3_perc_TRG = cdr3_freq_TRG/sum(cdr3_freq_TRG)*100) %T>% print()

gd_tcr %<>% left_join(freq_trd, by = c ('cdr3_aa_TRD'),  suffix = c('', '')) %>% 
  left_join(freq_trg, by = c ('func_TRG', 'cdr3_aa_TRG'),  suffix = c('', '')) 

gd_tcr %<>% mutate(cdr3_perc_func_TRG = case_when(gd_tcr$func_TRG == 'functional' ~ cdr3_perc_TRG))

colnames(gd_tcr)

# gd_tcr 4491
nrow(gd_tcr)

write.csv(gd_tcr, file = 'all_clean_gdTCR_rerun.csv', sep = ',')

# process abTCR
ab_tcr <- subset(tcr_bc, cell_type == "abT")

ab_tcr %>% colnames()

colnames(ab_tcr) <- gsub("chain1", replacement = "TRB", colnames(ab_tcr))

colnames(ab_tcr) <- gsub("chain2", replacement = "TRA", colnames(ab_tcr))

colnames(ab_tcr)

#remove TRDV abT cells
#TRDV1*01 TRDD3*01 TRAJ52*01 TRAC
ab_tcr %<>% filter(!v_TRA %in% c("TRDV1*01", "TRDV2*01", "TRDV3*01"))

#remove pseudo tra gene
tra_gene <- ab_tcr$v_TRA %>% table() %>% as_data_frame()

colnames(tra_gene) <- c('gene', 'freq')

tra_gene$gene_s <- gsub(pattern = "[*]\\d\\d$", replacement = "", tra_gene$gene)

tra_gene$gene_s 

# Funtional:
# non-Funtioal: TRAV8-7
# Pseudogene: TRAV8-5, TRAV11, TRAV15, TRAV28, TRAV31, TRAV32, TRAV33, and TRAV37 

tra_gene$func <- case_when(tra_gene$gene_s %in% c("TRAV8-5", "TRAV8-7", "TRAV11", "TRAV15", "TRAV28", "TRAV29/DV5", "TRAV31", "TRAV32", "TRAV33", "TRAV37") ~ 'pseudogene',
                           T ~ 'functional')

colnames(tra_gene) <- c("v_TRA", "freq", "gene_s", "func_TRA")

ab_tcr %<>% left_join(tra_gene[, c('v_TRA', 'func_TRA')], by = 'v_TRA')

###########
#remove pseudo tra gene
trb_gene <- ab_tcr$v_TRB %>% table() %>% as_data_frame()

colnames(trb_gene) <- c('gene', 'freq')

trb_gene$gene_s <- gsub(pattern = "[*]\\d\\d$", replacement = "", trb_gene$gene)

trb_gene$gene_s

# Funtional:
# non-Funtioal: TRAV8-7
# Pseudogene: TRAV8-5, TRAV11, TRAV15, TRAV28, TRAV31, TRAV32, TRAV33, and TRAV37 
hgnc_trb <- read.delim("hgnc_trbv.txt")

hgnc_trb$func_TRB <- str_extract(hgnc_trb$Name, "(?<=\\().+?(?=\\))")

colnames(hgnc_trb)[2] <- "gene_s"

trb_gene %<>% left_join(hgnc_trb[, c('gene_s', 'func_TRB')], by = 'gene_s')

trb_gene %<>% mutate(func_TRB = case_when(is.na(trb_gene$func_TRB) ~ "functional",
                                          T ~ trb_gene$func_TRB))

colnames(trb_gene) <- c("v_TRB", "freq", "gene_s", "func_TRB")

ab_tcr %<>% left_join(trb_gene[, c('v_TRB', 'func_TRB')], by = 'v_TRB')

# freq / pct
freq_trb <- ab_tcr %>% filter(!is.na(cdr3_aa_TRB) & cdr3_aa_TRB!= '') %>% group_by(func_TRB) %>% count(cdr3_aa_TRB, name = 'cdr3_freq_TRB' ) %>%  arrange(desc(cdr3_freq_TRB)) %>%  
  mutate(cdr3_perc_TRB = cdr3_freq_TRB/sum(cdr3_freq_TRB)*100) %T>% print()

freq_tra <- ab_tcr %>% filter(!is.na(cdr3_aa_TRA) & cdr3_aa_TRA!= '') %>% group_by(func_TRA) %>% count(cdr3_aa_TRA, name = 'cdr3_freq_TRA' ) %>%  arrange(desc(cdr3_freq_TRA)) %>%  
  mutate(cdr3_perc_TRA = cdr3_freq_TRA/sum(cdr3_freq_TRA)*100) %T>% print()

ab_tcr %<>% left_join(freq_trb, by = c ('func_TRB', 'cdr3_aa_TRB'),  suffix = c('', '')) %>% 
  left_join(freq_tra, by = c ('func_TRA', 'cdr3_aa_TRA'),  suffix = c('', '')) 

ab_tcr %<>% mutate(cdr3_perc_func_TRB = case_when(ab_tcr$func_TRB == 'functional' ~ cdr3_perc_TRB))

ab_tcr %<>% mutate(cdr3_perc_func_TRA = case_when(ab_tcr$func_TRA == 'functional' ~ cdr3_perc_TRA))

colnames(ab_tcr)

# ab_tcr 6474
nrow(ab_tcr)

write.csv(ab_tcr, "all_clean_abTCR_rerun.csv")

all_cells$bc_wells <- all_cells$bc

all_cells@meta.data  %<>%  left_join(ab_tcr, by = 'bc_wells', suffix = c('', ''))%>% 
  `rownames<-`(all_cells$bc)

# all_tcr 10965
all_cells$cell_type %>% table()

all_cells$bc_wells <- all_cells$bc

all_tcr <- full_join(ab_tcr, gd_tcr, suffix = c("", ""))

all_cells@meta.data  %<>%  left_join(all_tcr, by = 'bc_wells', suffix = c('', ''))%>% 
  `rownames<-`(all_cells$bc)

#  abT  gdT 
# 5787 4012 
DimPlot(all_cells, "predicted.celltype.l2", )

Feature_rast5(all_cells, "predicted.celltype.l2", do.label = T, colorset = "gg")

Feature_rast(all_cells, "cdr3_perc_func_TRB", color_grd = "grd", sz = 0.1, other = (alpha = 0.1))

saveRDS(
  object = all_cells,
  file = "aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno_tcr.Rds",
  destdir = paste0(wd_path, '/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno_tcr')
)

all_cells <- readRDS(paste0(wd_path, "/comb_s15/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno_tcr/aml_all_cells_cca_rpca_harmony_scvi_mrd008_anno_tcr.Rds"))
