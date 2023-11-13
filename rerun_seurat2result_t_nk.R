#characterize t_nk cluster 
#Seurat V5
#zheng

# packages and functions ----------------------------------------------------------------
library(Seurat) #v4.9.9.9058
library(SeuratObject)
library(BPCells)
library(dplyr)
library(purrr)
library(DESeq2)
library(Azimuth)

library(cowplot) # ggsave2
library(ggrastr) # geom_point_rast
library(ggplot2)
library(ggrepel)
library(patchwork)

#install.packages("viridis")
library(viridis)

library(stringr) # string manipulation
library(magrittr) # %>% 

library(future)
# functions
source('/home/big/zheng/rscripts/funcs.r')
source('/home/big/zheng/rscripts/funcs_umapfull.r')
source('/home/big/zheng/rscripts/themes.R')



# set this option when analyzing large datasets
options(future.globals.maxSize = 256*1024^3)
options(Seurat.object.assay.version = "v5")

future::plan(strategy = "multicore", workers = 32)

setwd("/home/big/zheng_song/aml")
wd_path <- getwd()

all_cells <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_all_cells_cca_harmony_mrd008_tcr_gsea/aml_all_cells_cca_harmony_mrd008_tcr_gsea.Rds"))

# Harmony integration -----------------------------------------------------
# we rerun harmony to increase the resolution
DefaultAssay(all_cells) <- 'sketch'

all_cells <- IntegrateLayers(all_cells, method = HarmonyIntegration, orig = "pca", new.reduction = "harmony",
                             dims = 1:14, verbose = T)
# cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "harmony", dims = 1:14)
all_cells <- FindClusters(all_cells, resolution = 2.5, cluster.name = "harmony_clusters_reso2.5")
all_cells <- RunUMAP(all_cells, reduction = "harmony", dims = 1:14, verbose = T)

#rejoin the layers in the sketched assay for DEG
# all_cells[["sketch"]] <- JoinLayers(all_cells[["sketch"]])
# c10_markers <- FindMarkers(object = all_cells, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
# head(c10_markers)

DimPlot(all_cells, group.by = "donor", reduction = "umap", pt.size = 0.1, raster = T) %T>% figsave('/comb_s15_rerun/m8_sketch_harmony_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap", pt.size = 0.1, raster = T) %T>% figsave('/comb_s15_rerun/m8_sketch_harmony_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "harmony_clusters_reso2.5", reduction = "umap", pt.size = 0.1, raster = T) %T>% figsave('/comb_s15_rerun/m8_sketch_harmony_umap_cluster_reso2.5.pdf', w = 150, h = 100, dpi = 300)

# integrate the full datasets 
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
all_cells[["sketch"]] <- split(all_cells[["sketch"]], f = all_cells$donor)
all_cells[["RNA"]] <- split(all_cells[["RNA"]], f = all_cells$donor)

all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", assay = "RNA", reduction = "harmony", ratio = 1)

all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                         full.reduction = "harmony.full", dims = 1:10, 
                         refdata = list(harmony_clusters_reso2.5_full = "harmony_clusters_reso2.5"),
                         normalization.method = "SCT")

#old plot m8_full_cca_umap_donor.pdf/ m8_full_cca_umap_cluster.pdf/ m8_full_cca_umap_transplant.pdf
#all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.cca.full",
#full.reduction = "integrated.cca.full", dims = 1:14, refdata = list(cluster_full = "seurat_clusters"))

all_cells <- RunUMAP(all_cells, reduction = "harmony.full", dims = 1:10, reduction.name = "umap.full", reduction.key = "UMAPfull_")

all_cells[["RNA"]] <- JoinLayers(all_cells[["RNA"]])
DefaultAssay(all_cells) <- 'RNA'

DimPlot(all_cells, group.by = "donor", reduction = "umap.full", alpha = 0.1, raster = T) %T>% 
  figsave('/comb_s15_rerun/m8_full_harmony_umap_donor.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "transplant", reduction = "umap.full", alpha = 0.1, raster = T) %T>% 
  figsave('/comb_s15_rerun/m8_full_harmony_umap_transplant.pdf', w = 150, h = 100, dpi = 300)
DimPlot(all_cells, group.by = "sample", reduction = "umap.full", alpha = 0.1, raster = T) %T>% 
  figsave('/comb_s15_rerun/m8_full_harmony_umap_sample.pdf', w = 150, h = 100, dpi = 300)

DimPlot(all_cells, group.by = "harmony_clusters_reso2.5_full", reduction = "umap.full", alpha = 0.1, raster = T, label = T, ) %T>% 
  figsave('/comb_s15_rerun/reso2.5/m8_full_harmony_umap_cluster.pdf', w = 150, h = 100, dpi = 300)

# all cells ---------------------------------------------------------------


all_cells@meta.data %<>% mutate(anno_2nd = case_when(harmony_clusters_reso2.5_full %in% c("39", "37", "42", "15", "27", "6", "26") ~ "B&Plasma cells",
                                                     harmony_clusters_reso2.5_full %in% c("22", "16", "21", "18", "33", "7", "8", "30",
                                                                                          "31", "35", "1", "11", "13", "14", "10", "3", "19",
                                                                                          "24", "5") ~ "Mono&DC cells",
                                                     harmony_clusters_reso2.5_full %in% c("32", "43", "41", "36", "44", "34") ~ "Mixed cells",
                                                     harmony_clusters_reso2.5_full %in% c("25", "0", "17", "23", "2", "4", "28", "20", "12",
                                                                                          "9", "40", "38", "29") ~ "T&NK cells"))

#umap.colors <- c("#5A5F6D80", "#6B5695", "#A29A7680",  "#a98a7b", "#95a674")
umap.colors <- c("#6B5695", "#5A5F6D80",  "#a98a7b", "#95a674")

(Feature_rast5(all_cells, "anno_2nd", do.label = T, colorset = "um", sz = 0.02, noaxis = F, titlesize = 0) + 
    xlab("UMAP_1") +
    ylab("UMAP_2")) %>% figsave("comb_s15_rerun/reso2.5/umap_anno_2nd.pdf", w = 150, h = 100) 

# dotplot of marker genes on four main populations
all_cells$anno_2nd <- factor(all_cells$anno_2nd, levels = c("T&NK cells", "B&Plasma cells", "Mono&DC cells", "Proliferating cells"))

Idents(all_cells) <- all_cells$anno_2nd

(Seurat::DotPlot(all_cells, c("CD3E", "CD3D", "CD3G", "NCAM1", "FCGR3A", "MS4A1", "CD19", 
                      "CD14", "VCAN", "LYZ", "CD163", "MCM2", "MKI67", "CD34"), 
                assay = "RNA",cols = c("#9fc5e8","#990000"), scale.max = 14)+
                guides(size = F,color = guide_colorbar(title = "Expression", title.position = "right"))+
                #scale_y_discrete()+
                #scale_x_discrete()+
                dotplot_theme + 
    theme(axis.text.x.bottom = element_text(angle = 30, hjust = 0.9, vjust = 0.9)))%T>% 
  figsave("/comb_s15_rerun/reso2.5/dotplot_marker_genes_main.pdf", 130, 50)

saveRDS(
  object = all_cells,
  file = "aml_all_cells_cca_harmony_mrd008_tcr_gsea_reso2.5.Rds",
  destdir = paste0(wd_path, '/comb_s15_rerun/aml_all_cells_cca_harmony_mrd008_tcr_gsea_reso2.5')
)

all_cells <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_all_cells_cca_harmony_mrd008_tcr_gsea_reso2.5/aml_all_cells_cca_harmony_mrd008_tcr_gsea_reso2.5.Rds"))

# T&NK compartment --------------------------------------------------------
t_nk <- subset(all_cells, anno_2nd == "T&NK cells")
Feature_rast5(t_nk, colorset = "gg")

t_nk <- FindNeighbors(t_nk, reduction = "harmony.full", dims = 1:24)
t_nk <- FindClusters(t_nk, resolution = 1.5, cluster.name = "harmony_clusters_reso1.5_full")
t_nk <- RunUMAP(t_nk, reduction = "harmony.full", dims = 1:24, reduction.name = "umap.full", reduction.key = "UMAPfull_")

t_nk <- Azimuth::RunAzimuth(t_nk, reference = "pbmcref", assay = "RNA", verbose = T)

Feature_rast5(t_nk, colorset = "gg", sz = 0.1, g = "predicted.celltype.l3") %T>% 
  figsave("/comb_s15_rerun/t_nk/azimuth_l3.pdf", w = 150, h = 100)

Feature_rast(t_nk, colorset = "gg", sz = 0.1, g = "predicted.celltype.l2") 

Feature_rast5(t_nk, colorset = "gg", sz = 0.3, g = "predicted.celltype.l2") %T>% 
  figsave("/comb_s15_rerun/t_nk/azimuth_l2.pdf", w = 150, h = 100)

Feature_rast5(t_nk, colorset = "gg", sz = 0.1, g = "harmony_clusters_reso1.5_full")%T>% 
  figsave("/comb_s15_rerun/t_nk/umap_cluster_2.pdf", w = 150, h = 100)

Feature_rast5(t_nk, c("CD3E", "CD3D", "CD3G", "CD4", "FOXP3",
                      "CD8A", "CD8B", "KLRB1", "SLC4A10", "TRDC", 
                      "TRGC1", "NCAM1", "FCGR3A", "B3GAT1"), 
              color_grd = "grd", sz = 0.2, ncol = 5) %T>% 
  figsave("/comb_s15_rerun/t_nk/umap_marker_genes_2.pdf", w = 500, h = 300)

VlnPlot(t_nk, features = c("CD4", "FOXP3", "CD8A", "CD8B", "TRDC", "TRGC1"), pt.size = 0.05, raster = T, ncol = 3) %T>% 
  figsave("/comb_s15_rerun/t_nk/vln_marker_genes.pdf", w = 1000, h = 200)

VlnPlot(t_nk, features = c("CD4", "FOXP3", "IL2RA", 
                           "CTLA4", "MTNFRSF4", "BATF"), pt.size = 0.05, raster = T, ncol = 3) %T>% 
  figsave("/comb_s15_rerun/t_nk/vln_marker_genes_treg.pdf", w = 1000, h = 200)

VlnPlot(t_nk, features = c("CD4", "FOXP3", "IL2RA", 
                           "CTLA4", "MTNFRSF4", "BATF"), pt.size = 0.05, raster = T, ncol = 3) %T>% 
  figsave("/comb_s15_rerun/t_nk/vln_marker_genes_treg.pdf", w = 1000, h = 200)

# gene module from Tan_SI_2019 
gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_c_list_ent

gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)

gene_cluster$GM_D

for (i in gc_name) {
  t_nk %<>%  AddModuleScore(features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
}

Feature_rast5(t_nk, c("GM_A1", "GM_B1", "GM_C1", "GM_D1", "GM_E1", "GM_F1", "GM_G1", "GM_H1"), color_grd = "grd", sz = 0.2, ncol = 4) %T>% 
  figsave("/comb_s15_rerun/t_nk/umap_cluster_gm.pdf", w = 400, h = 200)

VlnPlot(t_nk, features = c("GM_A1", "GM_B1", "GM_C1", "GM_D1", "GM_E1", "GM_F1", "GM_G1", "GM_H1"), pt.size = 0, raster = T, ncol = 4) %T>% 
  figsave("/comb_s15_rerun/t_nk/vln_gm.pdf", w = 1000, h = 200)

VlnPlot(t_nk, features = c("TRDC", "TRDV1", "TRDV2", "TRGV9"), pt.size = 0.1, raster = T, ncol = 3) %T>% 
  figsave("/comb_s15_rerun/t_nk/vln_gdt_2.pdf", w = 1000, h = 200)

c("4") ~ "CD8_Tcm/naive" 
c("5") ~ "CD4_Tcm/naive"
#c("5") ~ "CD4_Tcm/naive"
c("14") ~ "Treg"
c("15") ~ "CD4_Th17"
c("9", "7", "27") ~ "CD8_Tem"
c("11") ~ "CTL"
c("6") ~ "MAIT"
c("8") ~ "gdT"
c("2") ~ "CD8_Temra"
c("10") ~ "CD8_CD56_TNF" 
c("16") ~ "CD8_CD56"

c("5", "6") ~ "NK_CD56bright"
c("1", "0", "21", "3") ~ "NK_CD56dim"

c("29") ~ "ILC"

c("24") ~ "Unidentified"

c("19", "25", "26", "18", "20", "28", "22", "17", "23") ~ "Mixed"

marker_1 <- FindMarkers(t_nk, ident.1 = "9", ident.2 = "4", group.by = "harmony_clusters_reso1.5_full", assay = "RNA")

marker_1 %>% arrange(-avg_log2FC)

marker_2 <- FindMarkers(t_nk, ident.1 = "7", ident.2 = "9", group.by = "harmony_clusters_reso1.5_full", assay = "RNA")

marker_2 %>% arrange(-avg_log2FC)

marker_3 <- FindMarkers(t_nk, ident.1 = "11", ident.2 = "2", group.by = "harmony_clusters_reso1.5_full", assay = "RNA")

marker_3 %>% arrange(avg_log2FC)

marker_4 <- FindMarkers(t_nk, ident.1 = "24", ident.2 = "8", group.by = "harmony_clusters_reso1.5_full", assay = "RNA")

marker_4 %>% arrange(-avg_log2FC)

marker_5 <- FindMarkers(t_nk, ident.1 = "10", ident.2 = "16", group.by = "harmony_clusters_reso1.5_full", assay = "RNA")

marker_5 %>% arrange(avg_log2FC)

marker_6 <- FindMarkers(t_nk, ident.1 = "16", ident.2 = "11", group.by = "harmony_clusters_reso1.5_full", assay = "RNA")

marker_6 %>% arrange(avg_log2FC)

ctl_gm <- c("FAS","CD3G","ICAM1","ITGB2","FASLG","CD3D","CD3E","CD247","ITGB2","HLA-A","ITGAL","B2M","GZMB","PRF1","CD247")
t_nk <- AddModuleScore(t_nk, 
                       features = list(c(ctl_gm)),
                       name = "ctl_gm", assay ="RNA")

abTR <- grep("^TRAV|^TRBV|^TRAC|^TRBC" , rownames(t_nk@assays$RNA), value = T)
TRD <- grep("^TRDV|^TRDC" , rownames(t_nk@assays$RNA), value = T)
t_nk <- AddModuleScore(t_nk, 
                            features = list(c(abTR)),
                            name = "alpha_beta_score_tr", assay = "RNA")

t_nk <- AddModuleScore(t_nk, 
                            features = list(c(TRD)),
                            name = "gamma_delta_score_trd", assay ="RNA")

Feature_rast5(t_nk, c("alpha_beta_score_tr1", "gamma_delta_score_trd1", "ctl_gm1"), assay = "RNA", sz = 0.4, ncol = 1, color_grd = "grd")

Feature_rast5(t_nk, c("CCR7", "SELL", "LEF1",
                      "GZMB", "GNLY", "GZMK",
                      "PDCD1", "TIGIT", "HAVCR2",
                      "RORC", "CTLA4", "CCR6"), assay = "RNA", sz = 0.4, ncol = 3, color_grd = "grd")


VlnPlot(t_nk, c("alpha_beta_score_tr1", "gamma_delta_score_trd1", "ctl_gm1"), assay = "RNA", pt.size = 0, ncol = 1)

# focus on unknown/mix cluster

Feature_rast5(t_nk, colorset = "gg", sz = 0.1, g = "percent_mt", color_grd = "grd")

VlnPlot(t_nk, "nCount_RNA")

lin <- c("CD2", "CD3E", "CD14", "FCGR3A", "CD19", "MS4A1", "NCAM1", "GYPA")

t_nk <- AddModuleScore(t_nk, features = list(lin), assay = "RNA", name = "lin")

VlnPlot(t_nk, c("lin1", "IL7R"), assay = "RNA", pt.size = 0, ncol = 1)

Feature_rast5(t_nk, c("cell_type", "TRDV1"), assay = "RNA", sz = 2, ncol = 3, color_grd = "grd", colorset = "gg")

# focus on NK compartment
Feature_rast5(t_nk, c("GZMK", "XCL1", "IL7R", "TCF7", "GPR183", 
                      "GZMB", "PRF1", "CX3CR1", 
                      "CD7", "FCER1G", "KLRB1",
                      "KLRC2", "CD3E", "PATL2", "ZBTB38"), assay = "RNA", sz = 0.4, ncol = 3, color_grd = "grd")

Feature_rast5(t_nk, c("TRAC", "TRBC1", "TRBC2"), assay = "RNA", sz = 0.4, ncol = 3, color_grd = "grd")

Feature_rast5(t_nk, c("JAKMIP1", "CD3D", "IL32", "ZBTB16"), assay = "RNA", sz = 0.4, ncol = 3, color_grd = "grd")

Feature_rast5(t_nk, c("CD69", "IFNG", "TNF", "EOMES", "TBX21", "RUNX2"), assay = "RNA", sz = 0.4, ncol = 3, color_grd = "grd")

# anno t_nk compartment ---------------------------------------------------
 
t_nk@meta.data %<>% mutate(anno_3rd= case_when(harmony_clusters_reso1.5_full %in% c("4") ~ "CD8_Tcm/naive",
                                               harmony_clusters_reso1.5_full %in% c("13") ~ "CD4_Tcm/naive",
                                               harmony_clusters_reso1.5_full %in% c("14") ~ "Treg",
                                               harmony_clusters_reso1.5_full %in% c("15") ~ "CD4_Th17",
                                               harmony_clusters_reso1.5_full %in% c("9", "7", "27") ~ "CD8_Tem",
                                               harmony_clusters_reso1.5_full %in% c("11") ~ "CTL",
                                               harmony_clusters_reso1.5_full %in% c("12") ~ "MAIT",
                                               harmony_clusters_reso1.5_full %in% c("8") ~ "gdT",
                                               harmony_clusters_reso1.5_full %in% c("2") ~ "CD8_Temra",
                                               harmony_clusters_reso1.5_full %in% c("10") ~ "CD8_CD56_TNF",
                                               harmony_clusters_reso1.5_full %in% c("16") ~ "CD8_CD56", 
                                               harmony_clusters_reso1.5_full %in% c("5", "6") ~ "NK_CD56bright",
                                               harmony_clusters_reso1.5_full %in% c("1", "0", "21", "3") ~ "NK_CD56dim",
                                               harmony_clusters_reso1.5_full %in% c("29") ~ "ILC",
                                               harmony_clusters_reso1.5_full %in% c("24") ~ "Unidentified",
                                               harmony_clusters_reso1.5_full %in% c("19", "25", "26", "18", "20", "28", "22", "17", "23") ~ "Mixed")) 

Idents(t_nk) <- t_nk$anno_3rd

saveRDS(
  object = t_nk,
  file = "aml_t_nk_anno.Rds",
  destdir = paste0(wd_path, '/comb_s15_rerun/aml_t_nk_anno')
)

t_nk <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_t_nk_anno/aml_t_nk_anno.Rds"))

Feature_rast5(t_nk, colorset = "um", g = "harmony_clusters_reso2.5_full")

# we increase reso to separate Treg cells (foxp3+ foxp3-)
t_nk <- FindNeighbors(t_nk, reduction = "harmony.full", dims = 1:24)
t_nk <- FindClusters(t_nk, resolution = 2, cluster.name = "harmony_clusters_reso2_full")
t_nk <- RunUMAP(t_nk, reduction = "harmony.full", dims = 1:24, reduction.name = "umap.full", reduction.key = "UMAPfull_")

Feature_rast5(t_nk, colorset = "um", g = "harmony_clusters_reso2_full")

# we filed simple by simply increasing reso to separate Treg cells (foxp3+ foxp3-)
# we try to include more PCs(24-->30)
t_nk <- FindNeighbors(t_nk, reduction = "harmony.full", dims = 1:30)
t_nk <- FindClusters(t_nk, resolution = 3, cluster.name = "harmony_clusters_reso3_full")
t_nk <- RunUMAP(t_nk, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAPfull_")

Feature_rast5(t_nk, colorset = "gg", g = "harmony_clusters_reso3_full", noaxis = F, do.label = T, sz = 0.03) %T>% 
  figsave("/comb_s15_rerun/t_nk/umap_clusters_3.pdf", w = 120, h = 75)

Feature_rast5(t_nk, c("CD3E", "CD3D", "CD3G", "CD4", "FOXP3",
                      "CD8A", "CD8B", "KLRB1", "SLC4A10", "TRDC", 
                      "TRGC1", "NCAM1", "FCGR3A", "B3GAT1"), 
              color_grd = "grd", sz = 0.2, ncol = 5) %T>% 
  figsave("/comb_s15_rerun/t_nk/umap_marker_genes_reso3.pdf", w = 500, h = 300)

# gene module from Tan_SI_2019 
gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)

gene_cluster$GM_D

for (i in gc_name) {
  t_nk %<>%  AddModuleScore(features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
}

t_nk@meta.data %<>% mutate(anno_3rd= case_when(harmony_clusters_reso3_full %in% c("19", "21") ~ "CD8_Tcm/naive",
                                               harmony_clusters_reso3_full %in% c("25") ~ "CD4_Tcm/naive",
                                               harmony_clusters_reso3_full %in% c("35") ~ "Treg",
                                               harmony_clusters_reso3_full %in% c("20") ~ "Th17",
                                               harmony_clusters_reso3_full %in% c("9", "18", "48", "16") ~ "Th1/Th2",
                                               harmony_clusters_reso3_full %in% c("41", '46', "29", "15", "23", "34") ~ "CD8_Tem",
                                               harmony_clusters_reso3_full %in% c("8", '40') ~ "CTL",
                                               harmony_clusters_reso3_full %in% c("12") ~ "MAIT",
                                               harmony_clusters_reso3_full %in% c("7", "32") ~ "gdT",
                                               harmony_clusters_reso3_full %in% c("1", "5") ~ "CD8_Temra",
                                               harmony_clusters_reso3_full %in% c("11") ~ "CD8_CD56_TNF",
                                               harmony_clusters_reso3_full %in% c("26", "17") ~ "CD8_CD56", 
                                               harmony_clusters_reso3_full %in% c("10", "4", '45', '27') ~ "NK_CD56bright",
                                               harmony_clusters_reso3_full %in% c("13", "28", "0", '24', '14', '44', '6', '2', "3") ~ "NK_CD56dim",
                                               harmony_clusters_reso3_full %in% c("49") ~ "ILC",
                                               harmony_clusters_reso3_full %in% c("38") ~ "Doublets_2",
                                               harmony_clusters_reso3_full %in% c("36", "22", "39", "47", "33", "42", "43", "30", "31", "37") ~ "Doublets_1")) 

#double check cluster_32 MAIT or gdT
t_nk$anno_3rd <- factor(t_nk$anno_3rd, levels = c("CD4_Tcm/naive", "Th1/Th2", "Th17", "Treg", 
                                                  "CD8_Tcm/naive", "CD8_Tem", "CTL", "CD8_Temra", "CD8_CD56", "CD8_CD56_TNF", 
                                                  "gdT", "MAIT",
                                                  "NK_CD56bright", "NK_CD56dim", "ILC",
                                                  "Doublets_1", "Doublets_2"), ordered = T)

saveRDS(
  object = t_nk,
  file = "aml_t_nk_reso3_anno.Rds",
  destdir = paste0(wd_path, '/comb_s15_rerun/aml_t_nk_reso3_anno')
)

t_nk <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_t_nk_reso3_anno/aml_t_nk_reso3_anno.Rds"))

Idents(t_nk) <- t_nk$anno_3rd


t_nk@assays$sketch <- NULL

t_nk@reductions$integrated.cca <- NULL

t_nk@reductions$harmony <- NULL

t_nk@reductions$umap <- NULL

t_nk[['RNA']]$counts <- as(object = t_nk[['RNA']]$counts, Class = 'dgCMatrix')
t_nk[['RNA']]$data <- as(object = t_nk[['RNA']]$data, Class = 'dgCMatrix')
t_nk[['RNA']]$scale.data <- as(object = t_nk[['RNA']]$scale.data, Class = 'dgCMatrix')

SeuratObject::SaveSeuratRds(t_nk, file = '/home/big/zheng_song/aml/comb_s15_rerun/aml_t_nk_reso3_anno.Rds')

t_nk <- readRDS('/home/big/zheng_song/aml/comb_s15_rerun/aml_t_nk_reso3_anno.Rds')
