#characterize t_nk cluster 2(mainly on t cells) 
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


# t_cells -----------------------------------------------------------------

t_nk <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_t_nk_reso3_anno/aml_t_nk_reso3_anno.Rds"))

Idents(t_nk) <- t_nk$anno_3rd

clusters <- c("CD4_Tcm/naive", "Th1/Th2", "Th17", "Treg", 
  "CD8_Tcm/naive", "CD8_Tem", "CTL", "CD8_Temra", "CD8_CD56", "CD8_CD56_TNF", 
  "gdT", "MAIT",
  "NK_CD56bright", "NK_CD56dim", "ILC",
  "Doublets_1", "Doublets_2")

clusters <- c("Th1/Th2", "Th17", "Treg", 
              "CD8_Tcm/naive", "CD8_Tem", "CTL", "CD8_Temra", "CD8_CD56", "CD8_CD56_TNF", 
              "gdT", "MAIT",
              "NK_CD56bright", "NK_CD56dim", "ILC",
              "Doublets_1", "Doublets_2")

markers_list <- list()
markers_list[["CD4_Tcm/naive"]] <- markers

for (i in clusters) {
  print(i)
  markers <- FindMarkers(t_nk, ident.1 = i, group.by = "anno_3rd", assay = "RNA", verbose = T, max.cells.per.ident = 200)
  markers %<>% arrange(-avg_log2FC)
  markers_list[[i]] = markers
  #write.csv(markers, file = paste0("/home/big/zheng_song/aml/t_nk/", i, "_deg.csv"))
}

saveRDS(object = markers_list, file = paste0(wd_path, "/comb_s15_rerun/t_nk_DEGs.Rds"))

markers_list <- readRDS(paste0(wd_path, "/comb_s15_rerun/t_nk_DEGs.Rds"))

gene <- c("CD4", "TSHZ2","TCF7", "CCR7", "SELL", "LEF1",
          "CD28", "ICOS", "TBX21", "GATA3", 
          "RORC","CCR6", "PTPN13",  
          "FOXP3", "IL2RA", "CTLA4", "SLC16A10", "TIGIT",
          "CD8A", "CD8B",
          "IFNG-AS1","NELL2", "DUSP4",
          "CD5", "CD6", "MCOLN2", "ATG2A",
          "GZMB", "GZMH", "KLRG1", "TOX", "IFNG",
          "GNLY", "NCAM1", "IKZF2", "KLRC3", "KLRC4", "KLRC2", 
          "FCGR3A","TNF", "JUN",
          "TRGV9", "TRDV2", "TRDC", "TRGC1",
          "SLC4A10", "KLRB1",  
          "CD3E", "CD3G","NCAM1", "ZBTB16","NCR1", "PRF1", "KLRD1", "NKG7",
          "TNFRSF11A", "IL18R1", "KIR2DL4", "KLRC2",
          "IL7R","IL1R1", 
          "CD14", "CD163", "S100A8", "S100A9", "S100Z", "HLA-DRB1", "HLA-DRA")#"CSF3R", "TNFAIP2",, "CD86", "CSF1R")

(Seurat::DotPlot(t_nk, unique(gene), 
                 assay = "RNA",cols = c("#9fc5e8","#990000"), dot.scale = 2.5)+
    guides(size = F,color = guide_colorbar(title = "Expression", title.position = "right"))+
    #scale_y_discrete()+
    #scale_x_discrete()+
    dotplot_theme + 
    coord_flip()+
    theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5))
    #       axis.text.y.left = element_text(angle = 180) #hjust = 0.9, vjust = 0.9)
          )%>% 
  figsave("/comb_s15_rerun/t_nk/dotplot_marker_genes_main.pdf", w = 100, h = 180)


# (Feature_rast5(t_nk, "anno_3rd", do.label = T, colorset = "um", sz = 0.02, noaxis = F, titlesize = 0, navalue = "black") + 
#     xlab("UMAP_1") +
#     ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/t_nk/umap_anno_3rd.pdf", w = 150, h = 100) 



umap.colors <- c("#FFC133", "#990000", "#336699", "#e296ad", "#009999", "#DAA520", "#96ade2", "#a98a7b", "#6495ED", "#FFA07A", "#9370DB", "#669999", "#993300", "#666666", "#333366", "#95a674", "#BA55D3")

(Feature_rast5(t_nk, "anno_3rd", do.label = F, colorset = "um", sz = 0.04, noaxis = F, titlesize = 0, sort = T) + 
    umap_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/t_nk/umap_anno_3rd_2.pdf", w = 100, h = 75) 

Idents(t_nk) <- t_nk$anno_3rd

t_nk_markers <- FindAllMarkers(t_nk, assay = "RNA", verbose = T, logfc.threshold = 0.1, 
                               min.pct = 0.02, max.cells.per.ident = 500, test.use = "t")

t_nk_markers %>% write.csv("t_nk_all_markers_.csv")

# focus on gdT ------------------------------------------------------------
t_cells <- subset(t_nk, subset = anno_3rd %in% c("CD4_Tcm/naive", "Th1/Th2", "Th17", "Treg", 
                                                 "CD8_Tcm/naive", "CD8_Tem", "CTL", "CD8_Temra", "CD8_CD56", "CD8_CD56_TNF", 
                                                 "gdT", "MAIT"))

abTR <- grep("^TRAV|^TRBV|^TRAC|^TRBC" , rownames(t_cells@assays$RNA), value = T)

TRD <- grep("^TRDV|^TRDC" , rownames(t_cells@assays$RNA), value = T)

#gene module cal
t_cells <- AddModuleScore(t_cells, 
                       features = list(c(abTR)),
                       name = "alpha_beta_score_tr", assay = "RNA")

t_cells <- AddModuleScore(t_cells, 
                       features = list(c(TRD)),
                       name = "gamma_delta_score_trd", assay ="RNA")                                              

Feature_rast5(t_cells, "gamma_delta_score_trd1", do.label = F, color_grd = "grd", sz = 0.3, noaxis = F, titlesize = 0, sort = T)

t_cells@meta.data %<>% mutate(gdt_cal = case_when(anno_3rd!="gdT"&gamma_delta_score_trd1>=0.1&alpha_beta_score_tr1< -0.0045 ~ "gdT",
                                                  anno_3rd!="gdT"&gamma_delta_score_trd1<0.1 ~ "abT",
                                                  anno_3rd!="gdT"&gamma_delta_score_trd1>=0.1&alpha_beta_score_tr1>= -0.0045 ~ "dpT",
                                                  anno_3rd == "gdT" ~ "gdT"))

Feature_rast5(t_cells, "gdt_cal", do.label = F, colorset = "gg", sz = 1, noaxis = F, titlesize = 0, sort = T, facets = "gdt_cal")


Feature_rast(t_cells, d1 = "alpha_beta_score_tr1", d2 = "gamma_delta_score_trd1", g = "nCount_RNA", color_grd = "grd",noaxis = F, axis.number = T) +
  geom_hline(yintercept = 0.1)+
  geom_vline(xintercept = -0.0045)

# here we employed a loose cutoff to gain more gdT cells
t_cells@meta.data %<>% mutate(gdt_cal = case_when(anno_3rd!="gdT"&gamma_delta_score_trd1>=0.1&alpha_beta_score_tr1< -0.0045 ~ "gdT",
                                                  anno_3rd!="gdT"&gamma_delta_score_trd1<0.1 ~ "abT",
                                                  anno_3rd!="gdT"&gamma_delta_score_trd1>=0.1&alpha_beta_score_tr1>= -0.0045 ~ "gdT",
                                                  anno_3rd == "gdT" ~ "gdT"))

(Feature_rast5(t_cells, "gamma_delta_score_trd1", do.label = F, color_grd = "grd", sz = 0.1, noaxis = F, titlesize = 0) + 
    umap_theme +
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/t_nk/umap_t_cells_trd.pdf", w = 80, h = 75) 

t_cells$gdt_cal <- factor(t_cells$gdt_cal, levels = c("gdT", "abT"))

umap.colors <- c("#9370DB","grey85")

library(gghighlight)


(Feature_rast5(t_cells,g = "gdt_cal" ,sz = 0.1, do.label = F, labelsize = 8,noaxis = F, titlesize = 0, colorset = "gg")+
    gghighlight(gdt_cal == "gdT", label_key = gdt_cal, use_direct_label = F)+
    theme(plot.margin = margin(5, 3, 3, 5, "pt"),
          axis.text.y.left = element_blank(),
          axis.title.x = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt")),
          axis.text.x.bottom = element_blank(),
          axis.title.y = element_text(size = 8, color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt")))+
    scale_color_manual(labels = c("γδ T", "αβ T"), values = c("#9370DB","grey85"))+
    theme(legend.position = "right",
          legend.title = element_blank(), 
          legend.text = element_text(size = 8, color = "black" ,face = "plain",family = "Arial"),
          panel.background=element_rect(fill="transparent",colour=NA),
          #plot.background=element_rect(fill="transparent",colour=NA),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          axis.line.y.left = element_line(linewidth = 0.4, colour = "black"),
          axis.line.x.bottom = element_line(linewidth = 0.4, colour = "black"),
          legend.key.width = unit(6, "pt"),
          legend.spacing.y = unit(4, "pt"),
          legend.box.spacing = unit(0, "pt"))) %T>% figsave("comb_s15_rerun/t_nk/umap_t_cells_gdt_cal.pdf", w = 90, h = 75) 

cln = 2
ViolinPlot(t_cells, c("alpha_beta_score_tr1", "gamma_delta_score_trd1"), group.by = "anno_3rd", colors = umap.colors)

saveRDS(
  object = t_cells,
  file = "aml_t_cells_reso3_anno.Rds",
  destdir = paste0(wd_path, '/comb_s15_rerun/aml_t_cells_reso3_anno')
)

t_cells <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_t_cells_reso3_anno/aml_t_cells_reso3_anno.Rds"))

t_cells@assays$sketch <- NULL

t_cells@reductions$integrated.cca <- NULL

t_cells@reductions$harmony <- NULL

t_cells@reductions$umap <- NULL

t_cells[['RNA']]$counts <- as(object = t_cells[['RNA']]$counts, Class = 'dgCMatrix')
t_cells[['RNA']]$data <- as(object = t_cells[['RNA']]$data, Class = 'dgCMatrix')
t_cells[['RNA']]$scale.data <- as(object = t_cells[['RNA']]$scale.data, Class = 'dgCMatrix')

SeuratObject::SaveSeuratRds(t_cells, file = '/home/big/zheng_song/aml/comb_s15_rerun/aml_t_cells_reso3_anno.Rds')
# we delete the "destdir = paste0(wd_path, '/comb_s15_rerun/aml_t_cells_reso3_anno')"

# subset gdT based on cal
gd_t <- subset(t_cells, gdt_cal == "gdT")
# 13256 cells

gd_t <- FindNeighbors(gd_t, reduction = "harmony.full", dims = 1:30)
gd_t <- FindClusters(gd_t, resolution = 3, cluster.name = "harmony_clusters_reso3_full")
gd_t <- RunUMAP(gd_t, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAPfull_")

Feature_rast5(gd_t, colorset = "gg", g = "harmony_clusters_reso3_full", noaxis = F, do.label = T, sz = 0.3) %T>% 
  figsave("/comb_s15_rerun/gdt/umap_clusters_3.pdf", w = 120, h = 75)

# gene module from Tan_SI_2019 
gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)

gene_cluster$GM_D

for (i in gc_name) {
  gd_t %<>%  AddModuleScore(features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
}

Feature_rast5(gd_t, "GM_A1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_A1")
c("22") ~ "naive"

Feature_rast5(gd_t, "GM_B1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_B1")
c("26") ~ "gdT_proliferating" # 26 is more proflirating

Feature_rast5(gd_t, "GM_C1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_C1")
c("11", "20", "29") ~ "gdT_immature"

Feature_rast5(gd_t, "GM_D1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_D1")
c("5", "8", "12", "16", "18", "19", "28") ~ "gdT_type3"

Feature_rast5(gd_t, "GM_E1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_E1")
# not found

Feature_rast5(gd_t, "TRDV2", color_grd = "grd")
Feature_rast5(gd_t, "GM_F1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_F1")
c("23", "6", "0", "1", "4", "9", "10", "24", "15", "13", "14", "27", "2", "3", "25")

Feature_rast5(gd_t, "TRDV1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_G1")
c("6", "7", "21") ~ "Vd1"

Feature_rast5(gd_t, "GM_H1", color_grd = "grd")
ViolinPlot(gd_t, g = "GM_H1")
c("11", "20", "24", "29") ~ "gdT_acute_activation"

Feature_rast5(gd_t, g = c("TRDV1", "TRDV2", "TRDV3", "TRDC", "TRGC1", "TRGV9"), ncol = 3, assay = "RNA", color_grd = "grd")
c("7", "17", "22", "26") ~ "TRDV1"

gd_t@meta.data %<>% mutate(anno_4th= case_when(harmony_clusters_reso3_full %in% c("22") ~ "gdT_naive",
                                               harmony_clusters_reso3_full %in% c("26") ~ "gdT_proliferating",
                                               harmony_clusters_reso3_full %in% c("11") ~ "Vd2_eff_early",
                                               harmony_clusters_reso3_full %in% c("5", "8", "12", "16", "18", "19", "28") ~ "Vd2_eff_type3",
                                               harmony_clusters_reso3_full %in% c("29") ~ "gdT_immature",
                                               harmony_clusters_reso3_full %in% c("23", "6", "0", "1", "4", "10", "2", "25", "20") ~ "Vd2_eff_IFNG",
                                               harmony_clusters_reso3_full %in% c("9", "24", "15", "13", "14", "27","3") ~ "Vd2_eff_GZMK",
                                               harmony_clusters_reso3_full %in% c("17") ~ "Vd3_eff",
                                               harmony_clusters_reso3_full %in% c("7") ~ "Vd1_eff",
                                               harmony_clusters_reso3_full %in% c("21") ~ "Non_gdT")) 

gd_t$anno_4th <- factor(gd_t$anno_4th, levels = c("gdT_naive", "gdT_proliferating", "gdT_immature", 
                                                  "Vd2_eff_early", "Vd2_eff_IFNG", "Vd2_eff_GZMK", "Vd2_eff_type3",
                                                  "Vd1_eff", "Vd3_eff", "Non_gdT"))

umap.colors <- c("#990000", "#333366", "#336699", "#e296ad", "#009999", "#DAA520", "#95a674", "#a98a7b", "#6495ED", "#FFA07A", "#9370DB", "#669999", "#993300", "#666666",   "#BA55D3")

(Feature_rast5(gd_t, "anno_4th", do.label = F, colorset = "um", sz = 0.4, noaxis = F, titlesize = 0, sort = T) + 
    umap_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/gdt/umap_anno_4th.pdf", w = 100, h = 75) 

gd_t_marker <- FindMarkers(object = gd_t, ident.1 = "Vd2_eff", ident.2 = "Vd2_eff_IFNG", group.by = "anno_4th", assay = "RNA")

gd_t_marker %>% arrange(-avg_log2FC)

saveRDS(
  object = gd_t,
  file = "aml_gd_t_reso3_anno.Rds",
  destdir = paste0(wd_path, '/comb_s15_rerun/aml_gd_t_reso3_anno')
)

gd_t <- readRDS(file = '/home/big/zheng_song/aml/comb_s15_rerun/aml_gd_t_reso3_anno/aml_gd_t_reso3_anno.Rds')

gd_t@assays$sketch <- NULL

gd_t@reductions$integrated.cca <- NULL

gd_t@reductions$harmony <- NULL

gd_t@reductions$umap <- NULL

gd_t[['RNA']]$counts <- as(object = gd_t[['RNA']]$counts, Class = 'dgCMatrix')
gd_t[['RNA']]$data <- as(object = gd_t[['RNA']]$data, Class = 'dgCMatrix')
gd_t[['RNA']]$scale.data <- as(object = gd_t[['RNA']]$scale.data, Class = 'dgCMatrix')

SeuratObject::SaveSeuratRds(gd_t, file = '/home/big/zheng_song/aml/comb_s15_rerun/aml_gd_t_reso3_anno.Rds')

gd_t <- readRDS(file = '/home/big/zheng_song/aml/comb_s15_rerun/aml_gd_t_reso3_anno.Rds')

gd_t@meta.data

# gdt cell pct of total cell number
gdt_count <- gd_t@meta.data %>% group_by(donor) %>% dplyr::select(sample, tp) %>% dplyr::count(sample)
colnames(gdt_count) <- c("donor", "sample", "gdT_counts")

t_count <- t_cells@meta.data %>% group_by(donor) %>% dplyr::select(sample) %>% dplyr::count(sample)
colnames(t_count) <- c("donor", "sample", "total_T_counts")

gdt_pct <- left_join(gdt_count, t_count, by = c("sample", "donor"))

gdt_pct$gdT_pct <- gdt_pct$gdT_counts/gdt_pct$total_T_counts
colnames(gdt_pct) <- c("donor", "sample", "gdT_counts", "total_T_counts", "gdT_pct")

gdt_pct$tp <- str_extract(gdt_pct$sample, "Transplant|D30|D100|Before_SCT")

gdt_pct$tp <- factor(gdt_pct$tp, levels = c("Before_SCT", "Transplant", "D30", "D100"))

(ggplot(data = gdt_pct, aes(x = donor, y = gdT_pct, fill = tp))+
    geom_bar(stat="identity",position='dodge', width = 0.7) + 
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
    scale_fill_manual(values = c("#800000",  "#b5a699", "#67755c","#058789")) + violin_theme +
    ylab("% of gdT cells (gated on total T cells)") + xlab("Donor") + theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1)) + labs(fill = "Timepoint")) %>% 
  figsave("comb_s15_rerun/gdt/gdt_pct_on_t.pdf", w = 125, h = 75)

sum <- gdt_pct %>% group_by(tp) %>% summarise(mean_pct = mean(gdT_pct), sd_pct = sd(gdT_pct), n_pct = n(), se = sd_pct/sqrt(n_pct))

gdT_pct <- left_join()

(ggplot(data = gdt_pct,aes(x = factor(tp), y = gdT_pct, fill = factor(tp))) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
  geom_line(aes(group = donor), color = "grey") +  
  scale_fill_manual(values = c("#800000",  "#b5a699", "#67755c","#058789")) + violin_theme +
  theme(legend.position = "none") + 
  ylab("% of gdT cells (gated on total T cells)") + xlab("Timepoint") + theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1)) + labs(fill = "Timepoint")) %>% 
  figsave("comb_s15_rerun/gdt/gdt_pct_on_t_tp.pdf", w = 75, h = 75)

#gdt subsets pct
# doublets cell pct of total cell number
gdt_count <- gd_t@meta.data %>% group_by(tp) %>% dplyr::select(sample, anno_4th, tp) %>% dplyr::count(anno_4th)
colnames(gdt_count) <- c("tp", "anno_4th", "subset_counts")

tp_count <- gd_t@meta.data %>% group_by(tp) %>% dplyr::select(sample, anno_4th, tp) %>% dplyr::count(tp)

colnames(tp_count) <- c("tp", "total_gdt_counts")

gdt_pct <- left_join(gdt_count, tp_count, by = c("tp" ))

gdt_pct$gdt_subset_pct <- gdt_pct$subset_counts/gdt_pct$total_gdt_counts*100
colnames(gdt_pct) <- c("tp", "anno_4th", "subset_counts", "total_gdt_counts", "gdt_subset_pct")
gdt_pct$tp <- factor(gdt_pct$tp, levels = c("Before_SCT", "Transplant", "D30", "D100"))

(ggplot(data = gdt_pct, aes(x = tp, y = gdt_subset_pct, fill = anno_4th))+
    geom_bar(stat="identity",position = "stack", width = 0.7) + 
    scale_fill_manual(values = umap.colors) + violin_theme +
    theme(legend.title = element_blank())+
    ylab("% of subsets(gated on total gdT)") + xlab("Timepoint") + theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1)) + labs(fill = "Timepoint")) %>% 
  figsave("comb_s15_rerun/gdt/subset_pct_on_gdt.pdf", w = 100, h = 80)


# tcr_sharning ------------------------------------------------------------
library(viridis)
library(ggalluvial)
alluvium_theme <- theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
                        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                        panel.grid.major = element_blank(), #remove major gridlines
                        panel.grid.minor = element_blank(),
                        axis.line.x.bottom = element_line(),
                        axis.title = element_text(size = 10),
                        axis.text = element_text(size = 10))

trans_gdt <- subset(gd_t, donor != "MRD003")

# TCR sharing strata is donor
TCRG <- trans_gdt@meta.data %>% filter(func_TRG == "functional")%>% filter(!is.na(v_TRG) & cdr3_freq_TRG > 0) %>%
  group_by(cdr3_aa_TRG, tp, donor) %>%  summarise(cdr3_freq_TRG = n(), .groups = "keep") %>% arrange(cdr3_freq_TRG, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRG = factor(cdr3_aa_TRG, unique(cdr3_aa_TRG)))

TCRG$donor_strata <- paste0(TCRG$donor, "_", TCRG$cdr3_aa_TRG)
TCRG$tp <- factor(TCRG$tp, levels = c("Transplant", "D30", "D100"))

custom_colors <- viridis_pal(alpha = 0.5, option = "E", begin = 0.3, end = 0.7)(length(unique_donors))
unique_donors <- unique(TCR$donor)

donor_color_mapping <- data.frame(donor = unique_donors, custom_color = custom_colors)

donor_color_mapping$donor <- factor(donor_color_mapping$donor, levels = c("MRD002", "MRD004", "MRD007", "MRD008", "MRD009", "MRD010", "MRD012", "MRD013",
                                                                          "MRD016"))

(ggplot(TCRG,
        aes(x = tp, stratum = donor, alluvium = donor_strata,
            fill = cdr3_aa_TRG, label = donor)) +
    geom_alluvium(show.legend = F) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "grey", show.legend = F) +
    geom_stratum(aes(fill = donor), show.legend = F, 
                 fill = c(rev(c("#B0A47380", "#A29A7680", "#958F7880", "#88857980", "#7C7B7880", "#72727480", "#66697080", "#5A5F6D80", "#4E576C80")),
                          rev(c("#A29A7680", "#958F7880", "#88857980", "#7C7B7880", "#72727480", "#66697080", "#5A5F6D80", "#4E576C80")),
                          rev(c("#A29A7680", "#958F7880", "#88857980", "#5A5F6D80"))), 
                 alpha = 1)+  # Remove show.legend = TRUE
    #scale_color_manual(values = donor_color_mapping$custom_color) +
    #theme(legend.position = "right")+
    ggtitle("TCR Gamma sharing between timepoints") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8))+
    geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    #theme(text = element_rect(fill = "transparent", colour = "transparent"))+
    NoLegend() + xlab("") + ylab("frequence")+
    alluvium_theme
)%T>% 
  figsave("comb_s15_rerun/gdt/trg_sharing_donor.pdf", w = 100, h = 100)

c("#B0A47380", "#A29A7680", "#958F7880", "#88857980", "#7C7B7880", "#72727480", "#66697080", "#5A5F6D80", "#4E576C80")

# TCR sharing strata is donor
TCRD <- trans_gdt@meta.data %>% filter(!is.na(v_TRD) & cdr3_freq_TRD > 0) %>%
  group_by(cdr3_aa_TRD, tp, donor) %>%  summarise(cdr3_freq_TRD = n(), .groups = "keep") %>% arrange(cdr3_freq_TRD, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRD = factor(cdr3_aa_TRD, unique(cdr3_aa_TRD)))

TCRD$donor_strata <- paste0(TCRD$donor, "_", TCRD$cdr3_aa_TRD)
TCRD$tp <- factor(TCRD$tp, levels = c("Transplant", "D30", "D100"))

custom_colors <- viridis_pal(alpha = 0.5, option = "E", begin = 0.3, end = 0.7)(length(unique_donors))
unique_donors <- unique(TCRD$donor)

donor_color_mapping <- data.frame(donor = unique_donors, custom_color = custom_colors)

donor_color_mapping$donor <- factor(donor_color_mapping$donor, levels = c("MRD002", "MRD004", "MRD007", "MRD008", "MRD009", "MRD010", "MRD012", "MRD013",
                                                                          "MRD016"))

donor_color_mapping %<>% arrange(donor)


(ggplot(TCRD,
        aes(x = tp, stratum = donor, alluvium = donor_strata,
            fill = cdr3_aa_TRD, label = donor)) +
    geom_alluvium(show.legend = F) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "grey", show.legend = F) +
    geom_stratum(aes(fill = donor), show.legend = F, 
                 fill = c(rev(c("#B0A47380", "#A29A7680", "#958F7880", "#88857980", "#7C7B7880", "#72727480", "#5A5F6D80", "#4E576C80")),
                          rev(c("#A29A7680", "#958F7880", "#88857980", "#7C7B7880", "#72727480", "#66697080", "#4E576C80")),
                          rev(c("#A29A7680", "#958F7880", "#88857980", "#5A5F6D80"))), 
                 alpha = 1)+  # Remove show.legend = TRUE
    #scale_color_manual(values = donor_color_mapping$custom_color) +
    #theme(legend.position = "right")+
    ggtitle("TCR Delta sharing between timepoints") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8))+
    geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    #theme(text = element_rect(fill = "transparent", colour = "transparent"))+
    NoLegend() + xlab("") +ylab("frequence")+
    alluvium_theme
)%T>% 
  figsave("comb_s15_rerun/gdt/trd_sharing_donor.pdf", w = 100, h = 100)


# cross donor
# TCR sharing strata is donor
TCRD <- trans_gdt@meta.data %>% filter(!is.na(v_TRD) & cdr3_freq_TRD > 0) %>% dplyr::select(cdr3_aa_TRD, tp, donor) %>% 
  group_by(cdr3_aa_TRD, donor) %>%  summarise(cdr3_freq_TRD = n(), .groups = "keep") %>% arrange(cdr3_freq_TRD, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRD = factor(cdr3_aa_TRD, unique(cdr3_aa_TRD)))

TCRD$

TCRD <- trans_gdt@meta.data %>% filter(!is.na(v_TRD) & cdr3_freq_TRD > 0) %>% 
  group_by(cdr3_aa_TRD, donor) %>%  summarise(cdr3_freq_TRD = n(), .groups = "keep") %>% arrange(cdr3_freq_TRD, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(cdr3_aa_TRD = factor(cdr3_aa_TRD, unique(cdr3_aa_TRD)))


TCRD <- trans_gdt@meta.data %>% filter(!is.na(v_TRD) & cdr3_freq_TRD > 0) %>% 
  group_by(cdr3_aa_TRD, donor) %>%  summarise(donor = n(), .groups = "keep") %>% arrange(donor, .by_group = TRUE) %>%
  ungroup() %>% 
  mutate(donor = factor(donor, unique(donor)))


TCRD$donor_strata <- paste0(TCRD$donor, "_", TCRD$cdr3_aa_TRD)
TCRD$tp <- factor(TCRD$tp, levels = c("Transplant", "D30", "D100"))

custom_colors <- viridis_pal(alpha = 0.5, option = "E", begin = 0.3, end = 0.7)(length(unique_donors))
unique_donors <- unique(TCRD$donor)

donor_color_mapping <- data.frame(donor = unique_donors, custom_color = custom_colors)

donor_color_mapping$donor <- factor(donor_color_mapping$donor, levels = c("MRD002", "MRD004", "MRD007", "MRD008", "MRD009", "MRD010", "MRD012", "MRD013",
                                                                          "MRD016"))

donor_color_mapping %<>% arrange(donor)


(ggplot(TCRD,
        aes(x = donor, stratum = donor, alluvium = donor_strata,
            fill = cdr3_aa_TRD, label = donor)) +
    geom_alluvium(show.legend = F) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "grey", show.legend = F) +
    geom_stratum(aes(fill = donor), show.legend = F, alpha = 1) #+  # Remove show.legend = TRUE
    #scale_color_manual(values = donor_color_mapping$custom_color) +
    #theme(legend.position = "right")+
    ggtitle("TCR Delta sharing between timepoints") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8))+
    geom_label(stat = "stratum", aes(label = donor), label.size = 0, size = 2, fill = NA)+
    #theme(text = element_rect(fill = "transparent", colour = "transparent"))+
    NoLegend() + xlab("") +ylab("frequence")+
    alluvium_theme
)%T>% 
  figsave("comb_s15_rerun/gdt/trd_sharing_donor.pdf", w = 100, h = 100)

install.packages(divo)


# shared TRDV  ------------------------------------------------------------
trd_df <- select(gd_t@meta.data, c("tp", "donor", "cdr3_aa_TRD", "v_TRD"))

trd_df_bsct <- trd_df %>% filter(!is.na(cdr3_aa_TRD)) %>% filter(!is.na(v_TRD))%>% filter(tp == "Before_SCT") 

trd_df_tran <- trd_df %>% filter(!is.na(cdr3_aa_TRD)) %>% filter(!is.na(v_TRD))%>% filter(tp == "Transplant") 

trd_df_d30 <- trd_df  %>% filter(!is.na(cdr3_aa_TRD)) %>% filter(!is.na(v_TRD))%>% filter(tp == "D30")

trd_df_d100 <- trd_df  %>% filter(!is.na(cdr3_aa_TRD)) %>% filter(!is.na(v_TRD))%>% filter(tp == "D100")

inner_join(trd_df_tran, y = trd_df_d30, by="cdr3_aa_TRD")
inner_join(trd_df_d30, y = trd_df_d100, by="cdr3_aa_TRD")

inner_join(trd_df_bsct, y = trd_df_d30, by="cdr3_aa_TRD")
inner_join(trd_df_bsct, y = trd_df_d100, by="cdr3_aa_TRD")

install.packages("immunarch")           # Install the package
library(immunarch)

