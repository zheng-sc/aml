#characterize unknown/doublets cluster 
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

t_nk <- readRDS('/home/big/zheng_song/aml/comb_s15_rerun/aml_t_nk_reso3_anno.Rds')

Idents(t_nk) <- t_nk$anno_3rd

# doublets ----------------------------------------------------------------
umap.colors <- c("#FFC133", "#990000", "#336699", "#e296ad", "#009999", "#DAA520", "#96ade2", "#a98a7b", "#6495ED", "#FFA07A", "#9370DB", "#669999", "#993300", "#666666", "#333366", "#95a674", "#BA55D3")

(ViolinPlot(t_nk, g = "nFeature_RNA", group.by = "anno_3rd", colors = umap.colors, sz = 0) +
    theme(plot.margin = margin(5, 3, 1, 5, "pt"),
          #plot.title = element_blank(),
          axis.title.y = element_text(size = 8, color = "black" ,face = "plain", family = "Arial", margin = margin(0,0,0,0,"pt")),
          axis.text.y.left = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt")),
          axis.title.x = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(0,0,0,0,"pt")),
          axis.text.x.bottom = element_text(size = 8, color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt"), 
                                            angle = 30, hjust = 0.9, vjust = 0.9),
          legend.position = "none",
          panel.background=element_rect(fill="transparent",colour=NA),
          panel.grid = element_blank(),
          axis.line.y.left = element_line(linewidth = 0.4, colour = "black"),
          axis.line.x.bottom = element_line(linewidth = 0.4, colour = "black"))+
    ylab("nFeature_RNA")) %>% 
  figsave("/comb_s15_rerun/doublets/vln_nfeature_counts.pdf", w = 120, h = 75)

t_nk <- AddModuleScore(t_nk, 
                          features = list(c("CD3E", "CD3D", "CD3G", "TRBC1", "TRBC2", "TRAC", "TRDC", "CD4", "CD8A", "CD8B", "SLC4A10", "NCAM1", "CD28")),
                          name = "t_gm", assay = "RNA", slot = "counts", nbin = 10)

ViolinPlot(t_nk, g = "t_gm1")

t_nk <- AddModuleScore(t_nk, 
                          features = list(c("CD14", "CD163", "LYZ", "VCAN", "S100A8", "S100A9")),
                          name = "mono_gm", assay ="RNA", slot = "counts", nbin = 10)                                              

ViolinPlot(t_nk, g = "mono_gm1")

t_nk@meta.data %<>% mutate(Doublets = case_when(anno_3rd == "Doublets_1" ~ "Doublets_1",
                                                  anno_3rd == "Doublets_2" ~ "Doublets_2",
                                                  TRUE ~ NA))

Feature_rast(t_nk, d1 = "t_gm1", d2 = "mono_gm1", g = "Doublets", colorset = "um",noaxis = F, axis.number = T) +
  geom_hline(yintercept = 0.1)+
  geom_vline(xintercept = -0.0045)

# to clearly show the doublets T:Mono complex, we go back to all cells
all_cells <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_all_cells_cca_harmony_mrd008_tcr_gsea_reso2.5/aml_all_cells_cca_harmony_mrd008_tcr_gsea_reso2.5.Rds"))
  
all_cells <- AddModuleScore(all_cells, 
                       features = list(c("CD3E", "CD3D", "CD3G", "TRBC1", "TRBC2", "TRAC", "TRDC", "CD4", "CD8A", "CD8B", "SLC4A10", "NCAM1", "CD28")),
                       name = "t_gm", assay = "RNA")

#ViolinPlot(t_nk, g = "t_gm1")

all_cells <- AddModuleScore(all_cells, 
                       features = list(c("CD14", "CD163", "LYZ", "VCAN", "S100A8", "S100A9")),
                       name = "mono_gm", assay ="RNA")                                              

#ViolinPlot(t_nk, g = "mono_gm1")

all_cells@meta.data %<>% mutate(Doublets = case_when(harmony_clusters_reso2.5_full %in% c("22", "16", "21", "18", "33", "7", "8", "30",
                                                                                     "31", "35", "1", "11", "13", "14", "10", "3", "19",
                                                                                     "24", "5") ~ "Mono&DC cells",
                                                harmony_clusters_reso2.5_full %in% c("25", "0", "17", "23", "2", "4", "28", "20", "12",
                                                                                     "9", "40") ~ "T&NK cells",
                                                harmony_clusters_reso2.5_full %in% c("38", "29") ~ "Doublets",
                                                TRUE ~ NA))
umap.colors <- c("#990000", "#a98a7b","#6B5695")

library(gghighlight)

(Feature_rast(all_cells, d1 = "t_gm1", d2 = "mono_gm1", g = "Doublets", colorset = "um",noaxis = F, axis.number = T, 
              do.label = F, sz = 0.03, titlesize = 0, navalue = "transparent") + 
    #gghighlight(Doublets == "Doublets", label_key = Doublets, use_direct_label = F, )+
    xlab("T_cell_gene_module")+
    ylab("Myeloid_cell_gene_module")+
    theme(plot.margin = margin(5, 3, 3, 5, "pt"),
          axis.text.y.left = element_text(size = 8,color = "black" ,face = "plain", family = "Arial"),
          axis.title.x = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt")),
          axis.text.x.bottom = element_text(size = 8,color = "black" ,face = "plain", family = "Arial"),
          axis.title.y = element_text(size = 8, color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt"), 
                                      angle = 90),
          legend.position = "right",
          legend.title = element_blank(), 
          legend.text = element_text(size = 8, color = "black" ,face = "plain",family = "Arial"),
          panel.background=element_rect(fill="transparent",colour=NA),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          axis.line.y.left = element_line(linewidth = 0.4, colour = "black"),
          axis.line.x.bottom = element_line(linewidth = 0.4, colour = "black"),
          legend.spacing.y = unit(0, "pt"),
          legend.key.width = unit(6, "pt"),
          legend.box.spacing = unit(0, "pt")) ) %>% 
  figsave("/comb_s15_rerun/doublets/t_mono_gm.pdf", w = 100, h = 75)

rm(all_cells)






doublets <- subset(t_nk, subset = anno_3rd %in% c("Doublets_1", "Doublets_2"))
#19349 cells


#doublets <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_doublets_reso3_anno/aml_doublets_reso3_anno.Rds"))

doublets <- FindNeighbors(doublets, reduction = "harmony.full", dims = 1:30)
doublets <- FindClusters(doublets, resolution = 2, cluster.name = "harmony_clusters_reso2_full")
doublets <- RunUMAP(doublets, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAPfull_")

Feature_rast5(doublets, colorset = "gg", g = "harmony_clusters_reso2_full", noaxis = F, do.label = T, sz = 0.3) %T>% 
  figsave("/comb_s15_rerun/doublets/umap_clusters_reso2.pdf", w = 120, h = 75)

doublets@assays$sketch <- NULL

doublets@reductions$integrated.cca <- NULL

doublets@reductions$harmony <- NULL

doublets@reductions$umap <- NULL

doublets[['RNA']]$counts <- as(object = doublets[['RNA']]$counts, Class = 'dgCMatrix')
doublets[['RNA']]$data <- as(object = doublets[['RNA']]$data, Class = 'dgCMatrix')
doublets[['RNA']]$scale.data <- as(object = doublets[['RNA']]$scale.data, Class = 'dgCMatrix')

SeuratObject::SaveSeuratRds(doublets, destdir = paste0(wd_path, "/comb_s15_rerun/aml_doublets_reso2_anno.Rds"))

# saveRDS(
#   object = doublets,
#   file = "aml_doublets_reso2_anno.Rds",
#   destdir = paste0(wd_path, '/comb_s15_rerun/aml_doublets_reso2_anno')
# )

save.image("~/aml/doublets.RData")

# anno doublets -----------------------------------------------------------
Feature_rast5(doublets, colorset = "gg", g = "sample", noaxis = F, do.label = T, sz = 0.3) %T>% 
  figsave("/comb_s15_rerun/doublets/umap_clusters_reso2_2.pdf", w = 120, h = 75)

Feature_rast5(doublets, colorset = "gg", g = c("CD4", "CD8A", "RORC", "CD14", "CD3E", "SLC4A10", "TSHZ2", "NCAM1", "TRDC", "TRDV1", "TRDV2", "FOXP3"), 
              noaxis = F, do.label = T, sz = 0.6, color_grd = "grd", ncol = 4) #%T>% 
  
doublets@meta.data %<>% mutate(anno_4th= case_when(harmony_clusters_reso2_full %in% c("14", "8") ~ "CD4_CD8_2",
                                               harmony_clusters_reso2_full %in% c("9", "23") ~ "CD4_CD8_1",
                                               harmony_clusters_reso2_full %in% c("24") ~ "CD4_MAIT",
                                               harmony_clusters_reso2_full %in% c("17") ~ "CD4_CD8_CD56bright",
                                               harmony_clusters_reso2_full %in% c("4", "18", "22") ~ "NK_Mono_CD56bright",
                                               harmony_clusters_reso2_full %in% c("0", "13", "6", "10" , "20") ~ "NK_Mono_CD56dim",
                                               harmony_clusters_reso2_full %in% c("16") ~ "CD4_NK_Mono_CD56dim",
                                               harmony_clusters_reso2_full %in% c("12", "25") ~ "CD4_MAIT_Mono",
                                               harmony_clusters_reso2_full %in% c("2") ~ "CD4_Vd2_Mono",
                                               harmony_clusters_reso2_full %in% c("5") ~ "CD4_CD8_Mono",
                                               harmony_clusters_reso2_full %in% c("21") ~ "Treg_Mono",
                                               harmony_clusters_reso2_full %in% c("1", "11", "7", "3", "15", "26", "19") ~ "CD4_Mono")) 

doublets$anno_4th <- factor(doublets$anno_4th, levels = c("CD4_CD8_1", "CD4_CD8_2", "CD4_CD8_CD56bright", "CD4_NK_Mono_CD56dim", 
                                                  "CD4_MAIT", "CD4_MAIT_Mono", "CD4_CD8_Mono", "CD4_Vd2_Mono", "CD4_Mono", "Treg_Mono", "NK_Mono_CD56bright",
                                                  "NK_Mono_CD56dim"))

umap.colors <- c("#990000", "#333366", "#336699", "#e296ad", "#009999", "#DAA520", "#95a674", "#a98a7b", "#6495ED", "#FFA07A", "#9370DB", "#669999", "#993300", "#666666",   "#BA55D3")

(Feature_rast5(doublets, "anno_4th", do.label = F, colorset = "um", sz = 0.4, noaxis = F, titlesize = 0, sort = T) + 
    umap_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/doublets/umap_anno_4th.pdf", w = 120, h = 75) 

# doublets cell pct of total cell number
doublets_count <- doublets@meta.data %>% group_by(tp) %>% dplyr::select(sample, anno_4th) %>% dplyr::count(anno_4th)
colnames(doublets_count) <- c("tp", "anno_4th", "subset_counts")


tp_count <- doublets@meta.data %>% group_by(tp) %>% dplyr::select(sample, anno_4th, tp) %>% dplyr::count(tp)

colnames(tp_count) <- c("tp", "total_doublets_counts")

doublets_pct <- left_join(doublets_count, tp_count, by = c("tp" ))

doublets_pct$doublets_subset_pct <- doublets_pct$subset_counts/doublets_pct$total_doublets_counts*100
colnames(doublets_pct) <- c("tp", "anno_4th", "subset_counts", "total_doublets_counts", "doublets_subset_pct")

doublets_pct$tp <- str_extract(gdt_pct$sample, "Transplant|D30|D100|Before_SCT")

doublets_pct$tp <- factor(doublets_pct$tp, levels = c("Before_SCT", "Transplant", "D30", "D100"))

(ggplot(data = doublets_pct, aes(x = tp, y = doublets_subset_pct, fill = anno_4th))+
    geom_bar(stat="identity",position = "stack", width = 0.7) + 
    scale_fill_manual(values = umap.colors) + violin_theme +
    theme(legend.title = element_blank())+
    ylab("% of subsets(gated on total Doublets)") + xlab("Timepoint") + theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1)) + labs(fill = "Timepoint")) %>% 
  figsave("comb_s15_rerun/doublets/subset_pct_on_doublets.pdf", w = 125, h = 125)

(Feature_rast5(doublets, c("ICAM1", "ITGAL"), do.label = F, colorset = "um", sz = 0.4, noaxis = F, titlesize = 10, sort = T, color_grd = "grd") + 
    umap_split_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/doublets/umap_lfa1.pdf", w = 150, h = 75) 

(Feature_rast5(doublets, c("CSF1R", "CD86"), do.label = F, colorset = "um", sz = 0.4, noaxis = F, 
               titlesize = 10, sort = T, color_grd = "grd", facetcol = "tp") +
    #facet_wrap(doublets$tp,ncol = 2, dir = "h", as.table = T, nrow = 2)+
    umap_split_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/doublets/umap_csf1r.pdf", w = 150, h = 75) 

(Feature_rast5(doublets, c("CD14", "CSF3R"), do.label = F, colorset = "um", sz = 0.4, noaxis = F, titlesize = 10, sort = T, color_grd = "grd") + 
    umap_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/doublets/umap_csf3r.pdf", w = 150, h = 75) 

(Feature_rast5(doublets, c("CD34"), do.label = F, colorset = "um", sz = 0.4, noaxis = F, titlesize = 10, sort = T, color_grd = "grd") + 
    umap_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/doublets/umap_cd34.pdf", w = 75, h = 75) 

(Feature_rast5(doublets, c("ITGAM", "ARG1"), do.label = F, colorset = "um", sz = 0.4, noaxis = F, titlesize = 10, sort = T, color_grd = "grd") + 
    umap_theme+
    #geom_label_repel(aes(label = "anno_3rd")) +
    xlab("UMAP_1") +
    ylab("UMAP_2")) %T>% figsave("comb_s15_rerun/doublets/umap_itgam.pdf", w = 150, h = 75) 




doublets$tp <- factor(doublets$tp, levels = c("Before_SCT", "Transplant", "D30",  "D100"))
umap.colors <- c("#990000", "#333366", "#336699", "#e296ad", "#009999", "#DAA520", "#95a674", "#a98a7b", "#6495ED", "#FFA07A", "#9370DB", "#669999", "#993300", "#666666",   "#BA55D3")

umap.colors <- c("#990000", "#DAA520", "#669999", "#336699")

ViolinPlot(doublets, c("CSF1R", "CD86"), group.by = "tp", colors = umap.colors)


(Feature_rast5(doublets, "tp", do.label = F, colorset = "um", sz = 0.4, noaxis = F, 
               titlesize = 10, sort = T, color_grd = "grd", facetcol = "tp") + 
  facet_wrap(~tp,ncol = 2, dir = "h", as.table = T, nrow = 2) + 
  theme(plot.margin = margin(5, 3, 3, 7, "pt"), strip.background = element_rect(linewidth = 0), strip.text.x = element_blank())+
  ylab("UMAP_2")+
  xlab("UMAP_1")+
  ggtitle(label = "")+
  theme(plot.title = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(0,0,1.5,0,"pt"), hjust = 0.5),
        axis.title.y = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(0,0,0,0,"pt")),
        axis.text.y.left = element_text(size = 8,color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt")),
        axis.title.x = element_text(size = 8, color = "black" ,face = "plain",family = "Arial",margin = margin(0,0,0,0,"pt")),
        axis.text.x.bottom = element_text(size = 8, color = "black" ,face = "plain", family = "Arial", margin = margin(1.5,1,1.5,1.5,"pt")))+
  theme(legend.position = "right", legend.margin = margin(0,0,0,0,"pt"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 8,color = "black" , face = "plain",family = "Arial", margin = margin(1.5,0,0,0,"pt")),
        legend.spacing.y = unit(4, "pt"))) %>% figsave("comb_s15_rerun/doublets/umap_facet_tp.pdf", w = 120, h = 100)

# GSEA --------------------------------------------------------------------
#aggregation-based (pseudobulk) workflow
# bulk <- AggregateExpression(doublets, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("transplant_group", "donor"))
# bulk$transplant <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
# bulk$donor <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)


de_markers <- FindMarkers(doublets, ident.1 = "CD4_Mono", group.by = "anno_4th", verbose = T)

de_markers %>% write.csv(file = "/home/big/zheng_song/aml/comb_s15_rerun/doublets/CD4_Mono_DEGs.csv", row.names = T)

de_markers$gene <- rownames(de_markers)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05&avg_log2FC > 0.25, gene, "")), colour = "red", size = 3)

de_markers %>% arrange(-avg_log2FC)

# GSEA on de_markers from unknown clusters
# non-transplant vs transplant
# GSEA
library(clusterProfiler)
library(enrichplot)
library(dbplyr)
library("AnnotationHub")
#BiocManager::install("KEGGREST")

# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

organism = 'org.Hs.eg.db'
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

de_markers <- de_markers %>% dplyr::filter(avg_log2FC > 0.25) %>%  top_n(n = 500, wt = avg_log2FC)

#all_cells_marker_2 <- subset(all_cells_marker_2, rowSums(all_cells_marker_2[5] < 0.05) > 0)

de_markers$gene <- rownames(de_markers)

de_markers <-  bitr(de_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_tg <- enrichGO(de_markers$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                  maxGSSize = 500, readable = T)

(barplot(GO_tg, label_format = "text", showCategory = 10) + 
    #scale_fill_gradient() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("GO_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/go_doublets_CD4_Mono.pdf", w = 140, h = 75)

KEGG_tg <- enrichKEGG(de_markers$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10,
                      maxGSSize = 500, qvalueCutoff = 0.05, use_internal_data = F)

barplot(KEGG_tg, showCategory = 5, title = "KEGG Pathway", x = "GeneRatio")

(barplot(KEGG_tg, label_format = "text", showCategory = 10, color = ) + 
    #scale_fill_gradient(high = "slateblue4", low = "coral") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("KEGG_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/kegg_doublets_CD4_Mono.pdf", w = 140, h = 75)

# NK_Mono_CD56dim
de_markers <- FindMarkers(doublets, ident.1 = "NK_Mono_CD56dim", group.by = "anno_4th", verbose = T)

de_markers %>% write.csv(file = "/home/big/zheng_song/aml/comb_s15_rerun/doublets/NK_Mono_CD56dim_DEGs.csv", row.names = T)
de_markers <-  read.csv("/home/big/zheng_song/aml/comb_s15_rerun/doublets/NK_Mono_CD56dim_DEGs.csv", row.names = 1)
de_markers$gene <- rownames(de_markers)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05&avg_log2FC > 0.25, gene, "")), colour = "red", size = 3)

de_markers %>% arrange(-avg_log2FC)

# GSEA on de_markers from unknown clusters
# non-transplant vs transplant
# GSEA
library(clusterProfiler)
library(enrichplot)
library(dbplyr)
library("AnnotationHub")
#BiocManager::install("KEGGREST")

# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

organism = 'org.Hs.eg.db'
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

de_markers <- de_markers %>% dplyr::filter(avg_log2FC > 0.25) %>%  top_n(n = 500, wt = avg_log2FC)

#all_cells_marker_2 <- subset(all_cells_marker_2, rowSums(all_cells_marker_2[5] < 0.05) > 0)

de_markers$gene <- rownames(de_markers)

de_markers <-  bitr(de_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_tg <- enrichGO(de_markers$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                  maxGSSize = 500, readable = T)

(barplot(GO_tg, label_format = "text", showCategory = 10 ) + 
    #scale_fill_gradient() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("GO_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/go_doublets_NK_Mono_CD56dim.pdf", w = 140, h = 75)

KEGG_tg <- enrichKEGG(de_markers$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10,
                      maxGSSize = 500, qvalueCutoff = 0.05, use_internal_data = F)

barplot(KEGG_tg, showCategory = 5, title = "KEGG Pathway", x = "GeneRatio")

(barplot(KEGG_tg, label_format = "text", showCategory = 10 ) + 
    #scale_fill_gradient(high = "slateblue4", low = "coral") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("KEGG_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/kegg_doublets_NK_Mono_CD56dim.pdf", w = 140, h = 75)

#CD4_Vd2_Mono
de_markers <- FindMarkers(doublets, ident.1 = "CD4_Vd2_Mono", group.by = "anno_4th", verbose = T)

de_markers %>% write.csv(file = "/home/big/zheng_song/aml/comb_s15_rerun/doublets/CD4_Vd2_Mono_DEGs.csv", row.names = T)

de_markers$gene <- rownames(de_markers)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05&avg_log2FC > 0.25, gene, "")), colour = "red", size = 3)

de_markers %>% arrange(-avg_log2FC)

# GSEA on de_markers from unknown clusters
# non-transplant vs transplant
# GSEA
library(clusterProfiler)
library(enrichplot)
library(dbplyr)
library("AnnotationHub")
#BiocManager::install("KEGGREST")

# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

organism = 'org.Hs.eg.db'
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

de_markers <- de_markers %>% dplyr::filter(avg_log2FC > 0.25) %>%  top_n(n = 500, wt = avg_log2FC)

#all_cells_marker_2 <- subset(all_cells_marker_2, rowSums(all_cells_marker_2[5] < 0.05) > 0)

de_markers$gene <- rownames(de_markers)

de_markers <-  bitr(de_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_tg <- enrichGO(de_markers$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                  maxGSSize = 500, readable = T)

(barplot(GO_tg, label_format = "text", showCategory = 10 ) + 
    #scale_fill_gradient() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("GO_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/go_doublets_CD4_Vd2_Mono.pdf", w = 140, h = 75)

KEGG_tg <- enrichKEGG(de_markers$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10,
                      maxGSSize = 500, qvalueCutoff = 0.05, use_internal_data = F)

barplot(KEGG_tg, showCategory = 5, title = "KEGG Pathway", x = "GeneRatio")

(barplot(KEGG_tg, label_format = "text", showCategory = 10, color = ) + 
    #scale_fill_gradient(high = "slateblue4", low = "coral") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("KEGG_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/kegg_doublets_CD4_Vd2_Mono.pdf", w = 140, h = 75)


#transplant vs non
de_markers <- FindMarkers(doublets, ident.1 = "CD4_CD8_Mono", group.by = "anno_4th", verbose = T)

de_markers %>% write.csv(file = "/home/big/zheng_song/aml/comb_s15_rerun/doublets/CD4_CD8_Mono_DEGs.csv", row.names = T)

de_markers$gene <- rownames(de_markers)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05&avg_log2FC > 0.25, gene, "")), colour = "red", size = 3)

de_markers %>% arrange(-avg_log2FC)

# GSEA on de_markers from unknown clusters
# non-transplant vs transplant
# GSEA
library(clusterProfiler)
library(enrichplot)
library(dbplyr)
library("AnnotationHub")
#BiocManager::install("KEGGREST")

# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

organism = 'org.Hs.eg.db'
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

de_markers <- de_markers %>% dplyr::filter(avg_log2FC > 0.25) %>%  top_n(n = 500, wt = avg_log2FC)

#all_cells_marker_2 <- subset(all_cells_marker_2, rowSums(all_cells_marker_2[5] < 0.05) > 0)

de_markers$gene <- rownames(de_markers)

de_markers <-  bitr(de_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_tg <- enrichGO(de_markers$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                  maxGSSize = 500, readable = T)

(barplot(GO_tg, label_format = "text", showCategory = 10 ) + 
    #scale_fill_gradient() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
    dotplot_theme + 
    ggtitle("GO_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/go_doublets_CD4_CD8_Mono.pdf", w = 160, h = 75)

KEGG_tg <- enrichKEGG(de_markers$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10,
                      maxGSSize = 500, qvalueCutoff = 0.05, use_internal_data = F)

barplot(KEGG_tg, showCategory = 5, title = "KEGG Pathway", x = "GeneRatio")

(barplot(KEGG_tg, label_format = "text", showCategory = 10, color = ) + 
    #scale_fill_gradient(high = "slateblue4", low = "coral") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("KEGG_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/kegg_doublets_CD4_CD8_Mono.pdf", w = 140, h = 75)

#CD4_CD8_Mono
de_markers <- FindMarkers(doublets, ident.1 = "CD4_CD8_Mono", group.by = "anno_4th", verbose = T)

de_markers %>% write.csv(file = "/home/big/zheng_song/aml/comb_s15_rerun/doublets/CD4_CD8_Mono_DEGs.csv", row.names = T)

de_markers$gene <- rownames(de_markers)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05&avg_log2FC > 0.25, gene, "")), colour = "red", size = 3)

de_markers %>% arrange(-avg_log2FC)

# GSEA on de_markers from unknown clusters
# non-transplant vs transplant
# GSEA
library(clusterProfiler)
library(enrichplot)
library(dbplyr)
library("AnnotationHub")
#BiocManager::install("KEGGREST")

# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

organism = 'org.Hs.eg.db'
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

de_markers <- de_markers %>% dplyr::filter(avg_log2FC > 0.25) %>%  top_n(n = 500, wt = avg_log2FC)

#all_cells_marker_2 <- subset(all_cells_marker_2, rowSums(all_cells_marker_2[5] < 0.05) > 0)

de_markers$gene <- rownames(de_markers)

de_markers <-  bitr(de_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_tg <- enrichGO(de_markers$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, minGSSize = 10,
                  maxGSSize = 500, readable = T)

(barplot(GO_tg, label_format = "text", showCategory = 10 ) + 
    #scale_fill_gradient() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
    dotplot_theme + 
    ggtitle("GO_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/go_doublets_CD4_CD8_Mono.pdf", w = 160, h = 75)

KEGG_tg <- enrichKEGG(de_markers$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10,
                      maxGSSize = 500, qvalueCutoff = 0.05, use_internal_data = F)

barplot(KEGG_tg, showCategory = 5, title = "KEGG Pathway", x = "GeneRatio")

(barplot(KEGG_tg, label_format = "text", showCategory = 10, color = ) + 
    #scale_fill_gradient(high = "slateblue4", low = "coral") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    dotplot_theme + 
    ggtitle("KEGG_enrichment_analysis") +
    theme(axis.text.x = element_text(size = 6, family = "Arial", angle = 0, margin = margin(0,0,0,0)),
          axis.title.x = element_text(size = 6, family = "Arial"),
          plot.title = element_text(size = 6, family = "Arial", margin = margin(1.5,1.5,3,1.5,"pt")),
          legend.title = element_text(angle = 0),
          panel.border = element_blank())) %>% figsave("comb_s15_rerun/doublets/kegg_doublets_CD4_CD8_Mono.pdf", w = 140, h = 75)



