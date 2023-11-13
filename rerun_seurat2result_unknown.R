#characterize unknown cluster 
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

# focus on Unknown compartment --------------------------------------------
unknown <- subset(all_cells, anno_1st == "Unknown")

Feature_rast5(unknown, g = c("FCGR3A", "CD247", "IL1B", "NCAM1"), facets = "tp", colorset = "gg", assay = "RNA", color_grd = "grd")

Feature_rast5(all_cells, g = c("CD14", "FCGR3A", "CD247", "IL1B", "NCAM1", "PDCD1", "ITGAX", "TYROBP"), colorset = "gg", assay = "RNA", 
              color_grd = "grd", ncol = 4) %>% 
  figsave("comb_s15_rerun/umap_unknown_test.pdf", w = 400, h = 200)

all_cells_marker_2 %>% arrange(-avg_log2FC)

#options(Seurat.object.assay.version = "v3")

library(Seurat)

library(patchwork)

#unknown <- SplitObject(unknown, split.by = "stim")
as.sparse(unknown@assays$RNA$counts)

new_unknown <- Seurat::CreateSeuratObject(as.sparse(unknown@assays$RNA$counts))

new_unknown <- SCTransform(new_unknown, assay = "RNA")

new_unknown <- RunPCA(new_unknown) %>% RunUMAP(dims = 1:14)

new_unknown <- FindNeighbors(new_unknown, dims = 1:14, assay = "SCT")

new_unknown <- FindClusters(new_unknown, resolution = 1)

Feature_rast(new_unknown)

## pK Identification (no ground-truth) 
sweep.res.list_unknown <- DoubletFinder::paramSweep_v3(new_unknown, PCs = 1:14, sct = T, num.cores = 32)
sweep.stats_unknown <- DoubletFinder::summarizeSweep(sweep.res.list_unknown, GT = FALSE)
bcmvn_unknown <- DoubletFinder::find.pK(sweep.stats_unknown)

## pK Identification (ground-truth)
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE, assay = "SCT")
# 
# gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- new_unknown$seurat_clusters
homotypic.prop <- DoubletFinder::modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(new_unknown@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
new_unknown <- DoubletFinder::doubletFinder_v3(new_unknown, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)

#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)


# join doublets column into unknown object --------------------------------

left_join(new_unknown@meta.data, select(unknown@meta.data, ))

dbl <- new_unknown@meta.data %>% dplyr::select("bc_wells", "DF.classifications_0.25_0.09_1328")

colnames(dbl) <- c("bc_wells", "doublets")

unknown@meta.data <- left_join(unknown@meta.data, dbl, by = "bc_wells")

rownames(unknown@meta.data) <- unknown@meta.data$bc_wells

unknown <- RunAzimuth(query = unknown, reference = "bonemarrowref", verbose = T, assay = "RNA")

unknown@meta.data %<>% dplyr::mutate(transplant_group = case_when(tp == "Transplant" ~ "Transplant",
                                                                  tp == "D30" ~ "Non-transplant",
                                                                  tp == "D100" ~ "Non-transplant"))

Idents(unknown) <- unknown$transplant_group

Feature_rast(unknown, "predicted.celltype.l1", colorset = "gg", facets = "transplant")

#aggregation-based (pseudobulk) workflow
bulk <- AggregateExpression(unknown, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("transplant_group", "donor"))
bulk$transplant <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$donor <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
#bulk$tp <- sapply(strsplit(bulk$donor, split = "-"), "[", 2)
#bulk$donor <- sapply(strsplit(bulk$donor, split = "-"), "[", 1)

de_markers <- FindMarkers(bulk, ident.1 = "Non-transplant", ident.2 = "Transplant", group.by = "transplant", slot = "counts", test.use = "DESeq2",
                          verbose = T)

de_markers %>% write.csv(file = "/home/big/zheng_song/aml/comb_s15_rerun/unknown_non_trans_against_trans.csv", row.names = T)

de_markers$gene <- rownames(de_markers)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) +
  ylab("-log10(unadjusted p-value)") + 
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05&avg_log2FC > 1, gene, "")), colour = "red", size = 3)

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

de_markers <- de_markers %>% dplyr::filter(avg_log2FC > 1) %>%  top_n(n = 500, wt = avg_log2FC)

#all_cells_marker_2 <- subset(all_cells_marker_2, rowSums(all_cells_marker_2[5] < 0.05) > 0)

de_markers$gene <- rownames(de_markers)

de_markers <-  bitr(de_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_tg <- enrichGO(de_markers$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, minGSSize = 10,
                  maxGSSize = 500, readable = T)

barplot(GO_tg, showCategory = 100) #+ facet_grid(ONTOLOGY~., scales = "free") + 
#  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 100))

KEGG_tg <- enrichKEGG(de_markers$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10,
                      maxGSSize = 500, qvalueCutoff = 0.05, use_internal_data = F)

barplot(KEGG_tg, showCategory = 30, title = "KEGG Pathway", x = "GeneRatio")


saveRDS(
  object = unknown,
  file = "aml_unknown_neutrophil.Rds",
  destdir = paste0(wd_path, '/comb_s15_rerun/aml_unknown_neutrophil')
)

unknown <- readRDS("/home/big/zheng_song/aml/comb_s15_rerun/aml_unknown_neutrophil/aml_unknown_neutrophil.Rds")

unknown <- FindNeighbors(unknown, reduction = "harmony.full", dims = 1:14, assay = "RNA")
unknown <- FindClusters(unknown, resolution = 0.5, cluster.name = 'harmony_clusters')
unknown <- RunUMAP(unknown, reduction = "harmony.full", dims = 1:14, reduction.name = "umap.full", reduction.key = "UMAPfull_")

DefaultAssay(unknown) <- "RNA"

unknown$tp <- factor(unknown$tp, levels = c("Before_SCT", "Transplant", "D30", "D100"))


ViolinPlot(unknown, c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", 
                      "CD14", "NCAM1", "S100A8", "S100A9", "ELANE",
                      "SLAMF7", "CD163", "SIGLEC1", "IL7R", "IL1R2"), sz = 0.3, group.by = "tp") %>% 
  figsave("comb_s15_rerun/unknown_lin_myeloid_markers.pdf", w = 500, h = 300)


RidgePlot(unknown, c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", 
                     "CD14", "NCAM1", "S100A8", "S100A9", "ELANE",
                     "SLAMF7", "CD163", "SIGLEC1", "IL7R", "TBX21",
                     "EOMES", "RORC", "IL6", "IL1R2", "IL1RL1"), group.by = "tp") %>% 
  figsave("comb_s15_rerun/unknown_lin_myeloid_markers_ridge.pdf", w = 500, h = 300)

Feature_rast5(unknown, colorset = "gg", c("S100A8", "S100A9", "S100A12", "CEACAM8", "CD24", 
                                          "LTF", "CAMP", "FCGR3B", "CXCR2", "CMTM2", 
                                          "CD3E", "CD3D", "CD247", "TRAC", "TRDC", "NCAM1", 
                                          "TRGC1", "CD4", "CD8A", "CD8B", "IL7R"), 
              color_grd = "grd", ncol = 5) %>% 
  figsave("comb_s15_rerun/umap_neutrohil_test.pdf", w = 500, h = 400)


# unknown cell pct of total cell number
unknown_df <- unknown@meta.data %>% group_by(donor) %>% dplyr::select(sample, tp) %>% dplyr::count(sample)

sample_count <- all_cells@meta.data %>% group_by(donor) %>% dplyr::select(sample) %>% dplyr::count(sample)

unknown_pct <- left_join(unknown_df, sample_count, by = "sample")

unknown_pct$pct <- unknown_pct$n.x/unknown_pct$n.y

unknown_pct$tp <- str_extract(unknown_pct$sample, "Transplant|D30|D100|Before_SCT")

unknown_pct$tp <- factor(unknown_pct$tp, levels = c("Before_SCT", "Transplant", "D30", "D100"))

(ggplot(data = unknown_pct, aes(x = donor.x, y = pct, fill = tp))+
    geom_bar(stat="identity",position='dodge') + 
    scale_fill_manual(values = c("#800000",  "#b5a699", "#67755c","#058789")) + violin_theme +
    ylab("Cells") + xlab("Donor") + theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1)) + labs(fill = "Timepoint")) %>% 
  figsave("comb_s15_rerun/unknown_pct.pdf", w = 200, h = 100)

ggplot(data = unknown_pct, aes(x = tp, y = pct, stratum = donor.x)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
  scale_fill_manual(values = c("#800000", "#b5a699", "#67755c", "#058789")) +
  ylab("Cells") +
  xlab("Timepoint") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(fill = "Donor")

library(viridis)

ggplot(unknown_pct, aes(x = tp, y = pct)) + 
  geom_point(aes(color = factor(donor.x))) +
  geom_line(aes(group = donor.x, color = factor(donor.x))) +
  scale_x_discrete(limits = c("Before_SCT", "Transplant", "D30", "D100"))

# check ribo gene

ribo_genes <- grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(unknown@assays$RNA), value = T)

unknown[["percent_ribo"]] <- PercentageFeatureSet(unknown, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", assay = "RNA")

RidgePlot(unknown, c("percent_ribo"), group.by = "tp") %>% 
  figsave("comb_s15_rerun/unknown_pct_ribo_ridge.pdf", w = 300, h = 100)










Feature_rast5(all_cells, colorset = "gg", c("S100A8", "S100A9", "S100A12", "CEACAM8", "CD24", 
                                            "LTF", "CAMP", "FCGR3B", "CXCR2", "CMTM2", 
                                            "CD3E", "CD3D", "CD247", "TRBC1", "TRDC", "NCAM1",
                                            "FCGR3A", "MME", "FOXP3", "IL2RA"), 
              color_grd = "grd", ncol = 5) %>% 
  figsave("comb_s15_rerun/umap_all_cells_test.pdf", w = 500, h = 400)

Feature_rast5(all_cells, colorset = "gg", c("CD3E", "CD3D", "CD247", "TRAC", "TRBC1", 
                                            "TRBC2", "TRDC", "TRGC1", "TRGC2", "ITGAM", 
                                            "FUT4", "CEACAM8", "CD4", "CD8A", "CD8B", 
                                            "TBX21", "EOMES"), 
              color_grd = "grd", ncol = 5) %>% 
  figsave("comb_s15_rerun/umap_all_cells_tcr_neu.pdf", w = 500, h = 400)

ViolinPlot(unknown, c("S100A8", "S100A9", "ELANE", "SLAMF7", "CD163", "SIGLEC1"), sz = 0, )









# cd127 markers
ViolinPlot(unknown, c("IL7R"), sz = 0.5) %>% figsave("comb_s15_rerun/unknown_IL7R.pdf", w = 100, h = 100)
ViolinPlot(all_cells, c("ITGAM", "ITGAX"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/cd11b_cd11c.pdf", w = 200, h = 100)

ViolinPlot(all_cells, c("IL7R"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/IL7R.pdf", w = 100, h = 100)



ViolinPlot(all_cells, c("PTGDR2", "KIT"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/ckit.pdf", w = 200, h = 100)


ViolinPlot(all_cells, c("GNLY", "PRF1", "GZMB", "IFNG"), sz = 0.1)

Feature_rast(all_cells, d1 = "lin1", d2 = "IL7R", colorset = "gg", do.label = F)

Feature_rast(unknown, d1 = "lin1", d2 = "IL7R", colorset = "gg", do.label = F)

Feature_rast5(all_cells, c("cell_type", "ITGAM", "ITGAX", "CD69"), colorset = "gg", color_grd = "grd")

Feature_rast(unknown, d1 = "nFeature_RNA", d2 = "nCount_RNA", g = "transplant", colorset = "gg", do.label = F, noaxis = F)








# neutrophil gene list from https://doi.org/10.1016/j.cell.2020.08.001 
# cohort1 and cohort2 cluster 5 named neutrohpil
library("readxl")

ss_co_1 <- read_xlsx("/home/big/zheng_song/aml/schulte_schrepping_2020_cell/cohort_1.xlsx", col_names = T)

ss_co_1$gene...8 <- NULL

colnames(ss_co_1)[1] <- "gene"

ss_co_1 <- ss_co_1[-1,]

ss_co_1 %>% mutate(cell_type = case_when(cluster == "0" ~ "classical_mono",
                                         cluster == "1" ~ "cd83high_mono",
                                         cluster == "2" ~ "cd163high_mono",
                                         cluster == "3" ~ "s100ahigh_mono",
                                         cluster == "4" ~ "non_classical_mono",
                                         cluster == "5" ~ "neutrophils",
                                         cluster == "6" ~ "immature_neutrophils",
                                         cluster == "7" ~ "mDCs",
                                         cluster == "8" ~ "pDCs",
                                         cluster %in%  c("9", "10", "11") ~ "CD4_t",
                                         cluster %in%  c("12", "13", "14") ~ "CD8_t",
                                         cluster == "15" ~ "nk",
                                         cluster %in%  c("16", "17", "18") ~ "b",
                                         cluster ==  "19" ~ "Plasmablasts",
                                         cluster ==  "20" ~ "megakaryocyte",
                                         cluster ==  "21" ~ "mixed",
                                         cluster ==  "22" ~ "unidenfined")) -> ss_co_2

gc_name <- ss_co_2$cell_type %>% unique() %>% as.vector()

gene_cluster <- map(gc_name, function(x) ss_co_2 %>% filter(cell_type == x) %>% dplyr::select(gene) %>% pull()) %>% setNames(gc_name)

for (i in gc_name) {
  all_cells %<>%  AddModuleScore(features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
  colnames(all_cells@meta.data) %<>% str_replace("(?<=\\w)1", '')
}


ViolinPlot(all_cells, c("classical_mono", "s00ahigh_mono", "cd83high_mono", "cd63high_mono", 
                        "non_classical_mono", "immature_neutrophils", "neutrophils", "megakaryocyte"), sz = 0, ncol = 4)

ViolinPlot(all_cells, c("CD4_t", "CD8_t", "nk", 
                        "b", "mDCs", "pDCs"), sz = 0, ncol = 4)

ViolinPlot(all_cells, c("classical_mono", "s00ahigh_mono", "cd83high_mono", "cd63high_mono", 
                        "non_classical_mono", "immature_neutrophils", "neutrophils", "megakaryocyte", "nk"), 
           sz = 0, ncol = 3, facet = "transplant") %>% figsave("comb_s15_rerun/umap_all_cells_cell_gms.pdf", w = 300, h = 300)

neu_1 <- ss_co_1 %>% filter(cluster == "5")

neu_1 %<>% arrange(-avg_logFC) %>% top_n(n = 20, wt = avg_logFC)

all_cells <- AddModuleScore(all_cells, features = list(neu_1$gene), assay = "RNA", name = "cohort1_neu")

Feature_rast5(all_cells, g = "cohort1_neu1", color_grd = "grd")

ViolinPlot(all_cells, c("cohort1_neu1"))


gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_c_list_ent

gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)

gene_cluster$GM_D

for (i in gc_name) {
  all_cells %<>%  AddModuleScore(features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
  colnames(all_cells@meta.data) %<>% str_replace("(?<=\\w)1", '')
}



ViolinPlot(all_cells, "GATA3", sz = 0, group.by = "anno_st")




# Neutrophil-myeloid progenitor gene list from https://doi.org/10.1038/s41586-019-1652-y
pop_2019 <- read_xlsx("/home/big/zheng_song/aml/popescu_2019_nature/degs.xlsx", col_names = T)

pop_2019$`Neutrophil-myeloid progenitor`

all_cells <- AddModuleScore(all_cells, features = list(pop_2019$`Neutrophil-myeloid progenitor`), assay = "RNA", name = "pop_neu_prog")

ViolinPlot(all_cells, c("pop_neu_prog"))


Feature_rast5(all_cells, "anno_1st", colorset = "um", facets = "transplant", sz = 0.02)

Feature_rast5(unknown, colorset = "um", facets = "transplant", sz = 0.8)

# marker genes from DOI: 10.1126/sciimmunol.abn6429
# 
huo_2023 <- read_xlsx("/home/big/zheng_song/aml/huo_2023_science_immunology/marker_genes.xlsx", col_names = T, sheet = "Marker genes of TNCs")

colnames(huo_2023) <- huo_2023[1,]

huo_2023 <- huo_2023[-1,]

filter(huo_2023, TNCs == "MatureNeu1")$gene

all_cells <- AddModuleScore(all_cells, features = list(filter(huo_2023, TNCs == "MatureNeu1")$gene), assay = "RNA", name = "MatureNeu1_")

ViolinPlot(all_cells, c("MatureNeu1_1"))

ViolinPlot(all_cells, c("RNASE2", "CXCL8", "PLBD1", "VCAN", "MAFB", "CD14"), sz = 0)

ViolinPlot(all_cells, c("HLA-DRA", "HLA-DRB1", "S100A8", "S100A9", "S100A12"), sz = 0)

ViolinPlot(all_cells, c("CD68", "ADGRE1"), sz = 0)

# check ribo gene

ribo_genes <- grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(all_cells@assays$RNA), value = T)

all_cells[["percent_ribo"]] <- PercentageFeatureSet(all_cells, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", assay = "RNA")

Feature_rast5(all_cells, "percent_ribo", color_grd = "grd")

ViolinPlot(all_cells, "percent_mt")

ViolinPlot(all_cells, "lin1")

#check cell cycle gene
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all_cells <- CellCycleScoring(all_cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

ViolinPlot(all_cells, "G2M.Score", group.by = "anno_1st")



# annotation
unknown_markers <- FindAllMarkers(unknown, assay = "RNA", verbose = T, slot = "counts")

lin <- c("CD2", "CD3E", "CD14", "FCGR3A", "CD19", "MS4A1", "NCAM1", "GYPA")

all_cells <- AddModuleScore(all_cells, features = list(lin), assay = "RNA", name = "lin")

unknown <- AddModuleScore(unknown, features = list(lin), assay = "RNA", name = "lin")

Feature_rast5(all_cells, "lin1", color_grd = "grd")

ViolinPlot(all_cells, "lin1", sz = 0, group.by = "anno_1st") %T>% figsave("comb_s15_rerun/lin.pdf", w = 100, h = 100)

ViolinPlot(all_cells, lin, sz = 0, group.by = "anno_1st", ncol = 4) %T>% figsave("comb_s15_rerun/lin_sep.pdf", w = 400, h = 100)


cln <- 6
# myeloid cell markers

ViolinPlot(unknown, c("S100A8", "S100A9", "ELANE", "SLAMF7", "CD163", "SIGLEC1"), sz = 0, )

# cd127 markers
ViolinPlot(all_cells, c("IL7R"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/IL7R.pdf", w = 100, h = 100)
ViolinPlot(all_cells, c("ITGAM", "ITGAX"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/cd11b_cd11c.pdf", w = 200, h = 100)

ViolinPlot(all_cells, c("IL7R"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/IL7R.pdf", w = 100, h = 100)



ViolinPlot(all_cells, c("PTGDR2", "KIT"), sz = 0, group.by = "anno_1st") %>% figsave("comb_s15_rerun/ckit.pdf", w = 200, h = 100)


ViolinPlot(all_cells, c("GNLY", "PRF1", "GZMB", "IFNG"), sz = 0.1)

Feature_rast(all_cells, d1 = "lin1", d2 = "IL7R", colorset = "gg", do.label = F)

Feature_rast(unknown, d1 = "lin1", d2 = "IL7R", colorset = "gg", do.label = F)

Feature_rast5(all_cells, c("cell_type", "ITGAM", "ITGAX", "CD69"), colorset = "gg", color_grd = "grd")

Feature_rast(unknown, d1 = "nFeature_RNA", d2 = "nCount_RNA", g = "transplant", colorset = "gg", do.label = F, noaxis = F)


