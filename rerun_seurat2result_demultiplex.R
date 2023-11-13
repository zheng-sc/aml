# incorporate souporcell results
#Sketch integration rerun tcr analysis
#Seurat V5 
#zheng
#Tue Sep 13 14:02:50 2023
#15 sub-libs 777,530 cells

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

# SNV on sublib1(souporcell) ----------------------------------------------
all_cells$sublib <- str_extract(all_cells$bc_wells, "s\\d{1,2}")

sublib_1 <- subset(all_cells, sublib == "s1")

soup_1 <- read.table("/home/big/zheng/souporcell/aml_souporcell/aml_s1/aml_s1_clusters.tsv", sep = "\t", header = T)

soup_1_clean <- soup_1 %>% filter(barcode %in% grep("s1$", soup_1$barcode, value = T)) %>% dplyr::select(c("barcode", "status", "assignment")) %>% 
  filter(status == "singlet")

