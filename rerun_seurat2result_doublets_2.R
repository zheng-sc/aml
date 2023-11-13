#characterize unknown/doublets cluster include patient clinical info
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

# read rds ----------------------------------------------------------------
doublets <- readRDS(paste0(wd_path, "/comb_s15_rerun/aml_doublets_reso2.Rds"))


