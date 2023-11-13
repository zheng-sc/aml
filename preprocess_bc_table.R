library(dplyr)
library(stringr)

library(magrittr)

bcs <- read.csv('~/big/zheng/parse_bio/metastasis_liver_blood/analysis_novoseq/seurat_analysis/bc_table.csv', sep = '\t')
bcs[97:192,] <- bcs[1:96,]

bcs$bc2_well %<>% str_pad(width = 2, pad = '0', side = 'left')

bcs$bc3_well <- bcs$bc2_well

bcs <- bcs[,c(3:6)]

bcs %<>% as_data_frame()

bcs_2 <- read.csv('/home/big/zheng/parse_bio/aml/novaseq/analysis/comb_s1_s2/process/barcode_data.csv', header = 5, sep = "\t")

bcs_2 <- bcs_2[5:197,] %>% as_data_frame()


bcs_2 <- str_split_fixed(bcs_2$value, pattern = ',', n = 7)

bcs_2 %>% tail()

colnames(bcs_2) <- bcs_2[1,]

bcs_2 <- bcs_2[-1,]
bcs_2 %<>% as_data_frame() #%>% select(c('bci', 'sequence'))

bcs_2 <- bcs_2[,c(2,6)]

bcs_2 %<>% relocate(well_int, .after = sequence)

colnames(bcs_2) <- c('bc1', 'bc1_well')

bcs_2 %<>% as_data_frame()

bcs_2$bc1_well %<>% str_pad(width = 2, pad = '0', side = 'left')

mega_bc <- cbind(bcs_2, bcs)

write.csv(mega_bc, file = 'mega_barcode.csv', quote = F, row.names = F)
