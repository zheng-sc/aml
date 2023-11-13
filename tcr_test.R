tcr_raw <- c(1:15) %>% map(~read.csv(paste0('~/big/zheng/parse_bio/aml/novaseq/analysis/trust4_rerun_2/S',.x,'_barcode_report.tsv'),sep = '\t') %>% 
                             distinct(X.barcode, .keep_all = T) %>% 
                             mutate(bc1 = str_sub(X.barcode, 69, 76 ), 
                                    bc2 = str_sub(X.barcode, 39, 46 ),
                                    bc3 = str_sub(X.barcode, 1, 8 )) %>%
                             left_join(select(mega_bcs, 1,2), by = 'bc1') %>%
                             left_join(select(mega_bcs, 3,4), by = 'bc2') %>% 
                             left_join(select(mega_bcs, 5,6), by = 'bc3')  %>% 
                             mutate(bc_wells = paste0(bc1_well, '_', bc2_well,'_', bc3_well,'__s', .x)) %>% 
                             filter(!grepl('NA', bc_wells))
)


tcr_raw[[3]]

tcr_raw[[3]]$raw_bc <- paste0(tcr_raw[[3]]$bc1_well, '_', tcr_raw[[3]]$bc2_well,'_', tcr_raw[[3]]$bc3_well)

all_cells$raw_bc <- paste0(str_pad(all_cells$bc1_wind , 2, pad = "0"), '_', str_pad(all_cells$bc2_wind , 2, pad = "0"),'_', str_pad(all_cells$bc3_wind , 2, pad = "0"))

all_cells$sublib <- str_extract(string = all_cells$bc_wells, pattern = "s\\d{1,2}$")

filter(all_cells@meta.data, all_cells$raw_bc %in% tcr_raw[[3]]$raw_bc) %>% select(sublib) %>% table()

s3 <- read.csv("~/big/zheng_song/aml/comb_s15/s3_barcode.csv", row.names = 1)

tcr_raw[[3]]$raw_bc %in% s3$X0 %>% table()
