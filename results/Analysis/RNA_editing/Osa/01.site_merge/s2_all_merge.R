rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

# input editing sites
load("DSN_editing_sites.Rdata")
load("PolyA_editing_sites.Rdata")
load("Ribo-off_editing_sites.Rdata")
canonical_sites <- read.table("~/projects/DSN-seq/data/ref/Osa/Osa_canonical_sites.txt", sep = "\t", header = T) %>% 
  dplyr::rename(edit_id = gene_posi.edit)



all_sites <- full_join(DSN_sites,Ribo_sites, by = "edit_id") %>% 
  full_join(PolyA_sites) %>% 
  full_join(canonical_sites) %>% 
  select(edit_id, 
         DSN_det = Osa_DSN_merge, DSN_num,
         Ribo_det = `Osa_Ribo-off_merge`, `Ribo-off_num`,
         PolyA_det = Osa_PolyA_merge, PolyA_num) %>% 
    # NA -> FALSE/0
  mutate(DSN_det = if_else(is.na(DSN_det), FALSE, DSN_det),
         Ribo_det = if_else(is.na(Ribo_det), FALSE, Ribo_det),
         PolyA_det = if_else(is.na(PolyA_det), FALSE, PolyA_det),
         
         DSN_num = if_else(is.na(DSN_num), 0, DSN_num),
         `Ribo-off_num` = if_else(is.na(`Ribo-off_num`), 0, `Ribo-off_num`),
         PolyA_num = if_else(is.na(PolyA_num), 0, PolyA_num)
         ) 
  
  

save(all_sites,
     file = "all_sites_merge.Rdata")
  
