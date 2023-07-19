# If you have canonical sites, annotate the canonical info

rm(list=ls())
options(stringsAsFactors = F)

library(dplyr)

##---------------- input -----------------------

load("tmp2_aa_change_annotated.Rdata")
canonical <- pull(read.table("~/projects/DSN-seq/data/ref/Osa/Osa_canonical_sites.txt", sep = "\t"),V1)


##---------------- data process -----------------------

all_sites_annotation$canonical <- all_sites_annotation$edit_id %in% canonical

group_by(all_sites_annotation, canonical) %>% 
  summarise(count = n())


##---------------- output -----------------------

save(all_sites_annotation,
     file = "tmp3_canonical_sites_annotated.Rdata")
