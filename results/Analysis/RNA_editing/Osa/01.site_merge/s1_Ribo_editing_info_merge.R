rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

# input editing tables
`Osa_Ribo-off_rep1` <- read.table("../../../../P1_call_editing/02.call_editing/Osa_Ribo-off_rep1.filtered.table", header = T, sep = "\t")
`Osa_Ribo-off_rep2` <- read.table("../../../../P1_call_editing/02.call_editing/Osa_Ribo-off_rep2.filtered.table", header = T, sep = "\t")
`Osa_Ribo-off_rep3` <- read.table("../../../../P1_call_editing/02.call_editing/Osa_Ribo-off_rep3.filtered.table", header = T, sep = "\t")
`Osa_Ribo-off_merge` <- read.table("../../../../P1_call_editing/02.call_editing/Osa_Ribo-off_merge.filtered.table", header = T, sep = "\t")

Ribo <- list(`Osa_Ribo-off_rep1` = `Osa_Ribo-off_rep1`$edit_id,
            `Osa_Ribo-off_rep2` = `Osa_Ribo-off_rep2`$edit_id,
            `Osa_Ribo-off_rep3` = `Osa_Ribo-off_rep3`$edit_id,
            `Osa_Ribo-off_merge` = `Osa_Ribo-off_merge`$edit_id)


# merge editing sites
Ribo_sites <- c()
for (i in Ribo){
  Ribo_sites <- c(Ribo_sites, i)
}
Ribo_sites <- distinct(tibble(Ribo_sites))


# distribution of editing sites
Ribo_sites$`Osa_Ribo-off_rep1` <- Ribo_sites$Ribo_sites %in% Ribo$`Osa_Ribo-off_rep1`
Ribo_sites$`Osa_Ribo-off_rep2` <- Ribo_sites$Ribo_sites %in% Ribo$`Osa_Ribo-off_rep2`
Ribo_sites$`Osa_Ribo-off_rep3` <- Ribo_sites$Ribo_sites %in% Ribo$`Osa_Ribo-off_rep3`
Ribo_sites$`Osa_Ribo-off_merge` <- Ribo_sites$Ribo_sites %in% Ribo$`Osa_Ribo-off_merge`


Ribo_sites <- column_to_rownames(Ribo_sites, var = "Ribo_sites") 
Ribo_sites$`Ribo-off_num` <- rowSums(Ribo_sites) - Ribo_sites$`Osa_Ribo-off_merge`
Ribo_sites <- rownames_to_column(Ribo_sites, var = "edit_id")

save(Ribo_sites,
     file = "Ribo-off_editing_sites.Rdata")