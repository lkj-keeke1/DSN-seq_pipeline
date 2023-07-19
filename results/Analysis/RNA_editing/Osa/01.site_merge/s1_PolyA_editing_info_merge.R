rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

# input editing tables
Osa_PolyA_rep1 <- read.table("../../../../P1_call_editing/02.call_editing/Osa_PolyA_rep1.filtered.table", header = T, sep = "\t")
Osa_PolyA_rep2 <- read.table("../../../../P1_call_editing/02.call_editing/Osa_PolyA_rep2.filtered.table", header = T, sep = "\t")
Osa_PolyA_rep3 <- read.table("../../../../P1_call_editing/02.call_editing/Osa_PolyA_rep3.filtered.table", header = T, sep = "\t")
Osa_PolyA_merge <- read.table("../../../../P1_call_editing/02.call_editing/Osa_PolyA_merge.filtered.table", header = T, sep = "\t")

PolyA <- list(Osa_PolyA_rep1 = Osa_PolyA_rep1$edit_id,
            Osa_PolyA_rep2 = Osa_PolyA_rep2$edit_id,
            Osa_PolyA_rep3 = Osa_PolyA_rep3$edit_id,
            Osa_PolyA_merge = Osa_PolyA_merge$edit_id)


# merge editing sites
PolyA_sites <- c()
for (i in PolyA){
  PolyA_sites <- c(PolyA_sites, i)
}
PolyA_sites <- distinct(tibble(PolyA_sites))


# distribution of editing sites
PolyA_sites$Osa_PolyA_rep1 <- PolyA_sites$PolyA_sites %in% PolyA$Osa_PolyA_rep1
PolyA_sites$Osa_PolyA_rep2 <- PolyA_sites$PolyA_sites %in% PolyA$Osa_PolyA_rep2
PolyA_sites$Osa_PolyA_rep3 <- PolyA_sites$PolyA_sites %in% PolyA$Osa_PolyA_rep3
PolyA_sites$Osa_PolyA_merge <- PolyA_sites$PolyA_sites %in% PolyA$Osa_PolyA_merge


PolyA_sites <- column_to_rownames(PolyA_sites, var = "PolyA_sites") 
PolyA_sites$PolyA_num <- rowSums(PolyA_sites) - PolyA_sites$Osa_PolyA_merge
PolyA_sites <- rownames_to_column(PolyA_sites, var = "edit_id")

save(PolyA_sites,
     file = "PolyA_editing_sites.Rdata")