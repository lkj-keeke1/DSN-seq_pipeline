rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

# input editing tables
Osa_DSN_rep1 <- read.table("../../../../P1_call_editing/02.call_editing/Osa_DSN_rep1.filtered.table", header = T, sep = "\t")
Osa_DSN_rep2 <- read.table("../../../../P1_call_editing/02.call_editing/Osa_DSN_rep2.filtered.table", header = T, sep = "\t")
Osa_DSN_rep3 <- read.table("../../../../P1_call_editing/02.call_editing/Osa_DSN_rep3.filtered.table", header = T, sep = "\t")
Osa_DSN_merge <- read.table("../../../../P1_call_editing/02.call_editing/Osa_DSN_merge.filtered.table", header = T, sep = "\t")

DSN <- list(Osa_DSN_rep1 = Osa_DSN_rep1$edit_id,
            Osa_DSN_rep2 = Osa_DSN_rep2$edit_id,
            Osa_DSN_rep3 = Osa_DSN_rep3$edit_id,
            Osa_DSN_merge = Osa_DSN_merge$edit_id)


# merge editing sites
DSN_sites <- c()
for (i in DSN){
  DSN_sites <- c(DSN_sites, i)
}
DSN_sites <- distinct(tibble(DSN_sites))


# distribution of editing sites
DSN_sites$Osa_DSN_rep1 <- DSN_sites$DSN_sites %in% DSN$Osa_DSN_rep1
DSN_sites$Osa_DSN_rep2 <- DSN_sites$DSN_sites %in% DSN$Osa_DSN_rep2
DSN_sites$Osa_DSN_rep3 <- DSN_sites$DSN_sites %in% DSN$Osa_DSN_rep3
DSN_sites$Osa_DSN_merge <- DSN_sites$DSN_sites %in% DSN$Osa_DSN_merge


DSN_sites <- column_to_rownames(DSN_sites, var = "DSN_sites") 
DSN_sites$DSN_num <- rowSums(DSN_sites) - DSN_sites$Osa_DSN_merge
DSN_sites <- rownames_to_column(DSN_sites, var = "edit_id")

save(DSN_sites,
     file = "DSN_editing_sites.Rdata")