rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

# input editing tables
Ath_DSN_rep1 <- read.table("../../../../P1_call_editing/02.call_editing/Ath_DSN_rep1.filtered.table", header = T, sep = "\t")
Ath_DSN_rep2 <- read.table("../../../../P1_call_editing/02.call_editing/Ath_DSN_rep2.filtered.table", header = T, sep = "\t")
Ath_DSN_rep3 <- read.table("../../../../P1_call_editing/02.call_editing/Ath_DSN_rep3.filtered.table", header = T, sep = "\t")
Ath_DSN_merge <- read.table("../../../../P1_call_editing/02.call_editing/Ath_DSN_merge.filtered.table", header = T, sep = "\t")
Ath_STS <- read.table("../../../../sup_STS/02.call_editing/Ath_STS.filtered.table", header = T, sep = "\t")

Ath <- list(Ath_DSN_rep1 = Ath_DSN_rep1$edit_id,
            Ath_DSN_rep2 = Ath_DSN_rep2$edit_id,
            Ath_DSN_rep3 = Ath_DSN_rep3$edit_id,
            Ath_DSN_merge = Ath_DSN_merge$edit_id,
            Ath_STS = Ath_STS$edit_id)


# merge editing sites
all_sites <- c()
for (i in Ath){
  all_sites <- c(all_sites, i)
}
all_sites <- distinct(tibble(all_sites))
colnames(all_sites) <- "edit_id"


# distribution of editing sites
all_sites$Ath_DSN_rep1 <- all_sites$edit_id %in% Ath$Ath_DSN_rep1
all_sites$Ath_DSN_rep2 <- all_sites$edit_id %in% Ath$Ath_DSN_rep2
all_sites$Ath_DSN_rep3 <- all_sites$edit_id %in% Ath$Ath_DSN_rep3
all_sites$Ath_DSN_merge <- all_sites$edit_id %in% Ath$Ath_DSN_merge
all_sites$Ath_STS <- all_sites$edit_id %in% Ath$Ath_STS


save(all_sites,
     file = "plot.Rdata")

all_sites <- mutate(all_sites,
                    DSN_num = Ath_DSN_rep1 + Ath_DSN_rep2 + Ath_DSN_rep3) %>% 
  select(edit_id,
         DSN_det = Ath_DSN_merge, DSN_num,
         STS_det = Ath_STS)

save(all_sites,
     file = "all_sites_merge.Rdata")