# annotate editing extent and sequenceing depth information (only candidate RNA editing sites)

rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

##---------------- input -----------------------

table_list <- "sample.call.list"
table_list <- read.table(table_list)

#load(file = "tmp2_aa_change_annotated.Rdata")  # if no canonical information annotated
load(file = "tmp3_canonical_sites_annotated.Rdata")


##---------------- data process -----------------------

# annotate editing extent and sequenceing depth information
for (i in 1:nrow(table_list)) {
  sample <- basename(table_list$V1[i])
  tmp_table <- paste(table_list$V1[i], ".filtered.table", sep = "")
  
  tmp_table <- read.table(tmp_table, sep = "\t",header = T) %>% 
    select(edit_id, freq, depth)
  colnames(tmp_table) <- c("edit_id", paste(sample, "_extent", sep = ""), paste(sample, "_depth", sep = ""))

  all_sites_annotation <- left_join(all_sites_annotation, tmp_table, by  = "edit_id")
}


##---------------- output -----------------------

save(all_sites_annotation,
     file = "annotated_sites_filtered.Rdata")
write.table(all_sites_annotation, file = "annotated_sites_filtered.table",
            col.names = T, row.names = F,
            quote = F, sep = "\t")

