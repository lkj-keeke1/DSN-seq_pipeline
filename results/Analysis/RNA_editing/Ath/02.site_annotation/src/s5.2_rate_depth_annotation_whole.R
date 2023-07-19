# annotate editing extent and sequenceing depth information (all the sites)

rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)

##---------------- input -----------------------

table_list <- "sample.call.list"
table_list <- read.table(table_list)

load(file = "tmp2_aa_change_annotated.Rdata")  # if no canonical information annotated
#load(file = "tmp3_canonical_sites_annotated.Rdata")


##---------------- data process -----------------------

# posi in depth file = posi in all_sites_annotation + 100
all_sites_annotation <- mutate(all_sites_annotation,
                               seq_id = paste(seq, posi + 100, sep = " ")) # seq_id is used to annotate depth information

# annotate editing extent and sequenceing depth information
for (i in 1:nrow(table_list)) {
  sample <- basename(table_list$V1[i])                                 # sample name 
  tmp_filter <- paste(table_list$V1[i], ".filtered.table", sep = "")   # filtered table
  tmp_out <- paste(table_list$V1[i], ".out.table", sep = "")           # filtered out table
  tmp_depth <- paste(table_list$V1[i], ".depth", sep = "")             # depth table
  
  tmp_freq <- read.table(tmp_filter, sep = "\t",header = T) %>% 
    select(edit_id, freq) %>% 
    full_join(read.table(tmp_out, sep = "\t",header = T) %>% 
                select(edit_id, freq))
  tmp_depth <- read.table(tmp_depth, sep = "\t",header = F) %>% 
    select(seq_id = V1, depth = V2)
  
  colnames(tmp_freq) <- c("edit_id", paste(sample, "_extent", sep = ""))
  colnames(tmp_depth) <- c("seq_id",  paste(sample, "_depth", sep = ""))
  
  all_sites_annotation <- left_join(all_sites_annotation, tmp_freq, by  = "edit_id") %>% 
    left_join(tmp_depth, by = "seq_id")
  
  # modify the editing extent (particularlly NA, representing 0 editing extent, therefore not included in the vcf output of Varscan2) information
  for (i in 1:nrow(all_sites_annotation)) {
    depth <- as.numeric(all_sites_annotation[i,ncol(all_sites_annotation)])
    freq <- all_sites_annotation[i,ncol(all_sites_annotation)-1]
    
    # is.na(depth) => depth <- 0: no reads mapped to this sites
    if (is.na(depth)) {
      all_sites_annotation[i,ncol(all_sites_annotation)] <- 0
      next
    }
      
    # depth > 0 & is.na(editing_extent) => editing_extent <- 0: no editing detected in this sites
    if (depth > 0 & is.na(freq)){
      all_sites_annotation[i,ncol(all_sites_annotation)-1] <- 0
    }
    
  }
  
}

all_sites_annotation <- select(all_sites_annotation, -seq_id)


##---------------- output -----------------------

save(all_sites_annotation,
     file = "annotated_sites_nofiltered.Rdata")
write.table(all_sites_annotation, file = "annotated_sites_nofiltered.table",
            col.names = T, row.names = F,
            quote = F, sep = "\t")


