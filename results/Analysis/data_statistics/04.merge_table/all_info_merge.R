rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)

# input info tables
sample_info <- read.table("../../../../data/seq_data/RNA_editing/RNA_editing_sample_info.txt", header = T)
fq_info <- read.table("../01.fq_info/read_base_number.table",header = T)
rRNA_info <- read.table("../02.rRNA_content/rRNA_content.table",header = T)
mapping_info <- read.table("../03.organelle_mapping/RNA_editing_organelle_mapping.table", header = T) 


# merge info tables
merge_table <- left_join(sample_info, fq_info, by = 'sample') %>%
  left_join(rRNA_info, by = 'sample') %>% 
  left_join(mapping_info, by = 'sample')

# output table
write.table(merge_table,
            file = "all_info_merge.table",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)
