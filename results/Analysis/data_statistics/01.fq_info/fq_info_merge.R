rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)

# input seqkit stat table
clean <- read.table("../../../P0_data_qc/03.seqkit/RNA_editing_UID_seqkit_stat.table", header = T) %>% 
  select(c(1,4,5))
raw <- read.table("../../../P0_data_qc/03.seqkit/RNA_editing_raw_seqkit_stat.table", header = T) %>% 
  select(c(1,4,5))


# calculation of fragment number and base number
raw <- mutate(raw,
              sample=unlist(str_split(raw$file, "\\."))[seq(1,nrow(raw)*5,5)]) %>% 
  group_by(sample) %>% 
  summarise(Raw_fragments=sum(num_seqs)/2,
            Raw_bases=sum(sum_len))

clean <- mutate(clean,
               sample=unlist(str_split(clean$file, "\\."))[seq(1,nrow(clean)*5,5)]) %>% 
  group_by(sample) %>% 
  summarise(Clean_fragments=sum(num_seqs)/2,
            Clean_bases=sum(sum_len))

#
number_base_tbl <- left_join(raw,clean, by = "sample") 
write.table(number_base_tbl,
            file = "read_base_number.table",
            quote = F, sep = "\t",
            row.names = F)

