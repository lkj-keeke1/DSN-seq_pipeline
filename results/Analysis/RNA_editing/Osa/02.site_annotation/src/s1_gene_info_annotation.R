# annotate gene information


rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)


##---------------- input -----------------------

load("../01.site_merge/all_sites_merge.Rdata") 
gene_info_file <- "~/projects/DSN-seq/data/ref/Osa/Osa_organelle_gene_info.txt"

# input sites and gene info 
all_sites <- mutate(all_sites, edit_id_tmp = edit_id) %>% 
  separate(edit_id_tmp, c("seq", "posi", "editing_type"), sep = " ") %>% 
  mutate(posi = as.numeric(posi)) %>% 
  arrange(seq, posi) 

gene_info <-  read.table(gene_info_file, header = T, sep = "\t", quote = "") 


##---------------- data process -----------------------

# annotate gene info 
all_sites <- left_join(all_sites,gene_info,by = c("seq" = "name_in_ref")) %>% 
  mutate(codon_num = NA,
         codon_posi = NA,
         codon_start_bed = NA,
         codon_end_bed = NA)


# calculate codon position
for (i in 1:nrow(all_sites)) {
  
  if (all_sites$gene_biotype[i] %in% c("pseudogene", "protein_coding")) {
    all_sites$codon_num[i] <- (all_sites$posi[i] - 1) %/% 3 + 1
    all_sites$codon_posi[i] <- (all_sites$posi[i] - 1) %% 3 + 1
    all_sites$codon_start_bed[i] <- (all_sites$codon_num[i] - 1) * 3 
    all_sites$codon_end_bed[i] <- all_sites$codon_num[i] * 3
    } 

}


##---------------- output -----------------------

# output codon table and annotated table

## codon table
codon_table <- select(all_sites, edit_id, seq, posi, codon_start_bed, codon_end_bed) %>%
  na.omit() %>% 
  # The sequences in ref have upstream 100 bp flanking sequences
  mutate(codon_start_bed = codon_start_bed + 100,
         codon_end_bed = codon_end_bed + 100) %>% 
  select(seq, codon_start_bed, codon_end_bed, edit_id)

write.table(file = "tmp_codon_tk.bed",
            x = codon_table,
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)

## annotated table
save(all_sites,
     file = "tmp1_gene_info_annotated.Rdata")




