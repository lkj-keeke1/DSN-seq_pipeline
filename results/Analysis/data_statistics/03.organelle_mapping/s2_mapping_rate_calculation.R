rm(list = ls())

library(tidyverse)

Osa_gene_info_file <- "~/projects/DSN-seq/data/ref/Osa/Osa_organelle_gene_info.txt"
Ath_gene_info_file <- "~/projects/DSN-seq/data/ref/Ath/Ath_organelle_gene_info.txt"
flagstat_file <- "RNA_editing_organelle_mapping.summary"

# input gene info table
Osa_gene_info <- read.table(file = Osa_gene_info_file, sep = "\t", header = T) %>% 
  select(gene = name_in_ref, chr)
Ath_gene_info <- read.table(file = Ath_gene_info_file, sep = "\t", header = T) %>% 
  select(gene = name_in_ref, chr)

# input merged flagstat table
colnames <- c('sample',
              'total_mapping_records',
              'secondary',
              'suplementary',
              'duplicates',
              'mapped_records',
              'paired_reads',
              'read1',
              'read2',
              'properly_paired_reads',
              'R1+R2_mapped',
              'singletons',
              'different_chr_mapping',
              'different_chr_mapping(mapQ>=5)')
flagstat <- as.data.frame(t(read.table(file = flagstat_file, sep = "\t", header = F)))
colnames(flagstat) <- colnames
flagstat <- select(flagstat, sample, paired_reads)


# rice
Osa_samples <- c("Osa_DSN_rep1", "Osa_DSN_rep2", "Osa_DSN_rep3",
             "Osa_Ribo-off_rep1", "Osa_Ribo-off_rep2", "Osa_Ribo-off_rep3",
             "Osa_PolyA_rep1", "Osa_PolyA_rep2", "Osa_PolyA_rep3")

for (sample in Osa_samples) {
  idxstats_file <- paste("../../../P1_call_editing/01.mapping/", sample, ".filtered.idxstats", sep = "")
  
  idxstats <- read.table(file = idxstats_file,sep = "\t") %>% 
    select(V1, V3)
  colnames(idxstats) <- c("gene", sample)
  
  Osa_gene_info <- left_join(Osa_gene_info, idxstats, by = "gene")
}

Osa_mapping_rate <- gather(Osa_gene_info, key = "sample", value = "read_num", -chr, -gene) %>% 
  na.omit() %>% 
  group_by(sample, chr) %>% 
  summarise(count = sum(read_num)) %>% 
  left_join(flagstat, by = "sample") %>% 
  mutate(paired_reads = as.numeric(paired_reads),
         mapping_rate = round(100 * count/paired_reads, 2)) %>% 
  select(mapping_rate, sample, chr, mapping_rate) %>%
  spread(key = "chr", value = "mapping_rate")


# Arabidopsis
Ath_samples <- c("Ath_DSN_rep1", "Ath_DSN_rep2", "Ath_DSN_rep3")

for (sample in Ath_samples) {
  idxstats_file <- paste("../../../P1_call_editing/01.mapping/", sample, ".filtered.idxstats", sep = "")
  
  idxstats <- read.table(file = idxstats_file,sep = "\t") %>% 
    select(V1, V3)
  colnames(idxstats) <- c("gene", sample)
  
  Ath_gene_info <- left_join(Ath_gene_info, idxstats, by = "gene")
}

Ath_mapping_rate <- gather(Ath_gene_info, key = "sample", value = "read_num", -chr, -gene) %>% 
  na.omit() %>% 
  group_by(sample, chr) %>% 
  summarise(count = sum(read_num)) %>% 
  left_join(flagstat, by = "sample") %>% 
  mutate(paired_reads = as.numeric(paired_reads),
         mapping_rate = round(100 * count/paired_reads, 2)) %>% 
  select(mapping_rate, sample, chr, mapping_rate) %>%
  spread(key = "chr", value = "mapping_rate")


## output tables

mapping_rate <- full_join(Osa_mapping_rate, Ath_mapping_rate)
write.table(x = mapping_rate,
            file = "RNA_editing_organelle_mapping.table",
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")