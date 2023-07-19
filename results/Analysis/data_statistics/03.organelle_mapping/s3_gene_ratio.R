rm(list = ls())

library(tidyverse)
library(patchwork)

gene_info_file <- "~/projects/DSN-seq/data/2022-11-02_organelle_genes/merged_gene_info_table.txt"
# 读取 gene 信息表，提取对应染色体信息

edited <- read.table("edited_Pt_genes.txt") %>% pull(V1)
gene_info <- read.table(file = gene_info_file, sep = "\t", header = T) %>% 
  select(gene = name_in_ref, chr)

Ribo_rep1 <- read.table("Ribo_rep1.filtered.idxstats") %>% 
  select(gene = V1,
         count = V3) %>% 
  left_join(gene_info, by = "gene") %>% 
  filter(chr == "Mt") %>% 
  arrange(desc(count)) %>% 
  mutate(per = 100 * count / sum(count)) %>% 
  mutate(edited = gene %in% edited)

levels = Ribo_rep1$gene

DSN_rep2 <- read.table("DSN_rep2.filtered.idxstats") %>% 
  select(gene = V1,
         count = V3) %>% 
  left_join(gene_info, by = "gene") %>% 
  filter(chr == "Mt") %>% 
  arrange(desc(count)) %>% 
  mutate(per = 100 * count / sum(count)) %>% 
  mutate(edited = gene %in% edited)



DSN_rep2 <- mutate(DSN_rep2, gene = factor(gene, levels = levels))
Ribo_rep1 <- mutate(Ribo_rep1, gene = factor(gene, levels = levels))
group_by(DSN_rep2, edited) %>% 
  summarise(per = sum(per))
group_by(Ribo_rep1, edited) %>% 
  summarise(per = sum(per))
p1 <- ggplot(data = DSN_rep2) + 
  geom_col(mapping = aes(x = gene, y = per, fill = edited)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p2 <- ggplot(data = Ribo_rep1) + 
  geom_col(mapping = aes(x = gene, y = per, fill = edited)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p1/p2
