rm(list = ls())

library(tidyverse)
library(patchwork)


edit_gene <- read.table("~/projects/DSN-seq/data/ref/Osa/Osa_new_canonical_sites.txt",
                        sep = "\t", header = T) %>% 
  mutate(refname = if_else(Organelle == "Mt", "mt", "pt"),
         refname = paste(Gene, refname, sep = "_")) %>% 
  distinct(refname) %>% 
  pull(refname)

sample_info <- read.table("~/projects/DSN-seq/data/seq_data/RNA_editing/RNA_editing_sample_info.txt",
                          sep = "\t", header = T) %>% 
  select(sample, method)


exp <- read.table("../RNA_editing/Osa/03.gene_ratio/merge_expression.table", sep = "\t", header = F) %>% 
  filter(V1 != "*") %>% 
  mutate(tmp = gene) %>% 
  separate(col = tmp, sep = "_", into = c("name", "chr")) %>% 
  filter(chr == "pt") %>% 
  select(-name, -chr) %>% 
  column_to_rownames(var = "gene") %>% 
  select(- length) 

read_sum <-  colSums(exp)
for (i in 1:ncol(exp)){
  exp[,i] <- 100 * exp[,i] /read_sum[i]
}

exp_l <- rownames_to_column(exp, var = "gene") %>% 
  mutate(edit = gene %in% edit_gene) %>% 
  gather(key = sample, value = frac, -gene, -edit)

exp_l$sample[which(exp_l$sample == "Osa_Ribo.off_rep1")] <- "Osa_Ribo-off_rep1"
exp_l$sample[which(exp_l$sample == "Osa_Ribo.off_rep2")] <- "Osa_Ribo-off_rep2"
exp_l$sample[which(exp_l$sample == "Osa_Ribo.off_rep3")] <- "Osa_Ribo-off_rep3"

merge_plot <- group_by(exp_l, sample, edit) %>% 
  summarise(frac = sum(frac)) %>% 
  filter(edit == T) %>% 
  left_join(sample_info, by = "sample") %>% 
  group_by(method) %>% 
  summarise(mean = mean(frac),
            sd = sd(frac)) %>% 
  filter(method != "mRNA-seq")


# reads mapped to edited genes
p_edit_reads <- ggplot(data = merge_plot) +
  geom_errorbar(mapping = aes(x = method, ymin = mean-sd, ymax = mean+sd),
                width = 0.3) +
  geom_col(mapping = aes(x = method, y = mean, fill = method),
           width = 0.7) +
  scale_fill_manual(values = c("#5b9ecf", "#f59364")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8,
                                   color = "black"),
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.y = element_text(size = 8,
                                    color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  size = 8))+
  # plot.margin = margin(t = 10, b = 10, 
  #                      l = 10, r = 10,
  #                      unit = "mm")) +
  labs(x = "",
       y = "reads mapped to edited Pt genes/\n reads mapped to Pt genes (%)")

# ggsave(p_edit_reads, filename = "reads_on_edited_gene.pdf",
#        width = 60, height = 60, units = "mm")

# reads mapped to psbA_pt
single_plot <- left_join(exp_l, sample_info, by = "sample") %>% 
  group_by(method, gene, edit) %>% 
  summarise(mean = mean(frac)) %>% 
  arrange(method, edit, desc(mean)) %>% 
  filter(method != "mRNA-seq")


# mapping ratio of genes
p_single_genes_ratio <- ggplot(data = single_plot) +
  geom_jitter(mapping = aes(x = method, y = mean),
              width = 0.2) +
  scale_y_continuous(limits = c(0,75)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8,
                                   color = "black"),
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.y = element_text(size = 8,
                                    color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  size = 8)) +
  labs(x = "",
       y = "read count ratio of Pt genes (%)")


psbA <- left_join(exp_l, sample_info, by = "sample") %>% 
  filter(gene == "psbA_pt") %>% 
  group_by(method, gene, edit) %>% 
  summarise(mean = mean(frac), sd = sd(frac)) %>% 
  filter(method != "mRNA-seq") 


# reads mapped to psbA
p_psba <- ggplot(psbA)+
  geom_errorbar(mapping = aes(x = method, ymin = mean-sd, ymax = mean+sd),
                width = 0.3) +
  geom_col(mapping = aes(x = method, y = mean, fill = method),
           width = 0.7) +
  scale_fill_manual(values = c("#5b9ecf", "#f59364")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8,
                                   color = "black"),
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.y = element_text(size = 8,
                                    color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  size = 8)) +
  # plot.margin = margin(t = 10, b = 10, 
  #                      l = 10, r = 10,
  #                      unit = "mm")) +
  labs(x = "",
       y = "reads mapped to psbA_pt/\n reads mapped to Pt genes (%)")



ggsave(p_single_genes_ratio, filename = "fs8_single_genes_ratio.pdf",
       width = 80, height = 60, units = "mm")

p_read_ratio <- p_edit_reads + p_psba
# save(p_edit_reads, p_psba, file = "ribo_pt_plot.Rdata")

ggsave(p_read_ratio, filename = "f5bc_rice_gene_ratio_pt.pdf",
       width = 120, height = 60, units = "mm")
