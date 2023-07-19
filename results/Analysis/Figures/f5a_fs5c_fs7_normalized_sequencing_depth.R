# normalized sequencing depth = sequencing depth per 10M clean fragments

rm(list = ls())

library(tidyverse)
library(patchwork)


# polyA influence
load("../RNA_editing/Osa/02.site_annotation//annotated_sites_filtered.Rdata")

PolyA_clean <- 101592102
Ribo_clean <- 51172446
DSN_clean <- 43011843

depth_polya <- select(all_sites_annotation,
                  edit_id, PolyA_det, 
                  DSN_all_d = Osa_DSN_merge_depth, 
                  Ribo_all_d = `Osa_Ribo-off_merge_depth`, 
                  PolyA_all_d = `Osa_PolyA_merge_depth`) %>% 
  filter(PolyA_det == T) %>% 
  mutate(DSN_all_d = 10000000 * DSN_all_d / DSN_clean,
         Ribo_all_d = 10000000 * Ribo_all_d / Ribo_clean,
         PolyA_all_d = 10000000 * PolyA_all_d / PolyA_clean) %>% 
  mutate(polyA_influ_DSN = 100 * PolyA_all_d / DSN_all_d,
         polyA_influ_Ribo = 100 *PolyA_all_d / Ribo_all_d) %>% 
  select(edit_id, polyA_influ_DSN, polyA_influ_Ribo) %>% 
  gather(key = "method", value = "influence", -edit_id)


p_box <- ggplot(depth_polya) +
  geom_boxplot(aes(x = method, y = influence, color = method),
               size = 0.5) +
  geom_jitter(aes(x = method, y = influence),
              width = 0.2, alpha = .3, size = 1) +
  scale_x_discrete(labels = c("PolyA_merge/\nDSN_merge", "PolyA_merge/\nRibo-off_merge")) +
  scale_color_manual(values = c("#5b9ecf", "#f59364")) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.position = "") +
  labs(x = "",
       y = "relative depth ratio (%)")

ggsave(p_box, filename = "fs5c_polyA_influence.pdf",
       width = 60, height = 60, units = "mm")





# normalized

# calculate clean fragment number of samples 
f_num <- read.table("../data_statistics/01.fq_info/read_base_number.table",
                    header = T) %>% 
  select(sample, Clean_fragments)
sample_info <- read.table("~/projects/DSN-seq/data/seq_data/RNA_editing/RNA_editing_sample_info.txt", header = T, sep = "\t")

f_num <- left_join(f_num, sample_info, by = "sample") %>% 
  filter(project == "DSN_comparision") %>% 
  select(sample, Clean_fragments, method)

merge_num <- group_by(f_num, method) %>% 
  summarise(Clean_fragments = sum(Clean_fragments)) %>% 
  dplyr::rename(sample = method)

f_num <- select(f_num, -method) %>% 
  rbind(merge_num)


# only analyze the new canonical sites
load("../RNA_editing/Osa/02.site_annotation/annotated_sites_nofiltered.Rdata")


classical_sites_file <- "~/projects/DSN-seq/data/ref/Osa/Osa_new_canonical_sites.txt"

classcial_sites <- read.table(classical_sites_file, header = T, sep = "\t") %>% 
  mutate(chr = if_else(Organelle == "Mt", "_mt", "_pt"),
         edit_id = paste(Gene, chr, sep = ""),
         edit_id = paste(edit_id, `Position.in.CDS`, "CT", sep = " ")) %>% pull(edit_id)

all_sites_annotation <- all_sites_annotation[all_sites_annotation$edit_id %in% classcial_sites,]



all_sites_annotation <- select(all_sites_annotation,
                               chr,
                               edit_id,
                               Osa_DSN_rep1 = Osa_DSN_rep1_depth,
                               Osa_DSN_rep2 = Osa_DSN_rep2_depth,
                               Osa_DSN_rep3 = Osa_DSN_rep3_depth,
                               `DSN-seq` = Osa_DSN_merge_depth,
                               `Osa_Ribo-off_rep1` = `Osa_Ribo-off_rep1_depth`,
                               `Osa_Ribo-off_rep2` = `Osa_Ribo-off_rep2_depth`,
                               `Osa_Ribo-off_rep3` = `Osa_Ribo-off_rep3_depth`,
                               `Ribo-off-seq` = `Osa_Ribo-off_merge_depth`) %>% 
  gather(key = "sample", value = "depth", -chr, -edit_id) %>% 
  left_join(f_num, by = "sample") %>% 
  mutate(normalized_depth = 10000000 * depth / Clean_fragments)



# rep1
rep1_sequencing_depth <- filter(all_sites_annotation, sample %in% c("Osa_DSN_rep1", "Osa_Ribo-off_rep1"))

rep1_depth_plot <- ggplot(data = rep1_sequencing_depth) + 
  geom_boxplot(mapping = aes(x = sample, y = log2(normalized_depth), color = chr)) +
  scale_x_discrete(labels = c("Osa_DSN_rep1", "Osa_Ribo-off_rep1")) +
  scale_color_manual(values = c("#c04851", "#3c9566")) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black",
                                   angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(x = "",
       color = "Organelle")

# rep2
rep2_sequencing_depth <- filter(all_sites_annotation, sample %in% c("Osa_DSN_rep2", "Osa_Ribo-off_rep2"))

rep2_depth_plot <- ggplot(data = rep2_sequencing_depth) + 
  geom_boxplot(mapping = aes(x = sample, y = log2(normalized_depth), color = chr)) +
  scale_x_discrete(labels = c("Osa_DSN_rep2", "Osa_Ribo-off_rep2")) +
  scale_color_manual(values = c("#c04851", "#3c9566")) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black",
                                   angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(x = "",
       color = "Organelle")

# rep3
rep3_sequencing_depth <-  filter(all_sites_annotation, sample %in% c("Osa_DSN_rep3", "Osa_Ribo-off_rep3"))

rep3_depth_plot <- ggplot(data = rep3_sequencing_depth) + 
  geom_boxplot(mapping = aes(x = sample, y = log2(normalized_depth), color = chr)) +
  scale_x_discrete(labels = c("Osa_DSN_rep3", "Osa_Ribo-off_rep3")) +
  scale_color_manual(values = c("#c04851", "#3c9566")) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black",
                                   angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(x = "",
       color = "Organelle")


depth_plot <- rep1_depth_plot + rep2_depth_plot + rep3_depth_plot + patchwork::plot_layout(guides = "collect")

ggsave(filename = "fs7_rep123_normalized_depth.pdf", depth_plot,
       width = 180, height = 60, unit = "mm")



## merge

merge_sequencing_depth <-  filter(all_sites_annotation, sample %in% c("DSN-seq", "Ribo-off-seq"))

merge_depth_plot <- ggplot(data = merge_sequencing_depth) + 
  geom_boxplot(mapping = aes(x = sample, y = log2(normalized_depth), color = chr)) +
  scale_x_discrete(labels = c("Osa_DSN_merge", "Osa_Ribo-off_merge")) +
  scale_color_manual(values = c("#c04851", "#3c9566")) +
  theme_bw() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black",
                                   angle = 20, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black")) +
  labs(x = "",
       color = "Organelle")

ggsave(filename = "f5a_merge_normalized_depth.pdf", merge_depth_plot,
       width = 75, height = 60, unit = "mm")


t_normalized_depth <- select(all_sites_annotation, -Clean_fragments, -depth) %>% 
  spread(key = "sample", value = "normalized_depth")

write.table(t_normalized_depth, file = "normalized_depth.table",
            row.names = F, col.names = T,
            sep = "\t", quote = F)
