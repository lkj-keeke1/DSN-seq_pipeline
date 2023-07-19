rm(list = ls())

library(tidyverse)

load(file = "../RNA_editing/Osa/02.site_annotation/annotated_sites_filtered.Rdata")

# Editing extent distribution of newly identified editing sites vs canonical sites

## Mt sites
Mt_all_sites <- filter(all_sites_annotation, chr == "Mt",
                       DSN_det == T, Ribo_det == T) %>%
  select(canonical,
         Osa_DSN_merge_extent,
         `Osa_Ribo-off_merge_extent`) %>% 
  gather(key = sample,
         value = extent,
         -canonical)


mt_box <- ggplot(Mt_all_sites) +
  geom_boxplot(mapping = aes(x = canonical, y = extent, color = sample)) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_discrete(labels = c("non-classical", "classical")) + 
  scale_color_manual(values =  c("#5b9ecf", "#f59364"),
                     labels = c("Osa_DSN_merge", "Osa_Ribo-off_merge")) +
  theme_bw() +
  theme(axis.text.x = element_text(#angle = 30,
    #vjust = 1,
    #hjust = 1, 
    size = 8,
    color = "black"),
    axis.text.y = element_text(size = 8,
                               color = "black"),
    axis.title.y = element_text(size = 8,
                                color = "black"),
    legend.title = element_text(size = 8,
                                color = "black"),
    legend.text = element_text(size = 8,
                               color = "black"),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5,
                              size = 8),
    plot.margin = margin(t = 10, b = 10, 
                         l = 10, r = 10,
                         unit = "mm")) +
  labs(x ="",
       y = "RNA editing extent (%)",
       title = "")

ggsave(mt_box, filename = "f3b_mt_new_vs_old_extent.pdf",
       width = 120, height = 80, units = "mm")

## Pt sites
Pt_all_sites <- filter(all_sites_annotation, chr == "Pt",
                       DSN_det == T, Ribo_det == T) %>%
  select(canonical,
         Osa_DSN_merge_extent,
         `Osa_Ribo-off_merge_extent`) %>% 
  gather(key = sample,
         value = extent,
         -canonical)

pt_box <- ggplot(Pt_all_sites) +
  geom_boxplot(mapping = aes(x = canonical, y = extent, color = sample)) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_discrete(labels = c("non-classical", "classical")) + 
  scale_color_manual(values =  c("#5b9ecf", "#f59364"),
                     labels = c("Osa_DSN_merge", "Osa_Ribo-off_merge")) +
  theme_bw() +
  theme(axis.text.x = element_text(#angle = 30,
    #vjust = 1,
    #hjust = 1, 
    size = 8,
    color = "black"),
    axis.text.y = element_text(size = 8,
                               color = "black"),
    axis.title.y = element_text(size = 8,
                                color = "black"),
    legend.title = element_text(size = 8,
                                color = "black"),
    legend.text = element_text(size = 8,
                               color = "black"),
    #legend.position = "none",
    plot.title = element_text(hjust = 0.5,
                              size = 8),
    plot.margin = margin(t = 10, b = 10, 
                         l = 10, r = 10,
                         unit = "mm")) +
  labs(x ="",
       y = "RNA editing extent (%)",
       title = "")

ggsave(pt_box, filename = "f4b_pt_new_vs_old_extent.pdf",
       width = 120, height = 80, units = "mm")



# nad9-111
load(file = "../RNA_editing/Osa/02.site_annotation/annotated_sites_nofiltered.Rdata")


nad9_no_detect <- filter(all_sites_annotation, edit_id == "nad9_mt 111 CT") %>% 
  select(Osa_DSN_rep1_extent, Osa_DSN_rep2_extent, Osa_DSN_rep3_extent, Osa_DSN_merge_extent,
         `Osa_Ribo-off_rep1_extent`, `Osa_Ribo-off_rep2_extent`, `Osa_Ribo-off_rep3_extent`, `Osa_Ribo-off_merge_extent`) %>% 
  gather(key = sample, value = extent) %>% 
  mutate(method = rep(c("DSN", "Ribo"), each = 4))

nad9_extent <- ggplot(nad9_no_detect) +
  geom_col(mapping = aes(x = sample, y = extent, fill = method)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0,40)) +
  scale_x_discrete(labels = c(paste("Osa_DSN_", c("merge", "rep1", "rep2", "rep3"), sep = ""),
                              paste("Osa_Ribo-off_", c("merge", "rep1", "rep2", "rep3"), sep = "")
  )
  ) + 
  scale_fill_manual(values =  c("#5b9ecf", "#f59364")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30,
                                   vjust = 1,
                                   hjust = 1, 
                                   size = 8,
                                   color = "black"),
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.y = element_text(size = 8,
                                    color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  size = 8),
        plot.margin = margin(t = 10, b = 10, 
                             l = 10, r = 10,
                             unit = "mm")) +
  labs(x ="",
       y = "RNA editing extent (%)",
       title = "nad9-111")


ggsave(nad9_extent, filename = "fs6a_nad9_extent.pdf",
       width = 110, height = 80, units = "mm")
