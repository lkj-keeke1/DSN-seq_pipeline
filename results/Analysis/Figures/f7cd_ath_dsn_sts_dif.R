rm(list = ls())


library(tidyverse)

load("../RNA_editing/Ath/02.site_annotation/annotated_sites_nofiltered.Rdata")


all_extent <- mutate(all_sites_annotation, det = if_else(STS_det == T & DSN_det == F, "STS_only",
                                                         if_else(STS_det == F & DSN_det == T, "DSN_only",
                                                                 if_else(STS_det == T & DSN_det == T, "Both", "None")))) %>% 
  filter(det != "None") %>% 
  select(Ath_DSN_merge_extent, Ath_STS_extent, det, chr) %>% 
  gather(key = "sample", value = "extent", -chr, -det)

mt_plot <- filter(all_extent, chr == "Mt")
pt_plot <- filter(all_extent, chr == "Pt")

mt_extent <- ggplot(data = mt_plot) +
  geom_violin(mapping = aes(x = det, y = extent, color = sample), scale = "width", width = 0.5,
              position = position_dodge(width = 0.7)) +
  geom_point(mapping = aes(x = det, y = extent, color = sample, group = sample), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),alpha = 0.2) +
  geom_hline(yintercept = 10, linetype = "dashed", alpha = 0.7) +
  scale_color_manual(values = c("#5b9ecf", "#7e56d9"),
                     labels = c("Ath_DSN-seq", "Ath_STS-PCRseq")) +
  theme_bw() +
  theme(plot.title = element_text(size = 8, color = "#c04851"),
        axis.title.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.position = "none") +
  labs(x = "",
       y = "Editing extent",
       title = "Mitochondiral",
       color = "Sample")

pt_extent <- ggplot(data = filter(pt_plot, det != "DSN_only")) +
  geom_violin(mapping = aes(x = det, y = extent, color = sample), scale = "width", width = 0.5,
              position = position_dodge(width = 0.7)) +
  geom_point(mapping = aes(x = det, y = extent, color = sample, group = sample), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),alpha = 0.3) +
  geom_hline(yintercept = 10, linetype = "dashed", alpha = 0.7) +
  scale_color_manual(values = c("#5b9ecf", "#7e56d9"),
                     labels = c("Ath_DSN-seq", "Ath_STS-PCRseq")) +
  theme_bw() +
  theme(plot.title = element_text(size = 8, color = "#3c9566"),
        axis.title.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.position = "none") +
  labs(x = "",
       y = "Editing extent",
       title = "Chloroplast",
       color = "Sample")


p_out <- mt_extent / pt_extent

ggsave(file = "f7cd_ath_dsn_sts_extent_distribution.pdf",p_out,
       width = 120, height = 120, units = "mm")