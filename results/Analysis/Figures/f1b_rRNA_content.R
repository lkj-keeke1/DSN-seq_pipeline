library(tidyverse)
library(ggsci)
library(patchwork)


info <- read.table("../data_statistics/04.merge_table/all_info_merge.table", sep = "\t", header = T) %>% 
  filter(species == "nip")


content_mean <- group_by(info, method) %>% 
  summarise(mean = mean(rRNA_content),
            sd = sd(rRNA_content)) %>% 
  mutate(method = factor(method, levels = c("DSN-seq", "Ribo-off-seq", "mRNA-seq")))

p_bar <- ggplot(data = content_mean) +
  geom_col(aes(x=method, y=mean, fill = method),
           width = 0.8) +
  geom_errorbar(aes(x=method, ymin=mean-sd, ymax=mean+sd),
                width = 0.3) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 7)) +
  theme_classic() +
  scale_fill_manual(values = c("#5b9ecf", "#f59364", "#f16065")) +
  theme(axis.text.x = element_text(angle = 30,
                                   vjust = 1,
                                   hjust = 1, 
                                   size = 8,
                                   color = "black"),
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.y = element_text(size = 8,
                                    color = "black"),
        legend.position = "none") +
  labs(x = "",
       y = "rRNA content (%)",
       title = "")

ggsave("f1b_rRNA_content.pdf",
       p_bar, width = 55, height = 65, units = "mm")
