library(tidyverse)
library(ggsci)
library(patchwork)


info <- read.table("../data_statistics/04.merge_table/all_info_merge.table", sep = "\t", header = T) %>% 
  filter(species == "nip")


mt_mapping <- group_by(info, method)%>% 
  summarise(mean = mean(Mt),
            sd = sd(Mt))

pt_mapping <- group_by(info, method)%>% 
  summarise(mean = mean(Pt),
            sd = sd(Pt))

p_box_mt <- ggplot(mt_mapping) +
  geom_col(aes(x=method, y=mean, fill = method),
           width = 0.8) +
  geom_errorbar(aes(x=method, ymin=mean-sd, ymax=mean+sd),
                width = 0.3) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1.7)) +
  scale_fill_manual(values = c("#5b9ecf", "#f59364", "#f16065")) +
  theme_classic() +
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
  labs(x ="",
       y = "Mt gene mapping rate (%)")


p_box_pt <- ggplot(pt_mapping) +
  geom_col(aes(x=method, y=mean, fill = method),
           width = 0.8) +
  geom_errorbar(aes(x=method, ymin=mean-sd, ymax=mean+sd),
                width = 0.3) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) +
  scale_fill_manual(values = c("#5b9ecf", "#f59364", "#f16065")) +
  theme_classic() +
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
  labs(x ="",
       y = "Pt gene mapping rate (%)")






box <- p_box_mt / p_box_pt
ggsave(box, filename = "f2ab_organelle_mapping_rate.pdf", width = 65, height = 120, units = "mm")
