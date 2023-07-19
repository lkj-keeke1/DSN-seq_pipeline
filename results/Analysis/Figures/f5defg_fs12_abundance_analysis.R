##### Sequencing volume abundance analysis #####


##-------------------- rice --------------------
library(tidyverse)
library(patchwork)

load("~/projects/DSN-seq/results/P4_abundance_analysis/4.merge_info/Osa_edit_all.Rdata")


DSN_20M <- 10 * (20000000 / 43011843)
Ribo_20M <- 10 * (20000000 / 51172446)

edit_all <- filter(edit_all, classical == T) %>% separate(col = edit_id, into = c("gene", "posi", "type"), sep = " ") %>% 
  separate(col = gene, into = c("gene", "chr"), sep = "_") %>% select(-gene, -posi, -type)

mt_all <- filter(edit_all, chr == "mt") %>% 
  select(-chr)
mt_all <- as.data.frame(t(mt_all))
mt_all <- mutate(mt_all, count = rowSums(mt_all)) %>% 
  select(count) 
mt_plot <- rownames_to_column(mt_all, var = "sample_id") %>% 
  separate(col = "sample_id", sep = "\\.", into = c("sample", "sample_way")) %>% 
  separate(col = "sample", sep = "_", into = c("species", "method", "rep")) %>% 
  separate(col = "sample_way", sep = "_", into = c("seed", "frac"))



pt_all <- filter(edit_all, chr == "pt") %>% 
  select(-chr)
pt_all <- as.data.frame(t(pt_all))
pt_all <- mutate(mt_all, count = rowSums(pt_all)) %>% 
  select(count)
pt_plot <- rownames_to_column(pt_all, var = "sample_id") %>% 
  separate(col = "sample_id", sep = "\\.", into = c("sample", "sample_way")) %>% 
  separate(col = "sample", sep = "_", into = c("species", "method","rep")) %>% 
  separate(col = "sample_way", sep = "_", into = c("seed", "frac"))


my_theme <-   theme(text = element_text(size = 8, color = "black"),
                    axis.text.x = element_text(size = 8, color = "black"),
                    axis.text.y = element_text(size = 8, color = "black"),
                    axis.title.x = element_text(size = 8, color = "black"),
                    axis.title.y = element_text(size = 8, color = "black"),
                    plot.title = element_text(size = 8, color = "black"),
                    legend.title = element_text(size = 8, color = "black"),
                    legend.text = element_text(size = 8, color = "black"))


# Mt
## DSN
p_dsn_mt_abundance <- ggplot(data = filter(mt_plot, method == "DSN", rep == "merge")) +
  geom_boxplot(mapping = aes(x = frac, y = count),
               color = "#c04851") +
  geom_vline(xintercept = DSN_20M, linetype = "dashed", alpha = 0.5) +
  # geom_text(x = DSN_20M + 0.7, y = 475, label = "~ 20M clean reads", alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(450, 575)) +
  scale_x_discrete(labels = c(paste(1:9, "0", sep = ""))) +
  theme_bw() +
  labs(x = "sample size (%)",
       y = "detected Mt classical sites",
       title = "Osa_DSN_merge") + my_theme +
  theme(plot.title = element_text(colour = "#5b9ecf"))

## Ribo
p_ribo_mt_abundance <- ggplot(data = filter(mt_plot, method == "Ribo-off", rep == "merge")) +
  geom_boxplot(mapping = aes(x = frac, y = count),
               color = "#c04851") +
  geom_vline(xintercept = Ribo_20M, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(450, 575)) +
  scale_x_discrete(labels = c(paste(1:9, "0", sep = ""))) +
  theme_bw() +
  labs(x = "sample size (%)",
       y = "detected Mt classical sites",
       title = "Osa_Ribo-off_merge") + my_theme +
  theme(plot.title = element_text(colour = "#f59364"))

# Pt
## DSN
p_dsn_pt_abundance <- ggplot(data = filter(pt_plot, method == "DSN", rep == "merge")) +
  geom_boxplot(mapping = aes(x = frac, y = count),
               color = "#3c9566") +
  geom_vline(xintercept = DSN_20M, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(26, 33)) +
  scale_x_discrete(labels = c(paste(1:9, "0", sep = ""))) +
  theme_bw() +
  labs(x = "sample size (%)",
       y = "detected Pt classical sites",
       title = "Osa_DSN_merge") + my_theme +
  theme(plot.title = element_text(colour = "#5b9ecf"))
## Ribo
p_ribo_pt_abundance <- ggplot(data = filter(pt_plot, method == "Ribo-off", rep == "merge")) +
  geom_boxplot(mapping = aes(x = frac, y = count),
               color = "#3c9566") +
  geom_vline(xintercept = Ribo_20M, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(26, 33)) +
  scale_x_discrete(labels = c(paste(1:9, "0", sep = ""))) +
  theme_bw() +
  labs(x = "sample size (%)",
       y = "detected Pt classical sites",
       title = "Osa_Ribo-off_merge") + my_theme +
  theme(plot.title = element_text(colour = "#f59364"))

layout <- "
AB
CD"

p_abundance <- p_dsn_mt_abundance + p_dsn_pt_abundance + p_ribo_mt_abundance + p_ribo_pt_abundance + patchwork::plot_layout(design = layout)


ggsave(p_abundance, filename = "f5defg_rice_abundance.pdf",
       width = 180, height = 120, units = "mm")

write.table(mt_plot,file = "rice_mt_abundance.table",
            row.names = F, col.names = F,
            sep = "\t", quote = F)

write.table(pt_plot,file = "rice_pt_abundance.table",
            row.names = F, col.names = F,
            sep = "\t", quote = F)



##-------------------- Arabidopsis --------------------
rm(list = ls())
library(tidyverse)
library(patchwork)

load("~/projects/DSN-seq/results/P4_abundance_analysis/4.merge_info/Ath_edit_all.Rdata")


edit_all <- edit_all %>% separate(col = edit_id, into = c("gene", "posi", "type"), sep = " ") %>% 
  separate(col = gene, into = c("gene", "chr"), sep = "_") %>% select(-gene, -posi, -type)

mt_all <- filter(edit_all, chr == "mt") %>% 
  select(-chr)
mt_all <- as.data.frame(t(mt_all))
mt_all <- mutate(mt_all, count = rowSums(mt_all)) %>% 
  select(count) 
mt_plot <- rownames_to_column(mt_all, var = "sample_id") %>% 
  separate(col = "sample_id", sep = "\\.", into = c("sample", "sample_way")) %>% 
  separate(col = "sample", sep = "_", into = c("species", "method", "rep")) %>% 
  separate(col = "sample_way", sep = "_", into = c("seed", "frac"))



pt_all <- filter(edit_all, chr == "pt") %>% 
  select(-chr)
pt_all <- as.data.frame(t(pt_all))
pt_all <- mutate(mt_all, count = rowSums(pt_all)) %>% 
  select(count)
pt_plot <- rownames_to_column(pt_all, var = "sample_id") %>% 
  separate(col = "sample_id", sep = "\\.", into = c("sample", "sample_way")) %>% 
  separate(col = "sample", sep = "_", into = c("species", "method", "rep")) %>% 
  separate(col = "sample_way", sep = "_", into = c("seed", "frac"))


my_theme <-   theme(text = element_text(size = 8, color = "black"),
                    axis.text.x = element_text(size = 8, color = "black"),
                    axis.text.y = element_text(size = 8, color = "black"),
                    axis.title.x = element_text(size = 8, color = "black"),
                    axis.title.y = element_text(size = 8, color = "black"),
                    plot.title = element_text(size = 8, color = "black"),
                    legend.title = element_text(size = 8, color = "black"),
                    legend.text = element_text(size = 8, color = "black"))


# Mt
## DSN
p_dsn_mt_abundance <- ggplot(data = filter(mt_plot, method == "DSN", rep == "merge")) +
  geom_boxplot(mapping = aes(x = frac, y = count),
               color = "#c04851") +
  #geom_vline(xintercept = DSN_20M, linetype = "dashed", alpha = 0.5) +
  # geom_text(x = DSN_20M + 0.7, y = 475, label = "~ 20M clean reads", alpha = 0.7, size = 3) +
  scale_y_continuous(limits = c(250, 450)) +
  scale_x_discrete(labels = c(paste(1:9, "0", sep = ""))) +
  theme_bw() +
  labs(x = "sample size (%)",
       y = "detected Mt editing sites",
       title = "Ath_DSN_merge") + my_theme +
  theme(plot.title = element_text(colour = "#5b9ecf"))

# Pt
## DSN
p_dsn_pt_abundance <- ggplot(data = filter(pt_plot, method == "DSN", rep == "merge")) +
  geom_boxplot(mapping = aes(x = frac, y = count),
               color = "#3c9566") +
  #geom_vline(xintercept = DSN_20M, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(24, 35)) +
  scale_x_discrete(labels = c(paste(1:9, "0", sep = ""))) +
  theme_bw() +
  labs(x = "sample size (%)",
       y = "detected Pt editing sites",
       title = "Ath_DSN_merge") + my_theme +
  theme(plot.title = element_text(colour = "#5b9ecf"))

layout <- "
A
B
"

p_abundance <- p_dsn_mt_abundance + p_dsn_pt_abundance + patchwork::plot_layout(design = layout)

ggsave(p_abundance, filename = "fs12_ath_abundance.pdf",
       width = 90, height = 120, units = "mm")

write.table(mt_plot,file = "ath_mt_abundance.table",
            row.names = F, col.names = F,
            sep = "\t", quote = F)

write.table(pt_plot,file = "ath_pt_abundance.table",
            row.names = F, col.names = F,
            sep = "\t", quote = F)
