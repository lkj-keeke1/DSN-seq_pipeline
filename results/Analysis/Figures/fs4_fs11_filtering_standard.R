rm(list = ls())

library(tidyverse)
library(patchwork)

##### rice #####
DSN_rep1 <- read.table("../../P1_call_editing/02.call_editing/Osa_DSN_rep1.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
DSN_rep2 <- read.table("../../P1_call_editing/02.call_editing/Osa_DSN_rep2.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
DSN_rep3 <- read.table("../../P1_call_editing/02.call_editing/Osa_DSN_rep3.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
DSN_all <-  read.table("../../P1_call_editing/02.call_editing/Osa_DSN_merge.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)

Ribo_rep1 <- read.table("../../P1_call_editing/02.call_editing/Osa_Ribo-off_rep1.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
Ribo_rep2 <- read.table("../../P1_call_editing/02.call_editing/Osa_Ribo-off_rep2.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
Ribo_rep3 <- read.table("../../P1_call_editing/02.call_editing/Osa_Ribo-off_rep3.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
Ribo_all <- read.table("../../P1_call_editing/02.call_editing/Osa_Ribo-off_merge.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)

PolyA_rep1 <- read.table("../../P1_call_editing/02.call_editing/Osa_PolyA_rep1.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
PolyA_rep2 <- read.table("../../P1_call_editing/02.call_editing/Osa_PolyA_rep2.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
PolyA_rep3 <- read.table("../../P1_call_editing/02.call_editing/Osa_PolyA_rep3.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
PolyA_all <- read.table("../../P1_call_editing/02.call_editing/Osa_PolyA_merge.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)

## plot

#DSN
DSN_rep1_plot <- ggplot(data = DSN_rep1) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "density",
       title = "Osa_DSN_rep1") +
  theme_bw() 

DSN_rep2_plot <- ggplot(data = DSN_rep2) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "",
       title = "Osa_DSN_rep2") +
  theme_bw()

DSN_rep3_plot <- ggplot(data = DSN_rep3) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "",
       title = "Osa_DSN_rep3") +
  theme_bw()

DSN_all_plot <- ggplot(data = DSN_all) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "",
       title = "Osa_DSN_merge") +
  theme_bw()

# Ribo-off
Ribo_rep1_plot <- ggplot(data = Ribo_rep1) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "density",
       title = "Osa_Ribo-off_rep1") +
  theme_bw()

Ribo_rep2_plot <- ggplot(data = Ribo_rep2) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "",
       title = "Osa_Ribo-off_rep2") +
  theme_bw()

Ribo_rep3_plot <- ggplot(data = Ribo_rep3) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "",
       title = "Osa_Ribo-off_rep3") +
  theme_bw()

Ribo_all_plot <- ggplot(data = Ribo_all) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "",
       y = "",
       title = "Osa_Ribo-off_merge") +
  theme_bw()

# PolyA
PolyA_rep1_plot <- ggplot(data = PolyA_rep1) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "density",
       title = "Osa_PolyA_rep1") +
  theme_bw()

PolyA_rep2_plot <- ggplot(data = PolyA_rep2) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "",
       title = "Osa_PolyA_rep2") +
  theme_bw()

PolyA_rep3_plot <- ggplot(data = PolyA_rep3) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "",
       title = "Osa_PolyA_rep3") +
  theme_bw()

PolyA_all_plot <- ggplot(data = PolyA_all) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "",
       title = "Osa_PolyA_merge") +
  theme_bw()


my_theme <-  theme(axis.text.x = element_text(size = 8,
                                              color = "black"),
                   axis.text.y = element_text(size = 8,
                                              color = "black"),
                   axis.title.x = element_text(size = 8,
                                               color = "black"),
                   axis.title.y = element_text(size = 8,
                                               color = "black"),
                   plot.title = element_text(size = 8,
                                             color = "black"),
                   legend.title = element_text(size = 8,
                                               color = "black"),
                   legend.text = element_text(size = 8,
                                              color = "black"))



DSN_rep1_plot <- DSN_rep1_plot + my_theme
DSN_rep2_plot <- DSN_rep2_plot + my_theme
DSN_rep3_plot <- DSN_rep3_plot + my_theme
DSN_all_plot <- DSN_all_plot + my_theme
Ribo_rep1_plot <- Ribo_rep1_plot + my_theme
Ribo_rep2_plot <- Ribo_rep2_plot + my_theme
Ribo_rep3_plot <- Ribo_rep3_plot + my_theme
Ribo_all_plot <- Ribo_all_plot + my_theme
PolyA_rep1_plot <- PolyA_rep1_plot + my_theme
PolyA_rep2_plot <- PolyA_rep2_plot + my_theme
PolyA_rep3_plot <- PolyA_rep3_plot + my_theme
PolyA_all_plot <- PolyA_all_plot + my_theme




layout <- 
  "
ABCD
EFGH
IJKL
"

pdf(file = "fs4_rice_non-CT_editing_extent.pdf", width = 7, height = 5.5)  
DSN_rep1_plot + DSN_rep2_plot + DSN_rep3_plot + DSN_all_plot +
  Ribo_rep1_plot + Ribo_rep2_plot + Ribo_rep3_plot + Ribo_all_plot +
  PolyA_rep1_plot + PolyA_rep2_plot + PolyA_rep3_plot + PolyA_all_plot +
  plot_layout(design = layout)

dev.off()



##### Col-0 #####
rm(list = ls())

library(tidyverse)
library(patchwork)


DSN_rep1 <- read.table("../../P1_call_editing/02.call_editing/Ath_DSN_rep1.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
DSN_rep2 <- read.table("../../P1_call_editing/02.call_editing/Ath_DSN_rep2.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
DSN_rep3 <- read.table("../../P1_call_editing/02.call_editing/Ath_DSN_rep3.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)
DSN_merge <- read.table("../../P1_call_editing/02.call_editing/Ath_DSN_merge.out.table", sep = "\t", header = T) %>% 
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)


STS <- read.table("../../sup_STS/02.call_editing/Ath_STS.out.table", sep = "\t", header = T) %>%
  filter(edit != "CT", 
         posi >0, posi <= length,
         pvalue < 0.01)






p_dsn_1 <- ggplot(data = DSN_rep1) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "density",
       title = "Ath_DSN_rep1") +
  theme_bw()

p_dsn_2 <- ggplot(data = DSN_rep2) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "density",
       title = "Ath_DSN_rep2") +
  theme_bw()

p_dsn_3 <- ggplot(data = DSN_rep3) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "density",
       title = "Ath_DSN_rep3") +
  theme_bw()

p_dsn_all <- ggplot(data = DSN_merge) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "density",
       title = "Ath_DSN_merge") +
  theme_bw()

p_sts <- ggplot(data = STS) +
  geom_density(mapping = aes(x = freq),
               size = 0.5,
               color = "red") +
  geom_vline(xintercept = 10,
             linetype = "dashed") +
  scale_x_continuous(limits = c(-0.5,25),
                     breaks = c(0, 5, 10, 15, 20, 25)) +
  labs(x = "editing extent (%)",
       y = "density",
       title = "Ath_STS") +
  theme_bw()

my_theme <-  theme(axis.text.x = element_text(size = 8,
                                              color = "black"),
                   axis.text.y = element_text(size = 8,
                                              color = "black"),
                   axis.title.x = element_text(size = 8,
                                               color = "black"),
                   axis.title.y = element_text(size = 8,
                                               color = "black"),
                   plot.title = element_text(size = 8,
                                             color = "black"),
                   legend.title = element_text(size = 8,
                                               color = "black"),
                   legend.text = element_text(size = 8,
                                              color = "black"))


p_dsn_1 <- p_dsn_1 + my_theme
p_dsn_2 <- p_dsn_2 + my_theme
p_dsn_3 <- p_dsn_3 + my_theme
p_dsn_all <- p_dsn_all + my_theme
p_sts <- p_sts + my_theme

layout <- "
AABBCC
DDEE##
"

p_all <- p_dsn_1 + p_dsn_2 + p_dsn_3 + p_dsn_all + p_sts + patchwork::plot_layout(design = layout)

ggsave(p_all, filename = "fs11_Ath_non-CT_editing_extent.pdf", width = 18, height = 12, unit = "cm")

