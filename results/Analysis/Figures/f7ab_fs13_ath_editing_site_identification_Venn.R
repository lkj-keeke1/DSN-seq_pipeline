library(tidyverse)
library(ggvenn)

load("../RNA_editing/Ath/01.site_merge/plot.Rdata")

all_sites <- separate(all_sites, edit_id, into = c("seq", "posi", "edit_type"), sep = " ") %>% 
  separate(seq, into = c("gene", "chr"), sep = "_")

mt_plot <- filter(all_sites, chr == "mt")

pt_plot <- filter(all_sites, chr == "pt")


## DSN vs STS
p_mt_venn1 <- ggplot(mt_plot) +
  geom_venn(mapping = aes(
    A = `Ath_DSN_merge`,
    B = `Ath_STS`
  ),
  show_percentage = F,
  fill_color = c("#5b9ecf", "#7e56d9"),
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous() +
  theme_void() +
  labs(title = "Mitochondrial editing sites") +
  theme(plot.title = element_text(hjust = 0.5,
                                  #vjust = -10,
                                  size = 8,
                                  color = "#c04851"))


p_pt_venn1 <- ggplot(pt_plot) +
  geom_venn(mapping = aes(
    A = `Ath_DSN_merge`,
    B = `Ath_STS`
  ),
  show_percentage = F,
  fill_color = c("#5b9ecf", "#7e56d9"),
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous() +
  theme_void() +
  labs(title = "Chloroplast editing sites") +
  theme(plot.title = element_text(hjust = 0.5,
                                  #vjust = -10, 
                                  size = 8,
                                  color = "#3c9566"))
#save(file = "Ath_venn.Rdata", p_mt_venn1, p_pt_venn1)

ggsave(p_mt_venn1,filename = "f7a_DSN_STS_mt_editing_sites.pdf",width = 6,height =6.5,units = "cm")
ggsave(p_pt_venn1,filename = "f7b_DSN_STS_pt_editing_sites.pdf",width = 6,height =6.5,units = "cm")
# group_by(all_sites, seq) %>% 
#   summarise(count = n())



## DSN Venn

p_venn_dsn_mt <- ggplot(mt_plot) +
  geom_venn(mapping = aes(
    A = Ath_DSN_rep1,
    B = Ath_DSN_rep2,
    C = Ath_DSN_rep3,
    D = Ath_DSN_merge
  ),
  show_percentage = F,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()

p_venn_dsn_pt <- ggplot(pt_plot) +
  geom_venn(mapping = aes(
    A = Ath_DSN_rep1,
    B = Ath_DSN_rep2,
    C = Ath_DSN_rep3,
    D = Ath_DSN_merge
  ),
  show_percentage = F,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()


ggsave(plot = p_venn_dsn_mt,
       filename = "fs13a_DSN_ath_mt_editing_sites.pdf",
       width = 130, height = 75, units = "mm")
ggsave(plot = p_venn_dsn_pt,
       filename = "fs13b_DSN_ath_pt_editing_sites.pdf",
       width = 130, height = 75, units = "mm")




