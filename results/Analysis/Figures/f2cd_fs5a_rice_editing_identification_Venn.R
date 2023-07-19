library(tidyverse)

load("../RNA_editing/Osa/01.site_merge/DSN_editing_sites.Rdata")
load("../RNA_editing/Osa/01.site_merge/Ribo-off_editing_sites.Rdata")
load("../RNA_editing/Osa/01.site_merge/PolyA_editing_sites.Rdata")
load("../RNA_editing/Osa/01.site_merge/all_sites_merge.Rdata")


# DSN-seq
p_venn_dsn <- ggplot(DSN_sites) +
  geom_venn(mapping = aes(
    A = Osa_DSN_rep1,
    B = Osa_DSN_rep2,
    C = Osa_DSN_rep3,
    D = Osa_DSN_merge
  ),
  show_percentage = F,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()



ggsave(plot = p_venn_dsn,
       filename = "f2c_rice_DSN_identification_Venn.pdf",
       width = 130, height = 75, units = "mm")


# Ribo-off-seq
p_venn_ribo <- ggplot(Ribo_sites) +
  geom_venn(mapping = aes(
    A = `Osa_Ribo-off_rep1`,
    B = `Osa_Ribo-off_rep2`,
    C = `Osa_Ribo-off_rep3`,
    D = `Osa_Ribo-off_merge`
  ),
  show_percentage = F,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()

ggsave(plot = p_venn_ribo,
       filename = "f2d_rice_Ribo-off_identification_Venn.pdf",
       width = 130, height = 75, units = "mm")


# mRNA-seq
p_venn_polya <- ggplot(PolyA_sites) +
  geom_venn(mapping = aes(
    A = Osa_PolyA_rep1,
    B = Osa_PolyA_rep2,
    C = Osa_PolyA_rep3,
    D = Osa_PolyA_merge
  ),
  show_percentage = F,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()

ggsave(plot = p_venn_polya,
       filename = "fs5a_rice_PolyA_identification_Venn.pdf",
       width = 130, height = 75, units = "mm")

# among all three methods

p_venn_polya_influ <- ggplot(detection_plot2) +
  geom_venn(mapping = aes(
    A = `DSN-seq`,
    B = `Ribo-off-seq`,
    C = `mRNA-seq`
  ),
  show_percentage = F,
  fill_color = c("#5b9ecf", "#f59364", "#f16065"),
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()

ggsave(p_venn2,file = "PolyA_influence_venn.pdf",width = 9, height = 6, units = "cm")