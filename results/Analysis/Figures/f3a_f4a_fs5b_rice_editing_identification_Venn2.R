rm(list = ls())

library(tidyverse)
library(ggvenn)

load("../RNA_editing/Osa/02.site_annotation/annotated_sites_nofiltered.Rdata")

venn_plot <- select(all_sites_annotation,
                    chr,
                    `Osa_DSN_merge` = DSN_det,
                    `Osa_Ribo-off_merge` = Ribo_det,
                    `Osa_PolyA_merge` = PolyA_det,
                    Canonical = canonical
                    )


### Mt classical sites
p_venn_mt <- ggplot(filter(venn_plot, chr == "Mt")) +
  geom_venn(mapping = aes(
    A = `Osa_DSN_merge`,
    B = `Osa_Ribo-off_merge`,
    C = `Canonical`
  ),
  show_percentage = F,
  fill_color = c("#5b9ecf", "#f59364", "#ddc55a"),
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()

ggsave(filename = "f3a_Mt_classical_detection.pdf",
       p_venn_mt,
       width = 85, height = 60, units = "mm")

### Pt classical sites
p_venn_pt <- ggplot(filter(venn_plot, chr == "Pt")) +
  geom_venn(mapping = aes(
    A = `Osa_DSN_merge`,
    B = `Osa_Ribo-off_merge`,
    C = `Canonical`
  ),
  show_percentage = F,
  fill_color = c("#5b9ecf", "#f59364", "#ddc55a"),
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_size = 0.5,
  stroke_alpha = 1,
  set_name_size = 3,
  text_size = 3) +
  scale_x_continuous(limits = c(-3,3)) +
  theme_void()

ggsave(filename = "f4a_Pt_classical_detection.pdf",
       p_venn_pt,
       width = 85, height = 60, units = "mm")

## mRNA-seq influence
p_venn_polya <- ggplot(venn_plot) +
  geom_venn(mapping = aes(
    A = `Osa_DSN_merge`,
    B = `Osa_Ribo-off_merge`,
    C = `Osa_PolyA_merge`
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

ggsave(p_venn_polya,file = "fs5b_PolyA_influence_venn.pdf",width = 9, height = 6, units = "cm")