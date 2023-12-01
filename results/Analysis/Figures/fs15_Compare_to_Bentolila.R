library(readxl)
library(tidyverse)
library(ggvenn)
library(patchwork)

## raw data comparison

us <- readxl::read_xlsx("../../../data/2023-12-01_Compare_to_Bentolila_raw_data/STS-PCRseq_Bentolila.xlsx", sheet = "sites in dataset 4", col_names = T) %>% 
  mutate(posi = paste(Gene , `Position in CDS`))
be1 <- readxl::read_xlsx("../../../data/2023-12-01_Compare_to_Bentolila_raw_data/STS-PCRseq_Bentolila.xlsx", sheet = "only CDS+10%", col_names = T) %>% 
  mutate(posi = paste(Gene , `position ATG`)) %>% 
  mutate(inus = posi %in% us$posi)

benot1 <- filter(be1, inus == F)

venn1 <- list(
  us = us$posi,
  be1 = be1$posi
)

ggvenn(data = venn1)
ggplot(data = be1) +
  geom_bar(mapping = aes(x = Gene, fill = inus)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## name changed in Bentolila
# Based on the results of last step, editing sites of some genes are totally missing
# It turns out to be the differences in gene ID of the two sheets.

be2 <- readxl::read_xlsx("../../../data/2023-12-01_Compare_to_Bentolila_raw_data/STS-PCRseq_Bentolila.xlsx", sheet = "step1-name-change", col_names = T) %>% 
  mutate(posi = paste(`Changed name` , `position ATG`)) %>% 
  mutate(inus = posi %in% us$posi)

benot2 <- filter(be2, inus == F)

venn2 <- list(
  `This paper` = us$posi,
  `Bentolila et al.` = be2$posi
)

p_venn2 <- ggvenn(data = venn2,
       show_percentage = F,
       fill_color = c("#9592e0", "#B3E092"),
       stroke_color = "white",
       stroke_size = 0.5,
       stroke_alpha = 1,
       set_name_size = 3,
       text_size = 3) +
  scale_x_continuous() +
  theme_void() +
  theme() +
  labs(title = "before position change")

ggplot(data = be2) +
  geom_bar(mapping = aes(x = Gene, fill = inus)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## position changed in Bentolila
# After change the names of some genes, editing sites of several genes are still totally missing
# It turns out to be the annotation differences between the reference genome used
# This was explained in the Figure S15 of the paper
be3 <- readxl::read_xlsx("../../../data/2023-12-01_Compare_to_Bentolila_raw_data/STS-PCRseq_Bentolila.xlsx", sheet = "step2-position-change", col_names = T) %>% 
  mutate(posi = paste(`Changed name` , `Changed position`)) %>% 
  mutate(inus = posi %in% us$posi)

benot3 <- filter(be3, inus == F)

venn3 <- list(
  `This paper` = us$posi,
  `Bentolila et al.` = be3$posi
)

p_venn3 <- ggvenn(data = venn3,
                  show_percentage = F,
                  fill_color = c("#9592e0", "#B3E092"),
                  stroke_color = "white",
                  stroke_size = 0.5,
                  stroke_alpha = 1,
                  set_name_size = 3,
                  text_size = 3) +
  scale_x_continuous() +
  theme_void() +
  theme() +
  labs(title = "after position change")

ggplot(data = be3) +
  geom_bar(mapping = aes(x = Gene, fill = inus)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

write.table(be3, file = "fs15_sites.txt",
            col.names = T, row.names = F,
            quote = F, sep = "\t")

write.table(benot3, file = "fs15_missing_sites.txt",
            col.names = T, row.names = F,
            quote = F, sep = "\t")


p_venn <- p_venn2 /p_venn3
ggsave(p_venn,filename = "fs15_compare_to_Bentolila.pdf",width = 6,height = 14,units = "cm")
