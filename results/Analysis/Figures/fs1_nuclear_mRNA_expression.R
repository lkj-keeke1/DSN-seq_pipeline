##### featureCounts assignment #####
rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)
library(ggsci)

count_log <- as.data.frame(t(read.table("../../P3_nuclear_transcriptome/2.count/all.table.summary",
                                        row.names = 1,
                                        header = T)))

rownames(count_log) <- c("Osa_DSN_rep1", "Osa_DSN_rep2", "Osa_DSN_rep3",
                         "Osa_PolyA_rep1", "Osa_PolyA_rep2", "Osa_PolyA_rep3",
                         "Osa_Ribo-off_rep1", "Osa_Ribo-off_rep2", "Osa_Ribo-off_rep3")

count_plot <- rownames_to_column(count_log, var = "sample") %>% 
  gather(key = "assignment", value = "frag_num", -sample) %>% 
  filter(frag_num != 0)

p_featurecount <- ggplot() +
  geom_col(data = count_plot,
           mapping = aes(x = sample, y = frag_num, fill = assignment),
           position = "fill") +
  scale_fill_simpsons() +
  theme_bw() +
  labs(x = "",
       y = "ratio",
       fill = "Assignment") +
  coord_flip() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"))

ggsave(p_featurecount, file = "fs1a_featureCounts_assignment_evaluation.pdf",
       width = 150, height = 100, units = "mm")


##### nuclear expression correlation #####

rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)
library(ComplexHeatmap)
library(rtracklayer)

# input gtf info
gtf <- as.data.frame(
  rtracklayer::import("../../../data/ref/Osa/Osa_known_mRNA.gtf")) %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(chr = seqnames,
                transcript_id)

# input tpm matrix
all_exp <- read.table("../../P3_nuclear_transcriptome/2.count/all.tpm_matrix",
                      header = T, row.names = 1)


# delete 0, log2(tpm + 1)
all_exp_no0_log2 <- log2(all_exp[rowSums(all_exp) >0,] + 1)

colnames(all_exp_no0_log2) <- c("Osa_DSN_rep1", "Osa_DSN_rep2", "Osa_DSN_rep3",
                                "Osa_PolyA_rep1", "Osa_PolyA_rep2", "Osa_PolyA_rep3",
                                "Osa_Ribo-off_rep1", "Osa_Ribo-off_rep2", "Osa_Ribo-off_rep3")


# keep only nuclear mRNA
all_exp_no0_log2 <- rownames_to_column(all_exp_no0_log2, var = "transcript_id") %>%
  left_join(gtf, by = "transcript_id")
chr_matrix <- filter(all_exp_no0_log2, chr %in% 1:12) %>%
  select(-chr) %>%
  column_to_rownames(var = "transcript_id")


# correlation analysis
library(Hmisc)
chr_cor <- rcorr(as.matrix(chr_matrix))


col1 <- colorRampPalette(c("white", "white","white","white","white",
                           "white", "white","white","white","white",
                           "white", "white","blue",
                           "white", "red"))


# plot
pdf("fs1b_correlation_heatmap_chr.pdf", width = 5, height = 5)
corrplot::corrplot(corr = chr_cor$r,
                   method = "number",
                   col = col1(200),
                   diag = F,
                   outline = F,
                   order = "hclust",
                   #title = "Chr",
                   
                   # labels
                   tl.pos = "lt",
                   tl.cex = 0.7,
                   tl.col = "black",
                   
                   # color
                   # colors
                   cl.pos = "b",
                   cl.lim = c(0.7,1),
                   cl.length = 7,
                   cl.cex = 0.7,
                   cl.ratio = 0.15,
                   cl.align.text = "c",
                   
                   # number
                   number.cex = 0.00001
) +
  corrplot::corrplot(chr_cor$r, 
                     add = T,
                     method = "square", # 小格形式
                     type = "lower",   # 绘图位置
                     col = col1(200),
                     diag = F, # 是否显示中间为 1 的格子,
                     outline = T,
                     
                     #addgrid.col = "white",
                     #addCoef.col = "black",
                     addCoefasPercent = F,
                     order = "hclust",
                     
                     # label
                     tl.pos = "n",
                     # color
                     cl.pos = "n",
                     
                     # p-value
                     p.mat = chr_cor$P,
                     sig.level = c(0.001, 0.01, 0.05),
                     insig = "label_sig",
                     pch.cex = 0.7
  ) +
  corrplot::corrplot(corr = chr_cor$r,
                     add = T,
                     method = "number",
                     type = "upper",
                     col = "black",
                     diag = F,
                     outline = F,
                     order = "hclust",
                     
                     # labels
                     tl.pos = "n",
                     # color
                     cl.pos = "n",
                     # number
                     number.cex = 0.7
  )

dev.off()


