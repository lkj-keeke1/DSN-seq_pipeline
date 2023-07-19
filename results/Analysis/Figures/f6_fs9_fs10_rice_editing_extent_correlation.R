rm(list = ls())
options(stringsAsFactors = F)

load(file = "../RNA_editing/Osa/02.site_annotation/annotated_sites_filtered.Rdata")

# only new canonical sites

classical_sites_file <- "~/projects/DSN-seq/data/ref/Osa/Osa_new_canonical_sites.txt"

classcial_sites <- read.table(classical_sites_file, header = T, sep = "\t") %>% 
  mutate(chr = if_else(Organelle == "Mt", "_mt", "_pt"),
         edit_id = paste(Gene, chr, sep = ""),
         edit_id = paste(edit_id, `Position.in.CDS`, "CT", sep = " ")) %>% pull(edit_id)

all_sites_annotation <- all_sites_annotation[all_sites_annotation$edit_id %in% classcial_sites,] 
all_sites_annotation <- dplyr::rename(all_sites_annotation,
                               DSN_rep1 = Osa_DSN_rep1_extent,
                               DSN_rep2 = Osa_DSN_rep2_extent,
                               DSN_rep3 = Osa_DSN_rep3_extent,
                               DSN_all = Osa_DSN_merge_extent,
                               Ribo_rep1 = `Osa_Ribo-off_rep1_extent`,
                               Ribo_rep2 = `Osa_Ribo-off_rep2_extent`,
                               Ribo_rep3 = `Osa_Ribo-off_rep3_extent`,
                               Ribo_all = `Osa_Ribo-off_merge_extent`
                               )

selected_samples <- c('DSN_rep1',
                      'DSN_rep2',
                      'DSN_rep3',
                      'DSN_all',
                      'Ribo_rep1',
                      'Ribo_rep2',
                      'Ribo_rep3',
                      'Ribo_all')


#################### 方法内生物学重复间相关性分析 ##########################

# 相关性分析
library(Hmisc)
all_cor <- rcorr(as.matrix(all_sites_annotation[selected_samples]))


# 限定 xy 轴界限，方便进行文字 annotation
axis_min = 0
axis_max = 100


library(ggplot2)
library(patchwork)


##-------------------------------------------------------------------------
## DSN

### DSN_rep1 ~ DSN_rep2
lm_tmp <- round(lm(data = all_sites_annotation, DSN_rep2~DSN_rep1)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_rep1","DSN_rep2"],3)

DSN_inner1 <- ggplot(all_sites_annotation, aes(x= DSN_rep1, y=DSN_rep2)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_rep1",
       y = "Osa_DSN_rep2")



### DSN_rep1 ~ DSN_rep3
lm_tmp <- round(lm(data = all_sites_annotation, DSN_rep3~DSN_rep1)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_rep1","DSN_rep3"],3)

DSN_inner2 <- ggplot(all_sites_annotation, aes(x= DSN_rep1, y=DSN_rep3)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_rep1",
       y = "Osa_DSN_rep3")



### DSN_rep2 ~ DSN_rep3
lm_tmp <- round(lm(data = all_sites_annotation, DSN_rep3~DSN_rep2)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_rep2","DSN_rep3"],3)

DSN_inner3 <- ggplot(all_sites_annotation, aes(x= DSN_rep2, y=DSN_rep3)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_rep2",
       y = "Osa_DSN_rep3")

##-------------------------------------------------------------------------
## Ribo

### Ribo_rep1 ~ Ribo_rep2
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_rep2~Ribo_rep1)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["Ribo_rep1","Ribo_rep2"],3)

Ribo_inner1 <- ggplot(all_sites_annotation, aes(x= Ribo_rep1, y=Ribo_rep2)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_Ribo-off_rep1",
       y = "Osa_Ribo-off_rep2")



### Ribo_rep1 ~ Ribo_rep3
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_rep3~Ribo_rep1)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["Ribo_rep1","Ribo_rep3"],3)

Ribo_inner2 <- ggplot(all_sites_annotation, aes(x= Ribo_rep1, y=Ribo_rep3)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_Ribo-off_rep1",
       y = "Osa_Ribo-off_rep3")



### Ribo_rep2 ~ Ribo_rep3
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_rep3~Ribo_rep2)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["Ribo_rep2","Ribo_rep3"],3)

Ribo_inner3 <- ggplot(all_sites_annotation, aes(x= Ribo_rep2, y=Ribo_rep3)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_Ribo-off_rep2",
       y = "Osa_Ribo-off_rep3")


##------------------------------------------------------------------------
## DSN vs Ribo

### DSN_rep1 ~ Ribo_rep1
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_rep1~DSN_rep1)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_rep1","Ribo_rep1"],3)

DSN_Ribo_1 <- ggplot(all_sites_annotation, aes(x= DSN_rep1, y=Ribo_rep1)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_rep1",
       y = "Osa_Ribo-off_rep1")

### DSN_rep2 ~ Ribo_rep2
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_rep2~DSN_rep2)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_rep2","Ribo_rep2"],3)

DSN_Ribo_2 <- ggplot(all_sites_annotation, aes(x= DSN_rep2, y=Ribo_rep2)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_rep2",
       y = "Osa_Ribo-off_rep2")

### DSN_rep3 ~ Ribo_rep3
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_rep3~DSN_rep3)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_rep3","Ribo_rep3"],3)

DSN_Ribo_3 <- ggplot(all_sites_annotation, aes(x= DSN_rep3, y=Ribo_rep3)) +
  geom_point(size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  # 非线性拟合曲线
  #geom_smooth(method = "loess",color = "blue",
  #            size = 0.5) +  
  # slope = 1
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_rep3",
       y = "Osa_Ribo-off_rep3")


##### 

## DSN_merge vs Ribo_merge
lm_tmp <- round(lm(data = all_sites_annotation, Ribo_all~DSN_all)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["DSN_all","Ribo_all"],3)


#DSN_Ribo_mean
DSN_Ribo_merge <- 
  ggplot(all_sites_annotation, aes(x= DSN_all, y=Ribo_all)) +
  geom_point(aes(),size = 1, alpha = 0.7) +
  # 线性回归曲线
  geom_smooth(method = "lm",color = "red",
              size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "green",
              size = 0.5) +
  annotate("text", x = 75, y = 10, parse = F,
           label = paste("y=", lm_tmp[2], "x+", lm_tmp[1], sep = ""),
           size = 3) +
  annotate("text", x = 75, y = 3, parse = T,
           label = paste("R^2", "==", cor_tmp, sep = ""),
           size = 3) +  
  scale_x_continuous(limits = c(axis_min,axis_max)) +
  scale_y_continuous(limits = c(axis_min,axis_max)) +
  labs(x = "Osa_DSN_merge",
       y = "Osa_Ribo-off_merge",
       color = "Organelle")



##---------------------------------------------------------
## 汇总



theme <- theme_bw() +
  theme(text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"))

DSN_inner1 <- DSN_inner1 + theme
DSN_inner2 <- DSN_inner2 + theme
DSN_inner3 <- DSN_inner3 + theme
Ribo_inner1 <- Ribo_inner1 + theme
Ribo_inner2 <- Ribo_inner2 + theme
Ribo_inner3 <- Ribo_inner3 + theme
DSN_Ribo_1 <- DSN_Ribo_1 + theme
DSN_Ribo_2 <- DSN_Ribo_2 + theme
DSN_Ribo_3 <- DSN_Ribo_3 + theme

DSN_Ribo_merge <- DSN_Ribo_merge + theme

f6 <- DSN_inner1 + Ribo_inner1 + DSN_Ribo_merge + plot_annotation(tag_levels = "a")
f6_sup1 <- DSN_inner2 + DSN_inner3 + Ribo_inner2 + Ribo_inner3 + plot_layout(nrow = 2) 
f6_sup2 <- DSN_Ribo_1 + DSN_Ribo_2 + DSN_Ribo_3 

ggsave(f6, filename = "f6_rice_editing_extent_correlation.pdf",
       width = 180, height = 65, units = "mm")
ggsave(f6_sup1, filename = "fs9_rice_editing_extent_correlation.pdf",
       width = 120, height = 120, units = "mm")
ggsave(f6_sup2, filename = "fs10_rice_editing_extent_correlation.pdf",
       width = 180, height = 60, units = "mm")


# DSN_intra <- DSN_inner1 + DSN_inner2 + DSN_inner3
# Ribo_intra <- Ribo_inner1 + Ribo_inner2 + Ribo_inner3
# DSN_Ribo <- DSN_Ribo_1 + DSN_Ribo_2 +DSN_Ribo_3
# All <- DSN_intra / Ribo_intra / DSN_Ribo
# 
# ggsave(DSN_intra, filename = "fediting_extent_correlation_DSN.pdf", width = 18, height = 6 , units = "cm")
# ggsave(Ribo_intra, filename = "editing_extent_correlation_Ribo.pdf", width = 18, height = 6 , units = "cm")
# ggsave(DSN_Ribo, filename = "editing_extent_correlation_DSN-Ribo.pdf", width = 18, height = 6, units = "cm")
# 
# 
# ggsave(All, filename = "editing_extent_correlation_all.png", width = 18, height = 18 , units = "cm", dpi = 1000)
# ggsave(DSN_Ribo_merge, filename = "editing_extent_correlation_merge.pdf", width = 6, height = 6, units = "cm")

