rm(list = ls())
options(stringsAsFactors = F)

load(file = "../RNA_editing/Ath/02.site_annotation/annotated_sites_filtered.Rdata")

all_sites_annotation <- filter(all_sites_annotation, DSN_det == T)

#################### 方法内生物学重复间相关性分析 ##########################

# 相关性分析
library(Hmisc)
all_cor <- rcorr(as.matrix(all_sites_annotation[,c("Ath_DSN_rep1_extent", "Ath_DSN_rep2_extent", "Ath_DSN_rep3_extent")]))


# 限定 xy 轴界限，方便进行文字 annotation
axis_min = 0
axis_max = 100


library(ggplot2)
library(patchwork)


##-------------------------------------------------------------------------
## DSN

### Ath_DSN_rep1 ~ Ath_DSN_rep2
lm_tmp <- round(lm(data = all_sites_annotation, Ath_DSN_rep2_extent~Ath_DSN_rep1_extent)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["Ath_DSN_rep1_extent","Ath_DSN_rep2_extent"],3)

DSN_inner1 <- ggplot(all_sites_annotation, aes(x= Ath_DSN_rep1_extent, y=Ath_DSN_rep2_extent)) +
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
  labs(x = "Ath_DSN_rep1",
       y = "Ath_DSN_rep2")



### Ath_DSN_rep1_extent ~ Ath_DSN_rep3_extent
lm_tmp <- round(lm(data = all_sites_annotation, Ath_DSN_rep3_extent~Ath_DSN_rep1_extent)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["Ath_DSN_rep1_extent","Ath_DSN_rep3_extent"],3)

DSN_inner2 <- ggplot(all_sites_annotation, aes(x= Ath_DSN_rep1_extent, y=Ath_DSN_rep3_extent)) +
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
  labs(x = "Ath_DSN_rep1",
       y = "Ath_DSN_rep3")



### Ath_DSN_rep2 ~ Ath_DSN_rep3
lm_tmp <- round(lm(data = all_sites_annotation, Ath_DSN_rep3_extent~Ath_DSN_rep2_extent)$coefficients,3)
cor_tmp <- round(all_cor[["r"]]["Ath_DSN_rep2_extent","Ath_DSN_rep3_extent"],3)

DSN_inner3 <- ggplot(all_sites_annotation, aes(x= Ath_DSN_rep2_extent, y=Ath_DSN_rep3_extent)) +
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
  labs(x = "Ath_DSN_rep2",
       y = "Ath_DSN_rep3")



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



DSN_intra <- DSN_inner1 + DSN_inner2 + DSN_inner3


ggsave(DSN_intra, filename = "fs14_ath_editing_extent_correlation.pdf", width = 18, height = 6 , units = "cm")
