library(tidyverse)

### 如果数据是单细胞稀疏矩阵 -----
exprSet <- as.matrix(tumor@assays[["RNA"]]@data)
exprSet<-as.data.frame(t(exprSet)) #转置
write.csv(exprSet,file = 'expr_tumor.csv')
exprSet <- read.csv("~/expr_tumor.csv")

### 
y <- as.numeric(exprSet[,"KRAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

# 导入数据
write.csv(data,file = 'tumor_cor.csv')
# 找出正相关性最强的前5个基因和负相关性最强的前5个基因
top5_pos <- head(data[order(-data$correlation),], 20)
top5_neg <- head(data[order(data$correlation),], 20)

# 绘制火山图
library(ggplot2)
library(ggrepel)
ggplot(data, aes(x = correlation, y = -log10(pvalue), color = correlation)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw() +
  labs(x = "Correlation", y = "-log10(P-value)", color = "Correlation", 
       title = "Correlation of KRAS with other genes") +
  geom_text_repel(data = top5_pos, aes(x = correlation, y = -log10(pvalue), label = symbol), 
                  size = 3, color = "red", hjust = 0, vjust = 0, 
                  segment.color = NA, point.padding = 1, box.padding = 0.5) +
  geom_text_repel(data = top5_neg, aes(x = correlation, y = -log10(pvalue), label = symbol), 
                  size = 3, color = "blue", hjust = 1, vjust = 0, 
                  segment.color = NA, point.padding = 1, box.padding = 0.5) +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(size = 0.5))

ggsave(filename = 'cor.pdf',width = 10,height = 10)
write.csv(cor_data_df,file = 'cor.csv')


cor_data_sig_pos <- cor_data_df %>%
  filter(pvalue < 0.01) %>% filter(correlation > 0)%>%
  arrange(desc(correlation))
cor_data_sig_neg <- cor_data_df %>%
  filter(pvalue < 0.01) %>% filter(correlation < 0)%>%
  arrange(desc(abs(correlation)))

library(ggstatsplot)
ggscatterstats(data = exprSet,
               y = KRAS,
               x = IGSF10,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram
               title = "Relationship between KRAS and COL17A1")
