library(ggplot2) ##绘图
library(ggsignif) ##计算显著性差异
library(ggpubr) ##合并errorbar和显著性差异
df <- ToothGrowth

df$dose <- as.factor(df$dose)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

## 计算显著性
compare_means(len ~ supp, data = df, 
              group.by = "dose")


df2 <- data_summary(df, varname="len",
                    groupnames=c("supp", "dose"))

df2$dose=as.factor(df2$dose)

#为柱状图添加显著性标记先根据显著性分析的结果为每个柱状图添加醒醒的数目：

df2$label <- c("*","**","ns","*","**","ns")

ggplot(df2, aes(x=dose, y=len, fill=supp)) +
  ## 添加柱状图
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  ## 添加errorbar
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9))+
  ## 添加主标题
  labs(title="Effects of dose in different groups", x="Dose", y = "Results")+
  theme_classic() +
  ## 设置分组的颜色
  scale_fill_manual(values=c('#999999','black'))+
  ## 设置显著性标记的位置
  geom_text(
    aes(label = label, y= len+sd+0.5),
    position = position_dodge(0.9),
    vjust = 0,size = 6)+
  ## 对标题大小，字体代销，进行设置
  theme(plot.title = element_text(size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(colour="black", size=20),
        legend.text = element_text(colour="darkgrey", size = 16))

