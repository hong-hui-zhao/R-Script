# 读入数据
library(ggplot2)
library(ggsignif)
library(ggsci)
remotes::install_github("csdaw/ggprism")
library(ggprism)
library(tidyverse)
data <- openxlsx::read.xlsx("C:/Users/ZHH/Desktop/data.xlsx")

data$group <- stringr::str_split(data$group, "-", simplify = T)[, 1]
data$group <- factor(data$group, levels = c("KRAS Ctrl", "KRAS WT", "KRAS G12D"))

# 对expression列取对数
data$expression <- log(data$expression)

# Specify colors for each group
group_colors <- c("KRAS Ctrl" = "red", "KRAS WT" = "blue", "KRAS G12D" = "green")
data$group_color <- group_colors[data$group]

# 根据分组和基因求下表达量的均值和标准差
sum_data <- data %>%
  group_by(group) %>%
  summarise(across(starts_with("expression"), list(mean = mean, sd = sd)))

# 如果是两个分组用t.test函数就可以
## 定义差异分析函数，因为我这里是三个分组，组间两两比较，所以用下面的函数
func_HSD <- function(data, formula) {
  aov_re <- aov(formula, data)
  hsd_re <- TukeyHSD(aov_re)
  pvalue <- as.data.frame(hsd_re$group)
  return(pvalue)
}

formula <- expression ~ group

# 求出差异的p值，用于后面添加星号
pvalue <- func_HSD(data, formula)
pvalue$aster <- ifelse(pvalue$`p adj` <= 0.01, "**",
                       ifelse(pvalue$`p adj` <= 0.05, "*", "NS"))

# dir.create("img")
# openxlsx::write.xlsx(pvalue, file = file.path("img", paste0("pvalue_TP53.xlsx")), rowNames = T)
sum_data$group_color <- unique(data$group_color)
names(sum_data)
ggplot(sum_data, aes(group, expression_mean, fill = group_color)) +
  geom_bar(stat = "identity", width = 0.5) +
  # 添加误差棒
  geom_errorbar(aes(ymin = expression_mean - expression_sd, ymax = expression_mean + expression_sd), width = 0.2) +
  # x,y轴名
  xlab("") + ylab(paste0("Expression of KRAS")) +
  # stat_signif()和geom_signif()函数 
  # 添加差异是否显著星号，也可以添加p值
  geom_signif(comparisons = list(c("KRAS Ctrl", "KRAS WT"),
                                 c("KRAS Ctrl", "KRAS G12D"),
                                 c("KRAS WT", "KRAS G12D")),
              y_position = c(12, 14, 11),
              tip_length = 0.05, vjust = 0.2,
              annotation = pvalue$aster) +
  # 修改y轴区间
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3)) +
  # theme_prism修改主题为prism主题风格
  theme_prism(base_size = 14) + theme(legend.position = "none")
