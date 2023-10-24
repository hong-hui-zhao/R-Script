library(ggpubr)
library(rstatix)
library(tidyverse)
library(ggsci)
library(ggplot2)
library(Rmisc)
data(ToothGrowth)

## stat_summary+stat_compare_means差异比较 --------
ggplot(ToothGrowth, aes(x = factor(dose), y = len, fill = supp)) + 
  stat_summary(fun = mean, geom = "bar", width = 0.6, 
               position = position_dodge(0.6)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.2, width = 0.2, 
               position = position_dodge(0.6)) + 
  scale_fill_aaas() + 
  labs(x = "dose", y = "len") + 
  theme_bw() +
  stat_compare_means(aes(group = supp), method = "wilcox.test", label = "p.signif")

## summarySE求均值+compare_means差异比较------------------
myt_test = compare_means(len ~ supp, data = ToothGrowth, group.by = "dose", method = "wilcox.test")

ToothGrowth %>% 
  summarySE(., measurevar = "len", groupvars = c("dose", "supp")) %>% 
  ggplot(ToothGrowth, mapping = aes(x = factor(dose), y = len)) + 
  geom_bar(aes(fill = supp), stat = "identity", position = position_dodge(0.6), width = 0.6) + 
  geom_errorbar(aes(ymin = len - se, ymax = len + se, group = supp), position = position_dodge(0.6), size = 0.2, width = 0.2) + 
  scale_fill_aaas() + 
  labs(x = "dose", y = "len") + 
  theme_bw() + 
  geom_text(data = myt_test, mapping = aes(x = factor(dose), y = 30, label = p.signif), 
            position = position_dodge(0.6), size = 5, fontface = "bold")

##  summarySE求均值+rstatix差异比较--------------

myt_test <- ToothGrowth %>% group_by(dose) %>% 
  wilcox_test(., len ~ supp) %>% # 获取 p 值
  mutate(., p = round(.$p, 4)) %>% # 保留4位
  add_significance(., 'p') %>% # 根据 p 值添加显著性标记 * 符号
  add_xy_position(., x = 'dose', fun = "max", dodge = 0.6)

## ggbarplot求均值+stat_compare_means差异比较


ggbarplot(ToothGrowth, x = 'dose', y = 'len', 
          fill = 'supp', add = 'mean_se', color = 'black', 
          position = position_dodge(0.6), width = 0.6, size = 0.5, legend = 'top') +
  scale_fill_aaas() +
  labs(x = "dose", y = "len") + 
  theme_bw() + 
  stat_compare_means(aes(group = supp), method = "wilcox.test", label = "p.signif")

## 箱线图+stat_compare_means差异比较

ggplot(ToothGrowth, aes(x = factor(dose), y = len, fill = supp)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_aaas() +
  labs(x = "dose", y = "len") + 
  theme_bw()+
  stat_compare_means(aes(group = supp), method = "wilcox.test", label = "p.signif")
                     



ToothGrowth %>% 
  summarySE(., measurevar = "len", groupvars = c("dose", "supp")) %>% 
  ggplot(ToothGrowth, mapping = aes(x = factor(dose), y = len)) + 
  geom_bar(aes(fill = supp), stat = "identity", position = position_dodge(0.6), width = 0.6) + 
  geom_errorbar(aes(ymin = len - se, ymax = len + se, group = supp), position = position_dodge(0.6), size = 0.2, width = 0.2) + 
  scale_fill_aaas() + 
  labs(x = "dose", y = "len") + 
  theme_bw() + 
  geom_text(data = myt_test, mapping = aes(x = factor(dose), y = y.position, label = p.signif), 
            position = position_dodge(0.6), size = 4, fontface = "bold")