degsall <- readxl::read_xlsx("C:/Users/ZHH/Desktop/DEGs_siHIF1A-As2_vs_cntrl.xlsx")
rm(list = ls())
library(tidyverse)
library(patchwork)
library(ggrepel)

log2FC = 0.5
padj = 0.05 


degsall$threshold="ns";
degsall[which(degsall$log2FoldChange  > log2FC & degsall$padj <padj),]$threshold="up";
degsall[which(degsall$log2FoldChange  < (-log2FC) & degsall$padj < padj),]$threshold="down";
degsall$threshold=factor(degsall$threshold, levels=c('down','ns','up'))
table(degsall$threshold)


### gene of EnhancedVolcano
ggplot(degsall, aes(x=log2FoldChange, y=-log2(padj), color = threshold)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c( "blue","grey","red") ) + 
  geom_label_repel(data=subset(degsall, log2FoldChange >= 1 & padj <= 8.020140e-32), 
                   aes(label=symbol),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.1, 
                   #max.overlaps = 200,
                   segment.size = 0.3,  #框的大小
                   size=4)+
  geom_label_repel(data=subset(degsall, log2FoldChange <= -1.744253  & padj <= 4.836504e-87), 
                   aes(label=symbol), label.padding = 0.1, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4)+
  theme_classic()
ggsave("火山图.png",width = 20, height = 16, units = "cm")
