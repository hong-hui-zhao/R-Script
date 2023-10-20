library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggsci)
library(reshape2)
library(clustree)
library(data.table)
library(scRNAtoolVis)
sce <- readRDS("~/program/Human LUAD/EPI/epi_anno.rds")

### marker 基因平均值热图
gene <- read.csv("~/program/Human LUAD/EPI/top5_epi.csv")
gene <- c( "S100A9","WFDC2", "CXCL14","NDRG1", "IGFBP3","PIGR",  "HPGD",  "CIT",  
           "STEAP4","SCGB3A1",  "MALAT1","MT-ND4L" , "NEAT1", "LUCAT1","ANKRD36C", "TM4SF4"  ,
           "PSCA"   ,  "S100P", "LYZ",   "TPPP3", "TMEM190" , "LRRIQ1","DNAAF1"  ,
           "CAPS"   ,  "HSPA1A" ,  "IGHG4"  ,  "SCGB1A1" , "CXCL1" ,  
           "CXCL8"   , "SCGB3A2" , "BPIFB1"  , "SFTPC"   , "LAMP3"  ,  "SFTPD"   , "SFTPA1"  , "SFTPA2" , 
           "C4BPA"   , "LTF"   ,   "C3"   ,    "SERPINA1" ,"HP"   ,    "SERPINE2", "FN1"   ,   "LGALS1" , 
           "TGFBI"   , "TIMP1"  ,  "PTPRC" ,   "CD69"  ,   "CCL5" ,    "CD52"   ,  "SRGN"   ,  "AGER"  ,  
          "EMP2"   ,  "MT2A"   ,  "MT1E"  ,   "MT1X"  )
# remove cluster anno name
b <- AverageHeatmap(object =sce,
               markerGene = gene,
               clusterAnnoName = F)
b <- plot(b)
ggsave("average_marker_epi.png",b, path = "~",width = 10,height = 10)


### marker基因在各个亚群的分布
# self replot
pc <- Embeddings(sce,reduction = 'umap') %>% data.frame()
# add gene name
exp <- FetchData(sce,vars = gene)
# megrge
mer <- cbind(pc,exp)
# wide to long
dflong <- melt(mer,id.vars = colnames(mer)[1:2],
               value.name = 'expressiton',
               variable.name = 'gene')


# filter quantile in 0.02~0.98
quantile(dflong$expressiton,probs = c(0.02,0.98))


# filter expression
dflong <- dflong %>% filter(expressiton >= 0.000000 & expressiton <= 6.109739)
                  
# plot
ggplot(dflong,aes(x = UMAP_1,y = UMAP_2)) +
  geom_point(aes(color = expressiton),size = 0.1,alpha = 1) +
  scale_color_gradient(low = 'grey90',high = '#CC0033',
                       name = 'percentile \n expression',
                       breaks = c(0.000000,6.109739),
                       limits = c(0.000000,6.109739),
                       labels = c('< 2%','> 98%')
  ) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1,
        legend.position = 'right',
        strip.text.x = element_text(size = 16,face = 'italic'),
        strip.background = element_rect(fill = NA,color = NA)) +
  facet_wrap(~gene) +
  labs(x = '',y = '') +
  guides(color = guide_colorbar(title.position = 'left',
                                title.theme = element_text(angle = 90,hjust = 0.5)))

                  