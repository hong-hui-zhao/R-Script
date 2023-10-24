library(scRNAtoolVis)
library(tidyverse)
library(Seurat)
library(ggsci)
epi <- imm_anno

pc12 <- Embeddings(object = epi,reduction = 'umap') %>%
  data.frame()

# get botomn-left coord
lower <- floor(min(min(pc12$umap_1),min(pc12$umap_2)))

# get relative line length
linelen <- abs(0.3*lower) + lower

# mid point
mid <- abs(0.3*lower)/2 + lower

# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('umap2','umap1'),each = 2))

# axies label
label <- data.frame(lab = c('UMAP2','UMAP1'),angle = c(90,0),
                    x = c(lower ,mid),y = c(mid,lower))

library(tidydr)
p2 <- DimPlot(epi,group.by= 'tissue_type' ,reduction = 'umap',pt.size = 2,label = T) +
  theme_dr(xlength = 0.3,
           ylength = 0.3,
           arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  NoLegend() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

p1+p2+p3
ggsave('imm_patient.png',p3,width = 8,height = 8)
