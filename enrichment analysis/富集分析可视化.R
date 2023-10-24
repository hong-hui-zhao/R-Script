### load R package ----
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(pathview)
library(igraph)
library(cowplot)

### GO ----
ego <- pairwise_termsim(ego_ALL)
p1 <- emapplot(ego, showCategory = 20)

## cnetplot将基因和生物学概念(例如GO terms或KEGG pathways)之间的联系描述为一个网络
p2 <- cnetplot(ego_ALL,
               foldChange = genelist, 
               showCategory = 3,
               node_label = "all", # category | gene | all | none
               circular = TRUE, 
               colorEdge = TRUE)
ggsave('GO-pathway-gene_top5_pesu.png',p2 ,path = path2,width = 10,height=10)

p3 <- cnetplot(ego_ALL, showCategory = 5, circular = TRUE, colorEdge = TRUE)
ggsave('GP-pathway-gene_pesu.png',p3,path=path1,width = 15,height = 10)

terms <- ego_ALL$Description[1:10]
p4 <- pmcplot(terms, 2010:2023)
p5 <- pmcplot(terms, 2010:2023, proportion=FALSE)
plot_grid(p5, p4, ncol=2)
ggsave('go_top10_近10年研究趋势.png',path = path2,width = 40,height = 10 )

hsa04110 <- pathview(gene.data  = genelist,
                     pathway.id = "hsa04110",
                     species    = "hsa")

p6 <- dotplot(ego_ALL, showCategory=20) + ggtitle("dotplot for GO")
p7 <- barplot(ego_ALL, showCategory=20)
plot_grid(p6, p7, ncol=2)
ggsave('bar+dot_pesu.png',path = path1,width = 10,height = 10)

goplot(ego_BP)
ggsave('GO_从属关系_1_BP.png',path = path2,width = 8,height = 8)

goplot(ego_CC)
ggsave('GO_从属关系_1_CC.png',path = path2,width = 8,height = 8)

goplot(ego_MF)
ggsave('GO_从属关系_1_MF.png',path = path2,width = 8,height = 8)

ego_all <- data.frame(ego_ALL)
ego_all <- separate(data=ego_all, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
ego_all <- separate(data=ego_all, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
ego_all <- mutate(ego_all, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 

eGoBP <- ego_all %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- ego_all %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- ego_all %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)

p <- ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
p + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')
ggsave("GO_MF_CC_BP_TOP10.png",path = path1, width = 10,height = 10)

### KEGG 可视化
### 初步可视化
kegg_1 <- pairwise_termsim(KEGG)
p8 <- emapplot(kegg_1, showCategory = 20)
ggsave('kegg_path_path_pesu.png',p8,width = 8,height = 8,path = path2)

# cnetplot将基因和生物学概念(例如GO terms或KEGG pathways)之间的联系描述为一个网络
p9 <- cnetplot(KEGG_ges_result, categorySize="pvalue", foldChange=genelist)
ggsave('kegg_path_gene_pesu.png',p9,width = 8,height = 8,path = path2)

p10 <- cnetplot(KEGG_ges_result, foldChange=genelist, circular = TRUE, colorEdge = TRUE)
ggsave('kegg_path_gene_cir_pesu.png',path = path1,width = 8,height = 8,p10)

p11 <- dotplot(KEGG_ges_result, showCategory=20) + ggtitle("dotplot for kegg")
p12 <- barplot(KEGG_ges_result, showCategory=20)
plot_grid(p11, p12, ncol=2)
terms <- KEGG$Description[1:10]
p11 <- pmcplot(terms, 2010:2023)
p12 <- pmcplot(terms, 2010:2023, proportion=FALSE)
p13 <- plot_grid(p11 , p12, ncol=2)
ggsave('kegg_近10年变化.png',p13,width = 16,height = 8,path = path1)

hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))




