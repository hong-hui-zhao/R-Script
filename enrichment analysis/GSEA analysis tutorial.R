###  

# GSEA analysis tutorial --- R language

# author Honghui Zhao 2023/5/21

###
### Set the path
### load the R package for the analysis 
library(msigdbr)
library(clusterProfiler)
library(stringr)
library(tidyverse)
library(export)
library(readxl)
library(org.Hs.eg.db)
library(forcats)
library(ggstance)
library(GseaVis)
library(gggsea)
library(Hmisc)

### load the data for analysis
diff <- read.csv("",header = T)
colnames(diff)[1] <- 'symbol'

### GSEA --
set.seed(111) 
geneList <- diff$log2FoldChange
names(geneList) = diff$symbol
geneList = sort(geneList, decreasing = TRUE)
df <- msigdbr(species = "Homo sapiens")
df <- df %>% dplyr::select(gs_name, gene_symbol)

GSEA.res <- GSEA(geneList,minGSSize = 10,maxGSSize = 500,
                 pvalueCutoff = 0.15,pAdjustMethod = "BH",
                 verbose = FALSE,eps = 0,TERM2GENE = df)

GSEA.result <- GSEA.res@result
write.csv(GSEA.result ,file = "GSEA.csv")

### Visualization of GSEA results --
mygene <- c()
gseaNb(object = GSEA.res,geneSetID ='KEGG_ECM_RECEPTOR_INTERACTION',
       addPval = T,addGene = mygene,subPlot = 2)

### bacth plot
terms <- c('KEGG_VEGF_SIGNALING_PATHWAY',
           'KEGG_MAPK_SIGNALING_PATHWAY',
           'KEGG_NON_SMALL_CELL_LUNG_CANCER')

gseaNb(object = GSEA.res,geneSetID = terms,addPval = T,addGene = mygene)

# add segment line
dotplotGsea(data = GSEA.result,topn = 10,order.by = 'NES',add.seg = T)
ggsave('',width = 15,height = 12)

### fgsea------------------------------------------
library(fgsea)
library(data.table)

geneset = msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name,gene_symbol)
geneset$gs_name <- gsub('KEGG_','',geneset$gs_name)
geneset$gs_name <- tolower(geneset$gs_name)
geneset$gs_name <- gsub('_',' ',geneset$gs_name)

geneset$gs_name <- capitalize(geneset$gs_name)
GSEA_geneset <- geneset %>% split(x = .$gene_symbol, f = .$gs_name)

df <- diff[order(diff$avg_log2FC,decreasing = T),]
ranks <- df$avg_log2FC
names(ranks) <- df$gene

## analysis of fgsea
df <- fgsea(pathways = GSEA_geneset, stats = ranks,minSize=10,maxSize=500,eps=0.0)

### fgsea_KEGG Visualization of results
topPathwaysUp <- df[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- df[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GSEA_geneset[topPathways],ranks, df,gseaParam = 0.5)
#graph2ppt(file = 'GSEA-table.pptx',height = 7,width = 8.5)
fwrite(GSEA_df, file="fGSEA.csv")

### GSEA_KEGG --
genesymbol <- diff$symbol
entrezID <- bitr(genesymbol,fromType = "SYMBOL",
                 toType = "ENTREZID",OrgDb = "org.Hs.eg.db")# turn genesymbol into enetezID
genelist <- diff$log2FoldChange
names(genelist) <- diff$symbol
genelist <- genelist[names(genelist) %in% entrezID[,1]]
names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
genelist <- sort(genelist,decreasing = T)

KEGG <- gseKEGG(geneList = genelist,organism = "hsa",#https://www.genome.jp/kegg/catalog/org_list.html
                minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.15,
                pAdjustMethod = "BH",verbose = FALSE,eps = 0)
KEGG_result = setReadable(KEGG,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
KEGG_results <- KEGG_result@result
write.csv(KEGG_results,file = 'KEGG_results.csv')

### Visualization of KEGG results
dotplotGsea(data = KEGG,topn = 10,order.by = 'NES', add.seg = T)

gseaNb(object = KEGG,geneSetID ='hsa04010',addPval = T)

### GO analysis by clusterProfiler ----
GO <- gseGO(geneList = genelist,OrgDb = org.Hs.eg.db,
            ont = "ALL",minGSSize = 100,maxGSSize = 500,
            pvalueCutoff = 0.15,verbose = FALSE,seed = FALSE, by = "fgsea")

GO <- setReadable(GO,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
GO_result <- GO@result
GO@result<-GO[order(GO$NES,decreasing=T)]
# add segment line
dotplotGsea(data = GO,topn = 10,order.by = 'NES',add.seg = T)

gseaNb(object = GO,geneSetID ='GO:0001227',addPval = T)