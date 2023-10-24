library(Seurat)
pbmc <- readRDS("pbmc.rds")
table(pbmc$cell_type)
av <-AverageExpression(pbmc,
                       group.by = "cell_type",
                       assays = "RNA")
av=av[[1]]
head(av)

#选出标准差最大的1000个基因
cg=names(tail(sort(apply(av, 1, sd)),1000))
#查看这1000个基因在各细胞群中的表达矩阵
View(av[cg,])
#查看细胞群的相关性矩阵
View(cor(av[cg,],method = 'spearman'))
#pheatmap绘制热图
pheatmap::pheatmap(cor(av[cg,],method = 'spearman')) #默认是Pearson