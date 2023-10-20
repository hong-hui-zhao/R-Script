library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

anno_col <- brewer.pal(3, "Paired")
names(anno_col) <- paste('clsuter ',0:2,sep = '')

# get cells mean gene expression
mean_gene_exp <- AverageExpression(sce,
                                   features = mainmarkers,
                                   group.by = 'main_cell_type',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

# add colnames
colnames(mean_gene_exp) <- c("Immune","Epithelial","Stromal")

# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))

# color
col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))

# top annotation
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = anno_col))

Heatmap(htdf,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        column_title = "main_cell_type",
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_side = 'top',
        column_names_rot = 0,
        top_annotation = column_ha,
        col = col_fun,
        heatmap_legend_param = list(
          title = "Z-score", at = c(-2,1, 0,-1, -2)))
