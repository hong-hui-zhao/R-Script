# single-cell analysis package
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("C:/Users/ZHH/Desktop/program/Human LUAD/EPI/epi_anno.rds")

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "tutorial" # hdWGCNA实验名称
)

# 构造metacell
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type_epi", "sample_id"),  # 将来自同一Sample的同种cell_type构造成metacell
  reduction = "harmony", # 用于执行KNN算法，此处harmony为通过harmony校正过的PCA
  k = 25, # KNN参数
  max_shared = 10, # 两个metacell可共享的最大细胞数
  ident.group = "cell_type_epi"  # Idents名
)
seurat_obj <- NormalizeMetacells(seurat_obj)
# 标准化metacell表达矩阵
seurat_obj@misc$tutorial$wgcna_metacell_obj
head(seurat_obj@misc$tutorial$wgcna_metacell_obj, 2)
table(seurat_obj@meta.data$metacell_grouping)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = 'AT2', # 挑选感兴趣的细胞类型，可以c("INH", "EX")挑选多个感兴趣的细胞类型
  group.by = "cell_type_epi", # INH所在的metadata列名
  assay = "RNA",
  slot = "data"
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
ggsave('hdwgcna-wrap-at2.png',width = 10,height = 10)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=9,
  setDatExpr=FALSE,
  tom_name = 'INH' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='AT2 hdWGCNA Dendrogram')
ggsave('AT2-hdWGCNA-Dendrogram.png',width = 10,height = 10)

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="sample_id"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type_epi', group_name = 'AT2'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "AT2-L"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p

ggsave('KME-AT2.png',p,width = 25,height = 45)
