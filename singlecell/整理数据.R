###################################################################################
# 方法一：多个10X文件变成三个文件的
###################################################################################
rm(list = ls()); gc()
library(qs2)
library(tidyverse)
library(patchwork)
library(Seurat)
setwd("/home/yang/research/bioinformation/singlecell/database") # 根据实际情况修改工作路径

dir <- "GSE178318" # 根据实际情况修改参数
seurat <- ReadMtx(
  mtx = paste0(dir, "/GSE178318_matrix.mtx.gz"), # 根据实际情况修改文件的名字
  features = paste0(dir, "/GSE178318_genes.tsv.gz"), # 根据实际情况修改文件的名字
  cells = paste0(dir, "/GSE178318_barcodes.tsv.gz") # 根据实际情况修改文件的名字
)
seurat <- CreateSeuratObject(counts = seurat, min.cells = 0, min.features = 0)
seurat$orig.ident <- word(rownames(seurat@meta.data), start = 2, end = 3, sep = "_")
seurat <- split(seurat, seurat$orig.ident)
qs_save(seurat, paste0(dir, "/R/",dir, ".qs2"))

###################################################################################
# 质量控制
###################################################################################
# 计算质量控制指标
seurat[["percent.mt"]]   <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RPS|^RPL")
# 绘制散点图
for (sample in unique(seurat$orig.ident)) {
  seu <- subset(seurat, subset = orig.ident == sample)
  p1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F)
  p2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) 
  p <- p1 + p2
  ggsave(filename = paste0(dir, "/R/", "FeatureScatter_", sample, ".png"), plot = p, width = 8, height = 4)
}
# 绘制小提琴图
for (label in c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")) {
  p <- VlnPlot(seurat, features = label, group.by = "orig.ident", pt.size = 0.01, layer = "counts") + NoLegend()
  ggsave(filename = paste0(dir, "/R/", "VlnPlot_",label, "_", sample, ".png"), plot = p,
    width = 0.5*length(unique(seurat$orig.ident)), height = 5)
}
