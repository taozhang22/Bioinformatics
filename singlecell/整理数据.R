###################################################################################
# 方法一：多个10X文件变成三个文件的
###################################################################################
rm(list = ls()); gc()
library(Seurat)
library(tidyverse)
library(qs2)
setwd("/home/yang/research/bioinformation/singlecell/database")

dir <- "GSE178318"
seurat <- ReadMtx(
  mtx = paste0(dir, "/GSE178318_matrix.mtx.gz"),
  features = paste0(dir, "/GSE178318_genes.tsv.gz"),
  cells = paste0(dir, "/GSE178318_barcodes.tsv.gz")
)
seurat <- CreateSeuratObject(counts = seurat)
seurat$Sample <- word(rownames(seurat@meta.data), start = 2, end = 3, sep = "_")
seurat <- split(seurat, seurat$Sample)
qs_save(seurat, paste0(dir, "/R/",dir, ".qs2"))
