rm(list = ls()); gc()
library(tidyverse)
library(ggplot2)
library(Seurat)
setwd("E:/research/bioinformation/singlecell/huage")
load("anno.scRNA_harmony.rdata")
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$orig.ident=="s1"),"group.3"]="c"
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$orig.ident=="s2"),"group.3"]="c"
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$orig.ident=="s3"),"group.3"]="t"
scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$orig.ident=="s4"),"group.3"]="t"
cd4.naive=subset(scRNA_harmony,idents ="Tcells")
seurat <- cd4.naive

deg <- FindMarkers(seurat, ident.1 = "c", ident.2 = "t", group.by = "group.3",
                   logfc.threshold = 0.1, min.pct = 0.1) %>% 
  rownames_to_column("gene") %>% 
  mutate(label = factor(case_when(avg_log2FC > log2(1.5) & p_val_adj < 0.05 ~ "Up",
                                  avg_log2FC < -log2(1.5) & p_val_adj < 0.05 ~ "Down",
                                  TRUE ~ "Stable"),
                        levels = c("Up", "Stable", "Down")))

deg_10 <- deg[sample(nrow(deg), 10), ]
ggplot(data = deg, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj), color = label)) +
  geom_hline(yintercept = -log10(0.05), lty=4, col="black", lwd=0.8) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), lty = 4, col = "black", lwd = 0.8) +
  scale_color_manual(values=c("red", "grey", "blue")) +
  labs(x = "log2 fold change", y = "-log10(p_val_adj)", color = "Lable") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, NA)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white"),) +
  geom_label_repel(data = deg_10, aes(label = gene))
  
  






