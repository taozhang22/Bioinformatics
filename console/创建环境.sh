# R_4.4.1
conda activate base
conda create -n R_4.4.1 -y
conda activate R_4.4.1

conda search r-base
conda install r-base=4.4.1 -y
which R

# 使用conda安装package
install.packages("devtools")
install.packages("remotes")
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("data.table")
install.packages("qs2")
install.packages("Seurat")
install.packages("SeuratObject")
install.packages("SeuratData")
install.packages("hdf5r")
install.packages("ape")
install.packages("fields")
install.packages("maps")
install.packages("DoubletFinder")
install.packages("cowplot")
install.packages("lattice")
install.packages("Matrix")
BiocManager::install("glmGamPoi")
install.packages("Rfast2")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
