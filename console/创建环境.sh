# R_4.4.1
conda activate base
conda create -n R_4.4.1 -y
conda activate R_4.4.1

conda search r-base
conda install r-base=R_4.4.1 -y
which R

# 使用conda安装package
conda install r-devtools -y
onda install r-remotes -y
conda install r-BiocManager -y
conda install r-tidyverse -y
conda install r-ggplot2 -y
conda install r-patchwork
conda install r-data.table -y
conda install r-qs2 -y
conda install r-seurat -y
conda install r-seuratobject -y
conda install r-seurat-data -y
conda install r-hdf5r -y
conda install bioconductor-glmgampoi -y
conda install r-ape -y
conda install r-fields -y
conda install r-maps -y
conda install r-doubletfinder -y
conda install r-cowplot -y
conda install r-lattice -y
conda install r-matrix -y

# 使用R代码安装package
install.packages("Rfast2")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
