# scenv (win版本)
conda activate base
conda create -n scenv python=3.11 -y
conda activate scenv
pip install ipykernel
python -m ipykernel install --user --name=scenv --display-name "scenv"
pip install papermill ipython-autotime
pip install scanpy igraph leidenalg
pip install torch==2.7.1 torchvision==0.22.1 torchaudio==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pip install scvi-tools

# st
conda activate base
conda create -n st python=3.12 -y
conda activate st
pip install ipykernel
python -m ipykernel install --user --name=st --display-name "st"
pip install papermill
pip install scanpy igraph leidenalg
pip install squidpy
pip install stlearn






pip install datatable

conda install rpy2 -y


pip install openpyxl
pip install igraph
pip install leidenalg
pip install scikit-image
pip install harmonypy
pip install infercnvpy
pip install gseapy
pip install sccoda
pip install liana
pip install cellphonedb
pip install palantir
pip install scvelo



# scenv (linux版本)
conda activate base
conda create -n scenv python=3.10 -y
conda activate scenv
conda install r-base=4.4.1 -y
conda install scanpy -y
conda install rpy2 -y
pip install ipykernel
python -m ipykernel install --user --name=scenv --display-name "scenv"
pip install papermill
pip install ipython-autotime
pip install openpyxl
pip install igraph
pip install leidenalg
pip install scikit-image
pip install scvi-tools[cuda]
pip install harmonypy
pip install infercnvpy
pip install gseapy
pip install sccoda
pip install liana
pip install cellphonedb
pip install palantir
pip install scvelo
pip install datatable



















# R_4.4.1
conda activate base
conda create -n R_4.4.1 -y
conda activate R_4.4.1

conda search r-base
conda install r-base=4.4.1 -y
conda install r-pak r-remotes git gcc_linux-64 gxx_linux-64 gfortran_linux-64 make pkg-config -y

conda install r-devtools -y
conda install r-remotes -y
conda install r-biocManager -y
conda install r-tidyverse -y
conda install r-ggplot2 -y
conda install r-patchwork -y
conda install r-data.table -y
conda install r-openxlsx -y
conda install r-qs2 -y

conda install r-seurat -y
conda install r-seurat-data -y
conda install r-openxlsx -y
conda install r-qs2 -y


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

# scenv
$ conda install -c conda-forge scanpy python-igraph leidenalg




# scenv (win版本)
conda activate base
conda create -n scenv python=3.11 -y
conda activate scenv
conda install r-base=4.4.1 -y
conda install scanpy -y
conda install rpy2 -y
pip install ipykernel
python -m ipykernel install --user --name=scenv --display-name "scenv"
pip install papermill
pip install ipython-autotime
pip install openpyxl
pip install igraph
pip install leidenalg
pip install scikit-image
pip install harmonypy
pip install infercnvpy
pip install gseapy
pip install sccoda
pip install liana
pip install cellphonedb
pip install palantir
pip install scvelo
pip install datatable
pip install torch==2.7.1 torchvision==0.22.1 torchaudio==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pip install scvi-tools

# scenv (linux版本)
conda activate base
conda create -n scenv python=3.10 -y
conda activate scenv
conda install r-base=4.4.1 -y
conda install scanpy -y
conda install rpy2 -y
pip install ipykernel
python -m ipykernel install --user --name=scenv --display-name "scenv"
pip install papermill
pip install ipython-autotime
pip install openpyxl
pip install igraph
pip install leidenalg
pip install scikit-image
pip install scvi-tools[cuda]
pip install harmonypy
pip install infercnvpy
pip install gseapy
pip install sccoda
pip install liana
pip install cellphonedb
pip install palantir
pip install scvelo
pip install datatable
