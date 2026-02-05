# 创建环境st
conda activate base
conda create -n st python=3.10.14 -y
conda activate st
conda install scanpy -y
pip install stlearn

# 创建环境cell2location
conda activate base
conda create -n cell2location python=3.10.14 -y
conda activate cell2location
conda install scanpy -y
