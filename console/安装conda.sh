# 安装miniconda
cd ~/software
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_24.9.2-0-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ~/software/miniconda3
rm -f miniconda.sh
export PATH="$HOME/software/miniconda3/condabin:$HOME/software/miniconda3/bin:$PATH"
conda init
echo "$PATH" | tr ':' '\n'

#########################################################################################
# 接受协议（linux不需要）
#########################################################################################
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2

#########################################################################################
# 修改镜像
#########################################################################################
conda config --show channels
conda config --remove-key channels
conda config --add channels defaults
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --set channel_priority strict
conda config --show channels

pip config list
pip config unset global.index-url
pip config unset global.extra-index-url

pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
pip config set global.extra-index-url https://pypi.org/simple
pip config list

#########################################################################################
# 防止用conda在base中安装和修改
#########################################################################################
conda install -c conda-forge conda-protect -y
conda protect base
