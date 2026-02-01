%load_ext autotime
import time
start = time.time()

import warnings
warnings.filterwarnings("ignore")

import sys
from pathlib import Path
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scvi
os.chdir("/home/yang/research/bioinformation/singlecell/python/practice") # 需要根据实际的工作路径修改

print(f"Current python version: {sys.version}")
!conda list
pd.set_option("display.width", 1000)

dir = "result/base" # 根据需要修改，这个是你的sister下的目录，例如base、infercnv等
os.makedirs(dir, exist_ok=True)

######################################################################################################
# 读取文件
######################################################################################################
patterns = ["GSE132465", "GSE200997"] # 根据实际想要纳入的数据集名称修改
filenames = [f"../../database/{pattern}/python/{pattern}.h5ad" for pattern in patterns]
adata = sc.concat([sc.read_h5ad(f) for f in filenames], keys=patterns)

######################################################################################################
# qc
######################################################################################################
min_genes=100; min_cells=3; pct_counts_mt=25

# 计算质量控制指标
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)

# 绘制质控图片
fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
plt.tight_layout(); plt.show(); plt.close(fig)
    
fig, axes = plt.subplots(5, 1, figsize=(0.3 * adata.obs["Sample"].nunique(), 15))
for i, key in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]):
    sc.pl.violin(adata, keys=key, groupby="Sample", jitter=0.4, rotation=45, show=False, ax=axes[i])
plt.tight_layout(); plt.show(); plt.close(fig)

# 质控
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
adata = adata[adata.obs["pct_counts_mt"] < pct_counts_mt].copy()

# 去除双细胞
sc.pp.scrublet(adata, batch_key="Sample")
adata = adata[~adata.obs["predicted_doublet"]].copy()
adata.write_h5ad(f"{dir}/qc.h5ad")

######################################################################################################
# batch correction
######################################################################################################
# 查看数据是不是有批次效应
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata=adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata=adata, n_top_genes=3000, batch_key="Sample", subset=True)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata=adata, n_pcs=50, log=True)
sc.pp.neighbors(adata=adata, n_pcs=50)
sc.tl.umap(adata=adata)
sc.pl.umap(adata, color=["Sample", "Class"], title=["Sample", "Class"], legend_loc=None) # 绘制批次校正前的UMAP图
    
# scvi分析
scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata=adata, layer="counts", batch_key="Sample")
model = scvi.model.SCVI(adata=adata, n_layers=2, n_latent=30)
model.train(enable_progress_bar=False)
model.save(f"{dir}/scvi_model", overwrite=True)
    
# 查看批次效应去除后的效果
model = scvi.model.SCVI.load(f"{dir}/scvi_model", adata=adata)
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata=adata, use_rep="X_scVI")
sc.tl.umap(adata=adata)
sc.pl.umap(adata, color=["Sample", "Class"], title=["Sample", "Class"], legend_loc=None) # 绘制批次校正后的UMAP图

######################################################################################################
# celltype_annotation
######################################################################################################
# 参数，根据实际情况进行修改
res_list = np.arange(0.7, 0.8, 0.1)
markers = {
    "Epithelial": ["KRT18", "KRT8", "EPCAM"],
    "Endothelial": ["PLVAP", "VWF", "CLDN5"],
    "Fibroblast": ["COL1A1", "DCN", "COL3A1"],
    "Pericyte/SMC": ["RGS5", "PDGFRB", "TAGLN"],  # Smooth_muscle_cell
    "T_NK": ["CD3D", "CD3E", "KLRB1"],
    "B_cell": ["CD79A", "MS4A1", "BANK1"],
    "Plasma_B": ["MZB1", "DERL3", "IGHG2"],
    "McDC": ["LYZ", "CD68", "CD14"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C"],
    "Enteric_glial": ["PLP1", "CLU", "S100B"],
    "MAST": ["MS4A2", "KIT", "TPSAB1"],
    "Neutrophil": ["G0S2", "FCGR3B", "CSF3R"] 
}

for res in res_list:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res:.2f}")
    sc.pl.dotplot(adata, markers, groupby=f"leiden_{res:.2f}", standard_scale="var", title=f"Resolution {res:.2f}")
sc.pl.umap(adata=adata, color=[f"leiden_{res:.2f}" for res in res_list], legend_loc="on data")


# 参数，根据实际情况进行修改
my_res = 0.7
sc.tl.rank_genes_groups(adata, groupby=f"leiden_{my_res:.2f}", method="wilcoxon", pts=True)
deg = sc.get.rank_genes_groups_df(adata, group=None)
deg = deg.query("pvals_adj < 0.05 and logfoldchanges > 0.585 and pct_nz_group >= 0.25")
top20 = deg.groupby("group", as_index=False).head(20)
top20.to_csv(f"{dir}/top_deg.csv", index=False)


# 参数，根据实际情况进行修改
ct_map = {
    "0": "Epithelial",
    "1": "T_NK",
    "2": "T_NK",
    "3": "T_NK",
    "4": "B_cell",
    "5": "Epithelial",
    "6": "McDC",
    "7": "T_NK",
    "8": "T_NK",
    "9": "Plasma_B",
    "10": "Fibroblast",
    "11": "T_NK",
    "12": "Plasma_B",
    "13": "Epithelial",
    "14": "Endothelial",
    "15": "Plasma_B",
    "16": "Pericyte/SMC",
    "17": "Epithelial",
    "18": "Epithelial",
    "19": "B_cell",
    "20": "Enteric_glial",
    "21": "MAST",
    "22": "Epithelial",
    "23": "pDC",
    "24": "Neutrophil"
}

adata.obs["Celltype"] = adata.obs[f"leiden_{my_res:.2f}"].map(ct_map).astype("category") # 增加Celltype列
adata.obs["Celltype"] = adata.obs["Celltype"].cat.reorder_categories(list(markers.keys()), ordered=True) # 将Celltype类别的顺序设置为不区分大小写，按照字母顺序排序
celltypes = set(adata.obs["Celltype"].unique()) # 获取被注释到的细胞
filtered_genes = {ct: mk for ct, mk in markers.items() if ct in celltypes} # 保留被注释到的细胞和基因的字典

# 绘制注释后的气泡图
sc.pl.dotplot(adata, var_names=filtered_genes, groupby="Celltype", standard_scale="var", show=False)
plt.savefig(f"{dir}/dotplot_after_celltype_annotation.pdf", bbox_inches="tight")
plt.show(); plt.close()

# 绘制注释后的umap图
fig = sc.pl.umap(adata, color="Celltype", legend_loc="on data", show=False, return_fig=True, legend_fontsize=8)
fig.savefig(f"{dir}/umap_after_celltype_annotation.pdf", bbox_inches="tight")
plt.show(); plt.close(fig)

adata.write_h5ad(f"{dir}/celltype_annotated.h5ad")

######################################################################################################
# 获取qc之后带有注释的数据
######################################################################################################
obs_df = adata.obs["Celltype"].copy()
adata = sc.read_h5ad(f"{dir}/qc.h5ad")
adata.obs["Celltype"] = obs_df.reindex(adata.obs_names)
adata.write_h5ad(f"{dir}/adata.h5ad")
