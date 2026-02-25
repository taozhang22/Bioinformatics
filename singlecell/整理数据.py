###################################################################################################
# 代码开头
###################################################################################################
%load_ext autotime
import warnings
warnings.filterwarnings("ignore")

import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
os.chdir("/media/yang/9EBA46B1BA46862D/research/bioinformation/singlecell/database/") # 根据实际工作路径进行修改

###################################################################################################
# 方法一：表达矩阵和注释文件各为一个文件
###################################################################################################
# 以下参数根据实际情况进行修改
dir = "GSE132465"
file_count = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
file_meta1 = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"

# 制作adata
count = dt.fread(f"{dir}/{file_count}", sep="\t")
adata = count[:, 1:].to_pandas()
adata.index = count[:, 0].to_list()[0]
adata = adata.T

adata=sc.AnnData(adata)
adata.var_names_make_unique()

# 制作注释信息文件
meta1 = pd.read_csv(f"{dir}/{file_meta1}", index_col=0, sep="\t")
meta1 = meta1[["Sample", "Class"]]
meta2 = pd.read_csv(f"{dir}/{dir}.txt", sep="\t", index_col=0)

meta = meta1.join(meta2, on="Sample", how="inner")
meta["Source"] = dir
adata.obs = adata.obs.join(meta, how="left")
adata.write_h5ad(f"{dir}/python/{dir}.h5ad")

###################################################################################################
# 方法二: 存在barcodes、features和matrix，但是文件前面有前缀
###################################################################################################
dir = "GSE161277" # 根据实际情况进行修改

p = Path(f"{dir}/data")
seen = set()
for f in p.iterdir():
    name = f.name
    if "_barcodes" in name:
        seen.add(name.split("barcodes", 1)[0])
    elif "_features" in name:
        seen.add(name.split("features", 1)[0])
    elif "_matrix" in name:
        seen.add(name.split("matrix", 1)[0])
prefixes = sorted(seen)

adatas = []
for prefix in prefixes:
    adata = sc.read_10x_mtx(path=f"{dir}/data", var_names="gene_symbols", cache=False, prefix=prefix)
    adata.obs["Sample"] = "_".join(prefix.split("_")[1:3])
    adata.obs["Sample"] = adata.obs["Sample"].astype(str).str.replace("Patient", "P", regex=False)
    adata.obs["Sample"] = adata.obs["Sample"].astype(str).str.replace("normal", "N", regex=False)
    adata.obs["Sample"] = adata.obs["Sample"].astype(str).str.replace("carcinoma", "T", regex=False)
    adata.obs["Class"] = prefix.split("_")[2]
    adata.obs["Class"] = adata.obs["Class"].astype(str).str.replace("normal", "Normal", regex=False)
    adata.obs["Class"] = adata.obs["Class"].astype(str).str.replace("carcinoma", "Tumor", regex=False)
    adata.obs_names = adata.obs["Sample"].astype(str) + "_" + adata.obs_names.astype(str)
    adata.obs["Source"] = dir

    adata.var_names_make_unique()
    adatas.append(adata)
adata = sc.concat(adatas)
adata.var_names_make_unique()

meta = pd.read_csv(f"{dir}/{dir}.txt", sep="\t", index_col=0)
adata.obs = adata.obs.join(meta, on="Sample", how="left")

adata.write_h5ad(f"{dir}/python/{dir}.h5ad")

###################################################################################################
# 方法三：表达矩阵分为每个样本一个文件，注释信息只有一个文件
###################################################################################################
dir = "GSE166555" # 根据实际情况进行修改

filenames = os.listdir(f"{dir}/data")
adatas = []
for filename in filenames:
    count = dt.fread(f"{dir}/data/{filename}", sep="\t")
    adata = count[:, 1:].to_pandas()
    adata.index = count[:, 0].to_list()[0]
    adata = adata.T

    adata=sc.AnnData(adata)
    adata.var_names_make_unique()
    adatas.append(adata)
adata = sc.concat(adatas)
adata.var_names_make_unique()

meta1 = pd.read_csv(f"{dir}/GSE166555_meta_data.tsv.gz", sep="\t", index_col=0)
meta1["Sample"] = meta1.index.str.split(":").str[0]
meta1["Class"] = meta1["Sample"].str[-2:]
meta1.loc[meta1["Class"].str.contains("n"), "Class"] = "Normal"
meta1.loc[meta1["Class"].str.contains("t"), "Class"] = "Tumor"
meta1["Source"] = dir
meta1 = meta1.loc[:, ["Sample", "Class", "Source"]]
meta2 = pd.read_csv(f"{dir}/{dir}.txt", index_col=0, sep="\t")

meta = meta1.join(meta2, on="Sample", how="left")
adata.obs = adata.obs.join(meta, how="left")
adata.write_h5ad(f"{dir}/python/{dir}.h5ad")

###################################################################################################
# 质量控制
###################################################################################################
# 质量控制图
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)

for Sample in adata.obs["Sample"].astype(str).unique().tolist():
    fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))
    adat = adata[adata.obs['Sample'] == Sample].copy()
    sc.pl.scatter(adat, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
    sc.pl.scatter(adat, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
    plt.tight_layout(); plt.show(); plt.close(fig)

fig, axes = plt.subplots(5, 1, figsize=(0.3 * adata.obs["Sample"].nunique(), 15))
for i, key in enumerate(["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]):
    sc.pl.violin(adata, keys=key, groupby="Sample", jitter=0.4, rotation=45, show=False, ax=axes[i])
plt.tight_layout(); plt.show(); plt.close(fig)
print(adata.obs["Sample"].value_counts())
