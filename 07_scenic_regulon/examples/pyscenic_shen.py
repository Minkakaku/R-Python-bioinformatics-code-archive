## 加载模块
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import loompy
import os


# 获取当前目录
os.getcwd()
os.listdir(os.getcwd())
# 1. 设置工作目录
os.chdir("/home/hf/project/shen/pyscenic")
OUTPUT_DIR = "/home/hf/project/shen/pyscenic/data"
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

## 读取示例文件
adata = sc.read_h5ad("../cellphonedb/mfOA.h5ad")
## 检查示例文件

adata_raw = adata.raw.to_adata()
rownames = dict(Gene=adata_raw.var_names.tolist())
colnames = dict(CellID=adata_raw.obs_names.tolist())
loompy.create(filename=os.path.join(OUTPUT_DIR, "mf.loom"), layers=adata_raw.X.transpose(), row_attrs=rownames, col_attrs=colnames)