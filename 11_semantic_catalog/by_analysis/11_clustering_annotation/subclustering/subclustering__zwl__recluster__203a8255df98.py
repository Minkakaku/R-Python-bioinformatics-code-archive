import os
# 强制所有计算库使用单线程
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import font_manager as fm
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

# 加载亚群细分的函数
import sys
sys.path.append('/mnt/NexData_file/home/lst/scrna-code/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion
# os.chdir('/lvdata/wzb/scRNA/FW2023-656/2024.01.02/')
import pandas as pd

adata = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252792/result.h5ad")
print(adata.obs['cell_type'].unique())

tdTomato_h = (adata.raw.to_adata()[:,'tdTomato'].X.toarray().flatten() > 0)|(adata.raw.to_adata()[:,'Kit'].X.toarray().flatten() > 0)
adata.obs['tdTomato_cat'] = 'other'  # 创建副本
adata.obs.loc[tdTomato_h, 'tdTomato_cat'] = 'tdTomato+'

con_TST = adata[adata.obs['sample'].isin(["tst", "control"])]
### 对con_TST 重分群
con_TST = recluster(con_TST,species='mouse',neighbors=20,resolution=0.5,n_pcs = 10,batch_effect=2,output="/mnt/NexData_file/home/xzx/project/FWSC20252792/con_tst")
con_TST = con_TST.sub_cluster()