import os
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
# import sys
# sys.path.append('/mnt/NexData_file/home/wmx/Script/scanpy_fun/')
# from subcluster import recluster
# from find_marker import find_marker_gene
# from draw_group import draw_group
# from draw_bar import plot_celltype_proportion

# import pandas as pd
# outdir = '/mnt/NexData_file/home/wmx/Project/SC2601080130/new/subcluster/'

# adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/Project/SC2601080130/subcluster/B_and_Plasma_delete/result.h5ad')
# celltype = 'B_and_Plasma_delete'
# npcs_l = [10,15,20,25,30,35,40,45,50]
# batch_effect_l = [0,1,2]
# neighbors_l = [10,15,20,25,30,35,40]

# for npcs in npcs_l:
#     for batch_effect in batch_effect_l:
#         for neighbors in neighbors_l:
#             outdir1 = os.path.join(outdir, f"{celltype}_npcs_{npcs}_batch_effect_{batch_effect}_neighbors_{neighbors}")
#             adata_tmp = adata
#             os.makedirs(outdir1,exist_ok=True)
#             adata_tmp = recluster(adata_tmp,species='mouse',resolution = 0.1,neighbors = neighbors,n_pcs = npcs,batch_effect=batch_effect,output=outdir1)
#             adata_tmp = adata_tmp.sub_cluster()

# npcs_l = [6,8]
# batch_effect_l = [0,1,2]
# neighbors_l = [10,15]

# for npcs in npcs_l:
#     for batch_effect in batch_effect_l:
#         for neighbors in neighbors_l:
#             outdir1 = os.path.join(outdir, f"{celltype}_npcs_{npcs}_batch_effect_{batch_effect}_neighbors_{neighbors}")
#             adata_tmp = adata
#             os.makedirs(outdir1,exist_ok=True)
#             adata_tmp = recluster(adata_tmp,species='mouse',resolution = 0.1,neighbors = neighbors,n_pcs = npcs,batch_effect=batch_effect,output=outdir1)
#             adata_tmp = adata_tmp.sub_cluster()

# 加载亚群细分的函数
import sys
sys.path.append('/mnt/NexData_file/share/scRNA-tools/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion

import pandas as pd
outdir = '/mnt/NexData_file/home/lst/project/SC2601080130/售后2/GC_PB'
os.makedirs(outdir,exist_ok=True)

adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/Project/SC2601080130/new/Anno/B_and_Plasma_delete_2/adata_anno.h5ad')
os.makedirs(outdir,exist_ok=True)
adata_tmp = adata[adata.obs['cell_type']=='PB']
sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
adata_tmp = recluster(adata_tmp,species='mouse',resolution = 0.5,neighbors = 20,n_pcs = 10,batch_effect=2,output=outdir)
adata_tmp = adata_tmp.sub_cluster()
