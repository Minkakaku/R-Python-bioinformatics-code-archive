import anndata as ad
import diopy
import argparse
import warnings
import os
import scanpy as sc

warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(prog='trans.py', description="h5ad2h52rds")
parser.add_argument('-i', action='store', dest='h5ad', required=True, help='adata.h5ad')
parser.add_argument('-o', action='store', dest='out_dir',required=True, help='Result output directory./path/to/your/outdir/')
parser.add_argument('-n', action='store', dest='name',required=True, help='Result output h5name.')
parser.add_argument('-f', action='store', dest='RDSname',required=True, help='Result output RDSname.')
parser = parser.parse_args()


adata = ad.read_h5ad(parser.h5ad)
os.chdir(parser.out_dir)

adata.obs_names_make_unique()
# 这里有时候会报错:标准分析的h5ad中结构完整，有layers['counts],adata.raw,adata.X 但是recluster的h5ad中没有layers,adata.raw

# 此时结果中 raw(log1p)=seurat.obj['counts','data'],adata.X=seurat.obj['scale.data']
diopy.output.write_h5(adata,file = os.path.join(parser.out_dir,parser.name),save_X=False,save_graph = True)

os.system(f'sh /mnt/NexData_file/home/xzx/scripts/diopy/docker_run.sh "{os.path.join(parser.out_dir,parser.name)}" "{parser.out_dir}" "{parser.RDSname}"')

## 示例 
## conda activate /mnt/NexData_file/home/xzx/envs/scDIOR/
## python /mnt/NexData_file/home/xzx/scripts/diopy/trans_counts.py -i /mnt/NexData_file/home/xzx/project/SC2411133007-07/shouhou1/infer_adata.h5ad -o /mnt/NexData_file/home/xzx/project/SC2411133007-07/shouhou1/ -n "infer_adata.h5" -f "infer_adata.rds"