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
import sys
sys.path.append('/mnt/NexData_file/share/scRNA-tools/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion
import argparse
import pandas as pd

def write_run_parameters(args, outdir):
    """
    记录本次运行参数到输出目录
    """
    log_file = os.path.join(outdir, "run_parameters.txt")

    with open(log_file, "w", encoding="utf-8") as f:
        f.write("Subcluster Analysis Parameters\n")
        f.write("=" * 50 + "\n")
        #f.write(f"Run Time       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input File     : {args.input}\n")
        f.write(f"Output Folder  : {outdir}\n")
        f.write(f"Group Name     : {args.name}\n")
        f.write(f"Cell Types     : {', '.join(args.celltypes)}\n")
        f.write(f"Species        : {args.species}\n")
        f.write(f"Resolution     : {args.resolution}\n")
        f.write(f"Batch Effect   : {args.batch_effect}\n")
        f.write("=" * 50 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recluster subsets of an AnnData object.")
    parser.add_argument("-i", "--input", required=True, help="Input h5ad file path.")
    parser.add_argument("-o", "--outdir", required=True, help="Base output directory.")
    parser.add_argument("-n", "--name", required=True, help="Name of the sub-cluster group folder (e.g. NKT).")
    parser.add_argument("-c", "--celltypes", required=True, nargs='+', help="List of cell types to extract (e.g., 'T cells' 'NK cells').")
    parser.add_argument("-s", "--species", default="mouse", help="Species (human or mouse). Default: mouse.")
    parser.add_argument("-r", "--resolution", type=float, default=0.5, help="Clustering resolution. Default: 0.5.")
    parser.add_argument("-b", "--batch_effect", default=0, help="Batch effect column name, or 0 for none. Default: 0.")
    
    args = parser.parse_args()
    
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    outdir1 = os.path.join(args.outdir, args.name)
    os.makedirs(outdir1, exist_ok=True)
    
    write_run_parameters(args, outdir1)
    
    print(f"Subsetting data for cell types: {args.celltypes}")
    adata_tmp = adata[adata.obs['cell_type'].isin(args.celltypes)].copy()
    
    print("Plotting initial UMAP...")
    sc.pl.umap(adata_tmp, color=['cell_type', 'sample'], show=False)
    #save_fig(plt, out_dir=outdir1, file_name='pre_recluster_umap')
    
    print(f"Running recluster: resolution={args.resolution}, species={args.species}, batch_effect={args.batch_effect}")
    # format batch_effect back to int 0 if it is passing "0" as string from command line, otherwise keep string
    batch_effect_val = 0 if str(args.batch_effect) == "0" else args.batch_effect
    
    adata_tmp = recluster(adata_tmp, species=args.species, resolution=args.resolution, batch_effect=batch_effect_val, output=outdir1)
    adata_tmp = adata_tmp.sub_cluster()

    print("Reclustering complete.")
