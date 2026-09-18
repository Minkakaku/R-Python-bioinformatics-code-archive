#定义函数
import sys
sys.path.append('/mnt/NexData_file/home/yrc/scrna-code/scanpy_fun')
import os
import itertools
from find_marker import find_marker_gene

def pairwise_find_marker_genes(
    adata,
    group_key: str,
    species: str,
    outdir_base: str,
    exclude_groups: list = None
):
    """
    对 `group_key` 中的所有组进行两两差异表达分析
    
    参数:
        adata: AnnData 对象
        group_key: `.obs` 中的分组列名称
        species: 物种名称 (human/mouse)
        outdir_base: 输出结果的根目录
        exclude_groups: 可选，排除某些 group 不参与比较
    """
    groups = adata.obs[group_key].unique().tolist()
    
    if exclude_groups:
        groups = [g for g in groups if g not in exclude_groups]
    
    # 所有两两组合
    pair_combinations = list(itertools.combinations(groups, 2))

    for group1, group2 in pair_combinations:
        comp_name = f"{group1}_vs_{group2}"
        outdir = os.path.join(outdir_base, comp_name)
        os.makedirs(outdir, exist_ok=True)

        # 保留这两个组的数据
        adata_subset = adata[adata.obs[group_key].isin([group1, group2])].copy()

        print(f"▶ 正在处理: {group1} vs {group2}，细胞数: {adata_subset.n_obs}")
        find_marker_gene(
            adata=adata_subset,
            group=group_key,
            species=species,
            out_dir=outdir
        )






#运行
adata = sc.read_h5ad('/mnt/NexData_file/home/yrc/FWSC20251610/10.22/recluster/anno.h5ad')
adata_Teff = adata[adata.obs['cell_type']=='Teff',:]
sc.pl.umap(adata_Teff, color=['cell_type'], legend_loc='on data',title='cell_type', add_outline=False,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
pairwise_find_marker_genes(
    adata_Teff,
    group_key = 'group',
    species = 'mouse',
    outdir_base = '/mnt/NexData_file/home/yrc/FWSC20251610/11.12/group_deg/Teff'
)