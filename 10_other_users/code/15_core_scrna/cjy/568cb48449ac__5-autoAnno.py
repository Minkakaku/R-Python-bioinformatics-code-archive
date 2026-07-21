import os
import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

# def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
#     plt.tight_layout()
#     try:
#         plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
#     except:
#         pass
#     plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
#     plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
#     plt.clf()

def annotate_clusters_by_panglao(adata, cluster_key='cluster', obsm_key='ora_estimate', min_n=3,out_dir=None):
    """
    使用 PanglaoDB marker 基于 cluster 平均表达执行 ORA 注释。

    参数：
    - adata: AnnData 对象，包含已聚类的结果。
    - cluster_key: 用于聚类的 obs 列（默认为 'cluster'）。
    - obsm_key: enrichment score 的 obsm 键名。
    - min_n: ORA 的最小 marker 数目。

    返回：
    - adata: 增加了 'cell_type' 注释的 AnnData。
    - annotation_dict: cluster → cell_type 的映射字典。
    """
    # 保证 cluster 是字符串
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

    # 获取 PanglaoDB markers 并去重
    markers = dc.get_resource('PanglaoDB')
    markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]

    # 聚合为 cluster 平均表达矩阵
    cluster_expr = adata.to_df().groupby(adata.obs[cluster_key]).mean()
    adata_cluster = sc.AnnData(cluster_expr)
    adata_cluster.var_names = cluster_expr.columns.astype(str)
    adata_cluster.obs_names = cluster_expr.index.astype(str)

    # 执行 ORA
    dc.run_ora(
        mat=adata_cluster,
        net=markers,
        source='cell_type',
        target='genesymbol',
        min_n=min_n,
        use_raw=False,
        verbose=True
    )

    # enrichment 分数矩阵
    acts = dc.get_acts(adata_cluster, obsm_key=obsm_key)
    ora_scores = acts.obsm['ora_estimate']

    # 构造成 DataFrame，行为 cluster，列为 cell_type
    acts_df = pd.DataFrame(
    ora_scores,
    index=acts.obs_names,
    columns=acts.var_names
    )

    # 自动注释（最大 enrichment cell type）
    annotation_dict = dc.assign_groups(acts_df)

    # 映射到每个细胞
    adata.obs['cell_type'] = adata.obs[cluster_key].map(annotation_dict)
    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(adata, color=["cluster","cell_type"], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.pdf'))
    return adata, annotation_dict

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Finde marker genes")
    parser.add_argument("--input_file", required=True, help="输入文件 h5ad")
    parser.add_argument("--out_dir", required=True, help="输出文件夹路径")
    parser.add_argument("--species", default="human", help="物种：human 或 mouse")
    args = parser.parse_args()
    autoAnno_out_dir = os.path.join(args.out_dir, '06-Cell_type_annotation')
    if not os.path.exists(autoAnno_out_dir):
        os.makedirs(autoAnno_out_dir)
    else:
        os.system(f'rm -rf {autoAnno_out_dir}')
        os.makedirs(autoAnno_out_dir)
    adata = sc.read_h5ad(args.input_file)
    if args.species != 'human':
        adata.var['gene'] = adata.var.index
        adata.var.index = [i.upper() for i in adata.var.index]

    markers = dc.get_resource('PanglaoDB') 
    adata_autoAnno, annotation_dict = annotate_clusters_by_panglao(adata, cluster_key='cluster', obsm_key='ora_estimate', min_n=3,out_dir=autoAnno_out_dir)
    adata_autoAnno.write_h5ad(os.path.join(autoAnno_out_dir, 'adata_autoAnno.h5ad'))
    anno_df = pd.DataFrame(list(annotation_dict.items()), columns=['cluster', 'cell_type'])
    anno_df.to_csv(os.path.join(autoAnno_out_dir, 'annotation_cluster.xls'), sep='\t', index=False)