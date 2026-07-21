import os
import argparse
import pandas as pd
import json
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import scanpy as sc

def save_fig(plt, out_dir=os.getcwd(), file_name=None, fig_size={'w': 5, 'h': 5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches='tight', pad_inches=0.1, facecolor='white', dpi=300)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches='tight', pad_inches=0.1)
    plt.clf()

def plot_celltype_proportion(adata, out_dir, type_col='sample', fig_size=None):
    if type_col not in adata.obs.columns:
        print(f"Warning: column {type_col} not found in adata.obs. Skip proportion plot.")
        return
    # 获取细胞类型颜色映射
    keys = adata.obs['cell_type'].cat.categories.to_list()
    # 如果没有颜色则生成
    if 'cell_type_colors' not in adata.uns:
        sc.pl.umap(adata, color='cell_type', show=False)
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 计算比例表
    ctab = pd.crosstab(adata.obs[type_col], adata.obs['cell_type'], normalize='index')

    if fig_size is None:
        fig_width = max(len(ctab.index) * 1.5, 6)
        fig_size = (fig_width, 6)

    fig, ax = plt.subplots(figsize=fig_size)

    ctab.plot(kind='bar', stacked=True, ax=ax,
              color=[cell_type_colors[cell] for cell in ctab.columns],
              edgecolor='black', width=0.8)
    ax.grid(False)
    ax.set_xlabel(type_col.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    ax.tick_params(axis='x', labelsize=12, rotation=90)
    ax.tick_params(axis='y', labelsize=12)

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.,
        fontsize=10, title="Cell Types", title_fontsize=12,
        prop=legend_font, frameon=False,
        handleheight=1.5,
        handlelength=1.5,
        labelspacing=0.5,
        handletextpad=0.5,
        borderpad=0.8,
    )

    plt.subplots_adjust(right=0.75)
    file_name = f"celltype_proportion_{type_col}"
    save_fig(plt, out_dir, file_name=file_name, fig_size={'w': fig_size[0], 'h': fig_size[1]})


def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    print(f"Reading {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    # 1. 亚群注释
    if args.anno_file:
        print(f"Reading annotation file {args.anno_file}...")
        df_anno = pd.read_csv(args.anno_file, sep="\t", dtype=str)
        # 假设 anno_file 有两列：第一列是 cluster 编号，第二列是 cell_type
        cluster_col_name = df_anno.columns[0]
        celltype_col_name = df_anno.columns[1]
        anno_dict = dict(zip(df_anno[cluster_col_name], df_anno[celltype_col_name]))
        
        adata.obs['cell_type'] = [anno_dict.get(str(c), str(c)) for c in adata.obs[args.cluster_col]]
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    else:
        print(f"No anno_file provided. Using {args.cluster_col} as cell_type.")
        adata.obs['cell_type'] = adata.obs[args.cluster_col].astype('category')

    # 2. UMAP 展示 cell_type
    print("Plotting UMAP for cell_type...")
    sc.pl.umap(adata, color='cell_type', show=False, size=20, legend_loc='on data', legend_fontsize=10, legend_fontoutline=2)
    save_fig(plt, file_name='anno_umap', out_dir=args.outdir, fig_size={'w': 7, 'h': 5})
    
    # 3. Marker 基因可视化 (如果提供字典)
    if args.markers_dict:
        print(f"Reading markers dict {args.markers_dict}...")
        with open(args.markers_dict, 'r', encoding='utf-8') as f:
            markers = json.load(f)
        
        print("Plotting Marker Dotplot...")
        sc.pl.dotplot(adata, markers, groupby='cell_type', standard_scale='var', show=False, cmap='Reds')
        save_fig(plt, file_name='marker_list_bubble', out_dir=args.outdir, fig_size={'w': 10, 'h': 5})
        
        all_marker_genes = []
        for g_list in markers.values():
            all_marker_genes.extend(g_list)
        # 去重并保持顺序
        all_marker_genes = list(dict.fromkeys(all_marker_genes))
        valid_genes = [g for g in all_marker_genes if g in adata.var_names]
        
        if valid_genes:
            print("Plotting Marker UMAP...")
            sc.pl.umap(adata, color=valid_genes, cmap='jet', size=25, show=False)
            # 计算适合的画布大小
            n_genes = len(valid_genes)
            ncols = 4
            nrows = (n_genes - 1) // ncols + 1
            save_fig(plt, file_name='marker_list_umap', out_dir=args.outdir, fig_size={'w': ncols*4, 'h': nrows*4})

    # 4. Group 和 Sample 密度图
    for col in ['group', 'sample']:
        if col in adata.obs.columns:
            print(f"Plotting embedding density for {col}...")
            try:
                sc.tl.embedding_density(adata, basis='umap', groupby=col)
                sc.pl.embedding_density(adata, basis='umap', key=f'umap_density_{col}', show=False)
                save_fig(plt, out_dir=args.outdir, file_name=f'{col}_density', fig_size={'w': 20, 'h': 10})
            except Exception as e:
                print(f"Density plot for {col} failed: {e}")

    # 5. 统计与比例图
    if 'sample' in adata.obs.columns:
        result = adata.obs.groupby(['cell_type', 'sample']).size().reset_index(name='number').pivot(index='cell_type', columns='sample', values='number').reset_index()
        result.columns.name = None
        result.to_csv(os.path.join(args.outdir, 'celltype_statistics.xls'), sep='\t', index=None)
        
    print("Plotting celltype proportions...")
    plot_celltype_proportion(adata, out_dir=args.outdir, type_col='group')
    plot_celltype_proportion(adata, out_dir=args.outdir, type_col='sample')

    # 6. 按 Sample 拆分 UMAP
    if 'sample' in adata.obs.columns:
        samples = sorted(adata.obs['sample'].drop_duplicates().tolist())
        if len(samples) > 1:
            print("Plotting UMAP split by sample...")
            fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
            for i in range(len(samples)):
                sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", title=str(samples[i]), ax=axs[i] if len(samples)>1 else axs, wspace=0.1, legend_loc='on data' if i==len(samples)-1 else None, legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
            save_fig(plt, args.outdir, file_name='umap_split_sample_celltype', fig_size={'w': len(samples)*5, 'h': 5})

    # 7. 保存结果
    print("Saving results...")
    adata.obs.to_csv(os.path.join(args.outdir, 'all_cells_features.xls'), sep='\t')
    adata.write_h5ad(os.path.join(args.outdir, 'adata_anno.h5ad'))
    print("All done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Single-cell sub-cluster annotation and visualization")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--cluster_col", default="leiden", help="Column name of the cluster in adata.obs, default is 'leiden'")
    parser.add_argument("--anno_file", default=None, help="Tab-separated dictionary file for annotation (e.g. cluster_id \\t cell_type)")
    parser.add_argument("--markers_dict", default=None, help="JSON file containing marker dict for dotplot/umap (e.g. {'celltype': ['gene1', 'gene2']})")
    
    args = parser.parse_args()
    main(args)
