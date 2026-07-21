import os
import shutil
import random
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
from PIL import Image
import celltypist
from celltypist import models
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from matplotlib import font_manager as fm
import scipy.io
import anndata
from scipy.sparse import csr_matrix
import warnings
sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA-rawdata/script/abridged_report/tmp/.scanpy_cached"
def sample_cells(adata, n_cells=2000):
    if adata.n_obs < n_cells :
        raise ValueError(f"n_cells ({n_cells}) is greater than the number of cells in the dataset ({adata.n_obs}).")
    else:
        indices = random.sample(range(adata.n_obs), n_cells)
        
    return adata[indices, :]

# 多平台读取矩阵
def get_10x_mtx(sample_path, sample_name, out_dir):
    path = os.path.join(sample_path,  'outs', 'filtered_feature_bc_matrix')
    data = sc.read_10x_mtx(path, var_names='gene_symbols', cache=False)
    adata = data
    #adata = sample_cells(data, n_cells=2000)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs['group'] = sample_name

    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    # adata.write_h5ad(os.path.join(out_dir,f'{sample_name}adata.h5ad'))
    #adata.write_h5ad(os.path.join(out_dir,'adata2000.h5ad'))

    return adata

def get_xunyin_mtx(sample_path, sample_name, out_dir):
    sorce_path=''
    if os.path.exists(os.path.join(sample_path, sample_name)):
        sorce_path=os.path.join(sample_path, sample_name)
    else:
        sorce_path=sample_path
    path = os.path.join(sorce_path, 'Outs', 'filtered_feature_bc_matrix')
    data = sc.read_10x_mtx(path, var_names='gene_symbols', cache=False)
    adata = data
    #adata = sample_cells(data, n_cells=2000)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs['group'] = sample_name

    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    #adata.write_h5ad(os.path.join(out_dir,'adata2000.h5ad'))

    return adata

def read_anndata_C4(path):
    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat ,dtype="float32")
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[0] .values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    adata.obs_names = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.var_names_make_unique()
    return adata

def get_C4_mtx(sample_path, sample_name, out_dir):

    path = os.path.join(sample_path,'output', 'filter_matrix')
    data = read_anndata_C4(path)
    adata = data
    #adata = sample_cells(data, n_cells=2000)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs['group'] = sample_name    

    #adata.write_h5ad(os.path.join(out_dir,'adata2000.h5ad'))
    return adata

def read_anndata_BD(path):
    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat ,dtype="float32")
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[0] .values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values

    barcode = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.obs_names = np.char.add('CellLabel_', barcode.astype(str))
    
    adata.var_names_make_unique()
    return adata

def get_BD_mtx(sample_path, sample_name, out_dir):

    path = os.path.join(sample_path,'filtered_feature_bc_matrix')
    data = read_anndata_BD(path)
    adata = data
    #adata = sample_cells(data, n_cells=2000)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs['group'] = sample_name

    #adata.write_h5ad(os.path.join(out_dir,'adata2000.h5ad'))
    return adata

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def set_scanpy_color_palette(adata, color_key, color_palette=None):
    """
    为 Scanpy 对象设置统一配色方案。
    
    参数:
    - adata: AnnData 对象
    - color_key: str，obs 中的列名，用于分组着色（如 'cell_type'）
    - color_palette: list，可选，自定义的颜色列表（最多 50 个）

    默认色板支持最多 50 种颜色，超出会报错。
    """
    default_palette = [
        "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
        "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
        "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
        "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
        "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
        "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
        "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
        "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
        "#f17d80", "#737495"
    ]
    palette = color_palette if color_palette else default_palette

    if color_key not in adata.obs:
        raise KeyError(f"'{color_key}' 不存在于 adata.obs 中")

    categories = adata.obs[color_key].astype(str).unique().tolist()
    n_categories = len(categories)

    if n_categories > len(palette):
        raise ValueError(f"类别数为 {n_categories}，超过了可用颜色数 {len(palette)}，请扩展色板。")

    # 设置为 Categorical 类型并指定顺序
    adata.obs[color_key] = adata.obs[color_key].astype("category")
    adata.obs[color_key] = adata.obs[color_key].cat.set_categories(categories)

    # 设置颜色
    adata.uns[f"{color_key}_colors"] = palette[:n_categories]

import matplotlib as mpl

mpl.rcParams.update({
    'legend.handleheight': 1.2,    # 全局色块高度
    'legend.handlelength': 1.8,    # 全局色块长度
    'legend.labelspacing': 0.4,     # 全局标签间距
    'legend.handletextpad': 0.4,    # 全局色块文本间距
    'legend.borderpad': 0.6,        # 全局内边距
})

# 如果使用seaborn，需要重置
import seaborn as sns
sns.reset_orig()  # 清除seaborn样式影响
def plot_celltype_proportion(_adata, out_dir, type='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = _adata.obs['cell_type'].cat.categories.to_list()
    values = list(_adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 计算比例表
    ctab = pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index')

    # 设置图像尺寸
    if fig_size is None:
        fig_width = max(len(ctab.index) * 1.5, 6)
        fig_size = (fig_width, 6)

    # 创建图像和子图
    fig, ax = plt.subplots(figsize=fig_size)

    # 绘制堆积条形图
    ctab.plot(kind='bar', stacked=True, ax=ax,
              color=[cell_type_colors[cell] for cell in ctab.columns],
              edgecolor='black', width=0.8)
    ax.grid(False)
    # 设置坐标轴
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    ax.tick_params(axis='x', labelsize=12, rotation=90)
    ax.tick_params(axis='y', labelsize=12)

    # 加粗边框
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

    # 设置图例
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.,
        fontsize=10, title="Cell Types", title_fontsize=12,
        prop=legend_font, frameon=False,
    handleheight=1.5,      # 增加色块高度（默认0.7）
    handlelength=1.5,      # 增加色块长度
    labelspacing=0.5,      # 增加标签间距
    handletextpad=0.5,     # 增加色块与文本间距
    borderpad=0.8,         # 增加图例内边距
    )

    # 调整子图布局，确保 legend 不遮挡主图
    plt.subplots_adjust(right=0.75)
    file_name=f"celltype_proportion_{type}"
    # 保存图像
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

def annotate_clusters_by_panglao(adata, cluster_key='cluster', organism='human', obsm_key='ora_estimate', min_n=3):
    """
    使用 PanglaoDB marker 基于 cluster 平均表达执行 ORA 注释。

    参数：
    - adata: AnnData 对象，包含已聚类的结果。
    - cluster_key: 用于聚类的 obs 列（默认为 'cluster'）。
    - organism: 'human' 或 'mouse'。
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

    return adata

def generate_report(platform,pid,psample,ptitle,paddress,pcustmor,out_dir):
    rmd_template_path = f'/mnt/NexData_file/share/scRNA-rawdata/script/abridged_report/rmd/{platform}_simple_report.rmd'
    # 读取并替换 RMD 模板
    with open(rmd_template_path, 'r', encoding='utf-8') as f:
        rmd_content = f.read()
        replacements = {
            '%pid%': pid,
            '%psample%': psample,
            '%ptitle%': ptitle,
            '%paddress%': paddress,
            '%pcustmor%': pcustmor}

        for key, val in replacements.items():
            rmd_content = rmd_content.replace(key, val)
    # 写入替换后的 RMD 文件
    rmd_output = os.path.join(out_dir, 'report.rmd')
    with open(rmd_output, 'w', encoding='utf-8') as w:
        w.write(rmd_content)

    os.system(f"/mnt/NexData_file/share/scRNA-rawdata/envs/scRNA_p/bin/R -e \"rmarkdown::render('{rmd_output}', output_file='{out_dir}/{pid}_{psample}_report.html')\"")

def analysis(adata, organism, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-') 
    adata.var['ribo'] = adata.var_names.str.upper().str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
    adata.layers["counts"] = adata.X.copy()
    # data normalization and log10 transform
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    # select highly variable gene and data scale 
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, key_added='cluster')
    sc.tl.rank_genes_groups(adata, 'cluster', method='wilcoxon', pts=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, f'cluster{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        tmp = tmp[(tmp['logfoldchanges'].abs()>1) & (tmp['pvals_adj']<0.05) & (tmp['pct_nz_group'] >= 0.25)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)
    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=3,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=1,
        cmap='bwr',
        colorbar_title='log fold change',
        show=False
    )
    save_fig(plt, out_dir, file_name='dotplot_top3_marker_gene', fig_size={'width':20, 'height':8})
    # =========================
    # copy 一份用于 celltype 注释
    # =========================
    adata_anno = adata.copy()

    # 保存原始基因名
    adata_anno.var['original_gene'] = adata_anno.var_names

    # 去掉第一个 "_" 前缀，并转大写
    adata_anno.var_names = [
        i.split('_', 1)[1].upper() if '_' in i else i.upper()
        for i in adata_anno.var_names
    ]

    # 处理重复基因名
    adata_anno.var_names_make_unique()

    # PanglaoDB 注释
    adata_anno = annotate_clusters_by_panglao(
        adata_anno,
        cluster_key='cluster',
        organism=organism
    )

    # 将注释结果映射回原始 adata
    adata.obs['cell_type'] = adata_anno.obs['cell_type'].values
    set_scanpy_color_palette(adata, color_key="cell_type")
    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(adata, color='cell_type', show=False, title='', legend_loc='on data', legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2)
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
        #plt.savefig(os.path.join(out_dir, 'umap_dimplot_celltype.pdf'), facecolor='white', bbox_inches="tight")
    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(adata, color='cluster', show=False, title='', legend_loc='on data', legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2)
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)
        #plt.savefig(os.path.join(out_dir, 'umap_dimplot_cluster.pdf'), facecolor='white', bbox_inches="tight")
    #adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')
    plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
    #adata.write(os.path.join(out_dir, 'result.h5ad'))
    result = adata.obs.groupby(['cell_type', 'sample']).size().reset_index(name='number').pivot(index='cell_type', columns='sample', values='number').reset_index()
    result.columns.name = None
    #result.to_csv(os.path.join(out_dir, 'celltype_statistics.xls'), sep='\t', index=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="scRNA 数据分析自动化流程")

    parser.add_argument("--platform", required=True, help="平台名称，例如 10X / C4 / BD / xunyin")
    parser.add_argument("--pid", required=True, help="项目 ID")
    parser.add_argument("--psample", required=True, help="样本 ID")
    parser.add_argument("--paddress", required=True, help="客户地址")
    parser.add_argument("--pcustmor", required=True, help="客户名称")
    parser.add_argument("--input_dir", required=True, help="输入数据目录")
    parser.add_argument("--sample_name", required=True, help="样本名称")
    parser.add_argument("--output_dir", required=True, help="输出目录")
    parser.add_argument("--organism", default="human", choices=["human", "mouse"], help="物种类型（默认 human）")

    args = parser.parse_args()

    # 根据平台读取数据
    if args.platform.lower() == "10x":
        adata = get_10x_mtx(args.input_dir, args.sample_name, args.output_dir)
    elif args.platform.lower() == "xunyin":
        adata = get_xunyin_mtx(args.input_dir, args.sample_name, args.output_dir)
    elif args.platform.lower() == "c4":
        adata = get_C4_mtx(args.input_dir, args.sample_name, args.output_dir)
    elif args.platform.lower() == "bd":
        adata = get_BD_mtx(args.input_dir, args.sample_name, args.output_dir)
    else:
        raise ValueError(f"未知平台: {args.platform}")
    analysis(adata, organism=args.organism, out_dir=args.output_dir)
    # 生成报告
    generate_report(
        platform=args.platform,
        pid=args.pid,
        psample=args.psample,
        ptitle=args.pid+'-'+args.psample,
        paddress=args.paddress,
        pcustmor=args.pcustmor,
        out_dir=args.output_dir
    )

    
    
    