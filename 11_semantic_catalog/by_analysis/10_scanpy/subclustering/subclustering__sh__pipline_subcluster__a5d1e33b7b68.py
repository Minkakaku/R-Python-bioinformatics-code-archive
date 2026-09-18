import os, sys
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from matplotlib import font_manager as fm

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    '''
    out_dir : 输出文件夹地址
    file_name : 文件名称
    fig_size : 输出图片的大小
    '''
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None):
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    if fig_size is None:
        fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

# 映射所有细分 后面没映射
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2604012626/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-01/02.Mission/mapping_raw/"

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/01.Mission/Bcell/Bcell.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/01.Mission/Epithelial/Epithelial.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')


from scipy.io import mmwrite
import scipy.io
# scipy.io.mmwrite(os.path.join(out_dir,'expression_matrix.mtx'), adata_raw.X)
# # 同时输出基因名和细胞名
# pd.Series(adata_raw.var_names).to_csv(os.path.join(out_dir,'genes.csv'), index=False, header=False)
# pd.Series(adata_raw.obs_names).to_csv(os.path.join(out_dir,'cells.csv'), index=False, header=False)
# 获取表达矩阵
if hasattr(adata_raw.X, 'toarray'):
    expr_matrix = adata_raw.X.toarray()
else:
    expr_matrix = adata_raw.X

# 创建 DataFrame
df_expr = pd.DataFrame(
    expr_matrix,
    index=adata_raw.obs_names,
    columns=adata_raw.var_names
)

# 添加细胞元数据
meta_columns = ['cell_type', 'sub_cell_type', 'group', 'sample']  # 替换为你的列名
existing_meta = [col for col in meta_columns if col in adata_raw.obs.columns]

if existing_meta:
    for col in existing_meta:
        df_expr.insert(0, col, adata_raw.obs[col].values)

# 输出
df_expr.to_csv(os.path.join(out_dir, 'expression_matrix_with_metadata.csv'))
print(f"Saved: {df_expr.shape}")



# 获取所有细胞类型和样本
sub_cell_types = adata_raw.obs['sub_cell_type'].unique()
samples = adata_raw.obs['sample'].unique()

# 使用原始计数（如果存在）
if adata_raw.raw is not None:
    expr_matrix = adata_raw.raw.X
    var_names = adata_raw.raw.var_names
    print("Using raw counts")
else:
    expr_matrix = adata_raw.X
    var_names = adata_raw.var_names

# 转换为稠密矩阵
if hasattr(expr_matrix, 'toarray'):
    expr_matrix = expr_matrix.toarray()

# 创建完整 DataFrame
df_full = pd.DataFrame(
    expr_matrix,
    index=adata_raw.obs_names,
    columns=var_names
)

# 添加元数据
meta_columns = ['cell_type', 'sub_cell_type', 'group', 'sample']  # 替换为你的列名
existing_meta = [col for col in meta_columns if col in adata_raw.obs.columns]

if existing_meta:
    for col in existing_meta:
        df_full.insert(0, col, adata_raw.obs[col].values)

# 按细胞类型和样本分组输出
for ct in sub_cell_types:
    for sample in samples:
        mask = (df_full['sub_cell_type'] == ct) & (df_full['sample'] == sample)
        df_subset = df_full[mask]
        if len(df_subset) > 0:
            df_expr = df_subset.drop(columns=['cell_type', 'sub_cell_type', 'sample', 'group'])
            filename = f"{ct}_{sample}_expression.csv"
            filepath = os.path.join(out_dir, filename)
            df_expr.to_csv(filepath)
            print(f"Saved: {filename} ({len(df_subset)} cells)")