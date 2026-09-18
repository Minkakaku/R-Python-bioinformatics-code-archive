import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
import scanpy as sc
import os
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# ========== 绘图风格（Nature 系列） ==========
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.6
plt.rcParams['ytick.major.width'] = 0.6
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['ytick.major.size'] = 3
sns.set_context("paper", font_scale=1.2, rc={"lines.linewidth": 0.8})

ZONE_PALETTE = {
    'MPC_inner': '#E64B35',
    'MPC_outer': '#F39C12',
    'TNBC_outer': '#4DBBD5'
}

# ================== 数据准备（区域独立占比） ==================
def prepare_proportions(adata, sample_col='sample', zone_col='cancer_type_zone', celltype_col='annotation'):
    """计算每个样本内每个区域各自细胞类型占比（分母=该样本该区域总细胞数）"""
    obs = adata.obs[[sample_col, zone_col, celltype_col]].copy()
    obs = obs.dropna(subset=[celltype_col])
    target_zones = ['MPC_inner', 'MPC_outer', 'TNBC_outer']
    obs = obs[obs[zone_col].isin(target_zones)]

    # 计数：样本-区域-细胞类型
    count = obs.groupby([sample_col, zone_col, celltype_col]).size().reset_index(name='count')
    # 样本-区域总细胞数
    total = count.groupby([sample_col, zone_col])['count'].sum().reset_index(name='total')
    prop = count.merge(total, on=[sample_col, zone_col])
    prop['proportion'] = prop['count'] / prop['total']
    prop = prop.rename(columns={celltype_col: 'cell_type', zone_col: 'zone', sample_col: 'sample'})
    return prop[['sample', 'zone', 'cell_type', 'proportion']]

# ================== 两两显著性（原始 p 值，不校正） ==================
def pairwise_significance_raw(prop_df):
    """Mann-Whitney U 检验，返回原始 p 值及显著性星号（未校正）"""
    zones = sorted(prop_df['zone'].unique())
    pairs = list(combinations(zones, 2))
    celltypes = prop_df['cell_type'].unique()
    result = {}
    for ct in celltypes:
        for g1, g2 in pairs:
            v1 = prop_df.loc[(prop_df['cell_type']==ct) & (prop_df['zone']==g1), 'proportion']
            v2 = prop_df.loc[(prop_df['cell_type']==ct) & (prop_df['zone']==g2), 'proportion']
            if len(v1) < 3 or len(v2) < 3:
                p = np.nan
            else:
                _, p = mannwhitneyu(v1, v2, alternative='two-sided')
            signif = 'ns'
            if not np.isnan(p):
                if p < 0.001: signif = '***'
                elif p < 0.01: signif = '**'
                elif p < 0.05: signif = '*'
            result[(ct, g1, g2)] = {'p_val': p, 'signif': signif}
    return result

# ================== 箱线图 + 显著性星号 ==================
def plot_boxplot_significance(prop_df, signif_dict, save_path, max_cols=5):
    """每个细胞类型一个子图，三组箱线图，显著对用横线+星号标注（原始p值）"""
    celltypes = sorted(prop_df['cell_type'].unique())
    n = len(celltypes)
    n_cols = min(max_cols, n)
    n_rows = int(np.ceil(n / n_cols))
    figsize = (2.8 * n_cols, 2.5 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    zones = ['MPC_inner', 'MPC_outer', 'TNBC_outer']

    for idx, ct in enumerate(celltypes):
        ax = axes[idx]
        sub = prop_df[prop_df['cell_type'] == ct]
        # 箱线图
        sns.boxplot(data=sub, x='zone', y='proportion', order=zones,
                    palette=ZONE_PALETTE, width=0.6, linewidth=0.8, ax=ax)
        sns.stripplot(data=sub, x='zone', y='proportion', order=zones,
                      color='black', alpha=0.5, size=3, jitter=True, ax=ax)
        ax.set_title(ct, fontsize=9, pad=4)
        ax.set_xlabel('')
        if idx % n_cols == 0:
            ax.set_ylabel('Proportion')

        # 添加显著性星号
        zone_to_x = {z: i for i, z in enumerate(zones)}
        y_max = sub['proportion'].max() if len(sub) else 0
        y_offset = max(y_max, 0.2) * 0.08
        drawn = set()
        for i, (g1, g2) in enumerate(combinations(zones, 2)):
            info = signif_dict.get((ct, g1, g2), {})
            if info.get('signif', 'ns') == 'ns':
                continue
            x1, x2 = zone_to_x[g1], zone_to_x[g2]
            if (x1, x2) in drawn or (x2, x1) in drawn:
                continue
            drawn.add((x1, x2))
            y = y_max + y_offset * (i + 1)
            ax.plot([x1, x1, x2, x2], [y, y + y_offset*0.3, y + y_offset*0.3, y],
                    lw=0.8, color='black')
            ax.text((x1 + x2) / 2, y + y_offset*0.3, info['signif'],
                    ha='center', va='bottom', fontsize=9)
        ax.set_ylim(-0.02, y_max + y_offset * (len(drawn) + 1.5))

    for j in range(idx + 1, len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout(pad=1.2)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig 1] 箱线图已保存: {save_path}")

# ================== 跨区域相关性热图 ==================
def compute_cross_correlation(prop_df, zone1, zone2):
    """计算两个区域间样本的细胞组成 Spearman 相关矩阵"""
    df1 = prop_df[prop_df['zone'] == zone1]
    df2 = prop_df[prop_df['zone'] == zone2]
    mat1 = df1.pivot_table(index='sample', columns='cell_type', values='proportion', fill_value=0)
    mat2 = df2.pivot_table(index='sample', columns='cell_type', values='proportion', fill_value=0)
    common_cts = sorted(set(mat1.columns) & set(mat2.columns))
    mat1 = mat1[common_cts]
    mat2 = mat2[common_cts]
    corr = pd.DataFrame(index=mat1.index, columns=mat2.index, dtype=float)
    for s1 in mat1.index:
        for s2 in mat2.index:
            r, _ = spearmanr(mat1.loc[s1], mat2.loc[s2])
            corr.loc[s1, s2] = r
    return corr

def plot_correlation_heatmaps(prop_df, save_path):
    """并排热图：MPC_inner vs TNBC_outer 和 MPC_outer vs TNBC_outer"""
    comps = [('MPC_inner', 'TNBC_outer'), ('MPC_outer', 'TNBC_outer')]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True)
    for ax, (z1, z2) in zip(axes, comps):
        corr_mat = compute_cross_correlation(prop_df, z1, z2)
        annot = True if len(corr_mat) <= 10 else False
        sns.heatmap(corr_mat.astype(float), ax=ax, cmap='RdYlBu_r', center=0.5, vmin=0, vmax=1,
                    annot=annot, fmt='.2f', linewidths=0.5, cbar_kws={'label': 'Spearman r'})
        ax.set_title(f'{z1} vs {z2}', fontsize=11, pad=10)
        ax.set_xlabel(z2, fontsize=9)
        ax.set_ylabel(z1, fontsize=9)
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig 2] 相关性热图已保存: {save_path}")

# ================== 主程序 ==================
def main():
    adata_path = 'final_all.h5ad'
    save_dir = './figures'
    os.makedirs(save_dir, exist_ok=True)

    adata = sc.read(adata_path)
    print(f"细胞总数: {adata.n_obs}")

    prop_df = prepare_proportions(adata)
    print(f"样本数: {prop_df['sample'].nunique()}, 细胞类型数: {prop_df['cell_type'].nunique()}")

    signif_dict = pairwise_significance_raw(prop_df)

    # 打印显著差异（p<0.05）
    sig_list = []
    for (ct, g1, g2), info in signif_dict.items():
        if info['signif'] != 'ns':
            sig_list.append((ct, f'{g1} vs {g2}', info['p_val'], info['signif']))
    if sig_list:
        sig_df = pd.DataFrame(sig_list, columns=['cell_type', 'comparison', 'p_value', 'signif'])
        print("\n原始p值显著的细胞类型 (p<0.05, 未校正):")
        print(sig_df.to_string(index=False))
    else:
        print("\n未发现p<0.05的差异（原始p值）")

    # 箱线图
    plot_boxplot_significance(prop_df, signif_dict,
                              save_path=os.path.join(save_dir, 'Fig1_boxplot_significance_raw.pdf'))
    # 热图
    plot_correlation_heatmaps(prop_df,
                              save_path=os.path.join(save_dir, 'Fig2_correlation_heatmaps.pdf'))

    print(f"\n所有图形已保存至: {save_dir}")

if __name__ == "__main__":
    main()
