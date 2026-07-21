import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl   # 新增，用于临时修改 rcParams
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')
import matplotlib.font_manager as font_manager   # ★ 必须导入

# 字体路径
font_path = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'

# 添加字体并全局启用
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'Arial'          # 首字母大写！
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


def calculate_sample_celltype_proportions(adata, sample_col='sample'):
    print("正在计算样本维度的细胞类型比例...")
    required_cols = [sample_col, 'MPCpercent', 'annotation']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"在adata.obs中未找到'{col}'列")
    df = adata.obs[[sample_col, 'MPCpercent', 'annotation']].copy()
    sample_cell_counts = pd.crosstab(df[sample_col], df['annotation'])
    sample_totals = sample_cell_counts.sum(axis=1)
    sample_proportions = (sample_cell_counts.div(sample_totals, axis=0) * 100).round(4)
    sample_groups = df[[sample_col, 'MPCpercent']].drop_duplicates().set_index(sample_col)
    sample_proportions['MPCpercent'] = sample_proportions.index.map(sample_groups['MPCpercent'])
    sample_proportions = sample_proportions.reset_index()
    return sample_proportions

def kruskal_wallis_test_by_celltype(sample_proportions, sample_col='sample'):
    print("\n进行Kruskal-Wallis检验...")
    celltype_cols = [col for col in sample_proportions.columns 
                     if col not in [sample_col, 'MPCpercent']]
    results = {}
    for celltype in celltype_cols:
        groups = sample_proportions['MPCpercent'].unique()
        group_data = []
        for group in sorted(groups):
            group_values = sample_proportions[sample_proportions['MPCpercent'] == group][celltype].values
            group_values = group_values[~np.isnan(group_values)]
            if len(group_values) > 0:
                group_data.append(group_values)
        if len(group_data) >= 2:
            try:
                stat, p_value = stats.kruskal(*group_data)
                test_used = 'Kruskal-Wallis'
            except Exception as e:
                print(f"  细胞类型 {celltype}: 检验失败 - {e}")
                p_value = np.nan
                test_used = '检验失败'
        else:
            p_value = np.nan
            test_used = '数据不足'
        significant = p_value < 0.05 if not np.isnan(p_value) else False
        results[celltype] = {'test': test_used, 'p_value': p_value, 'significant': significant}
    return results

def plot_trend_boxplots(sample_proportions, significant_celltypes, 
                        sample_col='sample', output_dir='MPC_trend_plots'):
    print(f"\n绘制带折线趋势的箱线图...")
    os.makedirs(output_dir, exist_ok=True)

    # 增大字体设置
    original_rc = {k: mpl.rcParams[k] for k in ['font.size', 'axes.titlesize', 'axes.labelsize',
                                                'xtick.labelsize', 'ytick.labelsize', 'legend.fontsize',
                                                'figure.titlesize']}
    mpl.rcParams.update({'font.size': 12, 'axes.titlesize': 16, 'axes.labelsize': 14,
                         'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 12,
                         'figure.titlesize': 18})

    mpc_groups = sorted([g for g in sample_proportions['MPCpercent'].unique() 
                         if isinstance(g, str) and g.startswith('T')])
    if not mpc_groups:
        mpc_groups = sorted(sample_proportions['MPCpercent'].unique())
    print(f"分组顺序: {mpc_groups}")

    for i, celltype in enumerate(significant_celltypes, 1):
        print(f"  处理细胞类型 {i}/{len(significant_celltypes)}: {celltype}")
        fig, ax = plt.subplots(figsize=(10, 6))
        plot_data = []
        for group in mpc_groups:
            group_data = sample_proportions[sample_proportions['MPCpercent'] == group]
            values = group_data[celltype].dropna().values
            for value in values:
                plot_data.append({'MPCpercent': group, 'Proportion': value, 'CellType': celltype})
        plot_df = pd.DataFrame(plot_data)
        if plot_df.empty:
            print(f"    警告: {celltype} 没有有效数据")
            plt.close(fig)
            continue
        sns.boxplot(data=plot_df, x='MPCpercent', y='Proportion', 
                    order=mpc_groups, palette='Set2', width=0.6, ax=ax)
        sns.stripplot(data=plot_df, x='MPCpercent', y='Proportion', 
                      order=mpc_groups, color='black', size=4, alpha=0.6, jitter=True, ax=ax)
        median_values = plot_df.groupby('MPCpercent')['Proportion'].median().reindex(mpc_groups)
        ax.plot(range(len(mpc_groups)), median_values.values, 
                marker='o', markersize=8, linewidth=2, color='red', label='Median Trend')
        ax.set_xlabel('MPC Percent Group')
        ax.set_ylabel('Cell Proportion (%)')
        ax.set_title(f'{celltype} - Proportion by MPC Group')
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
        ax.legend()
        plt.tight_layout()

        safe_celltype = celltype.replace('/', '_').replace('\\', '_').replace(':', '_')
        png_path = os.path.join(output_dir, f"{safe_celltype}_trend_boxplot.png")
        pdf_path = os.path.join(output_dir, f"{safe_celltype}_trend_boxplot.pdf")
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
        plt.close(fig)

    mpl.rcParams.update(original_rc)
    print(f"\n所有图形已保存到目录: {output_dir}/")

def create_summary_plot(sample_proportions, significant_celltypes, 
                        sample_col='sample', output_dir='MPC_trend_plots'):
    print("\n创建显著细胞类型的汇总图...")
    original_rc = {k: mpl.rcParams[k] for k in ['font.size', 'axes.titlesize', 'axes.labelsize',
                                                'xtick.labelsize', 'ytick.labelsize', 'legend.fontsize',
                                                'figure.titlesize']}
    mpl.rcParams.update({'font.size': 12, 'axes.titlesize': 16, 'axes.labelsize': 14,
                         'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 12,
                         'figure.titlesize': 18})

    mpc_groups = sorted([g for g in sample_proportions['MPCpercent'].unique() 
                         if isinstance(g, str) and g.startswith('T')])
    if not mpc_groups:
        mpc_groups = sorted(sample_proportions['MPCpercent'].unique())

    n_celltypes = len(significant_celltypes)
    n_cols = 3
    n_rows = (n_celltypes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*4))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for idx, celltype in enumerate(significant_celltypes):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        plot_data = []
        for group in mpc_groups:
            group_data = sample_proportions[sample_proportions['MPCpercent'] == group]
            values = group_data[celltype].dropna().values
            for value in values:
                plot_data.append({'MPCpercent': group, 'Proportion': value, 'CellType': celltype})
        plot_df = pd.DataFrame(plot_data)
        if plot_df.empty:
            ax.text(0.5, 0.5, f'No data for\n{celltype}', ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        sns.boxplot(data=plot_df, x='MPCpercent', y='Proportion', 
                    order=mpc_groups, palette='Set2', width=0.6, ax=ax)
        median_values = plot_df.groupby('MPCpercent')['Proportion'].median().reindex(mpc_groups)
        ax.plot(range(len(mpc_groups)), median_values.values, 
                marker='o', markersize=6, linewidth=2, color='red')
        ax.set_title(f'{celltype}')
        ax.set_xlabel('')
        ax.tick_params(axis='x', rotation=45)
        if col == 0:
            ax.set_ylabel('Proportion (%)')
        else:
            ax.set_ylabel('')

    for idx in range(len(significant_celltypes), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)

    plt.suptitle('Significant Cell Types by MPC Group (with Median Trend)')
    plt.tight_layout()

    png_path = os.path.join(output_dir, "all_significant_celltypes_summary.png")
    pdf_path = os.path.join(output_dir, "all_significant_celltypes_summary.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()

    mpl.rcParams.update(original_rc)
    print(f"汇总图已保存到: {png_path} 和 {pdf_path}")

def main():
    h5ad_file = "adata_filtered_with_MPCpercent.h5ad"
    sample_col = 'sample'

    print("="*60)
    print("读取h5ad数据...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"数据形状: {adata.shape}")
    print(f"细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
    print(f"\n数据列: {adata.obs.columns.tolist()}")
    print(f"MPCpercent分组分布:")
    print(adata.obs['MPCpercent'].value_counts().sort_index())
    print(f"细胞类型数量: {adata.obs['annotation'].nunique()}")
    print(f"样本数量: {adata.obs[sample_col].nunique()}")

    sample_proportions = calculate_sample_celltype_proportions(adata, sample_col)
    sample_proportions.to_csv("sample_celltype_proportions.csv", index=False, encoding='utf-8-sig')
    print(f"\n样本细胞类型比例数据已保存到: sample_celltype_proportions.csv")

    results = kruskal_wallis_test_by_celltype(sample_proportions, sample_col)
    significant_celltypes = [celltype for celltype, result in results.items() 
                             if result.get('significant', False)]

    results_df = pd.DataFrame.from_dict(results, orient='index')
    results_df.index.name = 'celltype'
    results_df.to_csv("kruskal_wallis_results.csv", encoding='utf-8-sig')
    print(f"统计检验结果已保存到: kruskal_wallis_results.csv")

    print("\n" + "="*60)
    print("分析结果摘要")
    print("="*60)
    print(f"总细胞类型数: {len(results)}")
    print(f"有显著差异的细胞类型数: {len(significant_celltypes)}")

    if significant_celltypes:
        print(f"\n有显著差异的细胞类型 (p < 0.05):")
        for celltype in significant_celltypes:
            p_val = results[celltype].get('p_value', np.nan)
            print(f"  {celltype}: p = {p_val:.4e}")

        output_dir = "MPC_trend_plots"
        plot_trend_boxplots(sample_proportions, significant_celltypes, sample_col, output_dir)
        create_summary_plot(sample_proportions, significant_celltypes, sample_col, output_dir)

        print(f"\n生成的文件:")
        print(f"  1. sample_celltype_proportions.csv")
        print(f"  2. kruskal_wallis_results.csv")
        print(f"  3. {output_dir}/ 目录下包含所有图形的 PNG 和 PDF")
    else:
        print("\n由于没有显著差异的细胞类型，不生成图形。")

    print("\n分析完成!")

if __name__ == "__main__":
    main()
