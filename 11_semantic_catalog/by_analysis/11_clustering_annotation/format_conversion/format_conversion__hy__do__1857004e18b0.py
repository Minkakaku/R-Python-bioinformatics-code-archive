# -*- coding: utf-8 -*-
"""
单细胞数据差异分析及箱线图绘制脚本
功能：基于TREATMENT分组，以样本维度计算细胞类型差异，并可视化箱线图
依赖库：scanpy, pandas, numpy, matplotlib, seaborn, scipy
"""

# 1. 导入必要库
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图样式
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']  # 支持中文显示
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号
sns.set_style("whitegrid")  # 设置绘图风格

# 2. 加载数据
def load_data(file_path):
    """
    加载.h5ad格式的单细胞数据
    """
    try:
        adata = sc.read_h5ad(file_path)
        print("数据加载成功!")
        print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
        print(f"观察列(obs): {list(adata.obs.columns)}")
        return adata
    except Exception as e:
        print(f"数据加载失败: {e}")
        return None

# 3. 数据预处理和聚合
def preprocess_and_aggregate(adata, sample_col='sample', celltype_col='celltype', treatment_col='TREATMENT'):
    """
    预处理数据并聚合到样本水平
    返回每个样本中每种细胞类型的比例和表达量均值
    """
    # 检查必要的列是否存在
    required_cols = [sample_col, celltype_col, treatment_col]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        print(f"缺少必要列: {missing_cols}")
        return None
    
    # 提取样本元数据
    sample_meta = adata.obs[[sample_col, treatment_col]].drop_duplicates()
    
    # 计算每个样本中每种细胞类型的细胞数量
    celltype_counts = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col])
    
    # 计算每个样本的细胞总数
    total_cells_per_sample = celltype_counts.sum(axis=1)
    
    # 计算每种细胞类型在每个样本中的比例
    celltype_proportions = celltype_counts.div(total_cells_per_sample, axis=0)
    
    # 重置索引以便合并
    celltype_proportions = celltype_proportions.reset_index()
    celltype_proportions.rename(columns={'index': sample_col}, inplace=True)
    
    # 融合数据框，将细胞类型作为一列
    proportion_long = pd.melt(
        celltype_proportions, 
        id_vars=[sample_col], 
        value_name='Proportion', 
        var_name='CellType'
    )
    
    # 添加治疗组信息
    proportion_long = proportion_long.merge(sample_meta, on=sample_col)
    
    print(f"聚合数据完成，包含 {len(proportion_long)} 个样本-细胞类型组合")
    return proportion_long

# 4. 计算差异显著性
def calculate_differences(proportion_data, treatment_col='TREATMENT'):
    """
    计算每种细胞类型在不同治疗组间的统计差异
    """
    results = []
    celltypes = proportion_data['CellType'].unique()
    
    for celltype in celltypes:
        celltype_data = proportion_data[proportion_data['CellType'] == celltype]
        
        # 获取治疗组
        treatments = celltype_data[treatment_col].unique()
        if len(treatments) < 2:
            print(f"细胞类型 {celltype} 只有一组数据，无法计算差异")
            continue
        
        # 提取每组的数据
        groups = []
        for treatment in treatments:
            group_data = celltype_data[celltype_data[treatment_col] == treatment]['Proportion']
            groups.append(group_data)
        
        # 独立样本t检验（假设两组）
        if len(treatments) == 2:
            t_stat, p_value = stats.ttest_ind(groups[0], groups[1], nan_policy='omit')
            test_type = "t-test"
        else:
            # 多组使用ANOVA
            f_stat, p_value = stats.f_oneway(*groups)
            test_type = "ANOVA"
        
        # 计算每组的平均比例
        means = {f"{treatment}_mean": group.mean() for treatment, group in zip(treatments, groups)}
        
        results.append({
            'CellType': celltype,
            'Test': test_type,
            'P_Value': p_value,
            **means
        })
    
    results_df = pd.DataFrame(results)
    return results_df

# 5. 绘制箱线图
def plot_boxplots(proportion_data, results_df, treatment_col='TREATMENT', 
                  figsize=(15, 10), top_n=10):
    """
    绘制细胞类型比例的箱线图，突出显示显著差异
    """
    # 按P值排序，选择最显著的细胞类型
    if len(results_df) > 0:
        sig_results = results_df.nsmallest(top_n, 'P_Value')
        sig_celltypes = sig_results['CellType'].tolist()
        plot_data = proportion_data[proportion_data['CellType'].isin(sig_celltypes)]
    else:
        plot_data = proportion_data
        sig_celltypes = proportion_data['CellType'].unique()
    
    # 创建图形
    n_types = len(sig_celltypes)
    n_cols = min(3, n_types)
    n_rows = (n_types + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_types == 1:
        axes = np.array([axes])
    if n_rows * n_cols == 1:
        axes = axes.reshape(1, 1)
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # 为每个细胞类型绘制箱线图
    for idx, celltype in enumerate(sig_celltypes):
        row = idx // n_cols
        col = idx % n_cols
        
        if n_rows == 1 and n_cols == 1:
            ax = axes
        elif n_rows == 1:
            ax = axes[col]
        elif n_cols == 1:
            ax = axes[row]
        else:
            ax = axes[row, col]
        
        celltype_data = plot_data[plot_data['CellType'] == celltype]
        
        # 绘制箱线图
        sns.boxplot(data=celltype_data, x=treatment_col, y='Proportion', ax=ax, 
                   palette="Set2", showfliers=True)
        
        # 添加数据点
        sns.stripplot(data=celltype_data, x=treatment_col, y='Proportion', ax=ax,
                     color='black', size=4, jitter=0.2, alpha=0.6)
        
        # 添加标题和显著性标记
        p_value = results_df[results_df['CellType'] == celltype]['P_Value'].values[0] if len(results_df) > 0 else 1
        if p_value < 0.001:
            sig_symbol = '***'
        elif p_value < 0.01:
            sig_symbol = '**'
        elif p_value < 0.05:
            sig_symbol = '*'
        else:
            sig_symbol = 'NS'
        
        ax.set_title(f'{celltype}\nP-value: {p_value:.4f} {sig_symbol}', fontsize=12)
        ax.set_xlabel('')
        ax.set_ylabel('proportion', fontsize=10)
        
        # 旋转x轴标签
        ax.tick_params(axis='x', rotation=45)
    
    # 隐藏多余的子图
    for idx in range(len(sig_celltypes), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        if n_rows == 1 and n_cols == 1:
            fig.delaxes(axes)
        elif n_rows == 1:
            fig.delaxes(axes[col])
        elif n_cols == 1:
            fig.delaxes(axes[row])
        else:
            fig.delaxes(axes[row, col])
    
    plt.tight_layout()
    plt.suptitle(f'proportion {treatment_col}', y=1.02, fontsize=16)
    return fig

# 6. 主函数
def main():
    # 加载数据
    adata = load_data('extracted_specific_patients.h5ad')
    if adata is None:
        return
    
    # 显示数据基本信息
    print("\n数据基本信息:")
    print(adata.obs.head())
    
    # 预处理和聚合数据
    proportion_data = preprocess_and_aggregate(adata)
    if proportion_data is None:
        return
    
    # 计算统计差异
    results_df = calculate_differences(proportion_data)
    
    # 显示差异分析结果
    print("\n差异分析结果:")
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    print(results_df.round(4))
    
    # 绘制箱线图
    fig = plot_boxplots(proportion_data, results_df)
    
    # 保存结果
    proportion_data.to_csv('celltype_proportions.csv', index=False)
    results_df.to_csv('statistical_results.csv', index=False)
    plt.savefig('celltype_boxplots.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\n分析完成! 结果已保存到当前目录。")

# 7. 执行脚本
if __name__ == "__main__":
    main()
