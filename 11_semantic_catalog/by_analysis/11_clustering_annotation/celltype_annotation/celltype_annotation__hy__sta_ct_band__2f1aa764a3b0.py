import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import seaborn as sns
from scipy.stats import mannwhitneyu
from itertools import combinations
import os
import warnings
warnings.filterwarnings('ignore')

# ================= 字体设置 (CNS 发表级) =================
# 如果字体路径不对，请注释掉下面这行，改用系统默认
font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')

# ================== 全局设置 ==================
plt.rcParams['font.family'] = 'arial'
plt.rcParams['font.weight'] = 'heavy'        # 字体加粗
plt.rcParams['pdf.fonttype'] = 42            # PDF Type 42 (TrueType)
plt.rcParams['ps.fonttype'] = 42             # PS Type 42 (TrueType)

# ================== 放大所有字体 (至少增加 3 个字号) ==================
plt.rcParams['font.size'] = 18               # 基础字号 +4
plt.rcParams['axes.labelsize'] = 20          # X/Y轴标签 +4
plt.rcParams['axes.titlesize'] = 22          # 子图标题 +4
plt.rcParams['xtick.labelsize'] = 18         # X轴刻度标签 +4
plt.rcParams['ytick.labelsize'] = 18         # Y轴刻度标签 +4
plt.rcParams['legend.fontsize'] = 17         # 图例文字 +4
plt.rcParams['figure.titlesize'] = 24        # 总标题 +4

# ==================== 配置参数 ====================
# 确保 band_final 的分类严格按照此顺序排列
BAND_ORDER = ['N', 'AN1', 'AN2', 'A', 'CA1', 'CA2', 'C']
SAVE_DIR = 'sta_ct_band'
DATA_PATH = 'adata_with_bands.h5ad'

# ==================== 工具函数：显著性标注 ====================
def add_stat_annotation(ax, group_idx_map, p_values_dict, y_max):
    """
    在箱线图上绘制显著性差异连线和星号 (*, **, ***)
    """
    current_y = y_max * 1.05
    step = y_max * 0.05
    
    for (group1, group2), p_val in p_values_dict.items():
        if p_val < 0.05:
            x1 = group_idx_map[group1]
            x2 = group_idx_map[group2]
            
            if p_val < 0.001:
                star = '***'
            elif p_val < 0.01:
                star = '**'
            else:
                star = '*'
                
            ax.plot([x1, x1, x2, x2], [current_y, current_y + step, current_y + step, current_y], 
                    lw=1.2, c='black')
            # 星号字体同步加大 (改为 20，增加 4 个字号)
            ax.text((x1 + x2) * 0.5, current_y + step, star, ha='center', va='bottom', fontsize=20)
            
            current_y += step * 1.5 

# ==================== 第一部分：堆积柱状图 ====================
def plot_stacked_bar_chart(adata, save_dir):
    mask = adata.obs['band_final'].isin(BAND_ORDER)
    df = adata.obs[mask][['band_final', 'celltype']]
    
    cross_tab = pd.crosstab(df['band_final'], df['celltype'], normalize='index') * 100
    cross_tab = cross_tab.reindex(BAND_ORDER).fillna(0) 
    
    fig, ax = plt.subplots(figsize=(12, 6))
    cross_tab.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', alpha=0.85)
    
    ax.set_ylabel('Percentage (%)')
    ax.set_xlabel('Band Final Region')
    ax.set_title('Cell Type Composition across Different Band Regions', fontweight='heavy')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    base_name = os.path.join(save_dir, '01_stacked_barchart_celltype')
    fig.savefig(f'{base_name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{base_name}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"堆积柱状图已保存: {base_name}.png 和 {base_name}.pdf")

# ==================== 第二部分：多系列箱线图（按比例） ====================
def plot_boxplots_by_sample(adata, save_dir):
    # 检查列名
    if 'sample_new' not in adata.obs.columns:
        raise ValueError("未在 obs 中找到 'sample_new' 列，请先构建 sample_new 分组列！")
    
    mask = adata.obs['band_final'].isin(BAND_ORDER)
    df = adata.obs[mask][['sample_new', 'band_final', 'celltype']]
    
    # 计算每个样本、每个带区内，每种细胞类型的比例
    type_counts = df.groupby(['sample_new', 'band_final', 'celltype']).size().reset_index(name='count')
    total_counts = df.groupby(['sample_new', 'band_final']).size().reset_index(name='total')
    df_merged = pd.merge(type_counts, total_counts, on=['sample_new', 'band_final'])
    df_merged['proportion'] = df_merged['count'] / df_merged['total']

    cell_types = df_merged['celltype'].unique()
    group_idx_map = {group: i for i, group in enumerate(BAND_ORDER)}
    
    for cell_type in cell_types:
        subset = df_merged[df_merged['celltype'] == cell_type]
        subset = subset[subset['band_final'].isin(BAND_ORDER)]
        subset['band_final'] = pd.Categorical(subset['band_final'], categories=BAND_ORDER, ordered=True)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        sns.boxplot(x='band_final', y='proportion', data=subset, ax=ax, palette='Set3')
        sns.stripplot(x='band_final', y='proportion', data=subset, ax=ax, color='black', alpha=0.3, size=3)
        
        ax.set_title(f'Proportion of {cell_type} by Band Region', fontweight='heavy')
        ax.set_xlabel('Band Region')
        ax.set_ylabel(f'Proportion of {cell_type} Cells per Sample')
        
        # 计算显著性
        p_values_dict = {}
        group_data = {band: subset[subset['band_final'] == band]['proportion'].values for band in BAND_ORDER}
        
        for band1, band2 in combinations(BAND_ORDER, 2):
            data1 = group_data[band1]
            data2 = group_data[band2]
            if len(data1) < 2 or len(data2) < 2:
                continue
            stat, p_val = mannwhitneyu(data1, data2, alternative='two-sided')
            p_values_dict[(band1, band2)] = p_val
        
        y_max = subset['proportion'].max() if len(subset) > 0 else 1
        add_stat_annotation(ax, group_idx_map, p_values_dict, y_max)
        
        plt.tight_layout()
        
        base_name = os.path.join(save_dir, f'02_boxplot_{cell_type}_by_sample')
        fig.savefig(f'{base_name}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{base_name}.pdf', bbox_inches='tight')
        plt.close(fig)
        print(f"细胞类型 {cell_type} 比例箱线图已保存: {base_name}.png 和 {base_name}.pdf")

# ==================== 主函数 ====================
def run_statistical_analysis(adata_path, save_dir=SAVE_DIR):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        print(f"创建文件夹: {save_dir}")

    print("正在加载 AnnData 数据...")
    try:
        adata = sc.read_h5ad(adata_path)
    except Exception as e:
        print(f"读取文件错误: {e}，请确保文件路径正确！")
        return

    if 'band_final' not in adata.obs.columns:
        raise ValueError("数据中未找到 'band_final' 列，请先运行带状划分脚本！")
    if 'celltype' not in adata.obs.columns:
        raise ValueError("数据中未找到 'celltype' 列，请确保您的 obs 中包含细胞类型注释！")
    if 'sample_new' not in adata.obs.columns:
        print("警告：未找到 'sample_new' 列，请确保您的 obs 中有此分组列！")

    print("\n========== 开始绘制堆积柱状图 ==========")
    plot_stacked_bar_chart(adata, save_dir)
    
    print("\n========== 开始绘制多系列箱线图（基于样本比例） ==========")
    plot_boxplots_by_sample(adata, save_dir)
    
    print(f"\n分析结束！所有图表已保存至文件夹: {save_dir}")

if __name__ == "__main__":
    run_statistical_analysis(DATA_PATH)
