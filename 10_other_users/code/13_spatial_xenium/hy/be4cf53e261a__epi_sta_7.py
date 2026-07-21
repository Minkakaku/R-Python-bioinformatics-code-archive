import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib import rcParams

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def create_output_folder(folder_name='make_band_new_7_epi'):
    """
    创建输出文件夹
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

def load_and_filter_data(h5ad_path):
    """
    加载数据并按照条件过滤
    """
    print(f"加载数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据形状: {adata.shape}")
    print(f"原始数据列: {adata.obs.columns.tolist()}")
    
    # 检查必要的列是否存在
    if 'site_band' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'site_band'列")
    
    if 'annotation' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'annotation'列")
    
    # 定义我们要保留的site_band值
    target_bands = ['inner', 'in_d','in_c','in_b', 'in_a', 'on', 'out_a', 'out_b','out_c','out_d']
    
    # 第一步：提取指定site_band的数据
    band_mask = adata.obs['site_band'].isin(target_bands)
    adata_filtered = adata[band_mask].copy()
    
    print(f"提取指定site_band后的数据形状: {adata_filtered.shape}")
    print(f"各site_band的细胞数量:\n{adata_filtered.obs['site_band'].value_counts()}")
    
    # 第二步：提取包含"Epithelium"的细胞
    # 先找出所有包含"Epithelium"的annotation
    epithelium_annotations = adata_filtered.obs['annotation'][
        adata_filtered.obs['annotation'].astype(str).str.contains('Epithelium', case=False, na=False)
    ].unique()
    
    print(f"找到包含'Epithelium'的annotation类型: {len(epithelium_annotations)} 种")
    print(f"包含'Epithelium'的annotation列表:\n{epithelium_annotations}")
    
    # 过滤出包含"Epithelium"的细胞
    epithelium_mask = adata_filtered.obs['annotation'].astype(str).str.contains('Epithelium', case=False, na=False)
    adata_epi = adata_filtered[epithelium_mask].copy()
    
    print(f"提取包含'Epithelium'细胞后的数据形状: {adata_epi.shape}")
    print(f"各上皮细胞类型的数量:\n{adata_epi.obs['annotation'].value_counts()}")
    
    return adata_epi, target_bands

def calculate_proportions(adata_epi, target_bands):
    """
    计算每个site_band区域中特定上皮亚型的比例
    """
    # 定义我们要分析的三种上皮亚型
    target_epi_subtypes = [
        "Stress-responsive_glandular_Epithelium",
        "EMT_Epithelium", 
        "Proliferative_Epithelium"
    ]
    
    # 检查这些亚型是否存在于数据中
    available_annotations = adata_epi.obs['annotation'].unique()
    print(f"\n可用的annotation类型: {available_annotations}")
    
    # 找出实际存在于数据中的亚型
    existing_subtypes = [subtype for subtype in target_epi_subtypes 
                        if any(subtype in str(anno) for anno in available_annotations)]
    
    if not existing_subtypes:
        # 如果没有找到精确匹配，尝试部分匹配
        for subtype in target_epi_subtypes:
            for anno in available_annotations:
                if subtype.lower().replace('_', '').replace(' ', '') in str(anno).lower().replace('_', '').replace(' ', ''):
                    existing_subtypes.append(anno)
                    print(f"找到部分匹配: {subtype} -> {anno}")
    
    if not existing_subtypes:
        print("警告: 没有找到目标上皮亚型，尝试匹配所有包含这些关键词的细胞")
        # 创建新的列来标记这些亚型
        adata_epi.obs['epi_subtype'] = 'Other_Epithelium'
        for subtype in target_epi_subtypes:
            subtype_mask = adata_epi.obs['annotation'].astype(str).str.contains(
                subtype.split('_')[0], case=False, na=False
            )
            adata_epi.obs.loc[subtype_mask, 'epi_subtype'] = subtype
        existing_subtypes = target_epi_subtypes
        annotation_col = 'epi_subtype'
    else:
        annotation_col = 'annotation'
    
    print(f"将分析的上皮亚型: {existing_subtypes}")
    
    # 创建一个DataFrame来存储结果
    results = pd.DataFrame(index=target_bands, columns=existing_subtypes)
    
    # 计算每个site_band的总上皮细胞数
    total_per_band = adata_epi.obs['site_band'].value_counts()
    print(f"\n各site_band区域的总上皮细胞数:\n{total_per_band}")
    
    # 计算每个site_band区域中每种上皮亚型的比例
    for band in target_bands:
        # 获取该band的细胞
        band_mask = adata_epi.obs['site_band'] == band
        band_total = total_per_band.get(band, 0)
        
        if band_total > 0:
            for subtype in existing_subtypes:
                # 计算该亚型在该band中的细胞数
                if annotation_col == 'annotation':
                    # 精确匹配或包含匹配
                    subtype_mask = band_mask & adata_epi.obs[annotation_col].astype(str).str.contains(
                        subtype, case=False, na=False
                    )
                else:
                    subtype_mask = band_mask & (adata_epi.obs[annotation_col] == subtype)
                
                subtype_count = subtype_mask.sum()
                # 计算比例
                proportion = subtype_count / band_total * 100
                results.loc[band, subtype] = proportion
        else:
            # 如果没有细胞，比例为0
            for subtype in existing_subtypes:
                results.loc[band, subtype] = 0.0
    
    # 确保所有值为数值类型
    results = results.astype(float)
    
    print(f"\n各上皮亚型在不同site_band区域的比例 (%):")
    print(results)
    
    return results, existing_subtypes

def plot_line_chart(results, existing_subtypes, save_dir='make_band_new_7_epi'):
    """
    绘制折线图
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 定义颜色方案
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    marker_styles = ['o', 's', '^', 'D', 'v', '<']
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 左图：折线图
    ax1 = axes[0]
    
    x_positions = np.arange(len(results.index))
    bar_width = 0.8 / len(existing_subtypes)
    
    for i, subtype in enumerate(existing_subtypes):
        # 获取该亚型的数据
        y_values = results[subtype].values
        
        # 绘制折线
        ax1.plot(x_positions, y_values, marker=marker_styles[i % len(marker_styles)], 
                linewidth=2, markersize=8, label=subtype, color=colors[i % len(colors)])
    
    # 设置x轴
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(results.index, rotation=0, fontsize=10)
    ax1.set_xlabel('Site Band Region', fontsize=12)
    ax1.set_ylabel('Proportion of Epithelial Cells (%)', fontsize=12)
    ax1.set_title('Proportion of Epithelial Subtypes Across Site Bands', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # 设置y轴范围
    y_max = results.values.max() * 1.1
    ax1.set_ylim(0, max(y_max, 10))
    
    # 右图：堆叠柱状图
    ax2 = axes[1]
    
    # 准备堆叠柱状图数据
    bottom = np.zeros(len(results.index))
    for i, subtype in enumerate(existing_subtypes):
        y_values = results[subtype].values
        ax2.bar(x_positions, y_values, bottom=bottom, label=subtype, 
               color=colors[i % len(colors)], alpha=0.8, edgecolor='black', linewidth=0.5)
        bottom += y_values
    
    # 设置x轴
    ax2.set_xticks(x_positions)
    ax2.set_xticklabels(results.index, rotation=0, fontsize=10)
    ax2.set_xlabel('Site Band Region', fontsize=12)
    ax2.set_ylabel('Proportion of Epithelial Cells (%)', fontsize=12)
    ax2.set_title('Stacked Proportion of Epithelial Subtypes', fontsize=14, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    # 在柱状图上添加总计文本
    for i, band in enumerate(results.index):
        total = results.loc[band].sum()
        ax2.text(i, total + 0.5, f'{total:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.suptitle('Epithelial Cell Subtype Distribution Across Spatial Bands', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为suptitle留出空间
    
    # 保存图形
    fig_path_png = os.path.join(save_dir, 'epithelial_subtype_proportions.png')
    fig_path_pdf = os.path.join(save_dir, 'epithelial_subtype_proportions.pdf')
    fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    print(f"\n保存图形: {fig_path_png}")
    print(f"保存图形: {fig_path_pdf}")
    
    plt.show()
    
    return fig

def create_heatmap_visualization(adata_epi, results, existing_subtypes, save_dir='make_band_new_7_epi'):
    """
    创建热图可视化
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 创建热图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 创建热图数据
    heatmap_data = results.copy()
    
    # 创建热图
    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': 'Proportion (%)'}, ax=ax)
    
    ax.set_title('Heatmap of Epithelial Subtype Proportions by Site Band', fontsize=14, fontweight='bold')
    ax.set_xlabel('Epithelial Subtype', fontsize=12)
    ax.set_ylabel('Site Band', fontsize=12)
    
    plt.tight_layout()
    
    # 保存热图
    heatmap_path_png = os.path.join(save_dir, 'epithelial_subtype_heatmap.png')
    heatmap_path_pdf = os.path.join(save_dir, 'epithelial_subtype_heatmap.pdf')
    fig.savefig(heatmap_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(heatmap_path_pdf, bbox_inches='tight')
    print(f"保存热图: {heatmap_path_png}")
    print(f"保存热图: {heatmap_path_pdf}")
    
    plt.show()
    
    return fig

def save_results(results, adata_epi, save_dir='make_band_new_7_epi'):
    """
    保存结果到文件
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 保存比例结果
    results_csv = os.path.join(save_dir, 'epithelial_subtype_proportions.csv')
    results.to_csv(results_csv)
    print(f"\n保存比例结果: {results_csv}")
    
    # 保存每个site_band的详细统计
    stats_summary = os.path.join(save_dir, 'epithelial_cell_statistics.txt')
    
    with open(stats_summary, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("     Epithelial Cell Subtype Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("1. DATA OVERVIEW\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total epithelial cells analyzed: {len(adata_epi)}\n")
        f.write(f"Site bands analyzed: {list(results.index)}\n")
        f.write(f"Epithelial subtypes analyzed: {list(results.columns)}\n\n")
        
        f.write("2. CELL COUNTS BY SITE BAND\n")
        f.write("-" * 40 + "\n")
        band_counts = adata_epi.obs['site_band'].value_counts()
        for band, count in band_counts.items():
            f.write(f"{band}: {count} cells\n")
        
        f.write("\n3. EPITHELIAL SUBTYPE PROPORTIONS (%)\n")
        f.write("-" * 40 + "\n")
        
        # 写入表格格式的结果
        f.write("SiteBand")
        for subtype in results.columns:
            f.write(f"\t{subtype}")
        f.write("\n")
        
        for band in results.index:
            f.write(f"{band}")
            for subtype in results.columns:
                proportion = results.loc[band, subtype]
                f.write(f"\t{proportion:.1f}")
            f.write("\n")
        
        f.write("\n4. SPATIAL PATTERN ANALYSIS\n")
        f.write("-" * 40 + "\n")
        
        # 分析每种亚型的空间分布模式
        for subtype in results.columns:
            f.write(f"\n{subtype}:\n")
            proportions = results[subtype]
            
            # 找到最高和最低比例的区域
            max_band = proportions.idxmax()
            min_band = proportions.idxmin()
            max_val = proportions.max()
            min_val = proportions.min()
            
            f.write(f"  Highest proportion: {max_band} ({max_val:.1f}%)\n")
            f.write(f"  Lowest proportion: {min_band} ({min_val:.1f}%)\n")
            
            # 计算趋势（inner到outer）
            inner_val = results.loc['inner', subtype] if 'inner' in results.index else 0
            outer_val = results.loc['out_b', subtype] if 'out_b' in results.index else 0
            if inner_val > 0 or outer_val > 0:
                trend = (outer_val - inner_val) / max(inner_val, 1) * 100
                direction = "increases" if trend > 0 else "decreases"
                f.write(f"  Trend (inner→outer): {abs(trend):.1f}% {direction}\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Analysis completed successfully!\n")
        f.write("=" * 60 + "\n")
    
    print(f"保存统计摘要: {stats_summary}")
    
    return results_csv, stats_summary

def main():
    """
    主函数
    """
    # 数据路径
    h5ad_path = 'make_band_new_7/adata_with_band_by_sample.h5ad'
    
    # 检查文件是否存在
    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"未找到h5ad文件: {h5ad_path}")
    
    # 输出目录
    output_dir = 'make_band_new_7_epi'
    
    print("=" * 60)
    print("     Epithelial Cell Subtype Spatial Analysis")
    print("=" * 60)
    
    # 1. 加载并过滤数据
    adata_epi, target_bands = load_and_filter_data(h5ad_path)
    
    # 2. 计算比例
    results, existing_subtypes = calculate_proportions(adata_epi, target_bands)
    
    # 3. 绘制折线图
    fig_line = plot_line_chart(results, existing_subtypes, save_dir=output_dir)
    
    # 4. 创建热图
    fig_heatmap = create_heatmap_visualization(adata_epi, results, existing_subtypes, save_dir=output_dir)
    
    # 5. 保存结果
    results_csv, stats_summary = save_results(results, adata_epi, save_dir=output_dir)
    
    # 6. 显示关键结果
    print("\n" + "=" * 60)
    print("关键结果摘要")
    print("=" * 60)
    
    for subtype in existing_subtypes:
        print(f"\n{subtype} 的空间分布:")
        for band in target_bands:
            proportion = results.loc[band, subtype]
            print(f"  {band}: {proportion:.1f}%")
    
    print("\n" + "=" * 60)
    print("分析完成！所有结果已保存到目录:")
    print(f"  {output_dir}")
    print("=" * 60)
    
    return adata_epi, results

# 运行主函数
if __name__ == "__main__":
    adata_epi, results = main()
