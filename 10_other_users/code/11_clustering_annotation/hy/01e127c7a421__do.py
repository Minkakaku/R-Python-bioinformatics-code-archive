import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图风格
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def load_and_check_data(file_path):
    """加载数据并检查基本结构"""
    print("=== 加载数据 ===")
    adata = sc.read_h5ad(file_path)
    print(f"数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列是否存在
    required_cols = ['celltype', 'group_new']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"错误: 缺少必要的列: {missing_cols}")
        print(f"可用的列: {list(adata.obs.columns)}")
        return None
    
    print(f"\n=== 数据基本信息 ===")
    print(f"细胞类型数量: {adata.obs['celltype'].nunique()}")
    print(f"细胞类型分布:\n{adata.obs['celltype'].value_counts()}")
    print(f"\ngroup_new分组数量: {adata.obs['group_new'].nunique()}")
    print(f"group_new分布:\n{adata.obs['group_new'].value_counts()}")
    
    return adata

def check_celltype_coverage(adata):
    """检查每个group_new分组中的细胞类型覆盖情况"""
    print("\n=== 各分组细胞类型覆盖情况 ===")
    coverage_table = pd.crosstab(adata.obs['group_new'], adata.obs['celltype'])
    print(coverage_table)
    
    # 检查是否有足够细胞进行差异分析
    min_cells_per_group = 10  # 每组至少需要10个细胞
    problematic_combinations = []
    
    for group in adata.obs['group_new'].unique():
        for celltype in adata.obs['celltype'].unique():
            count = len(adata[(adata.obs['group_new'] == group) & 
                            (adata.obs['celltype'] == celltype)])
            if 0 < count < min_cells_per_group:
                problematic_combinations.append((group, celltype, count))
    
    if problematic_combinations:
        print(f"\n警告: 以下分组细胞数量较少（<{min_cells_per_group}），可能影响差异分析结果:")
        for group, celltype, count in problematic_combinations:
            print(f"  {group} - {celltype}: {count}个细胞")
    
    return coverage_table

def perform_global_differential_analysis(adata):
    """执行全局差异分析：比较所有group_new分组"""
    print("\n=== 执行全局差异分析 ===")
    
    # 确保数据格式正确
    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray()
    
    # 预处理数据
    print("1. 数据预处理...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 设置group_new为观察因子
    adata.obs['group_new'] = adata.obs['group_new'].astype('category')
    
    # 使用Scanpy的rank_genes_groups函数进行差异分析
    print("2. 执行差异表达分析...")
    sc.tl.rank_genes_groups(adata, groupby='group_new', method='wilcoxon', 
                           use_raw=False, pts=True)
    
    print("3. 差异分析完成!")
    return adata

def perform_celltype_specific_analysis(adata, min_cells=20):
    """对每种细胞类型执行分组间的差异分析"""
    print("\n=== 执行细胞类型特异性差异分析 ===")
    
    results = {}
    celltypes = adata.obs['celltype'].unique()
    
    for i, celltype in enumerate(celltypes):
        print(f"分析细胞类型 {i+1}/{len(celltypes)}: {celltype}")
        
        # 提取该细胞类型的子集
        celltype_mask = adata.obs['celltype'] == celltype
        adata_subset = adata[celltype_mask].copy()
        
        # 检查是否有足够的分组和细胞
        groups = adata_subset.obs['group_new'].unique()
        if len(groups) < 2:
            print(f"  ⚠️  {celltype} 只有1个分组，跳过分析")
            continue
            
        # 检查每个分组的细胞数量
        group_counts = adata_subset.obs['group_new'].value_counts()
        if group_counts.min() < min_cells:
            print(f"  ⚠️  {celltype} 有些分组细胞数 < {min_cells}，跳过分析")
            continue
        
        # 执行差异分析
        try:
            sc.tl.rank_genes_groups(adata_subset, groupby='group_new', 
                                  method='wilcoxon', use_raw=False)
            results[celltype] = adata_subset
            print(f"  ✅ {celltype} 分析完成")
        except Exception as e:
            print(f"  ❌ {celltype} 分析失败: {e}")
    
    return results

def create_volcano_plots(adata, top_n=10):
    """创建火山图可视化差异表达基因"""
    print("\n=== 生成火山图 ===")
    
    # 获取所有比较组
    groups = adata.obs['group_new'].cat.categories
    n_groups = len(groups)
    
    if n_groups < 2:
        print("分组数量不足，无法创建火山图")
        return
    
    # 计算子图布局
    ncols = min(3, n_groups - 1)
    nrows = (n_groups - 1 + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    if nrows * ncols == 1:
        axes = np.array([axes])
    axes = axes.flat
    
    # 为每个分组创建火山图（以第一个分组为参考）
    reference_group = groups[0]
    for i, target_group in enumerate(groups[1:]):
        ax = axes[i]
        
        try:
            # 获取差异分析结果
            deg_results = sc.get.rank_genes_groups_df(adata, group=target_group)
            
            # 设置显著性阈值
            deg_results['-log10_pval_adj'] = -np.log10(deg_results['pvals_adj'])
            deg_results['significant'] = (deg_results['pvals_adj'] < 0.05) & (np.abs(deg_results['scores']) > 0.5)
            
            # 创建火山图
            scatter = ax.scatter(deg_results['scores'], deg_results['-log10_pval_adj'],
                               c=deg_results['significant'], cmap='viridis', alpha=0.6, s=20)
            
            ax.set_xlabel('Log2 Fold Change')
            ax.set_ylabel('-Log10(Adj. P-value)')
            ax.set_title(f'{target_group} vs {reference_group}')
            ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.8)
            ax.axvline(x=0, color='black', linestyle='-', alpha=0.8)
            
            # 标记top N基因
            top_genes = deg_results.nlargest(top_n, 'scores')['names'].tolist()
            for _, row in deg_results.iterrows():
                if row['names'] in top_genes:
                    ax.annotate(row['names'], (row['scores'], row['-log10_pval_adj']),
                               xytext=(5, 5), textcoords='offset points', fontsize=8)
            
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {str(e)}', transform=ax.transAxes, 
                   ha='center', va='center')
            ax.set_title(f'{target_group} vs {reference_group}')
    
    # 隐藏多余的子图
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def create_heatmap_visualization(adata, top_genes=20):
    """创建热图展示差异表达基因"""
    print("\n=== 生成热图 ===")
    
    try:
        # 获取每个分组的前top_genes个差异基因
        sc.pl.rank_genes_groups_heatmap(adata, groupby='group_new', 
                                      n_genes=top_genes, show=False,
                                      cmap='RdBu_r', dendrogram=True,
                                      swap_axes=True, figsize=(12, 8))
        plt.title(f'Top {top_genes} Differential Genes per Group')
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"热图生成失败: {e}")
        # 尝试简化版本的热图
        try:
            # 使用所有细胞的前50个高变基因
            sc.pp.highly_variable_genes(adata, n_top_genes=50)
            adata_hvg = adata[:, adata.var.highly_variable]
            
            sc.pl.heatmap(adata_hvg, groupby='group_new', 
                         var_names=adata_hvg.var_names[:30], 
                         cmap='RdBu_r', figsize=(10, 6), show=False)
            plt.title('Top 30 Highly Variable Genes')
            plt.tight_layout()
            plt.show()
        except Exception as e2:
            print(f"简化热图也失败: {e2}")

def create_dotplot_visualization(adata, top_genes=10):
    """创建点图展示标记基因"""
    print("\n=== 生成点图 ===")
    
    try:
        sc.pl.rank_genes_groups_dotplot(adata, groupby='group_new', 
                                      n_genes=top_genes, show=False,
                                      cmap='RdBu_r', figsize=(12, 6))
        plt.title(f'Top {top_genes} Marker Genes per Group')
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"点图生成失败: {e}")

def create_stacked_barplot(adata):
    """创建堆叠柱状图展示细胞类型组成"""
    print("\n=== 生成细胞类型组成图 ===")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 1. 细胞类型在每个分组中的比例
    cross_tab = pd.crosstab(adata.obs['group_new'], adata.obs['celltype'], normalize='index')
    cross_tab.plot(kind='bar', stacked=True, ax=ax1, colormap='tab20')
    ax1.set_title('细胞类型在各分组中的比例')
    ax1.set_ylabel('比例')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. 每个细胞类型在不同分组中的分布
    cross_tab2 = pd.crosstab(adata.obs['celltype'], adata.obs['group_new'], normalize='index')
    cross_tab2.plot(kind='bar', stacked=True, ax=ax2, colormap='Set3')
    ax2.set_title('各细胞类型在不同分组中的分布')
    ax2.set_ylabel('比例')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def create_celltype_specific_plots(celltype_results):
    """为每种细胞类型创建特异性分析图"""
    print("\n=== 生成细胞类型特异性分析图 ===")
    
    for celltype, adata_subset in list(celltype_results.items())[:4]:  # 只显示前4个
        try:
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            axes = axes.flat
            
            # 1. 热图
            sc.pl.rank_genes_groups_heatmap(adata_subset, groupby='group_new', 
                                          n_genes=10, show=False, ax=axes[0])
            axes[0].set_title(f'{celltype} - Top 10 DEGs')
            
            # 2. 点图
            sc.pl.rank_genes_groups_dotplot(adata_subset, groupby='group_new', 
                                          n_genes=8, show=False, ax=axes[1])
            axes[1].set_title(f'{celltype} - Marker Genes')
            
            # 3. 矩阵图
            sc.pl.rank_genes_groups_matrixplot(adata_subset, groupby='group_new', 
                                             n_genes=5, show=False, ax=axes[2])
            axes[2].set_title(f'{celltype} - Gene Expression Matrix')
            
            # 4. 堆叠Violin图
            sc.pl.rank_genes_groups_stacked_violin(adata_subset, groupby='group_new', 
                                                 n_genes=5, show=False, ax=axes[3])
            axes[3].set_title(f'{celltype} - Expression Distribution')
            
            plt.suptitle(f'细胞类型: {celltype}', fontsize=16, y=1.02)
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            print(f"生成 {celltype} 图表失败: {e}")
            continue

def export_results(adata, celltype_results, output_dir='diff_analysis_results'):
    """导出分析结果"""
    print(f"\n=== 导出结果到目录: {output_dir} ===")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 导出差异基因结果
    try:
        deg_results = []
        groups = adata.obs['group_new'].cat.categories
        
        for group in groups:
            group_df = sc.get.rank_genes_groups_df(adata, group=group)
            group_df['comparison_group'] = group
            deg_results.append(group_df)
        
        all_deg = pd.concat(deg_results, ignore_index=True)
        all_deg.to_csv(f'{output_dir}/all_differential_genes.csv', index=False)
        print("1. 差异基因结果已导出")
    except Exception as e:
        print(f"导出差异基因失败: {e}")
    
    # 2. 导出细胞类型组成表
    try:
        composition_table = pd.crosstab(adata.obs['group_new'], adata.obs['celltype'])
        composition_table.to_csv(f'{output_dir}/celltype_composition.csv')
        print("2. 细胞类型组成表已导出")
    except Exception as e:
        print(f"导出组成表失败: {e}")
    
    # 3. 导出细胞类型特异性结果
    try:
        for celltype, adata_subset in celltype_results.items():
            celltype_deg = []
            groups = adata_subset.obs['group_new'].cat.categories
            
            for group in groups:
                group_df = sc.get.rank_genes_groups_df(adata_subset, group=group)
                group_df['celltype'] = celltype
                group_df['comparison_group'] = group
                celltype_deg.append(group_df)
            
            if celltype_deg:
                celltype_all_deg = pd.concat(celltype_deg, ignore_index=True)
                safe_celltype_name = "".join(c for c in celltype if c.isalnum() or c in (' ', '-', '_')).rstrip()
                celltype_all_deg.to_csv(f'{output_dir}/deg_{safe_celltype_name}.csv', index=False)
        
        print("3. 细胞类型特异性结果已导出")
    except Exception as e:
        print(f"导出细胞类型结果失败: {e}")
    
    # 4. 保存汇总统计
    try:
        summary_stats = {
            'total_cells': adata.n_obs,
            'total_genes': adata.n_vars,
            'n_celltypes': adata.obs['celltype'].nunique(),
            'n_groups': adata.obs['group_new'].nunique(),
            'analysis_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        summary_df = pd.DataFrame(list(summary_stats.items()), columns=['Metric', 'Value'])
        summary_df.to_csv(f'{output_dir}/analysis_summary.csv', index=False)
        print("4. 分析摘要已导出")
    except Exception as e:
        print(f"导出摘要失败: {e}")
    
    print("导出完成!")

def main():
    """主执行函数"""
    # 设置文件路径
    h5ad_file = "merged_single_cell_data_with_group_new.h5ad"
    
    # 1. 加载和检查数据
    adata = load_and_check_data(h5ad_file)
    if adata is None:
        return
    
    # 2. 检查细胞类型覆盖情况
    coverage_table = check_celltype_coverage(adata)
    
    # 3. 执行全局差异分析
    adata = perform_global_differential_analysis(adata)
    
    # 4. 执行细胞类型特异性分析
    celltype_results = perform_celltype_specific_analysis(adata)
    
    # 5. 创建可视化
    create_volcano_plots(adata)
    create_heatmap_visualization(adata)
    create_dotplot_visualization(adata)
    create_stacked_barplot(adata)
    
    # 6. 创建细胞类型特异性图表
    if celltype_results:
        create_celltype_specific_plots(celltype_results)
    
    # 7. 导出结果
    export_results(adata, celltype_results)
    
    print("\n" + "="*60)
    print("分析完成!")
    print("="*60)
    print("关键发现:")
    print(f"- 分析了 {adata.n_obs} 个细胞中的 {adata.n_vars} 个基因")
    print(f"- 覆盖 {adata.obs['celltype'].nunique()} 种细胞类型")
    print(f"- 比较 {adata.obs['group_new'].nunique()} 个分组")
    
    # 显示显著的差异基因
    try:
        groups = adata.obs['group_new'].cat.categories
        if len(groups) > 1:
            print(f"\n显著的差异基因示例 (分组: {groups[1]} vs {groups[0]}):")
            deg_df = sc.get.rank_genes_groups_df(adata, group=groups[1])
            sig_genes = deg_df[(deg_df['pvals_adj'] < 0.05) & (np.abs(deg_df['scores']) > 1)]
            if len(sig_genes) > 0:
                print(sig_genes.head(10)[['names', 'scores', 'pvals_adj']])
    except:
        pass

if __name__ == "__main__":
    main()
