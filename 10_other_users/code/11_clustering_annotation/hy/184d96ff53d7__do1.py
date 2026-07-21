import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def load_and_process_data(file_path):
    """加载并处理单细胞数据"""
    print("正在加载数据...")
    adata = sc.read_h5ad(file_path)
    print(f"数据加载成功！细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    
    # 检查必要的列是否存在
    required_cols = ['sample', 'celltype', 'TIME']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"警告: 缺少以下列: {missing_cols}")
        print(f"可用列: {list(adata.obs.columns)}")
        # 尝试使用可能的替代列名
        if 'sample' in adata.obs.columns and 'sample' not in adata.obs.columns:
            adata.obs['sample'] = adata.obs['sample']
        if 'cell_type' in adata.obs.columns and 'celltype' not in adata.obs.columns:
            adata.obs['celltype'] = adata.obs['cell_type']
        if 'treatment' in adata.obs.columns and 'TIME' not in adata.obs.columns:
            adata.obs['TIME'] = adata.obs['treatment']
    
    return adata

def calculate_celltype_proportions(adata, sample_col='sample', 
                                 celltype_col='celltype', treatment_col='TIME'):
    """计算每种细胞类型在每个样本中的比例[1,3](@ref)"""
    print("正在计算细胞类型比例...")
    
    # 计算每个样本中每种细胞类型的细胞数量[1](@ref)
    celltype_counts = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col])
    
    # 计算每个样本的细胞总数
    total_cells_per_sample = celltype_counts.sum(axis=1)
    
    # 计算每种细胞类型在每个样本中的比例[3](@ref)
    celltype_proportions = celltype_counts.div(total_cells_per_sample, axis=0)
    
    # 重置索引并转换为长格式
    celltype_proportions = celltype_proportions.reset_index()
    celltype_proportions.rename(columns={'index': sample_col}, inplace=True)
    
    # 转换为长格式[3](@ref)
    proportion_long = pd.melt(
        celltype_proportions, 
        id_vars=[sample_col], 
        value_name='Proportion', 
        var_name='CellType'
    )
    
    # 添加治疗组信息
    sample_meta = adata.obs[[sample_col, treatment_col]].drop_duplicates()
    proportion_long = proportion_long.merge(sample_meta, on=sample_col)
    
    print(f"计算完成！共 {len(proportion_long)} 个样本-细胞类型组合")
    print(f"细胞类型: {proportion_long['CellType'].unique().tolist()}")
    print(f"治疗组: {proportion_long[treatment_col].unique().tolist()}")
    
    return proportion_long

def plot_all_celltype_boxplots(proportion_data, treatment_col='TIME', 
                             figsize=(20, 15), palette="Set2"):
    """绘制所有细胞类型的箱线图[6,7](@ref)"""
    
    celltypes = proportion_data['CellType'].unique()
    n_types = len(celltypes)
    
    print(f"开始绘制 {n_types} 种细胞类型的箱线图...")
    
    # 计算子图布局
    n_cols = min(4, n_types)
    n_rows = (n_types + n_cols - 1) // n_cols
    
    # 创建图形
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    # 处理子图数组格式
    if n_types == 1:
        axes = np.array([axes])
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # 为每个细胞类型绘制箱线图
    for idx, celltype in enumerate(celltypes):
        row = idx // n_cols
        col = idx % n_cols
        
        # 获取当前子图的轴
        if n_rows == 1 and n_cols == 1:
            ax = axes[0]
        elif n_rows == 1:
            ax = axes[col]
        elif n_cols == 1:
            ax = axes[row]
        else:
            ax = axes[row, col]
        
        # 筛选当前细胞类型的数据
        celltype_data = proportion_data[proportion_data['CellType'] == celltype]
        
        # 绘制箱线图[6,7](@ref)
        sns.boxplot(data=celltype_data, x=treatment_col, y='Proportion', 
                   ax=ax, palette=palette, showfliers=True, width=0.6)
        
        # 添加数据点[7](@ref)
        sns.stripplot(data=celltype_data, x=treatment_col, y='Proportion', 
                     ax=ax, color='black', size=3, jitter=0.2, alpha=0.6)
        
        # 设置标题和标签
        ax.set_title(f'{celltype}', fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel('Treatment Group', fontsize=10)
        ax.set_ylabel('Cell Proportion', fontsize=10)
        
        # 旋转x轴标签
        ax.tick_params(axis='x', rotation=45, labelsize=9)
        ax.tick_params(axis='y', labelsize=9)
        
        # 添加网格
        ax.grid(True, alpha=0.3)
    
    # 隐藏多余的子图
    for idx in range(n_types, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        if n_rows == 1 and n_cols == 1:
            fig.delaxes(axes[0])
        elif n_rows == 1:
            fig.delaxes(axes[col])
        elif n_cols == 1:
            fig.delaxes(axes[row])
        else:
            if row < axes.shape[0] and col < axes.shape[1]:
                fig.delaxes(axes[row, col])
    
    plt.tight_layout()
    plt.suptitle(f'Cell Type Proportions by {treatment_col} Group', 
                y=1.02, fontsize=16, fontweight='bold')
    
    return fig

def create_combined_boxplot(proportion_data, treatment_col='TIME', 
                          figsize=(15, 8)):
    """创建所有细胞类型的组合箱线图（一个图中显示）[4](@ref)"""
    plt.figure(figsize=figsize)
    
    # 按细胞类型和中位数排序
    median_proportions = proportion_data.groupby('CellType')['Proportion'].median().sort_values(ascending=False)
    proportion_data['CellType'] = pd.Categorical(proportion_data['CellType'], 
                                               categories=median_proportions.index)
    
    ax = sns.boxplot(data=proportion_data, x='CellType', y='Proportion', 
                    hue=treatment_col, palette="Set2", showfliers=True)
    
    # 添加数据点
    sns.stripplot(data=proportion_data, x='CellType', y='Proportion', 
                 hue=treatment_col, dodge=True, alpha=0.5, size=3)
    
    plt.title(f'All Cell Type Proportions by {treatment_col}', fontsize=14, fontweight='bold')
    plt.xlabel('Cell Type', fontsize=12)
    plt.ylabel('Proportion', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.legend(title=treatment_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    return plt.gcf()

def main():
    """主函数"""
    # 文件路径
    file_path = 'extracted_specific_patients.h5ad'
    
    try:
        # 1. 加载数据
        adata = load_and_process_data(file_path)
        
        # 2. 计算细胞类型比例
        proportion_data = calculate_celltype_proportions(adata)
        
        if proportion_data is None or len(proportion_data) == 0:
            print("错误: 无法计算细胞类型比例")
            return
        
        # 3. 绘制分面箱线图（所有细胞类型分开）
        print("\n绘制分面箱线图...")
        fig1 = plot_all_celltype_boxplots(proportion_data)
        plt.savefig('all_celltype_boxplots_faceted.png', dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        
        # 4. 绘制组合箱线图（所有细胞类型在一个图中）
        print("绘制组合箱线图...")
        fig2 = create_combined_boxplot(proportion_data)
        plt.savefig('all_celltype_boxplots_combined.png', dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        
        # 5. 保存数据
        proportion_data.to_csv('celltype_proportions_all.csv', index=False)
        
        # 显示图形
        plt.show()
        
        print("\n" + "="*50)
        print("分析完成！")
        print("生成的文件:")
        print("- all_celltype_boxplots_faceted.png: 分面箱线图（每个细胞类型单独显示）")
        print("- all_celltype_boxplots_combined.png: 组合箱线图（所有细胞类型在一个图中）")
        print("- celltype_proportions_all.csv: 细胞类型比例数据")
        print("="*50)
        
        # 显示数据摘要
        celltype_summary = proportion_data.groupby(['CellType', 'TIME'])['Proportion'].agg(['mean', 'std', 'count']).round(4)
        print("\n细胞类型比例摘要:")
        print(celltype_summary)
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
