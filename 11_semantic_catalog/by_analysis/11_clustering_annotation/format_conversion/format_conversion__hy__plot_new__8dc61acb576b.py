import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib样式，使图形更适合发表
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'Arial',  # 使用Arial字体，适合发表
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

def plot_state_spatial_distribution(adata, output_dir, figsize=(8, 8), 
                                    dot_size=1, color_palette='tab20'):
    """
    绘制每个样本的state空间分布图
    
    参数:
    adata: AnnData对象，包含空间转录组数据
    output_dir: 输出目录路径
    figsize: 图形大小
    dot_size: 点的大小
    color_palette: 颜色调色板
    """
    
    print("=" * 60)
    print("开始绘制每个样本的state空间分布图")
    print("=" * 60)
    
    # 检查必要的列
    required_cols = ['sample', 'state']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"adata.obs中缺少必需的列: {col}")
    
    # 检查是否有空间坐标
    if 'spatial' not in adata.obsm_keys():
        raise ValueError("adata.obsm中缺少'spatial'坐标")
    
    # 获取所有样本
    library_names = adata.obs['sample'].drop_duplicates().tolist()
    print(f"找到 {len(library_names)} 个样本: {library_names}")
    
    # 获取所有state
    all_states = adata.obs['state'].dropna().unique()
    print(f"找到 {len(all_states)} 个state: {all_states}")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 创建颜色映射
    if color_palette == 'tab20':
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.get_cmap(color_palette)
    
    # 为每个state创建固定的颜色映射
    state_colors = {}
    for i, state in enumerate(all_states):
        if pd.isna(state):
            state_colors['NaN'] = '#808080'  # 灰色表示NaN
        else:
            color_idx = i % cmap.N
            state_colors[state] = cmap(color_idx)
    
    print(f"为state创建了颜色映射")
    
    # 为每个样本绘制空间图
    for library in library_names:
        print(f"处理样本: {library}")
        
        # 筛选当前样本的数据
        sample_data = adata[adata.obs['sample'] == library, :].copy()
        
        if sample_data.n_obs == 0:
            print(f"  警告: 样本 {library} 没有细胞，跳过")
            continue
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        
        # 获取空间坐标
        spatial_coords = sample_data.obsm['spatial']
        
        # 获取当前样本的state
        states = sample_data.obs['state']
        
        # 为每个state绘制散点
        for state, color in state_colors.items():
            if state == 'NaN':
                # NaN值的处理
                mask = states.isna()
            else:
                mask = states == state
            
            if mask.sum() > 0:
                coords = spatial_coords[mask]
                ax.scatter(coords[:, 0], coords[:, 1], 
                          c=color, s=dot_size, label=state, alpha=0.8, edgecolor='none')
        
        # 设置图形属性
        ax.set_title(f'{library} - State Distribution', fontsize=14, fontweight='bold')
        ax.set_xlabel('X coordinate', fontsize=12)
        ax.set_ylabel('Y coordinate', fontsize=12)
        
        # 反转y轴，使原点在左上角（空间转录组常见）
        ax.invert_yaxis()
        
        # 添加图例
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                 fontsize='small', title='State', title_fontsize='small',
                 markerscale=5)  # 放大图例中的点
        
        # 添加网格
        ax.grid(alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # 保存图形
        png_file = os.path.join(output_dir, f'spatial_state_{library}.png')
        pdf_file = os.path.join(output_dir, f'spatial_state_{library}.pdf')
        
        plt.savefig(png_file, facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(pdf_file, facecolor='white', bbox_inches="tight")
        plt.close()
        
        print(f"  已保存: {png_file}")
        print(f"  已保存: {pdf_file}")
        
        # 为当前样本创建统计信息
        create_sample_statistics(sample_data, library, output_dir)
    
    # 创建所有样本的汇总图
    create_summary_plots(adata, output_dir, state_colors)
    
    print(f"\n所有图形已保存到: {output_dir}")

def create_sample_statistics(sample_data, sample_name, output_dir):
    """为单个样本创建统计信息"""
    
    stats_file = os.path.join(output_dir, f'state_statistics_{sample_name}.txt')
    
    with open(stats_file, 'w') as f:
        f.write(f"Sample: {sample_name}\n")
        f.write("=" * 40 + "\n\n")
        
        f.write(f"Total cells: {sample_data.n_obs}\n")
        f.write(f"Total genes: {sample_data.n_vars}\n\n")
        
        f.write("State distribution:\n")
        state_counts = sample_data.obs['state'].value_counts(dropna=False)
        
        for state, count in state_counts.items():
            percentage = count / sample_data.n_obs * 100
            state_display = 'NaN' if pd.isna(state) else state
            f.write(f"  {state_display}: {count} cells ({percentage:.2f}%)\n")

def create_summary_plots(adata, output_dir, state_colors):
    """创建汇总图形"""
    
    # 1. 所有样本的state分布条形图
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # 计算每个样本的state分布
    sample_state_counts = adata.obs.groupby(['sample', 'state']).size().unstack(fill_value=0)
    
    # 转换为比例
    sample_state_proportions = sample_state_counts.div(sample_state_counts.sum(axis=1), axis=0)
    
    # 绘制堆叠条形图
    bottom = np.zeros(len(sample_state_proportions))
    
    for state in sample_state_proportions.columns:
        if pd.isna(state):
            color = '#808080'
            label = 'NaN'
        else:
            color = state_colors.get(state, '#808080')
            label = state
        
        ax.bar(range(len(sample_state_proportions)), 
               sample_state_proportions[state], 
               bottom=bottom, 
               color=color, 
               label=label)
        bottom += sample_state_proportions[state]
    
    # 设置图形属性
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title('State Distribution Across Samples', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(sample_state_proportions)))
    ax.set_xticklabels(sample_state_proportions.index, rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='State')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # 保存图形
    summary_file = os.path.join(output_dir, 'state_distribution_summary.png')
    plt.savefig(summary_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  已保存汇总图: {summary_file}")
    
    # 2. 创建state的热力图
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # 创建热图数据
    heatmap_data = sample_state_proportions.T
    
    # 绘制热图
    sns.heatmap(heatmap_data, cmap='YlOrRd', annot=False, 
                cbar_kws={'label': 'Proportion'}, ax=ax)
    
    ax.set_title('State Proportions by Sample', fontsize=14, fontweight='bold')
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('State', fontsize=12)
    
    plt.tight_layout()
    
    # 保存热图
    heatmap_file = os.path.join(output_dir, 'state_heatmap.png')
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  已保存热图: {heatmap_file}")

def create_individual_spatial_plots(adata, output_dir, sample_list=None, 
                                    figsize=(8, 8), dot_size=1):
    """
    为指定样本创建单独的空间图
    """
    
    if sample_list is None:
        sample_list = adata.obs['sample'].drop_duplicates().tolist()
    
    for sample in sample_list:
        print(f"处理样本: {sample}")
        
        # 筛选样本
        sample_data = adata[adata.obs['sample'] == sample, :].copy()
        
        if sample_data.n_obs == 0:
            print(f"  警告: 样本 {sample} 没有细胞，跳过")
            continue
        
        # 创建图形
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        
        # 1. 按state着色的空间图
        ax1 = axes[0]
        
        # 获取空间坐标
        spatial_coords = sample_data.obsm['spatial']
        
        # 获取state
        states = sample_data.obs['state']
        
        # 创建颜色映射
        unique_states = states.dropna().unique()
        cmap = plt.cm.tab20
        state_colors = {}
        
        for i, state in enumerate(unique_states):
            color_idx = i % 20
            state_colors[state] = cmap(color_idx)
        
        # 绘制每个state
        for state, color in state_colors.items():
            mask = states == state
            if mask.sum() > 0:
                coords = spatial_coords[mask]
                ax1.scatter(coords[:, 0], coords[:, 1], 
                          c=color, s=dot_size, label=state, alpha=0.8, edgecolor='none')
        
        # NaN值
        nan_mask = states.isna()
        if nan_mask.sum() > 0:
            coords = spatial_coords[nan_mask]
            ax1.scatter(coords[:, 0], coords[:, 1], 
                      c='#808080', s=dot_size, label='NaN', alpha=0.8, edgecolor='none')
        
        ax1.set_title(f'{sample} - State Distribution', fontsize=14, fontweight='bold')
        ax1.set_xlabel('X coordinate', fontsize=12)
        ax1.set_ylabel('Y coordinate', fontsize=12)
        ax1.invert_yaxis()
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                  fontsize='small', title='State')
        ax1.grid(alpha=0.3, linestyle='--')
        
        # 2. 按state分面图
        ax2 = axes[1]
        
        # 获取前几个主要的state
        top_states = states.value_counts().head(5).index.tolist()
        
        for i, state in enumerate(top_states):
            if pd.isna(state):
                continue
            
            mask = states == state
            if mask.sum() > 0:
                coords = spatial_coords[mask]
                ax2.scatter(coords[:, 0], coords[:, 1], 
                          c=state_colors.get(state, '#808080'), 
                          s=dot_size, label=state, alpha=0.8, edgecolor='none')
        
        ax2.set_title(f'{sample} - Top States', fontsize=14, fontweight='bold')
        ax2.set_xlabel('X coordinate', fontsize=12)
        ax2.set_ylabel('Y coordinate', fontsize=12)
        ax2.invert_yaxis()
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                  fontsize='small', title='State')
        ax2.grid(alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # 保存图形
        output_file = os.path.join(output_dir, f'{sample}_state_spatial.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  已保存: {output_file}")

def main():
    # 设置路径
    input_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state.h5ad"
    output_dir = "./state_spatial_plots_new"
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("State空间分布图绘制脚本")
    print("=" * 60)
    
    # 1. 读取数据
    print("\n1. 读取数据...")
    adata = sc.read_h5ad(input_file)
    print(f"   数据形状: {adata.shape}")
    print(f"   观测变量(obs)列: {list(adata.obs.columns)}")
    
    # 检查必要的列
    required_cols = ['sample', 'state']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"adata.obs中缺少必需的列: {col}")
    
    # 检查state列是否有数据
    if adata.obs['state'].isna().all():
        print("警告: state列全部为NaN，请先运行state匹配脚本")
        return
    
    # 统计信息
    print(f"\n2. 数据统计:")
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   总基因数: {adata.n_vars}")
    print(f"   样本数: {adata.obs['sample'].nunique()}")
    print(f"   state类型数: {adata.obs['state'].nunique()}")
    
    # 显示state分布
    state_counts = adata.obs['state'].value_counts(dropna=False)
    print(f"\n   state分布:")
    for state, count in state_counts.items():
        percentage = count / adata.n_obs * 100
        state_display = 'NaN' if pd.isna(state) else state
        print(f"     {state_display}: {count} ({percentage:.2f}%)")
    
    # 3. 绘制空间图
    print(f"\n3. 绘制空间分布图...")
    
    # 方法1: 使用自定义函数
    plot_state_spatial_distribution(
        adata, 
        output_dir=output_dir,
        figsize=(10, 8),
        dot_size=1,
        color_palette='tab20'
    )
    
    # 方法2: 创建单独的子目录用于分面图
    print(f"\n4. 创建详细分面图...")
    detailed_dir = os.path.join(output_dir, "detailed_plots")
    os.makedirs(detailed_dir, exist_ok=True)
    
    create_individual_spatial_plots(
        adata,
        output_dir=detailed_dir,
        figsize=(12, 8),
        dot_size=1
    )
    
    # 5. 创建分析报告
    print(f"\n5. 创建分析报告...")
    create_analysis_report(adata, output_dir)
    
    print("\n" + "=" * 60)
    print("绘图完成!")
    print(f"所有结果已保存到: {output_dir}")
    print("=" * 60)

def create_analysis_report(adata, output_dir):
    """创建详细的分析报告"""
    
    report_file = os.path.join(output_dir, "state_spatial_analysis_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("State空间分布分析报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("1. 数据基本信息\n")
        f.write(f"   总细胞数: {adata.n_obs}\n")
        f.write(f"   总基因数: {adata.n_vars}\n")
        f.write(f"   样本数: {adata.obs['sample'].nunique()}\n")
        f.write(f"   样本列表: {list(adata.obs['sample'].unique())}\n\n")
        
        f.write("2. State分布统计\n")
        state_counts = adata.obs['state'].value_counts(dropna=False)
        for state, count in state_counts.items():
            percentage = count / adata.n_obs * 100
            state_display = 'NaN' if pd.isna(state) else state
            f.write(f"   {state_display}: {count} ({percentage:.2f}%)\n")
        f.write("\n")
        
        f.write("3. 各样本State分布\n")
        for sample in adata.obs['sample'].unique():
            sample_data = adata[adata.obs['sample'] == sample, :]
            sample_state_counts = sample_data.obs['state'].value_counts(dropna=False)
            
            f.write(f"   样本: {sample} ({sample_data.n_obs} 个细胞)\n")
            for state, count in sample_state_counts.items():
                percentage = count / sample_data.n_obs * 100
                state_display = 'NaN' if pd.isna(state) else state
                f.write(f"     {state_display}: {count} ({percentage:.2f}%)\n")
            f.write("\n")
    
    print(f"   分析报告已保存到: {report_file}")

if __name__ == "__main__":
    main()
