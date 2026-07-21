import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.path import Path
import matplotlib.patches as patches

# 加载h5ad数据
adata1 = sc.read_h5ad("integrated_data.h5ad")
target_celltypes = ["\u200bCTMCs","SOCS3+ CTMCs","FAMs","M2","CD4+ Tem", "CD4+ Treg", "CD8+ Tex", "M2", "iDC","SIGLEC1+ Mø","SOCS3+ M1","MACRO+ Mø","cDC1","cDC2","CD4+ CD247+ T"]

# 筛选目标细胞类型
mask = adata1.obs['celltype_final'].isin(target_celltypes)
adata = adata1[mask].copy()




def plot_celltype_proportions(adata, group_col='group', celltype_col='celltype_final'):
    """
    绘制细胞类型比例堆积图，添加带状连接线
    
    参数:
    adata - AnnData对象
    group_col - 分组列名 (obs中的列)
    celltype_col - 细胞类型列名 (obs中的列)
    """
    # 1. 计算各分组中不同细胞类型的比例
    celltype_counts = adata.obs.groupby([group_col, celltype_col]).size().unstack()
    proportions = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    proportions = proportions.fillna(0)
    
    # 对分组和细胞类型进行排序（如需要）
    groups = sorted(proportions.index)  # 分组排序
    celltypes = proportions.columns.tolist()  # 保持原始顺序或手动排序
    
    # 2. 准备绘图
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 为每个细胞类型生成颜色（最多20种）
    n_celltypes = len(celltypes)
    colors = plt.cm.tab20(np.linspace(0, 1, min(n_celltypes, 20)))
    if n_celltypes > 20:
        # 补充额外颜色
        extra_colors = plt.cm.jet(np.linspace(0, 1, n_celltypes - 20))
        colors = np.vstack([colors, extra_colors])
    
    # 3. 绘制堆积柱状图
    bar_width = 0.8  # 柱状图宽度
    bottoms = np.zeros(len(groups))
    
    for i, celltype in enumerate(celltypes):
        # 获取当前细胞类型在各组的比例
        props = proportions.loc[groups, celltype].values
        
        # 绘制堆积柱
        bars = ax.bar(
            np.arange(len(groups)),
            props,
            bottom=bottoms,
            color=colors[i],
            edgecolor='white',
            width=bar_width,
            label=celltype
        )
        
        # 计算带状连接点的位置
        bar_positions = np.arange(len(groups))
        bar_heights = bottoms + props/2  # 中点位置
        
        # 创建带状连接路径
        for j in range(len(groups)-1):
            # 创建两点之间的带状连接
            path = Path([
                (bar_positions[j], bar_heights[j]),
                (bar_positions[j+1], bar_heights[j+1]),
                (bar_positions[j+1], bar_heights[j+1] - props[j+1]/2),
                (bar_positions[j], bar_heights[j] - props[j]/2),
                (bar_positions[j], bar_heights[j])
            ])
            
            # 添加半透明填充
            patch = patches.PathPatch(
                path,
                facecolor=colors[i],
                alpha=0.25,  # 半透明效果
                lw=0
            )
            ax.add_patch(patch)
        
        # 更新底部位置
        bottoms += props
    
    # 4. 设置图表样式
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, fontsize=12)
    ax.set_ylabel('Cell Proportion(%)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # 添加百分比的y轴刻度
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=10)
    
    # 5. 添加图例
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles, 
        labels, 
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        frameon=False,
        fontsize=10,
        title=celltype_col,
        title_fontsize=12
    )
    
    # 6. 优化布局
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)  # 为图例留出空间
    
    plt.savefig('celltype_proportions_with_bands.png', dpi=300, bbox_inches='tight')
    plt.show()

# 使用函数绘制图表
plot_celltype_proportions(adata)
