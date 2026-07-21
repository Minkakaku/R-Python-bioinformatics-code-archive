import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

# 加载h5ad数据

adata1 = sc.read_h5ad("integrated_data.h5ad")
target_celltypes = ["\u200bCTMCs","SOCS3+ CTMCs","FAMs","M2","CD4+ Tem", "CD4+ Treg", "CD8+ Tex", "M2", "iDC","SIGLEC1+ Mø","SOCS3+ M1","MACRO+ Mø","cDC1","cDC2","CD4+ CD247+ T"]

# 筛选目标细胞类型
mask = adata1.obs['celltype_final'].isin(target_celltypes)
adata = adata1[mask].copy()



def plot_celltype_proportions(adata, group_col='group', celltype_col='celltype_final'):
    """
    绘制细胞类型比例堆积图，添加精确的端点连接线
    
    参数:
    adata - AnnData对象
    group_col - 分组列名 (obs中的列)
    celltype_col - 细胞类型列名 (obs中的列)
    """
    # 1. 计算各分组中不同细胞类型的比例
    celltype_counts = adata.obs.groupby([group_col, celltype_col]).size().unstack()
    proportions = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    proportions = proportions.fillna(0)
    
    # 对分组进行排序（确保疾病组在前，正常组在后）
    groups = sorted(proportions.index, reverse=False)  # 疾病组:0，正常组:1
    
    # 2. 准备绘图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 细胞类型和颜色映射
    celltypes = proportions.columns.tolist()
    n_celltypes = len(celltypes)
    colors = plt.cm.tab20(np.linspace(0, 1, min(n_celltypes, 20)))
    
    # 存储每个细胞类型在每组中的位置信息
    group_positions = {group: [] for group in groups}
    
    # 3. 绘制堆积柱状图并记录位置信息
    bar_width = 0.6  # 柱状图宽度
    group_spacing = 0.8  # 两组之间的间距
    group_centers = [0, group_spacing]  # 每组中心位置
    
    # 用于存储每组细胞类型的边界点（用于连接线）
    boundary_points = {group: {"top": [], "bottom": []} for group in groups}
    
    for group_idx, group in enumerate(groups):
        bottom = 0
        group_center = group_centers[group_idx]
        
        for celltype_idx, celltype in enumerate(celltypes):
            prop = proportions.loc[group, celltype]
            
            # 记录位置信息
            left = group_center - bar_width/2
            right = group_center + bar_width/2
            bottom_pos = bottom
            top_pos = bottom + prop
            
            # 记录边界点（在柱子的右侧或左侧）
            if group_idx == 0:  # 疾病组（记录右侧边界）
                boundary_points[group]["top"].append((right, top_pos))
                boundary_points[group]["bottom"].append((right, bottom_pos))
            else:  # 正常组（记录左侧边界）
                boundary_points[group]["top"].append((left, top_pos))
                boundary_points[group]["bottom"].append((left, bottom_pos))
            
            # 绘制矩形色块
            rect = patches.Rectangle(
                (left, bottom_pos), 
                bar_width, 
                prop,
                facecolor=colors[celltype_idx],
                edgecolor='white',
                linewidth=1
            )
            ax.add_patch(rect)
            
            # 更新底部位置
            bottom += prop
    
    # 4. 绘制端点连接线
    for celltype_idx in range(n_celltypes):
        # 获取疾病组色块右侧边界点（右上角和右下角）
        disease_right_top = boundary_points[groups[0]]["top"][celltype_idx]
        disease_right_bottom = boundary_points[groups[0]]["bottom"][celltype_idx]
        
        # 获取正常组色块左侧边界点（左上角和左下角）
        normal_left_top = boundary_points[groups[1]]["top"][celltype_idx]
        normal_left_bottom = boundary_points[groups[1]]["bottom"][celltype_idx]
        
        # 创建连接路径（四边形）
        path = Path([
            disease_right_top,      # 疾病组右上角
            disease_right_bottom,   # 疾病组右下角
            normal_left_bottom,     # 正常组左下角
            normal_left_top,        # 正常组左上角
            disease_right_top       # 闭合路径
        ])
        
        # 添加半透明填充
        patch = patches.PathPatch(
            path,
            facecolor=colors[celltype_idx],
            alpha=0.3,  # 半透明效果
            lw=1,
            edgecolor='none'
        )
        ax.add_patch(patch)
    
    # 5. 设置图表样式
    ax.set_xticks(group_centers)
    ax.set_xticklabels(groups, fontsize=12)
    ax.set_ylabel('Cell Proportion(%)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.set_xlim(group_centers[0]-0.5, group_centers[1]+0.5)
    
    # 添加百分比的y轴刻度
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=10)
    
    # 6. 添加图例
    handles = [patches.Patch(color=colors[i], label=celltypes[i]) for i in range(n_celltypes)]
    legend = ax.legend(
        handles=handles,
        loc='upper left',
        bbox_to_anchor=(1.02, 1),
        frameon=False,
        fontsize=10,
        title=celltype_col,
        title_fontsize=12
    )
    
    # 7. 优化布局
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)
    
    plt.savefig('celltype_proportions_with_precise_links.png', dpi=300, bbox_inches='tight')
    plt.show()

# 使用函数绘制图表
plot_celltype_proportions(adata)
