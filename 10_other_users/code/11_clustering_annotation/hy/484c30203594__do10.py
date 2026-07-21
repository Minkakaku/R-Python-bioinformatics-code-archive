import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import matplotlib as mpl

# 设置全局字体
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# 1. 加载数据并提取指定细胞类型
def load_and_filter_data():
    adata = sc.read_h5ad('integrated_data.h5ad')
    target_celltypes = ["CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", "iDC", 
                        "THY1+ FB", "Kera Spinous", "Kera Inflammatory", 
                        "CTMCs", "Disease Stem", "POSTN+ FB"]
    
    # 筛选目标细胞类型
    mask = adata.obs['celltype_final'].isin(target_celltypes)
    adata = adata[mask].copy()
    return adata

# 2. 计算不同距离范围内的细胞数量（使用group列）
def calculate_cell_counts(adata):
    # 获取POSTN+ FB细胞的坐标
    postn_mask = adata.obs['celltype_final'] == "POSTN+ FB"
    postn_coords = adata[postn_mask].obsm['spatial']
    
    # 获取其他细胞类型的坐标和元数据
    other_mask = ~postn_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype_final']
    other_group = adata[other_mask].obs['group'].map({'disease': 0, 'normal': 1})
    
    # 构建KDTree
    tree = KDTree(other_coords)
    
    # 定义距离阈值
    distances = [0, 10, 20, 50, 100, 150, 200]
    counts = np.zeros((len(distances)-1, len(target_celltypes)-1, 2))
    
    # 对每个POSTN+ FB细胞计算距离范围内的细胞数量
    for center in postn_coords:
        # 修复：使用query_ball_point替代query
        indices = tree.query_ball_point(center, r=200)  # 获取200μm内所有点的索引
        
        for idx in indices:
            pt = other_coords[idx]
            dist = np.linalg.norm(center - pt)  # 计算欧氏距离
            
            # 获取细胞类型和分组信息
            ct = other_celltypes.iloc[idx]
            grp = other_group.iloc[idx]
            ct_idx = target_celltypes.index(ct) if ct in target_celltypes else -1
            
            # 如果细胞类型有效且不是POSTN+ FB
            if 0 <= ct_idx < len(target_celltypes)-1:
                # 确定距离区间
                for i in range(len(distances)-1):
                    if distances[i] <= dist < distances[i+1]:
                        counts[i, ct_idx, grp] += 1
                        break
    
    # 计算平均数量
    mean_counts = counts / len(postn_coords)
    return mean_counts, distances



# 创建绘图
def create_radial_plot(mean_counts, distances, celltypes):
    # 设置图形
    plt.figure(figsize=(12, 12))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    
    # 删除网格和标签
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    
    # 3. 绘制同心圆
    ring_colors = ['#D6EAF8', '#AED6F1', '#85C1E9', '#5DADE2', '#3498DB', '#2E86C1']
    for i, r in enumerate(distances[1:]):
        ax.bar(np.linspace(0, 2*np.pi, 100), r, 
               width=2*np.pi, bottom=0, 
               color=ring_colors[i], alpha=0.6, zorder=1)
    
    # 4. 绘制分隔线 (10条扇形分隔线)
    for i in range(len(celltypes)):
        angle = i * 2*np.pi/len(celltypes)
        ax.plot([angle, angle], [0, 220], 
                color='#34495E', linestyle='-', linewidth=1, alpha=0.6)
    
    # 5. 添加细胞类型标签 (外圈水平放置)
    for i, ct in enumerate(celltypes):
        if ct == "POSTN+ FB":
            continue
            
        angle = (i + 0.5) * 2*np.pi/len(celltypes)
        x = angle
        y = 220  # 固定在200μm圆环外侧
        
        # 添加白色背景
        ax.text(x, y, ct, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', 
                         facecolor='white', 
                         edgecolor='none', 
                         alpha=0.7))
    
    # 6. 绘制气泡图（使用group分组）
    bubble_scales = [0.5, 1.0, 2.0, 3.5, 5.0, 7.0]  # 气泡大小比例
    offset_angles = [-0.08, 0.08]  # 疾病和正常气泡的偏移角度
    
    # 细胞类型和分组颜色映射
    cmap_dis = mpl.colormaps['Reds']
    cmap_norm = mpl.colormaps['Blues']
    
    # 遍历每个距离环、每种细胞类型
    for ring_idx in range(len(distances)-1):  # 6个距离环
        for ct_idx in range(len(celltypes)-1):  # 10种细胞类型
            if celltypes[ct_idx] == "POSTN+ FB":
                continue
                
            # 获取数量和角度
            counts = mean_counts[ring_idx, ct_idx]
            base_angle = (ct_idx + 0.5) * 2*np.pi/len(celltypes)
            radius = distances[ring_idx + 1] - 5  # 气泡位置在当前环的内侧
            
            # 分别绘制疾病和正常气泡
            for grp_idx, offset in enumerate(offset_angles):
                angle = base_angle + offset
                size = counts[grp_idx] * bubble_scales[ring_idx] * 10
                
                # 选择颜色映射
                cmap = cmap_dis if grp_idx == 0 else cmap_norm
                color = cmap(0.6 + grp_idx*0.2)
                
                # 气泡边框颜色
                edge_color = '#E74C3C' if grp_idx == 0 else '#3498DB'
                
                # 绘制气泡
                ax.scatter(angle, radius, s=size, 
                           facecolor=color, 
                           edgecolor=edge_color,
                           linewidth=1.5, 
                           alpha=0.8, 
                           zorder=10)
    
    # 设置显示范围
    ax.set_ylim(0, 250)
    
    # 添加距离标签
    label_positions = np.linspace(0, 2*np.pi, 8)[:-1]
    for dist, pos in zip(distances[1:], label_positions):
        ax.text(pos, dist+15, f"{dist}μm", 
                ha='center', va='center', 
                fontsize=8, fontweight='bold')
    
    # 添加分组图例
    ax.text(0, 240, "Disease", color='#E74C3C', ha='center', fontweight='bold')
    ax.text(np.pi, 240, "Normal", color='#3498DB', ha='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('radial_plot_group.png', dpi=300, bbox_inches='tight')
    plt.show()

# 主程序
if __name__ == "__main__":
    # 1. 加载并过滤数据
    adata = load_and_filter_data()
    
    # 2. 计算细胞数量统计（使用group列）
    target_celltypes = ["CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", "iDC", 
                        "THY1+ FB", "Kera Spinous", "Kera Inflammatory", 
                        "CTMCs", "Disease Stem", "POSTN+ FB"]
    mean_counts, distances = calculate_cell_counts(adata)
    
    # 3-6. 创建环形图
    create_radial_plot(mean_counts, distances, target_celltypes)
