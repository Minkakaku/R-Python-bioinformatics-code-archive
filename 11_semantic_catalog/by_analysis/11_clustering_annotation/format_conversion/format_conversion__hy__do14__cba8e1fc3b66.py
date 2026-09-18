import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import matplotlib as mpl
import math

# 设置全局字体
#plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# 1. 加载数据并提取指定细胞类型
def load_and_filter_data():
    adata = sc.read_h5ad('integrated_data.h5ad')
    target_celltypes = ["CD4+ Tem", "CD4+ Treg", "CD8+ Tex", "M2", "iDC", 
                        "THY1+ FB", "Kera Spinous", "Kera Inflammatory", 
                        "\u200bCTMCs", "Disease Stem", "POSTN+ FB",'CD4+ CD247+ T',"SBC Proliferation","SBC Spinous","SIGLEC1+ Mø","LAMP3+ DC","cDC1","cDC2"]
    
    # 筛选目标细胞类型
    mask = adata.obs['celltype_final'].isin(target_celltypes)
    adata = adata[mask].copy()
    return adata

# 2. 计算不同距离范围内的细胞数量和密度
def calculate_cell_counts(adata):
    # 获取POSTN+ FB细胞的坐标
    postn_mask = adata.obs['celltype_final'] == "CD4+ Tem"
    postn_coords = adata[postn_mask].obsm['spatial']
    num_centers = len(postn_coords)
    
    # 获取其他细胞类型的坐标和元数据
    other_mask = ~postn_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype_final']
    other_group = adata[other_mask].obs['group'].map({'disease': 0, 'normal': 1})
    
    # 构建KDTree
    tree = KDTree(other_coords)
    
    # 定义距离阈值
    distances = [0, 50, 100, 150, 200]
    counts = np.zeros((len(distances)-1, len(target_celltypes)-1, 2))
    
    # 对每个POSTN+ FB细胞计算距离范围内的细胞数量
    for center in postn_coords:
        # 使用query_ball_point获取200μm内所有点的索引
        indices = tree.query_ball_point(center, r=200)
        
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
    mean_counts = counts / num_centers
    
    # 计算扇形面积和细胞密度
    num_sectors = len(target_celltypes) - 1  # 扇形数量（排除POSTN+ FB）
    sector_areas = np.zeros((len(distances)-1, num_sectors))
    densities = np.zeros_like(mean_counts)
    
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        # 计算圆环面积
        ring_area = math.pi * (r2**2 - r1**2)
        # 计算单个扇形面积
        sector_area = ring_area / num_sectors
        
        for sector_idx in range(num_sectors):
            sector_areas[dist_idx, sector_idx] = sector_area
            # 计算密度（细胞数/面积）
            densities[dist_idx, sector_idx, 0] = mean_counts[dist_idx, sector_idx, 0] / sector_area
            densities[dist_idx, sector_idx, 1] = mean_counts[dist_idx, sector_idx, 1] / sector_area
    
    return mean_counts, densities, distances, sector_areas

# 导出计数结果到CSV文件（添加面积和密度列）
def export_counts_to_csv(mean_counts, densities, distances, sector_areas, celltypes, filename='cell_counts.csv'):
    # 创建距离范围标签
    dist_labels = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    
    # 准备数据
    data = []
    for ct_idx in range(len(celltypes)-1):  # 排除POSTN+ FB
        cell_type = celltypes[ct_idx]
        for dist_idx in range(len(dist_labels)):
            area = sector_areas[dist_idx, ct_idx]
            # 疾病组数据
            data.append([
                cell_type,
                dist_labels[dist_idx],
                'disease',
                mean_counts[dist_idx, ct_idx, 0],
                area,
                densities[dist_idx, ct_idx, 0]
            ])
            # 正常组数据
            data.append([
                cell_type,
                dist_labels[dist_idx],
                'normal',
                mean_counts[dist_idx, ct_idx, 1],
                area,
                densities[dist_idx, ct_idx, 1]
            ])
    
    # 创建DataFrame并保存为CSV
    df = pd.DataFrame(data, columns=[
        'celltype', 
        'distance_range', 
        'group', 
        'cell_count',
        'sector_area',
        'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"细胞计数已导出到 {filename}")

# 创建绘图
def create_radial_plot(mean_counts, densities, distances, celltypes):
    # 设置图形
    plt.figure(figsize=(12, 12))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    
    # 删除网格和标签
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    
    # 3. 绘制同心圆（无填充颜色）
    for i, r in enumerate(distances[1:]):  # 从第一个距离边界开始
        ax.plot(np.linspace(0, 2*np.pi, 100), 
                [r] * 100, 
                color='gray', 
                linestyle='-', 
                linewidth=1, 
                alpha=0.5)
    
    # 4. 绘制分隔线 (10条扇形分隔线)
    for i in range(len(celltypes)):
        angle = i * 2*np.pi/len(celltypes)
        ax.plot([angle, angle], [0, 230], 
                color='#34495E', linestyle='-', linewidth=1, alpha=0.6)
    
    # 5. 添加细胞类型标签 (外圈水平放置)
    for i, ct in enumerate(celltypes):
        if ct == "CD4+ Tem":
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
    # 关键修改1：使用细胞密度而非细胞数量
    # 关键修改2：根据距离层级使用不同的缩放因子
    offset_angles = [-0.08, 0.08]  # 疾病和正常气泡的偏移角度
    
    # 细胞类型和分组颜色映射
    cmap_dis = mpl.colormaps['Reds']
    cmap_norm = mpl.colormaps['Blues']
    
    # 遍历每个距离环、每种细胞类型
    for ring_idx in range(len(distances)-1):  # 6个距离环
        for ct_idx in range(len(celltypes)-1):  # 10种细胞类型
            if celltypes[ct_idx] == "CD4+ Tem":
                continue
                
            # 获取密度值（而非细胞数量）
            density_vals = densities[ring_idx, ct_idx]
            base_angle = (ct_idx + 0.5) * 2*np.pi/len(celltypes)
            # 气泡位置在当前环的内侧（距外边界5um）
            radius = distances[ring_idx + 1] - 5  
            
            # 分别绘制疾病和正常气泡
            for grp_idx, offset in enumerate(offset_angles):
                angle = base_angle + offset
                
                # 根据距离层级使用不同的缩放因子
                if ring_idx == 0:  # 0-10um层级
                    size = density_vals[grp_idx] * 20000
                elif ring_idx == 1:  # 10-20um层级
                    size = density_vals[grp_idx] * 20000
                else:  # 其他层级
                    size = density_vals[grp_idx] * 100000
                
                # 选择颜色映射
                cmap = cmap_dis if grp_idx == 0 else cmap_norm
                color = cmap(0.6 + grp_idx*0.2)
                
                # 关键修改3：绘制实心圆点（无边框）
                ax.scatter(angle, radius, s=size, 
                           facecolor=color, 
                           edgecolor='none',  # 关键修改：去除边框[1,2](@ref)
                           alpha=0.8, 
                           zorder=10,
                           marker='o')  # 明确指定为实心圆点[5](@ref)
    
    # 设置显示范围
    ax.set_ylim(0, 260)
    
    # 添加距离标签
    label_positions = np.linspace(0, 2*np.pi, 7)[:-1]
    for i, dist in enumerate(distances[1:]):
        # 使用标准标签格式
        label = f"{distances[i]}-{distances[i+1]}μm"
        ax.text(label_positions[i], dist+15, label, 
                ha='center', va='center', 
                fontsize=8, fontweight='bold')
    
    # 添加分组图例
    ax.text(0, 230, "Disease", color='#E74C3C', ha='center', fontweight='bold')
    ax.text(np.pi, 230, "Normal", color='#3498DB', ha='center', fontweight='bold')
    
    # 添加密度图例
    ax.text(0, 250, "Cell Density Distribution Map", color='black', ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('Tem_radial_plot_group.png', dpi=300, bbox_inches='tight')
    plt.show()

# 主程序
if __name__ == "__main__":
    # 1. 加载并过滤数据
    adata = load_and_filter_data()
    
    # 定义目标细胞类型
    target_celltypes = ["CD4+ Tem", "CD4+ Treg", "CD8+ Tex", "M2", "iDC",
                        "THY1+ FB", "Kera Spinous", "Kera Inflammatory",
                        "\u200bCTMCs", "Disease Stem", "POSTN+ FB",'CD4+ CD247+ T',"SBC Proliferation","SBC Spinous","SIGLEC1+ Mø","LAMP3+ DC","cDC1","cDC2"]
    #计算细胞数量统计和密度
    mean_counts, densities, distances, sector_areas = calculate_cell_counts(adata)
    
    # 导出计数结果到CSV
    export_counts_to_csv(mean_counts, densities, distances, sector_areas, target_celltypes)
    
    # 3-6. 创建环形图
    create_radial_plot(mean_counts, densities, distances, target_celltypes)
