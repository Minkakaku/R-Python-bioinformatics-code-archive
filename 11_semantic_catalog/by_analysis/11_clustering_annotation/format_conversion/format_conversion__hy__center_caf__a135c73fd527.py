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
    adata = sc.read_h5ad('extracted_samples.h5ad')
    target_celltypes = ["CXCL9+ Tumor","M1","Tissue repairing M2","CXCL9+ CAF","CD8+ T","CD4+ T","Treg","CXCL9+ M1","Proliferative Tumor","EMT Tumor"]
    
    # 筛选目标细胞类型
    mask = adata.obs['celltype'].isin(target_celltypes)
    adata = adata[mask].copy()
    
    # 验证实际存在的细胞类型
    existing_celltypes = adata.obs['celltype'].unique()
    missing_types = [ct for ct in target_celltypes if ct not in existing_celltypes]
    if missing_types:
        print(f"警告: 以下细胞类型在数据中不存在: {', '.join(missing_types)}")
    
    return adata, target_celltypes

# 2. 计算不同距离范围内的细胞数量和密度
def calculate_cell_counts(adata, target_celltypes):
    # 获取CD4+ Tem细胞的坐标
    center_mask = adata.obs['celltype'] == "CXCL9+ CAF"
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    # 如果没有中心细胞，使用默认位置
    if num_centers == 0:
        print("警告: 没有找到CD4+ Tem细胞，使用空间中心点作为参考")
        center_coords = np.mean(adata.obsm['spatial'], axis=0).reshape(1, -1)
        num_centers = 1
    
    # 获取其他细胞类型的坐标和元数据
    other_mask = ~center_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    other_group = adata[other_mask].obs['group'].map({'K_pCR_P': 0, 'K_pCR_A': 1})
    
    # 构建KDTree
    tree = KDTree(other_coords)
    
    # 定义距离阈值
    distances = [0, 50, 100, 150, 200]
    
    # 修复点1：使用完整的细胞类型数量（包括中心类型）
    non_center_celltypes = [ct for ct in target_celltypes if ct != "CXCL9+ CAF"]
    num_celltypes = len(non_center_celltypes)
    
    # 初始化计数数组 - 修复点2：使用正确的维度
    counts = np.zeros((len(distances)-1, num_celltypes, 2))
    
    # 对每个中心细胞计算距离范围内的细胞数量
    for center in center_coords:
        # 使用query_ball_point获取200μm内所有点的索引
        indices = tree.query_ball_point(center, r=200)
        
        for idx in indices:
            pt = other_coords[idx]
            dist = np.linalg.norm(center - pt)  # 计算欧氏距离
            
            # 获取细胞类型和分组信息
            ct = other_celltypes.iloc[idx]
            grp = other_group.iloc[idx]
            
            # 修复点3：更健壮的细胞类型匹配
            if ct not in non_center_celltypes:
                continue
                
            ct_idx = non_center_celltypes.index(ct)
            
            # 确定距离区间
            for i in range(len(distances)-1):
                if distances[i] <= dist < distances[i+1]:
                    counts[i, ct_idx, grp] += 1
                    break
    
    # 计算平均数量
    mean_counts = counts / num_centers if num_centers > 0 else counts
    
    # 计算扇形面积和细胞密度
    sector_areas = np.zeros((len(distances)-1, num_celltypes))
    densities = np.zeros_like(mean_counts)
    
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        # 计算圆环面积
        ring_area = math.pi * (r2**2 - r1**2)
        # 计算单个扇形面积
        sector_area = ring_area / num_celltypes
        
        for sector_idx in range(num_celltypes):
            sector_areas[dist_idx, sector_idx] = sector_area
            # 计算密度（细胞数/面积）
            densities[dist_idx, sector_idx, 0] = mean_counts[dist_idx, sector_idx, 0] / sector_area
            densities[dist_idx, sector_idx, 1] = mean_counts[dist_idx, sector_idx, 1] / sector_area
    
    return mean_counts, densities, distances, sector_areas, non_center_celltypes

# 导出计数结果到CSV文件
def export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, filename='cell_counts.csv'):
    # 创建距离范围标签
    dist_labels = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    
    # 准备数据
    data = []
    for ct_idx, cell_type in enumerate(non_center_celltypes):
        for dist_idx in range(len(dist_labels)):
            area = sector_areas[dist_idx, ct_idx]
            # 疾病组数据
            data.append([
                cell_type,
                dist_labels[dist_idx],
                'K_pCR_P',
                mean_counts[dist_idx, ct_idx, 0],
                area,
                densities[dist_idx, ct_idx, 0]
            ])
            # 正常组数据
            data.append([
                cell_type,
                dist_labels[dist_idx],
                'K_pCR_A',
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
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes, all_celltypes):
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
    
    # 4. 绘制分隔线 (扇形分隔线数量等于所有细胞类型数量)
    for i in range(len(all_celltypes)):
        angle = i * 2*np.pi/len(all_celltypes)
        ax.plot([angle, angle], [0, 230], 
                color='#34495E', linestyle='-', linewidth=1, alpha=0.6)
    
    # 5. 添加细胞类型标签 (外圈水平放置)
    for i, ct in enumerate(all_celltypes):
        if ct == "CXCL9+ CAF":  # 跳过中心细胞类型标签
            continue
            
        angle = (i + 0.5) * 2*np.pi/len(all_celltypes)
        x = angle
        y = 230  # 固定在200μm圆环外侧
        
        # 添加白色背景
        ax.text(x, y, ct, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', 
                         facecolor='white', 
                         edgecolor='none', 
                         alpha=0.7))
    
    # 6. 绘制气泡图（使用group分组）
    offset_angles = [-0.08, 0.08]  # 疾病和正常气泡的偏移角度
    
    # 细胞类型和分组颜色映射
    cmap_dis = mpl.colormaps['Reds']
    cmap_norm = mpl.colormaps['Blues']
    
    # 遍历每个距离环、每种细胞类型（仅非中心细胞）
    for ring_idx in range(len(distances)-1):
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            # 找到在所有细胞类型中的位置
            all_ct_idx = all_celltypes.index(cell_type)
            
            density_vals = densities[ring_idx, ct_idx]
            base_angle = (all_ct_idx + 0.5) * 2*np.pi/len(all_celltypes)
            radius = distances[ring_idx + 1] - 5
            
            # 分别绘制疾病和正常气泡
            for grp_idx, offset in enumerate(offset_angles):
                angle = base_angle + offset
                
                # 根据距离层级使用不同的缩放因子
                if ring_idx == 0:  # 0-50um层级
                    size = density_vals[grp_idx] * 20000
                elif ring_idx == 1:  # 50-100um层级
                    size = density_vals[grp_idx] * 40000
                else:  # 其他层级
                    size = density_vals[grp_idx] * 100000
                
                # 选择颜色映射
                cmap = cmap_dis if grp_idx == 0 else cmap_norm
                color = cmap(0.6 + grp_idx*0.2)
                
                # 绘制实心圆点（无边框）
                ax.scatter(angle, radius, s=size, 
                           facecolor=color, 
                           edgecolor='none',
                           alpha=0.8, 
                           zorder=10,
                           marker='o')
    
    # 设置显示范围
    ax.set_ylim(0, 260)
    
    # 添加距离标签
    label_positions = np.linspace(0, 2*np.pi, 7)[:-1]
    for i, dist in enumerate(distances[1:]):
        label = f"{distances[i]}-{distances[i+1]}μm"
        ax.text(label_positions[i], dist+15, label, 
                ha='center', va='center', 
                fontsize=8, fontweight='bold')
    
    # 添加分组图例
    ax.text(0, 240, "K_pCR_P", color='#E74C3C', ha='center', fontweight='bold')
    ax.text(np.pi, 240, "K_pCR_A", color='#3498DB', ha='center', fontweight='bold')
    
    # 添加密度图例
    ax.text(0, 250, "Cell Density Distribution Map", color='black', ha='center', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('CXCL9_CAF_plot_group.png', dpi=300, bbox_inches='tight')
    plt.savefig('CXCL9_CAF_plot_group.pdf', dpi=300, bbox_inches='tight')
    plt.show()

# 主程序
if __name__ == "__main__":
    # 1. 加载并过滤数据
    adata, all_celltypes = load_and_filter_data()
    
    # 2. 计算细胞数量统计和密度
    mean_counts, densities, distances, sector_areas, non_center_celltypes = calculate_cell_counts(adata, all_celltypes)
    
    # 导出计数结果到CSV
    export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes)
    
    # 3. 创建环形图
    create_radial_plot(mean_counts, densities, distances, non_center_celltypes, all_celltypes)
