import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import matplotlib as mpl
import math

# 设置全局字体
plt.rcParams['font.size'] = 10

# 1. 加载数据并提取指定细胞类型
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter.h5ad')
    target_celltypes = ['Wnt-active stem','CD8+ T','CXCL9+ CD8+ T','MKI67+ T','M2','CXCL9+ M1','M1','EndoMT','Activated Endo','Microvascular Endo','CD4+ T']
    
    # 筛选目标细胞类型
    mask = adata.obs['celltype'].isin(target_celltypes)
    adata = adata[mask].copy()
    
    # 验证实际存在的细胞类型
    existing_celltypes = adata.obs['celltype'].unique()
    missing_types = [ct for ct in target_celltypes if ct not in existing_celltypes]
    if missing_types:
        print(f"警告: 以下细胞类型在数据中不存在: {', '.join(missing_types)}")
    
    return adata, target_celltypes

# 2. 计算不同距离范围内的细胞数量和密度（支持多组）
def calculate_cell_counts(adata, target_celltypes):
    # 获取CCL19+ CAF细胞的坐标（作为中心）
    center_mask = adata.obs['celltype'] == "Wnt-active stem"
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    # 如果没有中心细胞，使用空间中心点作为参考
    if num_centers == 0:
        print("警告: 没有找到CCL19+ CAF细胞，使用空间中心点作为参考")
        center_coords = np.mean(adata.obsm['spatial'], axis=0).reshape(1, -1)
        num_centers = 1
    
    # 获取其他细胞类型的坐标和元数据
    other_mask = ~center_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    
    # --- 多组处理：获取所有唯一分组名称 ---
    all_groups = adata.obs['group_final'].unique()
    all_groups = sorted([g for g in all_groups if pd.notna(g)])  # 去除NaN并排序
    num_groups = len(all_groups)
    group_to_idx = {g: i for i, g in enumerate(all_groups)}
    other_group = adata[other_mask].obs['group_final'].map(group_to_idx).values
    
    # 构建KDTree
    tree = KDTree(other_coords)
    
    # 定义距离阈值 (μm)
    distances = [0, 50, 100, 150, 200]
    
    # 非中心细胞类型（所有目标细胞类型中排除CCL19+ CAF）
    non_center_celltypes = [ct for ct in target_celltypes if ct != "Wnt-active stem"]
    num_celltypes = len(non_center_celltypes)
    
    # 初始化计数数组：距离环数 × 细胞类型数 × 组数
    counts = np.zeros((len(distances)-1, num_celltypes, num_groups))
    
    # 对每个中心细胞计算距离范围内的细胞数量
    for center in center_coords:
        indices = tree.query_ball_point(center, r=200)
        for idx in indices:
            dist = np.linalg.norm(center - other_coords[idx])
            ct = other_celltypes.iloc[idx]
            grp_idx = other_group[idx]
            
            if ct not in non_center_celltypes:
                continue
            ct_idx = non_center_celltypes.index(ct)
            
            # 确定距离区间
            for i in range(len(distances)-1):
                if distances[i] <= dist < distances[i+1]:
                    counts[i, ct_idx, grp_idx] += 1
                    break
    
    # 计算平均数量（每个中心细胞的平均）
    mean_counts = counts / num_centers if num_centers > 0 else counts
    
    # 计算扇形面积和细胞密度
    sector_areas = np.zeros((len(distances)-1, num_celltypes))
    densities = np.zeros_like(mean_counts)  # 形状同 mean_counts
    
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / num_celltypes  # 每个细胞类型占一个扇形
        for sector_idx in range(num_celltypes):
            sector_areas[dist_idx, sector_idx] = sector_area
            for grp_idx in range(num_groups):
                densities[dist_idx, sector_idx, grp_idx] = mean_counts[dist_idx, sector_idx, grp_idx] / sector_area
    
    return mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups

def export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups, filename='cell_counts.csv'):
    dist_labels = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    
    data = []
    for ct_idx, cell_type in enumerate(non_center_celltypes):
        for dist_idx in range(len(dist_labels)):
            area = sector_areas[dist_idx, ct_idx]
            for grp_idx, grp_name in enumerate(all_groups):
                data.append([
                    cell_type,
                    dist_labels[dist_idx],
                    grp_name,
                    mean_counts[dist_idx, ct_idx, grp_idx],
                    area,
                    densities[dist_idx, ct_idx, grp_idx]
                ])
    
    df = pd.DataFrame(data, columns=[
        'celltype', 'distance_range', 'group_final', 'cell_count', 'sector_area', 'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"细胞计数已导出到 {filename}")

# 创建径向气泡图（支持多组）
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes, all_celltypes, all_groups):
    num_groups = len(all_groups)
    
    plt.figure(figsize=(12, 12))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    
    # 删除默认网格和标签
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    
    # 绘制同心圆（距离边界）
    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 100), [r]*100, color='gray', linestyle='-', linewidth=1, alpha=0.5)
    
    # 绘制分隔线（数量等于所有细胞类型数量）
    num_sectors = len(all_celltypes)
    for i in range(num_sectors):
        angle = i * 2*np.pi / num_sectors
        ax.plot([angle, angle], [0, 230], color='#34495E', linestyle='-', linewidth=1, alpha=0.6)
    
    # 添加细胞类型标签（外圈）
    for i, ct in enumerate(all_celltypes):
        if ct == "CCL19+ CAF":
            continue
        angle = (i + 0.5) * 2*np.pi / num_sectors
        y = 230
        ax.text(angle, y, ct, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.7))
    
    # --- 多组绘图参数 ---
    # 每个扇形内，各组气泡的偏移角度（均匀分布在扇形中心两侧）
    # 总偏移宽度（弧度）根据组数动态调整，避免重叠
    total_offset_span = 0.12  # 总跨度（弧度）
    if num_groups > 1:
        offsets = np.linspace(-total_offset_span/2, total_offset_span/2, num_groups)
    else:
        offsets = [0]
    
    # 颜色映射：使用 tab10 或自定义列表，确保足够颜色
    colormap = plt.cm.tab10  # 最多10种不同颜色
    colors = [colormap(i % 10) for i in range(num_groups)]
    
    # 气泡大小缩放因子（不同距离环需要不同倍数，使气泡大小视觉合适）
    scale_factors = [20000, 40000, 60000, 80000]  # 对应0-50,50-100,100-150,150-200
    # 如果环数超过4，则扩展（这里只有4环）
    
    # 遍历每个距离环、每种非中心细胞类型
    for ring_idx in range(len(distances)-1):
        r_center = distances[ring_idx+1] - 5  # 气泡绘制在圆环中部偏外侧
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            all_ct_idx = all_celltypes.index(cell_type)
            base_angle = (all_ct_idx + 0.5) * 2*np.pi / num_sectors
            
            # 获取该细胞类型在该环上的各组密度值
            density_vals = densities[ring_idx, ct_idx, :]  # 长度为 num_groups
            
            for grp_idx in range(num_groups):
                angle = base_angle + offsets[grp_idx]
                size = density_vals[grp_idx] * scale_factors[min(ring_idx, len(scale_factors)-1)]
                # 限制最小和最大气泡尺寸，避免过大或过小
                size = np.clip(size, 5, 500)
                color = colors[grp_idx]
                ax.scatter(angle, r_center, s=size, facecolor=color, edgecolor='none', alpha=0.8, zorder=10)
    
    # 设置显示范围
    ax.set_ylim(0, 260)
    
    # 添加距离范围标签（在圆环外侧）
    label_angles = np.linspace(0, 2*np.pi, len(distances))[:-1]
    for i, dist in enumerate(distances[1:]):
        label = f"{distances[i]}-{distances[i+1]}μm"
        ax.text(label_angles[i], distances[i+1]+12, label, ha='center', va='center', fontsize=8, fontweight='bold')
    
    # 添加分组图例（放在图的上方两侧）
    # 计算图例位置：左右两侧
    legend_x_positions = np.linspace(0, 2*np.pi, num_groups+2)[1:-1]  # 均匀分布在圆周上
    # 更简单：在底部左右两侧或顶部集中显示
    # 这里选择在顶部（角度0附近）和底部（角度pi附近）分布，如果组数太多则换行
    if num_groups <= 4:
        for i, (grp_name, color) in enumerate(zip(all_groups, colors)):
            angle = (i - (num_groups-1)/2) * 0.3  # 弧度偏移
            ax.text(angle, 245, grp_name, color=color, ha='center', fontweight='bold', fontsize=10)
    else:
        # 组数较多时，在图右侧加一个矩形图例
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=grp_name)
                           for grp_name, color in zip(all_groups, colors)]
        ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 1.0), title="Groups")
    
    ax.text(0, 255, "Cell Density Distribution Map", color='black', ha='center', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('Wnt_immu_plot_group.png', dpi=300, bbox_inches='tight')
    plt.savefig('Wnt_immu_plot_group.pdf', dpi=300, bbox_inches='tight')
    plt.show()

# 主程序
if __name__ == "__main__":
    # 1. 加载并过滤数据
    adata, all_celltypes = load_and_filter_data()
    
    # 2. 计算细胞数量统计和密度
    mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups = calculate_cell_counts(adata, all_celltypes)
    
    # 导出计数结果到CSV
    export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups)
    
    # 3. 创建环形图
    create_radial_plot(mean_counts, densities, distances, non_center_celltypes, all_celltypes, all_groups)
