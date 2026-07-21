import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import matplotlib as mpl
import math

# 设置全局字体
plt.rcParams['font.size'] = 18

# 1. 加载数据并提取指定细胞类型
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter_state_samp_new_with_CNs10-15.h5ad')
    target_celltypes = ["RSPO3+ GREM1+ Crypt Fib","Prol Metab Epi","EMT-like Inv Ent","CXCL12+ Fib","CD8+ T","Pericyte-like Fib","Early stress Epi","WNT+ ISC"]
    
    # 筛选目标细胞类型
    mask = adata.obs['celltype'].isin(target_celltypes)
    adata = adata[mask].copy()
    
    # 验证实际存在的细胞类型
    existing_celltypes = adata.obs['celltype'].unique()
    missing_types = [ct for ct in target_celltypes if ct not in existing_celltypes]
    if missing_types:
        print(f"警告: 以下细胞类型在数据中不存在: {', '.join(missing_types)}")
    
    return adata, target_celltypes

# 2. 计算不同距离范围内的细胞数量和密度（state顺序固定为N,A,C）
def calculate_cell_counts(adata, target_celltypes):
    # 获取中心细胞坐标（WNT+ ISC）
    center_mask = adata.obs['celltype'] == "WNT+ ISC"
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    # 如果没有中心细胞，使用空间中心点作为参考
    if num_centers == 0:
        print("警告: 没有找到 WNT+ ISC 细胞，使用空间中心点作为参考")
        center_coords = np.mean(adata.obsm['spatial'], axis=0).reshape(1, -1)
        num_centers = 1
    
    # 获取其他细胞类型的坐标和元数据
    other_mask = ~center_mask
    
    # 强制 state 顺序为 N, A, C，并过滤掉其他状态
    all_groups = ['N', 'A', 'C']
    
    # 先获取其他细胞的 obs，再按 state 过滤
    other_obs = adata[other_mask].obs
    valid_state = other_obs['state'].isin(all_groups)
    # 使用 valid_state 再次过滤 other_mask
    other_mask_indices = other_mask[other_mask].index[valid_state]
    other_mask = adata.obs.index.isin(other_mask_indices)
    
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    other_state = adata[other_mask].obs['state']
    
    # 重新计算有效中心数量（中心细胞本身不需要 state 过滤）
    num_centers = len(center_coords)
    
    # 构建映射
    group_to_idx = {g: i for i, g in enumerate(all_groups)}
    other_group = other_state.map(group_to_idx).values
    
    # 构建 KDTree
    tree = KDTree(other_coords)
    
    # 定义距离阈值 (μm)
    distances = [0, 50, 100, 150, 200]
    
    # 非中心细胞类型（所有目标细胞类型中排除 WNT+ ISC）
    non_center_celltypes = [ct for ct in target_celltypes if ct != "WNT+ ISC"]
    num_celltypes = len(non_center_celltypes)
    
    # 初始化计数数组：距离环数 × 细胞类型数 × 组数
    counts = np.zeros((len(distances)-1, num_celltypes, len(all_groups)))
    
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
    densities = np.zeros_like(mean_counts)
    
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / num_celltypes  # 每个细胞类型占一个扇形
        for sector_idx in range(num_celltypes):
            sector_areas[dist_idx, sector_idx] = sector_area
            for grp_idx in range(len(all_groups)):
                densities[dist_idx, sector_idx, grp_idx] = mean_counts[dist_idx, sector_idx, grp_idx] / sector_area
    
    return mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups

def export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups, filename='WNT+_ISC_cell_counts.csv'):
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
        'celltype', 'distance_range', 'state', 'cell_count', 'sector_area', 'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"细胞计数已导出到 {filename}")

# 创建径向气泡图（标签外移、图例右移、画布加宽）
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes, all_celltypes, all_groups):
    num_groups = len(all_groups)
    
    # 增大画布宽度
    plt.figure(figsize=(16, 12))
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
    
    # 添加细胞类型标签（外圈——向外移动）
    for i, ct in enumerate(all_celltypes):
        if ct == "CCL19+ CAF":
            continue
        angle = (i + 0.5) * 2*np.pi / num_sectors
        y = 255                        # 原本为230，现在向外移动到255
        ax.text(angle, y, ct, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.7))
    
    # 每个扇形内，各组气泡的偏移角度
    total_offset_span = 0.28
    if num_groups > 1:
        offsets = np.linspace(-total_offset_span/2, total_offset_span/2, num_groups)
    else:
        offsets = [0]
    
    # 修改分组配色：N为绿色，A为蓝色，C为橙色
    group_colors = {
        'N': '#2ca02c',   # 绿色
        'A': '#1f77b4',   # 蓝色
        'C': '#ff7f0e'    # 橙色
    }
    colors = [group_colors[g] for g in all_groups]
    
    # 气泡大小缩放因子
    scale_factors = [5000, 8000, 10000, 12000]
    
    # 绘制气泡
    for ring_idx in range(len(distances)-1):
        r_center = distances[ring_idx+1] - 5
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            all_ct_idx = all_celltypes.index(cell_type)
            base_angle = (all_ct_idx + 0.5) * 2*np.pi / num_sectors
            density_vals = densities[ring_idx, ct_idx, :]
            for grp_idx in range(num_groups):
                angle = base_angle + offsets[grp_idx]
                size = density_vals[grp_idx] * scale_factors[min(ring_idx, len(scale_factors)-1)]
                size = np.clip(size, 5, 500)
                color = colors[grp_idx]
                ax.scatter(angle, r_center, s=size, facecolor=color, edgecolor='none', alpha=0.8, zorder=10)
    
    # 扩大 y 轴上限，确保标签不被裁剪
    ax.set_ylim(0, 270)                # 原来为260，现在增大到270
    
    # 添加距离范围标签（在圆环外侧）
    label_angles = np.linspace(0, 2*np.pi, len(distances))[:-1]
    for i, dist in enumerate(distances[1:]):
        label = f"{distances[i]}-{distances[i+1]}μm"
        ax.text(label_angles[i], distances[i+1]+15, label, ha='center', va='center', fontsize=8, fontweight='bold')
    
    # 图例右移，设置更大的锚点
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=grp_name)
                       for grp_name, color in zip(all_groups, colors)]
    
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.6, 0.5),
              title="Groups", framealpha=0.9)
    
    # 压缩主图区域，给图例留出足够空间
    plt.subplots_adjust(right=0.55)
    
    ax.text(0, 265, "Cell Density Distribution Map", color='black', ha='center', fontsize=12)  # 标题也相应上移
    
    plt.tight_layout()
    plt.savefig('WNT+_ISC_plot_group.png', dpi=300, bbox_inches='tight')
    plt.savefig('WNT+_ISC_plot_group.pdf', dpi=300, bbox_inches='tight')
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
