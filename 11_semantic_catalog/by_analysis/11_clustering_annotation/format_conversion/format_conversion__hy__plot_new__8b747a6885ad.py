import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import math

# ================== 全局设置 ==================
plt.rcParams['font.family'] = 'arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ================== 1. 数据加载 ==================
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter_state_samp_new_with_CNs10-15.h5ad')
    target_celltypes = [
        "WNT+ ISC", "matri-prod CAF", "CD8+ T", "Treg", "MKI67+ T", "CTSL+ PLAUR+ Macro"
    ]
    mask = adata.obs['celltype'].isin(target_celltypes)
    adata = adata[mask].copy()
    existing_celltypes = adata.obs['celltype'].unique()
    missing_types = [ct for ct in target_celltypes if ct not in existing_celltypes]
    if missing_types:
        print(f"警告: 以下细胞类型不存在: {', '.join(missing_types)}")
    return adata, target_celltypes

# ================== 2. 细胞计数与密度计算（环边界 50-250）==================
def calculate_cell_counts(adata, target_celltypes, center_celltype):
    center_mask = adata.obs['celltype'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    if num_centers == 0:
        print(f"警告: 未找到 {center_celltype}，使用全部细胞重心")
        center_coords = np.mean(adata.obsm['spatial'], axis=0).reshape(1, -1)
        num_centers = 1
    
    other_mask = ~center_mask
    all_groups = ['N', 'A', 'C']
    
    other_obs = adata[other_mask].obs
    valid_state = other_obs['state'].isin(all_groups)
    other_mask_indices = other_mask[other_mask].index[valid_state]
    other_mask = adata.obs.index.isin(other_mask_indices)
    
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    other_state = adata[other_mask].obs['state']
    
    group_to_idx = {g: i for i, g in enumerate(all_groups)}
    other_group = other_state.map(group_to_idx).values
    
    tree = KDTree(other_coords)
    distances = [50, 100, 150, 200, 250]
    
    non_center_celltypes = [ct for ct in target_celltypes if ct != center_celltype]
    num_celltypes = len(non_center_celltypes)
    
    counts = np.zeros((len(distances)-1, num_celltypes, len(all_groups)))
    
    for center in center_coords:
        indices = tree.query_ball_point(center, r=250)
        for idx in indices:
            dist = np.linalg.norm(center - other_coords[idx])
            ct = other_celltypes.iloc[idx]
            grp_idx = other_group[idx]
            if ct not in non_center_celltypes:
                continue
            ct_idx = non_center_celltypes.index(ct)
            for i in range(len(distances)-1):
                if distances[i] <= dist < distances[i+1]:
                    counts[i, ct_idx, grp_idx] += 1
                    break
    
    mean_counts = counts / num_centers if num_centers > 0 else counts
    
    sector_areas = np.zeros((len(distances)-1, num_celltypes))
    densities = np.zeros_like(mean_counts)
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / len(target_celltypes)
        for ct_idx in range(num_celltypes):
            sector_areas[dist_idx, ct_idx] = sector_area
            for grp_idx in range(len(all_groups)):
                densities[dist_idx, ct_idx, grp_idx] = mean_counts[dist_idx, ct_idx, grp_idx] / sector_area
    
    return mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups

# ================== 3. 导出CSV ==================
def export_counts_to_csv(mean_counts, densities, distances, sector_areas, 
                         non_center_celltypes, all_groups, center_celltype):
    safe_name = center_celltype.replace(" ", "_").replace("/", "_").replace("+", "plus")
    filename = f'{safe_name}_cell_counts.csv'
    dist_labels = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    data = []
    for ct_idx, cell_type in enumerate(non_center_celltypes):
        for dist_idx in range(len(dist_labels)):
            area = sector_areas[dist_idx, ct_idx]
            for grp_idx, grp_name in enumerate(all_groups):
                data.append([
                    cell_type, dist_labels[dist_idx], grp_name,
                    mean_counts[dist_idx, ct_idx, grp_idx],
                    area, densities[dist_idx, ct_idx, grp_idx]
                ])
    df = pd.DataFrame(data, columns=[
        'celltype', 'distance_range', 'state', 'cell_count', 'sector_area', 'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"数据已导出: {filename}")

# ================== 4. 绘图函数 ==================
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                       all_celltypes, all_groups, center_celltype):
    """
    每个扇形内三个状态的气泡位于同一半径（环的中径），
    角度偏移量增大，几乎占满整个扇形。
    中心细胞文字严格置于圆心。
    外圈标签再向外移动，距离标签字体缩小。
    """
    num_groups = len(all_groups)
    num_sectors = len(all_celltypes)
    sector_angle = 2 * np.pi / num_sectors   # 每个扇形的角度宽度
    
    # ----- 固定配色 -----
    group_colors = {
        'N': '#4caf50',   # 绿色
        'A': '#2196f3',   # 蓝色
        'C': '#ff9800'    # 橙色
    }
    colors = [group_colors[g] for g in all_groups]
    
    # 每个状态的角度偏移量（相对于扇区中心）
    angle_offset_factor = 0.33   # 可以调节，0.35 表示占据 70% 的扇形宽度
    angle_offset = {
        'N': -angle_offset_factor * sector_angle,
        'A': 0,
        'C': angle_offset_factor * sector_angle
    }
    
    # ----- 画布设置 -----
    fig = plt.figure(figsize=(7, 7), dpi=300, facecolor='white')
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    
    # ----- 绘制中心圆形区域（半径50）-----
    theta_center = np.linspace(0, 2*np.pi, 100)
    r_center = 50
    ax.fill_between(theta_center, 0, r_center, color='#EEEEEE', alpha=0.6, zorder=2)
    ax.plot(theta_center, [r_center]*100, color='#B0B0B0', linewidth=0.8, linestyle='-', alpha=0.8, zorder=3)
    
    # 中心细胞类型标签：严格放在圆心 (0,0)
    ct_lines = center_celltype.split(' ')
    if len(ct_lines) > 1:
        ct_text = '\n'.join(ct_lines)
    else:
        ct_text = center_celltype
    ax.text(0, 0, ct_text, ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9), zorder=4)
    
    # ----- 同心圆环（边界线：50,100,150,200,250）-----
    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 360), [r]*360,
                color='#D3D3D3', linewidth=0.5, linestyle='-', alpha=0.6, zorder=1)
    
    # ----- 扇形分隔虚线（从半径50到360）-----
    for i in range(num_sectors):
        angle = i * 2 * np.pi / num_sectors
        ax.plot([angle, angle], [50, 360], color='#B0B0B0',
                linewidth=0.4, linestyle='--', alpha=0.5, zorder=1)
    
    # ----- 扇形外圈背景（半径360-390，为标签提供浅色背景）-----
    sector_bg_colors = plt.cm.Greys(np.linspace(0.1, 0.25, num_sectors))
    for i in range(num_sectors):
        theta_start = i * 2 * np.pi / num_sectors
        theta_end = (i+1) * 2 * np.pi / num_sectors
        ax.fill_between(np.linspace(theta_start, theta_end, 50),
                        360, 390, color=sector_bg_colors[i], alpha=0.15, zorder=0)
    
    # ----- 外圈细胞类型标签（半径370，远离气泡区域）-----
    for i, ct in enumerate(all_celltypes):
        if ct == center_celltype:
            continue
        angle = (i + 0.5) * 2 * np.pi / num_sectors
        ax.text(angle, 360, ct, ha='center', va='center',
                fontsize=12, weight='bold', color='#333333',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          edgecolor='none', alpha=0.7))
    
    # ----- 距离数字（放在中心细胞扇区内，字体放大）-----
    center_idx = all_celltypes.index(center_celltype)
    center_angle = (center_idx + 0.5) * 2 * np.pi / num_sectors
    for ring_idx in range(len(distances)-1):
        r_mid = (distances[ring_idx] + distances[ring_idx+1]) / 2
        label = f"{distances[ring_idx]}–{distances[ring_idx+1]}"   # 只显示数字，不带单位
        ax.text(center_angle, r_mid, label, ha='center', va='center',
                fontsize=9, weight='bold', color='#555555',         # 字体放大到9
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.8))

    # ----- 在中心细胞扇区的外圈添加单位 “μm” -----
    ax.text(center_angle, 370, '(μm)', ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.7)) 
    # ----- 气泡绘制（同一半径，大角度偏移）-----
    max_density = np.max(densities) if np.max(densities) > 0 else 1
    base_scale = 400 / max_density
    min_size, max_size = 15, 220
    
    for ring_idx in range(len(distances)-1):
        r_inner = distances[ring_idx]
        r_outer = distances[ring_idx+1]
        r_center = (r_inner + r_outer) / 2   # 环的中径
        
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            all_ct_idx = all_celltypes.index(cell_type)
            base_angle = (all_ct_idx + 0.5) * 2 * np.pi / num_sectors
            
            for grp_idx, grp in enumerate(all_groups):
                angle = base_angle + angle_offset[grp]
                density_val = densities[ring_idx, ct_idx, grp_idx]
                size = density_val * base_scale
                size = np.clip(size, min_size, max_size)
                if size > 0:
                    ax.scatter(angle, r_center, s=size,
                               facecolor=colors[grp_idx], edgecolor='white',
                               linewidth=0.5, alpha=0.85, zorder=10)
    
    # ----- 图例（靠近图形）-----
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=color, markersize=8, label=grp)
                       for grp, color in zip(all_groups, colors)]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.02, 1.0), fontsize=9, title='State',
              title_fontsize=10, framealpha=0.9, edgecolor='#CCCCCC')
    
    # ----- 范围调整（ylim 要足够大容纳外圈标签）-----
    ax.set_ylim(0, 400)   # 增加到400，确保标签可见
    plt.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)
    
    # ----- 保存 -----
    safe_name = center_celltype.replace(" ", "_").replace("/", "_").replace("+", "plus")
    plt.savefig(f'{safe_name}_publication.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(f'{safe_name}_publication.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()

# ================== 5. 主程序 ==================
if __name__ == "__main__":
    adata, all_celltypes = load_and_filter_data()
    
    for center in all_celltypes:
        print(f"\n===== 处理中心细胞: {center} =====")
        mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups = \
            calculate_cell_counts(adata, all_celltypes, center)
        export_counts_to_csv(mean_counts, densities, distances, sector_areas,
                             non_center_celltypes, all_groups, center)
        create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                           all_celltypes, all_groups, center)
    
    print("\n所有图表生成完毕！")
