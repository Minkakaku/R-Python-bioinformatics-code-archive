import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import math
from tqdm import tqdm
from matplotlib import font_manager

# ================== 全局字体设置 ==================
font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ================== 1. 数据加载（免疫细胞版本）==================
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter_state_samp_new_ctzq_with_CNs10-15.h5ad')
    
    center_celltype = "Epithelial"
    target_celltypes = [
        "CTSL+ PLAUR+ Macro",
        "STAB1+ Macro",
        "CD4+ T",
        "γδ T",
        "CD8+ T",
        "Treg",
        "MKI67+ T",
        "DC",
        "GCB",
        "MBC",
        "Mast",
        "Plasmas",
        "Neutrophils"
    ]
    
    if 'celltype_zq' not in adata.obs.columns:
        raise KeyError("adata.obs 中缺少 'celltype_zq' 列")
    if 'celltype' not in adata.obs.columns:
        raise KeyError("adata.obs 中缺少 'celltype' 列")
    
    center_mask = adata.obs['celltype_zq'] == center_celltype
    if center_mask.sum() == 0:
        raise ValueError("没有找到任何 Epithelial 细胞")
    
    target_mask = adata.obs['celltype'].isin(target_celltypes) & (~center_mask)
    final_mask = center_mask | target_mask
    adata = adata[final_mask].copy()
    
    existing_celltypes = adata.obs['celltype'].unique()
    missing_types = [ct for ct in target_celltypes if ct not in existing_celltypes]
    if missing_types:
        print(f"警告: 以下目标细胞类型不存在: {', '.join(missing_types)}")
    
    print(f"中心细胞 (Epithelial): {center_mask.sum()} 个")
    print(f"周边目标细胞: {target_mask.sum()} 个")
    
    return adata, center_celltype, target_celltypes

# ================== 2. 高性能细胞计数与密度计算（环边界50-250）==================
def calculate_cell_counts(adata, center_celltype, target_celltypes):
    center_mask = adata.obs['celltype_zq'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    print(f"使用全部 {num_centers} 个中心细胞进行计算")
    
    other_mask = ~center_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    
    # 处理分组信息（state 列）
    if 'state' in adata.obs.columns:
        other_state = adata[other_mask].obs['state']
        all_groups = sorted(set(other_state.dropna().unique()))
        all_groups = [str(g) for g in all_groups]
        group_to_idx = {g: i for i, g in enumerate(all_groups)}
        other_group = other_state.astype(str).map(group_to_idx).values
    else:
        all_groups = ['All']
        group_to_idx = {'All': 0}
        other_group = np.zeros(len(other_coords), dtype=int)
    
    tree = KDTree(other_coords)
    # 新环边界：从50 μm到250 μm
    distances = [50, 100, 150, 200, 250]
    n_rings = len(distances) - 1
    ring_bins = distances[1:]   # [100,150,200,250]
    
    # 只保留实际存在的目标细胞类型
    unique_other_celltypes = other_celltypes.unique()
    non_center_celltypes = [ct for ct in target_celltypes if ct in unique_other_celltypes]
    num_celltypes = len(non_center_celltypes)
    num_groups = len(all_groups)
    
    ct_to_idx = {ct: i for i, ct in enumerate(non_center_celltypes)}
    counts = np.zeros((n_rings, num_celltypes, num_groups))
    
    # 遍历每个中心细胞（带进度条）
    for center in tqdm(center_coords, desc="处理中心细胞"):
        indices = tree.query_ball_point(center, r=distances[-1])
        if len(indices) == 0:
            continue
        
        dists = np.linalg.norm(center - other_coords[indices], axis=1)
        cts = other_celltypes.iloc[indices].values
        grps = other_group[indices]
        
        mask = np.isin(cts, non_center_celltypes)
        if not np.any(mask):
            continue
        dists = dists[mask]
        cts = cts[mask]
        grps = grps[mask]
        
        ring_ids = np.digitize(dists, ring_bins, right=False)   # 0~n_rings-1
        ct_ids = np.array([ct_to_idx[ct] for ct in cts])
        
        # 累加计数
        for i in range(len(dists)):
            counts[ring_ids[i], ct_ids[i], grps[i]] += 1
    
    mean_counts = counts / num_centers if num_centers > 0 else counts
    
    # 计算扇形面积和密度
    sector_areas = np.zeros((n_rings, num_celltypes))
    densities = np.zeros_like(mean_counts)
    for dist_idx in range(n_rings):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / num_celltypes
        for ct_idx in range(num_celltypes):
            sector_areas[dist_idx, ct_idx] = sector_area
            for grp_idx in range(num_groups):
                densities[dist_idx, ct_idx, grp_idx] = mean_counts[dist_idx, ct_idx, grp_idx] / sector_area
    
    return mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups

# ================== 3. 导出CSV ==================
def export_counts_to_csv(mean_counts, densities, distances, sector_areas,
                         non_center_celltypes, all_groups, filename='Epithelial_immune_cell_counts.csv'):
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
    print(f"细胞计数已导出到 {filename}")

# ================== 4. 顶级期刊风格径向气泡图 ==================
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                       all_celltypes, all_groups):
    """
    绘制符合顶级期刊要求的径向气泡图：
    - 中心细胞表示为半径50的圆形区域
    - 环从50到250 μm，每个环内气泡位于环中径
    - 不同分组通过角度偏移均匀分布
    - 外圈细胞标签放在远距离，距离数字字体放大且单位单独标注
    """
    num_groups = len(all_groups)
    num_sectors = len(all_celltypes)
    sector_angle = 2 * np.pi / num_sectors
    
    # 配色（若分组为 N/A/C 使用原始配色，否则使用 tab10）
    predefined_colors = {'N': '#4caf50', 'A': '#2196f3', 'C': '#ff9800'}
    if set(all_groups).issubset(predefined_colors.keys()):
        colors = [predefined_colors[g] for g in all_groups]
    else:
        cmap = plt.cm.tab10
        colors = [cmap(i % 10) for i in range(num_groups)]
    
    # 角度偏移因子（覆盖扇形宽度的 70%）
    angle_offset_factor = 0.35
    offsets = np.linspace(-angle_offset_factor * sector_angle,
                          angle_offset_factor * sector_angle, num_groups)
    
    # 画布
    fig = plt.figure(figsize=(7, 7), dpi=300, facecolor='white')
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    
    # 中心圆形区域（半径50）
    theta_center = np.linspace(0, 2*np.pi, 100)
    r_center = 50
    ax.fill_between(theta_center, 0, r_center, color='#EEEEEE', alpha=0.6, zorder=2)
    ax.plot(theta_center, [r_center]*100, color='#B0B0B0', linewidth=0.8, linestyle='-', alpha=0.8, zorder=3)
    
    # 中心细胞标签（圆心）
    ax.text(0, 0, "Epithelial", ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9), zorder=4)
    
    # 同心圆环（边界线：50,100,150,200,250）
    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 360), [r]*360,
                color='#D3D3D3', linewidth=0.5, linestyle='-', alpha=0.6, zorder=1)
    
    # 扇形分隔虚线（从半径50到370）
    for i in range(num_sectors):
        angle = i * 2 * np.pi / num_sectors
        ax.plot([angle, angle], [50, 370], color='#B0B0B0',
                linewidth=0.4, linestyle='--', alpha=0.5, zorder=1)
    
    # 扇形外圈背景（半径370-390）
    sector_bg_colors = plt.cm.Greys(np.linspace(0.1, 0.25, num_sectors))
    for i in range(num_sectors):
        theta_start = i * 2 * np.pi / num_sectors
        theta_end = (i+1) * 2 * np.pi / num_sectors
        ax.fill_between(np.linspace(theta_start, theta_end, 50),
                        370, 390, color=sector_bg_colors[i], alpha=0.15, zorder=0)
    
    # 外圈细胞类型标签（半径380）
    for i, ct in enumerate(all_celltypes):
        angle = (i + 0.5) * 2 * np.pi / num_sectors
        ax.text(angle, 380, ct, ha='center', va='center',
                fontsize=9, weight='normal', color='#333333',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          edgecolor='none', alpha=0.7))
    
    # 距离数字（放在角度0处，字体放大）
    label_angle = 0.0
    for ring_idx in range(len(distances)-1):
        r_mid = (distances[ring_idx] + distances[ring_idx+1]) / 2
        label = f"{distances[ring_idx]}–{distances[ring_idx+1]}"
        ax.text(label_angle, r_mid, label, ha='center', va='center',
                fontsize=9, weight='bold', color='#555555',
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.8))
    
    # 单位 μm（放在外圈）
    ax.text(label_angle, 390, 'μm', ha='center', va='center',
            fontsize=9, weight='normal', color='#333333',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.7))
    
    # 气泡绘制（同一半径，角度偏移）
    max_density = np.max(densities) if np.max(densities) > 0 else 1
    base_scale = 400 / max_density
    min_size, max_size = 15, 320
    
    for ring_idx in range(len(distances)-1):
        r_center_pos = (distances[ring_idx] + distances[ring_idx+1]) / 2   # 环中径
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            all_ct_idx = all_celltypes.index(cell_type)
            base_angle = (all_ct_idx + 0.5) * 2 * np.pi / num_sectors
            for grp_idx in range(num_groups):
                angle = base_angle + offsets[grp_idx]
                density_val = densities[ring_idx, ct_idx, grp_idx]
                size = density_val * base_scale
                size = np.clip(size, min_size, max_size)
                if size > 0:
                    ax.scatter(angle, r_center_pos, s=size,
                               facecolor=colors[grp_idx], edgecolor='white',
                               linewidth=0.5, alpha=0.85, zorder=10)
    
    # 图例（右侧）
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=colors[i], markersize=8, label=all_groups[i])
                       for i in range(num_groups)]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.02, 1.0), fontsize=9, title='State' if num_groups>1 else 'Group',
              title_fontsize=10, framealpha=0.9, edgecolor='#CCCCCC')
    
    ax.set_ylim(0, 410)
    plt.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)
    
    # 保存图片
    plt.savefig('Epithelial_immune_plot_publication.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig('Epithelial_immune_plot_publication.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()

# ================== 5. 主程序 ==================
if __name__ == "__main__":
    adata, center_celltype, target_celltypes = load_and_filter_data()
    
    mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups = calculate_cell_counts(
        adata, center_celltype, target_celltypes
    )
    
    export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups,
                         filename='Epithelial_immune_cell_counts.csv')
    
    create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                       non_center_celltypes, all_groups)
