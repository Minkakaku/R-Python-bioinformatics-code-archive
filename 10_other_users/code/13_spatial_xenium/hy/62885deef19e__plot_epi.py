import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import math
import os
from matplotlib import font_manager

# ================== 字体加载（请根据实际路径修改） ==================
try:
    font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')
except:
    print("警告: 未能加载 Arial 字体，将使用系统默认字体。")

# ================== 全局绘图设置 ==================
plt.rcParams['font.family'] = 'arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# 输出目录
OUTPUT_DIR = 'Epi_immu'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ================== 1. 数据加载 ==================
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter.h5ad')
    
    center_celltypes = [
        "Proliferative Epithelials",
        "Hypoxic Epithelials",
        "IFN-Epithelials",
        "TA-like Epithelials",
        "Wnt-active stem",
        "Stem-like Tumor"
    ]
    neighbor_celltypes = [
        "CD8+ T",
        "CXCL9+ CD8+ T",
        "MKI67+ T",
        "CD4+ T",
        "Treg",
        "CXCL9+ M1",
        "M2",
        "M1"
    ]
    all_celltypes = center_celltypes + neighbor_celltypes
    
    mask = adata.obs['celltype'].isin(all_celltypes)
    adata = adata[mask].copy()
    
    existing = adata.obs['celltype'].unique()
    missing = [ct for ct in all_celltypes if ct not in existing]
    if missing:
        print(f"警告: 以下细胞类型不存在: {', '.join(missing)}")
    
    return adata, center_celltypes, neighbor_celltypes

# ================== 2. 细胞计数与密度计算 ==================
def calculate_cell_counts(adata, center_celltype, neighbor_celltypes, distances):
    center_mask = adata.obs['celltype'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    if num_centers == 0:
        print(f"警告: 未找到 {center_celltype}，使用全部细胞重心")
        center_coords = np.mean(adata.obsm['spatial'], axis=0).reshape(1, -1)
        num_centers = 1
    
    all_groups = ['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R']
    group_to_idx = {g: i for i, g in enumerate(all_groups)}
    
    other_mask = ~center_mask
    other_obs = adata[other_mask].obs
    valid_group = other_obs['group_final'].isin(all_groups)
    other_mask_indices = other_mask[other_mask].index[valid_group]
    other_mask = adata.obs.index.isin(other_mask_indices)
    
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    other_groups = adata[other_mask].obs['group_final']
    
    neighbor_mask = other_celltypes.isin(neighbor_celltypes)
    other_coords = other_coords[neighbor_mask]
    other_celltypes = other_celltypes[neighbor_mask]
    other_groups = other_groups[neighbor_mask]
    
    other_celltypes = other_celltypes.reset_index(drop=True)
    other_groups = other_groups.reset_index(drop=True)
    other_coords = np.array(other_coords)
    
    num_rings = len(distances) - 1
    num_celltypes = len(neighbor_celltypes)
    num_groups = len(all_groups)
    counts = np.zeros((num_rings, num_celltypes, num_groups))
    
    if len(other_coords) > 0:
        tree = KDTree(other_coords)
    else:
        print(f"警告: {center_celltype} 周围没有符合条件的细胞")
        tree = None
    
    if tree is not None:
        for center in center_coords:
            indices = tree.query_ball_point(center, r=250)
            for idx in indices:
                dist = np.linalg.norm(center - other_coords[idx])
                ct = other_celltypes.iloc[idx]
                grp = other_groups.iloc[idx]
                if ct not in neighbor_celltypes:
                    continue
                ct_idx = neighbor_celltypes.index(ct)
                grp_idx = group_to_idx[grp]
                for i in range(num_rings):
                    if distances[i] <= dist < distances[i+1]:
                        counts[i, ct_idx, grp_idx] += 1
                        break
    
    mean_counts = counts / num_centers if num_centers > 0 else counts
    
    sector_areas = np.zeros((num_rings, num_celltypes))
    densities = np.zeros_like(mean_counts)
    for dist_idx in range(num_rings):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / len(neighbor_celltypes)
        for ct_idx in range(num_celltypes):
            sector_areas[dist_idx, ct_idx] = sector_area
            for grp_idx in range(num_groups):
                densities[dist_idx, ct_idx, grp_idx] = mean_counts[dist_idx, ct_idx, grp_idx] / sector_area
    
    return mean_counts, densities, distances, sector_areas, neighbor_celltypes, all_groups

# ================== 3. 导出CSV ==================
def export_counts_to_csv(mean_counts, densities, distances, sector_areas, 
                         neighbor_celltypes, all_groups, center_celltype):
    safe_name = center_celltype.replace(" ", "_").replace("/", "_").replace("+", "plus")
    filename = os.path.join(OUTPUT_DIR, f'{safe_name}_cell_counts.csv')
    dist_labels = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    data = []
    for ct_idx, cell_type in enumerate(neighbor_celltypes):
        for dist_idx in range(len(dist_labels)):
            area = sector_areas[dist_idx, ct_idx]
            for grp_idx, grp_name in enumerate(all_groups):
                data.append([
                    cell_type, dist_labels[dist_idx], grp_name,
                    mean_counts[dist_idx, ct_idx, grp_idx],
                    area, densities[dist_idx, ct_idx, grp_idx]
                ])
    df = pd.DataFrame(data, columns=[
        'celltype', 'distance_range', 'group', 'cell_count', 'sector_area', 'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"数据已导出: {filename}")

# ================== 4. 绘图函数 ==================
def create_radial_plot(mean_counts, densities, distances, neighbor_celltypes,
                       all_groups, center_celltype):
    num_groups = len(all_groups)
    num_sectors = len(neighbor_celltypes)   # 扇区数 = 周围细胞类型数
    sector_angle = 2 * np.pi / num_sectors
    
    group_colors = {
        'POST-NR': '#e41a1c',
        'POST-R':  '#377eb8',
        'PRE-NR':  '#4daf4a',
        'PRE-R':   '#984ea3'
    }
    colors = [group_colors[g] for g in all_groups]
    
    spread_ratio = 0.7
    offsets = np.linspace(-spread_ratio/2, spread_ratio/2, num_groups) * sector_angle
    angle_offset = {grp: offsets[i] for i, grp in enumerate(all_groups)}
    
    fig = plt.figure(figsize=(10, 10), dpi=300, facecolor='white')
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
    
    ct_lines = center_celltype.split(' ')
    ct_text = '\n'.join(ct_lines) if len(ct_lines) > 1 else center_celltype
    ax.text(0, 0, ct_text, ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9), zorder=4)
    
    # 同心圆环（边界线）
    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 360), [r]*360,
                color='#D3D3D3', linewidth=0.5, linestyle='-', alpha=0.6, zorder=1)
    
    # 扇形分隔虚线（只针对周围细胞）
    for i in range(num_sectors):
        angle = i * 2 * np.pi / num_sectors
        ax.plot([angle, angle], [50, 360], color='#B0B0B0',
                linewidth=0.4, linestyle='--', alpha=0.5, zorder=1)
    
    # 扇形外圈背景（半径360-390）
    sector_bg_colors = plt.cm.Greys(np.linspace(0.1, 0.25, num_sectors))
    for i in range(num_sectors):
        theta_start = i * 2 * np.pi / num_sectors
        theta_end = (i+1) * 2 * np.pi / num_sectors
        ax.fill_between(np.linspace(theta_start, theta_end, 50),
                        360, 390, color=sector_bg_colors[i], alpha=0.15, zorder=0)
    
    # 外圈细胞类型标签（只显示周围细胞）
    for i, ct in enumerate(neighbor_celltypes):
        angle = (i + 0.5) * 2 * np.pi / num_sectors
        ax.text(angle, 360, ct, ha='center', va='center',
                fontsize=9, weight='bold', color='#333333',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          edgecolor='none', alpha=0.7))
    
    # ----- 距离标签（修改处：改成放在对应的圆形边界线上） -----
    label_angle = 0  # 12点钟方向
    for ring_idx in range(len(distances)-1):
        # 原：r_mid = (distances[ring_idx] + distances[ring_idx+1]) / 2
        # 改：使用外边界半径，让标签精确落在对应的圆圈线上
        r_line = distances[ring_idx+1]  
        
        label = f"{distances[ring_idx]}–{distances[ring_idx+1]}"
        ax.text(label_angle, r_line, label, ha='center', va='center',
                fontsize=8, weight='bold', color='#555555',
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.8))
    # 单位 μm 放在稍外（不变）
    ax.text(label_angle, 380, '(μm)', ha='center', va='center',
            fontsize=9, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.7))
    
    # 绘制气泡（只对 neighbor_celltypes）
    max_density = np.max(densities) if np.max(densities) > 0 else 1
    base_scale = 400 / max_density
    min_size, max_size = 15, 220
    
    for ring_idx in range(len(distances)-1):
        r_inner = distances[ring_idx]
        r_outer = distances[ring_idx+1]
        r_center = (r_inner + r_outer) / 2
        
        for ct_idx, cell_type in enumerate(neighbor_celltypes):
            base_angle = (ct_idx + 0.5) * 2 * np.pi / num_sectors
            for grp_idx, grp in enumerate(all_groups):
                angle = base_angle + angle_offset[grp]
                density_val = densities[ring_idx, ct_idx, grp_idx]
                size = density_val * base_scale
                size = np.clip(size, min_size, max_size)
                if size > 0:
                    ax.scatter(angle, r_center, s=size,
                               facecolor=colors[grp_idx], edgecolor='white',
                               linewidth=0.5, alpha=0.85, zorder=10)
    
    # 图例
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=color, markersize=8, label=grp)
                       for grp, color in zip(all_groups, colors)]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.02, 1.0), fontsize=9, title='Group',
              title_fontsize=10, framealpha=0.9, edgecolor='#CCCCCC')
    
    ax.set_ylim(0, 400)
    plt.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)
    
    safe_name = center_celltype.replace(" ", "_").replace("/", "_").replace("+", "plus")
    pdf_path = os.path.join(OUTPUT_DIR, f'{safe_name}_publication.pdf')
    png_path = os.path.join(OUTPUT_DIR, f'{safe_name}_publication.png')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close()
    print(f"图形已保存: {pdf_path}")

# ================== 5. 主程序 ==================
if __name__ == "__main__":
    adata, center_celltypes, neighbor_celltypes = load_and_filter_data()
    distances = [50, 100, 150, 200, 250]
    
    for center in center_celltypes:
        print(f"\n===== 处理中心细胞: {center} =====")
        mean_counts, densities, distances, sector_areas, neighbor_celltypes, all_groups = \
            calculate_cell_counts(adata, center, neighbor_celltypes, distances)
        export_counts_to_csv(mean_counts, densities, distances, sector_areas,
                             neighbor_celltypes, all_groups, center)
        create_radial_plot(mean_counts, densities, distances, neighbor_celltypes,
                           all_groups, center)
    
    print("\n所有图表生成完毕！")
