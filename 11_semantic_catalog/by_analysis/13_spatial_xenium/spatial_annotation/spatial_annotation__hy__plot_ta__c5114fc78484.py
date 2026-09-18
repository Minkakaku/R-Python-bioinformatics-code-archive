import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scanpy as sc
from scipy.spatial import KDTree
import math
import os
import itertools

from matplotlib import font_manager

font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')

# ================== 全局设置 ==================
plt.rcParams['font.family'] = 'arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ================== 1. 数据加载（不进行任何筛选，只读入）==================
def load_full_data(file_path):
    """加载完整h5ad文件（不做细胞类型筛选）"""
    adata = sc.read_h5ad(file_path)
    return adata

# ================== 2. 细胞计数与密度计算（完全基于传入的adata）==================
def calculate_cell_counts(adata, target_celltypes, center_celltype, group_col='group_final'):
    """
    计算中心细胞周围环形区域内的细胞计数与密度。
    要求传入的 adata 已经过预处理：只包含需要对比的两个分组的所有细胞（中心+周围）。
    """
    # 提取中心细胞坐标
    center_mask = adata.obs['celltype'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    if num_centers == 0:
        print(f"警告: 在子集中未找到 {center_celltype}，跳过该对比。")
        return None, None, None, None, None, None

    # 周围细胞（非中心细胞）
    other_mask = ~center_mask
    other_obs = adata[other_mask].obs
    # 从周围细胞中提取分组（应该只有两个分组）
    all_groups = sorted(other_obs[group_col].dropna().unique())
    if len(all_groups) != 2:
        print(f"警告: 周围细胞分组数量为 {len(all_groups)}，期望2组，跳过。")
        return None, None, None, None, None, None

    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    other_state = adata[other_mask].obs[group_col]

    group_to_idx = {g: i for i, g in enumerate(all_groups)}
    other_group = other_state.map(group_to_idx).values

    tree = KDTree(other_coords)
    distances = [50, 100, 150, 200, 250]

    non_center_celltypes = [ct for ct in target_celltypes if ct != center_celltype]
    num_celltypes = len(non_center_celltypes)
    num_groups = len(all_groups)

    counts = np.zeros((len(distances)-1, num_celltypes, num_groups))

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

    mean_counts = counts / num_centers

    total_types = len(target_celltypes)
    sector_areas = np.zeros((len(distances)-1, num_celltypes))
    densities = np.zeros_like(mean_counts)
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / total_types
        for ct_idx in range(num_celltypes):
            sector_areas[dist_idx, ct_idx] = sector_area
            for grp_idx in range(num_groups):
                densities[dist_idx, ct_idx, grp_idx] = mean_counts[dist_idx, ct_idx, grp_idx] / sector_area

    return mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups

# ================== 3. 导出CSV ==================
def export_counts_to_csv(mean_counts, densities, distances, sector_areas,
                         non_center_celltypes, all_groups, center_celltype, batch_name, output_dir):
    filename = os.path.join(output_dir, f'{batch_name}_cell_counts.csv')
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
        'celltype', 'distance_range', 'group_final', 'cell_count', 'sector_area', 'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"数据已导出: {filename}")

# ================== 4. 绘图函数（不变）==================
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                       all_celltypes, all_groups, center_celltype, batch_name, output_dir):
    num_groups = len(all_groups)
    num_sectors = len(all_celltypes)
    sector_angle = 2 * np.pi / num_sectors

    if num_groups <= 10:
        cmap = cm.get_cmap('tab10', num_groups)
    else:
        cmap = cm.get_cmap('tab20', num_groups)
    group_colors = {grp: cmap(i) for i, grp in enumerate(all_groups)}
    colors = [group_colors[g] for g in all_groups]

    if num_groups == 1:
        angle_offset = {all_groups[0]: 0}
    else:
        offset_range = 0.4
        offsets = np.linspace(-offset_range * sector_angle / 2,
                              offset_range * sector_angle / 2, num_groups)
        angle_offset = {grp: offsets[i] for i, grp in enumerate(all_groups)}

    fig = plt.figure(figsize=(7, 7), dpi=300, facecolor='white')
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)

    theta_center = np.linspace(0, 2*np.pi, 100)
    ax.fill_between(theta_center, 0, 50, color='#EEEEEE', alpha=0.6, zorder=2)
    ax.plot(theta_center, [50]*100, color='#B0B0B0', linewidth=0.8, linestyle='-', alpha=0.8, zorder=3)

    ct_lines = center_celltype.split(' ')
    ct_text = '\n'.join(ct_lines) if len(ct_lines) > 1 else center_celltype
    ax.text(0, 0, ct_text, ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9), zorder=4)

    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 360), [r]*360,
                color='#D3D3D3', linewidth=0.5, linestyle='-', alpha=0.6, zorder=1)

    for i in range(num_sectors):
        angle = i * 2 * np.pi / num_sectors
        ax.plot([angle, angle], [50, 360], color='#B0B0B0',
                linewidth=0.4, linestyle='--', alpha=0.5, zorder=1)

    sector_bg_colors = plt.cm.Greys(np.linspace(0.1, 0.25, num_sectors))
    for i in range(num_sectors):
        theta_start = i * 2 * np.pi / num_sectors
        theta_end = (i+1) * 2 * np.pi / num_sectors
        ax.fill_between(np.linspace(theta_start, theta_end, 50),
                        360, 390, color=sector_bg_colors[i], alpha=0.15, zorder=0)

    for i, ct in enumerate(all_celltypes):
        if ct == center_celltype:
            continue
        angle = (i + 0.5) * 2 * np.pi / num_sectors
        display_ct = ct if len(ct) <= 20 else ct[:17] + '...'
        ax.text(angle, 370, display_ct, ha='center', va='center',
                fontsize=7, weight='bold', color='#333333',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          edgecolor='none', alpha=0.7))

    center_idx = all_celltypes.index(center_celltype)
    center_angle = (center_idx + 0.5) * 2 * np.pi / num_sectors
    for ring_idx in range(len(distances)-1):
        r_mid = (distances[ring_idx] + distances[ring_idx+1]) / 2
        label = f"{distances[ring_idx]}–{distances[ring_idx+1]}"
        ax.text(center_angle, r_mid, label, ha='center', va='center',
                fontsize=9, weight='bold', color='#555555',
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.8))
    ax.text(center_angle, 370, '(μm)', ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.7))

    max_density = np.max(densities) if np.max(densities) > 0 else 1
    base_scale = 400 / max_density
    min_size, max_size = 15, 220

    for ring_idx in range(len(distances)-1):
        r_center_ring = (distances[ring_idx] + distances[ring_idx+1]) / 2
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            all_ct_idx = all_celltypes.index(cell_type)
            base_angle = (all_ct_idx + 0.5) * 2 * np.pi / num_sectors
            for grp_idx, grp in enumerate(all_groups):
                angle = base_angle + angle_offset[grp]
                density_val = densities[ring_idx, ct_idx, grp_idx]
                size = density_val * base_scale
                size = np.clip(size, min_size, max_size)
                if size > 0:
                    ax.scatter(angle, r_center_ring, s=size,
                               facecolor=colors[grp_idx], edgecolor='white',
                               linewidth=0.5, alpha=0.85, zorder=10)

    legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=group_colors[grp], markersize=8, label=grp)
                       for grp in all_groups]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.02, 1.0), fontsize=9, title='Group (group_final)',
              title_fontsize=10, framealpha=0.9, edgecolor='#CCCCCC')
    ax.set_ylim(0, 400)
    plt.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)

    pdf_path = os.path.join(output_dir, f'{batch_name}_publication.pdf')
    png_path = os.path.join(output_dir, f'{batch_name}_publication.png')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()

# ================== 5. 主程序：先按分组提取，再按细胞类型提取，然后取两个分组子集合并 ==================
if __name__ == "__main__":
    input_file = 'All_filter_with_region_CNs10.h5ad'
    output_dir = 'TA-like_Epithelials_dis'
    os.makedirs(output_dir, exist_ok=True)

    center_cell = 'TA-like Epithelials'

    # 定义批次和周围细胞类型
    batch1_surrounding = [
        ' myCAF','iCAF','Matrix Fibroblasts','Matrix remodeling CAFs','cycling CAF','Tip cells','MKI67+ Endo','EndoMT','Activated Endo','Microvascular Endo','Arterial Endo'
    ]
    batch2_surrounding = [
        'CD8+ T', 'CXCL9+ CD8+ T', 'MKI67+ T', 'Treg',
        'LAMP3+ DC', 'CD4+T', 'M2',
        'Monocyte-derived Macrophages', 'CXCL9+ M1', 'M1'
    ]
    batches = {
        'TA-like_stromal': batch1_surrounding,
        'TA-like_immu': batch2_surrounding
    }

    # 要对比的所有分组（四个）
    all_groups_to_compare = ["PRE-R", "PRE-NR", "POST-R", "POST-NR"]

    # 第一步：加载完整数据（不筛选细胞类型，保留所有细胞）
    print("加载完整数据...")
    adata_full = load_full_data(input_file)
    print(f"原始数据形状: {adata_full.shape}")
    print(f"可用的分组: {sorted(adata_full.obs['group_final'].unique())}")
    print(f"可用的细胞类型: {adata_full.obs['celltype'].nunique()} 种")

    # 第二步：对于每个批次，先按分组提取子集（每个分组子集包含该分组下的所有细胞）
    for batch_name, surrounding_list in batches.items():
        print(f"\n===== 处理批次: {batch_name} =====")
        
        # 该批次需要的所有细胞类型（中心+周围）
        target_celltypes = [center_cell] + surrounding_list
        
        # 预先按分组提取四个子集（字典：分组名 -> AnnData）
        group_subsets = {}
        for grp in all_groups_to_compare:
            # 提取属于该分组的所有细胞（不限制细胞类型）
            mask_grp = adata_full.obs['group_final'] == grp
            if not np.any(mask_grp):
                print(f"警告: 分组 {grp} 没有数据，跳过")
                continue
            adata_grp = adata_full[mask_grp].copy()
            # 再筛选出该批次需要的细胞类型（中心细胞+周围细胞）
            mask_ct = adata_grp.obs['celltype'].isin(target_celltypes)
            adata_grp = adata_grp[mask_ct].copy()
            if adata_grp.n_obs > 0:
                group_subsets[grp] = adata_grp
                print(f"  分组 {grp} 提取到 {adata_grp.n_obs} 个细胞（包含中心细胞和周围细胞）")
        
        # 获取实际存在的分组
        existing_groups = list(group_subsets.keys())
        if len(existing_groups) < 2:
            print(f"  有效分组不足2个，跳过该批次。")
            continue
        
        # 生成所有两两组合
        combos = list(itertools.combinations(existing_groups, 2))
        print(f"  将生成以下分组组合图: {combos}")
        
        for grp1, grp2 in combos:
            combo_name = f"{grp1}_vs_{grp2}"
            print(f"    --- 合并分组: {combo_name} ---")
            
            # 合并两个分组的子集（直接 concatenate）
            adata_combined = sc.concat([group_subsets[grp1], group_subsets[grp2]], join='outer', index_unique=None)
            # 确保索引唯一（AnnData 的 concat 可能生成重复索引，这里不管）
            # 此时 adata_combined 中包含 grp1 和 grp2 的所有细胞（中心细胞+周围细胞）
            # 注意：中心细胞 TA-like Epithelials 只来源于这两个分组，不会重复（因为每个细胞只属于一个分组）
            
            # 检查中心细胞是否存在
            if not np.any(adata_combined.obs['celltype'] == center_cell):
                print(f"      合并后的子集中没有中心细胞 {center_cell}，跳过")
                continue
            
            # 此时 all_celltypes 应该包括中心细胞和实际存在的周围细胞类型
            actual_celltypes = adata_combined.obs['celltype'].unique().tolist()
            # 确保按 target_celltypes 顺序，但只保留实际存在的
            ordered_targets = [ct for ct in target_celltypes if ct in actual_celltypes]
            if len(ordered_targets) <= 1:
                print("      只有中心细胞，无周围细胞，跳过")
                continue
            
            # 调用统计函数（传入已合并的子集）
            result = calculate_cell_counts(adata_combined, ordered_targets, center_cell, group_col='group_final')
            if result[0] is None:
                continue
            mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups = result
            
            # 导出和绘图
            full_batch_name = f"{batch_name}_{combo_name}"
            export_counts_to_csv(mean_counts, densities, distances, sector_areas,
                                 non_center_celltypes, all_groups, center_cell,
                                 full_batch_name, output_dir)
            create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                               ordered_targets, all_groups, center_cell,
                               full_batch_name, output_dir)
    
    print("\n所有分析完成，结果保存于 TA-like Epithelials_dis/ 文件夹。")
