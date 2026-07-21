import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
import scanpy as sc
from scipy.spatial import KDTree
import math
from tqdm import tqdm

# ================== 字体设置 ==================
font_path = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 18

# ================== 辅助：长标签自动换行 ==================
def wrap_label(label):
    words = label.split()
    if len(words) <= 1:
        return label
    best_diff = float('inf')
    best_split = 1
    for i in range(1, len(words)):
        line1 = ' '.join(words[:i])
        line2 = ' '.join(words[i:])
        diff = abs(len(line1) - len(line2))
        if diff < best_diff:
            best_diff = diff
            best_split = i
    return ' '.join(words[:best_split]) + '\n' + ' '.join(words[best_split:])

# ================== 1. 数据加载 ==================
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter_state_samp_new_ctzq_with_CNs10-15.h5ad')
    center_celltype = "Epithelial"
    target_celltypes = [
        "Angiogenic EC", "Venule EC", "Lymphatic EC", "Arteriole EC",
        "ADAMDEC1+ Homeo-Immuno Fib", "Act Strom CAF", "matri-prod CAF",
        "RSPO3+ GREM1+ Crypt Fib", "WNT5B+ VEGFA+ Villus Fib",
        "Myofib", "apCAF", "Com Fib", "CXCL12+ Fib", "Pericyte-like Fib"
    ]
    
    if 'celltype_zq' not in adata.obs.columns:
        raise KeyError("缺少 'celltype_zq' 列")
    if 'celltype' not in adata.obs.columns:
        raise KeyError("缺少 'celltype' 列")
    
    center_mask = adata.obs['celltype_zq'] == center_celltype
    if center_mask.sum() == 0:
        raise ValueError("未找到 Epithelial 细胞")
    
    target_mask = adata.obs['celltype'].isin(target_celltypes) & (~center_mask)
    final_mask = center_mask | target_mask
    adata = adata[final_mask].copy()
    
    existing = adata.obs['celltype'].unique()
    missing = [ct for ct in target_celltypes if ct not in existing]
    if missing:
        print(f"警告: 以下细胞类型缺失: {', '.join(missing)}")
    
    print(f"中心细胞 (Epithelial): {center_mask.sum()} 个")
    print(f"周边目标细胞: {target_mask.sum()} 个")
    return adata, center_celltype, target_celltypes

# ================== 2. 细胞计数与密度计算 ==================
def calculate_cell_counts(adata, center_celltype, target_celltypes):
    center_mask = adata.obs['celltype_zq'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    if num_centers == 0:
        raise ValueError("无中心细胞坐标")
    
    other_mask = ~center_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    
    # 分组信息
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
    
    # 移除中心细胞类型，只保留实际存在的周边细胞
    non_center_celltypes = [ct for ct in target_celltypes if ct in other_celltypes.unique()]
    num_celltypes = len(non_center_celltypes)
    num_groups = len(all_groups)
    
    tree = KDTree(other_coords)
    distances = [0, 50, 100, 150, 200]  # 与之前保持一致
    
    counts = np.zeros((len(distances)-1, num_celltypes, num_groups))
    
    # 进度条包裹中心循环
    for center in tqdm(center_coords, desc="Processing center cells"):
        indices = tree.query_ball_point(center, r=distances[-1])
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
    
    sector_areas = np.zeros((len(distances)-1, num_celltypes))
    densities = np.zeros_like(mean_counts)
    for dist_idx in range(len(distances)-1):
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

# ================== 4. 径向气泡图（参考脚本风格） ==================
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                       all_groups, center_celltype):
    """
    扇区总数 = 周边细胞类型数 + 1（中心细胞占一个扇区，用于显示距离标签）
    气泡位于各环中径，同组气泡通过角度偏移区分。
    """
    # 构建完整的细胞类型列表（中心+周边），顺序保持中心在第一位
    all_celltypes = [center_celltype] + non_center_celltypes
    num_sectors = len(all_celltypes)          # 扇形总数
    sector_angle = 2 * np.pi / num_sectors    # 单个扇形的角度
    
    num_groups = len(all_groups)
    
    # ----- 颜色配置 -----
    if num_groups <= 10:
        group_colors = {grp: plt.cm.tab10(i % 10) for i, grp in enumerate(all_groups)}
    else:
        group_colors = {grp: plt.cm.gist_rainbow(i / num_groups) for i, grp in enumerate(all_groups)}
    colors = [group_colors[g] for g in all_groups]
    
    # ----- 角度偏移（每个分组在扇形内的偏移量） -----
    if num_groups == 1:
        angle_offsets = [0]
    else:
        max_offset = sector_angle * 0.4   # 占扇形宽度的80%（两侧各40%）
        angle_offsets = np.linspace(-max_offset, max_offset, num_groups).tolist()
    
    # ----- 画布设置 -----
    fig = plt.figure(figsize=(8, 8), dpi=300, facecolor='white')
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    
    # ----- 中心圆形区域（半径50，对应距离0-50）-----
    theta_center = np.linspace(0, 2*np.pi, 100)
    ax.fill_between(theta_center, 0, distances[1], color='#EEEEEE', alpha=0.6, zorder=2)
    ax.plot(theta_center, [distances[1]]*100, color='#B0B0B0', linewidth=0.8, linestyle='-', alpha=0.8, zorder=3)
    
    # 中心细胞标签置于圆心
    ax.text(0, 0, center_celltype, ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9), zorder=4)
    
    # ----- 同心圆环（距离边界线）-----
    for r in distances[2:]:  # 跳过0和50，因为50已画为中心圆边界
        ax.plot(np.linspace(0, 2*np.pi, 360), [r]*360,
                color='#D3D3D3', linewidth=0.5, linestyle='-', alpha=0.6, zorder=1)
    
    # ----- 扇形分隔虚线（从半径50到385）-----
    for i in range(num_sectors):
        angle = i * sector_angle
        ax.plot([angle, angle], [distances[1], 385], color='#B0B0B0',
                linewidth=0.4, linestyle='--', alpha=0.5, zorder=1)
    
    # ----- 扇形背景（外圈）-----
    sector_bg_colors = plt.cm.Greys(np.linspace(0.1, 0.25, num_sectors))
    for i in range(num_sectors):
        theta_start = i * sector_angle
        theta_end = (i+1) * sector_angle
        ax.fill_between(np.linspace(theta_start, theta_end, 50),
                        370, 400, color=sector_bg_colors[i], alpha=0.15, zorder=0)
    
    # ----- 外圈细胞类型标签（周边细胞）-----
    label_radius = 385
    for i, ct in enumerate(all_celltypes):
        if ct == center_celltype:
            continue   # 中心细胞标签已在圆心
        angle = (i + 0.5) * sector_angle
        wrapped = wrap_label(ct)          # 长名称自动换行
        ax.text(angle, label_radius, wrapped, ha='center', va='center',
                fontsize=12, weight='bold', color='#333333',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          edgecolor='none', alpha=0.7))
    
    # ----- 距离标签（放在中心细胞所在扇区内）-----
    center_idx = all_celltypes.index(center_celltype)   # 应为0
    center_angle = (center_idx + 0.5) * sector_angle
    for ring_idx in range(len(distances)-1):
        r_inner = distances[ring_idx]
        r_outer = distances[ring_idx+1]
        r_mid = (r_inner + r_outer) / 2
        label = f"{r_inner}–{r_outer}"
        ax.text(center_angle, r_mid, label, ha='center', va='center',
                fontsize=9, weight='bold', color='#555555',
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.8))
    # 单位标注
    ax.text(center_angle, 375, '(μm)', ha='center', va='center',
            fontsize=10, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.7))
    
    # ----- 气泡绘制 -----
    max_density = np.max(densities) if np.max(densities) > 0 else 1
    base_scale = 400 / max_density   # 可根据效果调整
    min_size, max_size = 15, 220
    
    for ring_idx in range(len(distances)-1):
        r_inner = distances[ring_idx]
        r_outer = distances[ring_idx+1]
        r_center = (r_inner + r_outer) / 2   # 环的中径
        
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            # 在all_celltypes中的索引（因为前面加了center，索引需+1）
            all_ct_idx = ct_idx + 1
            base_angle = (all_ct_idx + 0.5) * sector_angle
            
            for grp_idx, grp in enumerate(all_groups):
                angle = base_angle + angle_offsets[grp_idx]
                density_val = densities[ring_idx, ct_idx, grp_idx]
                size = density_val * base_scale
                size = np.clip(size, min_size, max_size)
                if size > 0:
                    ax.scatter(angle, r_center, s=size,
                               facecolor=colors[grp_idx], edgecolor='white',
                               linewidth=0.5, alpha=0.85, zorder=10)
    
    # ----- 图例 -----
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=group_colors[grp], markersize=8, label=grp)
                       for grp in all_groups]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.02, 1.0), fontsize=9, title='State',
              title_fontsize=10, framealpha=0.9, edgecolor='#CCCCCC')
    
    # ----- 范围与保存 -----
    ax.set_ylim(0, 410)
    plt.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)
    
    safe_name = center_celltype.replace(" ", "_").replace("/", "_").replace("+", "plus")
    plt.savefig(f'{safe_name}_publication.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(f'{safe_name}_publication.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()

# ================== 5. 主程序 ==================
if __name__ == "__main__":
    adata, center_celltype, target_celltypes = load_and_filter_data()
    
    # 可以针对多个中心循环，此处仅一个中心
    for center in [center_celltype]:
        print(f"\n===== 处理中心细胞: {center} =====")
        mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups = \
            calculate_cell_counts(adata, center, target_celltypes)
        export_counts_to_csv(mean_counts, densities, distances, sector_areas,
                             non_center_celltypes, all_groups, center)
        create_radial_plot(mean_counts, densities, distances, non_center_celltypes,
                           all_groups, center)
    
    print("\n所有图表生成完毕！")
