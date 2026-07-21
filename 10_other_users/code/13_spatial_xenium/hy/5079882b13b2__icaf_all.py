import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
import scanpy as sc
from scipy.spatial import KDTree
import math
from tqdm import tqdm

# ================== 字体 ==================
font_path = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 18

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

# ================== 1. 数据加载（不挑区域） ==================
def load_and_filter_data(adata_path):
    adata = sc.read_h5ad(adata_path)
    required = ['site_band', 'celltype', 'group_final']
    for c in required:
        if c not in adata.obs.columns:
            raise KeyError(f"缺少列: {c}")

    print("=== 数据诊断 ===")
    print("site_band:", adata.obs['site_band'].unique())
    print("celltype 前20:", adata.obs['celltype'].unique()[:20])
    print("group_final:", adata.obs['group_final'].unique())

    # 中心细胞：所有 iCAF，不限制 site_band
    center_celltype = "iCAF"
    center_mask = adata.obs['celltype'] == center_celltype

    # 周边免疫细胞列表（根据真实名称调整）
    target_celltypes = [
        "CXCL9+ CD8+ T", "MKI67+ T", "Treg", "LAMP3+ DC",
        "CD4+ T", "M2", "Monocyte-derived Macrophages",
        "CXCL9+ M1", "M1","CD8+ T"
    ]

    # 周边：所有部位中属于目标细胞列表，且不是中心细胞
    peri_mask = adata.obs['celltype'].isin(target_celltypes) & (~center_mask)

    final_mask = center_mask | peri_mask
    adata_filt = adata[final_mask].copy()

    # 再次确认周边细胞（仅用于诊断）
    peri_mask_filt = adata_filt.obs['celltype'].isin(target_celltypes) & \
                     (adata_filt.obs['celltype'] != center_celltype)

    existing = adata_filt[peri_mask_filt].obs['celltype'].unique()
    missing = [ct for ct in target_celltypes if ct not in existing]
    if missing:
        print(f"⚠️ 缺失类型: {missing}")

    print(f"中心细胞 ({center_celltype}, 全区域): {center_mask.sum()}")
    print(f"周边细胞 (全区域): {peri_mask.sum()}")
    return adata_filt, center_celltype, [ct for ct in target_celltypes if ct in existing]

# ================== 2. 计数与密度（限定分组 PRE-R / PRE-NR） ==================
def calculate_cell_counts(adata, center_celltype, target_celltypes):
    # 中心细胞：所有 iCAF
    center_mask = adata.obs['celltype'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    n_center = len(center_coords)
    if n_center == 0:
        raise ValueError("无中心细胞")

    other_mask = ~center_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_ct = adata[other_mask].obs['celltype']

    # ===== 分组筛选：仅保留 PRE-R 和 PRE-NR =====
    valid_groups = ['PRE-R', 'PRE-NR']
    other_grp_all = adata[other_mask].obs['group_final'].astype(str)
    group_filter = other_grp_all.isin(valid_groups)

    other_coords = other_coords[group_filter]
    other_ct = other_ct[group_filter]
    other_grp = other_grp_all[group_filter]

    all_groups = [g for g in valid_groups if g in other_grp.unique()]
    if len(all_groups) == 0:
        raise ValueError("没有符合 PRE-R 或 PRE-NR 的周边细胞")
    grp2idx = {g: i for i, g in enumerate(all_groups)}
    other_grp_idx = other_grp.map(grp2idx).values

    non_center_ct = [ct for ct in target_celltypes if ct in other_ct.unique()]
    n_ct = len(non_center_ct)
    n_grp = len(all_groups)

    distances = [0, 50, 100, 150, 200]
    counts = np.zeros((len(distances)-1, n_ct, n_grp))
    tree = KDTree(other_coords)

    for center in tqdm(center_coords, desc="Processing iCAF centers"):
        idxs = tree.query_ball_point(center, r=distances[-1])
        for idx in idxs:
            d = np.linalg.norm(center - other_coords[idx])
            ct = other_ct.iloc[idx]
            if ct not in non_center_ct:
                continue
            ct_idx = non_center_ct.index(ct)
            g_idx = other_grp_idx[idx]
            for i in range(len(distances)-1):
                if distances[i] <= d < distances[i+1]:
                    counts[i, ct_idx, g_idx] += 1
                    break

    mean_counts = counts / n_center
    sector_areas = np.zeros((len(distances)-1, n_ct))
    densities = np.zeros_like(mean_counts)
    for i in range(len(distances)-1):
        r1, r2 = distances[i], distances[i+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sec_area = ring_area / n_ct
        for j in range(n_ct):
            sector_areas[i, j] = sec_area
            for k in range(n_grp):
                densities[i, j, k] = mean_counts[i, j, k] / sec_area if sec_area > 0 else 0

    return mean_counts, densities, distances, sector_areas, non_center_ct, all_groups

# ================== 3. 导出CSV ==================
def export_csv(mean_counts, densities, distances, sector_areas, non_center_ct, all_groups, center_label):
    safe = center_label.replace(" ", "_").replace("/", "_").replace("+", "plus")
    fname = f"{safe}_PRE-R_vs_PRE-NR.csv"   # 新文件名
    dl = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    rows = []
    for ct_i, ct in enumerate(non_center_ct):
        for d_i in range(len(dl)):
            area = sector_areas[d_i, ct_i]
            for g_i, g in enumerate(all_groups):
                rows.append([ct, dl[d_i], g, mean_counts[d_i, ct_i, g_i], area, densities[d_i, ct_i, g_i]])
    df = pd.DataFrame(rows, columns=['celltype','distance_range','group_final','cell_count','sector_area','cell_density'])
    df.to_csv(fname, sep='\t', index=False)
    print(f"CSV → {fname}")

# ================== 4. 径向图（缩小气泡 + PRE分组） ==================
def plot_radial(mean_counts, densities, distances, non_center_ct, all_groups, center_label):
    all_types = [center_label] + non_center_ct
    n_sec = len(all_types)
    sec_angle = 2 * np.pi / n_sec
    n_grp = len(all_groups)

    colors_dict = {
        'PRE-R': '#D62728',
        'PRE-NR': '#1F77B4'
    }
    for g in all_groups:
        if g not in colors_dict:
            colors_dict[g] = plt.cm.tab10(len(colors_dict) % 10)
    grp_colors = [colors_dict[g] for g in all_groups]

    if n_grp == 1:
        offsets = [0]
    else:
        max_off = sec_angle * 0.3
        offsets = np.linspace(-max_off, max_off, n_grp).tolist()

    fig = plt.figure(figsize=(8,8), dpi=300, facecolor='white')
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)

    # 圆心标签
    ax.text(0, 0, center_label, ha='center', va='center',
            fontsize=12, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9), zorder=4)

    # 圆环边界
    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 360), [r]*360,
                color='#D3D3D3', linewidth=0.5, linestyle='-', alpha=0.6, zorder=1)

    # 扇形分隔线
    for i in range(n_sec):
        angle = i * sec_angle
        ax.plot([angle, angle], [0, 200], color='#B0B0B0',
                linewidth=0.4, linestyle='--', alpha=0.5, zorder=1)

    # 扇形背景
    bg_colors = plt.cm.Greys(np.linspace(0.1, 0.25, n_sec))
    for i in range(n_sec):
        t0 = i * sec_angle
        t1 = (i+1) * sec_angle
        ax.fill_between(np.linspace(t0, t1, 50), 200, 220,
                        color=bg_colors[i], alpha=0.15, zorder=0)

    # 外圈细胞类型标签
    label_rad = 210
    for i, ct in enumerate(all_types):
        if ct == center_label:
            continue
        angle = (i + 0.5) * sec_angle
        wrapped = wrap_label(ct)
        ax.text(angle, label_rad, wrapped, ha='center', va='center',
                fontsize=11, weight='bold', color='#333333',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          edgecolor='none', alpha=0.7))

    # 距离标签（放在中心细胞扇区）
    center_idx = all_types.index(center_label)
    center_angle = (center_idx + 0.5) * sec_angle
    for d_i in range(len(distances)-1):
        r_in = distances[d_i]
        r_out = distances[d_i+1]
        r_mid = (r_in + r_out) / 2
        ax.text(center_angle, r_mid, f"{r_in}–{r_out}", ha='center', va='center',
                fontsize=8, weight='bold', color='#555555',
                bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                          edgecolor='none', alpha=0.8))
    ax.text(center_angle, 195, '(μm)', ha='center', va='center',
            fontsize=9, weight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.7))

    # ===== 气泡（更小，带约束） =====
    max_dens = np.max(densities) if np.max(densities) > 0 else 1
    base_area_factor = 200000
    min_area_threshold = 3

    for ring_i in range(len(distances)-1):
        r_in = distances[ring_i]
        r_out = distances[ring_i+1]
        r_mid = (r_in + r_out) / 2
        ring_width = r_out - r_in
        arc_length = r_mid * sec_angle
        max_diameter = min(ring_width * 0.7, arc_length * 0.25)
        max_area = math.pi * (max_diameter/2)**2

        for ct_i in range(len(non_center_ct)):
            all_ct_i = ct_i + 1
            base_angle = (all_ct_i + 0.5) * sec_angle
            for g_i in range(n_grp):
                dens = densities[ring_i, ct_i, g_i]
                area = dens * base_area_factor
                area = min(area, max_area)
                if area >= min_area_threshold:
                    angle = base_angle + offsets[g_i]
                    ax.scatter(angle, r_mid, s=area,
                               facecolor=grp_colors[g_i], edgecolor='white',
                               linewidth=0.3, alpha=0.85, zorder=10)

    # 图例
    legend_elems = [plt.Line2D([0],[0], marker='o', color='w',
                               markerfacecolor=colors_dict[g], markersize=8, label=g)
                    for g in all_groups]
    ax.legend(handles=legend_elems, loc='upper right',
              bbox_to_anchor=(1.02,1.0), fontsize=9, title='Group',
              title_fontsize=10, framealpha=0.9, edgecolor='#CCCCCC')

    ax.set_ylim(0, 230)
    plt.subplots_adjust(left=0.05, right=0.80, top=0.92, bottom=0.08)

    safe = center_label.replace(" ", "_").replace("/", "_").replace("+", "plus")
    plt.savefig(f'{safe}_PRE-R_vs_PRE-NR.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(f'{safe}_PRE-R_vs_PRE-NR.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()

# ================== 5. 主流程 ==================
if __name__ == "__main__":
    ADATA_PATH = "adata_with_band_by_sample.h5ad"
    adata, center_label, target_ct = load_and_filter_data(ADATA_PATH)

    if len(target_ct) == 0:
        print("❌ 周边细胞类型不存在，请检查 celltype 列的名称。")
    else:
        mean_cnt, dens, dists, sec_areas, exist_ct, groups = \
            calculate_cell_counts(adata, center_label, target_ct)
        export_csv(mean_cnt, dens, dists, sec_areas, exist_ct, groups, center_label)
        plot_radial(mean_cnt, dens, dists, exist_ct, groups, center_label)
        print("完成！")
