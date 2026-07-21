import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.ndimage import gaussian_filter1d
from scipy.stats import sem
from itertools import permutations
import warnings
import os
import re
from matplotlib import font_manager

warnings.filterwarnings('ignore')

# ========================
# 1. 读取数据
# ========================
adata = ad.read_h5ad('All_filter.h5ad')
print(adata)

# 提取坐标和obs信息（必须包含 'sample' 列）
coords = adata.obsm['spatial']
obs = adata.obs[['celltype', 'group_final', 'sample']].copy()
obs['x'] = coords[:, 0]
obs['y'] = coords[:, 1]

obs = obs.dropna(subset=['celltype', 'group_final', 'sample'])
obs['celltype'] = obs['celltype'].str.strip()

# ========================
# 细胞类型名称
# ========================
required_types = [
    'Matrix Fibroblasts',
    'IFN-Epithelials',
    'Wnt-active stem',
    'CD8+ T'
]

existing_types = [t for t in required_types if t in obs['celltype'].unique()]
missing_types = set(required_types) - set(existing_types)
if missing_types:
    print(f"警告：以下细胞类型不存在，将跳过相关组合：{missing_types}")

combinations = list(permutations(existing_types, 2))
if not combinations:
    raise ValueError("没有有效的组合，请检查细胞类型名称")
print(f"有效组合数量：{len(combinations)}")

# ========================
# 2. 参数设置
# ========================
MAX_RADIUS = 250.0
BIN_WIDTH = 25.0
N_BINS = int(MAX_RADIUS / BIN_WIDTH)
bin_edges = np.linspace(0, MAX_RADIUS, N_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
annular_area = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)

group_colors = {
    'POST-NR': '#1f77b4',
    'POST-R': '#ff7f0e',
    'PRE-NR': '#2ca02c',
    'PRE-R': '#d62728'
}
groups = ['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R']

# ========================
# 3. 辅助函数
# ========================
def compute_density_per_sample(df_sample, center_type, target_type,
                               max_radius, bin_edges, annular_area):
    center_mask = df_sample['celltype'] == center_type
    target_mask = df_sample['celltype'] == target_type

    if center_mask.sum() == 0 or target_mask.sum() == 0:
        return None

    centers = df_sample[center_mask][['x', 'y']].values
    targets = df_sample[target_mask][['x', 'y']].values

    tree = KDTree(targets)
    total_counts = np.zeros(N_BINS, dtype=np.float64)
    total_num_centers = 0

    for (cx, cy) in centers:
        indices = tree.query_ball_point([cx, cy], r=max_radius)
        if not indices:
            continue
        dists = np.linalg.norm(targets[indices] - [cx, cy], axis=1)
        dists = dists[dists > 0]
        if len(dists) == 0:
            continue
        bin_idx = np.digitize(dists, bin_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < N_BINS)
        bin_idx = bin_idx[valid]
        if len(bin_idx) == 0:
            continue
        add_counts = np.bincount(bin_idx, minlength=N_BINS)[:N_BINS]
        total_counts += add_counts
        total_num_centers += 1

    if total_num_centers == 0:
        return None
    total_area = total_num_centers * annular_area
    density = total_counts / total_area
    return density

# ========================
# 4. 主计算
# ========================
sample_group_map = obs.groupby('sample')['group_final'].first().to_dict()
sample_groups = obs.groupby('sample')

all_results = {}

for (center_type, target_type) in combinations:
    print(f"Processing: {center_type} → {target_type}")
    group_curves = {g: [] for g in groups}

    for sample_name, df_sample in sample_groups:
        if sample_name not in sample_group_map:
            continue
        cur_group = sample_group_map[sample_name]
        if cur_group not in groups:
            continue

        density_curve = compute_density_per_sample(
            df_sample, center_type, target_type, MAX_RADIUS, bin_edges, annular_area
        )
        if density_curve is not None:
            group_curves[cur_group].append(density_curve)

    # 平滑
    sigma_smooth = 0.5
    for g in groups:
        for i in range(len(group_curves[g])):
            group_curves[g][i] = gaussian_filter1d(
                group_curves[g][i], sigma=sigma_smooth, mode='reflect'
            )

    all_results[(center_type, target_type)] = group_curves

# ========================
# 5. 绘图（图例置于图外右侧）
# ========================
out_dir = 'spatial_density_plots_by_group'
os.makedirs(out_dir, exist_ok=True)

def safe_filename(center, target):
    name = f"{center}_to_{target}"
    name = re.sub(r'[\\/*?:"<>|]', '_', name)
    name = name.replace('→', 'to').replace('+', 'plus').replace(' ', '_')
    return name

# 字体设置
try:
    font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')
    plt.rcParams['font.family'] = 'arial'
except:
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']

plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 17
plt.rcParams['axes.titlesize'] = 17
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

for (center_type, target_type), group_curves in all_results.items():
    # 增加图形宽度，为图例留出空间（6×4英寸）
    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    all_y_vals = []

    for group in groups:
        curves = group_curves[group]
        if len(curves) == 0:
            continue
        curves_arr = np.vstack(curves)
        mean_curve = np.mean(curves_arr, axis=0) * 1e6
        sem_curve = sem(curves_arr, axis=0, nan_policy='omit') * 1e6
        x = bin_centers
        ax.plot(x, mean_curve, color=group_colors[group], lw=2.0, label=group)
        ax.fill_between(x, mean_curve - sem_curve, mean_curve + sem_curve,
                         color=group_colors[group], alpha=0.2, edgecolor='none')
        all_y_vals.extend(mean_curve)

    ax.set_title(f"{center_type} → {target_type}", fontweight='bold', pad=10)
    ax.set_xlabel('Distance (μm)')
    ax.set_ylabel('Cell density (cells/mm²)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, MAX_RADIUS)
    if all_y_vals:
        ax.set_ylim(bottom=0, top=np.percentile(all_y_vals, 95) * 1.1)

    # === 关键修改：图例移到图的右侧外部 ===
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),
              frameon=False, handlelength=1.5)

    plt.tight_layout()

    base_name = safe_filename(center_type, target_type)
    pdf_path = os.path.join(out_dir, f"{base_name}.pdf")
    png_path = os.path.join(out_dir, f"{base_name}.png")
    plt.savefig(pdf_path, bbox_inches='tight', dpi=300)
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

    print(f"Saved: {pdf_path} , {png_path}")

print(f"所有组合的图形已分别保存至目录: {out_dir}")
