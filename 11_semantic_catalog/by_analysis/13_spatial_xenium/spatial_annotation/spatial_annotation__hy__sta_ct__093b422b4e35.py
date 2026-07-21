import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.ndimage import gaussian_filter1d
from scipy.stats import sem, mannwhitneyu
from itertools import permutations
import warnings
import os
import re
from matplotlib import font_manager

warnings.filterwarnings('ignore')

# ========================
# 1. 读取数据
# ========================
adata = ad.read_h5ad('extracted_samples.h5ad')
print(adata)

# 提取坐标和obs信息
coords = adata.obsm['spatial']  # (n_cells, 2)
obs = adata.obs[['annotation', 'mucus_zone', 'sample']].copy()
obs['x'] = coords[:, 0]
obs['y'] = coords[:, 1]

# 过滤 mucus_zone 和 annotation 缺失值，只保留 inner/outer
obs = obs.dropna(subset=['mucus_zone', 'annotation'])
obs = obs[obs['mucus_zone'].isin(['inner', 'outer'])]

# 标准化细胞类型名称（去除首尾空格）
obs['annotation'] = obs['annotation'].str.strip()

# ========================
# 细胞类型名称（按新数据指定）
# ========================
required_types = [
    'Epi-CAF',
    'Tumor_Vascular_endothelium',
    'EMT_myCAF',
    'MCAM+_CAF',
    'EMT_Epithelium'
]

# 检查存在性
existing_types = [t for t in required_types if t in obs['annotation'].unique()]
missing_types = set(required_types) - set(existing_types)
if missing_types:
    print(f"警告：以下细胞类型不存在，将跳过相关组合：{missing_types}")

# 自动生成所有有序组合（中心细胞类型, 周围细胞类型），即 A→B 且 A≠B
combinations = list(permutations(existing_types, 2))

if not combinations:
    raise ValueError("没有有效的组合，请检查细胞类型名称")
print(f"有效组合数量：{len(combinations)}")

# ========================
# 2. 新增：细胞类型组成分析（按样本、分区统计比例）
# ========================
print("\n=== 开始细胞类型组成分析 ===")
sample_groups = obs.groupby('sample')
out_dir_comp = 'composition_analysis'
os.makedirs(out_dir_comp, exist_ok=True)

# 2.1 统计每个样本、每个区域的各类型细胞数量及比例
comp_records = []
for sample_name, df_sample in sample_groups:
    for zone in ['inner', 'outer']:
        df_zone = df_sample[df_sample['mucus_zone'] == zone]
        total_cells = len(df_zone)
        if total_cells == 0:
            continue
        counts = df_zone['annotation'].value_counts()
        for cell_type in existing_types:
            count = counts.get(cell_type, 0)
            proportion = count / total_cells
            comp_records.append({
                'sample': sample_name,
                'mucus_zone': zone,
                'annotation': cell_type,
                'count': count,
                'proportion': proportion,
                'total_cells': total_cells
            })
comp_df = pd.DataFrame(comp_records)

# 2.2 导出统计表
comp_df.to_csv(os.path.join(out_dir_comp, 'composition_summary.csv'), index=False)
print(f"组成统计表已保存至：{os.path.join(out_dir_comp, 'composition_summary.csv')}")

# 2.3 为每种细胞类型绘制箱线图，对比 inner vs outer
mucus_zone_colors = {'inner': '#1f77b4', 'outer': '#d62728'}  # 沿用后续配色

for cell_type in existing_types:
    df_type = comp_df[comp_df['annotation'] == cell_type]
    inner_props = df_type[df_type['mucus_zone'] == 'inner']['proportion'].values
    outer_props = df_type[df_type['mucus_zone'] == 'outer']['proportion'].values
    
    # 统计检验（Mann-Whitney U，非参数检验）
    if len(inner_props) > 0 and len(outer_props) > 0:
        stat, p_val = mannwhitneyu(inner_props, outer_props, alternative='two-sided')
    else:
        p_val = np.nan

    fig, ax = plt.subplots(figsize=(4, 5), dpi=300)
    data = [inner_props, outer_props]
    bp = ax.boxplot(data, labels=['inner', 'outer'], patch_artist=True,
                    boxprops=dict(facecolor='lightgray'),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2),
                    flierprops=dict(marker='o', markersize=4, alpha=0.5))
    # 上色
    for patch, color in zip(bp['boxes'], [mucus_zone_colors['inner'], mucus_zone_colors['outer']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    ax.set_ylabel('Proportion of cells', fontsize=14)
    ax.set_title(f"{cell_type}", fontweight='bold', pad=10)
    ax.set_ylim(0, 1.0)
    # 添加P值
    if not np.isnan(p_val):
        p_text = f"p = {p_val:.2e}" if p_val < 0.01 else f"p = {p_val:.3f}"
        ax.text(0.5, 0.95, p_text, ha='center', transform=ax.transAxes, fontsize=12,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    else:
        ax.text(0.5, 0.95, "No data", ha='center', transform=ax.transAxes, fontsize=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    
    safe_name = cell_type.replace('+', 'plus').replace(' ', '_').replace('/', '_')
    plt.savefig(os.path.join(out_dir_comp, f"{safe_name}_composition.pdf"), bbox_inches='tight')
    plt.savefig(os.path.join(out_dir_comp, f"{safe_name}_composition.png"), bbox_inches='tight')
    plt.close(fig)

print(f"组成分析图表已保存至：{out_dir_comp}")
print("=== 组成分析完成 ===\n")

# ========================
# 3. 参数设置（最大距离250 μm，间隔25 μm）
# ========================
MAX_RADIUS = 250.0      # μm
BIN_WIDTH = 25.0        # μm
N_BINS = int(MAX_RADIUS / BIN_WIDTH)   # =10
bin_edges = np.linspace(0, MAX_RADIUS, N_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
annular_area = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)   # 单位: μm²

# 颜色映射（仅两个组）
mucus_zone_colors = {
    'inner': '#1f77b4',   # 蓝色
    'outer': '#d62728'    # 红色
}

mucus_zone_labels = {
    'inner': 'inner',
    'outer': 'outer'
}

# ========================
# 4. 辅助函数：计算单个样本的密度曲线（原始密度：cells/μm²）
# ========================
def compute_density_per_sample(df_sample, center_type, target_type,
                               max_radius, bin_edges, annular_area):
    center_mask = df_sample['annotation'] == center_type
    target_mask = df_sample['annotation'] == target_type

    if center_mask.sum() == 0 or target_mask.sum() == 0:
        return {}

    centers = df_sample[center_mask][['x', 'y']].values
    center_mucus_zones = df_sample.loc[center_mask, 'mucus_zone'].values
    targets = df_sample[target_mask][['x', 'y']].values

    tree = KDTree(targets)
    mucus_zone_counts = {s: np.zeros(N_BINS, dtype=np.float64) for s in ['inner', 'outer']}
    mucus_zone_num_centers = {s: 0 for s in ['inner', 'outer']}

    for (cx, cy), mucus_zone in zip(centers, center_mucus_zones):
        if mucus_zone not in mucus_zone_counts:
            continue
        indices = tree.query_ball_point([cx, cy], r=max_radius)
        if not indices:
            continue
        dists = np.linalg.norm(targets[indices] - [cx, cy], axis=1)
        # 排除距离为0的细胞（理论上不会出现，因为中心细胞类型与周围不同）
        dists = dists[dists > 0]
        if len(dists) == 0:
            continue
        bin_idx = np.digitize(dists, bin_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < N_BINS)
        bin_idx = bin_idx[valid]
        if len(bin_idx) == 0:
            continue
        add_counts = np.bincount(bin_idx, minlength=N_BINS)[:N_BINS]
        mucus_zone_counts[mucus_zone] += add_counts
        mucus_zone_num_centers[mucus_zone] += 1

    result = {}
    for mucus_zone in ['inner', 'outer']:
        if mucus_zone_num_centers[mucus_zone] == 0:
            continue
        total_area = mucus_zone_num_centers[mucus_zone] * annular_area   # μm²
        density = mucus_zone_counts[mucus_zone] / total_area             # cells/μm²
        result[mucus_zone] = (bin_centers.copy(), density)
    return result

# ========================
# 5. 对所有样本和组合计算
# ========================
sample_groups = obs.groupby('sample')
all_results = {}

for (center_type, target_type) in combinations:
    print(f"Processing: {center_type} → {target_type}")
    mucus_zone_curves = {s: [] for s in ['inner', 'outer']}

    for sample_name, df_sample in sample_groups:
        sample_res = compute_density_per_sample(
            df_sample, center_type, target_type, MAX_RADIUS, bin_edges, annular_area
        )
        for mucus_zone, (_, curve) in sample_res.items():
            mucus_zone_curves[mucus_zone].append(curve)

    # 平滑（高斯滤波）
    sigma_smooth = 0.5
    for mucus_zone in mucus_zone_curves:
        for i in range(len(mucus_zone_curves[mucus_zone])):
            mucus_zone_curves[mucus_zone][i] = gaussian_filter1d(mucus_zone_curves[mucus_zone][i],
                                                       sigma=sigma_smooth, mode='reflect')
    all_results[(center_type, target_type)] = mucus_zone_curves

# ========================
# 6. 绘图：每个组合单独生成一张图，图例放在右侧外部
# ========================
# 创建输出目录
out_dir = 'distance'
os.makedirs(out_dir, exist_ok=True)

# 辅助函数：生成安全的文件名
def safe_filename(center, target):
    name = f"{center}_to_{target}"
    name = re.sub(r'[\\/*?:"<>|]', '_', name)
    name = name.replace('→', 'to').replace('+', 'plus').replace(' ', '_')
    return name

# 设置字体（尝试加载 Arial，失败则用 sans-serif）
try:
    font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')
    plt.rcParams['font.family'] = 'arial'
except:
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']

# 字体大小设置
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 17
plt.rcParams['axes.titlesize'] = 17
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

for (center_type, target_type), mucus_zone_curves in all_results.items():
    # 创建新图，尺寸 5×4 英寸
    fig, ax = plt.subplots(figsize=(5, 4), dpi=300)
    all_y_vals = []

    for mucus_zone in ['inner', 'outer']:
        curves = mucus_zone_curves[mucus_zone]
        if len(curves) == 0:
            continue
        curves_arr = np.vstack(curves)                    # cells/μm²
        mean_curve = np.mean(curves_arr, axis=0) * 1e6    # cells/mm²
        sem_curve = sem(curves_arr, axis=0, nan_policy='omit') * 1e6
        x = bin_centers
        ax.plot(x, mean_curve, color=mucus_zone_colors[mucus_zone], lw=2.0, label=mucus_zone_labels[mucus_zone])
        ax.fill_between(x, mean_curve - sem_curve, mean_curve + sem_curve,
                         color=mucus_zone_colors[mucus_zone], alpha=0.2, edgecolor='none')
        all_y_vals.extend(mean_curve)

    # 装饰图形
    ax.set_title(f"{center_type} → {target_type}", fontweight='bold', pad=10)
    ax.set_xlabel('Distance (μm)')
    ax.set_ylabel('Cell density (cells/mm²)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, MAX_RADIUS)
    if all_y_vals:
        ax.set_ylim(bottom=0, top=np.percentile(all_y_vals, 95) * 1.1)
    
    # 将图例放到图形右侧外部
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, handlelength=1.5)

    plt.tight_layout()

    # 保存文件（bbox_inches='tight' 会自动包含轴外图例）
    base_name = safe_filename(center_type, target_type)
    pdf_path = os.path.join(out_dir, f"{base_name}.pdf")
    png_path = os.path.join(out_dir, f"{base_name}.png")
    plt.savefig(pdf_path, bbox_inches='tight', dpi=300)
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

    print(f"Saved: {pdf_path} , {png_path}")

print(f"所有组合的图形已分别保存至目录: {out_dir}")
