import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from scipy.ndimage import gaussian_filter1d
from scipy.stats import sem
import warnings
import re
import os
from matplotlib import font_manager

warnings.filterwarnings('ignore')

# ========================
# 1. 读取数据
# ========================
adata = ad.read_h5ad('All_filter_state_samp_new_with_CNs10-15.h5ad')
print(adata)

# 提取坐标和obs信息
coords = adata.obsm['spatial']  # (n_cells, 2)
obs = adata.obs[['celltype', 'state', 'sample']].copy()
obs['x'] = coords[:, 0]
obs['y'] = coords[:, 1]

# 过滤 state 和 celltype 缺失值，只保留 N/A/C
obs = obs.dropna(subset=['state', 'celltype'])
obs = obs[obs['state'].isin(['N', 'A', 'C'])]

# 标准化细胞类型名称（去除首尾空格）
obs['celltype'] = obs['celltype'].str.strip()

# ========================
# 细胞类型名称（请按数据中实际名称填写）
# ========================
required_types = [
    'WNT+ ISC',
    'matri-prod CAF',        # 小写 m
    'CTSL+ PLAUR+ Macro',
    'CD8+ T',
    'Treg',
    'MKI67+ T'
]

# 检查存在性
existing_types = [t for t in required_types if t in obs['celltype'].unique()]
missing_types = set(required_types) - set(existing_types)
if missing_types:
    print(f"警告：以下细胞类型不存在，将跳过相关组合：{missing_types}")

# 定义所有要分析的组合 (中心细胞类型, 周围细胞类型)
combinations = [
    ('WNT+ ISC', 'matri-prod CAF'),
    ('matri-prod CAF', 'WNT+ ISC'),
    ('CTSL+ PLAUR+ Macro', 'WNT+ ISC'),
    ('CD8+ T', 'WNT+ ISC'),
    ('Treg', 'WNT+ ISC'),
    ('MKI67+ T', 'WNT+ ISC'),
    ('matri-prod CAF', 'CTSL+ PLAUR+ Macro'),
    ('CTSL+ PLAUR+ Macro', 'matri-prod CAF')
]

# 过滤掉包含缺失细胞类型的组合
combinations = [(c, t) for (c, t) in combinations 
                if c in existing_types and t in existing_types]
if not combinations:
    raise ValueError("没有有效的组合，请检查细胞类型名称")
print(f"有效组合数量：{len(combinations)}")

# ========================
# 2. 参数设置（最大距离200 μm，间隔20 μm）
# ========================
MAX_RADIUS = 200.0      # μm
BIN_WIDTH = 20.0        # μm
N_BINS = int(MAX_RADIUS / BIN_WIDTH)   # =10
bin_edges = np.linspace(0, MAX_RADIUS, N_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
annular_area = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)

# 颜色映射
state_colors = {'N': '#4daf4a', 'A': '#377eb8', 'C': '#ff7f00'}
state_labels = {'N': 'N', 'A': 'A', 'C': 'C'}

# ========================
# 3. 辅助函数：计算单个样本的密度曲线
# ========================
def compute_density_per_sample(df_sample, center_type, target_type, 
                               max_radius, bin_edges, annular_area):
    """
    返回 dict: {state: (bin_centers, density_curve)}，density_curve长度为N_BINS
    """
    center_mask = df_sample['celltype'] == center_type
    target_mask = df_sample['celltype'] == target_type
    
    if center_mask.sum() == 0 or target_mask.sum() == 0:
        return {}
    
    centers = df_sample[center_mask][['x', 'y']].values
    center_states = df_sample.loc[center_mask, 'state'].values
    targets = df_sample[target_mask][['x', 'y']].values
    
    tree = KDTree(targets)
    state_counts = {s: np.zeros(N_BINS, dtype=np.float64) for s in ['N', 'A', 'C']}
    state_num_centers = {s: 0 for s in ['N', 'A', 'C']}
    
    for (cx, cy), state in zip(centers, center_states):
        if state not in state_counts:
            continue
        indices = tree.query_ball_point([cx, cy], r=max_radius)
        if not indices:
            continue
        dists = np.linalg.norm(targets[indices] - [cx, cy], axis=1)
        bin_idx = np.digitize(dists, bin_edges) - 1
        valid = (bin_idx >= 0) & (bin_idx < N_BINS)
        bin_idx = bin_idx[valid]
        if len(bin_idx) == 0:
            continue
        add_counts = np.bincount(bin_idx, minlength=N_BINS)[:N_BINS]
        state_counts[state] += add_counts
        state_num_centers[state] += 1
    
    result = {}
    for state in ['N', 'A', 'C']:
        if state_num_centers[state] == 0:
            continue
        total_area = state_num_centers[state] * annular_area
        density = state_counts[state] / total_area
        result[state] = (bin_centers.copy(), density)
    return result

# ========================
# 4. 对所有样本和组合计算
# ========================
sample_groups = obs.groupby('sample')
all_results = {}   # {(center, target): {state: list of density curves}}

for (center_type, target_type) in combinations:
    print(f"Processing: {center_type} → {target_type}")
    state_curves = {s: [] for s in ['N', 'A', 'C']}
    
    for sample_name, df_sample in sample_groups:
        sample_res = compute_density_per_sample(
            df_sample, center_type, target_type, MAX_RADIUS, bin_edges, annular_area
        )
        for state, (_, curve) in sample_res.items():
            state_curves[state].append(curve)
    
    # 平滑（高斯滤波，sigma适当减小因为bin变少）
    sigma_smooth = 0.5
    for state in state_curves:
        for i in range(len(state_curves[state])):
            state_curves[state][i] = gaussian_filter1d(state_curves[state][i], 
                                                       sigma=sigma_smooth, mode='reflect')
    all_results[(center_type, target_type)] = state_curves

# ========================
# 5. 绘图：每个组合单独输出 PDF 和 PNG
# ========================
# 设置字体（尝试加载 Arial，如果失败则使用 sans-serif）
try:
    font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')
    plt.rcParams['font.family'] = 'arial'
except:
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'Liberation Sans']

plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8

# 创建输出目录（可选）
out_dir = 'spatial_density_plots'
os.makedirs(out_dir, exist_ok=True)

# 辅助函数：生成安全的文件名
def safe_filename(center, target):
    name = f"{center}_to_{target}"
    # 替换可能导致问题的字符
    name = re.sub(r'[\\/*?:"<>|]', '_', name)
    name = name.replace('→', 'to').replace('+', 'plus').replace(' ', '_')
    return name

for (center_type, target_type), state_curves in all_results.items():
    # 创建新图形
    fig, ax = plt.subplots(figsize=(5, 4), dpi=300)
    all_y_vals = []
    
    for state in ['N', 'A', 'C']:
        curves = state_curves[state]
        if len(curves) == 0:
            continue
        curves_arr = np.vstack(curves)   # (n_samples, n_bins)
        mean_curve = np.mean(curves_arr, axis=0)
        sem_curve = sem(curves_arr, axis=0, nan_policy='omit')
        x = bin_centers
        ax.plot(x, mean_curve, color=state_colors[state], lw=1.8, label=state_labels[state])
        ax.fill_between(x, mean_curve - sem_curve, mean_curve + sem_curve,
                         color=state_colors[state], alpha=0.2, edgecolor='none')
        all_y_vals.extend(mean_curve)
    
    # 装饰图形
    ax.set_title(f"{center_type} → {target_type}", fontsize=10, fontweight='bold', pad=8)
    ax.set_xlabel('Distance (μm)', fontsize=8)
    ax.set_ylabel('Cell density (cells/μm²)', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=7)
    ax.set_xlim(0, MAX_RADIUS)
    if all_y_vals:
        ax.set_ylim(bottom=0, top=np.percentile(all_y_vals, 95) * 1.1)
    ax.legend(loc='upper right', frameon=False, fontsize=7, handlelength=1.5)
    
    plt.tight_layout()
    
    # 保存文件
    base_name = safe_filename(center_type, target_type)
    pdf_path = os.path.join(out_dir, f"{base_name}.pdf")
    png_path = os.path.join(out_dir, f"{base_name}.png")
    plt.savefig(pdf_path, bbox_inches='tight', dpi=300)
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    plt.close(fig)   # 关闭当前图形，释放内存
    
    print(f"Saved: {pdf_path} , {png_path}")

print("所有组合的图形已分别保存至目录:", out_dir)
