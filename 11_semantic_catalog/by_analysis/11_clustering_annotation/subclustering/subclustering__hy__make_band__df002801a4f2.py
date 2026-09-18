import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree, ConvexHull
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
import os
import warnings
warnings.filterwarnings('ignore')

# ==================== 工具函数 ====================
def create_output_folder(folder_name='make_band'):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

# ==================== 边界计算（含线条切断） ====================
def calculate_boundary_for_pair(adata, sample, sample_adata, coords, labels,
                                cat_a, cat_b, cat_focus,
                                n_neighbors=10, smooth_sigma=3,
                                min_vertex_count=20):
    # 1. 识别边界细胞
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    boundary_mask = np.zeros(len(sample_adata), dtype=bool)
    for i in range(len(sample_adata)):
        current_label = labels[i]
        if current_label not in [cat_a, cat_b]:
            continue
        neighbor_labels = labels[indices[i]]
        if current_label == cat_a and cat_b in neighbor_labels:
            boundary_mask[i] = True
        elif current_label == cat_b and cat_a in neighbor_labels:
            boundary_mask[i] = True

    # 2. 只使用参与边界的两类细胞训练分类器
    sub_mask = (labels == cat_a) | (labels == cat_b)
    if np.sum(sub_mask) < 10:
        return None, None, None, "none"
    
    coords_sub = coords[sub_mask]
    labels_sub = labels[sub_mask]
    y_sub = (labels_sub == cat_b).astype(int) # 以 cat_b（即A）为正类
    clf = KNeighborsClassifier(n_neighbors=n_neighbors)
    clf.fit(coords_sub, y_sub)

    # 3. 生成网格并提取轮廓线
    padding_percent = 0.1
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    x_pad = (x_max - x_min) * padding_percent
    y_pad = (y_max - y_min) * padding_percent
    grid_resolution = 300
    xx, yy = np.meshgrid(np.linspace(x_min - x_pad, x_max + x_pad, grid_resolution),
                         np.linspace(y_min - y_pad, y_max + y_pad, grid_resolution))
    grid_points = np.c_[xx.ravel(), yy.ravel()]

    nbrs_geometry = NearestNeighbors(n_neighbors=1).fit(coords)
    dists, _ = nbrs_geometry.kneighbors(grid_points)
    dists = dists.reshape(xx.shape)

    Z = clf.predict_proba(grid_points)[:, 1]
    Z = Z.reshape(xx.shape)
    Z_smooth = gaussian_filter(Z, sigma=smooth_sigma)

    # ==================== 切断多余的延伸线 ====================
    # 找出每个网格点距离最近的真实细胞，判断其属于哪一类
    nbrs_all = NearestNeighbors(n_neighbors=1).fit(coords)
    dists_all, idxs_all = nbrs_all.kneighbors(grid_points)
    nearest_labels = labels[idxs_all.flatten()]

    # 创建掩膜：只有网格点最近的真实细胞是 cat_a 或 cat_b 的地方，才允许画边界线
    valid_region_mask = (nearest_labels == cat_a) | (nearest_labels == cat_b)
    # 将无效区域（如 N-A 穿过 C 区域的点）设为 np.nan，使轮廓线断开
    Z_smooth[~valid_region_mask.reshape(xx.shape)] = np.nan
    # ==============================================================

    valid_boundary_coords = []
    method_used = "none"

    try:
        contours = plt.contour(xx, yy, Z_smooth, levels=[0.5], alpha=0)
        plt.close()
        if len(contours.allsegs) > 0 and len(contours.allsegs[0]) > 0:
            for seg in contours.allsegs[0]:
                if len(seg) >= min_vertex_count:
                    valid_boundary_coords.append(seg)
            if valid_boundary_coords:
                method_used = "contour"
    except Exception:
        pass

    # 后备方案：凸包 / 包围盒
    if not valid_boundary_coords:
        try:
            focus_coords = coords[labels == cat_focus]
            if len(focus_coords) > 5:
                hull = ConvexHull(focus_coords)
                hull_points = focus_coords[hull.vertices]
                hull_closed = np.vstack([hull_points, hull_points[0]])
                valid_boundary_coords.append(hull_closed)
                method_used = "convex_hull"
        except Exception:
            pass

    if not valid_boundary_coords:
        try:
            focus_coords = coords[labels == cat_focus]
            if len(focus_coords) > 0:
                x_min_inner, x_max_inner = focus_coords[:, 0].min(), focus_coords[:, 0].max()
                y_min_inner, y_max_inner = focus_coords[:, 1].min(), focus_coords[:, 1].max()
                padding = 50
                box_points = np.array([
                    [x_min_inner - padding, y_min_inner - padding],
                    [x_max_inner + padding, y_min_inner - padding],
                    [x_max_inner + padding, y_max_inner + padding],
                    [x_min_inner - padding, y_max_inner + padding],
                    [x_min_inner - padding, y_min_inner - padding]
                ])
                valid_boundary_coords.append(box_points)
                method_used = "bounding_box"
        except Exception:
            pass

    if len(valid_boundary_coords) == 0:
        return None, None, None, method_used

    all_boundary_points = np.vstack(valid_boundary_coords)
    all_boundary_points = np.unique(all_boundary_points, axis=0)
    tree = cKDTree(all_boundary_points)
    dist_abs, _ = tree.query(coords, k=1)

    signed_dist = dist_abs.copy()
    other_mask = (labels != cat_focus)
    signed_dist[other_mask] = -signed_dist[other_mask]

    return signed_dist, valid_boundary_coords, boundary_mask, method_used

def calculate_dual_boundaries_by_sample(adata, obs_col='state',
                                        n_neighbors=10, smooth_sigma=3,
                                        min_vertex_count=20,
                                        save_dir='make_band'):
    save_dir = create_output_folder(save_dir)

    if obs_col not in adata.obs.columns:
        raise ValueError(f"adata.obs 中未找到 '{obs_col}' 列")
    required_states = {'N', 'A', 'C'}
    if not required_states.issubset(adata.obs[obs_col].unique()):
        raise ValueError(f"{obs_col} 列必须包含 N, A, C 三种状态")

    sample_col = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col = col
            break
    if sample_col is None:
        adata.obs['sample'] = 'sample_1'
        sample_col = 'sample'

    adata.obs['dist_NA'] = np.nan
    adata.obs['dist_AC'] = np.nan
    adata.obs['boundary_method_NA'] = 'none'
    adata.obs['boundary_method_AC'] = 'none'

    samples = adata.obs[sample_col].unique()
    print(f"发现 {len(samples)} 个样本，开始计算双边界...")

    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask].copy()
        coords = sample_adata.obsm['spatial']
        labels = sample_adata.obs[obs_col].values

        print(f"\n处理样本: {sample} ({len(sample_adata)} 个细胞)")

        states_present = set(labels)
        if not {'N','A'}.issubset(states_present) and not {'A','C'}.issubset(states_present):
            print(f"  警告: 样本 {sample} 缺少必要的状态组合，跳过")
            continue

        segs_na_val = None
        segs_ac_val = None
        method_na_val = 'none'
        method_ac_val = 'none'

        if 'N' in states_present and 'A' in states_present:
            dist_na, segs_na_val, _, method_na_val = calculate_boundary_for_pair(
                adata, sample, sample_adata, coords, labels,
                cat_a='N', cat_b='A', cat_focus='A',
                n_neighbors=n_neighbors, smooth_sigma=smooth_sigma,
                min_vertex_count=min_vertex_count
            )
            if dist_na is not None:
                adata.obs.loc[sample_mask, 'dist_NA'] = dist_na
                adata.obs.loc[sample_mask, 'boundary_method_NA'] = method_na_val
                print(f"    N-A 边界计算方法: {method_na_val}, 距离范围 [{dist_na.min():.1f}, {dist_na.max():.1f}]")

        if 'A' in states_present and 'C' in states_present:
            dist_ac, segs_ac_val, _, method_ac_val = calculate_boundary_for_pair(
                adata, sample, sample_adata, coords, labels,
                cat_a='C', cat_b='A', cat_focus='A',
                n_neighbors=n_neighbors, smooth_sigma=smooth_sigma,
                min_vertex_count=min_vertex_count
            )
            if dist_ac is not None:
                adata.obs.loc[sample_mask, 'dist_AC'] = dist_ac
                adata.obs.loc[sample_mask, 'boundary_method_AC'] = method_ac_val
                print(f"    C-A 边界计算方法: {method_ac_val}, 距离范围 [{dist_ac.min():.1f}, {dist_ac.max():.1f}]")

        # 可视化边界
        fig, ax = plt.subplots(figsize=(8, 8))
        colors_state = {'N': 'blue', 'A': 'yellow', 'C': 'red'}
        for st, col in colors_state.items():
            mask = labels == st
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], c=col, s=5, alpha=0.6, label=st)

        if segs_na_val:
            for i, seg in enumerate(segs_na_val):
                ax.plot(seg[:, 0], seg[:, 1], 'black', linewidth=2,
                        label='N-A boundary' if i == 0 else "")
        if segs_ac_val:
            for i, seg in enumerate(segs_ac_val):
                ax.plot(seg[:, 0], seg[:, 1], 'green', linewidth=2, linestyle='--',
                        label='C-A boundary' if i == 0 else "")

        ax.set_title(f'Sample: {sample}\nN-A method: {method_na_val}, C-A method: {method_ac_val}')
        ax.legend()
        plt.tight_layout()

        fig_path_png = os.path.join(save_dir, f'sample_{sample}_boundaries.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_boundaries.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        plt.close(fig)
        print(f"  边界图已保存: {fig_path_png} 和 {fig_path_pdf}")

    return adata

# ==================== 带状划分与 final 列合并 ====================
def create_bands_by_sample(adata, obs_col='state',
                           width_na=200, width_ca=200,
                           save_dir='make_band'):
    if 'dist_NA' not in adata.obs.columns or 'dist_AC' not in adata.obs.columns:
        raise ValueError("请先运行 calculate_dual_boundaries_by_sample 计算两个边界的距离")

    save_dir = create_output_folder(save_dir)

    adata.obs['AN'] = adata.obs[obs_col].astype(str)
    adata.obs['CA'] = adata.obs[obs_col].astype(str)

    states = adata.obs[obs_col].values
    dist_na = adata.obs['dist_NA'].values
    dist_ac = adata.obs['dist_AC'].values

    # AN 列 (0~200um 范围)
    cond_an1 = (states == 'N') & (dist_na >= -width_na) & (dist_na < 0)
    cond_an2 = (states == 'A') & (dist_na > 0) & (dist_na <= width_na)
    adata.obs.loc[cond_an1, 'AN'] = 'AN1'
    adata.obs.loc[cond_an2, 'AN'] = 'AN2'

    # CA 列 (0~200um 范围)
    cond_ca1 = (states == 'A') & (dist_ac > 0) & (dist_ac <= width_ca)
    cond_ca2 = (states == 'C') & (dist_ac >= -width_ca) & (dist_ac < 0)
    adata.obs.loc[cond_ca1, 'CA'] = 'CA1'
    adata.obs.loc[cond_ca2, 'CA'] = 'CA2'

    # ==================== 新增：整合 band_final 列 ====================
    adata.obs['band_final'] = adata.obs[obs_col].astype(str)
    
    an_mask = adata.obs['AN'].isin(['AN1', 'AN2'])
    ca_mask = adata.obs['CA'].isin(['CA1', 'CA2'])
    
    # 1. 重叠区域
    overlap_mask = an_mask & ca_mask
    adata.obs.loc[overlap_mask, 'band_final'] = 'overlap'
    
    # 2. 纯 AN 区域（仅属于 AN，不属于 CA）
    pure_an_mask = an_mask & (~ca_mask)
    adata.obs.loc[pure_an_mask, 'band_final'] = adata.obs.loc[pure_an_mask, 'AN']
    
    # 3. 纯 CA 区域（仅属于 CA，不属于 AN）
    pure_ca_mask = ca_mask & (~an_mask)
    adata.obs.loc[pure_ca_mask, 'band_final'] = adata.obs.loc[pure_ca_mask, 'CA']
    # ================================================================

    print("\n========== AN 列分布 ==========")
    print(adata.obs['AN'].value_counts())
    print("\n========== CA 列分布 ==========")
    print(adata.obs['CA'].value_counts())
    print("\n========== band_final 列分布 ==========")
    print(adata.obs['band_final'].value_counts())
    print(f"\n重叠区域细胞数 (overlap): {overlap_mask.sum()}")
    
    sample_col = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col = col
            break
    if sample_col is None:
        adata.obs['sample'] = 'sample_1'
        sample_col = 'sample'

    samples = adata.obs[sample_col].unique()
    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask]
        coords = sample_adata.obsm['spatial']

        # ==================== 绘图调整：精简为 2 个子图 ====================
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))

        # 图1：原始 State（用于对照位置）
        ax0 = axes[0]
        colors_state = {'N': 'blue', 'A': 'yellow', 'C': 'red'}
        for st, col in colors_state.items():
            mask = sample_adata.obs[obs_col] == st
            if mask.any():
                ax0.scatter(coords[mask, 0], coords[mask, 1], c=col, s=5, alpha=0.6, label=st)
        ax0.set_title(f'Sample {sample} - State')
        ax0.legend()

        # 图2：Band Final（已过滤 Overlap 区域）
        ax1 = axes[1]
        # 过滤掉 overlap 标记的细胞
        clean_adata = sample_adata[~(sample_adata.obs['band_final'] == 'overlap')]
        clean_coords = clean_adata.obsm['spatial']
        
        colors_final = {
            'N': 'blue', 'A': 'yellow', 'C': 'red',
            'AN1': '#00BFFF', 'AN2': '#FFD700',
            'CA1': '#FF8C00', 'CA2': '#DC143C'
        }
        for label in clean_adata.obs['band_final'].unique():
            mask = clean_adata.obs['band_final'] == label
            if mask.any():
                ax1.scatter(clean_coords[mask, 0], clean_coords[mask, 1],
                            c=colors_final.get(label, 'gray'), s=5, alpha=0.6, label=label)
        ax1.set_title(f'Sample {sample} - Final Bands (overlap removed)')
        ax1.legend()

        plt.suptitle(f'Sample {sample} Band Analysis', fontsize=16)
        plt.tight_layout()

        fig_path_png = os.path.join(save_dir, f'sample_{sample}_bands.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_bands.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        plt.close(fig)
        print(f"  样本 {sample} band 图已保存: {fig_path_png} 和 {fig_path_pdf}")

    return adata

# ==================== 综合可视化 ====================
def create_comprehensive_visualization(adata, obs_col='state', save_dir='make_band',
                                       width_na=200, width_ca=200):
    save_dir = create_output_folder(save_dir)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    ax1 = axes[0, 0]
    if 'band_final' in adata.obs.columns:
        final_counts = adata.obs['band_final'].value_counts()
        colors_final = {
            'N': 'blue', 'A': 'yellow', 'C': 'red',
            'AN1': '#00BFFF', 'AN2': '#FFD700',
            'CA1': '#FF8C00', 'CA2': '#DC143C', 'overlap': 'gray'
        }
        bar_colors = [colors_final.get(b, 'gray') for b in final_counts.index]
        ax1.bar(final_counts.index.astype(str), final_counts.values, color=bar_colors)
        ax1.set_title('Final Band Distribution')
        ax1.set_ylabel('Cell Count')
        ax1.tick_params(axis='x', rotation=45)

    ax2 = axes[0, 1]
    if 'AN' in adata.obs.columns:
        an_counts = adata.obs['AN'].value_counts()
        colors_an = {'AN1': '#00BFFF', 'AN2': '#FFD700', 'N': 'lightblue', 'A': 'lightyellow', 'C': 'lightcoral'}
        bar_colors = [colors_an.get(b, 'gray') for b in an_counts.index]
        ax2.bar(an_counts.index.astype(str), an_counts.values, color=bar_colors)
        ax2.set_title('AN Column Distribution')
        ax2.set_ylabel('Cell Count')
        ax2.tick_params(axis='x', rotation=45)

    ax3 = axes[1, 0]
    if 'dist_NA' in adata.obs.columns:
        dist_na = adata.obs['dist_NA'].dropna()
        if len(dist_na) > 0:
            ax3.hist(dist_na, bins=50, alpha=0.5, label='dist_NA', color='blue')
    if 'dist_AC' in adata.obs.columns:
        dist_ac = adata.obs['dist_AC'].dropna()
        if len(dist_ac) > 0:
            ax3.hist(dist_ac, bins=50, alpha=0.5, label='dist_AC', color='green')
    ax3.axvline(0, color='black', linestyle='--')
    ax3.set_title('Distance to Boundaries')
    ax3.legend()

    ax4 = axes[1, 1]
    methods = []
    for mcol in ['boundary_method_NA', 'boundary_method_AC']:
        if mcol in adata.obs.columns:
            methods.append(adata.obs[mcol])
    if methods:
        method_combined = pd.concat(methods).value_counts()
        ax4.bar(method_combined.index.astype(str), method_combined.values)
        ax4.set_title('Boundary Methods Used')
        ax4.tick_params(axis='x', rotation=45)

    plt.suptitle('Comprehensive Dual-Boundary Analysis (AN / CA)', fontsize=16)
    plt.tight_layout()

    fig_path_png = os.path.join(save_dir, 'comprehensive_analysis.png')
    fig_path_pdf = os.path.join(save_dir, 'comprehensive_analysis.pdf')
    fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f"综合图已保存: {fig_path_png} 和 {fig_path_pdf}")

# ==================== 保存与主流程 ====================
def save_adata_safely(adata, filename, save_dir='make_band'):
    save_dir = create_output_folder(save_dir)
    filepath = os.path.join(save_dir, filename)
    adata_save = adata.copy()

    keys_to_remove = []
    for key in adata_save.uns.keys():
        try:
            val = adata_save.uns[key]
            if isinstance(val, list) and len(val) > 0:
                if all(isinstance(x, np.ndarray) for x in val):
                    shapes = [x.shape for x in val]
                    if len(set(shapes)) > 1:
                        keys_to_remove.append(key)
        except:
            keys_to_remove.append(key)

    for key in keys_to_remove:
        print(f"移除 uns 键: {key}")
        del adata_save.uns[key]

    adata_save.write(filepath)
    print(f"数据已保存: {filepath}")

    summary_cols = ['state', 'AN', 'CA', 'band_final', 'dist_NA', 'dist_AC',
                    'boundary_method_NA', 'boundary_method_AC']
    if 'sample' in adata.obs.columns:
        summary_cols.append('sample')
    existing = [c for c in summary_cols if c in adata.obs.columns]
    adata.obs[existing].to_csv(os.path.join(save_dir, 'band_summary.csv'), index=False)
    print("CSV 摘要已保存")

def process_regions_by_sample(adata, sample_col='sample', save_dir='make_band',
                             width_na=200, width_ca=200):
    save_dir = create_output_folder(save_dir)

    print("步骤1: 计算 N-A 和 C-A 双边界...")
    adata = calculate_dual_boundaries_by_sample(
        adata, obs_col='state', n_neighbors=10, smooth_sigma=3,
        min_vertex_count=20, save_dir=save_dir
    )

    print("\n步骤2: 根据边界向两侧扩展创建 AN 和 CA 列...")
    adata = create_bands_by_sample(
        adata, obs_col='state',
        width_na=width_na, width_ca=width_ca,
        save_dir=save_dir
    )

    print("\n步骤3: 创建综合可视化...")
    create_comprehensive_visualization(
        adata, obs_col='state', save_dir=save_dir,
        width_na=width_na, width_ca=width_ca
    )

    print("\n步骤4: 保存数据...")
    save_adata_safely(adata, 'adata_with_bands.h5ad', save_dir=save_dir)

    print("\n完成！最终 obs 列预览：")
    print(adata.obs[['state', 'AN', 'CA', 'band_final', 'dist_NA', 'dist_AC']].head(10))

    return adata

def main():
    adata = sc.read_h5ad('data_4_band.h5ad')

    if 'state' not in adata.obs.columns:
        raise ValueError("未找到 'state' 列")
    if 'spatial' not in adata.obsm:
        raise ValueError("未找到空间坐标 'spatial'")

    sample_col = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col = col
            break
    if sample_col is None:
        adata.obs['sample'] = 'sample_1'
        sample_col = 'sample'

    adata_result = process_regions_by_sample(
        adata,
        sample_col=sample_col,
        width_na=200,
        width_ca=200,
        save_dir='make_band'
    )
    return adata_result

if __name__ == "__main__":
    adata_result = main()
