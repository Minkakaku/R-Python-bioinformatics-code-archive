import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree, ConvexHull
from scipy.ndimage import gaussian_filter
from scipy.stats import mannwhitneyu
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
import os
import warnings
import pickle
warnings.filterwarnings('ignore')

# ====================== 基础功能（边界计算与条带划分） ======================
def create_output_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

def calculate_boundary_by_sample(adata, obs_col, cat_a, cat_b, sample_col='sample',
                                 n_neighbors=10, cat_focus='inner',
                                 smooth_sigma=3, min_vertex_count=20,
                                 max_dist_threshold=500, save_dir='make_band'):
    save_dir = create_output_folder(save_dir)
    if sample_col not in adata.obs.columns:
        raise ValueError(f"adata.obs中未找到'{sample_col}'列")
    adata.obs['is_boundary'] = False
    adata.obs['boundary_distance'] = np.nan
    adata.obs['boundary_method'] = 'none'
    samples = adata.obs[sample_col].unique()
    print(f"发现 {len(samples)} 个样本")
    sample_boundaries = {}
    processed_samples = []
    skipped_samples = []
    for idx, sample in enumerate(samples):
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask].copy()
        if len(sample_adata) == 0:
            continue
        print(f"\n处理样本: {sample} ({len(sample_adata)} 个细胞)")
        zones_in_sample = sample_adata.obs[obs_col].unique()
        if cat_a not in zones_in_sample or cat_b not in zones_in_sample:
            print(f"  警告: 样本 {sample} 不包含 {cat_a} 和 {cat_b} 两种区域，跳过")
            skipped_samples.append(sample)
            continue
        processed_samples.append(sample)
        coords = sample_adata.obsm['spatial']
        labels = sample_adata.obs[obs_col].values
        nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        boundary_mask = np.zeros(len(sample_adata), dtype=bool)
        for i in range(len(sample_adata)):
            current_label = labels[i]
            if current_label not in [cat_a, cat_b]:
                continue
            neighbor_labels = labels[indices[i]]
            if (current_label == cat_a and cat_b in neighbor_labels) or \
               (current_label == cat_b and cat_a in neighbor_labels):
                boundary_mask[i] = True
        print(f"  识别出 {sum(boundary_mask)} 个边界细胞")
        X = coords
        y = (sample_adata.obs[obs_col] == cat_focus).astype(int).values
        clf = KNeighborsClassifier(n_neighbors=n_neighbors)
        clf.fit(X, y)
        padding_percent = 0.1
        x_min, x_max = X[:, 0].min(), X[:, 0].max()
        y_min, y_max = X[:, 1].min(), X[:, 1].max()
        x_pad = (x_max - x_min) * padding_percent
        y_pad = (y_max - y_min) * padding_percent
        grid_resolution = 300
        xx, yy = np.meshgrid(np.linspace(x_min - x_pad, x_max + x_pad, grid_resolution),
                             np.linspace(y_min - y_pad, y_max + y_pad, grid_resolution))
        grid_points = np.c_[xx.ravel(), yy.ravel()]
        nbrs_geometry = NearestNeighbors(n_neighbors=1).fit(X)
        dists, _ = nbrs_geometry.kneighbors(grid_points)
        dists = dists.reshape(xx.shape)
        Z = clf.predict_proba(grid_points)[:, 1]
        Z = Z.reshape(xx.shape)
        Z_smooth = gaussian_filter(Z, sigma=smooth_sigma)
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
                    print(f"  使用标准等高线法找到边界")
        except Exception as e:
            print(f"  等高线法失败: {e}")
        if not valid_boundary_coords:
            print(f"  标准等高线法未找到边界，尝试凸包法...")
            try:
                inner_coords = coords[sample_adata.obs[obs_col] == cat_a]
                if len(inner_coords) > 5:
                    hull = ConvexHull(inner_coords)
                    hull_points = inner_coords[hull.vertices]
                    hull_closed = np.vstack([hull_points, hull_points[0]])
                    valid_boundary_coords.append(hull_closed)
                    method_used = "convex_hull"
                    print(f"  使用凸包法找到边界")
            except Exception as e:
                print(f"  凸包法失败: {e}")
        if not valid_boundary_coords:
            print(f"  凸包法也失败，尝试边界框法...")
            try:
                inner_coords = coords[sample_adata.obs[obs_col] == cat_a]
                if len(inner_coords) > 0:
                    x_min_inner, x_max_inner = inner_coords[:, 0].min(), inner_coords[:, 0].max()
                    y_min_inner, y_max_inner = inner_coords[:, 1].min(), inner_coords[:, 1].max()
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
                    print(f"  使用边界框法找到边界")
            except Exception as e:
                print(f"  边界框法失败: {e}")
        if valid_boundary_coords:
            all_boundary_points = np.vstack(valid_boundary_coords)
            all_boundary_points = np.unique(all_boundary_points, axis=0)
            tree = cKDTree(all_boundary_points)
            dist_values, _ = tree.query(X, k=1)
            signed_distances = dist_values.copy()
            outer_mask = (sample_adata.obs[obs_col] == cat_b).values
            signed_distances[outer_mask] = -signed_distances[outer_mask]
            sample_indices = adata.obs[sample_col] == sample
            adata.obs.loc[sample_indices, 'boundary_distance'] = signed_distances
            adata.obs.loc[sample_indices, 'is_boundary'] = boundary_mask
            adata.obs.loc[sample_indices, 'boundary_method'] = method_used
            sample_boundaries[sample] = {'coords': valid_boundary_coords, 'method': method_used}
            print(f"  距离计算完成：范围 [{signed_distances.min():.2f}, {signed_distances.max():.2f}]，方法: {method_used}")
        else:
            sample_indices = adata.obs[sample_col] == sample
            adata.obs.loc[sample_indices, 'is_boundary'] = boundary_mask
            adata.obs.loc[sample_indices, 'boundary_method'] = 'none'
        # 绘制样本边界图
        fig, ax = plt.subplots(figsize=(8, 8))
        for zone, color in [(cat_a, 'blue'), (cat_b, 'red')]:
            mask = sample_adata.obs[obs_col] == zone
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], c=color, s=5, alpha=0.6, label=zone)
        for seg in valid_boundary_coords:
            ax.plot(seg[:, 0], seg[:, 1], color='black', linewidth=2, label=f'Boundary ({method_used})')
        ax.set_title(f'Sample: {sample} - Boundary Fitting\nMethod: {method_used}')
        ax.legend()
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_boundary.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_boundary.pdf'), bbox_inches='tight')
        plt.close()
    with open(os.path.join(save_dir, 'sample_boundaries.pkl'), 'wb') as f:
        pickle.dump(sample_boundaries, f)
    return adata

def create_region_bands_by_sample(adata, sample_col='sample',
                                  on_band_width=100, out_band_width=200, in_band_width=200,
                                  save_dir='make_band'):
    if 'boundary_distance' not in adata.obs:
        raise ValueError("请先运行 calculate_boundary_by_sample")
    save_dir = create_output_folder(save_dir)
    adata.obs['site_band'] = np.nan
    samples = adata.obs[sample_col].unique()
    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        if not sample_mask.any():
            continue
        sample_distances = adata.obs.loc[sample_mask, 'boundary_distance'].values
        sample_zones = adata.obs.loc[sample_mask, 'region'].values
        if np.isnan(sample_distances).all():
            continue
        site_band_list = []
        for i, dist in enumerate(sample_distances):
            zone = sample_zones[i]
            if pd.isna(dist):
                site_band_list.append(zone)
            else:
                if abs(dist) <= on_band_width:
                    site_band_list.append('on')
                elif dist < -on_band_width and dist >= -(on_band_width + out_band_width):
                    site_band_list.append('out')
                elif dist > on_band_width and dist <= (on_band_width + in_band_width):
                    site_band_list.append('in')
                else:
                    site_band_list.append(zone)  # tumor/normal
        adata.obs.loc[sample_mask, 'site_band'] = site_band_list
        # 绘图
        sample_adata = adata[sample_mask]
        coords = sample_adata.obsm['spatial']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        colors_zone = {'tumor': 'blue', 'normal': 'red'}
        for zone in ['tumor', 'normal']:
            mask = sample_adata.obs['region'] == zone
            if mask.any():
                ax1.scatter(coords[mask, 0], coords[mask, 1], c=colors_zone[zone], s=5, alpha=0.6, label=zone)
        ax1.set_title(f'{sample} - Original')
        ax1.legend()
        colors_band = {'on': 'green', 'out': 'orange', 'in': 'purple', 'tumor': 'lightblue', 'normal': 'lightcoral'}
        bands = sample_adata.obs['site_band']
        for band in bands.dropna().unique():
            mask = bands == band
            ax2.scatter(coords[mask, 0], coords[mask, 1], c=colors_band.get(band, 'gray'), s=5, alpha=0.6, label=band)
        ax2.set_title(f'{sample} - Bands')
        ax2.legend()
        plt.suptitle(f'Sample: {sample}')
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_bands.png'), dpi=300)
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_bands.pdf'))
        plt.close()
    return adata

def create_comprehensive_visualization(adata, sample_col='sample', save_dir='make_band',
                                       on_band_width=100, out_band_width=200, in_band_width=200):
    save_dir = create_output_folder(save_dir)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax1 = axes[0, 0]
    if 'site_band' in adata.obs.columns:
        band_counts = adata.obs['site_band'].value_counts()
        colors = {'on': 'green', 'out': 'orange', 'in': 'purple', 'tumor': 'lightblue', 'normal': 'lightcoral'}
        bar_colors = [colors.get(b, 'gray') for b in band_counts.index]
        ax1.bar(band_counts.index.astype(str), band_counts.values, color=bar_colors)
        ax1.set_title('Overall Band Distribution')
    ax2 = axes[0, 1]
    if 'site_band' in adata.obs.columns and sample_col in adata.obs.columns:
        samples = adata.obs[sample_col].unique()[:10]
        all_bands = adata.obs['site_band'].dropna().unique()
        df_plot = pd.DataFrame(index=samples, columns=all_bands, data=0)
        for s in samples:
            mask = adata.obs[sample_col] == s
            counts = adata.obs.loc[mask, 'site_band'].value_counts()
            for b, c in counts.items():
                if b in df_plot.columns:
                    df_plot.loc[s, b] = c
        df_plot.plot(kind='bar', stacked=True, ax=ax2, color=[colors.get(b, 'gray') for b in df_plot.columns])
        ax2.set_title('Band Distribution by Sample')
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
    ax3 = axes[1, 0]
    if 'boundary_distance' in adata.obs.columns:
        dists = adata.obs['boundary_distance'].dropna()
        ax3.hist(dists, bins=50, alpha=0.7)
        ax3.axvline(0, color='red', linestyle='--')
        ax3.axvline(-on_band_width, color='orange', linestyle=':')
        ax3.axvline(on_band_width, color='orange', linestyle=':')
        ax3.set_title('Distance to Boundary')
    ax4 = axes[1, 1]
    if 'boundary_method' in adata.obs.columns:
        method_counts = adata.obs['boundary_method'].value_counts()
        ax4.bar(method_counts.index.astype(str), method_counts.values, color=['blue', 'green', 'orange', 'red'])
        ax4.set_title('Boundary Calculation Methods')
    plt.suptitle('Comprehensive Analysis')
    plt.tight_layout()
    fig.savefig(os.path.join(save_dir, 'comprehensive_analysis.png'), dpi=300)
    fig.savefig(os.path.join(save_dir, 'comprehensive_analysis.pdf'))
    plt.close()

def save_adata_safely(adata, filename, sample_col, save_dir='make_band'):
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
        if key in adata_save.uns:
            del adata_save.uns[key]
    adata_save.write(filepath)
    print(f"数据已保存: {filepath}")
    csv_path = os.path.join(save_dir, 'sample_band_summary.csv')
    summary_cols = ['region', 'site_band', 'boundary_distance', 'is_boundary', 'boundary_method']
    if sample_col in adata.obs.columns:
        summary_cols.append(sample_col)
    existing_cols = [c for c in summary_cols if c in adata.obs.columns]
    adata.obs[existing_cols].to_csv(csv_path)

# ====================== 新增功能1：全样本条带细胞组成统计（celltype + CNs_10） ======================
def band_celltype_statistics(adata, save_dir='band_sta'):
    """
    统计 out/on/in 三个条带中所有样本的：
    1. 细胞类型（celltype列）百分比，绘制堆叠柱状图。
    2. CNs_10 各分类（CN0~CN9）的百分比，绘制堆叠柱状图。
    """
    create_output_folder(save_dir)
    if 'site_band' not in adata.obs.columns:
        print("缺少 'site_band' 列，跳过统计。")
        return
    if 'celltype' not in adata.obs.columns:
        print("缺少 'celltype' 列，跳过细胞类型统计。")
    if 'CNs_10' not in adata.obs.columns:
        print("缺少 'CNs_10' 列，跳过CN统计。")

    bands = ['out', 'on', 'in']
    band_data = adata[adata.obs['site_band'].isin(bands)].obs.copy()
    if band_data.empty:
        print("没有条带细胞，跳过。")
        return

    # ----- 1. 细胞类型堆叠柱状图 -----
    if 'celltype' in band_data.columns:
        ct_counts = band_data.groupby(['site_band', 'celltype']).size().unstack(fill_value=0)
        ct_perc = ct_counts.div(ct_counts.sum(axis=1), axis=0) * 100
        fig, ax = plt.subplots(figsize=(10, 6))
        ct_perc.T.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
        ax.set_title('Cell Type Composition in Each Band (All Samples)')
        ax.set_ylabel('Percentage')
        ax.set_xlabel('Band')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, 'band_celltype_stacked_bar.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_dir, 'band_celltype_stacked_bar.pdf'), bbox_inches='tight')
        plt.close()
        ct_perc.to_csv(os.path.join(save_dir, 'band_celltype_percentage.csv'))

    # ----- 2. CNs_10 各分类堆叠柱状图 -----
    if 'CNs_10' in band_data.columns:
        cn_counts = band_data.groupby(['site_band', 'CNs_10']).size().unstack(fill_value=0)
        cn_perc = cn_counts.div(cn_counts.sum(axis=1), axis=0) * 100
        fig, ax = plt.subplots(figsize=(10, 6))
        cn_perc.T.plot(kind='bar', stacked=True, ax=ax, colormap='tab10')
        ax.set_title('CNs_10 Composition in Each Band (All Samples)')
        ax.set_ylabel('Percentage')
        ax.set_xlabel('Band')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, 'band_CNs10_stacked_bar.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_dir, 'band_CNs10_stacked_bar.pdf'), bbox_inches='tight')
        plt.close()
        cn_perc.to_csv(os.path.join(save_dir, 'band_CNs10_percentage.csv'))

    print(f"全样本条带统计完成，结果保存在 {save_dir}")

# ====================== 新增功能2：分样本分组比较（celltype + CNs_10） ======================
def band_celltype_group_comparison(adata, sample_col='sample', group_col='group_final',
                                   save_dir='band_sta_group'):
    """
    按样本统计 out/on/in 条带的：
    1. 各 celltype 百分比，在 PRE-NR vs PRE-R 两组间用箱线图比较。
    2. 各 CNs_10 分类百分比，在 PRE-NR vs PRE-R 两组间用箱线图比较。
    并标注 Mann-Whitney U 检验 p 值。
    """
    create_output_folder(save_dir)
    if 'site_band' not in adata.obs.columns:
        print("缺少 'site_band' 列，跳过。")
        return
    if group_col not in adata.obs.columns:
        print(f"缺少分组列 '{group_col}'，无法进行分组比较。")
        return
    if sample_col not in adata.obs.columns:
        print(f"缺少样本列 '{sample_col}'，无法按样本统计。")
        return

    bands = ['out', 'on', 'in']
    # 只取条带细胞
    band_mask = adata.obs['site_band'].isin(bands)
    adata_band = adata[band_mask].copy()

    # 获取每个样本的 group_final
    sample_group = adata.obs[[sample_col, group_col]].drop_duplicates()
    compare_groups = ['PRE-NR', 'PRE-R']
    sample_group = sample_group[sample_group[group_col].isin(compare_groups)]

    # ---- 函数：按样本计算某分类列的百分比 ----
    def calc_perc(adata_band, col_name):
        """返回 DataFrame: sample, site_band, category, percentage"""
        df = adata_band.obs.groupby([sample_col, 'site_band', col_name]).size().reset_index(name='count')
        total = adata_band.obs.groupby([sample_col, 'site_band']).size().reset_index(name='total')
        df = df.merge(total, on=[sample_col, 'site_band'])
        df['percentage'] = df['count'] / df['total'] * 100
        return df[[sample_col, 'site_band', col_name, 'percentage']].rename(columns={col_name: 'category'})

    # ---- 通用箱线图绘图函数 ----
    def plot_boxplot(perc_df, feature_name, band, save_prefix):
        """
        perc_df 含有 sample_col, 'category', 'percentage', group_col 列。
        feature_name: 'celltype' 或 'CNs_10'
        """
        categories = perc_df['category'].unique()
        for cat in categories:
            sub = perc_df[perc_df['category'] == cat]
            if sub.empty or sub['percentage'].nunique() < 2:
                continue
            group1 = sub[sub[group_col] == 'PRE-NR']['percentage']
            group2 = sub[sub[group_col] == 'PRE-R']['percentage']
            if len(group1) < 3 or len(group2) < 3:
                continue
            stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
            fig, ax = plt.subplots(figsize=(5, 5))
            data = [group1, group2]
            bp = ax.boxplot(data, labels=['PRE-NR', 'PRE-R'], patch_artist=True)
            for patch, color in zip(bp['boxes'], ['lightblue', 'lightcoral']):
                patch.set_facecolor(color)
            for i, d in enumerate(data):
                x = np.random.normal(i+1, 0.04, size=len(d))
                ax.scatter(x, d, alpha=0.6, color='black', s=20)
            ymax = max(group1.max(), group2.max())
            ymin = min(group1.min(), group2.min())
            y_line = ymax + (ymax - ymin) * 0.05
            ax.plot([1, 2], [y_line, y_line], color='black', linewidth=1)
            sig_text = f'p = {p:.3f}' if p >= 0.001 else 'p < 0.001'
            ax.text(1.5, y_line + (ymax-ymin)*0.01, sig_text, ha='center')
            ax.set_title(f'{feature_name} - {cat} ({band} band)')
            ax.set_ylabel('Percentage')
            plt.tight_layout()
            plt.savefig(os.path.join(save_prefix, f'{cat}_{band}_{feature_name}.png'), dpi=300)
            plt.savefig(os.path.join(save_prefix, f'{cat}_{band}_{feature_name}.pdf'), bbox_inches='tight')
            plt.close()

    # ---- 1. celltype 比较 ----
    if 'celltype' in adata_band.obs.columns:
        perc_celltype = calc_perc(adata_band, 'celltype')
        # 合并分组信息
        perc_celltype = perc_celltype.merge(sample_group, on=sample_col)
        for band in bands:
            band_dir = os.path.join(save_dir, f'celltype_band_{band}')
            create_output_folder(band_dir)
            sub_band = perc_celltype[perc_celltype['site_band'] == band]
            if not sub_band.empty:
                plot_boxplot(sub_band, 'celltype', band, band_dir)

    # ---- 2. CNs_10 比较 ----
    if 'CNs_10' in adata_band.obs.columns:
        perc_cn = calc_perc(adata_band, 'CNs_10')
        perc_cn = perc_cn.merge(sample_group, on=sample_col)
        for band in bands:
            band_dir = os.path.join(save_dir, f'CNs10_band_{band}')
            create_output_folder(band_dir)
            sub_band = perc_cn[perc_cn['site_band'] == band]
            if not sub_band.empty:
                plot_boxplot(sub_band, 'CNs_10', band, band_dir)

    print(f"分组比较完成，结果保存在 {save_dir}")

# ====================== 主流程 ======================
def process_regions_by_sample(adata, sample_col='sample', save_dir='make_band',
                             on_band_width=100, out_band_width=200, in_band_width=200):
    save_dir = create_output_folder(save_dir)
    print("步骤1: 按样本识别tumor-normal分界线...")
    adata = calculate_boundary_by_sample(
        adata, obs_col='region', cat_a='tumor', cat_b='normal',
        sample_col=sample_col, n_neighbors=10, cat_focus='tumor',
        smooth_sigma=3, min_vertex_count=20, max_dist_threshold=500,
        save_dir=save_dir
    )
    print("\n步骤2: 按样本创建条带标签...")
    adata = create_region_bands_by_sample(
        adata, sample_col=sample_col,
        on_band_width=on_band_width, out_band_width=out_band_width, in_band_width=in_band_width,
        save_dir=save_dir
    )
    print("\n步骤3: 创建综合分析图...")
    create_comprehensive_visualization(adata, sample_col=sample_col, save_dir=save_dir,
                                       on_band_width=on_band_width, out_band_width=out_band_width, in_band_width=in_band_width)
    print("\n步骤4: 保存数据...")
    save_adata_safely(adata, 'adata_with_band_by_sample.h5ad', sample_col, save_dir=save_dir)

    # 新增的统计功能
    print("\n步骤5: 全样本条带细胞组成统计...")
    band_celltype_statistics(adata, save_dir='band_sta')
    print("\n步骤6: 分组比较 (PRE-NR vs PRE-R)...")
    band_celltype_group_comparison(adata, sample_col=sample_col, group_col='group_final', save_dir='band_sta_group')

    return adata

def main():
    # 加载数据
    adata = sc.read_h5ad('All_filter_with_region.h5ad')
    # 检查必要列
    if 'region' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'region'列")
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm中未找到'spatial'坐标")

    # 自动识别样本列
    sample_col = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col = col
            print(f"使用 '{col}' 作为样本列")
            break
    if sample_col is None:
        print("未找到样本列，将所有数据视为一个样本")
        sample_col = 'sample'
        adata.obs['sample'] = 'sample_1'

    # 运行主流程（条带宽度已改为 on=100, out=200, in=200）
    adata_result = process_regions_by_sample(
        adata,
        sample_col=sample_col,
        on_band_width=50,    # 边界两侧各100μm
        out_band_width=100,   # 向外延伸200μm
        in_band_width=100,    # 向内延伸200μm
        save_dir='make_band'
    )
    return adata_result

if __name__ == "__main__":
    adata_result = main()
