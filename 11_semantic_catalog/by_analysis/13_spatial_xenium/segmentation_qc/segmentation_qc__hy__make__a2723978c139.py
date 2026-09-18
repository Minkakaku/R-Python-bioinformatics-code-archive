import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree, ConvexHull
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
import os
import warnings
import pickle
warnings.filterwarnings('ignore')

def create_output_folder(folder_name='make_band'):
    """
    创建输出文件夹
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

def calculate_boundary_by_sample(adata, obs_col, cat_a, cat_b, sample_col='sample',
                               n_neighbors=10, cat_focus='inner', 
                               smooth_sigma=3,
                               min_vertex_count=20,
                               max_dist_threshold=500,
                               save_dir='make_band'):
    """
    按样本独立计算tumor和normal的分界线
    """
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
            neighbor_indices = indices[i]
            neighbor_labels = labels[neighbor_indices]
            if current_label == cat_a and cat_b in neighbor_labels:
                boundary_mask[i] = True
            elif current_label == cat_b and cat_a in neighbor_labels:
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
        
        # 方法1: 等高线
        try:
            contours = plt.contour(xx, yy, Z_smooth, levels=[0.5], alpha=0)
            plt.close()
            if len(contours.allsegs) > 0 and len(contours.allsegs[0]) > 0:
                segments = contours.allsegs[0]
                for coords_segment in segments:
                    if len(coords_segment) >= min_vertex_count:
                        valid_boundary_coords.append(coords_segment)
                if valid_boundary_coords:
                    method_used = "contour"
                    print(f"  使用标准等高线法找到边界")
        except Exception as e:
            print(f"  等高线法失败: {e}")
        
        # 方法2: 凸包
        if not valid_boundary_coords:
            print(f"  标准等高线法未找到边界，尝试凸包法...")
            try:
                inner_coords = coords[sample_adata.obs[obs_col] == cat_a]
                if len(inner_coords) > 5:
                    hull = ConvexHull(inner_coords)
                    hull_points = inner_coords[hull.vertices]
                    hull_points_sorted = hull_points[hull.vertices]
                    hull_closed = np.vstack([hull_points_sorted, hull_points_sorted[0]])
                    valid_boundary_coords.append(hull_closed)
                    method_used = "convex_hull"
                    print(f"  使用凸包法找到边界，凸包顶点数: {len(hull_points_sorted)}")
                else:
                    print(f"  {cat_a}细胞数过少({len(inner_coords)})，无法计算凸包")
            except Exception as e:
                print(f"  凸包法失败: {e}")
        
        # 方法3: 边界框
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
        
        if len(valid_boundary_coords) > 0:
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
            
            sample_boundaries[sample] = {
                'coords': valid_boundary_coords,
                'method': method_used
            }
            print(f"  距离计算完成：范围 [{signed_distances.min():.2f}, {signed_distances.max():.2f}]，使用方法: {method_used}")
        else:
            print(f"  警告: 样本 {sample} 无法找到任何有效的分界线")
            sample_indices = adata.obs[sample_col] == sample
            adata.obs.loc[sample_indices, 'is_boundary'] = boundary_mask
            adata.obs.loc[sample_indices, 'boundary_method'] = 'none'
        
        # 绘制当前样本的边界图
        fig, ax = plt.subplots(figsize=(8, 8))
        # 区域颜色：tumor蓝色，normal红色
        for zone, color in [(cat_a, 'blue'), (cat_b, 'red')]:
            mask = sample_adata.obs[obs_col] == zone
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], c=color, 
                          s=5, alpha=0.6, label=zone)
        
        if len(valid_boundary_coords) > 0:
            for coords_segment in valid_boundary_coords:
                if method_used == "contour":
                    line_style = '-'
                    line_label = 'Fitted Boundary (Contour)'
                    line_color = 'black'
                elif method_used == "convex_hull":
                    line_style = '--'
                    line_label = 'Fitted Boundary (Convex Hull)'
                    line_color = 'green'
                else:
                    line_style = ':'
                    line_label = 'Fitted Boundary (Bounding Box)'
                    line_color = 'orange'
                ax.plot(coords_segment[:, 0], coords_segment[:, 1], 
                       line_style, color=line_color, 
                       linewidth=2, label=line_label)
        
        ax.set_title(f'Sample: {sample} - Boundary Fitting\nMethod: {method_used}')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), loc='upper right')
        plt.tight_layout()
        
        fig_path_png = os.path.join(save_dir, f'sample_{sample}_boundary.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_boundary.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        print(f"  保存边界图: {fig_path_png}")
        print(f"  保存边界图: {fig_path_pdf}")
        plt.show()
        plt.close(fig)
    
    print(f"\n处理完成统计:")
    print(f"  成功处理的样本: {len(processed_samples)} 个")
    if processed_samples:
        print(f"  成功处理的样本列表: {', '.join(processed_samples)}")
    print(f"  跳过的样本: {len(skipped_samples)} 个")
    if skipped_samples:
        print(f"  跳过的样本列表: {', '.join(skipped_samples)}")
    
    if 'boundary_method' in adata.obs.columns:
        method_counts = adata.obs['boundary_method'].value_counts()
        print(f"\n边界计算方法统计:")
        for method, count in method_counts.items():
            print(f"  {method}: {count} 个细胞")
    
    boundary_file = os.path.join(save_dir, 'sample_boundaries.pkl')
    with open(boundary_file, 'wb') as f:
        pickle.dump(sample_boundaries, f)
    print(f"\n样本边界坐标已保存为: {boundary_file}")
    
    return adata

def create_region_bands_by_sample(adata, sample_col='sample', 
                                on_band_width=50, out_band_width=100, in_band_width=100,
                                save_dir='make_band'):
    """
    按样本独立创建条带标签，并将region列改为region（tumor/normal）
    """
    if 'boundary_distance' not in adata.obs:
        raise ValueError("请先运行 calculate_boundary_by_sample 函数计算分界线距离")
    
    save_dir = create_output_folder(save_dir)
    adata.obs['site_band'] = np.nan
    
    samples = adata.obs[sample_col].unique()
    processed_samples = []
    skipped_samples = []
    
    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        if not sample_mask.any():
            continue
            
        sample_distances = adata.obs.loc[sample_mask, 'boundary_distance'].values
        sample_zones = adata.obs.loc[sample_mask, 'region'].values  # 改为region列
        
        if np.isnan(sample_distances).all():
            print(f"样本 {sample} 未计算距离，跳过")
            skipped_samples.append(sample)
            continue
            
        processed_samples.append(sample)
        
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
                    site_band_list.append(zone)
        
        adata.obs.loc[sample_mask, 'site_band'] = site_band_list
        
        band_counts = pd.Series(site_band_list).value_counts()
        print(f"\n样本 {sample} 条带细胞分布统计:")
        for band, count in band_counts.items():
            print(f"  {band}: {count} 个细胞")
        
        # 绘图
        sample_adata = adata[sample_mask]
        coords = sample_adata.obsm['spatial']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # 左图：原始区域分布 (tumor/normal)
        colors_zone = {'tumor': 'blue', 'normal': 'red'}
        for zone in ['tumor', 'normal']:
            mask = sample_adata.obs['region'] == zone
            if mask.any():
                ax1.scatter(coords[mask, 0], coords[mask, 1], c=colors_zone[zone], 
                           s=5, alpha=0.6, label=zone)
        ax1.set_title(f'Sample: {sample} - Original Regions')
        ax1.set_xlabel('X coordinate')
        ax1.set_ylabel('Y coordinate')
        ax1.legend()
        
        # 右图：条带分布
        colors_band = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                      'tumor': 'lightblue', 'normal': 'lightcoral'}
        bands = sample_adata.obs['site_band']
        if bands.notna().any():
            unique_bands = bands.dropna().unique()
            for band in unique_bands:
                mask = bands == band
                if mask.any():
                    ax2.scatter(coords[mask, 0], coords[mask, 1], 
                               c=colors_band.get(band, 'gray'), 
                               s=5, alpha=0.6, label=band)
        ax2.set_title(f'Sample: {sample} - Band Zones')
        ax2.set_xlabel('X coordinate')
        ax2.set_ylabel('Y coordinate')
        
        handles, labels = ax2.get_legend_handles_labels()
        unique_labels = []
        unique_handles = []
        for handle, label in zip(handles, labels):
            if label not in unique_labels:
                unique_labels.append(label)
                unique_handles.append(handle)
        if unique_labels:
            ax2.legend(unique_handles, unique_labels)
        
        plt.suptitle(f'Sample: {sample} - Region Analysis', fontsize=16)
        plt.tight_layout()
        
        fig_path_png = os.path.join(save_dir, f'sample_{sample}_bands.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_bands.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        print(f"  保存条带区域图: {fig_path_png}")
        print(f"  保存条带区域图: {fig_path_pdf}")
        plt.show()
        plt.close(fig)
    
    print(f"\n条带创建完成统计:")
    print(f"  成功处理的样本: {len(processed_samples)} 个")
    if processed_samples:
        print(f"  成功处理的样本列表: {', '.join(processed_samples)}")
    print(f"  跳过的样本: {len(skipped_samples)} 个")
    if skipped_samples:
        print(f"  跳过的样本列表: {', '.join(skipped_samples)}")
    
    return adata

def create_comprehensive_visualization(adata, sample_col='sample', save_dir='make_band',
                                      on_band_width=50, out_band_width=100, in_band_width=100):
    """
    创建综合可视化图形，包括所有样本的条带统计
    """
    save_dir = create_output_folder(save_dir)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. 所有样本的条带分布统计
    ax1 = axes[0, 0]
    if 'site_band' in adata.obs.columns:
        band_counts = adata.obs['site_band'].value_counts()
        colors = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                  'tumor': 'lightblue', 'normal': 'lightcoral'}
        bar_colors = [colors.get(band, 'gray') for band in band_counts.index]
        bars = ax1.bar(band_counts.index.astype(str), band_counts.values, color=bar_colors)
        ax1.set_title('Overall Band Distribution')
        ax1.set_xlabel('Band')
        ax1.set_ylabel('Cell Count')
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 5,
                    f'{int(height)}', ha='center', va='bottom')
    
    # 2. 每个样本的条带分布
    ax2 = axes[0, 1]
    if 'site_band' in adata.obs.columns and sample_col in adata.obs.columns:
        samples = adata.obs[sample_col].unique()
        band_colors = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                       'tumor': 'lightblue', 'normal': 'lightcoral'}
        all_bands = adata.obs['site_band'].dropna().unique()
        sample_data = []
        sample_names = []
        for sample in samples[:10]:
            sample_mask = adata.obs[sample_col] == sample
            if sample_mask.any():
                sample_bands = adata.obs.loc[sample_mask, 'site_band'].value_counts()
                sample_data.append(sample_bands)
                sample_names.append(str(sample))
        if sample_data:
            df = pd.DataFrame(sample_data).fillna(0)
            bottom = np.zeros(len(df))
            for band in all_bands:
                if band in df.columns:
                    ax2.bar(range(len(df)), df[band], bottom=bottom, 
                           label=band, color=band_colors.get(band, 'gray'))
                    bottom += df[band].values
            ax2.set_title('Band Distribution by Sample')
            ax2.set_xlabel('Sample')
            ax2.set_ylabel('Cell Count')
            ax2.set_xticks(range(len(sample_names)))
            ax2.set_xticklabels(sample_names, rotation=45, ha='right')
            if len(all_bands) <= 10:
                ax2.legend()
    
    # 3. 距离分布直方图
    ax3 = axes[1, 0]
    if 'boundary_distance' in adata.obs.columns:
        distances = adata.obs['boundary_distance'].dropna()
        if len(distances) > 0:
            ax3.hist(distances, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Boundary (0)')
            ax3.axvline(x=-on_band_width, color='orange', linestyle=':', linewidth=1, label=f'Outer edge (-{on_band_width})')
            ax3.axvline(x=on_band_width, color='orange', linestyle=':', linewidth=1, label=f'Inner edge (+{on_band_width})')
            ax3.set_title('Distribution of Distance to Boundary')
            ax3.set_xlabel('Distance to Boundary (um)')
            ax3.set_ylabel('Frequency')
            ax3.legend()
    
    # 4. 边界计算方法统计
    ax4 = axes[1, 1]
    if 'boundary_method' in adata.obs.columns and sample_col in adata.obs.columns:
        method_counts = adata.obs['boundary_method'].value_counts()
        if len(method_counts) > 0:
            colors = {'contour': 'blue', 'convex_hull': 'green', 
                     'bounding_box': 'orange', 'none': 'red'}
            bar_colors = [colors.get(method, 'gray') for method in method_counts.index]
            bars = ax4.bar(method_counts.index.astype(str), method_counts.values, color=bar_colors)
            ax4.set_title('Boundary Calculation Methods')
            ax4.set_xlabel('Method')
            ax4.set_ylabel('Number of Cells')
            ax4.set_xticklabels(method_counts.index, rotation=45, ha='right')
            for bar, count in zip(bars, method_counts.values):
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2., height + 5,
                        f'{int(count)}', ha='center', va='bottom')
    
    plt.suptitle('Comprehensive Analysis Results', fontsize=16)
    plt.tight_layout()
    
    fig3_path_png = os.path.join(save_dir, 'comprehensive_analysis.png')
    fig3_path_pdf = os.path.join(save_dir, 'comprehensive_analysis.pdf')
    fig.savefig(fig3_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(fig3_path_pdf, bbox_inches='tight')
    print(f"保存综合分析图: {fig3_path_png}")
    print(f"保存综合分析图: {fig3_path_pdf}")
    plt.show()
    plt.close(fig)

def save_adata_safely(adata, filename, save_dir='make_band'):
    """
    安全保存adata，避免不规则数组导致的错误
    """
    save_dir = create_output_folder(save_dir)
    filepath = os.path.join(save_dir, filename)
    adata_save = adata.copy()
    
    keys_to_remove = []
    for key in adata_save.uns.keys():
        try:
            test_value = adata_save.uns[key]
            if isinstance(test_value, (list, np.ndarray)):
                if isinstance(test_value, list) and len(test_value) > 0:
                    if all(isinstance(x, np.ndarray) for x in test_value):
                        shapes = [x.shape for x in test_value]
                        if len(set(shapes)) > 1:
                            keys_to_remove.append(key)
        except:
            keys_to_remove.append(key)
    
    if 'sample_boundaries' in adata_save.uns:
        keys_to_remove.append('sample_boundaries')
    
    for key in keys_to_remove:
        print(f"警告：移除不规则的uns键: {key}")
        if key in adata_save.uns:
            del adata_save.uns[key]
    
    adata_save.write(filepath)
    print(f"数据已安全保存为 '{filepath}'")
    
    csv_path = os.path.join(save_dir, 'sample_band_summary.csv')
    summary_cols = ['region', 'site_band', 'boundary_distance', 'is_boundary', 'boundary_method']
    if sample_col in adata.obs.columns:
        summary_cols.append(sample_col)
    existing_cols = [col for col in summary_cols if col in adata.obs.columns]
    summary_df = adata.obs[existing_cols].copy()
    summary_df.to_csv(csv_path, index=False)
    print(f"数据摘要已保存为 '{csv_path}'")

def process_regions_by_sample(adata, sample_col='sample', save_dir='make_band',
                             on_band_width=50, out_band_width=100, in_band_width=100):
    """
    一站式处理函数，按样本计算tumor-normal边界及条带
    """
    save_dir = create_output_folder(save_dir)
    
    print("步骤1: 按样本识别tumor-normal分界线...")
    adata = calculate_boundary_by_sample(
        adata, 
        obs_col='region',
        sample_col=sample_col,
        cat_a='tumor', 
        cat_b='normal',
        n_neighbors=10,
        cat_focus='tumor',
        smooth_sigma=3,
        min_vertex_count=20,
        max_dist_threshold=500,
        save_dir=save_dir
    )
    
    print("\n步骤2: 按样本创建条带标签...")
    adata = create_region_bands_by_sample(
        adata,
        sample_col=sample_col,
        on_band_width=on_band_width,
        out_band_width=out_band_width,
        in_band_width=in_band_width,
        save_dir=save_dir
    )
    
    print("\n步骤3: 创建综合分析图...")
    create_comprehensive_visualization(
        adata, 
        sample_col=sample_col, 
        save_dir=save_dir,
        on_band_width=on_band_width,
        out_band_width=out_band_width,
        in_band_width=in_band_width
    )
    
    print("\n步骤4: 保存数据...")
    save_adata_safely(adata, 'adata_with_band_by_sample.h5ad', save_dir=save_dir)
    
    print("\n最终的adata.obs列:")
    summary_cols = ['region', 'site_band', 'boundary_distance', 'is_boundary', 'boundary_method']
    if sample_col in adata.obs.columns:
        summary_cols.append(sample_col)
    existing_cols = [col for col in summary_cols if col in adata.obs.columns]
    print(adata.obs[existing_cols].head())
    
    save_statistics_summary(adata, sample_col, save_dir)
    
    return adata

def save_statistics_summary(adata, sample_col, save_dir):
    """
    保存统计摘要
    """
    save_dir = create_output_folder(save_dir)
    summary_path = os.path.join(save_dir, 'analysis_summary.txt')
    
    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("         Region (Tumor/Normal) Boundary Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("1. OVERALL STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total cells: {len(adata)}\n")
        
        if 'region' in adata.obs.columns:
            zone_counts = adata.obs['region'].value_counts()
            f.write("\nRegion Distribution:\n")
            for zone, count in zone_counts.items():
                percentage = count / len(adata) * 100
                f.write(f"  {zone}: {count} ({percentage:.1f}%)\n")
        
        if 'site_band' in adata.obs.columns:
            band_counts = adata.obs['site_band'].value_counts()
            f.write("\nSite Band Distribution:\n")
            for band, count in band_counts.items():
                if pd.notna(band):
                    percentage = count / len(adata) * 100
                    f.write(f"  {band}: {count} ({percentage:.1f}%)\n")
        
        if 'is_boundary' in adata.obs.columns:
            boundary_count = adata.obs['is_boundary'].sum()
            f.write(f"\nBoundary cells: {boundary_count} ({boundary_count/len(adata)*100:.1f}%)\n")
        
        if 'boundary_method' in adata.obs.columns:
            method_counts = adata.obs['boundary_method'].value_counts()
            f.write("\nBoundary Calculation Methods:\n")
            for method, count in method_counts.items():
                percentage = count / len(adata) * 100
                f.write(f"  {method}: {count} ({percentage:.1f}%)\n")
        
        if sample_col in adata.obs.columns:
            f.write("\n\n2. SAMPLE-WISE STATISTICS\n")
            f.write("-" * 40 + "\n")
            samples = adata.obs[sample_col].unique()
            for sample in samples:
                f.write(f"\nSample: {sample}\n")
                sample_mask = adata.obs[sample_col] == sample
                sample_cells = adata[sample_mask]
                f.write(f"  Total cells: {len(sample_cells)}\n")
                if 'region' in sample_cells.obs.columns:
                    zone_counts = sample_cells.obs['region'].value_counts()
                    for zone, count in zone_counts.items():
                        percentage = count / len(sample_cells) * 100
                        f.write(f"    {zone}: {count} ({percentage:.1f}%)\n")
                if 'site_band' in sample_cells.obs.columns:
                    band_counts = sample_cells.obs['site_band'].value_counts()
                    for band, count in band_counts.items():
                        if pd.notna(band):
                            percentage = count / len(sample_cells) * 100
                            f.write(f"    {band}: {count} ({percentage:.1f}%)\n")
                if 'boundary_method' in sample_cells.obs.columns:
                    method_counts = sample_cells.obs['boundary_method'].value_counts()
                    for method, count in method_counts.items():
                        f.write(f"    Boundary method: {method}: {count} cells\n")
        
        if 'boundary_distance' in adata.obs.columns:
            distances = adata.obs['boundary_distance'].dropna()
            if len(distances) > 0:
                f.write("\n\n3. DISTANCE STATISTICS\n")
                f.write("-" * 40 + "\n")
                f.write(f"Min distance: {distances.min():.2f} um\n")
                f.write(f"Max distance: {distances.max():.2f} um\n")
                f.write(f"Mean distance: {distances.mean():.2f} um\n")
                f.write(f"Std distance: {distances.std():.2f} um\n")
                f.write(f"Median distance: {distances.median():.2f} um\n")
                f.write("\nDistance quantiles:\n")
                for q in [0.1, 0.25, 0.5, 0.75, 0.9]:
                    f.write(f"  {q*100:.0f}%: {distances.quantile(q):.2f} um\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Analysis completed successfully!\n")
        f.write("=" * 60 + "\n")
    
    print(f"分析摘要已保存为: {summary_path}")

def main():
    # 加载数据（修改为您的输入文件）
    adata = sc.read_h5ad('All_filter_with_region.h5ad')
    
    # 检查必要列
    if 'region' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'region'列，请确保包含tumor/normal信息")
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm中未找到'spatial'坐标信息")
    
    # 自动识别样本列
    sample_col = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col = col
            print(f"使用 '{col}' 作为样本列")
            break
    if sample_col is None:
        print("警告: 未找到样本列，将使用整个数据集作为一个样本处理")
        sample_col = 'sample'
        adata.obs['sample'] = 'sample_1'
    
    # 一站式处理，输出到 make_band 文件夹
    adata_result = process_regions_by_sample(
        adata, 
        sample_col=sample_col,
        on_band_width=50,
        out_band_width=100,
        in_band_width=100,
        save_dir='make_band'
    )
    
    return adata_result

if __name__ == "__main__":
    adata_result = main()
