import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree, ConvexHull
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
from sklearn.cluster import DBSCAN
import os
import warnings
import pickle

warnings.filterwarnings('ignore')

def create_output_folder(folder_name='make_band_s_den'):
    """创建输出文件夹"""
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

# -------------------- 边界识别、条带创建、面积计算等函数（保持不变） --------------------
def calculate_boundary_by_sample(adata, obs_col, cat_a, cat_b, sample_col='sample',
                               n_neighbors=10, cat_focus='inner', 
                               smooth_sigma=3,
                               min_vertex_count=20,
                               max_dist_threshold=500,
                               save_dir='make_band_new_2'):
    """按样本独立计算inner和outer的分界线（同原文）"""
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
            print(f"  尝试凸包法...")
            try:
                inner_coords = coords[sample_adata.obs[obs_col] == cat_a]
                if len(inner_coords) > 5:
                    hull = ConvexHull(inner_coords)
                    hull_points = inner_coords[hull.vertices]
                    hull_closed = np.vstack([hull_points, hull_points[0]])
                    valid_boundary_coords.append(hull_closed)
                    method_used = "convex_hull"
                    print(f"  使用凸包法找到边界，顶点数: {len(hull_points)}")
            except Exception as e:
                print(f"  凸包法失败: {e}")
        
        # 方法3: 边界框
        if not valid_boundary_coords:
            print(f"  尝试边界框法...")
            try:
                inner_coords = coords[sample_adata.obs[obs_col] == cat_a]
                if len(inner_coords) > 0:
                    x_min_i, x_max_i = inner_coords[:, 0].min(), inner_coords[:, 0].max()
                    y_min_i, y_max_i = inner_coords[:, 1].min(), inner_coords[:, 1].max()
                    padding = 50
                    box_points = np.array([
                        [x_min_i - padding, y_min_i - padding],
                        [x_max_i + padding, y_min_i - padding],
                        [x_max_i + padding, y_max_i + padding],
                        [x_min_i - padding, y_max_i + padding],
                        [x_min_i - padding, y_min_i - padding]
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
            sample_boundaries[sample] = {'coords': valid_boundary_coords, 'method': method_used}
            print(f"  距离计算完成: [{signed_distances.min():.2f}, {signed_distances.max():.2f}]，方法: {method_used}")
        else:
            print(f"  警告: 样本 {sample} 无法找到任何有效的分界线")
            sample_indices = adata.obs[sample_col] == sample
            adata.obs.loc[sample_indices, 'is_boundary'] = boundary_mask
            adata.obs.loc[sample_indices, 'boundary_method'] = 'none'
        
        # 绘制当前样本边界图
        fig, ax = plt.subplots(figsize=(8, 8))
        for zone, color in [(cat_a, 'blue'), (cat_b, 'red')]:
            mask = sample_adata.obs[obs_col] == zone
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], c=color, s=5, alpha=0.6, label=zone)
        if len(valid_boundary_coords) > 0:
            for seg in valid_boundary_coords:
                if method_used == "contour":
                    ls, col = '-', 'black'
                elif method_used == "convex_hull":
                    ls, col = '--', 'green'
                else:
                    ls, col = ':', 'orange'
                ax.plot(seg[:, 0], seg[:, 1], ls, color=col, linewidth=2, label=f'Boundary ({method_used})')
        ax.set_title(f'Sample: {sample} - Boundary Fitting\nMethod: {method_used}')
        ax.set_xlabel('X'); ax.set_ylabel('Y')
        if ax.get_legend_handles_labels()[1]:
            ax.legend(loc='upper right')
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_boundary.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_boundary.pdf'), bbox_inches='tight')
        plt.show(); plt.close(fig)
    
    print(f"\n处理完成: 成功 {len(processed_samples)} 个, 跳过 {len(skipped_samples)} 个")
    boundary_file = os.path.join(save_dir, 'sample_boundaries.pkl')
    with open(boundary_file, 'wb') as f:
        pickle.dump(sample_boundaries, f)
    print(f"样本边界坐标已保存: {boundary_file}")
    return adata

def compute_area_from_coords(coords, method='cluster_convex', eps=100, min_samples=3):
    """根据二维坐标点集合计算面积（单位：坐标单位的平方）"""
    if len(coords) < 3:
        return 0.0
    try:
        if method == 'convex_hull':
            hull = ConvexHull(coords)
            return hull.volume
        elif method == 'cluster_convex':
            db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
            labels = db.labels_
            unique_labels = set(labels)
            total_area = 0.0
            for lab in unique_labels:
                if lab == -1:
                    continue
                cluster_pts = coords[labels == lab]
                if len(cluster_pts) >= 3:
                    try:
                        hull = ConvexHull(cluster_pts)
                        total_area += hull.volume
                    except:
                        pass
            noise_pts = coords[labels == -1]
            if len(noise_pts) >= 3:
                try:
                    hull_noise = ConvexHull(noise_pts)
                    total_area += hull_noise.volume
                except:
                    pass
            return total_area
        elif method == 'bounding_box':
            xmin, xmax = coords[:,0].min(), coords[:,0].max()
            ymin, ymax = coords[:,1].min(), coords[:,1].max()
            return (xmax - xmin) * (ymax - ymin)
        else:
            return 0.0
    except Exception as e:
        try:
            xmin, xmax = coords[:,0].min(), coords[:,0].max()
            ymin, ymax = coords[:,1].min(), coords[:,1].max()
            return (xmax - xmin) * (ymax - ymin)
        except:
            return 0.0

def create_mucus_bands_by_sample(adata, sample_col='sample',
                                on_band_width=50, out_band_width=100, in_band_width=100,
                                area_method='cluster_convex',
                                save_dir='make_band_new_2'):
    """按样本独立创建条带标签，并计算每个条带的面积"""
    if 'boundary_distance' not in adata.obs:
        raise ValueError("请先运行 calculate_boundary_by_sample 函数计算分界线距离")
    
    save_dir = create_output_folder(save_dir)
    adata.obs['site_band'] = np.nan
    samples = adata.obs[sample_col].unique()
    
    processed_samples = []
    skipped_samples = []
    
    if 'band_areas' not in adata.uns:
        adata.uns['band_areas'] = {}
    
    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        if not sample_mask.any():
            continue
        sample_distances = adata.obs.loc[sample_mask, 'boundary_distance'].values
        sample_zones = adata.obs.loc[sample_mask, 'mucus_zone'].values
        
        if np.isnan(sample_distances).all():
            print(f"样本 {sample} 未计算距离，跳过")
            skipped_samples.append(sample)
            continue
        processed_samples.append(sample)
        
        # 分配条带标签
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
        
        # 计算各条带面积
        sample_adata = adata[sample_mask].copy()
        coords = sample_adata.obsm['spatial']
        bands_in_sample = sample_adata.obs['site_band'].unique()
        band_areas = {}
        
        print(f"\n  计算样本 {sample} 各条带面积（方法: {area_method}）:")
        for band in bands_in_sample:
            if pd.isna(band):
                continue
            mask_band = sample_adata.obs['site_band'] == band
            band_coords = coords[mask_band]
            if len(band_coords) >= 3:
                area = compute_area_from_coords(band_coords, method=area_method, eps=100, min_samples=3)
            else:
                area = 0.0
            band_areas[band] = area
            area_mm2 = area / 1e6
            print(f"    {band}: {area:.2f} μm² ({area_mm2:.4f} mm²)")
        
        adata.uns['band_areas'][sample] = band_areas
        
        # 绘制条带区域图
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        colors_zone = {'inner': 'blue', 'outer': 'red'}
        for zone in ['inner', 'outer']:
            mask = sample_adata.obs['mucus_zone'] == zone
            if mask.any():
                ax1.scatter(coords[mask, 0], coords[mask, 1], c=colors_zone[zone],
                           s=5, alpha=0.6, label=zone)
        ax1.set_title(f'Sample: {sample} - Original Zones')
        ax1.set_xlabel('X'); ax1.set_ylabel('Y')
        ax1.legend()
        
        colors_band = {'on': 'green', 'out': 'orange', 'in': 'purple',
                      'inner': 'lightblue', 'outer': 'lightcoral'}
        bands = sample_adata.obs['site_band']
        if bands.notna().any():
            unique_bands = bands.dropna().unique()
            for band in unique_bands:
                mask = bands == band
                if mask.any():
                    ax2.scatter(coords[mask, 0], coords[mask, 1],
                               c=colors_band.get(band, 'gray'),
                               s=5, alpha=0.6, label=band)
        ax2.set_title(f'Sample: {sample} - Band Zones\nAreas (mm²): ' + 
                     ', '.join([f"{b}: {band_areas.get(b,0)/1e6:.2f}" for b in unique_bands if not pd.isna(b)]),
                     fontsize=10)
        ax2.set_xlabel('X'); ax2.set_ylabel('Y')
        handles, labels = ax2.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax2.legend(by_label.values(), by_label.keys(), loc='upper right')
        
        plt.suptitle(f'Sample: {sample} - Zone Analysis', fontsize=16)
        plt.tight_layout()
        
        fig_path_png = os.path.join(save_dir, f'sample_{sample}_bands.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_bands.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        print(f"  保存条带区域图: {fig_path_png}")
        plt.show(); plt.close(fig)
    
    print(f"\n条带创建完成: 成功 {len(processed_samples)} 个, 跳过 {len(skipped_samples)} 个")
    
    # 导出面积CSV
    area_records = []
    for sample, areas in adata.uns['band_areas'].items():
        for band, area in areas.items():
            area_records.append({'sample': sample, 'band': band, 'area_um2': area, 'area_mm2': area/1e6})
    area_df = pd.DataFrame(area_records)
    area_csv_path = os.path.join(save_dir, 'band_areas.csv')
    area_df.to_csv(area_csv_path, index=False)
    print(f"条带面积统计已保存为: {area_csv_path}")
    
    return adata

# ==================== 细胞密度计算与可视化 ====================
def calculate_cell_density_by_band(adata, sample_col, cell_type_col, save_dir):
    """
    基于条带划分和面积，计算每个样本、每个条带、每种细胞类型的细胞密度。
    密度单位：细胞数/mm²（1 mm² = 1e6 μm²）。
    结果存储在 adata.uns['cell_density'] 中，并导出CSV、绘制热图。
    """
    if cell_type_col not in adata.obs.columns:
        print(f"警告: 细胞类型列 '{cell_type_col}' 不存在，跳过密度计算。")
        return
    
    if 'band_areas' not in adata.uns:
        print("警告: 未找到条带面积信息，跳过密度计算。")
        return
    
    # 确保 site_band 列存在
    if 'site_band' not in adata.obs.columns:
        print("警告: 未找到 site_band 列，跳过密度计算。")
        return
    
    # 收集所有密度记录
    density_records = []
    
    # 遍历每个样本
    for sample in adata.obs[sample_col].unique():
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask]
        
        # 获取该样本的条带面积
        if sample not in adata.uns['band_areas']:
            print(f"警告: 样本 {sample} 缺少面积信息，跳过该样本密度计算。")
            continue
        band_areas = adata.uns['band_areas'][sample]
        
        # 按条带分组统计细胞类型数量
        grouped = sample_adata.obs.groupby(['site_band', cell_type_col]).size().reset_index(name='count')
        
        for _, row in grouped.iterrows():
            band = row['site_band']
            cell_type = row[cell_type_col]
            count = row['count']
            # 获取该条带的面积（μm²），若不存在则跳过
            area_um2 = band_areas.get(band, np.nan)
            if pd.isna(area_um2) or area_um2 == 0:
                continue
            # 密度 = 细胞数 / (面积/1e6) = 细胞数 * 1e6 / 面积
            density = count * 1e6 / area_um2
            density_records.append({
                'sample': sample,
                'band': band,
                'cell_type': cell_type,
                'count': count,
                'area_um2': area_um2,
                'density_mm2': density
            })
    
    if not density_records:
        print("没有生成任何密度记录，跳过。")
        return
    
    density_df = pd.DataFrame(density_records)
    
    # 保存到 uns
    adata.uns['cell_density'] = density_df
    
    # 导出为 CSV
    csv_path = os.path.join(save_dir, 'cell_density_by_band.csv')
    density_df.to_csv(csv_path, index=False)
    print(f"细胞密度统计已保存: {csv_path}")
    
    # 可视化：热图（样本 × 条带-细胞类型 密度矩阵）
    # 创建用于热图的透视表：行 = 样本，列 = 条带+细胞类型，值 = 密度
    density_df['band_cell'] = density_df['band'] + '_' + density_df['cell_type']
    pivot = density_df.pivot_table(index='sample', columns='band_cell', values='density_mm2', fill_value=0)
    
    # 如果样本或列太多，热图可能太大，可以适当调整
    if pivot.shape[0] > 0 and pivot.shape[1] > 0:
        # 设置绘图尺寸
        fig_width = max(10, pivot.shape[1] * 0.4)
        fig_height = max(6, pivot.shape[0] * 0.4)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        sns.heatmap(pivot, annot=False, cmap='viridis', cbar_kws={'label': 'Density (cells/mm²)'},
                   ax=ax)
        ax.set_title('Cell Density per Band and Cell Type')
        plt.tight_layout()
        
        heatmap_path_png = os.path.join(save_dir, 'cell_density_heatmap.png')
        heatmap_path_pdf = os.path.join(save_dir, 'cell_density_heatmap.pdf')
        fig.savefig(heatmap_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(heatmap_path_pdf, bbox_inches='tight')
        print(f"密度热图已保存: {heatmap_path_png}")
        plt.show()
        plt.close(fig)
    else:
        print("密度数据透视表为空，无法绘制热图。")

# ==================== 新增：细胞类型箱线图函数 ====================
def plot_cell_type_boxplots(adata, save_dir, metric='density_mm2'):
    """
    为每个细胞类型绘制箱线图，比较不同条带间的指标（默认密度）。
    数据来自 adata.uns['cell_density']。
    """
    if 'cell_density' not in adata.uns:
        print("未找到细胞密度数据，请先运行 calculate_cell_density_by_band")
        return
    density_df = adata.uns['cell_density']
    required = ['cell_type', 'band', metric]
    if not all(col in density_df.columns for col in required):
        print(f"密度数据框缺少列: {required}")
        return

    cell_types = density_df['cell_type'].unique()
    n_types = len(cell_types)
    if n_types == 0:
        print("无细胞类型数据")
        return

    ncols = min(4, n_types)
    nrows = (n_types + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes = axes.flatten() if nrows*ncols > 1 else [axes]

    for i, ct in enumerate(cell_types):
        ax = axes[i]
        ct_data = density_df[density_df['cell_type'] == ct]
        bands = ct_data['band'].unique()
        data_to_plot = [ct_data[ct_data['band'] == b][metric].values for b in bands]

        bp = ax.boxplot(data_to_plot, labels=bands, patch_artist=True)
        # 为盒子设置不同颜色
        colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow', 'lightpink']
        for patch, color in zip(bp['boxes'], colors * (len(bands)//len(colors)+1)):
            patch.set_facecolor(color)
        ax.set_title(f'Cell Type: {ct}')
        ax.set_ylabel(metric)
        ax.set_xlabel('Band')

    # 隐藏多余子图
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f'Boxplots of {metric} by Band per Cell Type')
    plt.tight_layout()
    fig_path_png = os.path.join(save_dir, f'boxplots_{metric}_by_band_per_celltype.png')
    fig_path_pdf = os.path.join(save_dir, f'boxplots_{metric}_by_band_per_celltype.pdf')
    fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    print(f"箱线图已保存: {fig_path_png}")
    plt.show()
    plt.close(fig)

# ==================== 修改综合可视化 ====================
def create_comprehensive_visualization(adata, sample_col='sample', save_dir='make_band_new_2',
                                      on_band_width=50, out_band_width=100, in_band_width=100):
    """创建综合可视化，包括细胞密度热图（如果存在）"""
    save_dir = create_output_folder(save_dir)
    
    # 判断是否存在密度数据，若存在则使用 2x3 布局，否则保持 2x2
    has_density = 'cell_density' in adata.uns
    if has_density:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
    
    ax_idx = 0
    
    # 1. 条带细胞数量分布
    if has_density:
        ax1 = axes[0, 0]
    else:
        ax1 = axes[ax_idx]; ax_idx += 1
    if 'site_band' in adata.obs.columns:
        band_counts = adata.obs['site_band'].value_counts()
        colors = {'on': 'green', 'out': 'orange', 'in': 'purple',
                  'inner': 'lightblue', 'outer': 'lightcoral'}
        bar_colors = [colors.get(b, 'gray') for b in band_counts.index]
        bars = ax1.bar(band_counts.index.astype(str), band_counts.values, color=bar_colors)
        ax1.set_title('Overall Band Distribution (Cell Count)')
        ax1.set_xlabel('Band'); ax1.set_ylabel('Cell Count')
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x()+bar.get_width()/2, height+5, f'{int(height)}', ha='center', va='bottom')
    
    # 2. 每个样本的条带细胞数量堆叠
    if has_density:
        ax2 = axes[0, 1]
    else:
        ax2 = axes[ax_idx]; ax_idx += 1
    if 'site_band' in adata.obs.columns and sample_col in adata.obs.columns:
        samples = adata.obs[sample_col].unique()[:10]
        band_colors = {'on': 'green', 'out': 'orange', 'in': 'purple',
                       'inner': 'lightblue', 'outer': 'lightcoral'}
        all_bands = adata.obs['site_band'].dropna().unique()
        sample_data = []
        sample_names = []
        for sample in samples:
            mask = adata.obs[sample_col] == sample
            if mask.any():
                sample_bands = adata.obs.loc[mask, 'site_band'].value_counts()
                sample_data.append(sample_bands)
                sample_names.append(str(sample))
        if sample_data:
            df = pd.DataFrame(sample_data).fillna(0)
            bottom = np.zeros(len(df))
            for band in all_bands:
                if band in df.columns:
                    ax2.bar(range(len(df)), df[band], bottom=bottom, label=band,
                           color=band_colors.get(band, 'gray'))
                    bottom += df[band].values
            ax2.set_title('Cell Count per Band by Sample')
            ax2.set_xlabel('Sample'); ax2.set_ylabel('Cell Count')
            ax2.set_xticks(range(len(sample_names)))
            ax2.set_xticklabels(sample_names, rotation=45, ha='right')
            if len(all_bands) <= 10:
                ax2.legend()
    
    # 3. 距离分布直方图
    if has_density:
        ax3 = axes[0, 2]
    else:
        ax3 = axes[ax_idx]; ax_idx += 1
    if 'boundary_distance' in adata.obs.columns:
        distances = adata.obs['boundary_distance'].dropna()
        if len(distances) > 0:
            ax3.hist(distances, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Boundary (0)')
            ax3.axvline(x=-on_band_width, color='orange', linestyle=':', linewidth=1,
                       label=f'Outer edge (-{on_band_width})')
            ax3.axvline(x=on_band_width, color='orange', linestyle=':', linewidth=1,
                       label=f'Inner edge (+{on_band_width})')
            ax3.set_title('Distribution of Distance to Boundary')
            ax3.set_xlabel('Distance (μm)'); ax3.set_ylabel('Frequency')
            ax3.legend()
    
    # 4. 条带总面积分布
    if has_density:
        ax4 = axes[1, 0]
    else:
        ax4 = axes[ax_idx]; ax_idx += 1
    if 'band_areas' in adata.uns and adata.uns['band_areas']:
        area_summary = {}
        for sample, areas in adata.uns['band_areas'].items():
            for band, area in areas.items():
                if band not in area_summary:
                    area_summary[band] = 0
                area_summary[band] += area
        bands = list(area_summary.keys())
        areas_mm2 = [a/1e6 for a in area_summary.values()]
        colors_pie = [colors.get(b, 'gray') for b in bands] if 'colors' in locals() else None
        wedges, texts, autotexts = ax4.pie(areas_mm2, labels=bands, autopct='%1.1f%%',
                                           colors=colors_pie, startangle=90)
        ax4.set_title(f'Total Band Area (mm²): {sum(areas_mm2):.2f}')
    else:
        ax4.text(0.5, 0.5, 'No area data', ha='center', va='center')
        ax4.set_title('Band Area Distribution')
    
    # 5. 如果有密度数据，添加密度热图（或提示）
    if has_density:
        ax5 = axes[1, 1]  # 第2行第2列
        ax5.text(0.5, 0.5, 'Cell density heatmap\nsaved separately', ha='center', va='center', fontsize=12)
        ax5.set_title('Density Visualization')
        ax5.axis('off')
        
        ax6 = axes[1, 2]
        ax6.text(0.5, 0.5, 'See cell_density_heatmap.png\nand cell_density_by_band.csv', ha='center', va='center')
        ax6.set_title('Density Data')
        ax6.axis('off')
    
    plt.suptitle('Comprehensive Analysis Results', fontsize=16)
    plt.tight_layout()
    
    fig_path_png = os.path.join(save_dir, 'comprehensive_analysis.png')
    fig_path_pdf = os.path.join(save_dir, 'comprehensive_analysis.pdf')
    fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    print(f"保存综合分析图: {fig_path_png}")
    plt.show(); plt.close(fig)

def save_adata_safely(adata, filename, save_dir='make_band_new_2'):
    """安全保存adata，避免不规则数组错误"""
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
    for key in set(keys_to_remove):
        print(f"移除不规则的uns键: {key}")
        del adata_save.uns[key]
    
    adata_save.write(filepath)
    print(f"数据已安全保存: '{filepath}'")
    
    csv_path = os.path.join(save_dir, 'sample_band_summary.csv')
    summary_cols = ['mucus_zone', 'site_band', 'boundary_distance', 'is_boundary', 'boundary_method']
    sample_col_in_data = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col_in_data = col
            summary_cols.append(col)
            break
    existing_cols = [col for col in summary_cols if col in adata.obs.columns]
    summary_df = adata.obs[existing_cols].copy()
    summary_df.to_csv(csv_path, index=False)
    print(f"数据摘要已保存: '{csv_path}'")

def save_statistics_summary(adata, sample_col, save_dir):
    """保存统计摘要，包含条带面积和细胞密度信息"""
    save_dir = create_output_folder(save_dir)
    summary_path = os.path.join(save_dir, 'analysis_summary.txt')
    
    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("         Mucus Zone Boundary Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        
        # 总体统计
        f.write("1. OVERALL STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total cells: {len(adata)}\n")
        if 'mucus_zone' in adata.obs.columns:
            zone_counts = adata.obs['mucus_zone'].value_counts()
            f.write("\nMucus Zone Distribution:\n")
            for zone, count in zone_counts.items():
                f.write(f"  {zone}: {count} ({count/len(adata)*100:.1f}%)\n")
        if 'site_band' in adata.obs.columns:
            band_counts = adata.obs['site_band'].value_counts()
            f.write("\nSite Band Distribution:\n")
            for band, count in band_counts.items():
                if pd.notna(band):
                    f.write(f"  {band}: {count} ({count/len(adata)*100:.1f}%)\n")
        if 'is_boundary' in adata.obs.columns:
            boundary_count = adata.obs['is_boundary'].sum()
            f.write(f"\nBoundary cells: {boundary_count} ({boundary_count/len(adata)*100:.1f}%)\n")
        if 'boundary_method' in adata.obs.columns:
            method_counts = adata.obs['boundary_method'].value_counts()
            f.write("\nBoundary Calculation Methods:\n")
            for method, count in method_counts.items():
                f.write(f"  {method}: {count} ({count/len(adata)*100:.1f}%)\n")
        
        # 2. 面积统计
        f.write("\n\n2. BAND AREA STATISTICS\n")
        f.write("-" * 40 + "\n")
        if 'band_areas' in adata.uns:
            total_area_um2 = {}
            for sample, areas in adata.uns['band_areas'].items():
                for band, area in areas.items():
                    total_area_um2[band] = total_area_um2.get(band, 0) + area
            f.write("Overall Band Area:\n")
            for band, area_um2 in total_area_um2.items():
                area_mm2 = area_um2 / 1e6
                f.write(f"  {band}: {area_um2:.2f} μm² ({area_mm2:.4f} mm²)\n")
            f.write(f"\nTotal area (all bands): {sum(total_area_um2.values())/1e6:.2f} mm²\n")
        
        # 3. 细胞密度统计
        f.write("\n\n3. CELL DENSITY STATISTICS\n")
        f.write("-" * 40 + "\n")
        if 'cell_density' in adata.uns:
            density_df = adata.uns['cell_density']
            # 总体平均密度（按条带、细胞类型）
            avg_density = density_df.groupby(['band', 'cell_type'])['density_mm2'].mean().reset_index()
            f.write("Average density per band and cell type (cells/mm²):\n")
            for _, row in avg_density.iterrows():
                f.write(f"  {row['band']} - {row['cell_type']}: {row['density_mm2']:.2f}\n")
            
            # 也可以按样本列出主要密度
            f.write("\nSample-wise detailed density saved in 'cell_density_by_band.csv'.\n")
        else:
            f.write("  No cell density data available.\n")
        
        # 4. 按样本统计
        if sample_col in adata.obs.columns:
            f.write("\n\n4. SAMPLE-WISE STATISTICS\n")
            f.write("-" * 40 + "\n")
            samples = adata.obs[sample_col].unique()
            for sample in samples:
                f.write(f"\nSample: {sample}\n")
                sample_mask = adata.obs[sample_col] == sample
                sample_cells = adata[sample_mask]
                f.write(f"  Total cells: {len(sample_cells)}\n")
                if 'mucus_zone' in sample_cells.obs.columns:
                    zone_counts = sample_cells.obs['mucus_zone'].value_counts()
                    for zone, count in zone_counts.items():
                        f.write(f"    {zone}: {count} ({count/len(sample_cells)*100:.1f}%)\n")
                if 'site_band' in sample_cells.obs.columns:
                    band_counts = sample_cells.obs['site_band'].value_counts()
                    for band, count in band_counts.items():
                        if pd.notna(band):
                            f.write(f"    {band}: {count} ({count/len(sample_cells)*100:.1f}%)\n")
                if 'band_areas' in adata.uns and sample in adata.uns['band_areas']:
                    f.write(f"    Band Areas (mm²):\n")
                    for band, area_um2 in adata.uns['band_areas'][sample].items():
                        f.write(f"      {band}: {area_um2/1e6:.4f} mm²\n")
                if 'cell_density' in adata.uns:
                    sample_density = density_df[density_df['sample'] == sample]
                    if not sample_density.empty:
                        f.write("    Top densities per band:\n")
                        top_density = sample_density.sort_values('density_mm2', ascending=False).head(5)
                        for _, row in top_density.iterrows():
                            f.write(f"      {row['band']} - {row['cell_type']}: {row['density_mm2']:.2f} cells/mm²\n")
        
        # 5. 距离统计
        if 'boundary_distance' in adata.obs.columns:
            distances = adata.obs['boundary_distance'].dropna()
            if len(distances) > 0:
                f.write("\n\n5. DISTANCE STATISTICS\n")
                f.write("-" * 40 + "\n")
                f.write(f"Min: {distances.min():.2f} μm\n")
                f.write(f"Max: {distances.max():.2f} μm\n")
                f.write(f"Mean: {distances.mean():.2f} μm\n")
                f.write(f"Std: {distances.std():.2f} μm\n")
                f.write(f"Median: {distances.median():.2f} μm\n")
                f.write("Quantiles:\n")
                for q in [0.1, 0.25, 0.5, 0.75, 0.9]:
                    f.write(f"  {q*100:.0f}%: {distances.quantile(q):.2f} μm\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Analysis completed successfully!\n")
        f.write("=" * 60 + "\n")
    print(f"分析摘要已保存: {summary_path}")

# ==================== 修改主流程，集成密度计算和箱线图绘制 ====================
def process_mucus_zones_by_sample(adata, sample_col='sample', save_dir='make_band_new_2',
                                 on_band_width=50, out_band_width=100, in_band_width=100,
                                 area_method='cluster_convex',
                                 cell_type_col=None):
    """一站式处理：边界识别 + 条带标签 + 面积计算 + 细胞密度计算 + 箱线图 + 可视化 + 保存"""
    save_dir = create_output_folder(save_dir)
    
    # 确定细胞类型列
    if cell_type_col is None:
        # 尝试常见列名
        possible_cols = ['annotation', 'cell_type', 'Cell type', 'celltype', 'cluster']
        for col in possible_cols:
            if col in adata.obs.columns:
                cell_type_col = col
                print(f"自动识别细胞类型列: '{cell_type_col}'")
                break
        if cell_type_col is None:
            print("警告: 未找到细胞类型列，跳过细胞密度计算和箱线图绘制。")
    else:
        if cell_type_col not in adata.obs.columns:
            print(f"警告: 指定的细胞类型列 '{cell_type_col}' 不存在，跳过密度计算和箱线图绘制。")
            cell_type_col = None
    
    print("步骤1: 按样本识别inner-outer分界线...")
    adata = calculate_boundary_by_sample(
        adata, obs_col='mucus_zone', sample_col=sample_col,
        cat_a='inner', cat_b='outer', n_neighbors=10, cat_focus='inner',
        smooth_sigma=3, min_vertex_count=20, max_dist_threshold=500, save_dir=save_dir
    )
    
    print("\n步骤2: 按样本创建条带标签并计算面积...")
    adata = create_mucus_bands_by_sample(
        adata, sample_col=sample_col,
        on_band_width=on_band_width, out_band_width=out_band_width, in_band_width=in_band_width,
        area_method=area_method,
        save_dir=save_dir
    )
    
    # 新增步骤：计算细胞密度
    if cell_type_col is not None:
        print("\n步骤3: 计算各条带细胞密度...")
        calculate_cell_density_by_band(adata, sample_col, cell_type_col, save_dir)
        
        print("\n步骤3.5: 绘制细胞类型密度箱线图...")
        plot_cell_type_boxplots(adata, save_dir, metric='density_mm2')
    else:
        print("\n步骤3: 跳过细胞密度计算和箱线图绘制。")
    
    print("\n步骤4: 创建综合分析图...")
    create_comprehensive_visualization(
        adata, sample_col=sample_col, save_dir=save_dir,
        on_band_width=on_band_width, out_band_width=out_band_width, in_band_width=in_band_width
    )
    
    print("\n步骤5: 保存数据...")
    save_adata_safely(adata, 'adata_with_band_by_sample.h5ad', save_dir=save_dir)
    
    # 显示最终obs列
    print("\n最终的adata.obs列:")
    summary_cols = ['mucus_zone', 'site_band', 'boundary_distance', 'is_boundary', 'boundary_method']
    if sample_col in adata.obs.columns:
        summary_cols.append(sample_col)
    if cell_type_col is not None:
        summary_cols.append(cell_type_col)
    existing_cols = [col for col in summary_cols if col in adata.obs.columns]
    print(adata.obs[existing_cols].head())
    
    # 保存统计摘要（包含密度）
    save_statistics_summary(adata, sample_col, save_dir)
    
    return adata

def main():
    """主程序入口"""
    adata = sc.read_h5ad('extracted_samples.h5ad')
    
    if 'mucus_zone' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'mucus_zone'列")
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm中未找到'spatial'坐标信息")
    
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
    
    # 按样本处理，指定细胞类型列（例如 'annotation'）
    adata_result = process_mucus_zones_by_sample(
        adata,
        sample_col=sample_col,
        on_band_width=50,
        out_band_width=100,
        in_band_width=100,
        area_method='cluster_convex',
        cell_type_col='annotation'   # 请根据实际列名修改
    )
    return adata_result

if __name__ == "__main__":
    adata_result = main()
