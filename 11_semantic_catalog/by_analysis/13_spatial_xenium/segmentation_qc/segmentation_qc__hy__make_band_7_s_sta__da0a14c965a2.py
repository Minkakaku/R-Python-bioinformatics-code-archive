import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree, ConvexHull
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
from sklearn.cluster import DBSCAN
import os
import warnings
import pickle

warnings.filterwarnings('ignore')

def create_output_folder(folder_name='make_band_new_7'):
    """创建输出文件夹"""
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

# -------------------- 原有边界识别、条带创建函数（保持不变） --------------------
def calculate_boundary_by_sample(adata, obs_col, cat_a, cat_b, sample_col='sample',
                               n_neighbors=10, cat_focus='inner', 
                               smooth_sigma=3,
                               min_vertex_count=20,
                               max_dist_threshold=500,
                               save_dir='make_band_new_7'):
    """按样本独立计算inner和outer的分界线（同原脚本）"""
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
            print(f"  距离计算完成：范围 [{signed_distances.min():.2f}, {signed_distances.max():.2f}]，使用方法: {method_used}")
        else:
            print(f"  警告: 样本 {sample} 无法找到任何有效的分界线")
            sample_indices = adata.obs[sample_col] == sample
            adata.obs.loc[sample_indices, 'is_boundary'] = boundary_mask
            adata.obs.loc[sample_indices, 'boundary_method'] = 'none'
        
        # 绘制当前样本边界图（省略部分代码，保持原样）
        fig, ax = plt.subplots(figsize=(8, 8))
        for zone, color in [(cat_a, 'blue'), (cat_b, 'red')]:
            mask = sample_adata.obs[obs_col] == zone
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], c=color, s=5, alpha=0.6, label=zone)
        if len(valid_boundary_coords) > 0:
            for coords_segment in valid_boundary_coords:
                if method_used == "contour":
                    ls, col = '-', 'black'
                elif method_used == "convex_hull":
                    ls, col = '--', 'green'
                else:
                    ls, col = ':', 'orange'
                ax.plot(coords_segment[:, 0], coords_segment[:, 1], ls, color=col, linewidth=2, label=f'Boundary ({method_used})')
        ax.set_title(f'Sample: {sample} - Boundary Fitting\nMethod: {method_used}')
        ax.set_xlabel('X'); ax.set_ylabel('Y')
        if ax.get_legend_handles_labels()[1]:
            ax.legend(loc='upper right')
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_boundary.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_boundary.pdf'), bbox_inches='tight')
        plt.show(); plt.close(fig)
    
    print(f"\n处理完成统计: 成功 {len(processed_samples)} 个, 跳过 {len(skipped_samples)} 个")
    boundary_file = os.path.join(save_dir, 'sample_boundaries.pkl')
    with open(boundary_file, 'wb') as f:
        pickle.dump(sample_boundaries, f)
    print(f"样本边界坐标已保存为: {boundary_file}")
    return adata

def create_mucus_bands_by_sample(adata, sample_col='sample', 
                                on_band_width=50, 
                                in_a_band_width=100, in_b_band_width=100, in_c_band_width=100, in_d_band_width=100,
                                out_a_band_width=100, out_b_band_width=100, out_c_band_width=100, out_d_band_width=100,
                                save_dir='make_band_new_7'):
    """按样本独立创建条带标签（9个条带：on, in_a, in_b, in_c, in_d, out_a, out_b, out_c, out_d）"""
    if 'boundary_distance' not in adata.obs:
        raise ValueError("请先运行 calculate_boundary_by_sample 函数计算分界线距离")
    
    save_dir = create_output_folder(save_dir)
    adata.obs['site_band'] = np.nan
    samples = adata.obs[sample_col].unique()
    
    processed_samples = []
    skipped_samples = []
    
    # 预先计算各条带的边界距离阈值
    in_a_boundary = on_band_width + in_a_band_width
    in_b_boundary = in_a_boundary + in_b_band_width
    in_c_boundary = in_b_boundary + in_c_band_width
    in_d_boundary = in_c_boundary + in_d_band_width
    
    out_a_boundary = on_band_width + out_a_band_width
    out_b_boundary = out_a_boundary + out_b_band_width
    out_c_boundary = out_b_boundary + out_c_band_width
    out_d_boundary = out_c_boundary + out_d_band_width
    
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
        
        site_band_list = []
        for i, dist in enumerate(sample_distances):
            zone = sample_zones[i]
            if pd.isna(dist):
                site_band_list.append(zone)
            else:
                if abs(dist) <= on_band_width:
                    site_band_list.append('on')
                elif dist > on_band_width and dist <= in_a_boundary:
                    site_band_list.append('in_a')
                elif dist > in_a_boundary and dist <= in_b_boundary:
                    site_band_list.append('in_b')
                elif dist > in_b_boundary and dist <= in_c_boundary:
                    site_band_list.append('in_c')
                elif dist > in_c_boundary and dist <= in_d_boundary:
                    site_band_list.append('in_d')
                elif dist < -on_band_width and dist >= -out_a_boundary:
                    site_band_list.append('out_a')
                elif dist < -out_a_boundary and dist >= -out_b_boundary:
                    site_band_list.append('out_b')
                elif dist < -out_b_boundary and dist >= -out_c_boundary:
                    site_band_list.append('out_c')
                elif dist < -out_c_boundary and dist >= -out_d_boundary:
                    site_band_list.append('out_d')
                else:
                    site_band_list.append(zone)
        
        adata.obs.loc[sample_mask, 'site_band'] = site_band_list
        
        # 绘制当前样本的条带图（省略，保持原样）
        sample_adata = adata[sample_mask]
        coords = sample_adata.obsm['spatial']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        colors_zone = {'inner': 'blue', 'outer': 'red'}
        for zone in ['inner', 'outer']:
            mask = sample_adata.obs['mucus_zone'] == zone
            if mask.any():
                ax1.scatter(coords[mask, 0], coords[mask, 1], c=colors_zone[zone], s=5, alpha=0.6, label=zone)
        ax1.set_title(f'Sample: {sample} - Original Zones')
        ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.legend()
        
        colors_band = {'on': 'green', 'in_a': 'purple', 'in_b': 'violet', 'in_c': 'magenta', 'in_d': 'pink',
                      'out_a': 'orange', 'out_b': 'gold', 'out_c': 'darkorange', 'out_d': 'peru',
                      'inner': 'lightblue', 'outer': 'lightcoral'}
        bands = sample_adata.obs['site_band']
        if bands.notna().any():
            unique_bands = bands.dropna().unique()
            for band in unique_bands:
                mask = bands == band
                if mask.any():
                    ax2.scatter(coords[mask, 0], coords[mask, 1], c=colors_band.get(band, 'gray'), s=5, alpha=0.6, label=band)
        ax2.set_title(f'Sample: {sample} - Band Zones')
        ax2.set_xlabel('X'); ax2.set_ylabel('Y')
        handles, labels = ax2.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        if by_label:
            ax2.legend(by_label.values(), by_label.keys(), loc='upper right')
        
        plt.suptitle(f'Sample: {sample} - Zone Analysis', fontsize=16)
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_bands.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(save_dir, f'sample_{sample}_bands.pdf'), bbox_inches='tight')
        plt.show(); plt.close(fig)
    
    print(f"\n条带创建完成统计: 成功 {len(processed_samples)} 个, 跳过 {len(skipped_samples)} 个")
    return adata

# ==================== 新增：面积计算函数 ====================
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

def compute_band_areas(adata, sample_col='sample', area_method='cluster_convex', save_dir=None):
    """
    计算每个样本、每个条带的面积，并存储到 adata.uns['band_areas'] 中。
    同时导出 CSV 文件。
    """
    if 'site_band' not in adata.obs.columns:
        raise ValueError("请先运行 create_mucus_bands_by_sample 创建条带标签")
    
    save_dir = create_output_folder(save_dir) if save_dir else create_output_folder()
    
    if 'band_areas' not in adata.uns:
        adata.uns['band_areas'] = {}
    
    samples = adata.obs[sample_col].unique()
    all_records = []
    
    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask].copy()
        coords = sample_adata.obsm['spatial']
        bands_in_sample = sample_adata.obs['site_band'].dropna().unique()
        band_areas = {}
        
        print(f"\n计算样本 {sample} 各条带面积（方法: {area_method}）:")
        for band in bands_in_sample:
            mask_band = sample_adata.obs['site_band'] == band
            band_coords = coords[mask_band]
            if len(band_coords) >= 3:
                area = compute_area_from_coords(band_coords, method=area_method, eps=100, min_samples=3)
            else:
                area = 0.0
            band_areas[band] = area
            area_mm2 = area / 1e6
            print(f"    {band}: {area:.2f} μm² ({area_mm2:.4f} mm²)")
            all_records.append({'sample': sample, 'band': band, 'area_um2': area, 'area_mm2': area_mm2})
        
        adata.uns['band_areas'][sample] = band_areas
    
    # 导出面积CSV
    area_df = pd.DataFrame(all_records)
    csv_path = os.path.join(save_dir, 'band_areas.csv')
    area_df.to_csv(csv_path, index=False)
    print(f"\n条带面积统计已保存为: {csv_path}")
    
    return adata

# ==================== 新增：细胞密度计算 ====================
def calculate_cell_density_by_band(adata, sample_col='sample', cell_type_col='annotation', save_dir=None):
    """
    基于条带面积和细胞类型计数，计算每个样本、每个条带、每种细胞类型的密度。
    密度单位：细胞数/mm²（1 mm² = 1e6 μm²）。
    结果存储在 adata.uns['cell_density'] 中，并导出CSV。
    """
    if cell_type_col not in adata.obs.columns:
        print(f"警告: 细胞类型列 '{cell_type_col}' 不存在，跳过密度计算。")
        return adata
    if 'band_areas' not in adata.uns:
        print("警告: 未找到条带面积信息，请先运行 compute_band_areas。")
        return adata
    
    save_dir = create_output_folder(save_dir) if save_dir else create_output_folder()
    
    density_records = []
    samples = adata.obs[sample_col].unique()
    
    for sample in samples:
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask]
        
        if sample not in adata.uns['band_areas']:
            print(f"警告: 样本 {sample} 缺少面积信息，跳过该样本密度计算。")
            continue
        band_areas = adata.uns['band_areas'][sample]
        
        grouped = sample_adata.obs.groupby(['site_band', cell_type_col]).size().reset_index(name='count')
        
        for _, row in grouped.iterrows():
            band = row['site_band']
            cell_type = row[cell_type_col]
            count = row['count']
            area_um2 = band_areas.get(band, np.nan)
            if pd.isna(area_um2) or area_um2 == 0:
                continue
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
        return adata
    
    density_df = pd.DataFrame(density_records)
    adata.uns['cell_density'] = density_df
    
    csv_path = os.path.join(save_dir, 'cell_density_by_band.csv')
    density_df.to_csv(csv_path, index=False)
    print(f"细胞密度统计已保存: {csv_path}")
    
    return adata

# ==================== 修改后的密度箱线图（按指定顺序排列条带）====================
def plot_cell_type_density_individual(adata, save_dir, metric='density_mm2'):
    """
    为每个细胞类型绘制单独的密度箱线图，比较不同条带间的指标。
    条带按照指定顺序排列：inner, in_d, in_c, in_b, in_a, on, out_a, out_b, out_c, out_d, outer。
    每张图保存为 PDF 和 PNG。
    """
    if 'cell_density' not in adata.uns:
        print("未找到细胞密度数据，请先运行 calculate_cell_density_by_band")
        return
    
    density_df = adata.uns['cell_density']
    required = ['cell_type', 'band', metric]
    if not all(col in density_df.columns for col in required):
        print(f"密度数据框缺少列: {required}")
        return
    
    # 定义固定的条带顺序
    band_order = ['inner', 'in_d', 'in_c', 'in_b', 'in_a', 'on', 
                  'out_a', 'out_b', 'out_c', 'out_d', 'outer']
    
    cell_types = density_df['cell_type'].unique()
    if len(cell_types) == 0:
        print("无细胞类型数据")
        return
    
    # 为每个细胞类型单独绘图
    for ct in cell_types:
        ct_data = density_df[density_df['cell_type'] == ct]
        # 获取该细胞类型存在的条带
        present_bands = ct_data['band'].unique()
        # 按照 band_order 过滤并排序
        ordered_bands = [b for b in band_order if b in present_bands]
        
        if len(ordered_bands) == 0:
            print(f"细胞类型 {ct} 没有条带数据，跳过")
            continue
        
        data_to_plot = [ct_data[ct_data['band'] == b][metric].values for b in ordered_bands]
        
        fig, ax = plt.subplots(figsize=(8, 6))
        bp = ax.boxplot(data_to_plot, labels=ordered_bands, patch_artist=True)
        # 为盒子设置不同颜色
        colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow', 'lightpink',
                  'lightcyan', 'lightgray', 'lightsalmon', 'lightseagreen', 'lightsteelblue', 'lightgoldenrodyellow']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        ax.set_title(f'Cell Type: {ct} - {metric} by Band')
        ax.set_ylabel(metric)
        ax.set_xlabel('Band')
        # 旋转标签以防重叠
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        # 保存
        safe_ct = ct.replace('/', '_').replace(' ', '_')
        fig_path_png = os.path.join(save_dir, f'density_{metric}_{safe_ct}.png')
        fig_path_pdf = os.path.join(save_dir, f'density_{metric}_{safe_ct}.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        print(f"保存细胞类型 {ct} 密度图: {fig_path_png}")
        plt.show()
        plt.close(fig)

# ==================== 原有综合可视化（保持不变） ====================
def create_comprehensive_visualization(adata, sample_col='sample', save_dir='make_band_new_7',
                                      on_band_width=50, 
                                      in_a_band_width=100, in_b_band_width=100, in_c_band_width=100, in_d_band_width=100,
                                      out_a_band_width=100, out_b_band_width=100, out_c_band_width=100, out_d_band_width=100):
    """创建综合可视化图形（同原脚本，未修改）"""
    save_dir = create_output_folder(save_dir)
    
    in_a_boundary = on_band_width + in_a_band_width
    in_b_boundary = in_a_boundary + in_b_band_width
    in_c_boundary = in_b_boundary + in_c_band_width
    in_d_boundary = in_c_boundary + in_d_band_width
    out_a_boundary = on_band_width + out_a_band_width
    out_b_boundary = out_a_boundary + out_b_band_width
    out_c_boundary = out_b_boundary + out_c_band_width
    out_d_boundary = out_c_boundary + out_d_band_width
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. 所有样本的条带分布统计
    ax1 = axes[0, 0]
    if 'site_band' in adata.obs.columns:
        band_counts = adata.obs['site_band'].value_counts()
        colors = {'on': 'green', 
                 'in_a': 'purple', 'in_b': 'violet', 'in_c': 'magenta', 'in_d': 'pink',
                 'out_a': 'orange', 'out_b': 'gold', 'out_c': 'darkorange', 'out_d': 'peru',
                 'inner': 'lightblue', 'outer': 'lightcoral'}
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
        band_colors = {'on': 'green', 
                      'in_a': 'purple', 'in_b': 'violet', 'in_c': 'magenta', 'in_d': 'pink',
                      'out_a': 'orange', 'out_b': 'gold', 'out_c': 'darkorange', 'out_d': 'peru',
                      'inner': 'lightblue', 'outer': 'lightcoral'}
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
                    ax2.bar(range(len(df)), df[band], bottom=bottom, label=band, color=band_colors.get(band, 'gray'))
                    bottom += df[band].values
            ax2.set_title('Band Distribution by Sample')
            ax2.set_xlabel('Sample')
            ax2.set_ylabel('Cell Count')
            ax2.set_xticks(range(len(sample_names)))
            ax2.set_xticklabels(sample_names, rotation=45, ha='right')
            if len(all_bands) <= 12:
                ax2.legend(loc='upper right')
    
    # 3. 距离分布直方图
    ax3 = axes[1, 0]
    if 'boundary_distance' in adata.obs.columns:
        distances = adata.obs['boundary_distance'].dropna()
        if len(distances) > 0:
            ax3.hist(distances, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Boundary (0)')
            ax3.axvline(x=-on_band_width, color='orange', linestyle=':', linewidth=1, label=f'Outer edge (-{on_band_width})')
            ax3.axvline(x=on_band_width, color='orange', linestyle=':', linewidth=1, label=f'Inner edge (+{on_band_width})')
            ax3.axvline(x=in_a_boundary, color='purple', linestyle=':', linewidth=0.8, label='Inner a edge')
            ax3.axvline(x=in_b_boundary, color='violet', linestyle=':', linewidth=0.8, label='Inner b edge')
            ax3.axvline(x=in_c_boundary, color='magenta', linestyle=':', linewidth=0.8, label='Inner c edge')
            ax3.axvline(x=in_d_boundary, color='pink', linestyle=':', linewidth=0.8, label='Inner d edge')
            ax3.axvline(x=-out_a_boundary, color='gold', linestyle=':', linewidth=0.8, label='Outer a edge')
            ax3.axvline(x=-out_b_boundary, color='darkorange', linestyle=':', linewidth=0.8, label='Outer b edge')
            ax3.axvline(x=-out_c_boundary, color='peru', linestyle=':', linewidth=0.8, label='Outer c edge')
            ax3.axvline(x=-out_d_boundary, color='saddlebrown', linestyle=':', linewidth=0.8, label='Outer d edge')
            ax3.set_title('Distribution of Distance to Boundary')
            ax3.set_xlabel('Distance to Boundary (um)')
            ax3.set_ylabel('Frequency')
            ax3.legend(loc='upper right', fontsize='x-small')
    
    # 4. 边界计算方法统计
    ax4 = axes[1, 1]
    if 'boundary_method' in adata.obs.columns:
        method_counts = adata.obs['boundary_method'].value_counts()
        if len(method_counts) > 0:
            colors = {'contour': 'blue', 'convex_hull': 'green', 'bounding_box': 'orange', 'none': 'red'}
            bar_colors = [colors.get(method, 'gray') for method in method_counts.index]
            bars = ax4.bar(method_counts.index.astype(str), method_counts.values, color=bar_colors)
            ax4.set_title('Boundary Calculation Methods')
            ax4.set_xlabel('Method')
            ax4.set_ylabel('Number of Cells')
            ax4.set_xticklabels(method_counts.index, rotation=45, ha='right')
            for bar, count in zip(bars, method_counts.values):
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2., height + 5, f'{int(count)}', ha='center', va='bottom')
    
    plt.suptitle('Comprehensive Analysis Results (Extended Bands)', fontsize=16)
    plt.tight_layout()
    fig.savefig(os.path.join(save_dir, 'comprehensive_analysis.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(save_dir, 'comprehensive_analysis.pdf'), bbox_inches='tight')
    print(f"保存综合分析图: {os.path.join(save_dir, 'comprehensive_analysis.png')}")
    plt.show(); plt.close(fig)

# ==================== 数据保存函数（保持不变） ====================
def save_adata_safely(adata, filename, save_dir='make_band_new_7'):
    """安全保存adata，避免不规则数组错误（同原脚本）"""
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
        print(f"警告：移除不规则的uns键: {key}")
        if key in adata_save.uns:
            del adata_save.uns[key]
    adata_save.write(filepath)
    print(f"数据已安全保存为 '{filepath}'")
    
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
    print(f"数据摘要已保存为 '{csv_path}'")

def save_statistics_summary(adata, sample_col, save_dir):
    """保存统计摘要（同原脚本）"""
    save_dir = create_output_folder(save_dir)
    summary_path = os.path.join(save_dir, 'analysis_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("         Mucus Zone Boundary Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
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
        # 面积和密度信息可以加在最后（略）
        if 'band_areas' in adata.uns:
            f.write("\n2. BAND AREAS\n")
            f.write("-" * 40 + "\n")
            for sample, areas in adata.uns['band_areas'].items():
                f.write(f"Sample {sample}:\n")
                for band, area in areas.items():
                    f.write(f"  {band}: {area:.2f} μm² ({area/1e6:.4f} mm²)\n")
        if 'cell_density' in adata.uns:
            f.write("\n3. CELL DENSITY\n")
            f.write("-" * 40 + "\n")
            density_df = adata.uns['cell_density']
            avg_density = density_df.groupby(['band', 'cell_type'])['density_mm2'].mean().reset_index()
            f.write("Average density per band and cell type (cells/mm²):\n")
            for _, row in avg_density.iterrows():
                f.write(f"  {row['band']} - {row['cell_type']}: {row['density_mm2']:.2f}\n")
        # ... 其余统计与原脚本一致
    print(f"分析摘要已保存为: {summary_path}")

# ==================== 修改主流程，集成面积和密度计算 ====================
def process_mucus_zones_by_sample(adata, sample_col='sample', save_dir='make_band_new_7',
                                 on_band_width=50, 
                                 in_a_band_width=100, in_b_band_width=100, in_c_band_width=100, in_d_band_width=100,
                                 out_a_band_width=100, out_b_band_width=100, out_c_band_width=100, out_d_band_width=100,
                                 cell_type_col=None):
    """
    一站式处理：边界识别 + 条带标签 + 面积计算 + 细胞密度计算 + 单独细胞类型箱线图 + 综合可视化 + 保存
    """
    save_dir = create_output_folder(save_dir)
    
    # 确定细胞类型列
    if cell_type_col is None:
        possible_cols = ['annotation', 'cell_type', 'Cell type', 'celltype', 'cluster']
        for col in possible_cols:
            if col in adata.obs.columns:
                cell_type_col = col
                print(f"自动识别细胞类型列: '{cell_type_col}'")
                break
        if cell_type_col is None:
            print("警告: 未找到细胞类型列，将跳过细胞密度计算和箱线图绘制。")
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
    
    print("\n步骤2: 按样本创建条带标签...")
    adata = create_mucus_bands_by_sample(
        adata, sample_col=sample_col,
        on_band_width=on_band_width,
        in_a_band_width=in_a_band_width, in_b_band_width=in_b_band_width,
        in_c_band_width=in_c_band_width, in_d_band_width=in_d_band_width,
        out_a_band_width=out_a_band_width, out_b_band_width=out_b_band_width,
        out_c_band_width=out_c_band_width, out_d_band_width=out_d_band_width,
        save_dir=save_dir
    )
    
    print("\n步骤3: 计算各条带面积...")
    adata = compute_band_areas(adata, sample_col=sample_col, area_method='cluster_convex', save_dir=save_dir)
    
    if cell_type_col is not None:
        print("\n步骤4: 计算各条带细胞密度...")
        adata = calculate_cell_density_by_band(adata, sample_col=sample_col, cell_type_col=cell_type_col, save_dir=save_dir)
        
        print("\n步骤5: 为每个细胞类型绘制单独的密度箱线图...")
        plot_cell_type_density_individual(adata, save_dir, metric='density_mm2')
    else:
        print("\n步骤4: 跳过细胞密度计算和箱线图绘制。")
    
    print("\n步骤6: 创建综合分析图...")
    create_comprehensive_visualization(
        adata, sample_col=sample_col, save_dir=save_dir,
        on_band_width=on_band_width,
        in_a_band_width=in_a_band_width, in_b_band_width=in_b_band_width,
        in_c_band_width=in_c_band_width, in_d_band_width=in_d_band_width,
        out_a_band_width=out_a_band_width, out_b_band_width=out_b_band_width,
        out_c_band_width=out_c_band_width, out_d_band_width=out_d_band_width
    )
    
    print("\n步骤7: 保存数据...")
    save_adata_safely(adata, 'adata_with_band_by_sample.h5ad', save_dir=save_dir)
    
    # 保存统计摘要（已包含面积和密度信息）
    save_statistics_summary(adata, sample_col, save_dir)
    
    return adata

# ==================== 主程序 ====================
def main():
    # 加载数据
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
    
    # 按样本处理，指定细胞类型列（根据您的数据修改）
    adata_result = process_mucus_zones_by_sample(
        adata,
        sample_col=sample_col,
        on_band_width=50,
        in_a_band_width=100, in_b_band_width=100, in_c_band_width=100, in_d_band_width=100,
        out_a_band_width=100, out_b_band_width=100, out_c_band_width=100, out_d_band_width=100,
        cell_type_col='annotation'   # 请根据实际列名修改，例如 'cell_type'
    )
    return adata_result

if __name__ == "__main__":
    adata_result = main()
