import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier

def calculate_boundary(adata, obs_col, cat_a, cat_b, n_neighbors=10, cat_focus='inner', 
                     smooth_sigma=2, min_vertex_count=200, max_dist_threshold=300):
    """
    基于KNN方法识别inner和outer的分界线
    """
    # 1. 获取坐标和标签
    coords = adata.obsm['spatial']
    labels = adata.obs[obs_col].values
    
    # 2. 构建空间邻居图识别边界细胞
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    boundary_mask = np.zeros(len(adata), dtype=bool)
    
    # 3. 遍历每个细胞，检查其邻居
    for i in range(len(adata)):
        current_label = labels[i]
        if current_label not in [cat_a, cat_b]:
            continue
            
        neighbor_indices = indices[i]
        neighbor_labels = labels[neighbor_indices]
        
        # 检查邻居中是否存在"对立"的类别
        if current_label == cat_a and cat_b in neighbor_labels:
            boundary_mask[i] = True
        elif current_label == cat_b and cat_a in neighbor_labels:
            boundary_mask[i] = True
    
    adata.obs['is_boundary'] = boundary_mask
    print(f"识别出 {sum(boundary_mask)} 个边界细胞")
    
    # 4. 训练KNN分类器拟合分界线
    X = adata.obsm['spatial']
    y = (adata.obs[obs_col] == cat_focus).astype(int).values
    clf = KNeighborsClassifier(n_neighbors=n_neighbors)
    clf.fit(X, y)
    
    # 5. 创建网格进行概率预测
    padding_percent = 0.1
    x_min, x_max = X[:, 0].min(), X[:, 0].max()
    y_min, y_max = X[:, 1].min(), X[:, 1].max()
    x_pad = (x_max - x_min) * padding_percent
    y_pad = (y_max - y_min) * padding_percent
    
    grid_resolution = 500
    xx, yy = np.meshgrid(np.linspace(x_min - x_pad, x_max + x_pad, grid_resolution),
                         np.linspace(y_min - y_pad, y_max + y_pad, grid_resolution))
    grid_points = np.c_[xx.ravel(), yy.ravel()]
    
    # 计算网格点距离最近真实细胞的距离
    nbrs_geometry = NearestNeighbors(n_neighbors=1).fit(X)
    dists, _ = nbrs_geometry.kneighbors(grid_points)
    dists = dists.reshape(xx.shape)
    
    # 6. 预测概率 + 高斯模糊
    Z = clf.predict_proba(grid_points)[:, 1]
    Z = Z.reshape(xx.shape)
    Z_smooth = gaussian_filter(Z, sigma=smooth_sigma)
    Z_smooth[dists > max_dist_threshold] = np.nan
    
    # 7. 提取分界线（概率为0.5的等高线）
    fig, ax = plt.subplots(figsize=(8, 8))
    contours = ax.contour(xx, yy, Z_smooth, levels=[0.5], alpha=0)
    
    # 收集有效的分界线坐标
    valid_boundary_coords = []
    if len(contours.allsegs) > 0:
        segments = contours.allsegs[0]
        for coords_segment in segments:
            if len(coords_segment) > min_vertex_count:
                ax.plot(coords_segment[:, 0], coords_segment[:, 1], color='black', linewidth=2)
                valid_boundary_coords.append(coords_segment)
    
    # 8. 计算每个细胞到分界线的距离（带符号）
    if len(valid_boundary_coords) > 0:
        all_boundary_points = np.vstack(valid_boundary_coords)
        tree = cKDTree(all_boundary_points)
        dist_values, _ = tree.query(X, k=1)
        
        # 给距离加上方向：inner侧为正，outer侧为负
        signed_distances = dist_values.copy()
        outer_mask = (adata.obs[obs_col] == cat_b).values
        signed_distances[outer_mask] = -signed_distances[outer_mask]
        
        adata.obs['boundary_distance'] = signed_distances
        
        # 修复：不将边界坐标保存到uns中，避免h5ad保存错误
        # adata.uns['fitted_boundary_coords'] = valid_boundary_coords
        
        print(f"距离计算完成：范围 [{signed_distances.min():.2f}, {signed_distances.max():.2f}]")
    else:
        print("警告：未找到有效的分界线")
        adata.obs['boundary_distance'] = np.nan
    
    # 可视化结果
    plt.scatter(X[:, 0], X[:, 1], c=adata.obs[obs_col].map({cat_a: 'blue', cat_b: 'red'}), 
                s=1, alpha=0.6)
    if len(valid_boundary_coords) > 0:
        for coords_segment in valid_boundary_coords:
            plt.plot(coords_segment[:, 0], coords_segment[:, 1], 'k-', linewidth=2, label='Fitted Boundary')
    plt.title(f'{cat_a}-{cat_b} Boundary Fitting')
    plt.legend()
    plt.show()
    
    return adata

def create_mucus_bands(adata, on_band_width=100, out_band_width=200, in_band_width=200):
    """
    根据分界线创建条带标签
    """
    if 'boundary_distance' not in adata.obs:
        raise ValueError("请先运行 calculate_boundary 函数计算分界线距离")
    
    distances = adata.obs['boundary_distance'].values
    
    # 修复Categorical类型错误：创建新的字符串列表
    site_band_list = []
    
    # 将原始区域转换为字符串类型
    if hasattr(adata.obs['mucus_zone'], 'cat'):
        # 如果是Categorical类型
        original_zones = adata.obs['mucus_zone'].astype(str).values
    else:
        original_zones = adata.obs['mucus_zone'].values
    
    for i, dist in enumerate(distances):
        zone = original_zones[i]
        
        # 1. 中心条带：分界线两侧各100um，总宽度200um
        if abs(dist) <= on_band_width:
            site_band_list.append('on')
        # 2. 外侧条带：在中心条带外向outer方向推进200um
        elif dist < -on_band_width and dist >= -(on_band_width + out_band_width):
            site_band_list.append('out')
        # 3. 内侧条带：在中心条带内向inner方向推进200um
        elif dist > on_band_width and dist <= (on_band_width + in_band_width):
            site_band_list.append('in')
        # 4. 其他细胞保持原mucus_zone标签
        else:
            site_band_list.append(zone)
    
    # 创建新的字符串列，避免Categorical问题
    adata.obs['site_band'] = site_band_list
    
    # 统计各条带细胞数量
    band_counts = pd.Series(site_band_list).value_counts()
    print("条带细胞分布统计:")
    for band, count in band_counts.items():
        print(f"  {band}: {count} 个细胞")
    
    return adata

def visualize_bands(adata):
    """
    可视化条带分布
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # 原始mucus_zone分布
    colors_zone = {'inner': 'blue', 'outer': 'red'}
    coords = adata.obsm['spatial']
    
    # 确保使用字符串类型
    if hasattr(adata.obs['mucus_zone'], 'cat'):
        zones = adata.obs['mucus_zone'].astype(str)
    else:
        zones = adata.obs['mucus_zone']
    
    for zone in ['inner', 'outer']:
        mask = zones == zone
        if mask.any():
            ax1.scatter(coords[mask, 0], coords[mask, 1], c=colors_zone[zone], 
                       s=1, alpha=0.6, label=zone)
    ax1.set_title('Original mucus_zone')
    ax1.legend()
    
    # 条带分布
    colors_band = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                   'inner': 'lightblue', 'outer': 'lightcoral'}
    bands = adata.obs['site_band']
    
    for band in np.unique(bands):
        mask = bands == band
        if mask.any():
            ax2.scatter(coords[mask, 0], coords[mask, 1], c=colors_band.get(band, 'gray'), 
                       s=1, alpha=0.6, label=band)
    ax2.set_title('Site Band Labels')
    ax2.legend()
    
    # 距离分布热图
    if 'boundary_distance' in adata.obs:
        scatter = ax3.scatter(coords[:, 0], coords[:, 1], c=adata.obs['boundary_distance'], 
                             cmap='coolwarm', s=1, alpha=0.8)
        ax3.set_title('Distance to Boundary (um)')
        plt.colorbar(scatter, ax=ax3)
    
    plt.tight_layout()
    plt.show()

def save_adata_safely(adata, filename):
    """
    安全保存adata，避免不规则数组导致的错误
    """
    # 创建一个adata的副本，避免修改原始数据
    adata_save = adata.copy()
    
    # 检查uns中是否有不规则的数据结构
    keys_to_remove = []
    for key in adata_save.uns.keys():
        try:
            # 尝试将数据转换为numpy数组来检查是否规则
            test_value = adata_save.uns[key]
            if isinstance(test_value, (list, np.ndarray)):
                # 如果是列表或数组，检查是否规则
                if isinstance(test_value, list) and len(test_value) > 0:
                    # 检查列表中的元素是否都是相同形状的数组
                    if all(isinstance(x, np.ndarray) for x in test_value):
                        shapes = [x.shape for x in test_value]
                        if len(set(shapes)) > 1:
                            # 形状不一致，需要移除
                            keys_to_remove.append(key)
        except:
            # 如果转换失败，也移除
            keys_to_remove.append(key)
    
    # 移除不规则的键
    for key in keys_to_remove:
        print(f"警告：移除不规则的uns键: {key}")
        del adata_save.uns[key]
    
    # 保存清理后的数据
    adata_save.write(filename)
    print(f"数据已安全保存为 '{filename}'")

# 主程序
def main():
    # 加载您的h5ad数据（请替换为您的文件路径）
    # adata = sc.read_h5ad('your_data.h5ad')
    
    # 假设您的数据已经加载为 adata
    # 检查数据是否包含必要的列
    if 'mucus_zone' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'mucus_zone'列")
    
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm中未找到'spatial'坐标信息")
    
    print("步骤1: 识别inner-outer分界线...")
    adata = calculate_boundary(
        adata, 
        obs_col='mucus_zone',
        cat_a='inner', 
        cat_b='outer',
        n_neighbors=10,
        cat_focus='inner',
        smooth_sigma=2,
        max_dist_threshold=300
    )
    
    print("\n步骤2: 创建条带标签...")
    adata = create_mucus_bands(
        adata,
        on_band_width=100,    # 中心条带宽度100um（两侧各100um）
        out_band_width=200,   # 外侧条带宽度200um
        in_band_width=200     # 内侧条带宽度200um
    )
    
    print("\n步骤3: 可视化结果...")
    visualize_bands(adata)
    
    print("\n步骤4: 保存数据...")
    # 使用安全保存函数
    save_adata_safely(adata, 'adata_with_band.h5ad')
    
    # 显示最终的数据结构
    print("\n最终的adata.obs列:")
    print(adata.obs[['mucus_zone', 'site_band', 'boundary_distance', 'is_boundary']].head())

# 简化使用方法：
def process_mucus_zones(adata):
    """
    一站式处理函数
    """
    # 1. 计算分界线
    adata = calculate_boundary(adata, 'mucus_zone', 'inner', 'outer')
    
    # 2. 创建条带
    adata = create_mucus_bands(adata)
    
    # 3. 使用安全保存函数
    save_adata_safely(adata, 'adata_with_band.h5ad')
    
    return adata

# 使用示例：
adata = sc.read_h5ad('merged_TMA_data.h5ad')
adata_result = process_mucus_zones(adata)
