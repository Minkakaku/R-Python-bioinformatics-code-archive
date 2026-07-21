import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
import os
import warnings
import pickle
warnings.filterwarnings('ignore')

def create_output_folder(folder_name='make_band_new'):
    """
    创建输出文件夹
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

def calculate_boundary_by_sample(adata, obs_col, cat_a, cat_b, sample_col='sample',
                               n_neighbors=10, cat_focus='inner', 
                               smooth_sigma=2, min_vertex_count=200, max_dist_threshold=300,
                               save_dir='make_band_new'):
    """
    按样本独立计算inner和outer的分界线，并为每个样本单独绘制边界图
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 确保有样本列
    if sample_col not in adata.obs.columns:
        raise ValueError(f"adata.obs中未找到'{sample_col}'列")
    
    # 初始化结果列
    adata.obs['is_boundary'] = False
    adata.obs['boundary_distance'] = np.nan
    
    # 获取所有样本
    samples = adata.obs[sample_col].unique()
    print(f"发现 {len(samples)} 个样本")
    
    # 存储每个样本的边界坐标
    sample_boundaries = {}
    
    # 统计处理结果
    processed_samples = []
    skipped_samples = []
    
    for idx, sample in enumerate(samples):
        # 获取当前样本的数据
        sample_mask = adata.obs[sample_col] == sample
        sample_adata = adata[sample_mask].copy()
        
        if len(sample_adata) == 0:
            continue
            
        print(f"\n处理样本: {sample} ({len(sample_adata)} 个细胞)")
        
        # 检查样本是否包含两种类别
        zones_in_sample = sample_adata.obs[obs_col].unique()
        if cat_a not in zones_in_sample or cat_b not in zones_in_sample:
            print(f"  警告: 样本 {sample} 不包含 {cat_a} 和 {cat_b} 两种区域，跳过")
            skipped_samples.append(sample)
            continue
            
        processed_samples.append(sample)
        
        # 1. 获取坐标和标签
        coords = sample_adata.obsm['spatial']
        labels = sample_adata.obs[obs_col].values
        
        # 2. 构建空间邻居图识别边界细胞
        nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        boundary_mask = np.zeros(len(sample_adata), dtype=bool)
        
        # 3. 遍历每个细胞，检查其邻居
        for i in range(len(sample_adata)):
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
        
        print(f"  识别出 {sum(boundary_mask)} 个边界细胞")
        
        # 4. 训练KNN分类器拟合分界线
        X = coords
        y = (sample_adata.obs[obs_col] == cat_focus).astype(int).values
        clf = KNeighborsClassifier(n_neighbors=n_neighbors)
        clf.fit(X, y)
        
        # 5. 创建网格进行概率预测（只使用当前样本的范围）
        padding_percent = 0.1
        x_min, x_max = X[:, 0].min(), X[:, 0].max()
        y_min, y_max = X[:, 1].min(), X[:, 1].max()
        x_pad = (x_max - x_min) * padding_percent
        y_pad = (y_max - y_min) * padding_percent
        
        grid_resolution = 300
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
        try:
            contours = plt.contour(xx, yy, Z_smooth, levels=[0.5], alpha=0)
            
            # 收集有效的分界线坐标
            valid_boundary_coords = []
            if len(contours.allsegs) > 0:
                segments = contours.allsegs[0]
                for coords_segment in segments:
                    if len(coords_segment) > min_vertex_count:
                        valid_boundary_coords.append(coords_segment)
            
            # 8. 计算每个细胞到分界线的距离（带符号）
            if len(valid_boundary_coords) > 0:
                all_boundary_points = np.vstack(valid_boundary_coords)
                tree = cKDTree(all_boundary_points)
                dist_values, _ = tree.query(X, k=1)
                
                # 给距离加上方向：inner侧为正，outer侧为负
                signed_distances = dist_values.copy()
                outer_mask = (sample_adata.obs[obs_col] == cat_b).values
                signed_distances[outer_mask] = -signed_distances[outer_mask]
                
                # 更新原始adata中的距离和边界标记
                sample_indices = adata.obs[sample_col] == sample
                adata.obs.loc[sample_indices, 'boundary_distance'] = signed_distances
                adata.obs.loc[sample_indices, 'is_boundary'] = boundary_mask
                
                # 存储样本边界坐标
                sample_boundaries[sample] = valid_boundary_coords
                
                print(f"  距离计算完成：范围 [{signed_distances.min():.2f}, {signed_distances.max():.2f}]")
            else:
                print(f"  警告: 样本 {sample} 未找到有效的分界线")
                valid_boundary_coords = []
                
        except Exception as e:
            print(f"  警告: 样本 {sample} 计算分界线时出错: {e}")
            valid_boundary_coords = []
        
        # 为当前样本单独绘制边界图
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # 绘制细胞点
        for zone, color in [(cat_a, 'blue'), (cat_b, 'red')]:
            mask = sample_adata.obs[obs_col] == zone
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], c=color, 
                          s=5, alpha=0.6, label=zone)
        
        # 绘制边界线
        if len(valid_boundary_coords) > 0:
            for coords_segment in valid_boundary_coords:
                ax.plot(coords_segment[:, 0], coords_segment[:, 1], 'k-', 
                       linewidth=2, label='Fitted Boundary')
        
        ax.set_title(f'Sample: {sample} - Boundary Fitting')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        
        # 添加图例
        handles, labels = ax.get_legend_handles_labels()
        if handles:  # 只有当有图例项时才添加图例
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())
        
        plt.tight_layout()
        
        # 保存当前样本的边界图
        fig_path_png = os.path.join(save_dir, f'sample_{sample}_boundary.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_boundary.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        print(f"  保存边界图: {fig_path_png}")
        print(f"  保存边界图: {fig_path_pdf}")
        
        plt.show()
        plt.close(fig)
    
    # 输出统计信息
    print(f"\n处理完成统计:")
    print(f"  成功处理的样本: {len(processed_samples)} 个")
    if processed_samples:
        print(f"  成功处理的样本列表: {', '.join(processed_samples)}")
    print(f"  跳过的样本: {len(skipped_samples)} 个")
    if skipped_samples:
        print(f"  跳过的样本列表: {', '.join(skipped_samples)}")
    
    # 将样本边界坐标单独保存为pickle文件
    boundary_file = os.path.join(save_dir, 'sample_boundaries.pkl')
    with open(boundary_file, 'wb') as f:
        pickle.dump(sample_boundaries, f)
    print(f"\n样本边界坐标已保存为: {boundary_file}")
    
    return adata

def create_mucus_bands_by_sample(adata, sample_col='sample', 
                                on_band_width=50, out_band_width=100, in_band_width=100,
                                save_dir='make_band_new'):
    """
    按样本独立创建条带标签，并为每个样本绘制条带区域图
    """
    if 'boundary_distance' not in adata.obs:
        raise ValueError("请先运行 calculate_boundary_by_sample 函数计算分界线距离")
    
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 初始化site_band列
    adata.obs['site_band'] = np.nan
    
    # 获取所有样本
    samples = adata.obs[sample_col].unique()
    
    # 统计
    processed_samples = []
    skipped_samples = []
    
    for sample in samples:
        # 获取当前样本的数据
        sample_mask = adata.obs[sample_col] == sample
        if not sample_mask.any():
            continue
            
        sample_distances = adata.obs.loc[sample_mask, 'boundary_distance'].values
        sample_zones = adata.obs.loc[sample_mask, 'mucus_zone'].values
        
        # 检查当前样本是否计算了距离
        if np.isnan(sample_distances).all():
            print(f"样本 {sample} 未计算距离，跳过")
            skipped_samples.append(sample)
            continue
            
        processed_samples.append(sample)
        
        # 创建当前样本的条带标签
        site_band_list = []
        
        for i, dist in enumerate(sample_distances):
            zone = sample_zones[i]
            
            if pd.isna(dist):
                # 如果距离为NaN，保持原始zone
                site_band_list.append(zone)
            else:
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
        
        # 更新原始adata中的条带标签
        adata.obs.loc[sample_mask, 'site_band'] = site_band_list
        
        # 统计当前样本的条带分布
        band_counts = pd.Series(site_band_list).value_counts()
        print(f"\n样本 {sample} 条带细胞分布统计:")
        for band, count in band_counts.items():
            print(f"  {band}: {count} 个细胞")
        
        # 为当前样本绘制条带区域图
        sample_adata = adata[sample_mask]
        coords = sample_adata.obsm['spatial']
        
        # 创建条带区域图
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # 左图：原始区域分布
        colors_zone = {'inner': 'blue', 'outer': 'red'}
        for zone in ['inner', 'outer']:
            mask = sample_adata.obs['mucus_zone'] == zone
            if mask.any():
                ax1.scatter(coords[mask, 0], coords[mask, 1], c=colors_zone[zone], 
                           s=5, alpha=0.6, label=zone)
        ax1.set_title(f'Sample: {sample} - Original Zones')
        ax1.set_xlabel('X coordinate')
        ax1.set_ylabel('Y coordinate')
        ax1.legend()
        
        # 右图：条带分布
        colors_band = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                      'inner': 'lightblue', 'outer': 'lightcoral'}
        
        # 获取当前样本的条带标签
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
        ax2.legend()
        
        plt.suptitle(f'Sample: {sample} - Zone Analysis', fontsize=16)
        plt.tight_layout()
        
        # 保存当前样本的条带区域图
        fig_path_png = os.path.join(save_dir, f'sample_{sample}_bands.png')
        fig_path_pdf = os.path.join(save_dir, f'sample_{sample}_bands.pdf')
        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
        fig.savefig(fig_path_pdf, bbox_inches='tight')
        print(f"  保存条带区域图: {fig_path_png}")
        print(f"  保存条带区域图: {fig_path_pdf}")
        
        plt.show()
        plt.close(fig)
    
    # 输出统计信息
    print(f"\n条带创建完成统计:")
    print(f"  成功处理的样本: {len(processed_samples)} 个")
    if processed_samples:
        print(f"  成功处理的样本列表: {', '.join(processed_samples)}")
    print(f"  跳过的样本: {len(skipped_samples)} 个")
    if skipped_samples:
        print(f"  跳过的样本列表: {', '.join(skipped_samples)}")
    
    return adata

def create_comprehensive_visualization(adata, sample_col='sample', save_dir='make_band_new'):
    """
    创建综合可视化图形，包括所有样本的条带统计
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 创建综合图形
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. 所有样本的条带分布统计
    ax1 = axes[0, 0]
    if 'site_band' in adata.obs.columns:
        # 统计所有细胞的条带分布
        band_counts = adata.obs['site_band'].value_counts()
        colors = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                  'inner': 'lightblue', 'outer': 'lightcoral'}
        bar_colors = [colors.get(band, 'gray') for band in band_counts.index]
        bars = ax1.bar(band_counts.index.astype(str), band_counts.values, color=bar_colors)
        ax1.set_title('Overall Band Distribution')
        ax1.set_xlabel('Band')
        ax1.set_ylabel('Cell Count')
        
        # 在柱状图上添加数字
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 5,
                    f'{int(height)}', ha='center', va='bottom')
    
    # 2. 每个样本的条带分布
    ax2 = axes[0, 1]
    if 'site_band' in adata.obs.columns and sample_col in adata.obs.columns:
        samples = adata.obs[sample_col].unique()
        band_colors = {'on': 'green', 'out': 'orange', 'in': 'purple', 
                       'inner': 'lightblue', 'outer': 'lightcoral'}
        
        # 获取所有样本的所有条带类型
        all_bands = adata.obs['site_band'].dropna().unique()
        
        # 为每个样本创建堆叠柱状图
        sample_data = []
        sample_names = []
        
        for sample in samples[:10]:  # 最多显示10个样本
            sample_mask = adata.obs[sample_col] == sample
            if sample_mask.any():
                sample_bands = adata.obs.loc[sample_mask, 'site_band'].value_counts()
                sample_data.append(sample_bands)
                sample_names.append(str(sample))
        
        if sample_data:
            # 创建DataFrame用于堆叠柱状图
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
            if len(all_bands) <= 10:  # 只有条带类型不多时才显示图例
                ax2.legend()
    
    # 3. 距离分布直方图
    ax3 = axes[1, 0]
    if 'boundary_distance' in adata.obs.columns:
        distances = adata.obs['boundary_distance'].dropna()
        if len(distances) > 0:
            ax3.hist(distances, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
            ax3.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Boundary (0)')
            ax3.axvline(x=-100, color='orange', linestyle=':', linewidth=1, label='Outer band edge')
            ax3.axvline(x=100, color='orange', linestyle=':', linewidth=1, label='Inner band edge')
            ax3.set_title('Distribution of Distance to Boundary')
            ax3.set_xlabel('Distance to Boundary (um)')
            ax3.set_ylabel('Frequency')
            ax3.legend()
    
    # 4. 边界细胞比例
    ax4 = axes[1, 1]
    if 'is_boundary' in adata.obs.columns and sample_col in adata.obs.columns:
        samples = adata.obs[sample_col].unique()
        boundary_ratios = []
        sample_names = []
        
        for sample in samples[:10]:  # 最多显示10个样本
            sample_mask = adata.obs[sample_col] == sample
            if sample_mask.any():
                sample_cells = adata.obs.loc[sample_mask, 'is_boundary']
                boundary_count = sample_cells.sum()
                total_count = len(sample_cells)
                if total_count > 0:
                    boundary_ratio = boundary_count / total_count * 100
                    boundary_ratios.append(boundary_ratio)
                    sample_names.append(str(sample))
        
        if boundary_ratios:
            bars = ax4.bar(range(len(boundary_ratios)), boundary_ratios, color='salmon')
            ax4.set_xticks(range(len(sample_names)))
            ax4.set_xticklabels(sample_names, rotation=45, ha='right')
            ax4.set_title('Percentage of Boundary Cells by Sample')
            ax4.set_xlabel('Sample')
            ax4.set_ylabel('Boundary Cells (%)')
            ax4.set_ylim(0, max(boundary_ratios) * 1.2 if boundary_ratios else 100)
            
            # 在柱状图上添加百分比
            for bar, ratio in zip(bars, boundary_ratios):
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                        f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.suptitle('Comprehensive Analysis Results', fontsize=16)
    plt.tight_layout()
    
    # 保存综合图形
    fig3_path_png = os.path.join(save_dir, 'comprehensive_analysis.png')
    fig3_path_pdf = os.path.join(save_dir, 'comprehensive_analysis.pdf')
    fig.savefig(fig3_path_png, dpi=300, bbox_inches='tight')
    fig.savefig(fig3_path_pdf, bbox_inches='tight')
    print(f"保存综合分析图: {fig3_path_png}")
    print(f"保存综合分析图: {fig3_path_pdf}")
    
    plt.show()
    plt.close(fig)

def save_adata_safely(adata, filename, save_dir='make_band_new'):
    """
    安全保存adata，避免不规则数组导致的错误
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 创建完整的文件路径
    filepath = os.path.join(save_dir, filename)
    
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
    
    # 特别处理'sample_boundaries'键，它肯定是不规则的
    if 'sample_boundaries' in adata_save.uns:
        keys_to_remove.append('sample_boundaries')
    
    # 移除不规则的键
    for key in keys_to_remove:
        print(f"警告：移除不规则的uns键: {key}")
        if key in adata_save.uns:
            del adata_save.uns[key]
    
    # 保存清理后的数据
    adata_save.write(filepath)
    print(f"数据已安全保存为 '{filepath}'")
    
    # 另外保存一个CSV文件，包含关键信息
    csv_path = os.path.join(save_dir, 'sample_band_summary.csv')
    
    # 创建摘要数据框
    summary_cols = ['mucus_zone', 'site_band', 'boundary_distance', 'is_boundary']
    
    # 添加样本信息
    sample_col = None
    for col in ['sample', 'Sample', 'tissue', 'Tissue']:
        if col in adata.obs.columns:
            sample_col = col
            summary_cols.append(col)
            break
    
    # 只选择存在的列
    existing_cols = [col for col in summary_cols if col in adata.obs.columns]
    summary_df = adata.obs[existing_cols].copy()
    
    summary_df.to_csv(csv_path, index=False)
    print(f"数据摘要已保存为 '{csv_path}'")

# 一站式处理函数
def process_mucus_zones_by_sample(adata, sample_col='sample', save_dir='make_band_new'):
    """
    一站式处理函数，按样本计算
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    print("步骤1: 按样本识别inner-outer分界线...")
    adata = calculate_boundary_by_sample(
        adata, 
        obs_col='mucus_zone',
        sample_col=sample_col,
        cat_a='inner', 
        cat_b='outer',
        n_neighbors=10,
        cat_focus='inner',
        smooth_sigma=2,
        max_dist_threshold=300,
        save_dir=save_dir
    )
    
    print("\n步骤2: 按样本创建条带标签...")
    adata = create_mucus_bands_by_sample(
        adata,
        sample_col=sample_col,
        on_band_width=100,
        out_band_width=200,
        in_band_width=200,
        save_dir=save_dir
    )
    
    print("\n步骤3: 创建综合分析图...")
    create_comprehensive_visualization(adata, sample_col=sample_col, save_dir=save_dir)
    
    print("\n步骤4: 保存数据...")
    save_adata_safely(adata, 'adata_with_band_by_sample.h5ad', save_dir=save_dir)
    
    # 显示最终的数据结构
    print("\n最终的adata.obs列:")
    summary_cols = ['mucus_zone', 'site_band', 'boundary_distance', 'is_boundary']
    if sample_col in adata.obs.columns:
        summary_cols.append(sample_col)
    
    # 只显示存在的列
    existing_cols = [col for col in summary_cols if col in adata.obs.columns]
    print(adata.obs[existing_cols].head())
    
    # 保存统计摘要
    save_statistics_summary(adata, sample_col, save_dir)
    
    return adata

def save_statistics_summary(adata, sample_col, save_dir):
    """
    保存统计摘要
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 保存统计摘要
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
        
        # 按样本统计
        if sample_col in adata.obs.columns:
            f.write("\n\n2. SAMPLE-WISE STATISTICS\n")
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
                        percentage = count / len(sample_cells) * 100
                        f.write(f"    {zone}: {count} ({percentage:.1f}%)\n")
                
                if 'site_band' in sample_cells.obs.columns:
                    band_counts = sample_cells.obs['site_band'].value_counts()
                    for band, count in band_counts.items():
                        if pd.notna(band):
                            percentage = count / len(sample_cells) * 100
                            f.write(f"    {band}: {count} ({percentage:.1f}%)\n")
        
        # 距离统计
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
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Analysis completed successfully!\n")
        f.write("=" * 60 + "\n")
    
    print(f"分析摘要已保存为: {summary_path}")

# 主程序
def main():
    # 加载数据
    adata = sc.read_h5ad('merged_TMA_data.h5ad')
    
    # 检查数据是否包含必要的列
    if 'mucus_zone' not in adata.obs.columns:
        raise ValueError("adata.obs中未找到'mucus_zone'列")
    
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm中未找到'spatial'坐标信息")
    
    # 检查是否有样本列
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
    
    # 按样本处理
    adata_result = process_mucus_zones_by_sample(adata, sample_col=sample_col)
    
    return adata_result

# 使用示例
if __name__ == "__main__":
    adata_result = main()
