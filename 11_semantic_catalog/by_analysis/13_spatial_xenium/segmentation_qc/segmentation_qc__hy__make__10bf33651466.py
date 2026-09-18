#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
import matplotlib.colors as mcolors
import warnings
import os
warnings.filterwarnings('ignore')

def calculate_mucus_zone_boundary_per_sample(adata_sample, sample_name, 
                                            obs_col='mucus_zone', 
                                            cat_a='inner', cat_b='outer', 
                                            n_neighbors=10, cat_focus='inner',
                                            smooth_sigma=2, min_vertex_count=200, 
                                            max_dist_threshold=300, shift_distance=0, 
                                            shift_direction='inner'):
    """
    为单个样本计算inner/outer边界线，并可向指定方向推进
    
    参数:
    ----------
    adata_sample : AnnData
        单个样本的空间转录组数据
    sample_name : str
        样本名称
    obs_col : str
        包含细胞标签的列名
    cat_a, cat_b : str
        两类细胞标签，如inner和outer
    n_neighbors : int
        KNN邻居数
    cat_focus : str
        关注的类别，用于构建概率场
    smooth_sigma : float
        高斯平滑的sigma值
    min_vertex_count : int
        最小顶点数，用于过滤短边界
    max_dist_threshold : float
        最大距离阈值，控制插值范围
    shift_distance : float
        推进距离（微米）
    shift_direction : str
        推进方向：'inner'（向inner推进）或 'outer'（向outer推进）
    
    返回:
    -------
    adata_sample : AnnData
        添加了边界信息的单个样本AnnData对象
    fig : matplotlib Figure
        边界图
    boundary_coords : dict
        边界坐标信息
    """
    
    print(f"\n处理样本: {sample_name}")
    print(f"  细胞数: {len(adata_sample)}")
    
    # 1. 获取坐标和标签
    if 'spatial' not in adata_sample.obsm:
        print(f"  警告: 样本 {sample_name} 没有空间坐标，跳过")
        return adata_sample, None, None
    
    coords = adata_sample.obsm['spatial']
    labels = adata_sample.obs[obs_col].values
    
    # 检查是否有足够的两种细胞类型
    inner_count = (labels == cat_a).sum()
    outer_count = (labels == cat_b).sum()
    
    print(f"  {cat_a}细胞数: {inner_count}")
    print(f"  {cat_b}细胞数: {outer_count}")
    
    if inner_count == 0 or outer_count == 0:
        print(f"  警告: 样本 {sample_name} 缺少 {cat_a} 或 {cat_b} 细胞，跳过边界检测")
        adata_sample.obs['is_boundary'] = False
        adata_sample.obs['dist_to_original_boundary'] = np.nan
        adata_sample.obs['signed_dist_to_boundary'] = np.nan
        return adata_sample, None, {}
    
    # 2. 构建空间邻居图 (KNN)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    
    # 3. 识别边界细胞
    boundary_mask = np.zeros(len(adata_sample), dtype=bool)
    for i in range(len(adata_sample)):
        current_label = labels[i]
        
        # 只关心 inner 或 outer 细胞
        if current_label not in [cat_a, cat_b]:
            continue
        
        # 获取邻居的标签
        neighbor_indices = indices[i]
        neighbor_labels = labels[neighbor_indices]
        
        # 检查邻居中是否存在"对立"的类别
        if current_label == cat_a and cat_b in neighbor_labels:
            boundary_mask[i] = True
        elif current_label == cat_b and cat_a in neighbor_labels:
            boundary_mask[i] = True
    
    adata_sample.obs['is_boundary'] = boundary_mask
    print(f"  识别出 {sum(boundary_mask)} 个边界细胞")
    
    # 4. 训练KNN分类器
    X = adata_sample.obsm['spatial']
    y = (adata_sample.obs[obs_col] == cat_focus).astype(int).values
    
    clf = KNeighborsClassifier(n_neighbors=n_neighbors)
    clf.fit(X, y)
    
    # 5. 创建网格
    padding_percent = 0.1
    x_min, x_max = X[:, 0].min(), X[:, 0].max()
    y_min, y_max = X[:, 1].min(), X[:, 1].max()
    x_pad = (x_max - x_min) * padding_percent
    y_pad = (y_max - y_min) * padding_percent
    
    grid_resolution = 500
    xx, yy = np.meshgrid(
        np.linspace(x_min - x_pad, x_max + x_pad, grid_resolution),
        np.linspace(y_min - y_pad, y_max + y_pad, grid_resolution)
    )
    grid_points = np.c_[xx.ravel(), yy.ravel()]
    
    # 计算网格点距离最近真实细胞的距离
    nbrs_geometry = NearestNeighbors(n_neighbors=1).fit(X)
    dists, _ = nbrs_geometry.kneighbors(grid_points)
    dists = dists.reshape(xx.shape)
    
    # 6. 预测概率 + 高斯平滑
    Z = clf.predict_proba(grid_points)[:, 1]
    Z = Z.reshape(xx.shape)
    Z_smooth = gaussian_filter(Z, sigma=smooth_sigma)
    
    # 如果网格点距离真实细胞太远，将其概率设为 NaN
    Z_smooth[dists > max_dist_threshold] = np.nan
    
    # 7. 创建图形
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # 绘制背景散点
    dot_size = 2
    is_boundary = boundary_mask
    
    # 绘制outer细胞
    mask_outer = (y == 0) & (~is_boundary)
    if mask_outer.any():
        ax.scatter(X[mask_outer, 0], X[mask_outer, 1], c='#F0A830', 
                  s=dot_size, alpha=0.8, linewidth=0, label=f'{cat_b} cells')
    
    # 绘制inner细胞
    mask_inner = (y == 1) & (~is_boundary)
    if mask_inner.any():
        ax.scatter(X[mask_inner, 0], X[mask_inner, 1], c='#7F1CD3', 
                  s=dot_size, alpha=0.8, linewidth=0, label=f'{cat_a} cells')
    
    # 绘制边界细胞
    if is_boundary.any():
        ax.scatter(X[is_boundary, 0], X[is_boundary, 1], c='red', 
                  s=dot_size*2, alpha=0.8, label='Boundary cells')
    
    # 提取原始边界轮廓
    contours = ax.contour(xx, yy, Z_smooth, levels=[0.5], alpha=0, colors='black')
    
    # 存储原始边界坐标
    original_boundary_coords = []
    boundary_coords_dict = {}
    
    if len(contours.allsegs) > 0:
        segments = contours.allsegs[0]
        for coords_segment in segments:
            if len(coords_segment) > min_vertex_count:
                ax.plot(coords_segment[:, 0], coords_segment[:, 1], 
                       color='black', linewidth=2, linestyle='-', label='Original boundary')
                original_boundary_coords.append(coords_segment)
    
    boundary_coords_dict['original'] = original_boundary_coords
    
    # 8. 如果指定了推进距离，创建推进后的边界
    if shift_distance > 0:
        print(f"  正在向 {shift_direction} 方向推进 {shift_distance} 单位距离...")
        shifted_boundary_coords = shift_boundary(
            original_boundary_coords, X, y, 
            shift_distance, shift_direction, cat_focus
        )
        
        # 绘制推进后的边界
        for i, coords_shifted in enumerate(shifted_boundary_coords):
            if len(coords_shifted) > min_vertex_count:
                color = 'blue' if shift_direction == 'inner' else 'green'
                label = f'Shifted boundary ({shift_direction})'
                ax.plot(coords_shifted[:, 0], coords_shifted[:, 1], 
                       color=color, linewidth=2, linestyle='--', label=label)
        
        boundary_coords_dict['shifted'] = shifted_boundary_coords
        
        # 计算细胞到推进后边界的距离
        if shifted_boundary_coords:
            all_shifted_points = np.vstack(shifted_boundary_coords)
            tree_shifted = cKDTree(all_shifted_points)
            dist_to_shifted, _ = tree_shifted.query(X, k=1)
            adata_sample.obs[f'dist_to_{shift_direction}_shifted_boundary'] = dist_to_shifted
            print(f"  已计算细胞到推进后边界的距离")
    
    # 9. 计算细胞到原始边界的距离
    if original_boundary_coords:
        all_boundary_points = np.vstack(original_boundary_coords)
        tree = cKDTree(all_boundary_points)
        dist_values, _ = tree.query(X, k=1)
        adata_sample.obs['dist_to_original_boundary'] = dist_values
        
        # 为距离添加符号，表示在边界的哪一侧
        signed_dist = dist_values.copy()
        signed_dist[y == 0] *= -1  # outer细胞距离为负
        
        adata_sample.obs['signed_dist_to_boundary'] = signed_dist
        print(f"  已计算细胞到原始边界的距离")
        print(f"    最小距离: {dist_values.min():.2f}, 最大距离: {dist_values.max():.2f}")
    else:
        # 如果没有找到边界，填充NaN
        adata_sample.obs['dist_to_original_boundary'] = np.nan
        adata_sample.obs['signed_dist_to_boundary'] = np.nan
        print(f"  警告: 未找到有效边界，距离设为NaN")
    
    ax.set_xlabel('X coordinate', fontsize=12)
    ax.set_ylabel('Y coordinate', fontsize=12)
    ax.set_title(f'Sample {sample_name}: Mucus Zone Boundary', fontsize=14, fontweight='bold')
    
    # 简化图例
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = []
    unique_handles = []
    for handle, label in zip(handles, labels):
        if label not in unique_labels:
            unique_labels.append(label)
            unique_handles.append(handle)
    
    if unique_handles:
        ax.legend(unique_handles, unique_labels, loc='upper right', fontsize=8)
    
    # Y轴坐标原点调整到左上角
    ylim = ax.get_ylim()
    if ylim[0] < ylim[1]:
        ax.invert_yaxis()
    
    plt.tight_layout()
    
    return adata_sample, fig, boundary_coords_dict

def shift_boundary(original_boundary_coords, cell_coords, cell_labels, 
                   shift_distance, shift_direction, cat_focus='inner'):
    """
    将边界线向指定方向推进
    
    参数:
    ----------
    original_boundary_coords : list
        原始边界坐标列表
    cell_coords : ndarray
        细胞坐标
    cell_labels : ndarray
        细胞标签 (1=关注类别, 0=其他)
    shift_distance : float
        推进距离
    shift_direction : str
        推进方向: 'inner' 或 'outer'
    cat_focus : str
        关注的类别
    
    返回:
    -------
    shifted_coords : list
        推进后的边界坐标列表
    """
    
    if not original_boundary_coords:
        return []
    
    shifted_coords = []
    
    for segment_coords in original_boundary_coords:
        if len(segment_coords) < 2:
            continue
        
        # 计算每个边界点的法向量
        normals = compute_normals(segment_coords)
        
        # 确定推进方向
        if shift_direction == 'inner':
            # 向inner方向推进，法向量指向inner的反方向
            direction = -1
        else:  # 'outer'
            # 向outer方向推进，法向量指向inner方向的反方向
            direction = 1
        
        # 推进边界点
        shifted_segment = segment_coords.copy()
        for i in range(len(segment_coords)):
            if i < len(normals):
                # 沿着法向量方向推进
                shifted_segment[i] = segment_coords[i] + normals[i] * shift_distance * direction
        
        shifted_coords.append(shifted_segment)
    
    return shifted_coords

def compute_normals(coords):
    """
    计算边界线上各点的法向量
    
    参数:
    ----------
    coords : ndarray
        边界点坐标，形状 (n_points, 2)
    
    返回:
    -------
    normals : ndarray
        单位法向量，形状 (n_points, 2)
    """
    n_points = len(coords)
    normals = np.zeros((n_points, 2))
    
    for i in range(n_points):
        if i == 0:
            # 第一个点：使用第一个线段的方向
            tangent = coords[1] - coords[0]
        elif i == n_points - 1:
            # 最后一个点：使用最后一个线段的方向
            tangent = coords[-1] - coords[-2]
        else:
            # 中间点：使用前后两点的平均方向
            tangent = (coords[i+1] - coords[i-1]) / 2.0
        
        # 计算法向量（旋转90度）
        # 切线 (dx, dy) 的法向量是 (-dy, dx) 或 (dy, -dx)
        # 这里我们选择指向边界外侧的法向量
        normal = np.array([-tangent[1], tangent[0]])
        
        # 归一化
        norm = np.linalg.norm(normal)
        if norm > 0:
            normal = normal / norm
        
        normals[i] = normal
    
    return normals

def analyze_boundary_zones(adata_sample, zone_widths=[50, 100, 150, 200]):
    """
    分析边界区域
    
    参数:
    ----------
    adata_sample : AnnData
        包含边界距离信息的数据
    zone_widths : list
        各区域的宽度列表
    
    返回:
    -------
    adata_sample : AnnData
        添加了区域标记的数据
    """
    
    if 'signed_dist_to_boundary' not in adata_sample.obs.columns:
        print("  警告: 没有边界距离信息，跳过区域分析")
        return adata_sample
    
    signed_dist = adata_sample.obs['signed_dist_to_boundary'].values
    
    # 检查是否有有效的距离值
    valid_mask = ~np.isnan(signed_dist)
    valid_count = np.sum(valid_mask)
    
    if valid_count == 0:
        print("  警告: 所有距离值都是NaN，跳过区域分析")
        return adata_sample
    
    print(f"  边界区域分析: {valid_count} 个有效距离值")
    
    # 定义区域
    # 只对有效距离进行区域划分
    valid_dist = signed_dist[valid_mask]
    
    # inner侧区域
    for i, width in enumerate(reversed(zone_widths)):
        zone_name = f'inner_zone_{i+1}'
        if i == 0:
            mask = valid_dist > 0
        else:
            prev_width = zone_widths[-i] if i < len(zone_widths) else 0
            mask = (valid_dist > prev_width) & (valid_dist <= width)
        
        # 初始化全为False
        adata_sample.obs[zone_name] = False
        # 只为有效距离设置True
        adata_sample.obs.loc[valid_mask, zone_name] = mask
    
    # 边界区域（靠近边界）
    boundary_mask = (valid_dist >= -zone_widths[0]) & (valid_dist <= zone_widths[0])
    adata_sample.obs['boundary_zone'] = False
    adata_sample.obs.loc[valid_mask, 'boundary_zone'] = boundary_mask
    
    # outer侧区域
    for i, width in enumerate(zone_widths):
        zone_name = f'outer_zone_{i+1}'
        if i == 0:
            mask = valid_dist < 0
        else:
            prev_width = -zone_widths[i-1] if i > 0 else 0
            mask = (valid_dist < prev_width) & (valid_dist >= -width)
        
        adata_sample.obs[zone_name] = False
        adata_sample.obs.loc[valid_mask, zone_name] = mask
    
    # 统计各区域细胞数
    print("  各区域细胞统计:")
    zone_columns = [f'inner_zone_{i+1}' for i in range(len(zone_widths))] + \
                  ['boundary_zone'] + \
                  [f'outer_zone_{i+1}' for i in range(len(zone_widths))]
    
    for col in zone_columns:
        if col in adata_sample.obs.columns:
            count = adata_sample.obs[col].sum()
            percentage = count / len(adata_sample) * 100
            print(f"    {col}: {count} 细胞 ({percentage:.1f}%)")
    
    return adata_sample

def visualize_boundary_zones_per_sample(adata_sample, sample_name, figsize=(10, 8)):
    """
    为单个样本可视化边界区域
    
    参数:
    ----------
    adata_sample : AnnData
        单个样本的数据
    sample_name : str
        样本名称
    figsize : tuple
        图形大小
    
    返回:
    -------
    fig : matplotlib Figure or None
        图形对象，如果无法绘制则返回None
    """
    
    if 'signed_dist_to_boundary' not in adata_sample.obs.columns:
        print(f"  样本 {sample_name} 没有边界距离信息，跳过可视化")
        return None
    
    signed_dist = adata_sample.obs['signed_dist_to_boundary'].values
    
    # 检查是否有有效的距离值
    valid_mask = ~np.isnan(signed_dist)
    valid_count = np.sum(valid_mask)
    
    if valid_count == 0:
        print(f"  样本 {sample_name} 所有距离值都是NaN，跳过可视化")
        return None
    
    # 提取有效距离
    valid_signed_dist = signed_dist[valid_mask]
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # 1. 边界距离分布直方图
    ax1 = axes[0]
    
    # 确保有有效的距离值
    if len(valid_signed_dist) > 0:
        ax1.hist(valid_signed_dist, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax1.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Boundary')
        ax1.set_xlabel('Signed Distance to Boundary', fontsize=12)
        ax1.set_ylabel('Cell Count', fontsize=12)
        ax1.set_title(f'Sample {sample_name}: Distance Distribution\n({valid_count} valid values)', 
                     fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, 'No valid distance values', 
                horizontalalignment='center', verticalalignment='center',
                transform=ax1.transAxes, fontsize=12)
        ax1.set_title(f'Sample {sample_name}: No Distance Data', fontsize=14, fontweight='bold')
    
    # 2. 边界区域可视化
    ax2 = axes[1]
    if 'spatial' in adata_sample.obsm and valid_count > 0:
        coords = adata_sample.obsm['spatial']
        
        # 只显示有有效距离的细胞
        valid_coords = coords[valid_mask]
        valid_signed_dist = signed_dist[valid_mask]
        
        # 定义颜色映射
        cmap = plt.cm.coolwarm
        
        # 确保有最小值和最大值
        vmin = np.nanmin(valid_signed_dist)
        vmax = np.nanmax(valid_signed_dist)
        
        # 如果所有值相同，添加一个小的偏移
        if vmin == vmax:
            vmin -= 0.1
            vmax += 0.1
        
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        
        sc = ax2.scatter(valid_coords[:, 0], valid_coords[:, 1], c=valid_signed_dist, 
                         cmap=cmap, s=2, alpha=0.7, norm=norm)
        
        ax2.set_xlabel('X coordinate', fontsize=12)
        ax2.set_ylabel('Y coordinate', fontsize=12)
        ax2.set_title(f'Sample {sample_name}: Boundary Zones\n({valid_count} cells)', 
                     fontsize=14, fontweight='bold')
        
        # 添加颜色条
        cbar = plt.colorbar(sc, ax=ax2)
        cbar.set_label('Distance to Boundary (signed)', fontsize=10)
        
        # 反转y轴
        ylim = ax2.get_ylim()
        if ylim[0] < ylim[1]:
            ax2.invert_yaxis()
    else:
        ax2.text(0.5, 0.5, 'No spatial data or valid distances', 
                horizontalalignment='center', verticalalignment='center',
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title(f'Sample {sample_name}: No Spatial Data', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    return fig

def save_boundary_coords_to_file(boundary_coords_dict, sample_name, output_dir):
    """
    将边界坐标保存到文件
    
    参数:
    ----------
    boundary_coords_dict : dict
        边界坐标字典
    sample_name : str
        样本名称
    output_dir : str
        输出目录
    """
    
    if not boundary_coords_dict:
        return
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 为每个样本创建一个子目录
    sample_dir = os.path.join(output_dir, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    
    for boundary_type, coords_list in boundary_coords_dict.items():
        if not coords_list:
            continue
        
        # 将坐标保存为文本文件
        output_file = os.path.join(sample_dir, f'{boundary_type}_boundary_coords.txt')
        with open(output_file, 'w') as f:
            f.write(f"# {boundary_type.capitalize()} Boundary Coordinates for Sample {sample_name}\n")
            f.write("# Format: segment_index, point_index, x, y\n")
            
            for seg_idx, segment in enumerate(coords_list):
                for pt_idx, (x, y) in enumerate(segment):
                    f.write(f"{seg_idx},{pt_idx},{x:.2f},{y:.2f}\n")
        
        print(f"  已保存{boundary_type}边界坐标到: {output_file}")
        
        # 同时保存为numpy文件
        npz_file = os.path.join(sample_dir, f'{boundary_type}_boundary_coords.npz')
        
        # 将每个线段保存为独立的数组
        save_dict = {}
        for seg_idx, segment in enumerate(coords_list):
            save_dict[f'segment_{seg_idx}'] = segment
        
        np.savez(npz_file, **save_dict)
        print(f"  已保存{boundary_type}边界坐标到: {npz_file}")

def main():
    """
    主函数
    """
    # 配置参数
    h5ad_file = "merged_TMA_data.h5ad"  # 修改为您的文件路径
    
    # 边界检测参数
    params = {
        'obs_col': 'mucus_zone',
        'cat_a': 'inner',
        'cat_b': 'outer',
        'n_neighbors': 15,
        'cat_focus': 'inner',
        'smooth_sigma': 3,
        'min_vertex_count': 100,
        'max_dist_threshold': 200,
        'shift_distance': 50,  # 推进距离（单位与坐标相同）
        'shift_direction': 'outer'  # 'inner' 或 'outer'
    }
    
    # 输出目录
    output_dir = "boundary_analysis_results"
    
    print("="*60)
    print("Mucus Zone边界分析 (按样本处理)")
    print("="*60)
    
    try:
        # 1. 读取数据
        print("正在读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        
        # 检查必要的列
        if 'sample' not in adata.obs.columns:
            raise ValueError("在adata.obs中未找到'sample'列")
        
        if 'mucus_zone' not in adata.obs.columns:
            raise ValueError("在adata.obs中未找到'mucus_zone'列")
        
        if 'spatial' not in adata.obsm:
            print("警告: 在adata.obsm中未找到'spatial'坐标，跳过空间分析")
            return
        
        # 获取所有样本
        all_samples = adata.obs['sample'].unique()
        print(f"\n发现 {len(all_samples)} 个样本: {list(all_samples)}")
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 存储所有样本的结果
        all_results = []
        all_boundary_coords = {}
        
        # 统计信息
        stats_summary = []
        
        # 2. 为每个样本单独计算边界
        for sample_idx, sample_name in enumerate(all_samples):
            print(f"\n{'='*60}")
            print(f"处理样本 {sample_idx+1}/{len(all_samples)}: {sample_name}")
            print('='*60)
            
            # 提取当前样本的数据
            sample_mask = adata.obs['sample'] == sample_name
            adata_sample = adata[sample_mask].copy()
            
            if len(adata_sample) == 0:
                print(f"  警告: 样本 {sample_name} 没有数据，跳过")
                continue
            
            # 初始化样本统计
            sample_stats = {
                'sample': sample_name,
                'total_cells': len(adata_sample),
                'inner_cells': (adata_sample.obs['mucus_zone'] == 'inner').sum(),
                'outer_cells': (adata_sample.obs['mucus_zone'] == 'outer').sum(),
                'boundary_cells': 0,
                'has_valid_boundary': False,
                'valid_distance_cells': 0
            }
            
            # 计算边界
            adata_sample, fig_boundary, boundary_coords = calculate_mucus_zone_boundary_per_sample(
                adata_sample, sample_name, **params
            )
            
            # 更新统计
            sample_stats['boundary_cells'] = adata_sample.obs['is_boundary'].sum()
            if 'signed_dist_to_boundary' in adata_sample.obs.columns:
                valid_dist_count = adata_sample.obs['signed_dist_to_boundary'].notna().sum()
                sample_stats['valid_distance_cells'] = valid_dist_count
                sample_stats['has_valid_boundary'] = valid_dist_count > 0
            
            stats_summary.append(sample_stats)
            
            if fig_boundary is not None:
                # 保存边界图
                boundary_plot_path = os.path.join(output_dir, f"boundary_sample_{sample_name}.png")
                fig_boundary.savefig(boundary_plot_path, dpi=300, bbox_inches='tight')
                plt.savefig(os.path.join(output_dir, f"boundary_sample_{sample_name}.pdf"), bbox_inches='tight')
                plt.close(fig_boundary)
                print(f"  已保存边界图: boundary_sample_{sample_name}.png/pdf")
            
            # 保存边界坐标
            if boundary_coords:
                save_boundary_coords_to_file(boundary_coords, sample_name, 
                                           os.path.join(output_dir, 'boundary_coords'))
                all_boundary_coords[sample_name] = boundary_coords
            
            # 分析边界区域
            adata_sample = analyze_boundary_zones(adata_sample, zone_widths=[50, 100, 150])
            
            # 可视化边界区域
            fig_zones = visualize_boundary_zones_per_sample(adata_sample, sample_name)
            if fig_zones is not None:
                zones_plot_path = os.path.join(output_dir, f"zones_sample_{sample_name}.png")
                fig_zones.savefig(zones_plot_path, dpi=300, bbox_inches='tight')
                plt.savefig(os.path.join(output_dir, f"zones_sample_{sample_name}.pdf"), bbox_inches='tight')
                plt.close(fig_zones)
                print(f"  已保存区域图: zones_sample_{sample_name}.png/pdf")
            else:
                print(f"  跳过区域图: 没有有效的边界距离数据")
            
            # 存储结果
            all_results.append(adata_sample)
        
        # 3. 合并所有样本的结果
        if all_results:
            print(f"\n{'='*60}")
            print("合并所有样本的结果...")
            print('='*60)
            
            # 合并所有样本的AnnData对象
            adata_combined = sc.concat(all_results, join='outer', label='original_sample_index')
            
            # 保存合并后的结果
            output_file = os.path.join(output_dir, "adata_with_boundary.h5ad")
            adata_combined.write_h5ad(output_file)
            print(f"已保存合并结果到: {output_file}")
            
            # 4. 创建结果摘要
            create_summary_report(adata_combined, params, output_dir, all_samples, stats_summary)
            
        else:
            print("没有成功处理任何样本")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

def create_summary_report(adata, params, output_dir, all_samples, stats_summary=None):
    """
    创建结果摘要报告
    """
    report_path = os.path.join(output_dir, "boundary_analysis_summary.txt")
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("MUCUS ZONE边界分析结果摘要 (按样本处理)\n")
        f.write("="*80 + "\n\n")
        
        f.write("1. 分析参数\n")
        f.write("-"*40 + "\n")
        for key, value in params.items():
            f.write(f"  {key}: {value}\n")
        
        f.write("\n2. 数据概览\n")
        f.write("-"*40 + "\n")
        f.write(f"总细胞数: {adata.n_obs}\n")
        f.write(f"总样本数: {len(all_samples)}\n")
        f.write(f"mucus_zone分布:\n")
        for zone, count in adata.obs['mucus_zone'].value_counts().items():
            f.write(f"  {zone}: {count} 细胞 ({count/adata.n_obs*100:.1f}%)\n")
        
        if 'is_boundary' in adata.obs.columns:
            boundary_count = adata.obs['is_boundary'].sum()
            f.write(f"边界细胞总数: {boundary_count} ({boundary_count/adata.n_obs*100:.1f}%)\n")
        
        f.write("\n3. 边界距离统计\n")
        f.write("-"*40 + "\n")
        if 'dist_to_original_boundary' in adata.obs.columns:
            dist_data = adata.obs['dist_to_original_boundary'].dropna()
            if len(dist_data) > 0:
                f.write(f"平均距离: {dist_data.mean():.2f}\n")
                f.write(f"距离标准差: {dist_data.std():.2f}\n")
                f.write(f"最小距离: {dist_data.min():.2f}\n")
                f.write(f"最大距离: {dist_data.max():.2f}\n")
                f.write(f"距离中位数: {dist_data.median():.2f}\n")
            else:
                f.write("没有有效的距离数据\n")
        
        f.write("\n4. 样本统计详情\n")
        f.write("-"*40 + "\n")
        if stats_summary:
            for stats in stats_summary:
                f.write(f"\n样本: {stats['sample']}\n")
                f.write(f"  总细胞数: {stats['total_cells']}\n")
                f.write(f"  inner细胞: {stats['inner_cells']} ({stats['inner_cells']/stats['total_cells']*100:.1f}%)\n")
                f.write(f"  outer细胞: {stats['outer_cells']} ({stats['outer_cells']/stats['total_cells']*100:.1f}%)\n")
                f.write(f"  边界细胞: {stats['boundary_cells']} ({stats['boundary_cells']/stats['total_cells']*100:.1f}%)\n")
                f.write(f"  有效距离细胞: {stats['valid_distance_cells']} ({stats['valid_distance_cells']/stats['total_cells']*100:.1f}%)\n")
                f.write(f"  成功检测边界: {'是' if stats['has_valid_boundary'] else '否'}\n")
        
        f.write("\n5. 边界区域统计\n")
        f.write("-"*40 + "\n")
        zone_columns = [col for col in adata.obs.columns if 'zone' in col and col != 'mucus_zone']
        for col in zone_columns:
            count = adata.obs[col].sum()
            f.write(f"{col}: {count} 细胞 ({count/adata.n_obs*100:.1f}%)\n")
        
        f.write("\n6. 输出文件\n")
        f.write("-"*40 + "\n")
        f.write("边界图 (每个样本): boundary_sample_*.png/pdf\n")
        f.write("区域图 (每个样本): zones_sample_*.png/pdf\n")
        f.write("边界坐标文件: boundary_coords/目录下\n")
        f.write("adata_with_boundary.h5ad - 包含边界信息的数据文件\n")
        f.write("boundary_analysis_summary.txt - 本摘要文件\n")
        
        f.write("\n分析完成!\n")
    
    print(f"已保存摘要报告: {report_path}")
    
    # 同时将统计摘要保存为CSV
    if stats_summary:
        stats_df = pd.DataFrame(stats_summary)
        stats_csv_path = os.path.join(output_dir, "sample_statistics.csv")
        stats_df.to_csv(stats_csv_path, index=False)
        print(f"已保存样本统计到: {stats_csv_path}")

if __name__ == "__main__":
    main()
