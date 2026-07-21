import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path

def load_and_process_data(h5ad_path, boundaries_path):
    """
    加载并处理h5ad文件和边界数据
    """
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 检查细胞类型列
    if 'annotation' not in adata.obs.columns:
        raise ValueError("未找到annotation列，可用的列有: " + ", ".join(adata.obs.columns))
    
    # 安全处理重复的cell_id
    if 'cell_id' in adata.obs.columns:
        duplicate_count = adata.obs.duplicated(subset=['cell_id'], keep=False).sum()
        
        if duplicate_count > 0:
            print(f"发现 {duplicate_count} 个重复的cell_id，正在处理...")
            original_size = adata.n_obs
            keep_mask = ~adata.obs.duplicated(subset='cell_id', keep='first')
            adata = adata[keep_mask, :].copy()
            print(f"去重完成: {original_size} → {adata.n_obs} 个细胞")
        
        adata.obs.set_index('cell_id', inplace=True)
        print("已使用cell_id列作为索引")
    else:
        print("警告: 未找到cell_id列，使用默认索引")
    
    print(f"细胞数量: {adata.n_obs}, 基因数量: {adata.n_vars}")
    print("细胞类型分布:")
    print(adata.obs['annotation'].value_counts())
    
    # 读取边界文件
    print("正在读取细胞边界文件...")
    boundaries = pd.read_csv(boundaries_path)
    
    if 'cell_id' not in boundaries.columns:
        raise ValueError("边界文件必须包含cell_id列")
    
    # 确定坐标列
    x_col, y_col = None, None
    coordinate_columns = boundaries.columns.tolist()
    coordinate_columns.remove('cell_id')
    
    numeric_columns = boundaries[coordinate_columns].select_dtypes(include=[np.number]).columns
    if len(numeric_columns) >= 2:
        x_col, y_col = numeric_columns[0], numeric_columns[1]
        print(f"自动检测到坐标列: {x_col} 和 {y_col}")
    else:
        possible_x = ['x', 'X', 'vertex_x', 'coord_x', 'Coordinate_X']
        possible_y = ['y', 'Y', 'vertex_y', 'coord_y', 'Coordinate_Y']
        
        for col in possible_x:
            if col in boundaries.columns:
                x_col = col
                break
        for col in possible_y:
            if col in boundaries.columns:
                y_col = col
                break
    
    if x_col is None or y_col is None:
        raise ValueError("未找到坐标列，请确保文件包含x、y坐标")
    
    print(f"使用坐标列: {x_col} 和 {y_col}")
    print(f"边界文件包含 {boundaries['cell_id'].nunique()} 个独特的细胞ID")
    
    return adata, boundaries, x_col, y_col

def create_cell_polygons(boundaries, adata, x_col, y_col, sample_name):
    """
    为指定样本的细胞创建多边形对象
    """
    print(f"正在为样本 {sample_name} 创建细胞多边形...")
    
    sample_cell_ids = set(adata.obs[adata.obs['sample'] == sample_name].index)
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cell_ids)]
    
    if len(sample_boundaries) == 0:
        print(f"警告: 样本 {sample_name} 没有找到匹配的边界数据")
        return [], [], []
    
    cell_groups = sample_boundaries.groupby('cell_id')
    
    polygons = []
    poly_annotations = []
    poly_cellids = []
    
    matched_count = 0
    for cell_id, group in cell_groups:
        if cell_id in adata.obs.index:
            cell_type_value = adata.obs.loc[cell_id, 'annotation']
            
            if isinstance(cell_type_value, pd.Series):
                cell_type = cell_type_value.iloc[0]
            else:
                cell_type = cell_type_value
            
            points = group[[x_col, y_col]].values
            
            if len(points) >= 3:
                polygon = patches.Polygon(points, closed=True)
                polygons.append(polygon)
                poly_annotations.append(cell_type)
                poly_cellids.append(cell_id)
                matched_count += 1
    
    print(f"样本 {sample_name}: 成功匹配 {matched_count}/{len(sample_cell_ids)} 个细胞")
    return polygons, poly_annotations, poly_cellids

def create_annotation_colormap(annotations):
    """
    创建细胞类型颜色映射
    """
    # 使用更鲜明的颜色
    extended_colors = {
        "Stress-responsive_glandular_Epithelium": "#FF6B6B",  # 鲜艳红色
        "Proliferative_Epithelium": "#4ECDC4",                # 青绿色
        "Luminal_Epithelium": "#45B7D1",                     # 蓝色
        "FGFR2+_Basal_Epithelium": "#96CEB4",                # 浅绿色
        "Basal_Epithelium": "#FECA57",                       # 黄色
        "EMT_Epithelium": "#FF9FF3",                         # 粉色
        "Renal_Epithelium": "#54A0FF"                        # 天蓝色
    }
    
    all_annotations = set(annotations)
    target_annotations = set(["Stress-responsive_glandular_Epithelium","Proliferative_Epithelium",
                             "Luminal_Epithelium","FGFR2+_Basal_Epithelium","Basal_Epithelium",
                             "EMT_Epithelium","Renal_Epithelium"]) 
    
    color_dict = {}
    for annotation in all_annotations:
        if annotation in target_annotations and annotation in extended_colors:
            color_dict[annotation] = extended_colors[annotation]
        elif annotation in extended_colors:
            color_dict[annotation] = extended_colors[annotation]
        else:
            color_dict[annotation] = "#D3D3D3"
    
    print(f"颜色映射: 已为 {len([ct for ct in all_annotations if ct in extended_colors])} 种细胞类型分配颜色")
    
    return color_dict, target_annotations

def calculate_optimal_plot_area(polygons):
    """
    计算最优的绘图区域，确保细胞居中显示
    """
    if len(polygons) == 0:
        return 0, 0, 1, 1
    
    # 获取所有多边形的顶点
    all_vertices = np.vstack([poly.get_xy() for poly in polygons])
    x_coords = all_vertices[:, 0]
    y_coords = all_vertices[:, 1]
    
    # 计算细胞区域的边界
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    
    # 计算中心点
    center_x = (x_min + x_max) / 2
    center_y = (y_min + y_max) / 2
    
    # 计算区域大小
    width = x_max - x_min
    height = y_max - y_min
    
    # 增加10%的边距
    margin_factor = 0.1
    width_with_margin = width * (1 + margin_factor)
    height_with_margin = height * (1 + margin_factor)
    
    # 确保绘图区域是正方形
    plot_size = max(width_with_margin, height_with_margin)
    
    # 计算新的边界
    new_x_min = center_x - plot_size / 2
    new_x_max = center_x + plot_size / 2
    new_y_min = center_y - plot_size / 2
    new_y_max = center_y + plot_size / 2
    
    return new_x_min, new_x_max, new_y_min, new_y_max

def plot_optimized_spatial_maps(adata, boundaries, x_col, y_col, output_dir):
    """
    绘制优化的空间原位图
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    if 'sample' not in adata.obs.columns:
        print("警告: 未找到sample列，将所有细胞视为一个样本")
        samples = ['all_samples']
        adata.obs['sample'] = 'all_samples'
    else:
        samples = adata.obs['sample'].unique()
    
    print(f"开始处理 {len(samples)} 个样本")
    
    all_annotations = adata.obs['annotation'].unique()
    color_dict, target_annotations = create_annotation_colormap(all_annotations)
    
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        polygons, poly_annotations, poly_cellids = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample
        )
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            continue
        
        # 计算最优绘图区域
        x_min, x_max, y_min, y_max = calculate_optimal_plot_area(polygons)
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        
        # 处理细胞类型
        processed_annotations = []
        for ct in poly_annotations:
            if ct in target_annotations:
                processed_annotations.append(ct)
            else:
                processed_annotations.append("Others")
        
        # 创建颜色列表
        colors = []
        for ct in poly_annotations:
            if ct in target_annotations:
                colors.append(color_dict.get(ct, "#D3D3D3"))
            else:
                colors.append("#D3D3D3")
        
        # 绘制多边形
        p = PatchCollection(polygons, alpha=0.9, edgecolor='black', linewidth=0.1)
        p.set_facecolor(colors)
        ax.add_collection(p)
        
        # 设置坐标轴范围和比例
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect('equal')
        
        # 隐藏坐标轴
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # 添加标题
        ax.set_title(f'Sample: {sample}', fontsize=16, pad=20)
        
        # 创建图例
        legend_celltypes = list(target_annotations) + ["Others"]
        legend_elements = []
        for ct in legend_celltypes:
            if ct in set(processed_annotations):
                color = color_dict.get(ct, "#D3D3D3") if ct != "Others" else "#D3D3D3"
                legend_elements.append(patches.Patch(color=color, label=ct, alpha=0.9))
        
        if len(legend_elements) > 0:
            ax.legend(handles=legend_elements, loc='center left', 
                     bbox_to_anchor=(1.05, 0.5), fontsize=10, framealpha=0.9,
                     ncol=1, title="Cell Types")
        
        # 添加统计信息
        target_count = len([ct for ct in processed_annotations if ct in target_annotations])
        others_count = len([ct for ct in processed_annotations if ct == "Others"])
        
        stats_text = f'Total: {len(polygons):,} cells\nTarget: {target_count:,}\nOthers: {others_count:,}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        plt.tight_layout()
        
        # 保存图片
        sample_filename = f"spatial_optimized_{str(sample).replace(' ', '_').replace('/', '_')}.png"
        output_path = output_dir / sample_filename
        plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
        print(f"优化图片已保存: {output_path}")
        plt.close()
        
        # 记录统计信息
        sample_stats[sample] = {
            'total_cells': len(polygons),
            'target_cells': target_count,
            'others_cells': others_count,
            'output_path': output_path
        }
    
    return sample_stats

def main():
    """
    主函数
    """
    # 文件路径配置
    h5ad_path = "adata_celltype_RCTD_removed.h5ad"
    boundaries_path = "cell_boundaries.csv"
    output_dir = "spatial_optimized_maps"
    
    try:
        # 加载和处理数据
        adata, boundaries, x_col, y_col = load_and_process_data(h5ad_path, boundaries_path)
        
        # 检查细胞ID匹配
        h5ad_cells = set(adata.obs.index)
        boundary_cells = set(boundaries['cell_id'].unique())
        matched_cells = h5ad_cells & boundary_cells
        
        print(f"\n=== 细胞ID匹配检查 ===")
        print(f"h5ad文件中的细胞数量: {len(h5ad_cells)}")
        print(f"边界文件中的细胞数量: {len(boundary_cells)}")
        print(f"匹配的细胞数量: {len(matched_cells)}")
        print(f"匹配率: {len(matched_cells)/len(h5ad_cells)*100:.2f}%")
        
        if len(matched_cells) == 0:
            print("错误: 没有找到匹配的细胞ID！")
            return
        
        # 绘制优化的空间图
        sample_stats = plot_optimized_spatial_maps(adata, boundaries, x_col, y_col, output_dir)
        
        # 打印统计信息
        print("\n" + "="*50)
        print("处理完成！样本统计信息:")
        for sample, stats in sample_stats.items():
            print(f"样本 {sample}: {stats['total_cells']:,} 个细胞, "
                  f"目标细胞: {stats['target_cells']:,}, "
                  f"其他细胞: {stats['others_cells']:,}")
        
        print(f"\n所有图像已保存至: {output_dir}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
