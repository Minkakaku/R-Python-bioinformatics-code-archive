import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path

def load_and_process_data(h5ad_path, boundaries_path):
    """
    加载并处理h5ad文件和边界数据，安全处理重复的cell_id
    """
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 检查细胞类型列
    if 'annotation' not in adata.obs.columns:
        raise ValueError("未找到annotation列，可用的列有: " + ", ".join(adata.obs.columns))
    
    # 关键修复：安全处理重复的cell_id
    if 'cell_id' in adata.obs.columns:
        # 检查重复的cell_id
        duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep=False)
        duplicate_count = duplicate_mask.sum()
        
        if duplicate_count > 0:
            print(f"发现 {duplicate_count} 个重复的cell_id，正在处理...")
            
            # 显示重复的cell_id示例
            duplicate_cell_ids = adata.obs[duplicate_mask]['cell_id'].unique()
            print(f"重复cell_id示例: {list(duplicate_cell_ids[:5])}")
            
            # 修复：使用布尔索引而不是直接替换obs
            original_size = adata.n_obs
            
            # 创建去重掩码，保留每个cell_id的第一个出现
            keep_mask = ~adata.obs.duplicated(subset='cell_id', keep='first')
            
            # 使用掩码过滤整个AnnData对象
            adata = adata[keep_mask, :].copy()
            
            print(f"去重完成: {original_size} → {adata.n_obs} 个细胞")
        
        # 设置cell_id为索引
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
    
    # 检查必要的列
    if 'cell_id' not in boundaries.columns:
        raise ValueError("边界文件必须包含cell_id列")
    
    # 确定坐标列
    x_col, y_col = None, None
    coordinate_columns = boundaries.columns.tolist()
    coordinate_columns.remove('cell_id')  # 移除细胞ID列
    
    # 通常坐标列是数值型的前两列
    numeric_columns = boundaries[coordinate_columns].select_dtypes(include=[np.number]).columns
    if len(numeric_columns) >= 2:
        x_col, y_col = numeric_columns[0], numeric_columns[1]
        print(f"自动检测到坐标列: {x_col} 和 {y_col}")
    else:
        # 如果自动检测失败，尝试常见列名
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
    
    # 获取当前样本的细胞ID
    sample_cell_ids = set(adata.obs[adata.obs['sample'] == sample_name].index)
    
    # 过滤边界数据，只包含当前样本的细胞
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cell_ids)]
    
    if len(sample_boundaries) == 0:
        print(f"警告: 样本 {sample_name} 没有找到匹配的边界数据")
        return [], [], []
    
    # 按细胞ID分组边界点
    cell_groups = sample_boundaries.groupby('cell_id')
    
    polygons = []
    poly_annotations = []
    poly_cellids = []
    
    matched_count = 0
    for cell_id, group in cell_groups:
        # 检查该细胞是否有类型信息
        if cell_id in adata.obs.index:
            # 获取细胞类型 - 确保获取标量值
            cell_type_value = adata.obs.loc[cell_id, 'annotation']
            
            # 确保cell_type是字符串而不是Series
            if isinstance(cell_type_value, pd.Series):
                cell_type = cell_type_value.iloc[0]  # 取第一个值
            else:
                cell_type = cell_type_value  # 已经是标量
            
            # 获取坐标点
            points = group[[x_col, y_col]].values
            
            if len(points) >= 3:  # 至少需要3个点来创建多边形
                # 创建多边形对象
                polygon = patches.Polygon(points, closed=True)
                polygons.append(polygon)
                poly_annotations.append(cell_type)
                poly_cellids.append(cell_id)
                matched_count += 1
    
    print(f"样本 {sample_name}: 成功匹配 {matched_count}/{len(sample_cell_ids)} 个细胞")
    return polygons, poly_annotations, poly_cellids

def create_annotation_colormap(annotations):
    """
    创建优化的细胞类型颜色映射
    """
    # 扩展的颜色字典 - 包含更多颜色选择
    extended_colors = {
        "Stress-responsive_glandular_Epithelium": "#FFA500",
        "Proliferative_Epithelium": "#FF7F00",
        "Luminal_Epithelium": "#B2DF8A",
        "FGFR2+_Basal_Epithelium": "#6A3D9A",
        "Basal_Epithelium": "#FB9A99",
        "EMT_Epithelium": "#33A02C",
        "Renal_Epithelium": "#B15928"
    }
    
    # 获取所有细胞类型
    all_annotations = set(annotations)
    
    # 目标细胞类型列表
    target_annotations = set(["Stress-responsive_glandular_Epithelium","Proliferative_Epithelium","Luminal_Epithelium","FGFR2+_Basal_Epithelium","Basal_Epithelium","EMT_Epithelium","Renal_Epithelium"]) 
    
    # 为每个细胞类型分配颜色
    color_dict = {}
    for annotation in all_annotations:
        if annotation in target_annotations and annotation in extended_colors:
            color_dict[annotation] = extended_colors[annotation]
        elif annotation in extended_colors:
            # 非目标但已定义的细胞类型也使用指定颜色
            color_dict[annotation] = extended_colors[annotation]
        else:
            # 未指定的细胞类型使用灰白色
            color_dict[annotation] = "#D3D3D3"
    
    print(f"颜色映射: 已为 {len([ct for ct in all_annotations if ct in extended_colors])} 种细胞类型分配颜色，其余使用灰白色")
    
    return color_dict, target_annotations

def plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir):
    """
    为每个样本绘制空间原位图
    """
    # 创建输出目录
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # 获取所有样本
    if 'sample' not in adata.obs.columns:
        print("警告: 未找到sample列，将所有细胞视为一个样本")
        samples = ['all_samples']
        adata.obs['sample'] = 'all_samples'
    else:
        samples = adata.obs['sample'].unique()
    
    print(f"开始处理 {len(samples)} 个样本")
    
    # 获取所有细胞类型并创建颜色映射
    all_annotations = adata.obs['annotation'].unique()
    color_dict, target_annotations = create_annotation_colormap(all_annotations)
    
    print(f"总细胞类型数量: {len(all_annotations)}")
    print("目标细胞类型:", target_annotations)
    
    # 样本级统计
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        # 创建细胞多边形
        polygons, poly_annotations, poly_cellids = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample
        )
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            sample_stats[sample] = {'total_cells': 0, 'cell_types': 0}
            continue
        
        # 将非目标细胞类型重命名为"Others"
        processed_annotations = []
        for ct in poly_annotations:
            if ct in target_annotations:
                processed_annotations.append(ct)
            else:
                processed_annotations.append("Others")
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        
        # 创建颜色列表
        colors = []
        for ct in poly_annotations:
            if ct in target_annotations:
                colors.append(color_dict.get(ct, "#D3D3D3"))
            else:
                colors.append("#D3D3D3")  # 非目标细胞类型统一使用灰白色
        
        # 使用PatchCollection批量添加多边形
        p = PatchCollection(polygons, alpha=0.8, edgecolor='white', linewidth=0.3)
        p.set_facecolor(colors)
        ax.add_collection(p)
        
        # 设置图形属性
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        ax.set_title(f'Spatial Map - Sample: {sample}\nCell Types: {len(set(processed_annotations))}', 
                    fontsize=16, pad=20)
        
        # 自动调整坐标轴范围
        if len(polygons) > 0:
            all_vertices = np.vstack([poly.get_xy() for poly in polygons])
            x_min, x_max = all_vertices[:, 0].min(), all_vertices[:, 0].max()
            y_min, y_max = all_vertices[:, 1].min(), all_vertices[:, 1].max()
            
            # 添加5%的边距
            x_margin = (x_max - x_min) * 0.05
            y_margin = (y_max - y_min) * 0.05
            
            ax.set_xlim(x_min - x_margin, x_max + x_margin)
            ax.set_ylim(y_min - y_margin, y_max + y_margin)
            
            # 关键修复：设置相等的纵横比
            ax.set_aspect('equal')
        
        # 添加图例 - 只显示目标细胞类型和"Others"
        legend_annotations = list(target_annotations) + ["Others"]
        legend_colors = [color_dict.get(ct, "#D3D3D3") if ct != "Others" else "#D3D3D3" 
                        for ct in legend_annotations if ct in set(processed_annotations)]
        
        # 创建图例元素
        legend_elements = []
        for ct in legend_annotations:
            if ct in set(processed_annotations):
                color = color_dict.get(ct, "#D3D3D3") if ct != "Others" else "#D3D3D3"
                legend_elements.append(patches.Patch(color=color, label=ct))
        
        if len(legend_elements) > 0:
            ax.legend(handles=legend_elements, loc='center left', 
                     bbox_to_anchor=(1.05, 0.5), fontsize=10, framealpha=0.9,
                     ncol=1)
        
        # 添加统计信息
        target_count = len([ct for ct in processed_annotations if ct in target_annotations])
        others_count = len([ct for ct in processed_annotations if ct == "Others"])
        
        stats_text = f'Total cells: {len(polygons)}\nTarget types: {len(target_annotations & set(processed_annotations))}\nTarget cells: {target_count}\nOthers: {others_count}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存图片
        sample_filename = f"spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.png"
        sample_filename1 = f"spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.pdf"
        output_path = output_dir / sample_filename
        output_path1 = output_dir / sample_filename1
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path1, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {output_path}")
        plt.close(fig)
        
        # 记录统计信息
        sample_stats[sample] = {
            'total_cells': len(polygons),
            'target_cells': target_count,
            'others_cells': others_count,
            'output_path': output_path
        }
    
    # 打印总体统计信息
    print("\n=== 处理完成 ===")
    print("样本统计信息:")
    for sample, stats in sample_stats.items():
        print(f"样本 {sample}: {stats['total_cells']} 个细胞, 目标细胞: {stats['target_cells']}, 其他细胞: {stats['others_cells']}")
    
    return sample_stats

def main():
    # 文件路径配置
    h5ad_path = "TMA_4.h5ad"  # 替换为您的h5ad文件路径
    boundaries_path = "cell_boundaries.csv"  # 替换为您的边界文件路径
    output_dir = "spatial_Tumor_TMA_4"
    
    try:
        # 1. 加载和处理数据（安全去重）
        adata, boundaries, x_col, y_col = load_and_process_data(h5ad_path, boundaries_path)
        
        # 2. 检查细胞ID匹配
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
            print("h5ad中的前5个细胞ID:", list(h5ad_cells)[:5])
            print("边界文件中的前5个细胞ID:", list(boundary_cells)[:5])
            return
        
        # 3. 为每个样本绘制空间图
        sample_stats = plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir)
        
        print(f"\n所有处理完成！图像保存在: {output_dir}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
