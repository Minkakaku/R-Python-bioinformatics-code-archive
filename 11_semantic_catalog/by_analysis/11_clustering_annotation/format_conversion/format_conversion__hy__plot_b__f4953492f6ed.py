import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path

def load_and_process_data(h5ad_path, boundaries_path):
    """
    加载并处理h5ad文件和边界数据，自动处理重复的cell_id
    """
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 检查细胞类型列
    if 'celltype' not in adata.obs.columns:
        raise ValueError("未找到celltype列，可用的列有: " + ", ".join(adata.obs.columns))
    
    # 关键修改：处理重复的cell_id
    if 'cell_id' in adata.obs.columns:
        # 检查重复的cell_id
        duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep=False)
        duplicate_count = duplicate_mask.sum()
        
        if duplicate_count > 0:
            print(f"发现 {duplicate_count} 个重复的cell_id，正在处理...")
            
            # 显示重复的cell_id示例
            duplicate_cell_ids = adata.obs[duplicate_mask]['cell_id'].unique()
            print(f"重复cell_id示例: {list(duplicate_cell_ids[:5])}")
            
            # 去除重复，保留每个cell_id的第一个出现
            original_size = adata.n_obs
            adata.obs = adata.obs.drop_duplicates(subset='cell_id', keep='first')
            adata = adata[adata.obs.index, :]  # 同步更新整个AnnData对象
            
            print(f"去重完成: {original_size} → {adata.n_obs} 个细胞")
        
        # 设置cell_id为索引
        adata.obs.set_index('cell_id', inplace=True)
        print("已使用cell_id列作为索引")
    else:
        print("警告: 未找到cell_id列，使用默认索引")
    
    print(f"细胞数量: {adata.n_obs}, 基因数量: {adata.n_vars}")
    print("细胞类型分布:")
    print(adata.obs['celltype'].value_counts())
    
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
    poly_celltypes = []
    poly_cellids = []
    
    matched_count = 0
    for cell_id, group in cell_groups:
        # 检查该细胞是否有类型信息
        if cell_id in adata.obs.index:
            # 获取细胞类型 - 确保获取标量值
            cell_type_value = adata.obs.loc[cell_id, 'celltype']
            
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
                poly_celltypes.append(cell_type)
                poly_cellids.append(cell_id)
                matched_count += 1
    
    print(f"样本 {sample_name}: 成功匹配 {matched_count}/{len(sample_cell_ids)} 个细胞")
    return polygons, poly_celltypes, poly_cellids

def create_celltype_colormap(celltypes):
    """
    创建优化的细胞类型颜色映射
    关键修改：所有未指定的细胞类型均使用灰白色 (#D3D3D3)
    """
    # 定义主要细胞类型及其颜色（使用更鲜明的色彩）
    primary_colors = {
        # 肿瘤相关细胞类型 - 使用暖色调
        "Metabolical Tumor": "#E31A1C",  # 亮红色
        "Stress-Induced Epi": "#FB9A99",  # 浅红色
        "Metaplasia Tumor": "#33A02C",    # 绿色
        "Stem-like Tumor": "#B2DF8A",     # 浅绿色
        "Hypoxic Tumor": "#1F78B4",       # 蓝色
        "Proliferative Tumor": "#A6CEE3", # 浅蓝色
        
        # 正常/干细胞类型 - 使用中性色
        "Stem-like Normal": "#FF7F00",     # 橙色
        "Normal Epi": "#FDBF6F",          # 浅橙色
        
        # 特殊细胞类型
        "Paneth-like": "#6A3D9A",         # 紫色
        "Absorptive enterocyte": "#CAB2D6", # 浅紫色
        
        # 免疫细胞
        "B cells": "#FFFF99",         # 黄色
    }
    
    # 获取所有细胞类型
    all_celltypes = set(celltypes)
    
    # 关键修改：所有未指定的细胞类型统一使用灰白色
    grey_color = "#D3D3D3"  # 灰白色
    
    # 为每个细胞类型分配颜色
    color_dict = {}
    for celltype in all_celltypes:
        if celltype in primary_colors:
            color_dict[celltype] = primary_colors[celltype]
        else:
            color_dict[celltype] = grey_color  # 未指定的细胞类型使用灰白色
    
    print(f"颜色映射: 已为 {len(primary_colors)} 种指定细胞类型分配颜色，其余 {len(all_celltypes) - len(primary_colors)} 种使用灰白色")
    
    return color_dict

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
    all_celltypes = adata.obs['celltype'].unique()
    color_dict = create_celltype_colormap(all_celltypes)
    
    print(f"总细胞类型数量: {len(all_celltypes)}")
    print("颜色分配情况:")
    for celltype, color in color_dict.items():
        if celltype in all_celltypes:
            print(f"  {celltype}: {color}")
    
    # 样本级统计
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        # 创建细胞多边形
        polygons, poly_celltypes, poly_cellids = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample
        )
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            sample_stats[sample] = {'total_cells': 0, 'cell_types': 0}
            continue
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        
        # 创建颜色列表 - 使用get方法确保未指定的细胞类型使用灰白色
        colors = [color_dict.get(ct, "#D3D3D3") for ct in poly_celltypes]  # 默认灰白色
        
        # 使用PatchCollection批量添加多边形
        p = PatchCollection(polygons, alpha=0.8, edgecolor='white', linewidth=0.3)
        p.set_facecolor(colors)
        ax.add_collection(p)
        
        # 设置图形属性
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        ax.set_title(f'Spatial Map - Sample: {sample}\nCell Types: {len(set(poly_celltypes))}', 
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
        
        # 添加图例（显示所有细胞类型）
        sample_celltypes = list(set(poly_celltypes))
        sample_celltypes.sort()  # 按字母顺序排序
        
        # 计算图例列数（根据细胞类型数量动态调整）
        n_col = 1
        if len(sample_celltypes) > 15:
            n_col = 2
        if len(sample_celltypes) > 30:
            n_col = 3
        
        # 创建图例元素，区分指定细胞类型和未指定细胞类型
        legend_elements = []
        specified_celltypes = set([
            "Metabolical Tumor", "Stress-Induced Epi", "Metaplasia Tumor", 
            "Stem-like Tumor", "Hypoxic Tumor", "Proliferative Tumor",
            "Stem-like Normal", "Normal Epi", "Paneth-like", 
            "Absorptive enterocyte", "B cells"
        ])
        
        for ct in sample_celltypes:
            color = color_dict.get(ct, "#D3D3D3")
            # 为未指定细胞类型添加特殊标签
            label = ct
            if ct not in specified_celltypes:
                label = f"{ct}"
            legend_elements.append(patches.Patch(color=color, label=label))
        
        if len(legend_elements) > 0:
            ax.legend(handles=legend_elements, loc='center left', 
                     bbox_to_anchor=(1.05, 0.5), fontsize=8, framealpha=0.9,
                     ncol=n_col)
        
        # 添加统计信息
        cell_type_count = len(set(poly_celltypes))
        specified_count = len([ct for ct in set(poly_celltypes) if ct in specified_celltypes])
        unspecified_count = cell_type_count - specified_count
        
        stats_text = f'Total: {len(polygons)}\nCelltypes: {cell_type_count}\nTarget celltype: {specified_count}\nOthers: {unspecified_count}'
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
            'cell_types': cell_type_count,
            'specified_types': specified_count,
            'unspecified_types': unspecified_count,
            'output_path': output_path
        }
    
    # 打印总体统计信息
    print("\n=== 处理完成 ===")
    print("样本统计信息:")
    for sample, stats in sample_stats.items():
        print(f"样本 {sample}: {stats['total_cells']} 个细胞, {stats['cell_types']} 种细胞类型")
        print(f"  指定类型: {stats['specified_types']}, 未指定类型: {stats['unspecified_types']}")
    
    return sample_stats

def main():
    # 文件路径配置
    h5ad_path = "merged_celltypes_unique.h5ad"  # 替换为您的h5ad文件路径
    boundaries_path = "cell_boundaries.csv"  # 替换为您的边界文件路径
    output_dir = "spatial_Tumor_B"
    
    try:
        # 1. 加载和处理数据（自动去重）
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
