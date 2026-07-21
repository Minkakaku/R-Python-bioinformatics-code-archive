import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path
import os

def load_and_process_data(h5ad_path):
    """
    加载并处理h5ad文件，安全处理重复的cell_id
    """
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 检查必要的列
    required_cols = ['mucus_zone', 'sample', 'CNs_10']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"未找到{col}列，可用的列有: " + ", ".join(adata.obs.columns))
    
    # 安全处理重复的cell_id
    if 'cell_id' in adata.obs.columns:
        # 检查重复的cell_id
        duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep=False)
        duplicate_count = duplicate_mask.sum()
        
        if duplicate_count > 0:
            print(f"发现 {duplicate_count} 个重复的cell_id，正在处理...")
            
            # 创建去重掩码，保留每个cell_id的第一个出现
            keep_mask = ~adata.obs.duplicated(subset='cell_id', keep='first')
            
            # 使用掩码过滤整个AnnData对象
            adata = adata[keep_mask, :].copy()
            print(f"去重完成: {adata.n_obs} 个细胞")
        
        # 设置cell_id为索引
        adata.obs.set_index('cell_id', inplace=True)
        print("已使用cell_id列作为索引")
    else:
        print("警告: 未找到cell_id列，使用默认索引")
    
    print(f"细胞数量: {adata.n_obs}, 基因数量: {adata.n_vars}")
    print("样本分布:")
    print(adata.obs['sample'].value_counts())
    
    return adata

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
    poly_colors = []
    poly_cellids = []
    
    matched_count = 0
    for cell_id, group in cell_groups:
        # 检查该细胞是否有类型信息
        if cell_id in adata.obs.index:
            # 获取细胞的mucus_zone和CNs_10信息
            mucus_zone = adata.obs.loc[cell_id, 'mucus_zone']
            cns_10 = adata.obs.loc[cell_id, 'CNs_10']
            
            # 根据规则确定颜色
            if cns_10 == 'CN2':
                color = 'red'
            elif cns_10 == 'CN7':
                color = 'blue'
            elif cns_10 == 'CN8':
                color = 'yellow'
            elif mucus_zone == 'inner':
                color = '#333333'  # 深灰色
            elif mucus_zone == 'outer':
                color = '#CCCCCC'  # 浅灰色
            else:
                color = '#888888'  # 默认灰色
            
            # 获取坐标点
            points = group[[x_col, y_col]].values
            
            if len(points) >= 3:  # 至少需要3个点来创建多边形
                # 创建多边形对象
                polygon = patches.Polygon(points, closed=True)
                polygons.append(polygon)
                poly_colors.append(color)
                poly_cellids.append(cell_id)
                matched_count += 1
    
    print(f"样本 {sample_name}: 成功匹配 {matched_count}/{len(sample_cell_ids)} 个细胞")
    return polygons, poly_colors, poly_cellids

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
    
    # 样本级统计
    sample_stats = {}
    
    # 定义图例元素
    legend_elements = [
        patches.Patch(color='red', label='CN2'),
        patches.Patch(color='blue', label='CN7'),
        patches.Patch(color='yellow', label='CN8'),
        patches.Patch(color='#333333', label='mucus_zone: inner'),
        patches.Patch(color='#CCCCCC', label='mucus_zone: outer'),
        patches.Patch(color='#888888', label='Other')
    ]
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        # 创建细胞多边形
        polygons, poly_colors, poly_cellids = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample
        )
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            sample_stats[sample] = {'total_cells': 0}
            continue
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        
        # 使用PatchCollection批量添加多边形
        p = PatchCollection(polygons, alpha=0.8, edgecolor='white', linewidth=0.3)
        p.set_facecolor(poly_colors)
        ax.add_collection(p)
        
        # 设置图形属性
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        ax.set_title(f'Sample: {sample}\nColored by mucus_zone and CN2/7/8', 
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
            
            # 设置相等的纵横比
            ax.set_aspect('equal')
        
        # 添加图例
        ax.legend(handles=legend_elements, loc='center left', 
                 bbox_to_anchor=(1.05, 0.5), fontsize=10, framealpha=0.9,
                 ncol=1)
        
        # 统计各种类型的细胞数量
        cn2_count = sum(1 for color in poly_colors if color == 'red')
        cn7_count = sum(1 for color in poly_colors if color == 'blue')
        cn8_count = sum(1 for color in poly_colors if color == 'yellow')
        inner_count = sum(1 for color in poly_colors if color == '#333333')
        outer_count = sum(1 for color in poly_colors if color == '#CCCCCC')
        other_count = len(poly_colors) - cn2_count - cn7_count - cn8_count - inner_count - outer_count
        
        # 添加统计信息
        stats_text = f'Total cells: {len(polygons)}\n'
        stats_text += f'CN2: {cn2_count}\n'
        stats_text += f'CN7: {cn7_count}\n'
        stats_text += f'CN8: {cn8_count}\n'
        stats_text += f'Inner: {inner_count}\n'
        stats_text += f'Outer: {outer_count}\n'
        stats_text += f'Other: {other_count}'
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存图片为PDF和PNG
        safe_sample_name = str(sample).replace(' ', '_').replace('/', '_').replace(':', '_')
        
        # 保存为PDF
        pdf_path = output_dir / f"sample_{safe_sample_name}.pdf"
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        print(f"PDF已保存: {pdf_path}")
        
        # 保存为PNG
        png_path = output_dir / f"sample_{safe_sample_name}.png"
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        print(f"PNG已保存: {png_path}")
        
        plt.close(fig)
        
        # 记录统计信息
        sample_stats[sample] = {
            'total_cells': len(polygons),
            'CN2': cn2_count,
            'CN7': cn7_count,
            'CN8': cn8_count,
            'inner': inner_count,
            'outer': outer_count,
            'other': other_count
        }
    
    return sample_stats

def main():
    # 文件路径配置
    h5ad_path = "adata_with_CNs10_TMA-3_with_patient_MZ.h5ad"  # 替换为您的h5ad文件路径
    boundaries_path = "cell_boundaries.csv"  # 替换为您的边界文件路径
    output_dir = "CN_2_7_8"
    
    try:
        # 1. 加载和处理h5ad数据
        adata = load_and_process_data(h5ad_path)
        
        # 2. 读取边界文件
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
        
        # 3. 检查细胞ID匹配
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
        
        # 4. 为每个样本绘制空间图
        sample_stats = plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir)
        
        # 5. 保存汇总统计
        summary_df = pd.DataFrame.from_dict(sample_stats, orient='index')
        summary_path = os.path.join(output_dir, "sample_summary.csv")
        summary_df.to_csv(summary_path)
        print(f"\n汇总统计已保存: {summary_path}")
        
        print(f"\n所有处理完成！图像保存在: {output_dir}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
