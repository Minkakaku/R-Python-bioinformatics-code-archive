import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path
import argparse
import sys
import os

def load_and_process_data(h5ad_path, boundaries_path):
    """
    加载并处理h5ad文件和边界数据
    """
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 检查必要的列
    required_cols = ['annotation', 'CNs_10', 'sample']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"未找到{col}列，可用的列有: " + ", ".join(adata.obs.columns))
    
    # 处理重复的cell_id
    if 'cell_id' in adata.obs.columns:
        duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep=False)
        duplicate_count = duplicate_mask.sum()
        
        if duplicate_count > 0:
            print(f"发现 {duplicate_count} 个重复的cell_id，正在处理...")
            keep_mask = ~adata.obs.duplicated(subset='cell_id', keep='first')
            adata = adata[keep_mask, :].copy()
            print(f"去重完成: {adata.n_obs} 个细胞")
        
        adata.obs.set_index('cell_id', inplace=True)
        print("已使用cell_id列作为索引")
    else:
        print("警告: 未找到cell_id列，使用默认索引")
    
    print(f"细胞数量: {adata.n_obs}, 基因数量: {adata.n_vars}")
    print("CNs_10分布:")
    print(adata.obs['CNs_10'].value_counts(dropna=False))
    print("样本数量:", adata.obs['sample'].nunique())
    
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
    
    # 获取当前样本的细胞ID
    sample_cell_ids = set(adata.obs[adata.obs['sample'] == sample_name].index)
    
    # 过滤边界数据，只包含当前样本的细胞
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cell_ids)]
    
    if len(sample_boundaries) == 0:
        print(f"警告: 样本 {sample_name} 没有找到匹配的边界数据")
        return [], [], [], []
    
    # 按细胞ID分组边界点
    cell_groups = sample_boundaries.groupby('cell_id')
    
    polygons = []
    poly_annotations = []
    poly_cellids = []
    is_cn2_list = []
    
    matched_count = 0
    cn2_count = 0
    
    for cell_id, group in cell_groups:
        # 检查该细胞是否有信息
        if cell_id in adata.obs.index:
            # 获取细胞类型
            cell_type_value = adata.obs.loc[cell_id, 'annotation']
            if isinstance(cell_type_value, pd.Series):
                cell_type = cell_type_value.iloc[0]
            else:
                cell_type = cell_type_value
            
            # 获取CNs_10值
            cns_value = adata.obs.loc[cell_id, 'CNs_10']
            if isinstance(cns_value, pd.Series):
                cns_value = cns_value.iloc[0]
            
            # 检查是否为CN2
            is_cn2 = (str(cns_value) == 'CN2')
            
            # 获取坐标点
            points = group[[x_col, y_col]].values
            
            if len(points) >= 3:
                polygon = patches.Polygon(points, closed=True)
                polygons.append(polygon)
                poly_annotations.append(cell_type)
                poly_cellids.append(cell_id)
                is_cn2_list.append(is_cn2)
                matched_count += 1
                
                if is_cn2:
                    cn2_count += 1
    
    print(f"样本 {sample_name}: 成功匹配 {matched_count} 个细胞, 其中CN2细胞: {cn2_count} 个")
    return polygons, poly_annotations, poly_cellids, is_cn2_list

def create_annotation_colormap(annotations):
    """
    创建细胞类型颜色映射
    使用Tab20c调色板，提供20种不同的颜色
    """
    import matplotlib.cm as cm
    
    # 获取所有细胞类型
    unique_annotations = sorted(set(annotations))
    
    # 使用Tab20c颜色映射
    tab20c = cm.get_cmap('tab20c')
    
    # 为每个细胞类型分配颜色
    color_dict = {}
    for i, annotation in enumerate(unique_annotations):
        color_dict[annotation] = tab20c(i % 20)  # Tab20c有20种颜色，循环使用
    
    print(f"颜色映射: 已为 {len(unique_annotations)} 种细胞类型分配颜色")
    
    return color_dict, unique_annotations

def plot_cn2_spatial_maps(adata, boundaries, x_col, y_col, output_dir):
    """
    为每个样本绘制CN2细胞的空间分布图
    CN2细胞按annotation着色，非CN2细胞为灰色
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
        samples = sorted(adata.obs['sample'].unique())
    
    print(f"开始处理 {len(samples)} 个样本")
    
    # 样本级统计
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        # 创建细胞多边形
        polygons, poly_annotations, poly_cellids, is_cn2_list = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample
        )
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            sample_stats[sample] = {'total_cells': 0, 'cn2_cells': 0}
            continue
        
        # 获取CN2细胞的annotation类型
        cn2_annotations = [ann for ann, is_cn2 in zip(poly_annotations, is_cn2_list) if is_cn2]
        unique_cn2_annotations = sorted(set(cn2_annotations))
        
        if len(unique_cn2_annotations) > 0:
            # 创建颜色映射
            color_dict, _ = create_annotation_colormap(unique_cn2_annotations)
        else:
            color_dict = {}
            print(f"样本 {sample} 中没有CN2细胞")
        
        # 创建颜色列表
        colors = []
        for ann, is_cn2 in zip(poly_annotations, is_cn2_list):
            if is_cn2:
                colors.append(color_dict.get(ann, "#808080"))  # CN2细胞按类型着色
            else:
                colors.append("#D3D3D3")  # 非CN2细胞为灰色
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        
        # 使用PatchCollection批量添加多边形
        p = PatchCollection(polygons, alpha=0.8, edgecolor='white', linewidth=0.3)
        p.set_facecolor(colors)
        ax.add_collection(p)
        
        # 设置图形属性
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        ax.set_title(f'CN2细胞空间分布 - Sample: {sample}\nCN2细胞类型: {len(unique_cn2_annotations)}种', 
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
            ax.set_aspect('equal')
        
        # 添加图例 - 只显示CN2细胞类型
        legend_elements = []
        for ann in unique_cn2_annotations:
            legend_elements.append(patches.Patch(color=color_dict.get(ann, "#808080"), 
                                                label=ann, alpha=0.8))
        
        # 添加灰色块表示非CN2细胞
        legend_elements.append(patches.Patch(color="#D3D3D3", label="Non-CN2 cells", alpha=0.8))
        
        if len(legend_elements) > 0:
            ax.legend(handles=legend_elements, loc='center left', 
                     bbox_to_anchor=(1.05, 0.5), fontsize=9, framealpha=0.9,
                     ncol=1, title="Cell Types")
        
        # 添加统计信息
        cn2_count = sum(is_cn2_list)
        non_cn2_count = len(polygons) - cn2_count
        
        stats_text = f'总细胞数: {len(polygons)}\nCN2细胞: {cn2_count}\n非CN2细胞: {non_cn2_count}\nCN2细胞类型数: {len(unique_cn2_annotations)}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存图片
        sample_filename = f"cn2_spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.png"
        sample_filename_pdf = f"cn2_spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.pdf"
        output_path = output_dir / sample_filename
        output_path_pdf = output_dir / sample_filename_pdf
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {output_path}")
        plt.close(fig)
        
        # 记录统计信息
        sample_stats[sample] = {
            'total_cells': len(polygons),
            'cn2_cells': cn2_count,
            'non_cn2_cells': non_cn2_count,
            'cn2_cell_types': len(unique_cn2_annotations),
            'output_path': output_path
        }
    
    # 打印总体统计信息
    print("\n=== 处理完成 ===")
    print("样本统计信息:")
    for sample, stats in sample_stats.items():
        print(f"样本 {sample}: {stats['total_cells']} 个细胞, CN2细胞: {stats['cn2_cells']}, 非CN2细胞: {stats['non_cn2_cells']}")
    
    return sample_stats

def create_summary_statistics(adata, sample_stats, output_dir):
    """
    创建汇总统计信息和图表
    """
    output_dir = Path(output_dir)
    
    # 创建汇总数据框
    summary_data = []
    for sample, stats in sample_stats.items():
        if stats['total_cells'] > 0:
            cn2_proportion = stats['cn2_cells'] / stats['total_cells'] * 100
        else:
            cn2_proportion = 0
            
        summary_data.append({
            'sample': sample,
            'total_cells': stats['total_cells'],
            'cn2_cells': stats['cn2_cells'],
            'non_cn2_cells': stats['non_cn2_cells'],
            'cn2_proportion': cn2_proportion,
            'cn2_cell_types': stats.get('cn2_cell_types', 0)
        })
    
    if not summary_data:
        print("没有有效的统计数据")
        return
    
    summary_df = pd.DataFrame(summary_data)
    
    # 1. CN2细胞比例条形图
    plt.figure(figsize=(12, 8))
    
    # 按CN2比例排序
    summary_df = summary_df.sort_values('cn2_proportion', ascending=False)
    
    bars = plt.bar(range(len(summary_df)), summary_df['cn2_proportion'], 
                   color='skyblue', edgecolor='black')
    
    plt.title('CN2细胞比例分布 (按样本)', fontsize=16, pad=20)
    plt.xlabel('样本', fontsize=14)
    plt.ylabel('CN2细胞比例 (%)', fontsize=14)
    plt.xticks(range(len(summary_df)), summary_df['sample'], rotation=45, ha='right')
    
    # 在每个柱子上添加数值
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'cn2_proportion_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. CN2细胞数量散点图
    plt.figure(figsize=(10, 8))
    
    plt.scatter(summary_df['total_cells'], summary_df['cn2_cells'], 
                s=100, alpha=0.7, c='coral', edgecolors='black')
    
    # 添加样本标签
    for i, row in summary_df.iterrows():
        plt.annotate(row['sample'], (row['total_cells'], row['cn2_cells']),
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.title('CN2细胞数量 vs 总细胞数量', fontsize=16, pad=20)
    plt.xlabel('总细胞数量', fontsize=14)
    plt.ylabel('CN2细胞数量', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'cn2_cell_count_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 保存汇总数据
    summary_df.to_csv(output_dir / 'cn2_cell_statistics.csv', index=False)
    
    print(f"\n汇总统计已保存到: {output_dir}")
    print("\nCN2细胞统计摘要:")
    print(f"总样本数: {len(summary_df)}")
    print(f"总CN2细胞数: {summary_df['cn2_cells'].sum()}")
    print(f"平均CN2比例: {summary_df['cn2_proportion'].mean():.2f}%")
    print(f"CN2比例范围: {summary_df['cn2_proportion'].min():.2f}% - {summary_df['cn2_proportion'].max():.2f}%")

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='绘制CN2细胞空间分布图',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python cn2_spatial_plot.py -d TMA_5.h5ad -b cell_boundaries.csv -o cn2_spatial_maps
  python cn2_spatial_plot.py --data data.h5ad --boundaries boundaries.csv --output results/
        '''
    )
    
    parser.add_argument('-d', '--data', required=True, 
                       help='h5ad文件路径 (必须包含CNs_10, annotation, sample列)')
    parser.add_argument('-b', '--boundaries', required=True,
                       help='细胞边界文件路径 (CSV格式，包含cell_id和坐标)')
    parser.add_argument('-o', '--output', default='cn2_spatial_maps',
                       help='输出目录 (默认: cn2_spatial_maps)')
    
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()
    
    h5ad_path = args.data
    boundaries_path = args.boundaries
    output_dir = args.output
    
    # 检查输入文件是否存在
    for file_path, desc in [(h5ad_path, "h5ad文件"), (boundaries_path, "边界文件")]:
        if not Path(file_path).exists():
            print(f"错误: 未找到{desc}: {file_path}")
            sys.exit(1)
    
    print(f"=== CN2细胞空间分布图生成工具 ===")
    print(f"h5ad文件: {h5ad_path}")
    print(f"边界文件: {boundaries_path}")
    print(f"输出目录: {output_dir}")
    print("=" * 50)
    
    try:
        # 1. 加载和处理数据
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
            return
        
        # 3. 检查CNs_10列
        cn2_cells = adata.obs[adata.obs['CNs_10'] == 'CN2'].shape[0]
        print(f"\n=== CN2细胞统计 ===")
        print(f"CN2细胞数量: {cn2_cells}")
        print(f"CN2细胞比例: {cn2_cells/adata.n_obs*100:.2f}%")
        
        if cn2_cells == 0:
            print("警告: 没有找到CN2细胞！")
        
        # 4. 为每个样本绘制CN2细胞空间图
        sample_stats = plot_cn2_spatial_maps(adata, boundaries, x_col, y_col, output_dir)
        
        # 5. 创建汇总统计
        if sample_stats:
            create_summary_statistics(adata, sample_stats, output_dir)
        
        print(f"\n所有处理完成！图像保存在: {output_dir}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
