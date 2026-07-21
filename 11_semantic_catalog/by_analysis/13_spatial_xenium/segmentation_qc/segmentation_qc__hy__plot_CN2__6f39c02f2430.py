import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path
import matplotlib.cm as cm
import colorsys
from matplotlib.colors import to_hex

def load_and_process_data(h5ad_path, boundaries_path):
    """
    加载并处理h5ad文件和边界数据，安全处理重复的cell_id
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

def get_top_annotations_for_sample(adata, sample_name, top_n=10):
    """
    获取指定样本中CN2细胞的前N个最常见的annotation类型
    """
    # 筛选当前样本的CN2细胞
    sample_cn2_mask = (adata.obs['sample'] == sample_name) & (adata.obs['CNs_10'] == 'CN2')
    sample_cn2_data = adata.obs[sample_cn2_mask]
    
    if len(sample_cn2_data) == 0:
        print(f"样本 {sample_name} 中没有CN2细胞")
        return [], pd.Series(dtype=int)
    
    # 统计每个annotation类型的数量
    annotation_counts = sample_cn2_data['annotation'].value_counts()
    
    # 取前N个
    top_annotations = annotation_counts.head(top_n).index.tolist()
    
    print(f"样本 {sample_name} 中CN2细胞的前{top_n}个常见annotation类型:")
    for i, (annotation, count) in enumerate(annotation_counts.head(top_n).items(), 1):
        print(f"  {i}. {annotation}: {count} 个细胞")
    
    return top_annotations, annotation_counts

def create_colors_for_annotations(top_annotations, exclude_colors=None):
    """
    为annotation类型创建颜色映射，排除指定的颜色
    """
    if exclude_colors is None:
        exclude_colors = ["#ADD8E6", "#D3D3D3"]  # 默认排除淡蓝和淡灰
    
    # 使用tab20c颜色映射
    tab20c = cm.get_cmap('tab20c')
    
    # 生成颜色列表
    colors = []
    used_colors = set()
    
    for i in range(len(top_annotations)):
        # 从tab20c获取颜色
        color = tab20c(i % 20)
        hex_color = to_hex(color)
        
        # 如果颜色在被排除的颜色列表中，尝试不同的索引
        if hex_color.upper() in [c.upper() for c in exclude_colors]:
            # 尝试使用不同的索引
            for j in range(20):
                alt_color = tab20c((i + j + 5) % 20)
                alt_hex = to_hex(alt_color)
                if alt_hex.upper() not in [c.upper() for c in exclude_colors] and alt_hex not in used_colors:
                    color = alt_color
                    hex_color = alt_hex
                    break
        
        colors.append(hex_color)
        used_colors.add(hex_color)
    
    # 创建颜色字典
    color_dict = {}
    for i, annotation in enumerate(top_annotations):
        color_dict[annotation] = colors[i]
    
    # 为特殊类别添加颜色
    color_dict['others-CN2'] = "#D3D3D3"  # 淡灰色
    color_dict['others'] = "#ADD8E6"  # 淡蓝色
    
    return color_dict

def process_cell_annotations(adata, cell_ids, sample_name, top_annotations):
    """
    处理细胞的annotation标签
    """
    processed_annotations = []
    
    for cell_id in cell_ids:
        if cell_id in adata.obs.index:
            # 获取原始annotation
            orig_annotation = adata.obs.loc[cell_id, 'annotation']
            if isinstance(orig_annotation, pd.Series):
                orig_annotation = orig_annotation.iloc[0]
            
            # 获取CNs_10值
            cns_value = adata.obs.loc[cell_id, 'CNs_10']
            if isinstance(cns_value, pd.Series):
                cns_value = cns_value.iloc[0]
            
            # 处理annotation
            if str(cns_value) == 'CN2':
                if orig_annotation in top_annotations:
                    processed_annotations.append(orig_annotation)
                else:
                    processed_annotations.append('others-CN2')
            else:
                processed_annotations.append('others')
        else:
            processed_annotations.append('others')
    
    return processed_annotations

def create_cell_polygons(boundaries, adata, x_col, y_col, sample_name, top_annotations):
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
    poly_cellids = []
    
    matched_count = 0
    for cell_id, group in cell_groups:
        # 检查该细胞是否有信息
        if cell_id in adata.obs.index:
            # 获取坐标点
            points = group[[x_col, y_col]].values
            
            if len(points) >= 3:  # 至少需要3个点来创建多边形
                # 创建多边形对象
                polygon = patches.Polygon(points, closed=True)
                polygons.append(polygon)
                poly_cellids.append(cell_id)
                matched_count += 1
    
    print(f"样本 {sample_name}: 成功匹配 {matched_count} 个细胞")
    
    # 处理细胞的annotation标签
    processed_annotations = process_cell_annotations(adata, poly_cellids, sample_name, top_annotations)
    
    return polygons, processed_annotations, poly_cellids

def plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir, top_n=10):
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
        samples = sorted(adata.obs['sample'].unique())
    
    print(f"开始处理 {len(samples)} 个样本")
    
    # 样本级统计
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        # 获取当前样本CN2细胞的前N个annotation类型
        top_annotations, annotation_counts = get_top_annotations_for_sample(adata, sample, top_n)
        
        # 创建颜色映射
        color_dict = create_colors_for_annotations(top_annotations)
        
        # 创建细胞多边形
        polygons, processed_annotations, poly_cellids = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample, top_annotations
        )
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            sample_stats[sample] = {'total_cells': 0, 'cn2_cells': 0}
            continue
        
        # 统计各类别数量
        annotation_counts_in_plot = pd.Series(processed_annotations).value_counts()
        
        # 创建图形
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        
        # 创建颜色列表
        colors = [color_dict.get(ann, "#808080") for ann in processed_annotations]
        
        # 使用PatchCollection批量添加多边形
        p = PatchCollection(polygons, alpha=0.8, edgecolor='white', linewidth=0.3)
        p.set_facecolor(colors)
        ax.add_collection(p)
        
        # 设置图形属性
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        
        # 获取统计信息
        cn2_cells = sum(1 for ann in processed_annotations if ann != 'others')
        top_annotation_count = sum(annotation_counts_in_plot.get(ann, 0) for ann in top_annotations)
        others_cn2_count = annotation_counts_in_plot.get('others-CN2', 0)
        others_count = annotation_counts_in_plot.get('others', 0)
        
        title = f'CN2细胞空间分布 - Sample: {sample}\n'
        title += f'CN2细胞总数: {cn2_cells} | Top{top_n}类型: {top_annotation_count} | Others-CN2: {others_cn2_count} | Others: {others_count}'
        ax.set_title(title, fontsize=16, pad=20)
        
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
        
        # 创建图例元素
        legend_elements = []
        
        # 添加top annotations
        for annotation in top_annotations:
            if annotation in color_dict:
                legend_elements.append(patches.Patch(
                    color=color_dict[annotation], 
                    label=f'{annotation} ({annotation_counts_in_plot.get(annotation, 0)})',
                    alpha=0.8
                ))
        
        # 添加others-CN2
        if 'others-CN2' in annotation_counts_in_plot:
            legend_elements.append(patches.Patch(
                color=color_dict['others-CN2'],
                label=f'others-CN2 ({annotation_counts_in_plot.get("others-CN2", 0)})',
                alpha=0.8
            ))
        
        # 添加others
        if 'others' in annotation_counts_in_plot:
            legend_elements.append(patches.Patch(
                color=color_dict['others'],
                label=f'others ({annotation_counts_in_plot.get("others", 0)})',
                alpha=0.8
            ))
        
        if len(legend_elements) > 0:
            # 根据图例项数量调整列数
            ncols = 1
            if len(legend_elements) > 15:
                ncols = 2
            if len(legend_elements) > 30:
                ncols = 3
            
            ax.legend(handles=legend_elements, loc='center left', 
                     bbox_to_anchor=(1.05, 0.5), fontsize=9, framealpha=0.9,
                     ncol=ncols, title="Cell Categories")
        
        # 添加统计信息
        stats_text = f'总细胞数: {len(polygons)}\nCN2细胞: {cn2_cells}\nTop{top_n}类型: {top_annotation_count}\nOthers-CN2: {others_cn2_count}\nOthers: {others_count}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存图片
        sample_filename = f"cn2_top{top_n}_spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.png"
        sample_filename_pdf = f"cn2_top{top_n}_spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.pdf"
        output_path = output_dir / sample_filename
        output_path_pdf = output_dir / sample_filename_pdf
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {output_path}")
        plt.close(fig)
        
        # 记录统计信息
        sample_stats[sample] = {
            'total_cells': len(polygons),
            'cn2_cells': cn2_cells,
            'top_annotation_count': top_annotation_count,
            'others_cn2_count': others_cn2_count,
            'others_count': others_count,
            'top_annotations': top_annotations,
            'output_path': output_path
        }
    
    # 打印总体统计信息
    print("\n=== 处理完成 ===")
    print("样本统计信息:")
    for sample, stats in sample_stats.items():
        print(f"样本 {sample}: {stats['total_cells']} 个细胞, CN2细胞: {stats['cn2_cells']}, "
              f"Top{top_n}类型: {stats['top_annotation_count']}, Others-CN2: {stats['others_cn2_count']}, "
              f"Others: {stats['others_count']}")
    
    return sample_stats

def create_summary_statistics(adata, sample_stats, output_dir, top_n=10):
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
            'cn2_proportion': cn2_proportion,
            'top_annotation_count': stats.get('top_annotation_count', 0),
            'others_cn2_count': stats.get('others_cn2_count', 0),
            'others_count': stats.get('others_count', 0)
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
    
    plt.title(f'CN2细胞比例分布 (按样本, Top{top_n}分析)', fontsize=16, pad=20)
    plt.xlabel('样本', fontsize=14)
    plt.ylabel('CN2细胞比例 (%)', fontsize=14)
    plt.xticks(range(len(summary_df)), summary_df['sample'], rotation=45, ha='right')
    
    # 在每个柱子上添加数值
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'cn2_top{top_n}_proportion_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. 各类别细胞数量堆叠条形图
    plt.figure(figsize=(12, 8))
    
    categories = ['top_annotation_count', 'others_cn2_count', 'others_count']
    colors = ['#2ca02c', '#D3D3D3', '#ADD8E6']  # 绿色, 淡灰色, 淡蓝色
    labels = [f'Top{top_n} CN2', 'Others-CN2', 'Others']
    
    bottom = np.zeros(len(summary_df))
    for i, (category, color, label) in enumerate(zip(categories, colors, labels)):
        values = summary_df[category].values
        plt.bar(range(len(summary_df)), values, bottom=bottom, 
                color=color, edgecolor='black', label=label)
        bottom += values
    
    plt.title(f'细胞类别分布 (按样本, Top{top_n}分析)', fontsize=16, pad=20)
    plt.xlabel('样本', fontsize=14)
    plt.ylabel('细胞数量', fontsize=14)
    plt.xticks(range(len(summary_df)), summary_df['sample'], rotation=45, ha='right')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / f'cn2_top{top_n}_category_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 保存汇总数据
    summary_df.to_csv(output_dir / f'cn2_top{top_n}_statistics.csv', index=False)
    
    print(f"\n汇总统计已保存到: {output_dir}")
    print(f"\nCN2细胞统计摘要 (Top{top_n}分析):")
    print(f"总样本数: {len(summary_df)}")
    print(f"总CN2细胞数: {summary_df['cn2_cells'].sum()}")
    print(f"平均CN2比例: {summary_df['cn2_proportion'].mean():.2f}%")
    print(f"CN2比例范围: {summary_df['cn2_proportion'].min():.2f}% - {summary_df['cn2_proportion'].max():.2f}%")

def main():
    # 文件路径配置
    h5ad_path = "TMA_5.h5ad"  # 替换为您的h5ad文件路径
    boundaries_path = "cell_boundaries.csv"  # 替换为您的边界文件路径
    output_dir = f"spatial_CN2_top10"
    top_n = 10  # 选择前10个最常见的annotation类型
    
    print(f"=== CN2细胞空间分布图生成工具 (Top{top_n}分析) ===")
    print(f"h5ad文件: {h5ad_path}")
    print(f"边界文件: {boundaries_path}")
    print(f"输出目录: {output_dir}")
    print(f"分析Top{top_n}个最常见的CN2细胞annotation类型")
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
            print("h5ad中的前5个细胞ID:", list(h5ad_cells)[:5])
            print("边界文件中的前5个细胞ID:", list(boundary_cells)[:5])
            return
        
        # 3. 检查CN2细胞
        cn2_cells = adata.obs[adata.obs['CNs_10'] == 'CN2'].shape[0]
        print(f"\n=== CN2细胞统计 ===")
        print(f"CN2细胞总数: {cn2_cells}")
        print(f"CN2细胞比例: {cn2_cells/adata.n_obs*100:.2f}%")
        
        if cn2_cells == 0:
            print("警告: 没有找到CN2细胞！")
        
        # 4. 为每个样本绘制CN2细胞空间图
        sample_stats = plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir, top_n)
        
        # 5. 创建汇总统计
        if sample_stats:
            create_summary_statistics(adata, sample_stats, output_dir, top_n)
        
        print(f"\n所有处理完成！图像保存在: {output_dir}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
