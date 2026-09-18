import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import os

# 读取h5ad数据
adata = sc.read_h5ad("merged_single_cell_data_modified.h5ad")

# 读取细胞边界文件
cell_boundaries = pd.read_csv("cell_boundaries.csv")

# 确保有raw矩阵
if adata.raw is None:
    raise ValueError("h5ad数据中没有raw矩阵，请确保数据包含原始计数")

# 检查是否有样本列
if 'sample' not in adata.obs.columns:
    raise ValueError("h5ad数据的obs中未找到'sample'列，请确保数据包含样本信息")

# 获取所有唯一样本
unique_samples = adata.obs['sample'].unique()
print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")

# 获取基因在raw矩阵中的索引
if 'CXCL12' not in adata.raw.var_names:
    raise ValueError("raw矩阵中未找到CXCL12基因")
if 'CXCR4' not in adata.raw.var_names:
    raise ValueError("raw矩阵中未找到CXCR4基因")

cxcl12_idx = np.where(adata.raw.var_names == 'CXCL12')[0][0]
cxcr4_idx = np.where(adata.raw.var_names == 'CXCR4')[0][0]

# 提取基因表达量
cxcl12_expression = adata.raw.X[:, cxcl12_idx].toarray().flatten() if hasattr(adata.raw.X, 'toarray') else adata.raw.X[:, cxcl12_idx]
cxcr4_expression = adata.raw.X[:, cxcr4_idx].toarray().flatten() if hasattr(adata.raw.X, 'toarray') else adata.raw.X[:, cxcr4_idx]

# 创建输出目录
output_dir = "spatial_plots_by_sample"
os.makedirs(output_dir, exist_ok=True)

# 为每个样本创建图像[7](@ref)
for sample in unique_samples:
    print(f"处理样本: {sample}")
    
    # 筛选当前样本的细胞
    sample_mask = adata.obs['sample'] == sample
    sample_cell_ids = adata.obs[sample_mask].index
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    # 创建颜色映射
    caf_cmap = mcolors.LinearSegmentedColormap.from_list('white_to_red', ['white', 'red'])
    t_cmap = mcolors.LinearSegmentedColormap.from_list('white_to_blue', ['white', 'blue'])
    
    # 准备绘制多边形
    patches = []
    colors = []
    edge_colors = []
    
    # 遍历当前样本的所有细胞[7](@ref)
    for cell_id in sample_cell_ids:
        # 获取细胞在原始数据中的索引
        original_idx = list(adata.obs_names).index(cell_id)
        
        # 获取细胞边界信息
        cell_boundary = cell_boundaries[cell_boundaries['cell_id'] == cell_id]
        
        if cell_boundary.empty:
            continue
            
        # 获取细胞类型
        cell_type = adata.obs.loc[cell_id, 'celltype']
        
        # 根据细胞类型设置颜色和边框
        if cell_type == 'CCL19+ CAF':
            # 使用CXCL12表达量设置填充色
            # 使用样本内的归一化
            sample_caf_mask = sample_mask & (adata.obs['celltype'] == 'CCL19+ CAF')
            sample_caf_expressions = cxcl12_expression[sample_caf_mask]
            caf_norm = mcolors.Normalize(vmin=np.min(sample_caf_expressions), 
                                        vmax=np.max(sample_caf_expressions))
            fill_color = caf_cmap(caf_norm(cxcl12_expression[original_idx]))
            edge_color = 'red'
        elif cell_type == 'CD4+ T':
            # 使用CXCR4表达量设置填充色
            # 使用样本内的归一化
            sample_t_mask = sample_mask & (adata.obs['celltype'] == 'CD4+ T')
            sample_t_expressions = cxcr4_expression[sample_t_mask]
            t_norm = mcolors.Normalize(vmin=np.min(sample_t_expressions), 
                                     vmax=np.max(sample_t_expressions))
            fill_color = t_cmap(t_norm(cxcr4_expression[original_idx]))
            edge_color = 'blue'
        else:
            # 其他细胞类型
            fill_color = 'lightgray'
            edge_color = 'gray'
        
        # 根据实际边界文件格式调整下面的代码
        # 假设边界文件包含x,y坐标序列
        try:
            # 尝试不同的列名可能性
            if 'x_coordinates' in cell_boundary.columns and 'y_coordinates' in cell_boundary.columns:
                x_coords = eval(cell_boundary['x_coordinates'].values[0])  # 如果是字符串格式
                y_coords = eval(cell_boundary['y_coordinates'].values[0])
            elif 'x' in cell_boundary.columns and 'y' in cell_boundary.columns:
                x_coords = cell_boundary['x'].values
                y_coords = cell_boundary['y'].values
            else:
                # 尝试使用前两列作为坐标
                x_coords = cell_boundary.iloc[:, 0].values
                y_coords = cell_boundary.iloc[:, 1].values
        except:
            print(f"警告: 无法解析细胞 {cell_id} 的边界坐标，跳过该细胞")
            continue
        
        # 创建多边形
        polygon = Polygon(list(zip(x_coords, y_coords)), closed=True)
        patches.append(polygon)
        colors.append(fill_color)
        edge_colors.append(edge_color)
    
    # 批量绘制所有多边形以提高效率
    if patches:  # 确保有有效的多边形
        p = PatchCollection(patches, facecolors=colors, edgecolors=edge_colors, linewidths=1.2)
        ax.add_collection(p)
        
        # 设置坐标轴
        ax.autoscale_view()
        ax.set_aspect('equal')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')
        ax.set_title(f'Spatial Plot - Sample {sample}\nCCL19+ CAF and CD4+ T Cells')
        
        # 创建图例
        legend_elements = [
            mpatches.Patch(facecolor='red', edgecolor='red', label='CCL19+ CAF (CXCL12 expr)'),
            mpatches.Patch(facecolor='blue', edgecolor='blue', label='CD4+ T (CXCR4 expr)'),
            mpatches.Patch(facecolor='lightgray', edgecolor='gray', label='Other Cells')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # 添加颜色条
        if 'sample_caf_expressions' in locals() and len(sample_caf_expressions) > 0:
            caf_sm = plt.cm.ScalarMappable(cmap=caf_cmap, norm=caf_norm)
            caf_sm.set_array([])
            cbar_caf = plt.colorbar(caf_sm, ax=ax, fraction=0.046, pad=0.04)
            cbar_caf.set_label('CXCL12 Expression in CCL19+ CAF')
        
        if 'sample_t_expressions' in locals() and len(sample_t_expressions) > 0:
            t_sm = plt.cm.ScalarMappable(cmap=t_cmap, norm=t_norm)
            t_sm.set_array([])
            cbar_t = plt.colorbar(t_sm, ax=ax, fraction=0.046, pad=0.08)
            cbar_t.set_label('CXCR4 Expression in CD4+ T')
        
        plt.tight_layout()
        
        # 保存图像
        filename = f"spatial_plot_sample_{sample}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()  # 关闭图形以释放内存
        print(f"已保存: {filepath}")
    else:
        print(f"样本 {sample} 没有有效的细胞边界数据，跳过该样本")

print(f"\n所有样本的空间原位图已保存到目录: {output_dir}")
print(f"共处理了 {len(unique_samples)} 个样本")
