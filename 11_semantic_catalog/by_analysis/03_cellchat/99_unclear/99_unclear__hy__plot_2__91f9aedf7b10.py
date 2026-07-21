import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import os

def create_spatial_plots_per_sample(adata_path, boundaries_path, output_dir="spatial_plots_by_sample_CCL19CAF_CD8T_COL4A1"):
    """
    为每个样本创建空间原位图
    直接使用raw矩阵中已归一化的表达量，不进行额外归一化
    """
    # 读取h5ad数据
    adata = sc.read_h5ad(adata_path)
    print("h5ad数据加载完成")
    
    # 读取细胞边界文件
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    
    # 检查必要的列是否存在
    if 'sample' not in adata.obs.columns:
        raise ValueError("h5ad数据的obs中未找到'sample'列")
    if 'celltype' not in adata.obs.columns:
        raise ValueError("h5ad数据的obs中未找到'celltype'列")
    if 'cell_id' not in adata.obs.columns:
        raise ValueError("h5ad数据的obs中未找到'cell_id'列")
    
    # 确保有raw矩阵
    if adata.raw is None:
        raise ValueError("h5ad数据中没有raw矩阵")
    
    # 检查目标基因是否存在
    if 'COL4A1' not in adata.raw.var_names:
        raise ValueError("raw矩阵中未找到CCL19基因")
    if 'CD44' not in adata.raw.var_names:
        raise ValueError("raw矩阵中未找到CCR7基因")
    
    # 获取基因在raw矩阵中的索引
    cxcl12_idx = np.where(adata.raw.var_names == 'COL4A1')[0][0]
    cxcr4_idx = np.where(adata.raw.var_names == 'CD44')[0][0]
    
    # 提取基因表达量（直接使用已归一化的数据）
    if hasattr(adata.raw.X, 'toarray'):
        cxcl12_expression = adata.raw.X[:, cxcl12_idx].toarray().flatten()
        cxcr4_expression = adata.raw.X[:, cxcr4_idx].toarray().flatten()
    else:
        cxcl12_expression = adata.raw.X[:, cxcl12_idx]
        cxcr4_expression = adata.raw.X[:, cxcr4_idx]
    
    # 获取所有唯一样本
    unique_samples = adata.obs['sample'].unique()
    print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 预处理边界数据：按cell_id分组顶点
    print("正在预处理细胞边界数据...")
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices
    
    # 为每个样本创建图像
    for sample in unique_samples:
        print(f"\n处理样本: {sample}")
        
        try:
            # 创建图形
            fig, ax = plt.subplots(1, 1, figsize=(15, 12))
            
            # 筛选当前样本的细胞
            sample_mask = adata.obs['sample'] == sample
            sample_cell_ids = adata.obs[sample_mask]['cell_id']  # 使用cell_id列
            
            # 准备绘制多边形
            patches = []
            colors = []
            edge_colors = []
            
            # 获取当前样本中目标细胞类型
            sample_caf_mask = sample_mask & (adata.obs['celltype'] == 'CCL19+ CAF')
            sample_t_mask = sample_mask & (adata.obs['celltype'] == 'CD8+ T')
            
            # 创建颜色映射（白色到红色，白色到蓝色）
            caf_cmap = mcolors.LinearSegmentedColormap.from_list('white_to_red', ['white', 'red'])
            t_cmap = mcolors.LinearSegmentedColormap.from_list('white_to_blue', ['white', 'blue'])
            
            # 使用原始归一化数据的实际范围（不进行额外归一化）
            # 获取当前样本中目标细胞类型的表达量范围
            if sample_caf_mask.any():
                caf_expressions = cxcl12_expression[sample_caf_mask]
                caf_vmin, caf_vmax = np.min(caf_expressions), np.max(caf_expressions)
            else:
                caf_vmin, caf_vmax = 0, 1  # 默认范围
            
            if sample_t_mask.any():
                t_expressions = cxcr4_expression[sample_t_mask]
                t_vmin, t_vmax = np.min(t_expressions), np.max(t_expressions)
            else:
                t_vmin, t_vmax = 0, 1  # 默认范围
            
            # 遍历当前样本的所有细胞
            cells_plotted = 0
            for cell_id in sample_cell_ids:
                # 获取细胞在原始数据中的索引（通过cell_id匹配）
                cell_mask = adata.obs['cell_id'] == cell_id
                if not cell_mask.any():
                    continue
                
                original_idx = np.where(cell_mask)[0][0]
                
                # 获取细胞边界信息
                if cell_id not in cell_boundaries_dict:
                    continue
                
                vertices = cell_boundaries_dict[cell_id]
                
                if len(vertices) < 3:  # 多边形至少需要3个顶点
                    continue
                
                # 获取细胞类型
                cell_type = adata.obs.loc[cell_mask, 'celltype'].values[0]
                
                # 根据细胞类型设置颜色和边框
                if cell_type == 'CCL19+ CAF':
                    # 使用CCL19表达量设置填充色（直接使用归一化值）
                    expr_value = cxcl12_expression[original_idx]
                    # 处理可能的NaN值
                    if np.isnan(expr_value):
                        fill_color = 'white'
                    else:
                        # 将表达量直接映射到颜色，不进行0-1归一化
                        # 使用线性插值将实际表达量值映射到颜色
                        if caf_vmax > caf_vmin:  # 避免除零
                            # 直接使用实际值进行颜色映射
                            normalized_value = (expr_value - caf_vmin) / (caf_vmax - caf_vmin)
                        else:
                            normalized_value = 0
                        normalized_value = np.clip(normalized_value, 0, 1)  # 确保在有效范围内
                        fill_color = caf_cmap(normalized_value)
                    edge_color = 'red'
                    
                elif cell_type == 'CD8+ T':
                    # 使用CCR7表达量设置填充色（直接使用归一化值）
                    expr_value = cxcr4_expression[original_idx]
                    # 处理可能的NaN值
                    if np.isnan(expr_value):
                        fill_color = 'white'
                    else:
                        if t_vmax > t_vmin:  # 避免除零
                            normalized_value = (expr_value - t_vmin) / (t_vmax - t_vmin)
                        else:
                            normalized_value = 0
                        normalized_value = np.clip(normalized_value, 0, 1)
                        fill_color = t_cmap(normalized_value)
                    edge_color = 'blue'
                    
                else:
                    # 其他细胞类型
                    fill_color = 'lightgray'
                    edge_color = 'gray'
                
                # 创建多边形
                polygon = Polygon(vertices, closed=True)
                patches.append(polygon)
                colors.append(fill_color)
                edge_colors.append(edge_color)
                cells_plotted += 1
            
            # 批量绘制所有多边形以提高效率
            if patches:
                p = PatchCollection(patches, facecolors=colors, edgecolors=edge_colors,
                                  linewidths=0.3, alpha=0.9)
                ax.add_collection(p)
                
                # 设置坐标轴
                ax.autoscale_view()
                ax.set_aspect('equal')
                ax.set_xlabel('X Position', fontsize=12)
                ax.set_ylabel('Y Position', fontsize=12)
                ax.set_title(f'Spatial Plot - Sample {sample}\nCCL19+ CAF (COL4A1) and CD8+ T (CD44)', 
                           fontsize=14, pad=20)
                
                # 创建图例
                legend_elements = [
                    mpatches.Patch(facecolor='red', edgecolor='red', 
                                 label='CCL19+ CAF (COL4A11)'),
                    mpatches.Patch(facecolor='blue', edgecolor='blue', 
                                 label='CD8+ T (CD44 expr)'),
                    mpatches.Patch(facecolor='lightgray', edgecolor='gray', 
                                 label='Other Cells')
                ]
                ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
                
                # 添加颜色条（只在有对应细胞类型时添加）
                if sample_caf_mask.any() and caf_vmax > caf_vmin:
                    caf_norm = mcolors.Normalize(vmin=caf_vmin, vmax=caf_vmax)
                    caf_sm = plt.cm.ScalarMappable(cmap=caf_cmap, norm=caf_norm)
                    caf_sm.set_array([])
                    cbar_caf = plt.colorbar(caf_sm, ax=ax, fraction=0.046, pad=0.04)
                    cbar_caf.set_label('COL4A1 Expression in CCL19+ CAF', fontsize=5)
                    # 设置颜色条刻度为实际值
                    cbar_caf.set_ticks([caf_vmin, (caf_vmin + caf_vmax) / 2, caf_vmax])
                    cbar_caf.set_ticklabels([f'{caf_vmin:.2f}', f'{(caf_vmin + caf_vmax) / 2:.2f}', f'{caf_vmax:.2f}'])
                
                if sample_t_mask.any() and t_vmax > t_vmin:
                    t_norm = mcolors.Normalize(vmin=t_vmin, vmax=t_vmax)
                    t_sm = plt.cm.ScalarMappable(cmap=t_cmap, norm=t_norm)
                    t_sm.set_array([])
                    cbar_t = plt.colorbar(t_sm, ax=ax, fraction=0.046, pad=0.08)
                    cbar_t.set_label('CD44 Expression in CD8+ T', fontsize=5)
                    # 设置颜色条刻度为实际值
                    cbar_t.set_ticks([t_vmin, (t_vmin + t_vmax) / 2, t_vmax])
                    cbar_t.set_ticklabels([f'{t_vmin:.2f}', f'{(t_vmin + t_vmax) / 2:.2f}', f'{t_vmax:.2f}'])
                
                plt.tight_layout()
                
                # 保存图像
                filename = f"spatial_plot_sample_{str(sample).replace(' ', '_').replace('/', '_')}.png"
                filepath = os.path.join(output_dir, filename)
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                filename1 = f"spatial_plot_sample_{str(sample).replace(' ', '_').replace('/', '_')}.pdf"
                filepath1 = os.path.join(output_dir, filename1)
                plt.savefig(filepath1, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"已保存: {filepath}")
                print(f"样本 {sample} 统计: 绘制了 {cells_plotted} 个细胞")
                
                # 打印细胞类型统计信息
                n_caf = sample_caf_mask.sum()
                n_t = sample_t_mask.sum()
                n_total = sample_mask.sum()
                print(f"CCL19+ CAF: {n_caf}, CD8+ T: {n_t}, 总细胞数: {n_total}")
                
                # 打印表达量范围信息
                if sample_caf_mask.any():
                    print(f"COL4A1表达量范围: {caf_vmin:.4f} - {caf_vmax:.4f}")
                if sample_t_mask.any():
                    print(f"CD44表达量范围: {t_vmin:.4f} - {t_vmax:.4f}")
                
            else:
                print(f"样本 {sample} 没有有效的细胞边界数据，跳过该样本")
                plt.close()
                
        except Exception as e:
            print(f"处理样本 {sample} 时出错: {str(e)}")
            import traceback
            traceback.print_exc()
            plt.close()
    
    print(f"\n所有样本的空间原位图已保存到目录: {output_dir}")
    return output_dir

# 使用示例
if __name__ == "__main__":
    # 设置文件路径
    adata_path = "merged_single_cell_data_modified.h5ad"
    boundaries_path = "cell_boundaries.csv"
    
    # 创建空间图
    output_dir = create_spatial_plots_per_sample(adata_path, boundaries_path)
    
    print(f"处理完成！图像保存在: {output_dir}")
