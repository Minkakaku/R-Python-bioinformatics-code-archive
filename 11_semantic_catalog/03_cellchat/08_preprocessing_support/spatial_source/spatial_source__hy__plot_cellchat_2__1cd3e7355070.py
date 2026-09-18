import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib import cm
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
    
    # 安全处理重复的cell_id
    if 'cell_id' in adata.obs.columns:
        duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep=False)
        duplicate_count = duplicate_mask.sum()
        
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

def plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir):
    """
    为每个样本绘制空间原位图：
    - EMT_Epithelium: 橙色边框，COL4A1/COL4A2/COL2A1评分（白->红）
    - Epi-CAF: 深蓝色边框，CD44表达量（白->蓝）
    - 其他: 白色边框，根据mucus_zone涂灰（inner深灰，outer浅灰）
    """
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
    
    # ========== 1. 预处理基因表达数据 ==========
    # 获取原始/归一化矩阵（优先使用 raw）
    if adata.raw is not None:
        expr_mat = adata.raw.X
        var_names = adata.raw.var_names
    else:
        expr_mat = adata.X
        var_names = adata.var_names
    obs_names = adata.obs_names
    
    # --- 计算 COL 评分 (用于 EMT_Epithelium) ---
    col_genes = ['COL4A2', 'COL4A1', 'COL2A1']
    existing_cols = [g for g in col_genes if g in var_names]
    col_score = np.zeros(adata.n_obs)
    if len(existing_cols) > 0:
        idx = [var_names.get_loc(g) for g in existing_cols]
        col_mat = expr_mat[:, idx]
        if hasattr(col_mat, 'toarray'):
            col_mat = col_mat.toarray()
        col_score = col_mat.mean(axis=1)  # 取平均值
        print(f"检测到 COL 基因: {existing_cols}，评分范围: {col_score.min():.4f} - {col_score.max():.4f}")
    else:
        print("警告: 未找到任何 COL 基因 (COL4A2, COL4A1, COL2A1)，EMT_Epithelium 将全部显示为白色")
    
    col_score = np.nan_to_num(col_score, nan=0.0)
    col_min, col_max = np.min(col_score), np.max(col_score)
    if col_max == col_min:
        col_max += 1.0
    
    # --- CD44 表达 (用于 Epi-CAF) ---
    cd44_vals = np.zeros(adata.n_obs)
    if 'CD44' in var_names:
        idx = var_names.get_loc('CD44')
        cd44_vals = expr_mat[:, idx]
        if hasattr(cd44_vals, 'toarray'):
            cd44_vals = cd44_vals.toarray().flatten()
        else:
            cd44_vals = cd44_vals.flatten()
        print(f"检测到 CD44 基因，表达范围: {cd44_vals.min():.4f} - {cd44_vals.max():.4f}")
    else:
        print("警告: CD44 基因不在数据中，Epi-CAF 将全部显示为白色")
    
    cd44_vals = np.nan_to_num(cd44_vals, nan=0.0)
    cd44_min, cd44_max = np.min(cd44_vals), np.max(cd44_vals)
    if cd44_max == cd44_min:
        cd44_max += 1.0
    
    # 构建快速查询字典 (cell_id -> value)
    col_dict = dict(zip(obs_names, col_score))
    cd44_dict = dict(zip(obs_names, cd44_vals))
    mucus_dict = adata.obs['mucus_zone'].to_dict() if 'mucus_zone' in adata.obs.columns else {}
    annot_dict = adata.obs['annotation'].to_dict()
    
    # ========== 2. 样本循环绘图 ==========
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== 处理样本: {sample} ===")
        
        sample_cell_ids = set(adata.obs[adata.obs['sample'] == sample].index)
        sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cell_ids)]
        if len(sample_boundaries) == 0:
            print(f"样本 {sample} 没有找到匹配的边界数据，跳过")
            continue
        
        cell_groups = sample_boundaries.groupby('cell_id')
        
        polygons = []
        face_colors = []
        edge_colors = []
        
        # 统计计数
        emtepi_count = 0
        epicaf_count = 0
        has_inner = False
        has_outer = False
        
        for cell_id, group in cell_groups:
            if cell_id not in annot_dict:
                continue
            
            points = group[[x_col, y_col]].values
            if len(points) < 3:
                continue
            polygon = patches.Polygon(points, closed=True)
            polygons.append(polygon)
            
            cell_type = annot_dict[cell_id]
            mucus_zone = mucus_dict.get(cell_id, None)
            
            # --- 核心颜色逻辑（根据你的新要求） ---
            if cell_type == 'EMT_Epithelium':
                emtepi_count += 1
                edge_colors.append('orange')
                # 使用 COL 评分，白到红
                val = col_dict.get(cell_id, 0.0)
                norm_val = (val - col_min) / (col_max - col_min)
                norm_val = max(0.0, min(1.0, norm_val))
                face_colors.append(cm.Reds(norm_val))
                
            elif cell_type == 'Epi-CAF':
                epicaf_count += 1
                edge_colors.append('darkblue')
                # 使用 CD44 表达，白到蓝
                val = cd44_dict.get(cell_id, 0.0)
                norm_val = (val - cd44_min) / (cd44_max - cd44_min)
                norm_val = max(0.0, min(1.0, norm_val))
                face_colors.append(cm.Blues(norm_val))
                
            else:
                # 其他细胞：只看 mucus_zone
                edge_colors.append('white')
                if mucus_zone == 'inner':
                    has_inner = True
                    face_colors.append('#333333')  # 深灰
                elif mucus_zone == 'outer':
                    has_outer = True
                    face_colors.append('#D3D3D3')  # 浅灰
                else:
                    face_colors.append('#EEEEEE')  # 极浅灰（缺失值）
        
        if len(polygons) == 0:
            print(f"样本 {sample} 没有可绘制的细胞，跳过")
            continue
        
        # --- 绘制 ---
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        p = PatchCollection(polygons, alpha=0.8, linewidth=0.3)
        p.set_facecolor(face_colors)
        p.set_edgecolor(edge_colors)
        ax.add_collection(p)
        
        # 坐标轴范围与比例
        all_vertices = np.vstack([poly.get_xy() for poly in polygons])
        x_min, x_max = all_vertices[:, 0].min(), all_vertices[:, 0].max()
        y_min, y_max = all_vertices[:, 1].min(), all_vertices[:, 1].max()
        x_margin = (x_max - x_min) * 0.05
        y_margin = (y_max - y_min) * 0.05
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_aspect('equal')
        
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        ax.set_title(f'Spatial Map - Sample: {sample}\n(EMT_Epithelium & Epi-CAF Highlighted)', fontsize=16, pad=20)
        
        # --- 动态生成图例 ---
        legend_elements = []
        if emtepi_count > 0:
            legend_elements.append(patches.Patch(facecolor=cm.Reds(0.5), edgecolor='orange', label='EMT_Epithelium (COLs)'))
        if epicaf_count > 0:
            legend_elements.append(patches.Patch(facecolor=cm.Blues(0.5), edgecolor='darkblue', label='Epi-CAF (CD44)'))
        if has_inner:
            legend_elements.append(patches.Patch(facecolor='#333333', edgecolor='white', label='Other - Inner'))
        if has_outer:
            legend_elements.append(patches.Patch(facecolor='#D3D3D3', edgecolor='white', label='Other - Outer'))
        
        if legend_elements:
            ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)
        
        # --- 统计信息 ---
        stats_text = f'Total cells: {len(polygons)}\nEMT_Epithelium: {emtepi_count}\nEpi-CAF: {epicaf_count}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存
        sample_filename = f"spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.png"
        sample_filename1 = f"spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.pdf"
        output_path = output_dir / sample_filename
        output_path1 = output_dir / sample_filename1
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path1, dpi=300, bbox_inches='tight')
        print(f"图片已保存: {output_path}")
        plt.close(fig)
        
        sample_stats[sample] = {'total_cells': len(polygons)}
    
    print("\n=== 处理完成 ===")
    return sample_stats

def main():
    # 文件路径配置
    h5ad_path = "adata_with_CNs10_TMA-5_with_patient_MZ.h5ad"
    boundaries_path = "cell_boundaries.csv"
    output_dir = "emtepi--epicaf"  # 新的输出目录
    
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
        
        # 3. 绘制空间图
        plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir)
        
        print(f"\n所有处理完成！图像保存在: {output_dir}")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
