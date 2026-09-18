import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import os
import warnings

warnings.filterwarnings('ignore')

def create_spatial_plots_subset(
    adata_path="subset_data.h5ad",
    boundaries_path="cell_boundaries.csv",
    output_dir="matri-prod-isc-macro"
):
    """
    为每个样本创建空间原位图：
    - 展示 matri-prod CAF (COL4A1)、WNT+ ISC (CD44)、CTSL+ PLAUR+ Macro (CD44) 的表达量
    - X矩阵为原始计数，内部进行文库大小归一化+log1p
    """
    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample_new', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    if 'COL4A1' not in adata.var_names:
        raise ValueError("var_names中未找到'COL4A1'")
    if 'CD44' not in adata.var_names:
        raise ValueError("var_names中未找到'CD44'")

    # ========================= 标准化 =========================
    print("正在进行文库大小归一化和log1p转换...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    col4a1_expr = adata_norm[:, 'COL4A1'].X
    cd44_expr   = adata_norm[:, 'CD44'].X

    if hasattr(col4a1_expr, 'toarray'):
        col4a1_expr = col4a1_expr.toarray().flatten()
    else:
        col4a1_expr = col4a1_expr.flatten()
    if hasattr(cd44_expr, 'toarray'):
        cd44_expr = cd44_expr.toarray().flatten()
    else:
        cd44_expr = cd44_expr.flatten()

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    # 三种目标细胞类型及新边框颜色
    target_types = {
        'matri-prod CAF': '#D73027',           # 红色调
        'WNT+ ISC': '#3F48CC',                 # 蓝紫
        'CTSL+ PLAUR+ Macro': '#33A02C'        # 蓝色
    }
    # 表达量颜色映射（低->高）
    cmap_col4a1 = mcolors.LinearSegmentedColormap.from_list(
        'col4a1_cmap', ['#FFFFFF', '#FC8D59']   # 白 -> 橙
    )
    cmap_cd44 = mcolors.LinearSegmentedColormap.from_list(
        'cd44_cmap', ['#FFFFFF', '#4575B4']     # 白 -> 蓝
    )

    unique_samples = adata_norm.obs['sample_new'].unique()
    print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")
    os.makedirs(output_dir, exist_ok=True)

    # cell_id -> 整数位置索引
    cell_id_to_idx = pd.Series(range(adata_norm.n_obs), index=adata_norm.obs['cell_id']).to_dict()

    for sample in unique_samples:
        print(f"\n处理样本: {sample}")
        try:
            sample_mask = adata_norm.obs['sample_new'] == sample
            sample_cell_ids = adata_norm.obs.loc[sample_mask, 'cell_id']

            fig, ax = plt.subplots(1, 1, figsize=(15, 12))

            caf_mask = sample_mask & (adata_norm.obs['celltype'] == 'matri-prod CAF')
            isc_mask = sample_mask & (adata_norm.obs['celltype'] == 'WNT+ ISC')
            macro_mask = sample_mask & (adata_norm.obs['celltype'] == 'CTSL+ PLAUR+ Macro')

            if caf_mask.any():
                caf_vals = col4a1_expr[caf_mask]
                vmin_col4a1, vmax_col4a1 = np.min(caf_vals), np.max(caf_vals)
            else:
                vmin_col4a1, vmax_col4a1 = 0, 1

            cd44_mask = isc_mask | macro_mask
            if cd44_mask.any():
                cd44_vals = cd44_expr[cd44_mask]
                vmin_cd44, vmax_cd44 = np.min(cd44_vals), np.max(cd44_vals)
            else:
                vmin_cd44, vmax_cd44 = 0, 1

            norm_col4a1 = mcolors.Normalize(vmin=vmin_col4a1, vmax=vmax_col4a1)
            norm_cd44   = mcolors.Normalize(vmin=vmin_cd44, vmax=vmax_cd44)

            patches, facecolors, edgecolors = [], [], []
            cells_plotted = 0

            for cell_id in sample_cell_ids:
                if cell_id not in cell_boundaries_dict or cell_id not in cell_id_to_idx:
                    continue
                vertices = cell_boundaries_dict[cell_id]
                if len(vertices) < 3:
                    continue

                idx = cell_id_to_idx[cell_id]
                cell_type = adata_norm.obs['celltype'].iloc[idx]

                if cell_type == 'matri-prod CAF':
                    expr = col4a1_expr[idx]
                    if np.isnan(expr):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((expr - vmin_col4a1) / (vmax_col4a1 - vmin_col4a1)
                                      if vmax_col4a1 > vmin_col4a1 else 0, 0, 1)
                        facecolor = cmap_col4a1(val)
                    edgecolor = '#D73027'
                elif cell_type in ['WNT+ ISC', 'CTSL+ PLAUR+ Macro']:
                    expr = cd44_expr[idx]
                    if np.isnan(expr):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((expr - vmin_cd44) / (vmax_cd44 - vmin_cd44)
                                      if vmax_cd44 > vmin_cd44 else 0, 0, 1)
                        facecolor = cmap_cd44(val)
                    edgecolor = target_types[cell_type]
                else:
                    facecolor = 'white'
                    edgecolor = 'lightgray'

                polygon = Polygon(vertices, closed=True)
                patches.append(polygon)
                facecolors.append(facecolor)
                edgecolors.append(edgecolor)
                cells_plotted += 1

            if not patches:
                print(f"样本 {sample} 没有有效边界数据，跳过")
                plt.close()
                continue

            p = PatchCollection(patches, facecolors=facecolors, edgecolors=edgecolors,
                                linewidths=0.1, alpha=0.9)
            ax.add_collection(p)
            ax.autoscale_view()
            ax.set_aspect('equal')
            ax.set_xlabel('X Position', fontsize=12)
            ax.set_ylabel('Y Position', fontsize=12)
            ax.set_title(f'Spatial Plot - Sample {sample}\n'
                         'matri-prod CAF (COL4A1) | WNT+ ISC, CTSL+ PLAUR+ Macro (CD44)',
                         fontsize=14, pad=20)

            legend_elements = [
                mpatches.Patch(facecolor='white', edgecolor='#D73027',
                               label='matri-prod CAF (COL4A1)'),
                mpatches.Patch(facecolor='white', edgecolor='#3F48CC',
                               label='WNT+ ISC (CD44)'),
                mpatches.Patch(facecolor='white', edgecolor='#33A02C',
                               label='CTSL+ PLAUR+ Macro (CD44)'),
                mpatches.Patch(facecolor='white', edgecolor='lightgray',
                               label='Other Cells')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

            if caf_mask.any() and vmax_col4a1 > vmin_col4a1:
                sm_caf = plt.cm.ScalarMappable(cmap=cmap_col4a1, norm=norm_col4a1)
                sm_caf.set_array([])
                cbar_caf = plt.colorbar(sm_caf, ax=ax, fraction=0.046, pad=0.04)
                cbar_caf.set_label('COL4A1 (log1p norm)', fontsize=9)
            if cd44_mask.any() and vmax_cd44 > vmin_cd44:
                sm_cd44 = plt.cm.ScalarMappable(cmap=cmap_cd44, norm=norm_cd44)
                sm_cd44.set_array([])
                cbar_cd44 = plt.colorbar(sm_cd44, ax=ax, fraction=0.046, pad=0.12)
                cbar_cd44.set_label('CD44 (log1p norm)', fontsize=9)

            plt.tight_layout()
            sample_str = str(sample).replace(' ', '_').replace('/', '_')
            png_path = os.path.join(output_dir, f"spatial_plot_{sample_str}.png")
            pdf_path = os.path.join(output_dir, f"spatial_plot_{sample_str}.pdf")
            plt.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"已保存: {png_path}")
            print(f"样本统计: 绘制 {cells_plotted} 个细胞")
            print(f"  matri-prod CAF: {caf_mask.sum()}, WNT+ ISC: {isc_mask.sum()}, "
                  f"CTSL+ PLAUR+ Macro: {macro_mask.sum()}, 总细胞: {sample_mask.sum()}")
            if caf_mask.any():
                print(f"  COL4A1 表达范围: {vmin_col4a1:.4f} - {vmax_col4a1:.4f}")
            if cd44_mask.any():
                print(f"  CD44 表达范围: {vmin_cd44:.4f} - {vmax_cd44:.4f}")

        except Exception as e:
            print(f"处理样本 {sample} 时出错: {e}")
            import traceback
            traceback.print_exc()
            plt.close()

    print(f"\n所有图像已保存至: {output_dir}")
    return output_dir

if __name__ == "__main__":
    create_spatial_plots_subset(
        adata_path="subset_data.h5ad",
        boundaries_path="cell_boundaries.csv",
        output_dir="matri-prod-isc-macro"
    )
