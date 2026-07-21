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
    output_dir="matri-prod_wnt_isc"            # 新输出目录
):
    """
    为每个样本创建空间原位图：
    - 展示 matri-prod CAF (THBS1) 和 WNT+ ISC (CD47) 的表达量
    - X矩阵为原始计数，内部进行文库大小归一化+log1p
    """
    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample_new', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    # 检查基因
    if 'THBS1' not in adata.var_names:
        raise ValueError("var_names中未找到'THBS1'")
    if 'CD47' not in adata.var_names:
        raise ValueError("var_names中未找到'CD47'")

    # ========================= 标准化 =========================
    print("正在进行文库大小归一化和log1p转换...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    thbs1_expr = adata_norm[:, 'THBS1'].X
    cd47_expr  = adata_norm[:, 'CD47'].X

    if hasattr(thbs1_expr, 'toarray'):
        thbs1_expr = thbs1_expr.toarray().flatten()
    else:
        thbs1_expr = thbs1_expr.flatten()
    if hasattr(cd47_expr, 'toarray'):
        cd47_expr = cd47_expr.toarray().flatten()
    else:
        cd47_expr = cd47_expr.flatten()

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    # 边框颜色
    target_edgecolors = {
        'matri-prod CAF': '#D73027',            # 红色
        'WNT+ ISC': '#3F48CC'                   # 蓝紫色
    }
    # 表达量颜色映射（低->高）
    cmap_thbs1 = mcolors.LinearSegmentedColormap.from_list(
        'thbs1_cmap', ['#FFFFFF', '#FC8D59']     # 白 -> 橙
    )
    cmap_cd47 = mcolors.LinearSegmentedColormap.from_list(
        'cd47_cmap', ['#FFFFFF', '#4575B4']      # 白 -> 蓝
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

            # THBS1 范围（仅 matri-prod CAF）
            if caf_mask.any():
                caf_vals = thbs1_expr[caf_mask]
                vmin_thbs1, vmax_thbs1 = np.min(caf_vals), np.max(caf_vals)
            else:
                vmin_thbs1, vmax_thbs1 = 0, 1

            # CD47 范围（仅 WNT+ ISC）
            if isc_mask.any():
                isc_vals = cd47_expr[isc_mask]
                vmin_cd47, vmax_cd47 = np.min(isc_vals), np.max(isc_vals)
            else:
                vmin_cd47, vmax_cd47 = 0, 1

            norm_thbs1 = mcolors.Normalize(vmin=vmin_thbs1, vmax=vmax_thbs1)
            norm_cd47  = mcolors.Normalize(vmin=vmin_cd47, vmax=vmax_cd47)

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
                    expr = thbs1_expr[idx]
                    if np.isnan(expr):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((expr - vmin_thbs1) / (vmax_thbs1 - vmin_thbs1)
                                      if vmax_thbs1 > vmin_thbs1 else 0, 0, 1)
                        facecolor = cmap_thbs1(val)
                    edgecolor = '#D73027'
                elif cell_type == 'WNT+ ISC':
                    expr = cd47_expr[idx]
                    if np.isnan(expr):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((expr - vmin_cd47) / (vmax_cd47 - vmin_cd47)
                                      if vmax_cd47 > vmin_cd47 else 0, 0, 1)
                        facecolor = cmap_cd47(val)
                    edgecolor = '#3F48CC'
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
                         'matri-prod CAF (THBS1) | WNT+ ISC (CD47)',
                         fontsize=14, pad=20)

            legend_elements = [
                mpatches.Patch(facecolor='white', edgecolor='#D73027',
                               label='matri-prod CAF (THBS1)'),
                mpatches.Patch(facecolor='white', edgecolor='#3F48CC',
                               label='WNT+ ISC (CD47)'),
                mpatches.Patch(facecolor='white', edgecolor='lightgray',
                               label='Other Cells')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

            if caf_mask.any() and vmax_thbs1 > vmin_thbs1:
                sm_thbs1 = plt.cm.ScalarMappable(cmap=cmap_thbs1, norm=norm_thbs1)
                sm_thbs1.set_array([])
                cbar_thbs1 = plt.colorbar(sm_thbs1, ax=ax, fraction=0.046, pad=0.04)
                cbar_thbs1.set_label('THBS1 (log1p norm)', fontsize=9)
            if isc_mask.any() and vmax_cd47 > vmin_cd47:
                sm_cd47 = plt.cm.ScalarMappable(cmap=cmap_cd47, norm=norm_cd47)
                sm_cd47.set_array([])
                cbar_cd47 = plt.colorbar(sm_cd47, ax=ax, fraction=0.046, pad=0.12)
                cbar_cd47.set_label('CD47 (log1p norm)', fontsize=9)

            plt.tight_layout()
            sample_str = str(sample).replace(' ', '_').replace('/', '_')
            png_path = os.path.join(output_dir, f"spatial_plot_{sample_str}.png")
            pdf_path = os.path.join(output_dir, f"spatial_plot_{sample_str}.pdf")
            plt.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"已保存: {png_path}")
            print(f"样本统计: 绘制 {cells_plotted} 个细胞")
            print(f"  matri-prod CAF: {caf_mask.sum()}, "
                  f"WNT+ ISC: {isc_mask.sum()}, 总细胞: {sample_mask.sum()}")
            if caf_mask.any():
                print(f"  THBS1 表达范围: {vmin_thbs1:.4f} - {vmax_thbs1:.4f}")
            if isc_mask.any():
                print(f"  CD47 表达范围: {vmin_cd47:.4f} - {vmax_cd47:.4f}")

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
        output_dir="matri-prod_wnt_isc"          # 保存路径
    )
