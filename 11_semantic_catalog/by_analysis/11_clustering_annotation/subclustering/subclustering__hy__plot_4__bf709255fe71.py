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
    output_dir="matri-prod-macro_2"           # 新输出目录
):
    """
    为每个样本创建空间原位图：
    - 展示 matri-prod CAF (THY1) 和 CTSL+ PLAUR+ Macro (ITGB2) 的表达量
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
    if 'THY1' not in adata.var_names:
        raise ValueError("var_names中未找到'THY1'")
    if 'ITGB2' not in adata.var_names:
        raise ValueError("var_names中未找到'ITGB2'")

    # ========================= 标准化 =========================
    print("正在进行文库大小归一化和log1p转换...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    thy1_expr = adata_norm[:, 'THY1'].X
    cd11c_expr = adata_norm[:, 'ITGB2'].X

    if hasattr(thy1_expr, 'toarray'):
        thy1_expr = thy1_expr.toarray().flatten()
    else:
        thy1_expr = thy1_expr.flatten()
    if hasattr(cd11c_expr, 'toarray'):
        cd11c_expr = cd11c_expr.toarray().flatten()
    else:
        cd11c_expr = cd11c_expr.flatten()

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    # 两种目标细胞类型及边框颜色
    target_edgecolors = {
        'matri-prod CAF': '#D73027',            # 红色
        'CTSL+ PLAUR+ Macro': '#33A02C'         # 蓝色
    }
    # 表达量颜色映射（低->高）
    cmap_thy1 = mcolors.LinearSegmentedColormap.from_list(
        'thy1_cmap', ['#FFFFFF', '#FC8D59']      # 白 -> 橙 (原COL4A1颜色)
    )
    cmap_cd11c = mcolors.LinearSegmentedColormap.from_list(
        'cd11c_cmap', ['#FFFFFF', '#4575B4']     # 白 -> 蓝 (原CD44颜色)
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
            macro_mask = sample_mask & (adata_norm.obs['celltype'] == 'CTSL+ PLAUR+ Macro')

            # THY1 范围（仅 matri-prod CAF）
            if caf_mask.any():
                caf_vals = thy1_expr[caf_mask]
                vmin_thy1, vmax_thy1 = np.min(caf_vals), np.max(caf_vals)
            else:
                vmin_thy1, vmax_thy1 = 0, 1

            # ITGB2 范围（仅 CTSL+ PLAUR+ Macro）
            if macro_mask.any():
                macro_vals = cd11c_expr[macro_mask]
                vmin_cd11c, vmax_cd11c = np.min(macro_vals), np.max(macro_vals)
            else:
                vmin_cd11c, vmax_cd11c = 0, 1

            norm_thy1 = mcolors.Normalize(vmin=vmin_thy1, vmax=vmax_thy1)
            norm_cd11c = mcolors.Normalize(vmin=vmin_cd11c, vmax=vmax_cd11c)

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
                    expr = thy1_expr[idx]
                    if np.isnan(expr):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((expr - vmin_thy1) / (vmax_thy1 - vmin_thy1)
                                      if vmax_thy1 > vmin_thy1 else 0, 0, 1)
                        facecolor = cmap_thy1(val)
                    edgecolor = '#D73027'
                elif cell_type == 'CTSL+ PLAUR+ Macro':
                    expr = cd11c_expr[idx]
                    if np.isnan(expr):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((expr - vmin_cd11c) / (vmax_cd11c - vmin_cd11c)
                                      if vmax_cd11c > vmin_cd11c else 0, 0, 1)
                        facecolor = cmap_cd11c(val)
                    edgecolor = '#33A02C'
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
                         'matri-prod CAF (THY1) | CTSL+ PLAUR+ Macro (ITGB2)',
                         fontsize=14, pad=20)

            legend_elements = [
                mpatches.Patch(facecolor='white', edgecolor='#D73027',
                               label='matri-prod CAF (THY1)'),
                mpatches.Patch(facecolor='white', edgecolor='#33A02C',
                               label='CTSL+ PLAUR+ Macro (ITGB2)'),
                mpatches.Patch(facecolor='white', edgecolor='lightgray',
                               label='Other Cells')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

            if caf_mask.any() and vmax_thy1 > vmin_thy1:
                sm_thy1 = plt.cm.ScalarMappable(cmap=cmap_thy1, norm=norm_thy1)
                sm_thy1.set_array([])
                cbar_thy1 = plt.colorbar(sm_thy1, ax=ax, fraction=0.046, pad=0.04)
                cbar_thy1.set_label('THY1 (log1p norm)', fontsize=9)
            if macro_mask.any() and vmax_cd11c > vmin_cd11c:
                sm_cd11c = plt.cm.ScalarMappable(cmap=cmap_cd11c, norm=norm_cd11c)
                sm_cd11c.set_array([])
                cbar_cd11c = plt.colorbar(sm_cd11c, ax=ax, fraction=0.046, pad=0.12)
                cbar_cd11c.set_label('ITGB2 (log1p norm)', fontsize=9)

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
                  f"CTSL+ PLAUR+ Macro: {macro_mask.sum()}, 总细胞: {sample_mask.sum()}")
            if caf_mask.any():
                print(f"  THY1 表达范围: {vmin_thy1:.4f} - {vmax_thy1:.4f}")
            if macro_mask.any():
                print(f"  ITGB2 表达范围: {vmin_cd11c:.4f} - {vmax_cd11c:.4f}")

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
        output_dir="matri-prod_CTSL_macro_2"
    )
