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
    output_dir="matri-prod_CTSL_macro_mix"   # 新输出文件夹
):
    """
    为每个样本创建空间原位图：
    - 展示 matri-prod CAF (THY1) 和 CTSL+ PLAUR+ Macro (sqrt(ITGAX * ITGB2)) 的表达量
    - X矩阵为原始计数，内部进行文库大小归一化+log1p
    """
    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample_new', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    # 检查三个基因
    for gene in ['THY1', 'ITGAX', 'ITGB2']:
        if gene not in adata.var_names:
            raise ValueError(f"var_names中未找到'{gene}'")

    # ========================= 标准化 =========================
    print("正在进行文库大小归一化和log1p转换...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # 提取各基因表达（稠密数组）
    thy1_expr = adata_norm[:, 'THY1'].X
    itgax_expr = adata_norm[:, 'ITGAX'].X
    itgb2_expr = adata_norm[:, 'ITGB2'].X

    def to_dense(arr):
        if hasattr(arr, 'toarray'):
            return arr.toarray().flatten()
        return arr.flatten()

    thy1_expr = to_dense(thy1_expr)
    itgax_expr = to_dense(itgax_expr)
    itgb2_expr = to_dense(itgb2_expr)

    # 计算 CTSL+ PLAUR+ Macro 的受体得分：sqrt(ITGAX * ITGB2)
    receptor_score = np.sqrt(itgax_expr * itgb2_expr)

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    target_edgecolors = {
        'matri-prod CAF': '#D73027',
        'CTSL+ PLAUR+ Macro': '#33A02C'
    }
    cmap_thy1 = mcolors.LinearSegmentedColormap.from_list(
        'thy1_cmap', ['#FFFFFF', '#FC8D59']
    )
    cmap_receptor = mcolors.LinearSegmentedColormap.from_list(
        'receptor_cmap', ['#FFFFFF', '#4575B4']   # 使用原 CD44 蓝色渐变
    )

    unique_samples = adata_norm.obs['sample_new'].unique()
    print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")
    os.makedirs(output_dir, exist_ok=True)

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

            # 受体得分范围（仅 CTSL+ PLAUR+ Macro）
            if macro_mask.any():
                macro_vals = receptor_score[macro_mask]
                vmin_receptor, vmax_receptor = np.min(macro_vals), np.max(macro_vals)
            else:
                vmin_receptor, vmax_receptor = 0, 1

            norm_thy1 = mcolors.Normalize(vmin=vmin_thy1, vmax=vmax_thy1)
            norm_receptor = mcolors.Normalize(vmin=vmin_receptor, vmax=vmax_receptor)

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
                    score = receptor_score[idx]
                    if np.isnan(score):
                        facecolor = '#FFFFFF'
                    else:
                        val = np.clip((score - vmin_receptor) / (vmax_receptor - vmin_receptor)
                                      if vmax_receptor > vmin_receptor else 0, 0, 1)
                        facecolor = cmap_receptor(val)
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
                         'matri-prod CAF (THY1) | CTSL+ PLAUR+ Macro (Receptor Score)',
                         fontsize=14, pad=20)

            legend_elements = [
                mpatches.Patch(facecolor='white', edgecolor='#D73027',
                               label='matri-prod CAF (THY1)'),
                mpatches.Patch(facecolor='white', edgecolor='#33A02C',
                               label='CTSL+ PLAUR+ Macro (sqrt(ITGAX*ITGB2))'),
                mpatches.Patch(facecolor='white', edgecolor='lightgray',
                               label='Other Cells')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

            if caf_mask.any() and vmax_thy1 > vmin_thy1:
                sm_thy1 = plt.cm.ScalarMappable(cmap=cmap_thy1, norm=norm_thy1)
                sm_thy1.set_array([])
                cbar_thy1 = plt.colorbar(sm_thy1, ax=ax, fraction=0.046, pad=0.04)
                cbar_thy1.set_label('THY1 (log1p norm)', fontsize=9)
            if macro_mask.any() and vmax_receptor > vmin_receptor:
                sm_receptor = plt.cm.ScalarMappable(cmap=cmap_receptor, norm=norm_receptor)
                sm_receptor.set_array([])
                cbar_receptor = plt.colorbar(sm_receptor, ax=ax, fraction=0.046, pad=0.12)
                cbar_receptor.set_label('Receptor Score (sqrt(ITGAX*ITGB2))', fontsize=9)

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
                print(f"  Receptor Score 范围: {vmin_receptor:.4f} - {vmax_receptor:.4f}")

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
        output_dir="matri-prod_CTSL_macro_mix"
    )
