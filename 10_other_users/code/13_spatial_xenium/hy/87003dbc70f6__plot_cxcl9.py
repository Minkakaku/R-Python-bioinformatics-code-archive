import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import os
import warnings

warnings.filterwarnings('ignore')

def create_spatial_plot_CXCL9(
    adata_path="TMA2.h5ad",
    boundaries_path="cell_boundaries.csv",
    output_dir="sp_CXCL9"
):
    """
    为每个样本创建基因 CXCL9 的空间原位图：
    - 仅基因：CXCL9
    - 不对细胞类型筛选，所有细胞均绘制
    - 表达量经文库归一化 + log1p 后，用白→红渐变填充
    - 所有细胞边框均为浅灰色
    """
    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    # 检查必要的列（样本列改为'sample'）
    for col in ['sample', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    target_gene = 'CXCL9'
    if target_gene not in adata.var_names:
        raise ValueError(f"var_names中未找到基因'{target_gene}'")

    # ========================= 标准化 =========================
    print("正在进行文库大小归一化和log1p转换...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # 提取CXCL9的表达量（稠密数组）
    expr = adata_norm[:, target_gene].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    else:
        expr = expr.flatten()
    gene_expr = expr

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    cmap_red = mcolors.LinearSegmentedColormap.from_list(
        'white_red', ['#FFFFFF', '#FF0000']
    )

    unique_samples = adata_norm.obs['sample'].unique()
    print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")
    os.makedirs(output_dir, exist_ok=True)

    # cell_id -> 整数位置索引
    cell_id_to_idx = pd.Series(range(adata_norm.n_obs), index=adata_norm.obs['cell_id']).to_dict()

    for sample in unique_samples:
        sample_mask = adata_norm.obs['sample'] == sample
        sample_cell_ids = adata_norm.obs.loc[sample_mask, 'cell_id']
        sample_int_idx = np.where(sample_mask)[0]

        print(f"\n处理样本: {sample}, 基因: {target_gene}")
        try:
            fig, ax = plt.subplots(1, 1, figsize=(15, 12))

            # 获取该样本中该基因的表达值范围
            expr_vals = gene_expr[sample_int_idx]
            if len(expr_vals) > 0:
                vmin, vmax = np.min(expr_vals), np.max(expr_vals)
            else:
                vmin, vmax = 0, 1
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

            patches, facecolors, edgecolors = [], [], []
            cells_plotted = 0

            for cell_id in sample_cell_ids:
                if cell_id not in cell_boundaries_dict or cell_id not in cell_id_to_idx:
                    continue
                vertices = cell_boundaries_dict[cell_id]
                if len(vertices) < 3:
                    continue

                idx = cell_id_to_idx[cell_id]
                expr_val = gene_expr[idx]

                if np.isnan(expr_val):
                    facecolor = '#FFFFFF'
                else:
                    val = np.clip((expr_val - vmin) / (vmax - vmin) if vmax > vmin else 0, 0, 1)
                    facecolor = cmap_red(val)
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
            ax.set_title(f'Spatial Plot - Sample {sample} - {target_gene}\nAll cell types',
                         fontsize=14, pad=20)

            if vmax > vmin:
                sm = plt.cm.ScalarMappable(cmap=cmap_red, norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label(f'{target_gene} (log1p norm)', fontsize=9)

            plt.tight_layout()
            sample_str = str(sample).replace(' ', '_').replace('/', '_')
            gene_str = target_gene
            png_path = os.path.join(output_dir, f"spatial_plot_{sample_str}_{gene_str}.png")
            pdf_path = os.path.join(output_dir, f"spatial_plot_{sample_str}_{gene_str}.pdf")
            plt.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"已保存: {png_path}")
            print(f"样本统计: 绘制 {cells_plotted} 个细胞, 表达范围: {vmin:.4f} - {vmax:.4f}")

        except Exception as e:
            print(f"处理样本 {sample} 基因 {target_gene} 时出错: {e}")
            import traceback
            traceback.print_exc()
            plt.close()

    print(f"\n所有图像已保存至: {output_dir}")
    return output_dir

if __name__ == "__main__":
    create_spatial_plot_CXCL9(
        adata_path="TMA2.h5ad",
        boundaries_path="cell_boundaries.csv",
        output_dir="sp_CXCL9"
    )
