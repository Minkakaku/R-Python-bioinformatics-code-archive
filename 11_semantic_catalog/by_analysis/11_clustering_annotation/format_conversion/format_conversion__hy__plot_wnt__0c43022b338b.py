import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import os
import warnings

warnings.filterwarnings('ignore')

def plot_celltype_spatial(
    adata_path="subset_data.h5ad",
    boundaries_path="cell_boundaries.csv",
    output_dir="sp_CT6",
    target_celltype="WNT+ ISC",   # 目标细胞类型
    target_gene="CD274"           # 目标基因
):
    """
    所有细胞都绘制：非目标细胞为白填充+浅灰边框；
    目标细胞(WNT+ ISC)用蓝色边框，并按CD274表达量红色渐变填充。
    """
    # ==================== 颜色设置 ====================
    # 非目标细胞统一样式
    OTHER_FACE = 'white'
    OTHER_EDGE = 'lightgray'

    # 目标细胞边框颜色
    TARGET_EDGE = 'blue'

    # 目标细胞填充颜色映射：使用红色渐变
    cmap = plt.cm.Reds
    norm = None  # 根据目标细胞表达范围动态计算

    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample_new', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    # 检查目标基因
    if target_gene not in adata.var_names:
        raise ValueError(f"基因 '{target_gene}' 不在数据集的变量中")

    # 提取目标基因表达值
    gene_idx = adata.var_names.get_loc(target_gene)
    if hasattr(adata.X, "toarray"):
        expr = adata.X[:, gene_idx].toarray().flatten()
    else:
        expr = adata.X[:, gene_idx].flatten()

    # 筛选目标细胞类型
    is_target = (adata.obs['celltype'] == target_celltype).values

    # 准备表达值向量：非目标细胞设为 NaN，方便跳过填充颜色
    target_expr = np.where(is_target, expr, np.nan)

    # 基于有效值设定归一化范围
    valid_expr = target_expr[~np.isnan(target_expr)]
    if len(valid_expr) == 0:
        raise ValueError(f"目标细胞类型 '{target_celltype}' 中没有有效的表达数据")
    norm = Normalize(vmin=np.min(valid_expr), vmax=np.max(valid_expr))

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    unique_samples = adata.obs['sample_new'].unique()
    print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")
    os.makedirs(output_dir, exist_ok=True)

    # 构建 cell_id -> obs 索引映射
    cell_id_to_idx = pd.Series(range(adata.n_obs), index=adata.obs['cell_id']).to_dict()

    for sample in unique_samples:
        print(f"\n处理样本: {sample}")
        try:
            sample_mask = adata.obs['sample_new'] == sample
            sample_cell_ids = adata.obs.loc[sample_mask, 'cell_id']

            fig, ax = plt.subplots(1, 1, figsize=(15, 12))
            patches, facecolors, edgecolors = [], [], []
            cells_plotted = 0
            target_cells_plotted = 0

            for cell_id in sample_cell_ids:
                if cell_id not in cell_boundaries_dict or cell_id not in cell_id_to_idx:
                    continue
                vertices = cell_boundaries_dict[cell_id]
                if len(vertices) < 3:
                    continue

                idx = cell_id_to_idx[cell_id]
                polygon = Polygon(vertices, closed=True)
                patches.append(polygon)

                if is_target[idx] and not np.isnan(target_expr[idx]):
                    # 目标细胞：蓝色边框，红色填充（按表达量）
                    color = cmap(norm(target_expr[idx]))
                    facecolors.append(color)
                    edgecolors.append(TARGET_EDGE)
                    target_cells_plotted += 1
                else:
                    # 非目标细胞（或无表达的目标细胞，极少数情况）：白填充+灰边
                    facecolors.append(OTHER_FACE)
                    edgecolors.append(OTHER_EDGE)

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
            ax.set_title(f'{target_gene} expression in {target_celltype} - Sample {sample}',
                         fontsize=14, pad=20)

            # 颜色条
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label(f'{target_gene} expression', fontsize=11)

            plt.tight_layout()
            sample_str = str(sample).replace(' ', '_').replace('/', '_')
            png_path = os.path.join(output_dir, f"celltype_map_{sample_str}.png")
            pdf_path = os.path.join(output_dir, f"celltype_map_{sample_str}.pdf")
            plt.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"已保存: {png_path}")
            print(f"样本细胞总数绘制: {cells_plotted}, 其中目标细胞({target_celltype})绘制: {target_cells_plotted}")

        except Exception as e:
            print(f"处理样本 {sample} 时出错: {e}")
            import traceback
            traceback.print_exc()
            plt.close()

    print(f"\n所有图像已保存至: {output_dir}")
    return output_dir

if __name__ == "__main__":
    plot_celltype_spatial(
        adata_path="subset_data.h5ad",
        boundaries_path="cell_boundaries.csv",
        output_dir="sp_wnt_274"
    )
