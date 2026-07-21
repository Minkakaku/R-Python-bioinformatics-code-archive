import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LinearSegmentedColormap
import os
import warnings

warnings.filterwarnings('ignore')

def plot_celltype_spatial(
    adata_path="TMA1.h5ad",
    boundaries_path="cell_boundaries.csv",
    output_dir="sp_epi_cd274",
    target_gene="CD274"
):
    """
    所有细胞轮廓都绘制：
    - 目标细胞（15种）：蓝色边框
        * 若 CD274 > 0：浅红 → 深红渐变填充（无白色）
        * 若 CD274 = 0：纯白填充
    - 其余细胞：白底灰边
    """
    # ==================== 目标细胞类型 ====================
    target_celltypes = ['Proliferative Epithelials','Hypoxic Epithelials','IFN-Epithelials','TA-like Epithelials','Wnt-active stem、','stem-like Tumor']

    # ==================== 颜色定义 ====================
    OTHER_FACE = 'white'
    OTHER_EDGE = 'lightgray'
    TARGET_EDGE = 'blue'

    # 仅用于 >0 表达的渐变：从明显的浅红到深红（不含白）
    cmap = LinearSegmentedColormap.from_list('LightToDarkRed', ['lightcoral', 'darkred'])

    # ========================= 读数据 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    if target_gene not in adata.var_names:
        raise ValueError(f"基因 '{target_gene}' 不在数据集中")

    gene_idx = adata.var_names.get_loc(target_gene)
    if hasattr(adata.X, "toarray"):
        expr = adata.X[:, gene_idx].toarray().flatten()
    else:
        expr = adata.X[:, gene_idx].flatten()

    is_target = adata.obs['celltype'].isin(target_celltypes).values
    target_expr = np.where(is_target, expr, np.nan)

    # 非零表达的范围用于颜色映射
    non_zero_expr = target_expr[~np.isnan(target_expr) & (target_expr > 0)]
    if len(non_zero_expr) == 0:
        norm = Normalize(vmin=0, vmax=1)  # 所有表达都是0时，不会用到cmap
    else:
        norm = Normalize(vmin=np.min(non_zero_expr), vmax=np.max(non_zero_expr))

    # ========================= 读边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 绘图 =========================
    unique_samples = adata.obs['sample'].unique()
    print(f"样本数: {len(unique_samples)}")
    os.makedirs(output_dir, exist_ok=True)

    cell_id_to_idx = pd.Series(range(adata.n_obs), index=adata.obs['cell_id']).to_dict()

    for sample in unique_samples:
        print(f"\n处理样本: {sample}")
        try:
            sample_mask = adata.obs['sample'] == sample
            sample_cell_ids = adata.obs.loc[sample_mask, 'cell_id']

            fig, ax = plt.subplots(1, 1, figsize=(15, 12))
            patches, facecolors, edgecolors = [], [], []
            cells_plotted = 0
            target_plotted = 0

            for cell_id in sample_cell_ids:
                if cell_id not in cell_boundaries_dict or cell_id not in cell_id_to_idx:
                    continue
                vertices = cell_boundaries_dict[cell_id]
                if len(vertices) < 3:
                    continue

                idx = cell_id_to_idx[cell_id]
                polygon = Polygon(vertices, closed=True)
                patches.append(polygon)

                if is_target[idx]:
                    edgecolors.append(TARGET_EDGE)
                    if not np.isnan(target_expr[idx]) and target_expr[idx] > 0:
                        # 非零表达：浅红 → 深红
                        facecolors.append(cmap(norm(target_expr[idx])))
                    else:
                        # 零表达：纯白
                        facecolors.append('white')
                    target_plotted += 1
                else:
                    facecolors.append(OTHER_FACE)
                    edgecolors.append(OTHER_EDGE)
                cells_plotted += 1

            if not patches:
                plt.close()
                continue

            p = PatchCollection(patches, facecolors=facecolors, edgecolors=edgecolors,
                                linewidths=0.1, alpha=0.9)
            ax.add_collection(p)
            ax.autoscale_view()
            ax.set_aspect('equal')
            ax.set_xlabel('X', fontsize=12)
            ax.set_ylabel('Y', fontsize=12)
            ax.set_title(f'CD274 in target cells - Sample {sample}', fontsize=14, pad=20)

            # 颜色条仅显示非零渐变范围
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label('CD274 (>0)', fontsize=11)

            plt.tight_layout()
            sample_str = str(sample).replace(' ', '_').replace('/', '_')
            plt.savefig(os.path.join(output_dir, f"celltype_map_{sample_str}.png"), dpi=300)
            plt.savefig(os.path.join(output_dir, f"celltype_map_{sample_str}.pdf"), dpi=300)
            plt.close()

            print(f"绘制细胞: {cells_plotted}, 目标: {target_plotted}")

        except Exception as e:
            print(f"样本 {sample} 出错: {e}")
            plt.close()

    print(f"\n图像保存在: {output_dir}")
    return output_dir

if __name__ == "__main__":
    plot_celltype_spatial(
        adata_path="TMA1.h5ad",
        boundaries_path="cell_boundaries.csv",
        output_dir="sp_epi_cd274"
    )
