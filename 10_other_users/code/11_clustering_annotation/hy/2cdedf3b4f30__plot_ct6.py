import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import os
import warnings

warnings.filterwarnings('ignore')

def plot_celltype_spatial(
    adata_path="subset_data.h5ad",
    boundaries_path="cell_boundaries.csv",
    output_dir="sp_CT6"                 # 修改为 sp_CT6
):
    """
    为每个样本绘制指定细胞类型的空间分布图（按类别颜色填充，无基因表达映射）
    """
    # ==================== 细胞类型 → 填充颜色 ====================
    celltype_colors = {
        "WNT+ ISC": "#FF5733",
        "matri-prod CAF": "#00C9A7",
        "CTSL+ PLAUR+ Macro": "#D81B60",
        "Treg": "#8E24AA",
        "CD8+ T": "#8BC34A",
        "MKI67+ T": "#00B8D4",
    }

    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample_new', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    # 确认所有目标细胞类型存在（至少出现一次）
    existing_types = set(adata.obs['celltype'].unique())
    missing = set(celltype_colors.keys()) - existing_types
    if missing:
        print(f"警告：以下细胞类型在数据中未找到: {missing}")

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

    cell_id_to_idx = pd.Series(range(adata.n_obs), index=adata.obs['cell_id']).to_dict()

    for sample in unique_samples:
        print(f"\n处理样本: {sample}")
        try:
            sample_mask = adata.obs['sample_new'] == sample
            sample_cell_ids = adata.obs.loc[sample_mask, 'cell_id']

            fig, ax = plt.subplots(1, 1, figsize=(15, 12))

            patches, facecolors, edgecolors = [], [], []
            cells_plotted = 0

            for cell_id in sample_cell_ids:
                if cell_id not in cell_boundaries_dict or cell_id not in cell_id_to_idx:
                    continue
                vertices = cell_boundaries_dict[cell_id]
                if len(vertices) < 3:
                    continue

                idx = cell_id_to_idx[cell_id]
                cell_type = adata.obs['celltype'].iloc[idx]

                # 设置填充颜色和边框
                if cell_type in celltype_colors:
                    facecolor = celltype_colors[cell_type]
                else:
                    facecolor = 'white'
                edgecolor = 'lightgray'   # 所有细胞统一浅灰边框

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
            ax.set_title(f'Spatial Cell Type Map - Sample {sample}', fontsize=14, pad=20)

            # 图例（仅显示存在的细胞类型）
            legend_elements = []
            for ct, color in celltype_colors.items():
                if ct in adata.obs['celltype'].values:
                    legend_elements.append(
                        mpatches.Patch(facecolor=color, edgecolor='lightgray', label=ct)
                    )
            legend_elements.append(
                mpatches.Patch(facecolor='white', edgecolor='lightgray', label='Other')
            )
            ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

            plt.tight_layout()
            sample_str = str(sample).replace(' ', '_').replace('/', '_')
            png_path = os.path.join(output_dir, f"celltype_map_{sample_str}.png")
            pdf_path = os.path.join(output_dir, f"celltype_map_{sample_str}.pdf")
            plt.savefig(png_path, dpi=300, bbox_inches='tight')
            plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
            plt.close()

            # 统计输出
            type_counts = {ct: (sample_mask & (adata.obs['celltype'] == ct)).sum()
                           for ct in celltype_colors}
            print(f"已保存: {png_path}")
            print(f"样本细胞总数: {sample_mask.sum()}, 绘制: {cells_plotted}")
            for ct, n in type_counts.items():
                if n > 0:
                    print(f"  {ct}: {n}")

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
        output_dir="sp_CT6"            # 输出文件夹
    )
