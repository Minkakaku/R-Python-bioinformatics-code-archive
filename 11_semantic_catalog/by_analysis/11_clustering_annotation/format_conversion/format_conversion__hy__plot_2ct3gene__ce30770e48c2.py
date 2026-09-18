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

def create_spatial_plots(
    adata_path="subset_data.h5ad",
    boundaries_path="cell_boundaries.csv",
    output_dir="2CT_3gene"
):
    """
    为每个样本、每个基因创建空间原位图：
    - 目标细胞：matri-prod CAF（固定绿色）和 WNT+ ISC（阈值判定，高深红/低浅红）
    - 目标基因：MYC, SOX9, ASCL2
    - X矩阵为 raw counts，内部进行文库大小归一化+log1p
    """
    # ========================= 数据读取 =========================
    print("正在读取h5ad数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据维度: {adata.shape}")

    for col in ['sample_new', 'celltype', 'cell_id']:
        if col not in adata.obs.columns:
            raise ValueError(f"obs中缺少'{col}'列")

    target_genes = ['MYC', 'SOX9', 'ASCL2']
    for gene in target_genes:
        if gene not in adata.var_names:
            raise ValueError(f"var_names中未找到'{gene}'")

    # ========================= 标准化 =========================
    print("正在进行文库大小归一化和log1p转换...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # 提取各基因表达（稠密数组）
    gene_expr_dict = {}
    for gene in target_genes:
        expr = adata_norm[:, gene].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        gene_expr_dict[gene] = expr

    # ========================= 读取边界 =========================
    print("正在读取细胞边界文件...")
    boundaries_df = pd.read_csv(boundaries_path)
    cell_boundaries_dict = {}
    for cell_id, group in boundaries_df.groupby('cell_id'):
        vertices = list(zip(group['vertex_x'], group['vertex_y']))
        cell_boundaries_dict[cell_id] = vertices

    # ========================= 准备绘图 =========================
    # 颜色定义
    color_caf = '#33A02C'          # 绿色（matri-prod CAF 固定）
    color_other = '#E0F7FA'        # 浅蓝色（其他细胞）

    # 定义 WNT+ ISC 的阈值（根据您提供的折线图提取）
    gene_thresholds = {
        'MYC': 0.7,
        'SOX9': 0.5,
        'ASCL2': 0.35
    }

    # 高低表达的颜色
    color_high = '#E60000'         # 深红（高于等于阈值）
    color_low = '#FFCCCC'          # 浅红（低于阈值）

    unique_samples = adata_norm.obs['sample_new'].unique()
    print(f"找到 {len(unique_samples)} 个样本: {list(unique_samples)}")
    os.makedirs(output_dir, exist_ok=True)

    cell_id_to_idx = pd.Series(range(adata_norm.n_obs), index=adata_norm.obs['cell_id']).to_dict()

    for sample in unique_samples:
        print(f"\n处理样本: {sample}")
        try:
            sample_mask = adata_norm.obs['sample_new'] == sample
            sample_cell_ids = adata_norm.obs.loc[sample_mask, 'cell_id']

            # 获取该样本中 WNT+ ISC 的细胞索引
            isc_mask = sample_mask & (adata_norm.obs['celltype'] == 'WNT+ ISC')
            isc_indices = np.where(isc_mask)[0] if isc_mask.any() else []

            # 对每个基因生成一张图
            for gene in target_genes:
                print(f"  绘制基因: {gene}")

                fig, ax = plt.subplots(1, 1, figsize=(15, 12))

                expr_vec = gene_expr_dict[gene]
                threshold = gene_thresholds[gene]

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
                        facecolor = color_caf
                        edgecolor = '#1F7A2A'
                    elif cell_type == 'WNT+ ISC':
                        expr = expr_vec[idx]
                        if np.isnan(expr):
                            facecolor = '#FFFFFF'
                        else:
                            # 使用阈值二分化着色
                            if expr >= threshold:
                                facecolor = color_high
                            else:
                                facecolor = color_low
                        edgecolor = '#B22222'
                    else:
                        facecolor = color_other
                        edgecolor = '#B0C4DE'

                    polygon = Polygon(vertices, closed=True)
                    patches.append(polygon)
                    facecolors.append(facecolor)
                    edgecolors.append(edgecolor)
                    cells_plotted += 1

                if not patches:
                    print(f"样本 {sample} 基因 {gene} 没有有效边界数据，跳过")
                    plt.close()
                    continue

                # 将 linewidths 从 0.3 调整为 0.05，边框更细
                p = PatchCollection(patches, facecolors=facecolors, edgecolors=edgecolors,
                                    linewidths=0.05, alpha=0.9)
                ax.add_collection(p)
                ax.autoscale_view()
                ax.set_aspect('equal')
                ax.set_xlabel('X Position', fontsize=12)
                ax.set_ylabel('Y Position', fontsize=12)
                ax.set_title(f'Spatial Plot - Sample {sample}, Gene {gene}\n'
                             'matri-prod CAF (fixed green) | WNT+ ISC (High/Low)',
                             fontsize=14, pad=20)

                # 动态更新图例，加入高/低表达区分
                legend_elements = [
                    mpatches.Patch(facecolor=color_caf, edgecolor='#1F7A2A',
                                   label='matri-prod CAF'),
                    mpatches.Patch(facecolor=color_high, edgecolor='#B22222',
                                   label=f'WNT+ ISC ({gene} High)'),
                    mpatches.Patch(facecolor=color_low, edgecolor='#B22222',
                                   label=f'WNT+ ISC ({gene} Low)'),
                    mpatches.Patch(facecolor=color_other, edgecolor='#B0C4DE',
                                   label='Other Cells')
                ]
                ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

                # 二元着色，不再使用连续的 colorbar
                plt.tight_layout()
                sample_str = str(sample).replace(' ', '_').replace('/', '_')
                gene_str = gene.replace(' ', '_')
                png_path = os.path.join(output_dir, f"spatial_{sample_str}_{gene_str}.png")
                pdf_path = os.path.join(output_dir, f"spatial_{sample_str}_{gene_str}.pdf")
                plt.savefig(png_path, dpi=300, bbox_inches='tight')
                plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
                plt.close()

                print(f"  已保存: {png_path}")
                print(f"  统计: 绘制 {cells_plotted} 个细胞")
                print(f"    matri-prod CAF: {(sample_mask & (adata_norm.obs['celltype'] == 'matri-prod CAF')).sum()}, "
                      f"WNT+ ISC: {len(isc_indices)}, 总细胞: {sample_mask.sum()}")
                if len(isc_indices) > 0:
                    print(f"    阈值: {threshold} | 深红: 表达 ≥ {threshold}, 浅红: 表达 < {threshold}")

        except Exception as e:
            print(f"处理样本 {sample} 时出错: {e}")
            import traceback
            traceback.print_exc()
            plt.close()

    print(f"\n所有图像已保存至: {output_dir}")
    return output_dir

if __name__ == "__main__":
    create_spatial_plots(
        adata_path="subset_data.h5ad",
        boundaries_path="cell_boundaries.csv",
        output_dir="2CT_3gene"
    )
