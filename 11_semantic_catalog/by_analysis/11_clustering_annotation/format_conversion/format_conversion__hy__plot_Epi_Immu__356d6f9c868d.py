import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import KDTree
import matplotlib as mpl
import math

# 设置全局字体
plt.rcParams['font.size'] = 18

# 1. 加载数据并提取指定细胞类型
def load_and_filter_data():
    adata = sc.read_h5ad('All_filter_state_samp_new_ctzq_with_CNs10-15.h5ad')
    
    # 中心细胞类型：来自 celltype_zq 列
    center_celltype = "Epithelial"
    # 周边（目标）细胞类型列表（来自 celltype 列）
    target_celltypes = [
        "CTSL+ PLAUR+ Macro",
        "STAB1+ Macro",
        "CD4+ T",
        "γδ T",
        "CD8+ T",
        "Treg",
        "MKI67+ T",
        "DC",
        "GCB",
        "MBC",
        "Mast",
        "Plasmas",
        "Neutrophils"
    ]
    
    # 检查列是否存在
    if 'celltype_zq' not in adata.obs.columns:
        raise KeyError("adata.obs 中缺少 'celltype_zq' 列")
    if 'celltype' not in adata.obs.columns:
        raise KeyError("adata.obs 中缺少 'celltype' 列")
    
    # 获取中心细胞索引
    center_mask = adata.obs['celltype_zq'] == center_celltype
    if center_mask.sum() == 0:
        raise ValueError("没有找到任何 Epithelial 细胞（celltype_zq == 'Epithelial'）")
    
    # 获取周边细胞：celltype 在目标列表中，并且排除中心细胞本身（避免重复）
    target_mask = adata.obs['celltype'].isin(target_celltypes) & (~center_mask)
    
    # 合并需要保留的细胞：中心细胞 + 周边细胞
    final_mask = center_mask | target_mask
    adata = adata[final_mask].copy()
    
    # 验证实际存在的细胞类型
    existing_celltypes = adata.obs['celltype'].unique()
    missing_types = [ct for ct in target_celltypes if ct not in existing_celltypes]
    if missing_types:
        print(f"警告: 以下目标细胞类型在数据中不存在: {', '.join(missing_types)}")
    
    print(f"中心细胞 (Epithelial): {center_mask.sum()} 个")
    print(f"周边目标细胞: {target_mask.sum()} 个")
    
    return adata, center_celltype, target_celltypes

# 2. 计算不同距离范围内的细胞数量和密度
def calculate_cell_counts(adata, center_celltype, target_celltypes):
    # 获取中心细胞坐标
    center_mask = adata.obs['celltype_zq'] == center_celltype
    center_coords = adata[center_mask].obsm['spatial']
    num_centers = len(center_coords)
    
    if num_centers == 0:
        raise ValueError("没有找到中心细胞坐标")
    
    # 获取周边细胞（非中心）
    other_mask = ~center_mask
    other_coords = adata[other_mask].obsm['spatial']
    other_celltypes = adata[other_mask].obs['celltype']
    
    # 处理分组信息：如果存在 'state' 列则使用，否则所有细胞归为一组
    if 'state' in adata.obs.columns:
        other_state = adata[other_mask].obs['state']
        # 收集所有出现过的 state 值（仅保留 N, A, C 等常见分组，全部保留也可以）
        all_groups = sorted(set(other_state.dropna().unique()))
        # 如果全是数字或其他，直接转为字符串列表
        all_groups = [str(g) for g in all_groups]
        group_to_idx = {g: i for i, g in enumerate(all_groups)}
        other_group = other_state.astype(str).map(group_to_idx).values
    else:
        # 无 state 列，创建单一分组
        all_groups = ['All']
        group_to_idx = {'All': 0}
        other_group = np.zeros(len(other_coords), dtype=int)
    
    # 构建 KDTree
    tree = KDTree(other_coords)
    
    # 距离阈值 (μm)
    distances = [0, 50, 100, 150, 200]
    
    # 非中心细胞类型只包含我们分析的目标列表（排除可能存在的其他类型，但我们已经筛选过，这里再确保）
    non_center_celltypes = [ct for ct in target_celltypes if ct in other_celltypes.unique()]
    num_celltypes = len(non_center_celltypes)
    num_groups = len(all_groups)
    
    # 初始化计数数组：[距离环数, 细胞类型数, 分组数]
    counts = np.zeros((len(distances)-1, num_celltypes, num_groups))
    
    # 对每个中心细胞计算周围细胞分布
    for center in center_coords:
        indices = tree.query_ball_point(center, r=distances[-1])  # 最远距离
        for idx in indices:
            dist = np.linalg.norm(center - other_coords[idx])
            ct = other_celltypes.iloc[idx]
            grp_idx = other_group[idx]
            
            if ct not in non_center_celltypes:
                continue
            ct_idx = non_center_celltypes.index(ct)
            
            # 确定距离区间
            for i in range(len(distances)-1):
                if distances[i] <= dist < distances[i+1]:
                    counts[i, ct_idx, grp_idx] += 1
                    break
    
    # 计算平均每个中心细胞的计数
    mean_counts = counts / num_centers
    
    # 计算扇形面积和细胞密度
    sector_areas = np.zeros((len(distances)-1, num_celltypes))
    densities = np.zeros_like(mean_counts)
    
    for dist_idx in range(len(distances)-1):
        r1 = distances[dist_idx]
        r2 = distances[dist_idx+1]
        ring_area = math.pi * (r2**2 - r1**2)
        sector_area = ring_area / num_celltypes  # 每个细胞类型等分圆环
        for sector_idx in range(num_celltypes):
            sector_areas[dist_idx, sector_idx] = sector_area
            for grp_idx in range(num_groups):
                densities[dist_idx, sector_idx, grp_idx] = mean_counts[dist_idx, sector_idx, grp_idx] / sector_area
    
    return mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups

def export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups, filename='Epithelial_immune_cell_counts.csv'):
    dist_labels = [f"{distances[i]}-{distances[i+1]}" for i in range(len(distances)-1)]
    
    data = []
    for ct_idx, cell_type in enumerate(non_center_celltypes):
        for dist_idx in range(len(dist_labels)):
            area = sector_areas[dist_idx, ct_idx]
            for grp_idx, grp_name in enumerate(all_groups):
                data.append([
                    cell_type,
                    dist_labels[dist_idx],
                    grp_name,
                    mean_counts[dist_idx, ct_idx, grp_idx],
                    area,
                    densities[dist_idx, ct_idx, grp_idx]
                ])
    
    df = pd.DataFrame(data, columns=[
        'celltype', 'distance_range', 'state', 'cell_count', 'sector_area', 'cell_density'
    ])
    df.to_csv(filename, sep='\t', index=False)
    print(f"细胞计数已导出到 {filename}")

# 创建径向气泡图（标签外移、图例右移、画布加宽）
def create_radial_plot(mean_counts, densities, distances, non_center_celltypes, all_celltypes, all_groups):
    num_groups = len(all_groups)
    
    # 增大画布宽度
    plt.figure(figsize=(16, 12))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    
    # 删除默认网格和标签
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    
    # 绘制同心圆（距离边界）
    for r in distances[1:]:
        ax.plot(np.linspace(0, 2*np.pi, 100), [r]*100, color='gray', linestyle='-', linewidth=1, alpha=0.5)
    
    # 绘制分隔线（数量等于所有细胞类型数量）
    num_sectors = len(all_celltypes)
    for i in range(num_sectors):
        angle = i * 2*np.pi / num_sectors
        ax.plot([angle, angle], [0, 230], color='#34495E', linestyle='-', linewidth=1, alpha=0.6)
    
    # 添加细胞类型标签（外圈——向外移动）
    for i, ct in enumerate(all_celltypes):
        angle = (i + 0.5) * 2*np.pi / num_sectors
        y = 255                        # 向外移动
        ax.text(angle, y, ct, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.7))
    
    # 每个扇形内，各组气泡的偏移角度
    total_offset_span = 0.28
    if num_groups > 1:
        offsets = np.linspace(-total_offset_span/2, total_offset_span/2, num_groups)
    else:
        offsets = [0]
    
    # 颜色映射（如果有多个组）
    colormap = plt.cm.tab10
    colors = [colormap(i % 10) for i in range(num_groups)]
    
    # 气泡大小缩放因子
    scale_factors = [5000, 8000, 10000, 12000]
    
    # 绘制气泡
    for ring_idx in range(len(distances)-1):
        r_center = distances[ring_idx+1] - 5
        for ct_idx, cell_type in enumerate(non_center_celltypes):
            all_ct_idx = all_celltypes.index(cell_type)
            base_angle = (all_ct_idx + 0.5) * 2*np.pi / num_sectors
            density_vals = densities[ring_idx, ct_idx, :]
            for grp_idx in range(num_groups):
                angle = base_angle + offsets[grp_idx]
                size = density_vals[grp_idx] * scale_factors[min(ring_idx, len(scale_factors)-1)]
                size = np.clip(size, 5, 500)
                color = colors[grp_idx]
                ax.scatter(angle, r_center, s=size, facecolor=color, edgecolor='none', alpha=0.8, zorder=10)
    
    # 扩大 y 轴上限，确保标签不被裁剪
    ax.set_ylim(0, 270)
    
    # 添加距离范围标签（在圆环外侧）
    label_angles = np.linspace(0, 2*np.pi, len(distances))[:-1]
    for i, dist in enumerate(distances[1:]):
        label = f"{distances[i]}-{distances[i+1]}μm"
        ax.text(label_angles[i], distances[i+1]+15, label, ha='center', va='center', fontsize=8, fontweight='bold')
    
    # 图例右移
    legend_title = "State" if num_groups > 1 and 'All' not in all_groups else "Group"
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=grp_name)
                       for grp_name, color in zip(all_groups, colors)]
    
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.6, 0.5),
              title=legend_title, framealpha=0.9)
    
    # 压缩主图区域
    plt.subplots_adjust(right=0.55)
    
    ax.text(0, 265, "Immune Cell Density around Epithelial", color='black', ha='center', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('Epithelial_immune_plot_group.png', dpi=300, bbox_inches='tight')
    plt.savefig('Epithelial_immune_plot_group.pdf', dpi=300, bbox_inches='tight')
    plt.show()

# 主程序
if __name__ == "__main__":
    # 1. 加载并过滤数据
    adata, center_celltype, target_celltypes = load_and_filter_data()
    
    # 2. 计算细胞数量统计和密度
    mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups = calculate_cell_counts(
        adata, center_celltype, target_celltypes
    )
    
    # 导出计数结果到CSV
    export_counts_to_csv(mean_counts, densities, distances, sector_areas, non_center_celltypes, all_groups,
                         filename='Epithelial_immune_cell_counts.csv')
    
    # 3. 创建环形图（周边细胞类型列表需保持与 all_celltypes 一致）
    create_radial_plot(mean_counts, densities, distances, non_center_celltypes, target_celltypes, all_groups)
