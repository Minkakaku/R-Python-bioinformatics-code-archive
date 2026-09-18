import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

# 1. 加载h5ad数据
adata = sc.read_h5ad("../adata_celltype.h5ad")  # 替换为实际文件路径

# 2. 筛选Myeloids细胞
myeloids = adata[adata.obs['celltype'] == 'Myeloids'].copy()

# 3. 定义目标基因列表
#target_genes = ["CSF1R", "CEBPA", "LILRB2", "TYK2", "IL17RA", "IKBKB", "PYCARD", "CD4", "IRAK1",
#                "CSNK2B", "PTK2B", "PTPN6", "PLP2", "CCL19", "P4HB", "PRKACA", "JAK3", "CCR2"]

#target_genes = ["UPF1","CSF1R","CEBPA","LILRB2","TYK2","RHOA","IKBKB","PYCARD","CD4","SIRPA","PTPN6","ARHGEF2","PLP2","CCL19","PRKACA","GAS6","JAK3","CCR2"]
#target_genes = ["APP","CEBPA","GRN","IL10RA","PRKCD","TNC","AOAH","PARK7","IL17RA","PYCARD","ALOX5","SIRPA","NAIP","CCR2","BRD4"]
target_genes = ["PYCARD","APP","SIRPA","CYBA","ARRB1","PTPN6","ARHGEF2","LILRB2","GAS6","IL17RA"]


# 4. 检查并过滤不存在的基因
if myeloids.raw is not None:
    valid_genes = [g for g in target_genes if g in myeloids.raw.var_names]
    # 从raw属性提取原始计数数据
    raw_myeloids = myeloids.raw.to_adata()
    counts_matrix = raw_myeloids[:, valid_genes].X.copy()
else:
    print("Warning: No raw data found. Using .X instead.")
    valid_genes = [g for g in target_genes if g in myeloids.var_names]
    counts_matrix = myeloids[:, valid_genes].X.copy()

print(f"找到 {len(valid_genes)}/{len(target_genes)} 个有效基因")

# 5. 转换为DataFrame
if hasattr(counts_matrix, 'todense'):
    expr_df = pd.DataFrame(
        counts_matrix.todense(),
        columns=valid_genes,
        index=myeloids.obs_names
    )
else:
    expr_df = pd.DataFrame(
        counts_matrix,
        columns=valid_genes,
        index=myeloids.obs_names
    )

# 6. 添加样本和分组信息
expr_df['sample'] = myeloids.obs['sample'].values
expr_df['group'] = myeloids.obs['group'].values

# 7. 对数归一化
for gene in valid_genes:
    expr_df[gene] = np.log1p(expr_df[gene])

# 8. 按样本计算平均表达量
sample_means = expr_df.groupby('sample', observed=True)[valid_genes].mean().reset_index()

# 9. 合并分组信息
group_info = expr_df[['sample', 'group']].drop_duplicates()
sample_means = pd.merge(sample_means, group_info, on='sample')

# 10. Z-score标准化基因表达量 (按基因)
scaled_means = sample_means.copy()
for gene in valid_genes:
    scaled_means[gene] = (scaled_means[gene] - scaled_means[gene].mean()) / scaled_means[gene].std()

# 11. 排序样本（按分组）
scaled_means = scaled_means.sort_values('group')

# 12. 创建绘图数据
plot_data = scaled_means[valid_genes].T  # 转置为(基因x样本)
sample_labels = scaled_means['sample'].values
group_labels = scaled_means['group'].values

# 13. 创建自定义图例
unique_groups = pd.unique(group_labels)
# 使用第一张图的配色：绿色(disease)和黄色(normal)
group_color_map = {
    'disease': '#F0E68C',  # 绿色
    'normal': '#66CDAA'    # 黄色
}
# 如果分组名称不同，自动分配颜色
if not all(g in group_color_map for g in unique_groups):
    palette = sns.color_palette("Set3", len(unique_groups))
    group_color_map = {g: palette[i] for i, g in enumerate(unique_groups)}

# 14. 创建热图
fig = plt.figure(figsize=(16, 10))

# 创建网格布局：顶部分组条、主体热图和右侧颜色条
# 优化布局，预留更多空间给标题
gs = gridspec.GridSpec(3, 2, 
                      width_ratios=[0.95, 0.05],  # 减小右侧颜色条宽度
                      height_ratios=[0.06, 0.9, 0.04],  # 增加顶部空间
                      wspace=0.05, 
                      hspace=0.05)

# 绘制分组标签 - 更窄的指示条
ax0 = plt.subplot(gs[0, 0])
current_group = None
group_start = 0
for i, group in enumerate(group_labels):
    if group != current_group:
        if current_group is not None:
            # 绘制上一个组的矩形
            ax0.add_patch(mpatches.Rectangle(
                (group_start, 0), i - group_start, 1, 
                color=group_color_map[current_group], ec='k'
            ))
        current_group = group
        group_start = i
# 绘制最后一个组的矩形
ax0.add_patch(mpatches.Rectangle(
    (group_start, 0), len(group_labels) - group_start, 1, 
    color=group_color_map[current_group], ec='k'
))

ax0.set_xlim(0, len(group_labels))
ax0.set_ylim(0, 1)
ax0.set_xticks([])
ax0.set_yticks([])
# 移除ax0的标题，使用统一的标题
ax0.set_title("", fontsize=12)

# 绘制主热图
ax1 = plt.subplot(gs[1, 0])
# 创建蓝色到红色的颜色映射
blue_red_colors = ["#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", 
                   "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"]
cmap = LinearSegmentedColormap.from_list("blue_red", blue_red_colors)

# 绘制热图 - 增加色块间距
sns.heatmap(
    plot_data,
    cmap=cmap,
    ax=ax1,
    cbar=False,  # 稍后单独添加颜色条
    yticklabels=True,
    linewidths=0.8,  # 增加色块间距
    linecolor='white',  # 白色间隔线
    vmin=-0.4,
    vmax=1.0
)
ax1.set_ylabel("Genes", fontsize=12)
ax1.set_xlabel("Samples", fontsize=12)
# 移除热图子图的标题，使用统一的标题

# 设置样本标签 - 避免重叠问题
ax1.set_xticks(np.arange(len(sample_labels)) + 0.5)
ax1.set_xticklabels(
    sample_labels, 
    rotation=45, 
    ha="right", 
    rotation_mode="anchor",
    fontsize=9
)

# 设置基因标签
ax1.set_yticks(np.arange(len(valid_genes)) + 0.5)
ax1.set_yticklabels(valid_genes, fontsize=9, rotation=0)

# 添加右侧颜色条（更窄）
ax2 = plt.subplot(gs[1, 1])
cbar = plt.colorbar(ax1.collections[0], cax=ax2, aspect=15)  # aspect控制宽度
cbar.set_label('Z-scored Expression', fontsize=10)
cbar.ax.tick_params(labelsize=8)
cbar.set_ticks([-0.4, 0.0, 1.0])  # 设置与图片匹配的刻度

# 15. 添加分组图例
legend_handles = [mpatches.Patch(color=c, label=g) for g, c in group_color_map.items()]
ax0.legend(
    handles=legend_handles,
    loc='center left',
    bbox_to_anchor=(1.2, 0.5),
    frameon=False,
    title='Groups',
    fontsize=9
)

# 16. 添加统一的标题（放置在顶部中央）
fig.suptitle("Regulation Of Interleukin-6 Production Differentially Expressed Genes", fontsize=14, fontweight='normal', y=0.97)

# 17. 调整布局并保存
plt.tight_layout()
plt.subplots_adjust(top=0.92)  # 为顶部标题留出更多空间
plt.savefig("Regulation_Of_Interleukin-6_gene_expression_heatmap.pdf", dpi=300, bbox_inches='tight')
plt.show()
