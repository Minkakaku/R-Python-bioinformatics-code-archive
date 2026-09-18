import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import matplotlib as mpl
import warnings

# 忽略特定警告
warnings.filterwarnings("ignore", message="Variable names are not unique.")
warnings.filterwarnings("ignore", message="is_categorical_dtype")
warnings.filterwarnings("ignore", message="use_inf_as_na")

# 设置绘图全局参数
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# 1. 读取h5ad文件
adata = sc.read_h5ad("integrated_data.h5ad")  # 替换为你的h5ad文件路径

# 2. 检查空间坐标是否存在
if 'spatial' not in adata.obsm_keys():
    raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")
    
print(f"Found spatial coordinates: {adata.obsm['spatial'].shape}")

# 3. 提取POSTN+ FB细胞
postn_fb_mask = adata.obs['celltype_final'] == 'POSTN+ FB'
if postn_fb_mask.sum() == 0:
    raise ValueError("No POSTN+ FB cells found in the dataset")
    
adata_postn = adata[postn_fb_mask].copy()

# 4. 创建其他细胞的数据集（排除POSTN+ FB自身）
other_mask = adata.obs['celltype_final'] != 'POSTN+ FB'
adata_other = adata[other_mask].copy()

# 使用空间坐标进行距离计算
X_postn = adata_postn.obsm['spatial']
X_other = adata_other.obsm['spatial']

# 打印空间坐标范围以供参考
print(f"POSTN+ FB spatial coordinates range: x=[{X_postn[:,0].min():.2f}, {X_postn[:,0].max():.2f}], y=[{X_postn[:,1].min():.2f}, {X_postn[:,1].max():.2f}]")
print(f"Other cells spatial coordinates range: x=[{X_other[:,0].min():.2f}, {X_other[:,0].max():.2f}], y=[{X_other[:,1].min():.2f}, {X_other[:,1].max():.2f}]")

# 5. 构建KDTree进行最近邻搜索
print("Building KDTree for spatial coordinates...")
kdtree = cKDTree(X_other)
distances, indices = kdtree.query(X_postn, k=1)  # 每个POSTN+ FB细胞的最近邻

print(f"Distance stats: min={distances.min():.2f}, max={distances.max():.2f}, mean={distances.mean():.2f}")

# 6. 获取最近邻细胞信息
# 6.1 从adata_other获取最近邻的cellID和细胞类型
nn_cells = adata_other.obs.index[indices]
nn_celltypes = adata_other.obs['celltype_final'].values[indices]

# 6.2 创建结果DataFrame
results_df = pd.DataFrame({
    'query_cell': adata_postn.obs.index,
    'query_type': 'POSTN+ FB',
    'query_group': adata_postn.obs['group'],
    'nn_cell': nn_cells,
    'nn_celltype': nn_celltypes,
    'distance': distances
})

# 7. 输出最近邻细胞列表文件 (文件1)
nn_list_df = results_df[['nn_cell', 'nn_celltype']]
nn_list_df.to_csv("nearest_neighbor_cells_spatial.csv", index=False)

# 8. 统计最近邻细胞类型结果 (文件2)
# 统计每个query_group中各种nn_celltype的数量
counts_df = results_df.groupby(['query_group', 'nn_celltype']).size().unstack(fill_value=0)
# 计算频率
freq_df = counts_df.div(counts_df.sum(axis=1), axis=0)

# 输出完整的统计结果
stats_df = freq_df.reset_index().melt(id_vars='query_group', 
                                     var_name='celltype', 
                                     value_name='proportion')
stats_df = stats_df.sort_values(['query_group', 'proportion'], ascending=[True, False])

# 保存统计结果
with pd.ExcelWriter('nearest_neighbor_stats_spatial.xlsx') as writer:
    results_df.to_excel(writer, sheet_name='All Cells', index=False)
    counts_df.to_excel(writer, sheet_name='Cell Counts')
    freq_df.to_excel(writer, sheet_name='Cell Frequencies')
    stats_df.to_excel(writer, sheet_name='Proportions', index=False)

# 9. 绘制分组柱状图与趋势线
# 准备绘图数据 - 只取前15个最丰富的细胞类型
# 计算平均频率并排序
average_freq = freq_df.mean(axis=0)
top_types = average_freq.sort_values(ascending=False).index[:15]
freq_df_top = freq_df[top_types]

celltypes = freq_df_top.columns
disease_freq = freq_df_top.loc['disease'].values
normal_freq = freq_df_top.loc['normal'].values

# 设置柱状图位置
bar_width = 0.35
r1 = np.arange(len(celltypes))
r2 = [x + bar_width + 0.1 for x in r1]  # 增加0.1作为两组之间的间隔

# 创建图形
fig, ax = plt.subplots(figsize=(16, 8))

# 创建更丰富的颜色映射
cmap = plt.get_cmap('tab20')
colors_disease = [cmap(i % 20) for i in range(len(celltypes))]
colors_normal = [cmap((i + 10) % 20) for i in range(len(celltypes))]  # 偏移10种颜色

# 绘制分组柱状图
bars1 = ax.bar(r1, disease_freq, width=bar_width, 
               color=colors_disease, edgecolor='black', 
               alpha=0.9, label='Disease')

bars2 = ax.bar(r2, normal_freq, width=bar_width, 
               color=colors_normal, edgecolor='black', 
               alpha=0.9, label='Normal')

# 添加连接线 - 用比柱状图淡20%的颜色
for i in range(len(celltypes)):
    # 获取柱子的坐标信息
    disease_bar = bars1[i]
    normal_bar = bars2[i]
    
    # 提取Disease柱子右侧的两个端点
    disease_xr = disease_bar.get_x() + disease_bar.get_width()
    disease_top = disease_bar.get_y() + disease_bar.get_height()
    disease_bottom = disease_bar.get_y()
    
    # 提取Normal柱子左侧的两个端点
    normal_xl = normal_bar.get_x()
    normal_top = normal_bar.get_y() + normal_bar.get_height()
    normal_bottom = normal_bar.get_y()
    
    # 创建更淡的颜色
    light_color_disease = (*mpl.colors.to_rgb(colors_disease[i]), 0.4)  # RGBA格式
    light_color_normal = (*mpl.colors.to_rgb(colors_normal[i]), 0.4)  # RGBA格式
    
    # 添加连接线
    ax.plot([disease_xr, normal_xl], [disease_top, normal_top], 
            c=light_color_disease, lw=1.5, linestyle='-', alpha=0.7)
    ax.plot([disease_xr, normal_xl], [disease_bottom, normal_bottom], 
            c=light_color_normal, lw=1.5, linestyle='-', alpha=0.7)

    # 添加填充区域增强可视化效果
    ax.fill_between([disease_xr, normal_xl], 
                   [disease_bottom, normal_bottom],
                   [disease_top, normal_top],
                   color=light_color_disease, alpha=0.1)

# 添加标签和标题
ax.set_xlabel('Cell Types', fontsize=14)
ax.set_ylabel('Proportion', fontsize=14)
ax.set_title('Spatial Nearest Neighbors of POSTN+ FB Cells by Group', fontsize=16)
ax.set_xticks([r + bar_width/2 for r in r1])
ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=10)
ax.legend(fontsize=12)

# 添加网格以增强可读性
ax.grid(True, linestyle='--', alpha=0.3)

# 空间信息标注
ax.text(0.02, 0.95, 'Spatial Analysis', transform=ax.transAxes, 
        fontsize=14, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# 添加距离统计信息
dist_stats = f"Distance Stats (min/avg/max): {distances.min():.2f} / {distances.mean():.2f} / {distances.max():.2f}"
ax.text(0.02, 0.02, dist_stats, transform=ax.transAxes, 
        fontsize=11, verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.tight_layout()
plt.savefig('nearest_neighbor_analysis_spatial.pdf', dpi=300, bbox_inches='tight')

# 10. 附加：绘制空间分布图
if 'spatial' in adata.uns:
    print("Generating spatial plots...")
    # 创建分组图
    fig, axs = plt.subplots(1, 2, figsize=(16, 8), sharex=True, sharey=True)
    groups = ['disease', 'normal']
    titles = ['Disease Group', 'Normal Group']
    
    for i, group in enumerate(groups):
        # 选择当前组的POSTN+ FB细胞
        group_mask = (adata.obs['group'] == group) & postn_fb_mask
        adata_group = adata[group_mask].copy()
        
        # 绘制空间分布
        sc.pl.spatial(adata_group, color='celltype_final', 
                      title=titles[i], ax=axs[i],
                      show=False, size=0.5 if adata_group.n_obs < 5000 else 0.1)
        
        # 添加比例信息
        proportion = len(adata_group) / len(adata_postn)
        axs[i].text(0.05, 0.95, f"{len(adata_group)} cells ({proportion:.1%})", 
                   transform=axs[i].transAxes, fontsize=12,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))
    
    plt.tight_layout()
    plt.savefig('postn_fb_spatial_distribution.pdf', dpi=300, bbox_inches='tight')
    
    print("Analysis completed. Spatial plot saved.")
else:
    print("Spatial plot skipped: No spatial image information found in adata.uns['spatial']")

print("All outputs saved successfully.")
