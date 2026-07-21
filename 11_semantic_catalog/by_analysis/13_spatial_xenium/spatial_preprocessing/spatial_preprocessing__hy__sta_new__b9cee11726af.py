import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# 1. 读取h5ad数据
print("读取h5ad数据...")
# 请将"your_data.h5ad"替换为您的实际文件名
adata = sc.read_h5ad("adata_celltype.h5ad")
print(f"数据形状: {adata.shape}")
print(f"细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
print(f"obs列: {adata.obs.columns.tolist()}")

# 检查必需的列
if 'sample' not in adata.obs.columns:
    print("错误: 数据中缺少 'sample' 列")
    exit(1)
    
if 'celltype' not in adata.obs.columns:
    print("错误: 数据中缺少 'celltype' 列")
    exit(1)

# 2. 添加group_final列
print("\n添加group_final列...")

# 定义映射字典
sample_to_group_final = {
    "PRE_5-1": "PRE-R", "PRE_5-2": "PRE-R", "PRE_6-1": "PRE-R", 
    "PRE_6-2": "PRE-R", "PRE_6-3": "PRE-R", "PRE_8": "PRE-R", 
    "PRE_15-1": "PRE-R", "PRE_15-2": "PRE-R", "PRE_15-3": "PRE-R", 
    "PRE_15-4": "PRE-R", "PRE_16": "PRE-R", "PRE_17-1": "PRE-R", 
    "PRE_17-2": "PRE-R", "PRE_17-3": "PRE-R", "PRE_18-1": "PRE-R", 
    "PRE_18-2": "PRE-R", "PRE_19": "PRE-R", "PRE_3-1": "PRE-NR", 
    "PRE_3-2": "PRE-NR", "PRE_10-1": "PRE-NR", "PRE_10-2": "PRE-NR", 
    "PRE_20-1": "PRE-NR", "PRE_20-2": "PRE-NR", "PRE_20-3": "PRE-NR", 
    "PRE_20-4": "PRE-NR", "PRE_21": "PRE-NR", "PRE_22-1": "PRE-NR", 
    "PRE_22-2": "PRE-NR", "PRE_23-1": "PRE-NR", "PRE_23-2": "PRE-NR", 
    "PRE_25-1": "PRE-NR", "PRE_25-2": "PRE-NR", "PRE_25-3": "PRE-NR", 
    "PRE_25-4": "PRE-NR", "PRE_26-1": "PRE-NR", "PRE_26-3": "PRE-NR", 
    "PRE_27-1": "PRE-NR", "PRE_27-2": "PRE-NR", "POST_5": "POST-R", 
    "POST_6": "POST-R", "POST_8": "POST-R", "POST_15": "POST-R", 
    "POST_16": "POST-R", "POST_18": "POST-R", "POST_19": "POST-R", 
    "POST_3": "POST-NR", "POST_20": "POST-NR", "POST_21": "POST-NR", 
    "POST_22-1": "POST-NR", "POST_25": "POST-NR", "POST_26": "POST-NR", 
    "POST_27": "POST-NR"
}

# 应用映射
adata.obs['group_final'] = adata.obs['sample'].map(sample_to_group_final)

# 检查未映射的样本
unmapped_samples = adata.obs[adata.obs['group_final'].isna()]['sample'].unique()
if len(unmapped_samples) > 0:
    print(f"警告: 以下样本未找到group_final映射关系: {unmapped_samples.tolist()}")
    print(f"未映射的细胞数: {adata.obs['group_final'].isna().sum()}")
    # 移除未映射的样本
    adata = adata[~adata.obs['group_final'].isna(), :].copy()
    print(f"移除未映射样本后细胞数: {adata.shape[0]}")

# 统计group_final分布
print(f"\ngroup_final分布:")
group_counts = adata.obs['group_final'].value_counts()
for group, count in group_counts.items():
    print(f"  {group}: {count} 个细胞 ({count/len(adata)*100:.2f}%)")

# 3. 以样本为维度计算细胞类型比例
print("\n计算每个样本的细胞类型比例...")

# 提取需要的列
df = adata.obs[['sample', 'celltype', 'group_final']].copy()

# 计算每个样本中每个细胞类型的细胞数
sample_cell_counts = df.groupby(['sample', 'celltype']).size().unstack(fill_value=0)

# 计算每个样本的总细胞数
sample_totals = sample_cell_counts.sum(axis=1)

# 计算每个样本中每个细胞类型的比例（百分比）
sample_proportions = (sample_cell_counts.div(sample_totals, axis=0) * 100).round(4)

# 添加group_final信息
sample_groups = df[['sample', 'group_final']].drop_duplicates().set_index('sample')
sample_proportions['group_final'] = sample_proportions.index.map(sample_groups['group_final'])

# 重新组织数据为长格式，便于绘图
sample_long = sample_proportions.reset_index().melt(
    id_vars=['sample', 'group_final'], 
    value_vars=[col for col in sample_proportions.columns if col not in ['sample', 'group_final']],
    var_name='celltype', 
    value_name='proportion'
)

# 移除比例为0的行
sample_long = sample_long[sample_long['proportion'] > 0]

# 显示统计信息
print(f"样本数: {sample_proportions.shape[0]}")
print(f"细胞类型数: {len(sample_cell_counts.columns)}")
print(f"group_final分组: {sorted(sample_proportions['group_final'].unique().tolist())}")

# 保存样本维度的数据
sample_proportions.to_csv("sample_celltype_proportions.csv", encoding='utf-8-sig')
print(f"样本细胞类型比例数据已保存到: sample_celltype_proportions.csv")

# 4. 创建输出目录
output_dir = "celltype_proportion_analysis"
os.makedirs(output_dir, exist_ok=True)

# 5. 绘制多系列箱线图（一张图中包含所有细胞类型）
print(f"\n绘制多系列箱线图...")
celltypes = df['celltype'].unique()
print(f"总细胞类型数: {len(celltypes)}")

# 设置图形大小（根据细胞类型数量调整宽度）
fig_width = max(20, len(celltypes) * 0.8)  # 每个细胞类型0.8个单位宽度
fig_height = 8
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# 设置颜色调色板，为四个分组定义颜色
group_colors = {
    'PRE-R': '#1f77b4',  # 蓝色
    'PRE-NR': '#ff7f0e',  # 橙色
    'POST-R': '#2ca02c',  # 绿色
    'POST-NR': '#d62728'  # 红色
}

# 定义每个细胞类型的位置
x_positions = np.arange(len(celltypes))

# 计算每个箱线图的宽度
box_width = 0.15  # 每个箱线图的宽度

# 定义每个分组相对于中心位置的偏移
group_offsets = {
    'PRE-R': -1.5 * box_width,
    'PRE-NR': -0.5 * box_width,
    'POST-R': 0.5 * box_width,
    'POST-NR': 1.5 * box_width
}

# 绘制每个细胞类型、每个分组的箱线图
for i, celltype in enumerate(celltypes):
    for group_name, offset in group_offsets.items():
        # 获取当前细胞类型和分组的数据
        celltype_group_data = sample_long[
            (sample_long['celltype'] == celltype) & 
            (sample_long['group_final'] == group_name)
        ]['proportion'].values
        
        if len(celltype_group_data) > 0:
            # 计算箱线图的位置
            box_position = i + offset
            
            # 绘制箱线图
            boxprops = dict(facecolor=group_colors[group_name], color='black')
            whiskerprops = dict(color='black')
            capprops = dict(color='black')
            medianprops = dict(color='black')
            
            bp = ax.boxplot(
                celltype_group_data, 
                positions=[box_position],
                widths=box_width,
                patch_artist=True,
                boxprops=boxprops,
                whiskerprops=whiskerprops,
                capprops=capprops,
                medianprops=medianprops
            )
            
            # 添加散点显示每个样本
            scatter_x = np.random.normal(box_position, 0.02, size=len(celltype_group_data))
            ax.scatter(scatter_x, celltype_group_data, color='black', alpha=0.6, s=20, zorder=3)

# 设置x轴
ax.set_xticks(x_positions)
ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=9)

# 设置y轴
ax.set_ylabel('Cell Proportion (%)', fontsize=12)
ax.set_ylim(bottom=0)  # 确保y轴从0开始

# 添加网格
ax.grid(True, axis='y', alpha=0.3, linestyle='--')

# 添加图例
legend_handles = []
for group_name, color in group_colors.items():
    legend_handles.append(plt.Rectangle((0, 0), 1, 1, fc=color, label=group_name))
ax.legend(handles=legend_handles, title='Groups', fontsize=10, title_fontsize=11)

# 添加标题
ax.set_title('Cell Type Proportions by Group', fontsize=14, fontweight='bold', pad=20)

# 调整布局
plt.tight_layout()

# 保存图形
output_path = os.path.join(output_dir, "all_celltypes_group_comparison.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"图形已保存到: {output_path}")

# 6. 创建分组的汇总统计
print(f"\n创建分组汇总统计...")

# 计算每个细胞类型在每个分组中的平均比例
summary_stats = sample_long.groupby(['celltype', 'group_final'])['proportion'].agg([
    'mean', 'median', 'std', 'count', 'min', 'max'
]).reset_index()

# 重新排列列顺序
summary_stats.columns = ['celltype', 'group_final', 'mean_proportion', 'median_proportion', 
                         'std_proportion', 'sample_count', 'min_proportion', 'max_proportion']

# 保存汇总统计
summary_stats.to_csv(os.path.join(output_dir, "celltype_group_summary_stats.csv"), 
                     index=False, encoding='utf-8-sig')
print(f"汇总统计已保存到: {output_dir}/celltype_group_summary_stats.csv")

# 7. 生成汇总报告
print(f"\n汇总报告:")
print(f"  总细胞类型数: {len(celltypes)}")
print(f"  样本数: {sample_proportions.shape[0]}")
print(f"  分组: {sorted(sample_proportions['group_final'].unique().tolist())}")
print(f"  图形保存到: {output_path}")
print(f"  汇总统计保存到: {output_dir}/celltype_group_summary_stats.csv")

# 8. 保存更新后的h5ad数据
output_h5ad = "data_with_group_final.h5ad"
adata.write_h5ad(output_h5ad, compression='gzip')
print(f"\n更新后的h5ad数据已保存到: {output_h5ad}")

print("\n分析完成!")
