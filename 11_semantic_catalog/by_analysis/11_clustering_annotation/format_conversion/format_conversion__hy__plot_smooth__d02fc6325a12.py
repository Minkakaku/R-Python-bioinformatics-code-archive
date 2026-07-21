import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# 创建输出目录
output_dir = "LJD"
os.makedirs(output_dir, exist_ok=True)

# 1. 加载数据
adata = sc.read_h5ad("extracted_samples.h5ad")

# 2. 准备空间坐标
if 'spatial' in adata.obsm:
    coordinates = adata.obsm['spatial']
    x_coords, y_coords = coordinates[:, 0], coordinates[:, 1]
elif 'X' in adata.obs.columns and 'Y' in adata.obs.columns:
    x_coords = adata.obs['X'].values
    y_coords = adata.obs['Y'].values
else:
    raise ValueError("未找到空间坐标信息")

# 3. 获取所有样本和细胞类型
samples = adata.obs['sample'].unique()
all_celltypes = adata.obs['celltype'].unique()
center_celltype = "CXCL9+ CAF"

# 4. 设置新的距离区间：0-100, 100-200, 200-300, 300-400, 400-500
distance_bins = [0, 100, 200, 300, 400, 500]
distance_labels = ['0-100', '100-200', '200-300', '300-400', '400-500']
x_positions = [50, 150, 250, 350, 450]  # 区间中点作为x轴位置

# 5. 存储结果的数据结构
sample_results = {}

# 6. 添加调试信息
print(f"总样本数: {len(samples)}")
print(f"总细胞类型数: {len(all_celltypes)}")
print(f"中心细胞类型: {center_celltype}")
print(f"非中心细胞类型: {[ct for ct in all_celltypes if ct != center_celltype]}")

# 7. 对每个样本进行分析
for sample in samples:
    print(f"\n处理样本: {sample}")
    sample_mask = adata.obs['sample'] == sample
    sample_adata = adata[sample_mask].copy()
    sample_x = x_coords[sample_mask]
    sample_y = y_coords[sample_mask]
    
    # 获取中心细胞的坐标
    center_mask = sample_adata.obs['celltype'] == center_celltype
    if not any(center_mask):
        print(f"样本 {sample} 中没有找到 {center_celltype}，跳过")
        continue
    
    center_coords = np.vstack([sample_x[center_mask], sample_y[center_mask]]).T
    print(f"样本 {sample} 中的中心细胞数量: {len(center_coords)}")
    
    # 为样本中的每个非中心细胞计算距离
    other_cells_mask = sample_adata.obs['celltype'] != center_celltype
    other_coords = np.vstack([sample_x[other_cells_mask], sample_y[other_cells_mask]]).T
    other_celltypes = sample_adata.obs['celltype'].values[other_cells_mask]
    
    if len(other_coords) == 0:
        print(f"样本 {sample} 中没有非中心细胞，跳过")
        continue
    
    print(f"样本 {sample} 中的非中心细胞数量: {len(other_coords)}")
    
    # 获取样本中存在的所有非中心细胞类型
    unique_other_cts = np.unique(other_celltypes)
    print(f"样本 {sample} 中的非中心细胞类型: {unique_other_cts}")
    
    # 构建KDTree用于快速距离计算
    kdtree = KDTree(other_coords)
    
    # 初始化样本结果结构
    for ct in unique_other_cts:
        if ct not in sample_results:
            sample_results[ct] = {}
        if sample not in sample_results[ct]:
            sample_results[ct][sample] = []
    
    # 对每个中心细胞，找出周围非中心细胞并计算距离
    for idx, center_coord in enumerate(center_coords):
        if idx % 100 == 0 and idx > 0:  # 每100个中心细胞打印一次进度
            print(f"  处理第 {idx}/{len(center_coords)} 个中心细胞...")
            
        # 计算距离到最近的所有非中心细胞
        # 注意：查询半径设置为最大距离500，以节省计算时间
        indices = kdtree.query_ball_point(center_coord, r=500)
        
        if not indices:  # 如果没有找到任何细胞
            continue
            
        # 获取距离和细胞类型
        distances = np.linalg.norm(other_coords[indices] - center_coord, axis=1)
        celltypes = other_celltypes[indices]
        
        # 为每个细胞类型统计不同距离区间的细胞数
        ct_counts = {}
        for ct in unique_other_cts:
            ct_mask = celltypes == ct
            if not any(ct_mask):
                # 如果这个中心细胞周围没有该细胞类型，计数为0
                counts = np.zeros(len(distance_bins) - 1, dtype=int)
            else:
                # 获取该细胞类型的距离
                ct_distances = distances[ct_mask]
                
                # 统计不同距离区间的细胞数
                counts, _ = np.histogram(ct_distances, bins=distance_bins)
            
            # 确保中心细胞周围至少有一个该类型的细胞，才记录
            if np.any(counts > 0):
                ct_counts[ct] = counts
        
        # 将结果存储到样本结果中
        for ct, counts in ct_counts.items():
            sample_results[ct][sample].append(counts)

# 8. 计算每个样本每个距离区间的均值
sample_means = {}
sample_counts = {}  # 用于存储每个样本中每个细胞类型的中心细胞数量

for ct in sample_results:
    sample_means[ct] = {}
    sample_counts[ct] = {}
    
    for sample in sample_results[ct]:
        if sample_results[ct][sample]:  # 确保有数据
            # 计算均值
            data_array = np.array(sample_results[ct][sample])
            means = np.mean(data_array, axis=0)
            sample_means[ct][sample] = means
            sample_counts[ct][sample] = len(data_array)  # 记录有多少个中心细胞
            print(f"细胞类型 {ct}, 样本 {sample}: 均值={means}, 基于 {len(data_array)} 个中心细胞")
        else:
            print(f"细胞类型 {ct}, 样本 {sample}: 无数据")

# 9. 按组计算平均值和标准差
groups = adata.obs['group'].unique()
group_stats = {}

# 添加调试信息
print(f"\n=== 分组统计 ===")
print(f"所有组: {groups}")

for ct in sample_means:
    group_stats[ct] = {}
    
    for group in groups:
        # 获取该组的所有样本
        group_samples = adata.obs[adata.obs['group'] == group]['sample'].unique()
        
        # 收集该组所有样本的均值
        group_data = []
        valid_samples = []
        
        for sample in group_samples:
            if sample in sample_means[ct]:
                group_data.append(sample_means[ct][sample])
                valid_samples.append(sample)
        
        if group_data:
            # 计算组内平均值和标准差
            group_mean = np.mean(np.array(group_data), axis=0)
            group_std = np.std(np.array(group_data), axis=0)
            
            # 确保没有负值
            group_mean = np.maximum(group_mean, 0)
            
            group_stats[ct][group] = {
                'mean': group_mean, 
                'std': group_std,
                'n_samples': len(group_data)
            }
            print(f"细胞类型 {ct}, 组 {group}: 平均值={group_mean}, 标准差={group_std}, 基于 {len(group_data)} 个样本")
        else:
            print(f"细胞类型 {ct}, 组 {group}: 无数据")

# 10. 过滤掉没有数据的细胞类型
# 只保留在两个组中都有数据的细胞类型
filtered_celltypes = []
for ct in group_stats:
    if len(group_stats[ct]) >= 2:  # 至少有两个组有数据
        # 检查每个组是否有数据
        has_data = all('mean' in group_stats[ct][group] for group in group_stats[ct])
        if has_data:
            filtered_celltypes.append(ct)
        else:
            print(f"细胞类型 {ct} 不满足数据要求，跳过")
    else:
        print(f"细胞类型 {ct} 只有 {len(group_stats[ct])} 个组有数据，跳过")

print(f"\n有数据的细胞类型: {filtered_celltypes}")
print(f"有数据的细胞类型数量: {len(filtered_celltypes)}")

# 11. 根据参考图片调整颜色和样式
# 参考图片使用红色和蓝色，并且线条样式不同
colors = {}
line_styles = {}
markers = {}

group_list = list(adata.obs['group'].unique())
if len(group_list) >= 2:
    # 假设第一个组是红色实线，第二个组是蓝色虚线
    colors = {group_list[0]: '#E41A1C', group_list[1]: '#377EB8'}
    line_styles = {group_list[0]: '-', group_list[1]: '--'}  # 实线和虚线
    markers = {group_list[0]: 'o', group_list[1]: 's'}  # 圆形和方形标记
    if len(group_list) > 2:
        # 如果超过2个组，使用颜色映射
        cmap = plt.cm.tab10
        line_style_options = ['-', '--', '-.', ':']
        marker_options = ['o', 's', '^', 'v', '<', '>', 'p', '*', 'h', 'H']
        for i, group in enumerate(group_list[2:], start=2):
            colors[group] = cmap(i % 10)
            line_styles[group] = line_style_options[(i-2) % len(line_style_options)]
            markers[group] = marker_options[(i-2) % len(marker_options)]
else:
    # 如果只有一个组，使用红色
    colors = {group_list[0]: '#E41A1C'}
    line_styles = {group_list[0]: '-'}
    markers = {group_list[0]: 'o'}

# 12. 为每个有数据的细胞类型绘制单独图表（参考图片样式）
for ct in filtered_celltypes:
    if len(group_stats[ct]) == 0:
        continue
    
    # 创建图表，根据参考图片调整尺寸
    fig, ax = plt.subplots(figsize=(6, 4.5))
    
    # 绘制每个组的曲线（根据参考图片去掉误差棒）
    for group in group_stats[ct]:
        stats = group_stats[ct][group]
        mean_values = stats['mean']
        
        # 使用预定义颜色、线型和标记
        color = colors.get(group, '#E41A1C')
        line_style = line_styles.get(group, '-')
        marker = markers.get(group, 'o')
        
        # 绘制折线（去掉误差棒）
        ax.plot(x_positions, mean_values, 
                label=f'{group}',
                color=color,
                linestyle=line_style,
                linewidth=2.5,
                marker=marker,
                markersize=6,  # 根据参考图片调整大小
                markerfacecolor='white',
                markeredgewidth=1.5)
    
    # 设置图表属性 - 根据参考图片调整
    ax.set_xlabel('Distance (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean count', fontsize=12, fontweight='bold')
    
    # 根据参考图片，标题在左上角
    ax.set_title(ct, fontsize=12, fontweight='bold', loc='left', pad=10)
    
    # 设置x轴 - 根据参考图片调整
    ax.set_xticks(x_positions)
    ax.set_xticklabels(distance_labels, fontsize=10)
    ax.set_xlim(0, 500)
    
    # 设置y轴从0开始
    ax.set_ylim(bottom=0)
    
    # 自动调整y轴上限，留15%的边距
    y_max = ax.get_ylim()[1]
    if y_max > 0:
        ax.set_ylim(0, y_max * 1.15)
    
    # 添加网格 - 根据参考图片，只有水平网格线
    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    # 根据参考图片，移除上边框和右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 添加图例 - 根据参考图片调整位置和样式
    ax.legend(title='Group', fontsize=10, title_fontsize=11,
             frameon=True, fancybox=True, shadow=False,
             loc='upper left')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表
    safe_ct_name = ct.replace('/', '_').replace('+', 'pos').replace(' ', '_')
    fig.savefig(f'{output_dir}/{safe_ct_name}_distance_analysis.png', 
                dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_dir}/{safe_ct_name}_distance_analysis.pdf', 
                bbox_inches='tight')
    
    plt.close()

# 13. 创建组合图表（模仿参考图片的布局）
if len(filtered_celltypes) > 0:
    # 根据参考图片，我们可以选择前3个细胞类型创建组合图
    selected_celltypes = filtered_celltypes[:3] if len(filtered_celltypes) >= 3 else filtered_celltypes
    
    # 创建3个子图，垂直排列
    fig, axes = plt.subplots(len(selected_celltypes), 1, figsize=(6, 4.5*len(selected_celltypes)))
    
    # 如果只有一个子图，确保axes是数组
    if len(selected_celltypes) == 1:
        axes = [axes]
    
    # 为每个子图添加字母标签（模仿参考图片的E, F, G）
    subplot_labels = ['E', 'F', 'G', 'H', 'I', 'J'][:len(selected_celltypes)]
    
    for idx, (ct, ax) in enumerate(zip(selected_celltypes, axes)):
        if len(group_stats[ct]) == 0:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, 
                    ha='center', va='center', fontsize=12)
            ax.set_title(ct, fontsize=12, fontweight='bold', loc='left', pad=10)
            continue
        
        # 绘制每个组的曲线（去掉误差棒）
        for group in group_stats[ct]:
            stats = group_stats[ct][group]
            mean_values = stats['mean']
            
            # 使用预定义颜色、线型和标记
            color = colors.get(group, '#E41A1C')
            line_style = line_styles.get(group, '-')
            marker = markers.get(group, 'o')
            
            # 绘制折线（去掉误差棒）
            ax.plot(x_positions, mean_values, 
                    label=f'{group}',
                    color=color,
                    linestyle=line_style,
                    linewidth=2.5,
                    marker=marker,
                    markersize=6,
                    markerfacecolor='white',
                    markeredgewidth=1.5)
        
        # 设置子图属性
        ax.set_ylabel('Mean count', fontsize=11, fontweight='bold')
        if idx == len(selected_celltypes) - 1:  # 最后一个子图显示x轴标签
            ax.set_xlabel('Distance (μm)', fontsize=11, fontweight='bold')
        
        # 添加子图标签（字母）
        ax.text(-0.1, 1.05, subplot_labels[idx], transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='right')
        
        # 设置标题在左上角
        ax.set_title(ct, fontsize=12, fontweight='bold', loc='left', pad=10)
        
        # 设置x轴
        ax.set_xticks(x_positions)
        ax.set_xticklabels(distance_labels, fontsize=10)
        ax.set_xlim(0, 500)
        
        # 设置y轴从0开始
        ax.set_ylim(bottom=0)
        
        # 自动调整y轴上限
        y_max = ax.get_ylim()[1]
        if y_max > 0:
            ax.set_ylim(0, y_max * 1.15)
        
        # 添加网格
        ax.grid(True, linestyle='--', alpha=0.3, axis='y')
        
        # 移除上边框和右边框
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # 只在第一个子图添加图例
        if idx == 0:
            ax.legend(title='Group', fontsize=10, title_fontsize=11,
                     frameon=True, fancybox=True, shadow=False,
                     loc='upper left')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存组合图表
    fig.savefig(f'{output_dir}/combined_distance_analysis.png', 
                dpi=300, bbox_inches='tight')
    fig.savefig(f'{output_dir}/combined_distance_analysis.pdf', 
                bbox_inches='tight')
    
    plt.close()

# 14. 保存结果为CSV文件
summary_data = []
for ct in filtered_celltypes:
    for group, stats in group_stats[ct].items():
        for i, (dist_label, mean_val, std_val) in enumerate(zip(distance_labels, stats['mean'], stats['std'])):
            summary_data.append({
                'Celltype': ct,
                'Group': group,
                'Distance_range': dist_label,
                'Mean_count': mean_val,
                'Std_count': std_val,
                'N_samples': stats['n_samples']
            })

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(f'{output_dir}/distance_analysis_summary.csv', index=False)

# 15. 保存每个样本的均值结果
sample_data = []
for ct in filtered_celltypes:
    for sample in sample_means.get(ct, {}):
        group = adata.obs[adata.obs['sample'] == sample]['group'].iloc[0]
        n_centers = sample_counts.get(ct, {}).get(sample, 0)
        for i, (dist_label, mean_val) in enumerate(zip(distance_labels, sample_means[ct][sample])):
            sample_data.append({
                'Celltype': ct,
                'Sample': sample,
                'Group': group,
                'Distance_range': dist_label,
                'Mean_count': mean_val,
                'N_centers': n_centers
            })

sample_df = pd.DataFrame(sample_data)
sample_df.to_csv(f'{output_dir}/sample_mean_counts.csv', index=False)

# 16. 保存处理日志
from datetime import datetime
with open(f'{output_dir}/processing_log.txt', 'w') as f:
    f.write("=== 空间距离分析处理日志 ===\n")
    f.write(f"处理时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"总样本数: {len(samples)}\n")
    f.write(f"总细胞类型数: {len(all_celltypes)}\n")
    f.write(f"中心细胞类型: {center_celltype}\n")
    f.write(f"距离区间: {distance_labels}\n")
    f.write(f"有数据的细胞类型数: {len(filtered_celltypes)}\n")
    f.write(f"有数据的细胞类型: {filtered_celltypes}\n")
    f.write(f"组列表: {list(groups)}\n")

print(f"\n分析完成！结果已保存到 {output_dir} 目录")
print(f"分析了 {len(filtered_celltypes)} 个有数据的细胞类型")
print(f"已生成 {len(filtered_celltypes)} 张单个细胞类型图表")
if len(filtered_celltypes) > 0:
    print(f"已生成 1 张组合图表")
print(f"数据汇总已保存为 CSV 文件")
print(f"处理日志已保存")
