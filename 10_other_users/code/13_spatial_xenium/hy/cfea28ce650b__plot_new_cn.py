import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
import warnings
warnings.filterwarnings('ignore')

# 创建输出目录
output_dir = "plot_new_cn"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 创建单独箱线图输出目录
boxplot_dir = f"{output_dir}/individual_boxplots"
if not os.path.exists(boxplot_dir):
    os.makedirs(boxplot_dir)

# 1. 加载数据
print("Loading data...")
adata = sc.read_h5ad("extracted_samples.h5ad")
print(f"Data shape: {adata.shape}")
print(f"Original annotation column: {adata.obs['annotation'].unique()}")
print(f"Using CNs_10 column: {adata.obs['CNs_10'].unique()}")
print(f"Patient count: {adata.obs['patient'].nunique()}")
print(f"Sample count: {adata.obs['sample'].nunique()}")

# 检查CNs_10列是否存在
if 'CNs_10' not in adata.obs.columns:
    print("Error: CNs_10 column not found in adata.obs!")
    print(f"Available columns: {adata.obs.columns.tolist()}")
    exit()

# 2. 预处理和检查数据完整性
print("\nChecking data integrity...")
# 确保每个patient在inner和outer都有数据
patient_zone_counts = adata.obs.groupby(['patient', 'mucus_zone']).size().unstack(fill_value=0)
print("Cell counts per patient in inner/outer regions:")
print(patient_zone_counts)

# 检查是否有病人缺少inner或outer数据
patients_with_both = patient_zone_counts.dropna().index.tolist()
print(f"\nPatients with both inner and outer data: {len(patients_with_both)}")
print(f"Patient IDs: {patients_with_both}")

# 3. 计算每个patient在每个区域的细胞类型比例
print("\nCalculating cell type proportions (using CNs_10)...")

# 获取所有细胞类型 (基于CNs_10列)
all_cell_types = adata.obs['CNs_10'].unique()
print(f"Total cell types (CNs_10): {len(all_cell_types)}")
print(f"Cell types: {all_cell_types}")

# 为每个病人创建比例矩阵
patient_proportions = {}
patients_to_analyze = []  # 用于存储有足够数据的病人

for patient in patients_with_both:
    patient_data = adata[adata.obs['patient'] == patient]
    
    # 提取inner区域
    inner_data = patient_data[patient_data.obs['mucus_zone'] == 'inner']
    inner_counts = inner_data.obs['CNs_10'].value_counts()
    inner_total = len(inner_data)
    
    # 提取outer区域
    outer_data = patient_data[patient_data.obs['mucus_zone'] == 'outer']
    outer_counts = outer_data.obs['CNs_10'].value_counts()
    outer_total = len(outer_data)
    
    # 检查两个区域是否都有细胞
    if inner_total == 0 or outer_total == 0:
        print(f"Warning: Patient {patient} has no cells in inner or outer region, skipping")
        continue
    
    # 计算比例
    inner_props = {cell_type: inner_counts.get(cell_type, 0) / inner_total for cell_type in all_cell_types}
    outer_props = {cell_type: outer_counts.get(cell_type, 0) / outer_total for cell_type in all_cell_types}
    
    patient_proportions[patient] = {
        'inner': inner_props,
        'outer': outer_props
    }
    patients_to_analyze.append(patient)

print(f"\nPatients analyzed: {len(patients_to_analyze)}")
if len(patients_to_analyze) < 2:
    print("Error: At least 2 patients with complete data are required for paired analysis")
    exit()

# 4. 将数据转换为DataFrame以便分析
prop_data = []

for patient in patients_to_analyze:
    inner_props = patient_proportions[patient]['inner']
    outer_props = patient_proportions[patient]['outer']
    
    for cell_type in all_cell_types:
        prop_data.append({
            'patient': patient,
            'cell_type': cell_type,
            'zone': 'inner',
            'proportion': inner_props[cell_type]
        })
        
        prop_data.append({
            'patient': patient,
            'cell_type': cell_type,
            'zone': 'outer',
            'proportion': outer_props[cell_type]
        })

prop_df = pd.DataFrame(prop_data)

# 5. 检查每个细胞类型在所有病人中的存在情况
print("\nChecking cell type presence across patients...")
cell_type_presence = {}

for cell_type in all_cell_types:
    # 检查在多少病人中这个细胞类型在inner或outer区域存在
    present_in_patients = 0
    for patient in patients_to_analyze:
        inner_val = patient_proportions[patient]['inner'][cell_type]
        outer_val = patient_proportions[patient]['outer'][cell_type]
        if inner_val > 0 or outer_val > 0:
            present_in_patients += 1
    
    cell_type_presence[cell_type] = present_in_patients

# 只保留在至少3个病人中存在的细胞类型进行统计分析
min_patients = 3
valid_cell_types = [ct for ct, count in cell_type_presence.items() if count >= min_patients]
print(f"Cell types present in at least {min_patients} patients: {len(valid_cell_types)}")
print(f"Filtered cell types: {set(all_cell_types) - set(valid_cell_types)}")

# 6. 执行配对样本统计检验
print("\nPerforming paired statistical tests...")

test_results = []

for cell_type in valid_cell_types:
    # 提取该细胞类型的数据
    cell_data = prop_df[prop_df['cell_type'] == cell_type]
    
    # 转换为宽格式
    wide_data = cell_data.pivot(index='patient', columns='zone', values='proportion')
    
    # 确保有足够的数据点
    if len(wide_data) < 3:  # 至少需要3对数据
        print(f"Warning: Insufficient data for {cell_type}, skipping statistical test")
        continue
    
    # 提取inner和outer的比例
    inner_props = wide_data['inner'].values
    outer_props = wide_data['outer'].values
    
    # 计算描述性统计
    mean_inner = np.mean(inner_props)
    mean_outer = np.mean(outer_props)
    std_inner = np.std(inner_props, ddof=1)
    std_outer = np.std(outer_props, ddof=1)
    
    # 配对t检验
    t_stat, t_p = stats.ttest_rel(inner_props, outer_props)
    
    # Wilcoxon符号秩检验（非参数）
    w_stat, w_p = stats.wilcoxon(inner_props, outer_props)
    
    # 计算均值差异
    mean_diff = mean_inner - mean_outer
    rel_change = mean_diff / mean_outer if mean_outer > 0 else np.nan
    
    test_results.append({
        'cell_type': cell_type,
        'n_patients': len(wide_data),
        'mean_inner': mean_inner,
        'mean_outer': mean_outer,
        'std_inner': std_inner,
        'std_outer': std_outer,
        'mean_diff': mean_diff,
        'rel_change': rel_change,
        't_statistic': t_stat,
        't_pvalue': t_p,
        'wilcoxon_stat': w_stat,
        'wilcoxon_pvalue': w_p
    })

if not test_results:
    print("Error: No valid cell types for statistical analysis")
    exit()

# 转换为DataFrame
results_df = pd.DataFrame(test_results)

# 多重检验校正
results_df['t_pvalue_adj'] = multipletests(results_df['t_pvalue'], method='fdr_bh')[1]
results_df['wilcoxon_pvalue_adj'] = multipletests(results_df['wilcoxon_pvalue'], method='fdr_bh')[1]

# 排序结果
results_df = results_df.sort_values('wilcoxon_pvalue_adj')

# 7. 保存统计结果
print("\nSaving statistical results...")
results_df.to_csv(f"{output_dir}/paired_test_results.csv", index=False)
print(f"Statistical results saved to: {output_dir}/paired_test_results.csv")

# 8. 可视化
print("\nGenerating visualizations...")

# 8.1 总体细胞类型组成条形图
plt.figure(figsize=(14, 8))

# 计算总体平均比例
overall_props = prop_df.groupby(['cell_type', 'zone'])['proportion'].mean().unstack()

# 只包含有效细胞类型
overall_props = overall_props.loc[valid_cell_types]

# 排序
overall_props['diff'] = overall_props['inner'] - overall_props['outer']
overall_props = overall_props.sort_values('diff', ascending=False)

# 绘制条形图
x = np.arange(len(overall_props))
width = 0.35

fig, ax = plt.subplots(figsize=(max(12, len(overall_props)/2), 8))
bars1 = ax.bar(x - width/2, overall_props['inner'], width, label='inner', alpha=0.8, color='skyblue')
bars2 = ax.bar(x + width/2, overall_props['outer'], width, label='outer', alpha=0.8, color='lightcoral')

# 添加显著性标记
for i, cell_type in enumerate(overall_props.index):
    if cell_type in results_df['cell_type'].values:
        p_val = results_df[results_df['cell_type'] == cell_type]['wilcoxon_pvalue_adj'].values[0]
        
        if p_val < 0.001:
            sig = '***'
        elif p_val < 0.01:
            sig = '**'
        elif p_val < 0.05:
            sig = '*'
        else:
            sig = 'ns'
    else:
        sig = 'NA'
    
    y_pos = max(overall_props.loc[cell_type, 'inner'], overall_props.loc[cell_type, 'outer']) + 0.01
    
    ax.text(i, y_pos, sig, ha='center', fontsize=10, fontweight='bold')

ax.set_xlabel('Cell Type (CNs_10)', fontsize=12)
ax.set_ylabel('Average Proportion', fontsize=12)
ax.set_title('Cell Type Composition Comparison: Inner vs Outer Regions (CNs_10)', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(overall_props.index, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_dir}/overall_composition.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_dir}/overall_composition.pdf", bbox_inches='tight')
plt.close()

# 8.2 显著差异细胞类型的配对点线图
# 筛选显著差异的细胞类型
sig_cell_types = results_df[results_df['wilcoxon_pvalue_adj'] < 0.05]['cell_type'].tolist()

if sig_cell_types:
    print(f"\nFound {len(sig_cell_types)} significantly different cell types:")
    for ct in sig_cell_types:
        row = results_df[results_df['cell_type'] == ct].iloc[0]
        p_val = row['wilcoxon_pvalue_adj']
        diff = row['mean_diff']
        dir_text = "inner > outer" if diff > 0 else "inner < outer"
        print(f"  - {ct}: p={p_val:.4f}, direction: {dir_text}")
    
    # 创建子图
    n_sig = len(sig_cell_types)
    n_cols = 3
    n_rows = (n_sig + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_sig == 1:
        axes = np.array([axes])
    if n_rows * n_cols > 1 and n_sig > 1:
        axes = axes.flatten()
    
    for i, cell_type in enumerate(sig_cell_types[:n_rows*n_cols]):
        if n_sig > 1:
            ax = axes[i]
        else:
            ax = axes
        
        # 提取该细胞类型的数据
        cell_data = prop_df[prop_df['cell_type'] == cell_type]
        
        # 绘制配对点线图
        for patient in cell_data['patient'].unique():
            patient_data = cell_data[cell_data['patient'] == patient]
            inner_val = patient_data[patient_data['zone'] == 'inner']['proportion'].values[0]
            outer_val = patient_data[patient_data['zone'] == 'outer']['proportion'].values[0]
            
            # 绘制连线
            ax.plot([0, 1], [inner_val, outer_val], 'gray', alpha=0.5, linewidth=1)
            # 绘制点
            ax.scatter(0, inner_val, s=50, alpha=0.7, color='skyblue')
            ax.scatter(1, outer_val, s=50, alpha=0.7, color='lightcoral')
        
        # 添加均值线
        mean_inner = cell_data[cell_data['zone'] == 'inner']['proportion'].mean()
        mean_outer = cell_data[cell_data['zone'] == 'outer']['proportion'].mean()
        ax.axhline(mean_inner, color='blue', linestyle='--', alpha=0.7, linewidth=2, xmin=0, xmax=0.4)
        ax.axhline(mean_outer, color='red', linestyle='--', alpha=0.7, linewidth=2, xmin=0.6, xmax=1)
        
        # 获取p值
        p_val = results_df[results_df['cell_type'] == cell_type]['wilcoxon_pvalue_adj'].values[0]
        
        # 设置图形属性
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['inner', 'outer'])
        ax.set_ylabel('Proportion')
        ax.set_title(f'{cell_type}\np={p_val:.4f}')
        ax.grid(alpha=0.3)
    
    # 隐藏多余的子图
    if n_sig > 1:
        for i in range(len(sig_cell_types), len(axes)):
            axes[i].axis('off')
    
    plt.suptitle('Paired Comparison of Significantly Different Cell Types (CNs_10)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/significant_paired_comparisons.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/significant_paired_comparisons.pdf", bbox_inches='tight')
    plt.close()
else:
    print("No significantly different cell types found")

# 8.3 差异热图
# 创建差异矩阵
diff_matrix = pd.pivot_table(
    prop_df[prop_df['cell_type'].isin(valid_cell_types)], 
    values='proportion', 
    index='patient', 
    columns=['cell_type', 'zone']
)

# 计算每个病人每个细胞类型的差异
patient_diffs = {}
for patient in patients_to_analyze:
    diffs = []
    for cell_type in valid_cell_types:
        try:
            inner_val = diff_matrix.loc[patient, (cell_type, 'inner')]
            outer_val = diff_matrix.loc[patient, (cell_type, 'outer')]
            diffs.append(inner_val - outer_val)
        except:
            diffs.append(0)
    patient_diffs[patient] = diffs

diff_df = pd.DataFrame(patient_diffs, index=valid_cell_types).T

# 绘制热图
plt.figure(figsize=(max(12, len(valid_cell_types)), max(8, len(patients_to_analyze)/2)))

# 添加显著性标记
sig_mask = pd.DataFrame(False, index=valid_cell_types, columns=patients_to_analyze)
for cell_type in valid_cell_types:
    if cell_type in results_df['cell_type'].values:
        p_val = results_df[results_df['cell_type'] == cell_type]['wilcoxon_pvalue_adj'].values[0]
        if p_val < 0.05:
            sig_mask.loc[cell_type, :] = True

# 转置以便热图显示
diff_df_T = diff_df.T

# 创建热图
fig, ax = plt.subplots(figsize=(max(12, len(diff_df.columns)/3), max(8, len(diff_df.index)/2)))
sns.heatmap(diff_df_T, cmap='RdBu_r', center=0, 
            cbar_kws={'label': 'Inner - Outer Proportion Difference'},
            ax=ax)

# 添加显著性标记
for i, cell_type in enumerate(diff_df_T.index):
    for j, patient in enumerate(diff_df_T.columns):
        if sig_mask.loc[cell_type, patient]:
            # 在热图单元格中添加星号
            ax.text(j + 0.5, i + 0.5, '*', 
                   ha='center', va='center', 
                   fontsize=10, fontweight='bold',
                   color='white')

ax.set_xlabel('Patient')
ax.set_ylabel('Cell Type (CNs_10)')
ax.set_title('Cell Type Composition Differences by Patient (Inner - Outer, CNs_10)', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{output_dir}/difference_heatmap.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_dir}/difference_heatmap.pdf", bbox_inches='tight')
plt.close()

# 8.4 每个病人的细胞组成堆叠条形图
# 为每个病人创建子图
n_patients = len(patients_to_analyze)
n_cols = min(4, n_patients)
n_rows = (n_patients + n_cols - 1) // n_cols

if n_patients > 0:
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
    if n_patients == 1:
        axes = np.array([axes])
    if n_rows * n_cols > 1:
        axes = axes.flatten()
    
    for i, patient in enumerate(patients_to_analyze[:n_rows*n_cols]):
        if n_patients > 1:
            ax = axes[i]
        else:
            ax = axes
        
        # 提取该病人的数据
        patient_data = prop_df[prop_df['patient'] == patient]
        
        # 转换为宽格式用于堆叠图
        stack_data = patient_data.pivot(index='zone', columns='cell_type', values='proportion')
        
        # 绘制堆叠条形图
        stack_data.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
        
        ax.set_title(f'Patient {patient}')
        ax.set_xlabel('Mucus Zone')
        ax.set_ylabel('Cell Proportion')
        if i == 0:  # 只在第一个图中显示图例
            ax.legend(title='Cell Type (CNs_10)', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
        else:
            ax.legend().remove()
        ax.grid(axis='y', alpha=0.3)
        ax.set_ylim([0, 1])
    
    # 隐藏多余的子图
    if n_patients > 1:
        for i in range(len(patients_to_analyze), len(axes)):
            axes[i].axis('off')
    
    plt.suptitle('Cell Composition by Patient in Inner/Outer Regions (CNs_10)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/patient_composition_stacked.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/patient_composition_stacked.pdf", bbox_inches='tight')
    plt.close()

# 8.5 综合箱线图（所有细胞类型在一起）
print("\nGenerating combined boxplot...")

# 筛选数据用于箱线图
boxplot_data = prop_df[prop_df['cell_type'].isin(valid_cell_types)].copy()

# 为了更好的可视化，可以限制显示的细胞类型数量（按平均比例排序）
cell_type_means = boxplot_data.groupby('cell_type')['proportion'].mean().sort_values(ascending=False)
top_cell_types = cell_type_means.head(30).index.tolist()  # 只显示前30个最常见的细胞类型
boxplot_data_filtered = boxplot_data[boxplot_data['cell_type'].isin(top_cell_types)]

plt.figure(figsize=(max(16, len(top_cell_types) * 0.5), 8))
ax = sns.boxplot(data=boxplot_data_filtered, x='cell_type', y='proportion', hue='zone', 
                 palette={'inner': 'skyblue', 'outer': 'lightcoral'}, showfliers=False)
ax = sns.stripplot(data=boxplot_data_filtered, x='cell_type', y='proportion', hue='zone',
                   dodge=True, jitter=True, alpha=0.5, size=3, palette={'inner': 'blue', 'outer': 'red'})

# 添加显著性标记
for i, cell_type in enumerate(top_cell_types):
    if cell_type in results_df['cell_type'].values:
        p_val = results_df[results_df['cell_type'] == cell_type]['wilcoxon_pvalue_adj'].values[0]
        
        if p_val < 0.001:
            sig = '***'
        elif p_val < 0.01:
            sig = '**'
        elif p_val < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        
        # 在箱线图上添加显著性标记
        cell_data = boxplot_data_filtered[boxplot_data_filtered['cell_type'] == cell_type]
        y_max = cell_data['proportion'].max()
        ax.text(i, y_max + 0.01, sig, ha='center', fontsize=10, fontweight='bold')

ax.set_xlabel('Cell Type (CNs_10)', fontsize=12)
ax.set_ylabel('Proportion', fontsize=12)
ax.set_title('Cell Type Proportions in Inner vs Outer Regions (Boxplot, CNs_10)', fontsize=14, fontweight='bold')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.legend(title='Region')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_dir}/combined_boxplot.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_dir}/combined_boxplot.pdf", bbox_inches='tight')
plt.close()

# 8.6 为每个细胞类型绘制单独的箱线图
print(f"\nGenerating individual boxplots for each cell type...")
print(f"Saving to directory: {boxplot_dir}")

# 为每个细胞类型创建单独的箱线图
for cell_type in valid_cell_types:
    # 提取该细胞类型的数据
    cell_data = prop_df[prop_df['cell_type'] == cell_type]
    
    if len(cell_data) < 6:  # 至少需要3个病人（每个病人有inner和outer）
        continue
    
    plt.figure(figsize=(8, 6))
    
    # 创建箱线图
    ax = sns.boxplot(data=cell_data, x='zone', y='proportion', 
                     order=['inner', 'outer'],
                     palette={'inner': 'skyblue', 'outer': 'lightcoral'},
                     width=0.5)
    
    # 添加散点图显示个体数据点
    ax = sns.stripplot(data=cell_data, x='zone', y='proportion',
                       order=['inner', 'outer'],
                       color='black', size=6, alpha=0.6, jitter=True)
    
    # 获取统计结果
    if cell_type in results_df['cell_type'].values:
        row = results_df[results_df['cell_type'] == cell_type].iloc[0]
        p_val = row['wilcoxon_pvalue_adj']
        mean_inner = row['mean_inner']
        mean_outer = row['mean_outer']
        
        # 添加p值注释
        y_max = cell_data['proportion'].max()
        y_pos = y_max + 0.02
        
        if p_val < 0.001:
            p_text = f"p = {p_val:.2e} (***)"
        elif p_val < 0.01:
            p_text = f"p = {p_val:.2e} (**)"
        elif p_val < 0.05:
            p_text = f"p = {p_val:.3f} (*)"
        else:
            p_text = f"p = {p_val:.3f} (ns)"
        
        # 添加p值标注
        plt.text(0.5, y_pos, p_text, ha='center', fontsize=12, fontweight='bold')
        
        # 添加均值标注
        plt.text(0, mean_inner + 0.005, f'{mean_inner:.3f}', ha='center', fontsize=10, color='blue')
        plt.text(1, mean_outer + 0.005, f'{mean_outer:.3f}', ha='center', fontsize=10, color='red')
    
    # 设置图形属性
    plt.xlabel('Mucus Zone', fontsize=12)
    plt.ylabel('Proportion', fontsize=12)
    plt.title(f'{cell_type} Proportions in Inner vs Outer Regions (CNs_10)', fontsize=14, fontweight='bold')
    plt.grid(axis='y', alpha=0.3)
    
    # 调整y轴范围以容纳标注
    y_min, y_max = cell_data['proportion'].min(), cell_data['proportion'].max()
    plt.ylim(y_min - 0.02, y_max + 0.05)
    
    plt.tight_layout()
    
    # 保存图片
    # 使用安全的文件名（替换特殊字符）
    safe_cell_type = cell_type.replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
    plt.savefig(f"{boxplot_dir}/boxplot_{safe_cell_type}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{boxplot_dir}/boxplot_{safe_cell_type}.pdf", bbox_inches='tight')
    plt.close()

print(f"Generated {len([f for f in os.listdir(boxplot_dir) if f.startswith('boxplot_')])} individual boxplots")

# 8.7 创建显著的细胞类型综合箱线图
if sig_cell_types:
    print("\nGenerating combined boxplot for significant cell types only...")
    
    # 只使用显著差异的细胞类型
    sig_boxplot_data = prop_df[prop_df['cell_type'].isin(sig_cell_types)].copy()
    
    # 按差异大小排序
    sig_results = results_df[results_df['wilcoxon_pvalue_adj'] < 0.05].sort_values('mean_diff', ascending=False)
    ordered_cell_types = sig_results['cell_type'].tolist()
    
    # 保持顺序
    sig_boxplot_data['cell_type'] = pd.Categorical(sig_boxplot_data['cell_type'], categories=ordered_cell_types, ordered=True)
    sig_boxplot_data = sig_boxplot_data.sort_values('cell_type')
    
    plt.figure(figsize=(max(12, len(sig_cell_types) * 0.8), 8))
    
    # 创建分组箱线图
    ax = sns.boxplot(data=sig_boxplot_data, x='cell_type', y='proportion', hue='zone',
                     palette={'inner': 'skyblue', 'outer': 'lightcoral'}, 
                     showfliers=False, width=0.8)
    
    # 添加显著性标记
    for i, cell_type in enumerate(ordered_cell_types):
        row = sig_results[sig_results['cell_type'] == cell_type].iloc[0]
        p_val = row['wilcoxon_pvalue_adj']
        
        if p_val < 0.001:
            sig = '***'
        elif p_val < 0.01:
            sig = '**'
        elif p_val < 0.05:
            sig = '*'
        
        # 在箱线图上添加显著性标记
        cell_data = sig_boxplot_data[sig_boxplot_data['cell_type'] == cell_type]
        y_max = cell_data['proportion'].max()
        ax.text(i, y_max + 0.01, sig, ha='center', fontsize=10, fontweight='bold')
        
        # 添加差异方向箭头
        mean_diff = row['mean_diff']
        if mean_diff > 0:
            arrow_text = '↑'  # inner > outer
            arrow_color = 'blue'
        else:
            arrow_text = '↓'  # inner < outer
            arrow_color = 'red'
        
        ax.text(i, y_max + 0.03, arrow_text, ha='center', fontsize=12, fontweight='bold', color=arrow_color)
    
    ax.set_xlabel('Cell Type (CNs_10)', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title('Significantly Different Cell Types: Inner vs Outer Regions (CNs_10)', fontsize=14, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.legend(title='Region', loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/significant_celltypes_boxplot.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/significant_celltypes_boxplot.pdf", bbox_inches='tight')
    plt.close()

# 9. 生成统计摘要报告
print("\n" + "="*60)
print("Statistical Summary Report")
print("="*60)

# 显著差异总结
sig_results = results_df[results_df['wilcoxon_pvalue_adj'] < 0.05].sort_values('mean_diff', ascending=False)
n_sig = len(sig_results)

print(f"\n1. Overall Information:")
print(f"   - Number of patients: {len(patients_to_analyze)}")
print(f"   - Total cell types (CNs_10): {len(all_cell_types)}")
print(f"   - Analyzed cell types: {len(valid_cell_types)}")
print(f"   - Significantly different cell types: {n_sig}")

if n_sig > 0:
    print(f"\n2. Significantly different cell types (FDR-adjusted p<0.05):")
    for idx, row in sig_results.iterrows():
        direction = "inner > outer" if row['mean_diff'] > 0 else "inner < outer"
        print(f"   - {row['cell_type']}: {direction}")
        print(f"     Mean difference: {row['mean_diff']:.4f}, p={row['wilcoxon_pvalue_adj']:.4f}")
        print(f"     Mean proportion (inner): {row['mean_inner']:.4f}, (outer): {row['mean_outer']:.4f}")

# 差异最大的细胞类型
if len(results_df) > 0:
    max_diff_idx = results_df['mean_diff'].abs().idxmax()
    max_diff_row = results_df.loc[max_diff_idx]
    print(f"\n3. Cell type with largest absolute difference: {max_diff_row['cell_type']}")
    print(f"   Absolute difference: {abs(max_diff_row['mean_diff']):.4f}")
    print(f"   Direction: {'inner > outer' if max_diff_row['mean_diff'] > 0 else 'inner < outer'}")

# 保存摘要报告
with open(f"{output_dir}/analysis_summary.txt", "w") as f:
    f.write("Cell Composition Difference Analysis Summary (CNs_10)\n")
    f.write("="*50 + "\n\n")
    f.write(f"Analysis date: {pd.Timestamp.now()}\n")
    f.write(f"Number of patients: {len(patients_to_analyze)}\n")
    f.write(f"Total cell types (CNs_10): {len(all_cell_types)}\n")
    f.write(f"Analyzed cell types: {len(valid_cell_types)}\n")
    f.write(f"Significantly different cell types (FDR-adjusted p<0.05): {n_sig}\n\n")
    
    if n_sig > 0:
        f.write("Significant differences details:\n")
        for idx, row in sig_results.iterrows():
            direction = "inner > outer" if row['mean_diff'] > 0 else "inner < outer"
            f.write(f"- {row['cell_type']}: {direction}, p={row['wilcoxon_pvalue_adj']:.4f}\n")
            f.write(f"  Mean difference: {row['mean_diff']:.4f}, inner: {row['mean_inner']:.4f}, outer: {row['mean_outer']:.4f}\n\n")
    
    f.write("\nGenerated figures:\n")
    f.write(f"- overall_composition.png/pdf: Bar plot showing average proportions of all cell types\n")
    f.write(f"- significant_paired_comparisons.png/pdf: Paired line plots for significantly different cell types\n")
    f.write(f"- difference_heatmap.png/pdf: Heatmap showing differences by patient\n")
    f.write(f"- patient_composition_stacked.png/pdf: Stacked bar plots for each patient\n")
    f.write(f"- combined_boxplot.png/pdf: Combined boxplot for all cell types\n")
    f.write(f"- significant_celltypes_boxplot.png/pdf: Boxplot for significantly different cell types only\n")
    f.write(f"- individual_boxplots/*.png/pdf: Individual boxplots for each cell type\n")

print(f"\nAnalysis complete! All results saved to {output_dir} directory:")
print(f"  1. Statistical results table: {output_dir}/paired_test_results.csv")
print(f"  2. Overall composition bar plot: {output_dir}/overall_composition.png/pdf")
if sig_cell_types:
    print(f"  3. Significant paired comparisons: {output_dir}/significant_paired_comparisons.png/pdf")
print(f"  4. Difference heatmap: {output_dir}/difference_heatmap.png/pdf")
if n_patients > 0:
    print(f"  5. Patient composition stacked bars: {output_dir}/patient_composition_stacked.png/pdf")
print(f"  6. Combined boxplot: {output_dir}/combined_boxplot.png/pdf")
if sig_cell_types:
    print(f"  7. Significant cell types boxplot: {output_dir}/significant_celltypes_boxplot.png/pdf")
print(f"  8. Individual boxplots: {boxplot_dir}/*.png/pdf")
print(f"  9. Analysis summary: {output_dir}/analysis_summary.txt")
