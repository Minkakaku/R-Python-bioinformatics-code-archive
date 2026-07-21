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

# 1. 加载数据
print("Loading data...")
adata = sc.read_h5ad("extracted_samples.h5ad")
print(f"Data shape: {adata.shape}")

# 检查必要的列是否存在
required_cols = ['patient', 'mucus_zone', 'CNs_10']
missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    print(f"Error: Missing required columns: {missing_cols}")
    exit()

print(f"Patients: {adata.obs['patient'].nunique()}")
print(f"Samples: {adata.obs['sample'].nunique()}")
print(f"Cell types (CNs_10): {len(adata.obs['CNs_10'].unique())}")

# 2. 数据预处理
print("\nPreprocessing data...")
# 筛选inner和outer区域
adata_filtered = adata[adata.obs['mucus_zone'].isin(['inner', 'outer']), :].copy()
print(f"Filtered data shape: {adata_filtered.shape}")

# 检查每个病人在两个区域是否有数据
patient_zone_counts = adata_filtered.obs.groupby(['patient', 'mucus_zone']).size().unstack(fill_value=0)
patients_with_both = patient_zone_counts.dropna().index.tolist()
print(f"Patients with both inner and outer {len(patients_with_both)}")

if len(patients_with_both) < 2:
    print("Error: At least 2 patients with complete data are required")
    exit()

# 3. 计算细胞类型比例
print("\nCalculating cell type proportions...")

# 获取所有细胞类型
all_cell_types = adata_filtered.obs['CNs_10'].unique()
print(f"Total cell types: {len(all_cell_types)}")

# 计算每个病人在每个区域的细胞类型比例
prop_data = []
patients_to_analyze = []

for patient in patients_with_both:
    patient_data = adata_filtered[adata_filtered.obs['patient'] == patient]
    
    inner_data = patient_data[patient_data.obs['mucus_zone'] == 'inner']
    outer_data = patient_data[patient_data.obs['mucus_zone'] == 'outer']
    
    if len(inner_data) == 0 or len(outer_data) == 0:
        continue
    
    inner_total = len(inner_data)
    outer_total = len(outer_data)
    
    inner_counts = inner_data.obs['CNs_10'].value_counts()
    outer_counts = outer_data.obs['CNs_10'].value_counts()
    
    for cell_type in all_cell_types:
        inner_prop = inner_counts.get(cell_type, 0) / inner_total
        outer_prop = outer_counts.get(cell_type, 0) / outer_total
        
        prop_data.extend([
            {'patient': patient, 'cell_type': cell_type, 'zone': 'inner', 'proportion': inner_prop},
            {'patient': patient, 'cell_type': cell_type, 'zone': 'outer', 'proportion': outer_prop}
        ])
    
    patients_to_analyze.append(patient)

prop_df = pd.DataFrame(prop_data)
print(f"Final patients analyzed: {len(patients_to_analyze)}")

# 4. 统计检验
print("\nPerforming statistical tests...")

# 筛选在至少3个病人中存在的细胞类型
cell_type_presence = prop_df.groupby('cell_type').apply(
    lambda x: x['patient'].nunique() >= 3
)
valid_cell_types = cell_type_presence[cell_type_presence].index.tolist()
print(f"Cell types present in at least 3 patients: {len(valid_cell_types)}")

# 执行配对检验
test_results = []

for cell_type in valid_cell_types:
    cell_data = prop_df[prop_df['cell_type'] == cell_type]
    
    if len(cell_data) < 6:  # 至少需要3个病人
        continue
    
    # 转换为宽格式
    wide_data = cell_data.pivot(index='patient', columns='zone', values='proportion')
    
    inner_vals = wide_data['inner'].values
    outer_vals = wide_data['outer'].values
    
    # 描述性统计
    mean_inner = np.mean(inner_vals)
    mean_outer = np.mean(outer_vals)
    std_inner = np.std(inner_vals, ddof=1)
    std_outer = np.std(outer_vals, ddof=1)
    
    # 统计检验
    t_stat, t_p = stats.ttest_rel(inner_vals, outer_vals)
    w_stat, w_p = stats.wilcoxon(inner_vals, outer_vals)
    
    test_results.append({
        'cell_type': cell_type,
        'n_patients': len(wide_data),
        'mean_inner': mean_inner,
        'mean_outer': mean_outer,
        'std_inner': std_inner,
        'std_outer': std_outer,
        'mean_diff': mean_inner - mean_outer,
        't_statistic': t_stat,
        't_pvalue': t_p,
        'wilcoxon_stat': w_stat,
        'wilcoxon_pvalue': w_p
    })

if not test_results:
    print("Error: No valid cell types for statistical analysis")
    exit()

results_df = pd.DataFrame(test_results)

# 多重检验校正
results_df['t_pvalue_adj'] = multipletests(results_df['t_pvalue'], method='fdr_bh')[1]
results_df['wilcoxon_pvalue_adj'] = multipletests(results_df['wilcoxon_pvalue'], method='fdr_bh')[1]

# 排序并保存
results_df = results_df.sort_values('wilcoxon_pvalue_adj')
results_df.to_csv(f"{output_dir}/statistical_results.csv", index=False)
print("Statistical results saved")

# 5. 可视化
print("\nGenerating visualizations...")

# 设置样式
plt.style.use('default')
sns.set_palette("husl")

# 5.1 总体组成条形图
print("Creating overall composition plot...")

overall_props = prop_df.groupby(['cell_type', 'zone'])['proportion'].mean().unstack()
overall_props = overall_props.loc[valid_cell_types]
overall_props['diff'] = overall_props['inner'] - overall_props['outer']
overall_props = overall_props.sort_values('diff', ascending=False)

fig, ax = plt.subplots(figsize=(max(12, len(overall_props)/2), 8))
x = np.arange(len(overall_props))
width = 0.35

bars1 = ax.bar(x - width/2, overall_props['inner'], width, label='Inner', alpha=0.8, color='skyblue')
bars2 = ax.bar(x + width/2, overall_props['outer'], width, label='Outer', alpha=0.8, color='lightcoral')

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

ax.set_xlabel('Cell Type (CNs_10)')
ax.set_ylabel('Average Proportion')
ax.set_title('Cell Type Composition: Inner vs Outer Regions')
ax.set_xticks(x)
ax.set_xticklabels(overall_props.index, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_dir}/overall_composition.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_dir}/overall_composition.pdf", bbox_inches='tight')
plt.close()

# 5.2 组合箱线图
print("Creating combined boxplot...")

boxplot_data = prop_df[prop_df['cell_type'].isin(valid_cell_types)].copy()

# 按平均比例排序
cell_type_means = boxplot_data.groupby('cell_type')['proportion'].mean().sort_values(ascending=False)
ordered_cell_types = cell_type_means.index.tolist()

# 动态调整图形大小
n_cell_types = len(ordered_cell_types)
fig_width = max(16, n_cell_types * 0.8)

fig, ax = plt.subplots(figsize=(fig_width, 8))
sns.boxplot(data=boxplot_data, x='cell_type', y='proportion', hue='zone',
            order=ordered_cell_types,
            palette={'inner': 'skyblue', 'outer': 'lightcoral'},
            ax=ax)

ax.set_xlabel('Cell Type (CNs_10)')
ax.set_ylabel('Proportion')
ax.set_title('Cell Type Distribution: Inner vs Outer Regions')
ax.tick_params(axis='x', rotation=45)
ax.legend(title='Region')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_dir}/combined_boxplot.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_dir}/combined_boxplot.pdf", bbox_inches='tight')
plt.close()

# 5.3 单个细胞类型箱线图
print("Creating individual boxplots...")

individual_dir = f"{output_dir}/individual_boxplots"
if not os.path.exists(individual_dir):
    os.makedirs(individual_dir)

for cell_type in valid_cell_types:
    cell_data = prop_df[prop_df['cell_type'] == cell_type]
    
    if len(cell_data) < 6:
        continue
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sns.boxplot(data=cell_data, x='zone', y='proportion', 
                order=['inner', 'outer'],
                palette={'inner': 'skyblue', 'outer': 'lightcoral'},
                ax=ax)
    
    sns.stripplot(data=cell_data, x='zone', y='proportion',
                  order=['inner', 'outer'],
                  color='black', size=4, alpha=0.6, jitter=True,
                  ax=ax)
    
    # 添加统计信息
    if cell_type in results_df['cell_type'].values:
        row = results_df[results_df['cell_type'] == cell_type].iloc[0]
        p_val = row['wilcoxon_pvalue_adj']
        
        if p_val < 0.001:
            p_text = f"p = {p_val:.2e} (***)"
        elif p_val < 0.01:
            p_text = f"p = {p_val:.2e} (**)"
        elif p_val < 0.05:
            p_text = f"p = {p_val:.3f} (*)"
        else:
            p_text = f"p = {p_val:.3f} (ns)"
        
        # 添加p值标注
        y_max = cell_data['proportion'].max()
        ax.text(0.5, y_max + 0.02, p_text, ha='center', fontweight='bold')
    
    ax.set_xlabel('Region')
    ax.set_ylabel('Proportion')
    ax.set_title(f'{cell_type} Distribution')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # 安全文件名
    safe_name = cell_type.replace('/', '_').replace('\\', '_').replace(':', '_')
    safe_name = safe_name.replace('*', '_').replace('?', '_').replace('"', '_')
    safe_name = safe_name.replace('<', '_').replace('>', '_').replace('|', '_')
    
    plt.savefig(f"{individual_dir}/boxplot_{safe_name}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{individual_dir}/boxplot_{safe_name}.pdf", bbox_inches='tight')
    plt.close()

print(f"Generated {len([f for f in os.listdir(individual_dir) if f.startswith('boxplot_')])} individual boxplots")

# 5.4 热图 - 显示细胞类型区域偏好性
print("Creating enrichment heatmap...")

# 计算每个细胞类型在不同区域的Z-score
enrichment_data = []
for cell_type in valid_cell_types:
    cell_data = prop_df[prop_df['cell_type'] == cell_type]
    
    if len(cell_data) == 0:
        continue
    
    inner_mean = cell_data[cell_data['zone'] == 'inner']['proportion'].mean()
    outer_mean = cell_data[cell_data['zone'] == 'outer']['proportion'].mean()
    
    # 计算Z-score
    means = [inner_mean, outer_mean]
    if np.std(means) > 0:
        z_scores = (means - np.mean(means)) / np.std(means)
    else:
        z_scores = [0, 0]
    
    enrichment_data.append({
        'cell_type': cell_type,
        'inner_zscore': z_scores[0],
        'outer_zscore': z_scores[1]
    })

if enrichment_data:
    enrichment_df = pd.DataFrame(enrichment_data)
    enrichment_df = enrichment_df.set_index('cell_type')
    
    # 按inner z-score排序
    enrichment_df = enrichment_df.sort_values('inner_zscore', ascending=False)
    
    fig, ax = plt.subplots(figsize=(10, max(8, len(enrichment_df)/3)))
    sns.heatmap(enrichment_df, annot=True, cmap='RdBu_r', center=0,
                ax=ax, cbar_kws={'label': 'Z-score'})
    ax.set_title('Cell Type Regional Enrichment (Z-score)')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/enrichment_heatmap.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/enrichment_heatmap.pdf", bbox_inches='tight')
    plt.close()

# 5.5 病人级别的组成图
print("Creating patient-level composition plots...")

# 为每个病人创建堆叠条形图
n_patients = len(patients_to_analyze)
n_cols = min(3, n_patients)
n_rows = (n_patients + n_cols - 1) // n_cols

if n_patients > 0:
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    for i, patient in enumerate(patients_to_analyze):
        if i < len(axes):
            ax = axes[i]
            patient_data = prop_df[prop_df['patient'] == patient]
            
            # 转换为宽格式
            wide_data = patient_data.pivot(index='cell_type', columns='zone', values='proportion')
            wide_data = wide_data.loc[valid_cell_types]  # 只显示有效细胞类型
            
            wide_data.plot(kind='bar', stacked=True, ax=ax)
            ax.set_title(f'Patient {patient}')
            ax.set_ylabel('Proportion')
            ax.tick_params(axis='x', rotation=45)
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
    
    # 隐藏多余的子图
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)
    
    plt.suptitle('Cell Composition by Patient', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/patient_composition.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/patient_composition.pdf", bbox_inches='tight')
    plt.close()

# 6. 生成报告
print("\nGenerating analysis report...")

report_path = f"{output_dir}/analysis_report.txt"
with open(report_path, 'w') as f:
    f.write("Cell Composition Analysis Report\n")
    f.write("="*50 + "\n\n")
    
    f.write("DATA SUMMARY:\n")
    f.write(f"- Total patients: {len(patients_to_analyze)}\n")
    f.write(f"- Total cell types: {len(all_cell_types)}\n")
    f.write(f"- Analyzed cell types: {len(valid_cell_types)}\n")
    f.write(f"- Regions compared: inner vs outer\n\n")
    
    f.write("STATISTICAL RESULTS SUMMARY:\n")
    f.write(f"- Significant cell types (p < 0.05): {len(results_df[results_df['wilcoxon_pvalue_adj'] < 0.05])}\n")
    f.write(f"- Significant cell types (p < 0.01): {len(results_df[results_df['wilcoxon_pvalue_adj'] < 0.01])}\n")
    f.write(f"- Significant cell types (p < 0.001): {len(results_df[results_df['wilcoxon_pvalue_adj'] < 0.001])}\n\n")
    
    if len(results_df[results_df['wilcoxon_pvalue_adj'] < 0.05]) > 0:
        f.write("SIGNIFICANT CELL TYPES (p < 0.05):\n")
        sig_results = results_df[results_df['wilcoxon_pvalue_adj'] < 0.05].sort_values('mean_diff', ascending=False)
        
        for _, row in sig_results.iterrows():
            direction = "inner > outer" if row['mean_diff'] > 0 else "inner < outer"
            f.write(f"- {row['cell_type']}: {direction} (p={row['wilcoxon_pvalue_adj']:.4f}, diff={row['mean_diff']:.4f})\n")
    
    f.write(f"\nGENERATED FIGURES:\n")
    f.write(f"- overall_composition.png/pdf: Bar plot of average composition\n")
    f.write(f"- combined_boxplot.png/pdf: Combined boxplot of all cell types\n")
    f.write(f"- individual_boxplots/: Individual boxplots for each cell type\n")
    f.write(f"- enrichment_heatmap.png/pdf: Heatmap of regional enrichment\n")
    f.write(f"- patient_composition.png/pdf: Patient-level composition plots\n")
    f.write(f"- statistical_results.csv: Complete statistical results\n")

print("Analysis complete!")
print(f"\nAll results saved to: {output_dir}")
print(f"- Statistical results: {output_dir}/statistical_results.csv")
print(f"- Overall composition: {output_dir}/overall_composition.png/pdf")
print(f"- Combined boxplot: {output_dir}/combined_boxplot.png/pdf")
print(f"- Individual boxplots: {output_dir}/individual_boxplots/")
print(f"- Enrichment heatmap: {output_dir}/enrichment_heatmap.png/pdf")
print(f"- Patient composition: {output_dir}/patient_composition.png/pdf")
print(f"- Analysis report: {output_dir}/analysis_report.txt")
