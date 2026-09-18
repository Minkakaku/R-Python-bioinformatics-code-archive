import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
from datetime import datetime
warnings.filterwarnings('ignore')

def create_output_folder(folder_name='plot_new_band_300'):
    """创建输出文件夹"""
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Created output folder: {folder_name}")
    return folder_name

def save_figure(fig, filename, folder='plot_new_band_300'):
    """保存图形为PNG和PDF格式"""
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    # 保存PNG
    png_path = os.path.join(folder, f"{filename}.png")
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"Saved PNG: {png_path}")
    
    # 保存PDF
    pdf_path = os.path.join(folder, f"{filename}.pdf")
    fig.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved PDF: {pdf_path}")
    
    plt.close(fig)

def analyze_band_composition(adata):
    """分析site_band区域中细胞组成的差异"""
    band_mask = adata.obs['site_band'].isin(['in', 'out', 'on'])
    adata_bands = adata[band_mask, :].copy()
    
    print(f"Filtered cell count: {adata_bands.n_obs}")
    print(f"Distribution of cells in different bands:")
    print(adata_bands.obs['site_band'].value_counts())
    
    # 创建统计表格
    count_table = pd.crosstab(
        index=[adata_bands.obs['patient'], adata_bands.obs['site_band']],
        columns=adata_bands.obs['CNs_10']
    )
    
    # 计算百分比
    percentage_table = count_table.div(count_table.sum(axis=1), axis=0) * 100
    
    print("\nCell type count statistics by region:")
    print(count_table.head())
    
    return adata_bands, count_table, percentage_table

def perform_statistical_tests(adata_bands, count_table):
    """执行配对样本的显著性检验"""
    results = {}
    
    cell_types = adata_bands.obs['CNs_10'].unique()
    patients = adata_bands.obs['patient'].unique()
    bands = ['in', 'on', 'out']
    
    print(f"\nPerforming statistical analysis:")
    print(f"Number of cell types: {len(cell_types)}")
    print(f"Number of patients: {len(patients)}")
    
    for cell_type in cell_types:
        print(f"\n=== Cell type: {cell_type} ===")
        
        # 准备数据：每个患者在三个区域的比例
        data = []
        valid_patients = []
        for patient in patients:
            patient_data = []
            valid_bands = 0
            for band in bands:
                try:
                    cell_count = count_table.loc[(patient, band), cell_type]
                    total_count = count_table.loc[(patient, band)].sum()
                    proportion = cell_count / total_count if total_count > 0 else 0
                    patient_data.append(proportion)
                    if proportion > 0:
                        valid_bands += 1
                except:
                    patient_data.append(0)
            
            if valid_bands >= 2:
                data.append(patient_data)
                valid_patients.append(patient)
        
        if len(data) < 3:
            print(f"Insufficient data, skipping {cell_type}")
            results[cell_type] = None
            continue
            
        data = np.array(data)
        
        # Friedman检验
        try:
            valid_columns = [i for i in range(data.shape[1]) if np.any(data[:, i] > 0)]
            valid_data = [data[:, i] for i in valid_columns]
            
            if len(valid_data) >= 2:
                friedman_stat, friedman_p = stats.friedmanchisquare(*valid_data)
                print(f"Friedman test: statistic={friedman_stat:.3f}, p-value={friedman_p:.3f}")
                
                posthoc_results = None
                if friedman_p < 0.05 and len(valid_data) > 1:
                    posthoc_results = {}
                    valid_band_names = [bands[i] for i in valid_columns]
                    
                    for i in range(len(valid_band_names)):
                        for j in range(i+1, len(valid_band_names)):
                            band1 = valid_band_names[i]
                            band2 = valid_band_names[j]
                            
                            data1 = data[:, i]
                            data2 = data[:, j]
                            
                            valid_pairs = []
                            for k in range(len(data1)):
                                if data1[k] > 0 or data2[k] > 0:
                                    valid_pairs.append((data1[k], data2[k]))
                            
                            if len(valid_pairs) >= 2:
                                data1_vals = [p[0] for p in valid_pairs]
                                data2_vals = [p[1] for p in valid_pairs]
                                
                                if len(set(data1_vals)) > 1 or len(set(data2_vals)) > 1:
                                    try:
                                        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(data1_vals, data2_vals)
                                        posthoc_results[f"{band1}_vs_{band2}"] = {
                                            'statistic': wilcoxon_stat,
                                            'p_value': wilcoxon_p
                                        }
                                        print(f"  {band1} vs {band2}: Wilcoxon statistic={wilcoxon_stat:.3f}, p-value={wilcoxon_p:.3f}")
                                    except:
                                        print(f"  {band1} vs {band2}: Unable to perform Wilcoxon test")
                
                results[cell_type] = {
                    'friedman_stat': friedman_stat,
                    'friedman_p': friedman_p,
                    'posthoc': posthoc_results,
                    'data': data,
                    'patients': valid_patients,
                    'bands': [bands[i] for i in valid_columns]
                }
        except Exception as e:
            print(f"Statistical test error: {e}")
            results[cell_type] = None
    
    return results

def plot_individual_celltype_boxplots(percentage_table, output_folder):
    """为每个细胞类型单独绘制箱线图"""
    celltype_folder = os.path.join(output_folder, 'individual_celltype_boxplots')
    if not os.path.exists(celltype_folder):
        os.makedirs(celltype_folder)
        print(f"Created folder for individual cell type boxplots: {celltype_folder}")
    
    # 准备数据
    boxplot_data = []
    for cell_type in percentage_table.columns:
        for band in ['in', 'on', 'out']:
            try:
                values = percentage_table.xs(band, level='site_band')[cell_type]
                for value in values:
                    boxplot_data.append({
                        'cell_type': cell_type,
                        'band': band,
                        'percentage': value
                    })
            except:
                continue
    
    if not boxplot_data:
        print("No data available for individual cell type boxplots")
        return
    
    boxplot_df = pd.DataFrame(boxplot_data)
    
    # 为每个细胞类型创建单独的箱线图
    for cell_type in percentage_table.columns:
        celltype_data = boxplot_df[boxplot_df['cell_type'] == cell_type]
        
        if len(celltype_data) == 0:
            continue
            
        # 创建图形
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # 绘制箱线图
        sns.boxplot(data=celltype_data, x='band', y='percentage', ax=ax, 
                    order=['in', 'on', 'out'])
        sns.stripplot(data=celltype_data, x='band', y='percentage', ax=ax, 
                      color='black', alpha=0.5, size=3, jitter=True,
                      order=['in', 'on', 'out'])
        
        # 设置图形属性
        ax.set_title(f'Distribution of {cell_type} Across Regions', fontsize=14, fontweight='bold')
        ax.set_ylabel('Percentage (%)')
        ax.set_xlabel('Region')
        ax.tick_params(axis='x', rotation=45)  # 修改：X轴标签旋转45度
        
        plt.tight_layout()
        
        # 保存图形
        filename = f"boxplot_{cell_type.replace(' ', '_').replace('/', '_')}"
        save_figure(fig, filename, celltype_folder)
    
    print(f"Individual cell type boxplots saved to: {celltype_folder}")

def visualize_composition(adata_bands, percentage_table, results):
    """可视化细胞组成结果"""
    plt.style.use('default')
    sns.set_palette("husl")
    
    output_folder = create_output_folder('plot_new_band_300')
    
    # 1. 总体堆叠条形图
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    
    band_percentage = percentage_table.groupby(level='site_band').mean()
    
    if not band_percentage.empty:
        band_percentage.plot(kind='bar', stacked=True, ax=ax1)
        ax1.set_title('Average Cell Type Composition by Region', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Percentage (%)')
        ax1.set_xlabel('Region')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.tick_params(axis='x', rotation=45)  # 修改：X轴标签旋转45度
        plt.tight_layout()
        
        save_figure(fig1, 'average_composition_by_region', output_folder)
    else:
        print("Warning: No data available for stacked bar chart")
    
    # 2. 每个患者的堆叠条形图
    unique_patients = percentage_table.index.get_level_values('patient').unique()
    n_patients = len(unique_patients)
    
    if n_patients > 0:
        n_cols = min(3, n_patients)
        n_rows = (n_patients + n_cols - 1) // n_cols
        
        fig2, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        
        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        elif n_rows == 1 or n_cols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        for i, patient in enumerate(unique_patients):
            if i < len(axes):
                ax = axes[i]
                patient_data = []
                for band in ['in', 'on', 'out']:
                    try:
                        patient_band_data = percentage_table.loc[(patient, band)]
                        patient_band_data.name = band
                        patient_data.append(patient_band_data)
                    except:
                        continue
                
                if patient_data:
                    patient_df = pd.DataFrame(patient_data).T
                    patient_df.plot(kind='bar', stacked=True, ax=ax)
                    ax.set_title(f'Patient {patient}', fontweight='bold')
                    ax.set_ylabel('Percentage (%)')
                    ax.set_xlabel('Cell Type')
                    ax.tick_params(axis='x', rotation=45)  # 修改：X轴标签旋转45度
                else:
                    ax.set_visible(False)
        
        for j in range(i+1, len(axes)):
            axes[j].set_visible(False)
        
        plt.tight_layout()
        save_figure(fig2, 'composition_by_patient', output_folder)
    
    # 3. 箱线图和小提琴图
    if not percentage_table.empty:
        fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(20, 8))
        
        boxplot_data = []
        for cell_type in percentage_table.columns:
            for band in ['in', 'on', 'out']:
                try:
                    values = percentage_table.xs(band, level='site_band')[cell_type]
                    for value in values:
                        boxplot_data.append({
                            'cell_type': cell_type,
                            'band': band,
                            'percentage': value
                        })
                except:
                    continue
        
        if boxplot_data:  # 修复：添加冒号和正确的条件语句
            boxplot_df = pd.DataFrame(boxplot_data)
            
            top_cell_types = percentage_table.sum().sort_values(ascending=False).head(10).index
            
            if not boxplot_df.empty and len(top_cell_types) > 0:
                filtered_df = boxplot_df[boxplot_df['cell_type'].isin(top_cell_types)]
                if not filtered_df.empty:
                    sns.boxplot(data=filtered_df, x='cell_type', y='percentage', hue='band', ax=ax3)
                    ax3.set_title('Cell Type Distribution Across Regions (Boxplot)', fontsize=14, fontweight='bold')
                    ax3.set_ylabel('Percentage (%)')
                    ax3.set_xlabel('Cell Type')
                    ax3.tick_params(axis='x', rotation=45)  # 修改：X轴标签旋转45度
                    ax3.legend(title='Region')
                    
                    sns.violinplot(data=filtered_df, x='cell_type', y='percentage', hue='band', ax=ax4)
                    ax4.set_title('Cell Type Distribution Across Regions (Violin)', fontsize=14, fontweight='bold')
                    ax4.set_ylabel('Percentage (%)')
                    ax4.set_xlabel('Cell Type')
                    ax4.tick_params(axis='x', rotation=45)  # 修改：X轴标签旋转45度
                    ax4.legend(title='Region')
            
            plt.tight_layout()
            save_figure(fig3, 'distribution_boxplot_violin', output_folder)
    
    # 4. 为每个细胞类型单独绘制箱线图
    plot_individual_celltype_boxplots(percentage_table, output_folder)
    
    # 5. 热图
    if not percentage_table.empty:
        fig4, ax5 = plt.subplots(figsize=(12, 8))
        
        enrichment_matrix = []
        cell_type_list = []
        for cell_type in percentage_table.columns:
            band_means = []
            for band in ['in', 'on', 'out']:
                try:
                    mean_percentage = percentage_table.xs(band, level='site_band')[cell_type].mean()
                    band_means.append(mean_percentage)
                except:
                    band_means.append(0)
            enrichment_matrix.append(band_means)
            cell_type_list.append(cell_type)
        
        if enrichment_matrix:
            enrichment_df = pd.DataFrame(
                enrichment_matrix,
                index=cell_type_list,
                columns=['in', 'on', 'out']
            )
            
            enrichment_zscore = pd.DataFrame(index=enrichment_df.index, columns=enrichment_df.columns)
            
            for idx in enrichment_df.index:
                row = enrichment_df.loc[idx]
                if row.std() == 0:
                    z_row = pd.Series([0, 0, 0], index=enrichment_df.columns)
                else:
                    z_row = (row - row.mean()) / row.std()
                enrichment_zscore.loc[idx] = z_row
            
            enrichment_zscore = enrichment_zscore.astype(float)
            
            if not enrichment_zscore.empty and not enrichment_zscore.isnull().all().all():
                try:
                    sns.heatmap(enrichment_zscore, annot=True, cmap='RdBu_r', center=0, 
                                ax=ax5, cbar_kws={'label': 'Z-score'}, fmt='.2f')
                    ax5.set_title('Cell Type Regional Enrichment Heatmap\n(Positive values indicate enrichment)', 
                                  fontsize=14, fontweight='bold')
                    ax5.set_xlabel('Region')
                    ax5.set_ylabel('Cell Type')
                    # 热图的X轴标签旋转45度
                    ax5.tick_params(axis='x', rotation=45)
                    ax5.tick_params(axis='y', rotation=0)
                except Exception as e:
                    print(f"Heatmap plotting error: {e}")
                    ax5.text(0.5, 0.5, 'Unable to plot heatmap', ha='center', va='center', transform=ax5.transAxes)
                    ax5.set_title('Cell Type Regional Enrichment Heatmap')
            else:
                ax5.text(0.5, 0.5, 'Insufficient data for heatmap', ha='center', va='center', transform=ax5.transAxes)
                ax5.set_title('Cell Type Regional Enrichment Heatmap')
        
        plt.tight_layout()
        save_figure(fig4, 'regional_enrichment_heatmap', output_folder)
    
    # 6. 显著性检验结果可视化
    if results:
        significant_cell_types = []
        for cell_type, result in results.items():
            if result and result['friedman_p'] < 0.05:
                significant_cell_types.append(cell_type)
        
        if significant_cell_types:
            fig5, ax6 = plt.subplots(figsize=(12, 8))
            
            sig_data = []
            display_cell_types = significant_cell_types[:10]
            for cell_type in display_cell_types:
                result = results[cell_type]
                if result and 'bands' in result and 'data' in result:
                    for i, band in enumerate(result['bands']):
                        for j, proportion in enumerate(result['data'][:, i]):
                            if proportion > 0:
                                sig_data.append({
                                    'cell_type': cell_type,
                                    'band': band,
                                    'proportion': proportion
                                })
            
            if sig_data:
                sig_df = pd.DataFrame(sig_data)
                
                sns.boxplot(data=sig_df, x='cell_type', y='proportion', hue='band', ax=ax6)
                ax6.set_title('Significantly Different Cell Types Across Regions\n(Friedman test p < 0.05)', 
                             fontsize=14, fontweight='bold')
                ax6.set_ylabel('Proportion')
                ax6.set_xlabel('Cell Type')
                ax6.tick_params(axis='x', rotation=45)  # 修改：X轴标签旋转45度
                ax6.legend(title='Region')
                
                plt.tight_layout()
                save_figure(fig5, 'significant_cell_types', output_folder)

def generate_analysis_report(adata_bands, count_table, percentage_table, results, output_folder):
    """生成分析报告"""
    report_path = os.path.join(output_folder, 'analysis_report.txt')
    
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("CELL COMPOSITION ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        # 报告基本信息
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total cells analyzed: {adata_bands.n_obs}\n")
        f.write(f"Number of patients: {len(adata_bands.obs['patient'].unique())}\n")
        f.write(f"Number of cell types: {len(adata_bands.obs['CNs_10'].unique())}\n")
        f.write(f"Regions analyzed: in, on, out\n\n")
        
        # 细胞分布统计
        f.write("CELL DISTRIBUTION BY REGION:\n")
        f.write("-" * 40 + "\n")
        region_counts = adata_bands.obs['site_band'].value_counts()
        for region, count in region_counts.items():
            f.write(f"{region}: {count} cells ({count/adata_bands.n_obs*100:.1f}%)\n")
        f.write("\n")
        
        # 显著性检验结果
        f.write("STATISTICAL TEST RESULTS:\n")
        f.write("-" * 40 + "\n")
        
        significant_cell_types = []
        for cell_type, result in results.items():
            if result:
                f.write(f"\n{cell_type}:\n")
                f.write(f"  Friedman test: statistic={result['friedman_stat']:.3f}, p-value={result['friedman_p']:.3f}\n")
                if result['friedman_p'] < 0.05:
                    f.write("  **Significant difference between regions (p < 0.05)**\n")
                    significant_cell_types.append(cell_type)
                    
                    if result['posthoc']:
                        f.write("  Post-hoc comparisons:\n")
                        for comparison, posthoc_result in result['posthoc'].items():
                            f.write(f"    {comparison}: Wilcoxon statistic={posthoc_result['statistic']:.3f}, p-value={posthoc_result['p_value']:.3f}\n")
                else:
                    f.write("  No significant difference between regions\n")
        
        f.write(f"\nTotal significant cell types: {len(significant_cell_types)}\n")
        if significant_cell_types:
            f.write("Significant cell types: " + ", ".join(significant_cell_types) + "\n")
        
        # 细胞类型组成摘要
        f.write("\nCELL TYPE COMPOSITION SUMMARY:\n")
        f.write("-" * 40 + "\n")
        for region in ['in', 'on', 'out']:
            try:
                region_data = percentage_table.xs(region, level='site_band')
                f.write(f"\n{region.upper()} REGION:\n")
                for cell_type in percentage_table.columns:
                    mean_percentage = region_data[cell_type].mean()
                    std_percentage = region_data[cell_type].std()
                    f.write(f"  {cell_type}: {mean_percentage:.2f}% ± {std_percentage:.2f}%\n")
            except:
                f.write(f"\n{region.upper()} REGION: No data available\n")
    
    print(f"Analysis report saved to: {report_path}")
    return report_path

def save_statistical_results(adata_bands, count_table, percentage_table, results):
    """保存统计分析结果"""
    output_folder = create_output_folder('plot_new_band_300_cn')
    
    # 保存统计表格
    count_table.to_csv(os.path.join(output_folder, 'band_cell_count_statistics.csv'))
    percentage_table.to_csv(os.path.join(output_folder, 'band_cell_percentage_statistics.csv'))
    
    # 保存显著性检验结果
    results_df = []
    for cell_type, result in results.items():
        if result:
            results_df.append({
                'cell_type': cell_type,
                'friedman_statistic': result['friedman_stat'],
                'friedman_p_value': result['friedman_p'],
                'significant': result['friedman_p'] < 0.05
            })
    
    if results_df:
        results_df = pd.DataFrame(results_df)
        results_df.to_csv(os.path.join(output_folder, 'statistical_test_results.csv'), index=False)
    
    # 生成分析报告
    generate_analysis_report(adata_bands, count_table, percentage_table, results, output_folder)
    
    print(f"\nStatistical analysis results saved to folder: {output_folder}")
    print(f"- band_cell_count_statistics.csv: Cell count statistics")
    print(f"- band_cell_percentage_statistics.csv: Cell percentage statistics")
    print(f"- statistical_test_results.csv: Statistical test results")
    print(f"- analysis_report.txt: Comprehensive analysis report")

def analyze_band_composition_simple(adata):
    """简化的一站式分析函数"""
    required_cols = ['patient', 'site_band', 'CNs_10']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    print("Starting statistical analysis...")
    
    # 1. 数据筛选和基本统计
    adata_bands, count_table, percentage_table = analyze_band_composition(adata)
    
    # 2. 显著性检验
    results = perform_statistical_tests(adata_bands, count_table)
    
    # 3. 可视化
    visualize_composition(adata_bands, percentage_table, results)
    
    # 4. 保存结果
    save_statistical_results(adata_bands, count_table, percentage_table, results)
    
    print("\nAnalysis completed!")
    
    return adata_bands, count_table, percentage_table, results

# 主程序
if __name__ == "__main__":
    print("Loading data file...")
    try:
        adata = sc.read_h5ad("adata_with_band_300.h5ad")
        print(f"Data loaded successfully, cell count: {adata.n_obs}")
        print(f"Data columns: {list(adata.obs.columns)}")
    except FileNotFoundError:
        print("Error: File 'adata_with_band.h5ad' not found")
        print("Please ensure the file exists in the current directory or provide the correct file path")
        exit(1)
    except Exception as e:
        print(f"Error reading {e}")
        exit(1)
    
    # 运行分析
    adata_bands, count_table, percentage_table, results = analyze_band_composition_simple(adata)
    
    print("\nAll analysis completed! Check the 'plot_new_band_300' folder for output files.")
