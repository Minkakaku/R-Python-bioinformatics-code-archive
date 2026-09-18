import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def analyze_band_composition(adata):
    """
    分析site_band区域中细胞组成的差异
    """
    # 1. 筛选site_band为in, out, on的细胞
    band_mask = adata.obs['site_band'].isin(['in', 'out', 'on'])
    adata_bands = adata[band_mask, :].copy()
    
    print(f"筛选后细胞数量: {adata_bands.n_obs}")
    print(f"各区域细胞分布:")
    print(adata_bands.obs['site_band'].value_counts())
    
    # 2. 创建统计表格
    count_table = pd.crosstab(
        index=[adata_bands.obs['patient'], adata_bands.obs['site_band']],
        columns=adata_bands.obs['annotation']
    )
    
    # 计算百分比
    percentage_table = count_table.div(count_table.sum(axis=1), axis=0) * 100
    
    print("\n各区域细胞类型数量统计:")
    print(count_table)
    
    return adata_bands, count_table, percentage_table

def perform_statistical_tests(adata_bands, count_table):
    """
    执行配对样本的显著性检验 - 不使用scikit-posthocs
    """
    results = {}
    
    cell_types = adata_bands.obs['annotation'].unique()
    patients = adata_bands.obs['patient'].unique()
    bands = ['in', 'on', 'out']
    
    print(f"\n进行统计分析:")
    print(f"细胞类型数量: {len(cell_types)}")
    print(f"患者数量: {len(patients)}")
    
    for cell_type in cell_types:
        print(f"\n=== 细胞类型: {cell_type} ===")
        
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
            
            # 只有当患者在至少2个区域有数据时才包括
            if valid_bands >= 2:
                data.append(patient_data)
                valid_patients.append(patient)
        
        if len(data) < 3:  # 至少需要3个患者
            print(f"数据不足，跳过 {cell_type}")
            results[cell_type] = None
            continue
            
        data = np.array(data)
        
        # Friedman检验（配对非参数检验）
        try:
            # 只选择有数据的列
            valid_columns = [i for i in range(data.shape[1]) if np.any(data[:, i] > 0)]
            valid_data = [data[:, i] for i in valid_columns]
            
            if len(valid_data) >= 2:
                friedman_stat, friedman_p = stats.friedmanchisquare(*valid_data)
                print(f"Friedman检验: 统计量={friedman_stat:.3f}, p值={friedman_p:.3f}")
                
                # 如果Friedman检验显著，进行事后检验（使用Wilcoxon符号秩检验）
                posthoc_results = None
                if friedman_p < 0.05 and len(valid_data) > 1:
                    posthoc_results = {}
                    valid_band_names = [bands[i] for i in valid_columns]
                    
                    # 进行两两比较
                    for i in range(len(valid_band_names)):
                        for j in range(i+1, len(valid_band_names)):
                            band1 = valid_band_names[i]
                            band2 = valid_band_names[j]
                            
                            # 提取两个区域的数据
                            data1 = data[:, i]
                            data2 = data[:, j]
                            
                            # 只保留在两个区域都有数据的患者
                            valid_pairs = []
                            for k in range(len(data1)):
                                if data1[k] > 0 or data2[k] > 0:  # 至少在一个区域有数据
                                    valid_pairs.append((data1[k], data2[k]))
                            
                            if len(valid_pairs) >= 2:
                                # 使用Wilcoxon符号秩检验
                                data1_vals = [p[0] for p in valid_pairs]
                                data2_vals = [p[1] for p in valid_pairs]
                                
                                # 检查是否有差异
                                if len(set(data1_vals)) > 1 or len(set(data2_vals)) > 1:
                                    try:
                                        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(data1_vals, data2_vals)
                                        posthoc_results[f"{band1}_vs_{band2}"] = {
                                            'statistic': wilcoxon_stat,
                                            'p_value': wilcoxon_p
                                        }
                                        print(f"  {band1} vs {band2}: Wilcoxon统计量={wilcoxon_stat:.3f}, p值={wilcoxon_p:.3f}")
                                    except:
                                        print(f"  {band1} vs {band2}: 无法计算Wilcoxon检验")
                
                results[cell_type] = {
                    'friedman_stat': friedman_stat,
                    'friedman_p': friedman_p,
                    'posthoc': posthoc_results,
                    'data': data,
                    'patients': valid_patients,
                    'bands': [bands[i] for i in valid_columns]
                }
        except Exception as e:
            print(f"统计检验错误: {e}")
            results[cell_type] = None
    
    return results

def visualize_composition(adata_bands, percentage_table, results):
    """
    可视化细胞组成结果
    """
    # 设置图形样式
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. 总体堆叠条形图
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    
    band_percentage = percentage_table.groupby(level='site_band').mean()
    band_percentage.plot(kind='bar', stacked=True, ax=ax1)
    ax1.set_title('各区域细胞类型平均组成', fontsize=14, fontweight='bold')
    ax1.set_ylabel('百分比 (%)')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.show()
    
    # 2. 每个患者的堆叠条形图
    unique_patients = percentage_table.index.get_level_values('patient').unique()
    n_patients = len(unique_patients)
    
    # 动态计算子图布局
    n_cols = min(3, n_patients)
    n_rows = (n_patients + n_cols - 1) // n_cols
    
    if n_rows > 0 and n_cols > 0:
        fig2, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        elif n_rows == 1 or n_cols == 1:
            axes = axes.reshape(-1)
        
        for i, patient in enumerate(unique_patients):
            if i < len(axes):
                ax = axes[i] if n_patients > 1 else axes
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
                    ax.set_title(f'患者 {patient}', fontweight='bold')
                    ax.set_ylabel('百分比 (%)')
                    ax.tick_params(axis='x', rotation=45)
            
            # 隐藏多余的子图
            for j in range(i+1, len(axes)):
                if n_patients > 1:
                    axes[j].set_visible(False)
        
        plt.tight_layout()
        plt.show()
    
    # 3. 箱线图和小提琴图
    fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(20, 8))
    
    # 准备箱线图数据
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
    
    if boxplot_data:
        boxplot_df = pd.DataFrame(boxplot_data)
        
        # 选择前10个最常见的细胞类型进行可视化
        top_cell_types = percentage_table.sum().sort_values(ascending=False).head(10).index
        
        # 箱线图
        if not boxplot_df.empty and len(top_cell_types) > 0:
            sns.boxplot(data=boxplot_df[boxplot_df['cell_type'].isin(top_cell_types)], 
                        x='cell_type', y='percentage', hue='band', ax=ax3)
            ax3.set_title('各细胞类型在不同区域的分布', fontsize=14, fontweight='bold')
            ax3.set_ylabel('百分比 (%)')
            ax3.tick_params(axis='x', rotation=45)
            ax3.legend(title='区域')
            
            # 小提琴图
            sns.violinplot(data=boxplot_df[boxplot_df['cell_type'].isin(top_cell_types)], 
                           x='cell_type', y='percentage', hue='band', ax=ax4)
            ax4.set_title('各细胞类型在不同区域的分布（小提琴图）', fontsize=14, fontweight='bold')
            ax4.set_ylabel('百分比 (%)')
            ax4.tick_params(axis='x', rotation=45)
            ax4.legend(title='区域')
        
        plt.tight_layout()
        plt.show()
    
    # 4. 热图 - 细胞类型区域偏好性
    fig4, ax5 = plt.subplots(figsize=(12, 8))
    
    # 计算每个细胞类型在不同区域的相对富集程度
    enrichment_matrix = []
    for cell_type in percentage_table.columns:
        band_means = []
        for band in ['in', 'on', 'out']:
            try:
                mean_percentage = percentage_table.xs(band, level='site_band')[cell_type].mean()
                band_means.append(mean_percentage)
            except:
                band_means.append(0)
        enrichment_matrix.append(band_means)
    
    if enrichment_matrix:
        enrichment_df = pd.DataFrame(
            enrichment_matrix,
            index=percentage_table.columns,
            columns=['in', 'on', 'out']
        )
        
        # 归一化处理（Z-score）
        enrichment_zscore = enrichment_df.apply(stats.zscore, axis=1)
        sns.heatmap(enrichment_zscore, annot=True, cmap='RdBu_r', center=0, 
                    ax=ax5, cbar_kws={'label': 'Z-score'})
        ax5.set_title('细胞类型区域偏好性热图\n(正值表示在该区域富集)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    # 5. 显著性检验结果可视化
    significant_cell_types = []
    for cell_type, result in results.items():
        if result and result['friedman_p'] < 0.05:
            significant_cell_types.append(cell_type)
    
    if significant_cell_types:
        fig5, ax6 = plt.subplots(figsize=(12, 8))
        
        # 准备显著性检验结果数据
        sig_data = []
        for cell_type in significant_cell_types[:10]:  # 只显示前10个显著的类型
            result = results[cell_type]
            for i, band in enumerate(result['bands']):
                for j, proportion in enumerate(result['data'][:, i]):
                    if proportion > 0:  # 只包括有数据的点
                        sig_data.append({
                            'cell_type': cell_type,
                            'band': band,
                            'proportion': proportion
                        })
        
        if sig_data:
            sig_df = pd.DataFrame(sig_data)
            
            # 创建带有显著性标记的箱线图
            sns.boxplot(data=sig_df, x='cell_type', y='proportion', hue='band', ax=ax6)
            ax6.set_title('显著差异的细胞类型在不同区域的分布\n(Friedman检验 p < 0.05)', 
                         fontsize=14, fontweight='bold')
            ax6.set_ylabel('比例')
            ax6.tick_params(axis='x', rotation=45)
            ax6.legend(title='区域')
            
            plt.tight_layout()
            plt.show()

def save_statistical_results(adata_bands, count_table, percentage_table, results):
    """
    保存统计分析结果
    """
    # 保存统计表格
    count_table.to_csv('band_cell_count_statistics.csv')
    percentage_table.to_csv('band_cell_percentage_statistics.csv')
    
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
        results_df.to_csv('statistical_test_results.csv', index=False)
    
    print("\n统计分析结果已保存:")
    print("- band_cell_count_statistics.csv: 细胞数量统计")
    print("- band_cell_percentage_statistics.csv: 细胞百分比统计")
    print("- statistical_test_results.csv: 显著性检验结果")

# 主函数
def main():
    # 加载数据（请替换为您的文件路径）
    # adata = sc.read_h5ad('your_data.h5ad')
    
    # 检查必要的列是否存在
    required_cols = ['patient', 'site_band', 'annotation']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        raise ValueError(f"缺少必要的列: {missing_cols}")
    
    print("开始统计分析...")
    
    # 1. 数据筛选和基本统计
    adata_bands, count_table, percentage_table = analyze_band_composition(adata)
    
    # 2. 显著性检验
    results = perform_statistical_tests(adata_bands, count_table)
    
    # 3. 可视化
    visualize_composition(adata_bands, percentage_table, results)
    
    # 4. 保存结果
    save_statistical_results(adata_bands, count_table, percentage_table, results)
    
    print("\n分析完成!")

# 简化的一站式函数
def analyze_band_composition_simple(adata):
    """
    简化的一站式分析函数
    """
    # 数据筛选和基本统计
    adata_bands, count_table, percentage_table = analyze_band_composition(adata)
    
    # 显著性检验
    results = perform_statistical_tests(adata_bands, count_table)
    
    # 可视化
    visualize_composition(adata_bands, percentage_table, results)
    
    # 保存结果
    save_statistical_results(adata_bands, count_table, percentage_table, results)
    
    return adata_bands, count_table, percentage_table, results

# 使用示例
if __name__ == "__main__":
    # 加载您的数据
    adata = sc.read_h5ad("adata_with_band.h5ad")
    
    # 运行分析
    adata_bands, count_table, percentage_table, results = analyze_band_composition_simple(adata)
