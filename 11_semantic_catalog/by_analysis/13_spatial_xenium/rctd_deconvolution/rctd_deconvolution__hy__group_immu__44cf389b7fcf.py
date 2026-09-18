import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from statannotations import Annotator
import warnings
warnings.filterwarnings('ignore')

def analyze_celltype_composition(h5ad_path, before_patients, after_patients, 
                                target_celltypes, annotation_col='annotation', 
                                patient_col='patient', output_h5ad="Immu_unpaired.h5ad"):
    """
    分析before和after组之间指定细胞类型的占比差异
    
    参数:
    h5ad_path: h5ad文件路径
    before_patients: before组patient列表
    after_patients: after组patient列表  
    target_celltypes: 目标细胞类型列表
    annotation_col: 细胞类型列名
    patient_col: patient列名
    output_h5ad: 输出h5ad文件名
    """
    
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列是否存在
    if patient_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{patient_col}'列")
    
    if annotation_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{annotation_col}'列")
    
    # 合并所有目标patient
    all_patients = before_patients + after_patients
    existing_patients = [p for p in all_patients if p in adata.obs[patient_col].values]
    missing_patients = [p for p in all_patients if p not in adata.obs[patient_col].values]
    
    print(f"\n=== Patient匹配情况 ===")
    print(f"目标patient数量: {len(all_patients)}")
    print(f"找到的patient数量: {len(existing_patients)}")
    print(f"未找到的patient: {missing_patients}")
    
    if len(existing_patients) < 2:
        print("错误: 没有找到足够的patient!")
        return None
    
    # 提取目标patient数据
    adata_subset = adata[adata.obs[patient_col].isin(existing_patients), :].copy()
    
    print(f"\n=== 数据提取结果 ===")
    print(f"提取后的细胞数: {adata_subset.n_obs}")
    
    # 添加分组信息
    adata_subset.obs['group'] = 'unknown'
    adata_subset.obs.loc[adata_subset.obs[patient_col].isin(before_patients), 'group'] = 'before'
    adata_subset.obs.loc[adata_subset.obs[patient_col].isin(after_patients), 'group'] = 'after'
    
    # 统计每个patient的细胞数量
    patient_counts = adata_subset.obs[patient_col].value_counts()
    print(f"\n=== 各Patient细胞数量统计 ===")
    for patient in existing_patients:
        count = patient_counts[patient]
        group = 'before' if patient in before_patients else 'after'
        print(f"Patient {patient} ({group}组): {count} 个细胞")
    
    # 检查指定的细胞类型是否存在
    all_celltypes = adata_subset.obs[annotation_col].unique()
    existing_celltypes = [ct for ct in target_celltypes if ct in all_celltypes]
    missing_celltypes = [ct for ct in target_celltypes if ct not in all_celltypes]
    
    print(f"\n=== 细胞类型匹配情况 ===")
    print(f"指定的细胞类型: {len(target_celltypes)} 种")
    print(f"找到的细胞类型: {len(existing_celltypes)} 种")
    print(f"未找到的细胞类型: {missing_celltypes}")
    
    if len(existing_celltypes) == 0:
        print("错误: 没有找到任何指定的细胞类型!")
        return None
    
    # 计算每个patient中各种细胞类型的占比
    celltype_composition = []
    
    for patient in existing_patients:
        patient_data = adata_subset[adata_subset.obs[patient_col] == patient]
        total_cells = len(patient_data)
        group = 'before' if patient in before_patients else 'after'
        
        # 计算每个细胞类型的数量和占比
        celltype_counts = patient_data.obs[annotation_col].value_counts()
        
        for celltype in existing_celltypes:
            count = celltype_counts.get(celltype, 0)
            percentage = (count / total_cells) * 100 if total_cells > 0 else 0
            
            celltype_composition.append({
                'patient': patient,
                'group': group,
                'celltype': celltype,
                'count': count,
                'percentage': percentage,
                'total_cells': total_cells
            })
    
    # 转换为DataFrame
    comp_df = pd.DataFrame(celltype_composition)
    
    print(f"\n=== 细胞类型占比统计 ===")
    print(comp_df.groupby(['group', 'celltype'])['percentage'].describe())
    
    # 进行显著性检验并绘制箱线图
    plot_celltype_comparison(comp_df, existing_celltypes)
    
    # 保存提取后的数据
    if adata_subset.raw is None:
        adata_subset.raw = adata_subset
    
    adata_subset.write(output_h5ad)
    print(f"\n提取后的数据已保存为: {output_h5ad}")
    
    return adata_subset, comp_df

def plot_celltype_comparison(comp_df, celltypes):
    """
    绘制细胞类型占比的箱线图并进行显著性检验
    """
    # 设置图形样式
    plt.style.use('default')
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 10
    
    # 计算需要的行数和列数
    n_celltypes = len(celltypes)
    n_cols = min(4, n_celltypes)  # 最多4列
    n_rows = (n_celltypes + n_cols - 1) // n_cols  # 向上取整
    
    # 创建图形
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    
    # 如果只有一行，确保axes是二维数组
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # 颜色设置
    colors = ['#1f77b4', '#ff7f0e']  # before和after的颜色
    
    # 对每个细胞类型进行绘图和检验
    for idx, celltype in enumerate(celltypes):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # 提取当前细胞类型的数据
        celltype_data = comp_df[comp_df['celltype'] == celltype]
        
        # 准备绘图数据
        plot_data = []
        groups = ['before', 'after']
        
        for group in groups:
            group_data = celltype_data[celltype_data['group'] == group]['percentage']
            if len(group_data) > 0:
                plot_data.append(group_data)
        
        # 绘制箱线图和小提琴图
        if len(plot_data) == 2:  # 两组数据都存在
            # 小提琴图
            parts = ax.violinplot(plot_data, showmeans=False, showmedians=True)
            
            # 设置颜色
            for pc, color in zip(parts['bodies'], colors):
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
            
            for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
                if partname in parts:
                    parts[partname].set_color('black')
                    parts[partname].set_linewidth(1)
            
            # 添加散点图显示每个样本
            for i, group in enumerate(groups):
                group_data = celltype_data[celltype_data['group'] == group]['percentage']
                x = np.random.normal(i + 1, 0.04, size=len(group_data))
                ax.scatter(x, group_data, alpha=0.7, color=colors[i], s=60, 
                          edgecolors='white', linewidth=1)
            
            # 进行显著性检验（Mann-Whitney U检验）
            before_data = celltype_data[celltype_data['group'] == 'before']['percentage']
            after_data = celltype_data[celltype_data['group'] == 'after']['percentage']
            
            if len(before_data) > 1 and len(after_data) > 1:
                try:
                    # Mann-Whitney U检验
                    stat, p_value = stats.mannwhitneyu(before_data, after_data, 
                                                      alternative='two-sided')
                    
                    # 根据p值确定显著性标记
                    if p_value < 0.001:
                        sig_symbol = '***'
                    elif p_value < 0.01:
                        sig_symbol = '**'
                    elif p_value < 0.05:
                        sig_symbol = '*'
                    else:
                        sig_symbol = 'ns'
                    
                    # 添加显著性标记
                    max_val = max(before_data.max(), after_data.max())
                    y_pos = max_val + max_val * 0.1
                    
                    # 绘制显著性线
                    ax.plot([1, 1, 2, 2], [y_pos, y_pos + max_val * 0.05, 
                                          y_pos + max_val * 0.05, y_pos], 
                           lw=1.5, c='black')
                    
                    # 添加显著性符号
                    ax.text(1.5, y_pos + max_val * 0.07, sig_symbol, 
                           ha='center', va='bottom', color='black', fontsize=12, 
                           fontweight='bold')
                    
                    # 在标题中显示p值
                    ax.set_title(f'{celltype}\n(p={p_value:.4f})', fontsize=11)
                    
                except Exception as e:
                    print(f"对细胞类型 {celltype} 进行显著性检验时出错: {e}")
                    ax.set_title(celltype, fontsize=11)
            else:
                ax.set_title(celltype, fontsize=11)
        
        # 设置图表属性
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['I_II', 'III'], rotation=45)
        ax.set_ylabel('Percentage (%)')
        ax.grid(True, alpha=0.3)
    
    # 隐藏多余的子图
    for idx in range(len(celltypes), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        fig.delaxes(axes[row, col])
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors[0], alpha=0.6, label='I_II'),
        Patch(facecolor=colors[1], alpha=0.6, label='III')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2, 
               bbox_to_anchor=(0.5, 0.02), fontsize=12)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)  # 为图例留出空间
    
    # 保存图片
    output_filename = "Immu_comparison_I_II_III.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    output_filename1 = "Immu_comparison_I_II_III.pdf"
    plt.savefig(output_filename1, dpi=300, bbox_inches='tight')
    print(f"\n比较图已保存为: {output_filename}")
    
    plt.show()

# 使用示例
if __name__ == "__main__":
    # 文件路径
    h5ad_path = "adata_celltype_RCTD_removed_with_patient.h5ad"  # 请替换为您的h5ad文件路径
    
    # 定义分组
    before_patients = ["ZF2304103-8","ZF2208790-1","ZF2207170-8","ZF2202916-8","Z2101519-1","Z2171766-3","Z2262137-5","Z2271577-12","Z2224197-10","Z2332034-5","Z2336500-6","Z2345459-7","Z2323818-8","Q2236500-1","F1805049-2","Z1530765-11","Z2038176-15","Z2353385-10","Z1954693-9","Z1635989-12","Z1642635-14","Z1802235-10","Z1814675-2","Z1826868-6","Z1915172-7","Z1952876-9","F1706967-11","Z2123354-1","Z2123354-2","Z2417954-8","Z1430967-9"]
    after_patients = ["Z2057500-35","Z1736838-15","Z1944219-14","Q2406275-2","Q2324341-1"]    
# 指定要分析的细胞类型
    target_celltypes = ["CD8+_Tem","M2","MKI67+_T","B_follicular","CD4+_Tex","LYVE1+_TAM","CD8+_Tcm","NK","IFN-T","B_germinal_center","EMT_M2","Treg","CD8+_Tex","CD4+_Tem","B_memory","TAM"]
    
    try:
        result = analyze_celltype_composition(
            h5ad_path, before_patients, after_patients, target_celltypes
        )
        
        if result is not None:
            adata_subset, comp_df = result
            
            # 显示详细统计信息
            print(f"\n=== 详细统计信息 ===")
            
            # 按组统计
            for group in ['before', 'after']:
                group_data = comp_df[comp_df['group'] == group]
                group_patients = group_data['patient'].unique()
                
                print(f"\n{group.upper()} 组 ({len(group_patients)}个patient):")
                for patient in group_patients:
                    patient_data = group_data[group_data['patient'] == patient]
                    total_cells = patient_data.iloc[0]['total_cells'] if len(patient_data) > 0 else 0
                    
                    print(f"  {patient}: {total_cells} 个细胞")
                    for _, row in patient_data.iterrows():
                        if row['count'] > 0:  # 只显示存在的细胞类型
                            print(f"    {row['celltype']}: {row['count']} 个细胞 ({row['percentage']:.2f}%)")
            
            # 保存统计结果到CSV
            comp_df.to_csv("Immu_composition_statistics.csv", index=False)
            print(f"\n详细统计结果已保存为: celltype_composition_statistics.csv")
            
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
