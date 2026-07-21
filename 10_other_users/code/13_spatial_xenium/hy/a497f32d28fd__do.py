import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import warnings
warnings.filterwarnings('ignore')

def main():
    # 文件路径配置
    cns_file = "CNs_10_all_cell_features.xls"  # CNs数据文件
    releate_file = "releate"                   # 分组信息文件
    
    # 读取CNs数据文件（制表符分隔）
    print("正在读取CNs数据文件...")
    try:
        cns_df = pd.read_csv(cns_file, sep='\t')
        print(f"CNs数据形状: {cns_df.shape}")
        print(f"CNs数据列名: {list(cns_df.columns)}")
    except Exception as e:
        print(f"读取CNs文件错误: {e}")
        return
    
    # 检查必要的列是否存在
    if 'patient' not in cns_df.columns:
        print("错误: CNs数据中未找到'patient'列")
        return
    if 'CNs_10' not in cns_df.columns:
        print("错误: CNs数据中未找到'CNs_10'列")
        return
    
    # 读取releate文件（假设是制表符分隔）
    print("\n正在读取releate文件...")
    try:
        # 尝试不同的分隔符
        try:
            releate_df = pd.read_csv(releate_file, sep='\t', header=None)
        except:
            try:
                releate_df = pd.read_csv(releate_file, sep='\s+', header=None)
            except:
                releate_df = pd.read_csv(releate_file, header=None)
        
        # 设置列名
        if releate_df.shape[1] >= 2:
            releate_df.columns = ['patient', 'T_group'] if releate_df.shape[1] == 2 else ['patient', 'T_group'] + [f'extra_{i}' for i in range(2, releate_df.shape[1])]
        else:
            print("错误: releate文件列数不足")
            return
        
        print(f"releate数据形状: {releate_df.shape}")
        print(f"releate数据预览:\n{releate_df.head()}")
    except Exception as e:
        print(f"读取releate文件错误: {e}")
        return
    
    # 确保patient列类型一致
    cns_df['patient'] = cns_df['patient'].astype(str)
    releate_df['patient'] = releate_df['patient'].astype(str)
    
    # 合并数据：只保留releate文件中存在的patient
    print("\n正在合并数据...")
    merged_df = cns_df.merge(releate_df[['patient', 'T_group']], on='patient', how='inner')
    merged_df.rename(columns={'T_group': 'MPCpercent'}, inplace=True)
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"MPCpercent分组分布:\n{merged_df['MPCpercent'].value_counts()}")
    
    if merged_df.shape[0] == 0:
        print("错误: 合并后无数据，请检查patient列是否匹配")
        return
    
    # 计算每个patient中不同CNs_10的比例
    print("\n正在计算CNs_10比例...")
    
    # 获取所有唯一的CNs_10值
    cns_values = sorted(merged_df['CNs_10'].unique())
    print(f"CNs_10的唯一值: {cns_values}")
    
    # 按patient和MPCpercent分组，计算每个CNs_10值的比例
    patient_proportions = []
    
    for patient in merged_df['patient'].unique():
        patient_data = merged_df[merged_df['patient'] == patient]
        group = patient_data['MPCpercent'].iloc[0]
        
        # 计算该patient中每个CNs_10值的计数
        cns_counts = patient_data['CNs_10'].value_counts()
        total_cells = len(patient_data)
        
        # 计算比例
        proportions = {f'CNs_{cns}': cns_counts.get(cns, 0) / total_cells if total_cells > 0 else 0 
                      for cns in cns_values}
        proportions['patient'] = patient
        proportions['MPCpercent'] = group
        proportions['total_cells'] = total_cells
        
        patient_proportions.append(proportions)
    
    # 转换为DataFrame
    proportion_df = pd.DataFrame(patient_proportions)
    
    # 统计学检验 - 卡方检验比较不同分组间CNs_10分布的差异
    print("\n正在进行差异分析（卡方检验）...")
    
    # 构建列联表：每个MPCpercent分组中每个CNs_10的细胞计数
    contingency_table = pd.crosstab(merged_df['MPCpercent'], merged_df['CNs_10'])
    
    print("列联表（细胞计数）:")
    print(contingency_table)
    
    # 进行卡方检验
    try:
        chi2, p_value, dof, expected = chi2_contingency(contingency_table)
        print(f"\n卡方检验结果:")
        print(f"卡方值: {chi2:.4f}")
        print(f"P值: {p_value:.4f}")
        print(f"自由度: {dof}")
        
        if p_value < 0.05:
            print("结论: 不同MPCpercent分组间的CNs_10分布存在显著差异 (p < 0.05)")
        else:
            print("结论: 不同MPCpercent分组间的CNs_10分布无显著差异 (p ≥ 0.05)")
    except Exception as e:
        print(f"卡方检验错误: {e}")
        p_value = np.nan
    
    # 计算每个MPCpercent分组中各CNs_10的平均比例
    print("\n各MPCpercent分组中CNs_10的平均比例:")
    cns_columns = [col for col in proportion_df.columns if col.startswith('CNs_')]
    mean_proportions = proportion_df.groupby('MPCpercent')[cns_columns].mean()
    print(mean_proportions)
    
    # 保存结果到Excel文件
    output_file = "CNs_10_MPCpercent_analysis.xlsx"
    print(f"\n正在保存结果到 {output_file}...")
    
    with pd.ExcelWriter(output_file) as writer:
        merged_df.to_excel(writer, sheet_name='合并后数据', index=False)
        proportion_df.to_excel(writer, sheet_name='患者比例数据', index=False)
        mean_proportions.to_excel(writer, sheet_name='分组平均比例')
        
        # 保存检验结果
        stats_df = pd.DataFrame({
            '检验方法': ['卡方检验'],
            '统计量': [chi2],
            'P值': [p_value],
            '自由度': [dof],
            '显著性': [p_value < 0.05 if not np.isnan(p_value) else False]
        })
        stats_df.to_excel(writer, sheet_name='差异分析结果', index=False)
    
    print("分析完成！")
    print(f"结果已保存到: {output_file}")
    
    # 提供进一步分析建议
    print("\n进一步分析建议:")
    print("1. 查看Excel文件中的'分组平均比例'工作表，了解各分组的CNs_10分布")
    print("2. 如果卡方检验显著，可进行事后检验分析具体哪些CNs_10值在不同分组间差异显著")
    print("3. 可以考虑使用可视化方法（如堆叠柱状图）展示比例差异")

if __name__ == "__main__":
    main()
