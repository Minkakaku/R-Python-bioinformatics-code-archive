import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

def load_and_merge_data(main_file, relate_file):
    """
    加载主数据文件和分级文件，并合并数据
    """
    print("正在读取主数据文件...")
    # 读取主Excel文件
    df = pd.read_csv(main_file,sep="\t")
    print(f"主数据形状: {df.shape}")
    
    print("\n正在读取分级文件...")
    # 读取分级文件（假设为CSV格式）
    relate_df = pd.read_csv(relate_file,sep="\t")
    print(f"分级数据形状: {relate_df.shape}")
    print(f"分级文件列名: {list(relate_df.columns)}")
    
    # 检查分级文件列名
    if len(relate_df.columns) < 2:
        raise ValueError("分级文件应包含至少两列：patient和分级")
    
    # 假设第一列是patient，第二列是分级
    patient_col = relate_df.columns[0]
    grading_col = relate_df.columns[1]
    
    print(f"使用列: {patient_col} 作为患者ID, {grading_col} 作为分级")
    
    # 合并数据，只保留主数据中在分级文件中存在的patient
    merged_df = df.merge(relate_df[[patient_col, grading_col]], 
                        left_on='patient', right_on=patient_col, how='inner')
    
    # 重命名分级列为tie
    merged_df = merged_df.rename(columns={grading_col: 'tie'})
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"分级分布:\n{merged_df['tie'].value_counts()}")
    print(f"CNs_10的唯一值数量: {merged_df['CNs_10'].nunique()}")
    
    return merged_df

def create_boxplots_by_cns(df, output_dir="cns_boxplots"):
    """
    为每个CNs_10类别创建箱线图，按tie分组显示
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有唯一的CNs_10值
    unique_cns = df['CNs_10'].unique()
    print(f"\n将为 {len(unique_cns)} 个CNs_10类别创建箱线图...")
    
    # 为每个CNs_10类别创建箱线图
    for i, cns_value in enumerate(unique_cns):
        try:
            # 筛选当前CNs_10的数据
            cns_data = df[df['CNs_10'] == cns_value]
            
            # 检查是否有足够的数据
            if len(cns_data) < 5:
                print(f"CNs_10={cns_value} 的数据量不足({len(cns_data)}个)，跳过")
                continue
            
            # 创建图形
            plt.figure(figsize=(10, 6))
            
            # 使用Seaborn绘制箱线图
            sns.boxplot(data=cns_data, x='tie', y='CNs_10', palette='Set2')
            
            # 添加数据点（抖动散点图）
            sns.stripplot(data=cns_data, x='tie', y='CNs_10', 
                         color='black', alpha=0.5, size=3, jitter=True)
            
            # 设置图表属性
            plt.title(f'CNs_10 = {cns_value} (tie) Proportion\n(Sample number: {len(cns_data)})', 
                     fontsize=14, fontweight='bold', pad=20)
            plt.xlabel('MPC class', fontsize=12, fontweight='bold')
            plt.ylabel('CNs_10 ', fontsize=12, fontweight='bold')
            
            # 添加网格线
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            
            # 调整布局
            plt.tight_layout()
            
            # 保存图像
            filename = f"boxplot_CNs_10_{str(cns_value).replace(' ', '_').replace('/', '_')}.png"
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"已生成: {filename} (样本量: {len(cns_data)})")
            
        except Exception as e:
            print(f"生成CNs_10={cns_value}的箱线图时出错: {str(e)}")
            plt.close()
    
    return output_dir

def create_summary_boxplot(df, output_dir):
    """
    创建所有CNs_10类别的汇总箱线图
    """
    print("\n正在生成汇总箱线图...")
    
    try:
        # 创建图形
        plt.figure(figsize=(14, 8))
        
        # 使用Seaborn绘制箱线图
        sns.boxplot(data=df, x='CNs_10', y='CNs_10', hue='tie', palette='Set2')
        
        # 设置图表属性
        plt.title('CNs_10 proportions', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('CNs_10 ', fontsize=12, fontweight='bold')
        plt.ylabel('CNs_10 ', fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        
        # 添加图例
        plt.legend(title='MPC class', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # 添加网格线
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图像
        filepath = os.path.join(output_dir, "summary_boxplot_all_CNs_10.png")
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"已生成汇总箱线图: summary_boxplot_all_CNs_10.png")
        
    except Exception as e:
        print(f"生成汇总箱线图时出错: {str(e)}")
        plt.close()

def main():
    """
    主函数
    """
    # 配置参数
    main_file = "CNs_10_all_cell_features.xls"
    relate_file = "releate.csv"  # 假设为CSV文件
    
    try:
        # 1. 读取和合并数据
        df = load_and_merge_data(main_file, relate_file)
        
        # 检查必要的列是否存在
        if 'tie' not in df.columns:
            raise ValueError("合并后数据中未找到'tie'列")
        if 'CNs_10' not in df.columns:
            raise ValueError("数据中未找到'CNs_10'列")
        
        print(f"\n数据基本信息:")
        print(f"总样本数: {len(df)}")
        print(f"分级(tie)的唯一值: {sorted(df['tie'].unique())}")
        print(f"CNs_10的唯一值: {sorted(df['CNs_10'].unique())}")
        
        # 2. 为每个CNs_10类别创建箱线图
        output_dir = create_boxplots_by_cns(df)
        
        # 3. 创建汇总箱线图
        create_summary_boxplot(df, output_dir)
        
        # 4. 保存统计摘要
        stats_summary = df.groupby(['CNs_10', 'tie'])['CNs_10'].agg(['count', 'mean', 'std', 'min', 'max']).round(4)
        stats_filepath = os.path.join(output_dir, "cns_boxplots_statistics.xlsx")
        stats_summary.to_excel(stats_filepath)
        
        print(f"\n分析完成！")
        print(f"所有箱线图已保存到目录: {output_dir}")
        print(f"统计摘要已保存到: {output_dir}/cns_boxplots_statistics.xlsx")
        print(f"共处理了 {df['CNs_10'].nunique()} 个CNs_10类别")
        print(f"总样本数: {len(df)}")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
