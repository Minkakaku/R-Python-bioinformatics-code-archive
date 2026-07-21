import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

def process_cn_distribution(csv_path, output_dir="cn_distribution_tables"):
    """
    处理CN分布数据，为每个CN值生成细胞类型在样本中的分布表
    
    参数:
    csv_path: CSV文件路径
    output_dir: 输出目录
    """
    print(f"正在读取数据文件: {csv_path}")
    
    # 读取CSV文件
    df = pd.read_csv(csv_path, sep='\t')
    
    # 检查必要的列
    required_cols = ['CNs_10', 'sample', 'annotation', 'cell_id']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"错误: 文件中缺少必要的列: {missing_cols}")
        print(f"文件中可用的列: {df.columns.tolist()}")
        return None
    
    print(f"数据形状: {df.shape}")
    print(f"总细胞数: {len(df)}")
    
    # 获取所有唯一的CN值
    cn_values = df['CNs_10'].dropna().unique()
    print(f"\n发现 {len(cn_values)} 个不同的CN值:")
    for i, cn in enumerate(sorted(cn_values), 1):
        print(f"  {i}. {cn}")
    
    # 获取所有样本
    all_samples = sorted(df['sample'].dropna().unique())
    print(f"\n发现 {len(all_samples)} 个样本: {all_samples[:10]}{'...' if len(all_samples) > 10 else ''}")
    
    # 获取所有细胞类型
    all_celltypes = sorted(df['annotation'].dropna().unique())
    print(f"发现 {len(all_celltypes)} 种细胞类型: {all_celltypes[:10]}{'...' if len(all_celltypes) > 10 else ''}")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 为每个CN值生成分布表
    for cn_value in sorted(cn_values):
        print(f"\n=== 处理 CN: {cn_value} ===")
        
        # 筛选当前CN的数据
        cn_df = df[df['CNs_10'] == cn_value].copy()
        
        if len(cn_df) == 0:
            print(f"  CN {cn_value} 没有数据，跳过")
            continue
        
        print(f"  CN {cn_value} 的细胞总数: {len(cn_df)}")
        
        # 创建交叉表：样本 × 细胞类型
        crosstab = pd.crosstab(
            cn_df['sample'], 
            cn_df['annotation'],
            values=cn_df['cell_id'], 
            aggfunc='count',
            dropna=False
        )
        
        # 确保所有样本都在表中
        missing_samples = set(all_samples) - set(crosstab.index)
        for sample in missing_samples:
            crosstab.loc[sample] = 0
        
        # 确保所有细胞类型都在表中
        missing_celltypes = set(all_celltypes) - set(crosstab.columns)
        for celltype in missing_celltypes:
            crosstab[celltype] = 0
        
        # 重新排序行和列
        crosstab = crosstab.reindex(index=all_samples, columns=all_celltypes, fill_value=0)
        
        # 添加总计行和列
        crosstab['样本总计'] = crosstab.sum(axis=1)
        total_row = crosstab.sum(axis=0)
        total_row.name = '细胞类型总计'
        
        # 修复：使用pd.concat替代append
        total_row_df = total_row.to_frame().T
        crosstab = pd.concat([crosstab, total_row_df])
        
        # 添加百分比列（每个样本中各细胞类型的比例）
        for col in all_celltypes:
            if crosstab['样本总计'].iloc[:-1].sum() > 0:  # 排除总计行
                crosstab[f'{col}_百分比'] = crosstab[col] / crosstab['样本总计'] * 100
            else:
                crosstab[f'{col}_百分比'] = 0
        
        # 保存到CSV文件
        output_filename = f"CN_{str(cn_value).replace('/', '_').replace(' ', '_')}_celltype_distribution.csv"
        output_path = os.path.join(output_dir, output_filename)
        crosstab.to_csv(output_path, encoding='utf-8-sig')
        
        print(f"  已保存: {output_path}")
        print(f"  表格形状: {crosstab.shape}")
        print(f"  样本数量: {len(all_samples)}")
        print(f"  细胞类型数量: {len(all_celltypes)}")
        
        # 显示前几行和前几列
        print("\n  表格预览 (前5行，前5列):")
        preview_rows = min(5, len(crosstab))
        preview_cols = min(5, len(crosstab.columns))
        preview = crosstab.iloc[:preview_rows, :preview_cols]
        print(preview)
    
    return cn_values

def create_summary_statistics(df, output_dir="cn_distribution_tables"):
    """
    创建汇总统计信息
    """
    print(f"\n=== 创建汇总统计 ===")
    
    # 1. 按CN统计总体分布
    cn_summary = df['CNs_10'].value_counts().reset_index()
    cn_summary.columns = ['CN', '细胞数量']
    cn_summary['百分比'] = cn_summary['细胞数量'] / len(df) * 100
    cn_summary = cn_summary.sort_values('细胞数量', ascending=False)
    
    # 2. 按样本统计CN分布
    sample_cn_dist = pd.crosstab(df['sample'], df['CNs_10'])
    
    # 3. 按细胞类型统计CN分布
    celltype_cn_dist = pd.crosstab(df['annotation'], df['CNs_10'])
    
    # 保存汇总统计
    cn_summary_path = os.path.join(output_dir, "CN_overall_summary.csv")
    sample_cn_path = os.path.join(output_dir, "sample_CN_distribution.csv")
    celltype_cn_path = os.path.join(output_dir, "celltype_CN_distribution.csv")
    
    cn_summary.to_csv(cn_summary_path, encoding='utf-8-sig', index=False)
    sample_cn_dist.to_csv(sample_cn_path, encoding='utf-8-sig')
    celltype_cn_dist.to_csv(celltype_cn_path, encoding='utf-8-sig')
    
    print(f"CN总体统计已保存: {cn_summary_path}")
    print(f"样本-CN分布已保存: {sample_cn_path}")
    print(f"细胞类型-CN分布已保存: {celltype_cn_path}")
    
    # 打印简要统计
    print(f"\nCN分布统计:")
    for _, row in cn_summary.iterrows():
        print(f"  {row['CN']}: {row['细胞数量']} 个细胞 ({row['百分比']:.2f}%)")
    
    return cn_summary, sample_cn_dist, celltype_cn_dist

def create_visualizations(df, cn_summary, output_dir="cn_distribution_tables"):
    """
    创建可视化图表
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        print(f"\n=== 创建可视化图表 ===")
        
        # 设置中文字体
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
        
        # 1. CN分布条形图
        plt.figure(figsize=(12, 8))
        bars = plt.bar(range(len(cn_summary)), cn_summary['细胞数量'], color='skyblue', edgecolor='black')
        
        plt.title('CN分布统计', fontsize=16, pad=20)
        plt.xlabel('CN', fontsize=14)
        plt.ylabel('细胞数量', fontsize=14)
        plt.xticks(range(len(cn_summary)), cn_summary['CN'], rotation=45, ha='right')
        
        # 在每个柱子上添加数值
        for i, bar in enumerate(bars):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{int(height)}', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "CN_distribution_bar.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  CN分布条形图已保存")
        
        # 2. 样本-CN热图
        if len(df['sample'].unique()) > 1 and len(df['CNs_10'].unique()) > 1:
            plt.figure(figsize=(12, 8))
            
            # 创建样本-CN交叉表
            sample_cn_heatmap = pd.crosstab(df['sample'], df['CNs_10'])
            
            # 绘制热图
            sns.heatmap(sample_cn_heatmap, annot=True, fmt='d', cmap='YlOrRd', 
                       cbar_kws={'label': '细胞数量'})
            
            plt.title('各样本中CN分布热图', fontsize=16, pad=20)
            plt.xlabel('CN', fontsize=14)
            plt.ylabel('样本', fontsize=14)
            plt.tight_layout()
            
            plt.savefig(os.path.join(output_dir, "sample_CN_heatmap.png"), dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  样本-CN热图已保存")
        
        # 3. 细胞类型-CN堆叠条形图（前20种细胞类型）
        top_celltypes = df['annotation'].value_counts().head(20).index
        celltype_cn_data = df[df['annotation'].isin(top_celltypes)]
        
        if len(celltype_cn_data) > 0:
            plt.figure(figsize=(14, 10))
            
            # 创建交叉表
            celltype_cn_crosstab = pd.crosstab(celltype_cn_data['annotation'], celltype_cn_data['CNs_10'])
            
            # 绘制堆叠条形图
            celltype_cn_crosstab.plot(kind='bar', stacked=True, figsize=(14, 10), 
                                     colormap='tab20c')
            
            plt.title('前20种细胞类型的CN分布（堆叠条形图）', fontsize=16, pad=20)
            plt.xlabel('细胞类型', fontsize=14)
            plt.ylabel('细胞数量', fontsize=14)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='CN', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            plt.savefig(os.path.join(output_dir, "celltype_CN_stacked_bar.png"), dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  细胞类型-CN堆叠条形图已保存")
        
        print(f"  所有可视化图表已保存到: {output_dir}")
        
    except Exception as e:
        print(f"  创建可视化图表时出错: {e}")
        import traceback
        traceback.print_exc()

def main():
    # 文件路径配置
    csv_path = "CNs_10_all_cell_features.csv"  # 替换为您的CSV文件路径
    output_dir = "cn_distribution_tables"  # 输出目录
    
    print("=" * 60)
    print("CN分布统计工具")
    print("=" * 60)
    
    try:
        # 读取数据
        print(f"正在读取文件: {csv_path}")
        df = pd.read_csv(csv_path,sep="\t")
        print(f"数据读取成功: {df.shape[0]} 行, {df.shape[1]} 列")
        
        # 检查必要的列
        required_cols = ['CNs_10', 'sample', 'annotation', 'cell_id']
        for col in required_cols:
            if col not in df.columns:
                print(f"错误: 缺少必要列 '{col}'")
                print(f"可用的列: {list(df.columns)}")
                return
        
        # 处理每个CN的分布
        cn_values = process_cn_distribution(csv_path, output_dir)
        
        if cn_values is not None:
            # 创建汇总统计
            cn_summary, sample_cn_dist, celltype_cn_dist = create_summary_statistics(df, output_dir)
            
            # 创建可视化图表
            create_visualizations(df, cn_summary, output_dir)
            
            print(f"\n=== 处理完成 ===")
            print(f"所有结果已保存到目录: {output_dir}")
            print(f"生成了 {len(cn_values)} 个CN分布表")
            
            # 显示汇总信息
            print(f"\n汇总统计:")
            print(f"  总细胞数: {len(df)}")
            print(f"  样本数量: {df['sample'].nunique()}")
            print(f"  细胞类型数量: {df['annotation'].nunique()}")
            print(f"  CN值数量: {len(cn_values)}")
            
            # 显示最多的CN
            top_cn = df['CNs_10'].value_counts().head(5)
            print(f"\n最多的5个CN:")
            for cn, count in top_cn.items():
                percentage = count / len(df) * 100
                print(f"  {cn}: {count} 个细胞 ({percentage:.1f}%)")
    
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
