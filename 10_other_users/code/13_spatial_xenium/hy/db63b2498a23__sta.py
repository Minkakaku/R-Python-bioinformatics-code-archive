import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']  # 支持中文显示
plt.rcParams['axes.unicode_minus'] = False  # 正确显示负号
sns.set_style("whitegrid")  # 设置样式

# 读取Excel文件
df = pd.read_csv('CNs_10_all_cell_features.csv',sep="\t")

# 检查必要的列是否存在
if 'CNs_10' not in df.columns or 'annotation' not in df.columns:
    raise ValueError("数据文件中必须包含'CNs_10'和'annotation'列")

# 统计每个CNs_10中细胞类型的数量
annotation_counts = df.groupby(['CNs_10', 'annotation']).size().reset_index(name='cell_count')

# 计算每个CNs_10的总细胞数
total_cells_per_cns = df.groupby('CNs_10').size().reset_index(name='total_cells')

# 合并数据并计算占比
result = pd.merge(annotation_counts, total_cells_per_cns, on='CNs_10')
result['percentage'] = (result['cell_count'] / result['total_cells']) * 100

# 按CNs_10和占比排序
result = result.sort_values(['CNs_10', 'percentage'], ascending=[True, False])

# 添加排名列
result['rank'] = result.groupby('CNs_10')['percentage'].rank(method='dense', ascending=False).astype(int)

# 选择并重排列
result = result[['CNs_10', 'annotation', 'cell_count', 'total_cells', 'percentage', 'rank']]

# 保存结果到CSV文件
result.to_csv('stat.csv', index=False, float_format='%.2f')

print("统计完成！结果已保存到stat.csv")
print(f"共处理了 {len(result['CNs_10'].unique())} 个CNs_10类别")
print(f"共包含 {len(result['annotation'].unique())} 种细胞类型")
print(f"总细胞数: {len(df)}")

# 创建输出目录
os.makedirs('top10_visualization', exist_ok=True)

# 为每个CNs_10类别生成Top 10细胞类型可视化
unique_cns = result['CNs_10'].unique()

print("\n开始生成Top 10细胞类型可视化图...")

for cns in unique_cns:
    try:
        # 筛选当前CNs_10的数据
        cns_data = result[result['CNs_10'] == cns].copy()
        
        # 获取Top 10细胞类型（按占比降序）
        top10_data = cns_data.head(20)
        
        # 如果细胞类型超过10个，将剩余类型合并为"其他"
        if len(cns_data) > 20:
            other_count = cns_data.iloc[20:]['cell_count'].sum()
            other_percentage = cns_data.iloc[20:]['percentage'].sum()
            other_row = pd.DataFrame({
                'CNs_10': [cns],
                'annotation': ['Others'],
                'cell_count': [other_count],
                'total_cells': [top10_data.iloc[0]['total_cells']],
                'percentage': [other_percentage],
                'rank': [21]
            })
            top10_data = pd.concat([top10_data, other_row], ignore_index=True)
        
        # 创建图形
        plt.figure(figsize=(14, 10))
        
        # 使用seaborn绘制柱状图
        ax = sns.barplot(data=top10_data, x='annotation', y='percentage', 
                         palette='viridis', alpha=0.8)
        
        # 添加数值标签
        for i, (_, row) in enumerate(top10_data.iterrows()):
            ax.text(i, row['percentage'] + 0.5, f'{row["percentage"]:.2f}%', 
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # 设置图表属性
        plt.title(f'CNs_10 {cns} - Top 20 Celltypes Proportion', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Celltype', fontsize=12, fontweight='bold')
        plt.ylabel('Proportion(%)', fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha='right')  # 旋转x轴标签
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图表
        filename = f"top10_annotations_CNs_10_{str(cns).replace(' ', '_').replace('/', '_')}.png"
        filepath = os.path.join('top10_visualization', filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        filename1 = f"top10_annotations_CNs_10_{str(cns).replace(' ', '_').replace('/', '_')}.pdf"
        filepath1 = os.path.join('top10_visualization', filename1)
        plt.savefig(filepath1, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"已生成: {filename}")
        
    except Exception as e:
        print(f"生成CNs_10 {cns} 的可视化图时出错: {str(e)}")
        plt.close()

# 生成汇总可视化图（所有CNs_10的Top 5细胞类型）
print("\n生成汇总可视化图...")

try:
    # 获取每个CNs_10的Top 5细胞类型
    top_annotations_by_cns = []
    
    for cns in unique_cns:
        cns_data = result[result['CNs_10'] == cns]
        top5 = cns_data.head(20)
        for _, row in top5.iterrows():
            top_annotations_by_cns.append({
                'CNs_10': cns,
                'annotation': row['annotation'],
                'percentage': row['percentage'],
                'rank': row['rank']
            })
    
    top_annotations_df = pd.DataFrame(top_annotations_by_cns)
    
    # 创建热图显示Top细胞类型分布
    plt.figure(figsize=(max(12, len(unique_cns) * 1.5), max(8, len(top_annotations_df['annotation'].unique()) * 0.5)))
    
    # 创建数据透视表
    pivot_data = top_annotations_df.pivot_table(index='annotation', columns='CNs_10', 
                                            values='percentage', aggfunc='first', fill_value=0)
    
    # 绘制热图
    sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='YlOrRd', 
                cbar_kws={'label': '占比 (%)'}, linewidths=0.5)
    
    plt.title('CNs_10 Tpo10 Celltype Heatmap', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('CNs_10', fontsize=12, fontweight='bold')
    plt.ylabel('Celltype', fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    # 保存热图
    plt.savefig('top10_visualization/summary_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('top10_visualization/summary_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("已生成汇总热图: summary_heatmap.png")
    
except Exception as e:
    print(f"生成汇总可视化图时出错: {str(e)}")
    plt.close()

print(f"\n所有可视化图已保存到目录: top10_visualization")
print(f"共处理了 {len(unique_cns)} 个CNs_10类别")
print("脚本执行完成！")
