import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 设置中文字体（可选，解决中文显示问题）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']  # 使用黑体或Arial
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

def plot_cluster_cnv_boxplot(file_path, output_file='cnv_score_boxplot.png'):
    """
    绘制不同cluster的cnv_score分布箱线图
    
    参数:
    file_path: 输入文件路径
    output_file: 输出图片文件名
    """
    
    # 1. 读取数据[6,7](@ref)
    try:
        data = pd.read_csv(file_path)  # 默认以逗号分隔
        print(f"成功读取数据，共 {len(data)} 行")
        print(f"数据列名: {list(data.columns)}")
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        return
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return
    
    # 2. 数据验证和预处理
    # 检查必要的列是否存在[6](@ref)
    required_columns = ['cell', 'cluster', 'cnv_score']
    missing_columns = [col for col in required_columns if col not in data.columns]
    
    if missing_columns:
        print(f"错误: 缺少必要的列: {missing_columns}")
        print(f"现有列: {list(data.columns)}")
        return
    
    # 检查数据质量
    print(f"Cluster数量: {data['cluster'].nunique()}")
    print(f"CNV分数范围: {data['cnv_score'].min():.3f} ~ {data['cnv_score'].max():.3f}")
    print(f"缺失值数量: {data.isnull().sum().sum()}")
    
    # 去除缺失值[6](@ref)
    data_clean = data.dropna(subset=['cluster', 'cnv_score'])
    if len(data_clean) < len(data):
        print(f"已去除 {len(data) - len(data_clean)} 行包含缺失值的记录")
    
    # 3. 绘制箱线图[2,3](@ref)
    plt.figure(figsize=(12, 8))  # 设置图形大小
    
    # 方法1: 使用seaborn绘制箱线图（推荐，更美观）[2](@ref)
    sns.boxplot(
        x='cluster', 
        y='cnv_score', 
        data=data_clean,
        palette='viridis',  # 配色方案
        width=0.7,          # 箱体宽度
        showfliers=True,    # 显示离群点
        fliersize=3         # 离群点大小
    )
    
    # 4. 美化图表[3](@ref)
    plt.title('CNV Score Distribution by Cluster', fontsize=16, fontweight='bold')
    plt.xlabel('Cluster', fontsize=12)
    plt.ylabel('CNV Score', fontsize=12)
    plt.xticks(rotation=45)  # 如果cluster名称较长，旋转x轴标签
    
    # 添加网格线便于阅读
    plt.grid(axis='y', alpha=0.3)
    
    # 调整布局防止标签被裁剪
    plt.tight_layout()
    
    # 5. 保存和显示图表[3](@ref)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"箱线图已保存为: {output_file}")
    
    # 显示图表
    plt.show()
    
    # 6. 可选: 输出统计摘要
    print("\n各Cluster的CNV分数统计摘要:")
    stats_summary = data_clean.groupby('cluster')['cnv_score'].describe()
    print(stats_summary)
    
    return data_clean

# 使用示例
if __name__ == "__main__":
    # 假设您的文件名为 'a'（无扩展名）
    file_path = 'a'  # 如果文件有扩展名，如 'a.csv'，请相应修改
    
    # 运行函数
    data_result = plot_cluster_cnv_boxplot(file_path)
