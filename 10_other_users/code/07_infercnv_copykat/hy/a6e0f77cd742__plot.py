import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_cnv_boxplot(data_plot, output_file="cnv_score_boxplot_patient.png"):
    """
    绘制CNV Score箱线图，箱体只显示四分位数
    
    参数:
    data_plot: DataFrame，包含type和cnv_score两列
    output_file: 输出图片文件名
    """
    
    # 设置中文字体（如果需要显示中文）
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 1. 数据验证
    if data_plot is None or data_plot.empty:
        print("错误: 数据为空")
        return
    
    if 'type' not in data_plot.columns or 'cnv_score' not in data_plot.columns:
        print("错误: 数据必须包含'type'和'cnv_score'列")
        print(f"当前列名: {list(data_plot.columns)}")
        return
    
    # 2. 确定聚类顺序
    all_types = data_plot['type'].unique()
    numeric_types = []
    named_types = []
    
    for t in all_types:
        try:
            numeric_val = int(t)
            numeric_types.append(numeric_val)
        except (ValueError, TypeError):
            named_types.append(str(t))
    
    numeric_types.sort()
    ordered_types = [str(t) for t in numeric_types] + sorted(named_types)
    
    print(f"检测到 {len(ordered_types)} 个聚类类型")
    print(f"聚类顺序: {ordered_types}")
    
    # 3. 创建图形
    plt.figure(figsize=(16, 8))
    
    # 4. 使用seaborn绘制箱线图，只显示四分位数
    # 创建自定义箱线图属性
    boxprops = dict(linestyle='-', linewidth=1.5, color='darkblue')
    medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
    whiskerprops = dict(linestyle='--', linewidth=1, color='black')
    capprops = dict(linestyle='-', linewidth=1, color='black')
    
    sns.boxplot(
        data=data_plot,
        x='type',
        y='cnv_score',
        order=ordered_types,
        palette='Set3',
        width=0.7,
        fliersize=0,  # 不显示异常值点
        showmeans=False,  # 不显示均值
        showbox=True,  # 显示箱体
        showcaps=True,  # 显示须线帽
        showfliers=False,  # 不显示异常值
        boxprops=boxprops,
        medianprops=medianprops,
        whiskerprops=whiskerprops,
        capprops=capprops
    )
    
    # 5. 设置标题和标签
    plt.title("CNV SCORE", fontsize=16, fontweight='bold', pad=20)
    plt.xlabel("Seurat clusters", fontsize=12, labelpad=10)
    plt.ylabel("CNV_score", fontsize=12, labelpad=10)
    
    # 6. 设置Y轴
    plt.ylim(0, 0.10)
    plt.yticks([i * 0.02 for i in range(6)], 
               ['0.00', '0.02', '0.04', '0.06', '0.08', '0.10'])
    
    # 7. 设置X轴标签旋转
    if len(ordered_types) > 10:
        plt.xticks(rotation=45)
    
    # 8. 添加网格线
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    
    # 9. 添加四分位数标记
    # 计算每个聚类的四分位数
    quartiles = data_plot.groupby('type')['cnv_score'].quantile([0.25, 0.5, 0.75]).unstack()
    
    # 在图上标注四分位数
    for i, cluster in enumerate(ordered_types):
        if cluster in quartiles.index:
            q1 = quartiles.loc[cluster, 0.25]
            median = quartiles.loc[cluster, 0.5]
            q3 = quartiles.loc[cluster, 0.75]
            
            # 在箱体上方标注四分位数
            plt.text(
                i, 
                q3 + 0.005, 
                f'Q3: {q3:.3f}',
                ha='center',
                fontsize=8,
                color='darkblue'
            )
            
            # 在箱体下方标注中位数
            plt.text(
                i, 
                median - 0.005, 
                f'Med: {median:.3f}',
                ha='center',
                fontsize=8,
                color='firebrick'
            )
            
            # 在箱体下方标注Q1
            plt.text(
                i, 
                q1 - 0.005, 
                f'Q1: {q1:.3f}',
                ha='center',
                fontsize=8,
                color='darkblue'
            )
    
    # 10. 调整布局并保存
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"箱线图已保存为: {output_file}")
    
    # 11. 显示图表
    plt.show()
    
    # 12. 输出统计信息
    print("\n各聚类CNV分数统计:")
    stats = data_plot.groupby('type')['cnv_score'].describe()
    print(stats)

# 使用示例
if __name__ == "__main__":
    # 示例数据创建（您应该使用自己的data_plot）
    # 假设data_plot是您已有的DataFrame
    
    data_plot = pd.read_csv("plot_data")
    
    # 绘制箱线图
    plot_cnv_boxplot(data_plot, "my_cnv_boxplot_quartiles.png")
