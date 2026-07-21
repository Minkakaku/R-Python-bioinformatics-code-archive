import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_cnv_boxplot(data_plot, output_file="cnv_score_boxplot_0.4.png"):
    """
    绘制CNV Score箱线图，与描述中的图表一致
    
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
    
    # 2. 确定聚类顺序（与描述一致）
    # 先提取所有唯一的type值
    all_types = data_plot['type'].unique()
    
    # 分离数字类型和命名类型
    numeric_types = []
    named_types = []
    
    for t in all_types:
        try:
            # 尝试转换为数字
            numeric_val = int(t)
            numeric_types.append(numeric_val)
        except (ValueError, TypeError):
            # 如果不能转换为数字，则视为命名类型
            named_types.append(str(t))
    
    # 排序数字类型
    numeric_types.sort()
    
    # 确定最终顺序：数字类型 + 命名类型（T_NK和Myeloid）
    # 根据描述，数字类型是0-32，然后是T_NK和Myeloid
    ordered_types = [str(t) for t in numeric_types] + sorted(named_types)
    
    print(f"检测到 {len(ordered_types)} 个聚类类型")
    print(f"聚类顺序: {ordered_types}")
    
    # 3. 创建图形
    plt.figure(figsize=(16, 8))  # 与描述中的尺寸一致
    
    # 4. 使用seaborn绘制箱线图
    sns.boxplot(
        data=data_plot,
        x='type',
        y='cnv_score',
        order=ordered_types,  # 使用确定的顺序
        palette='Set3',       # 使用多彩色系
        width=0.7,            # 箱体宽度
        fliersize=2,          # 离群点大小
        linewidth=1           # 线条宽度
    )
    
    # 5. 设置标题和标签（与描述一致）
    plt.title("CNV SCORE", fontsize=16, fontweight='bold', pad=20)
    plt.xlabel("Seurat clusters", fontsize=12, labelpad=10)
    plt.ylabel("CNV_score", fontsize=12, labelpad=10)
    
    # 6. 设置Y轴（与描述一致）
    plt.ylim(0, 0.10)  # Y轴范围0-0.10
    plt.yticks([i * 0.02 for i in range(6)], 
               ['0.00', '0.02', '0.04', '0.06', '0.08', '0.10'])
    
    # 7. 设置X轴标签旋转（如果标签过多）
    if len(ordered_types) > 10:
        plt.xticks(rotation=45)
    
    # 8. 添加网格线
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    
    # 9. 调整布局并保存
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"箱线图已保存为: {output_file}")
    
    # 10. 显示图表
    plt.show()
    
    # 11. 输出统计信息
    print("\n各聚类CNV分数统计:")
    stats = data_plot.groupby('type')['cnv_score'].describe()
    print(stats)

# 使用示例
if __name__ == "__main__":
    # 示例数据创建（您应该使用自己的data_plot）
    # 假设data_plot是您已有的DataFrame
    
    
    data_plot =  pd.read_csv("plot_data_new")
    
    # 绘制箱线图
    plot_cnv_boxplot(data_plot, "my_cnv_boxplot_new.png")    

