import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description='统计单细胞数据中细胞类型数量')
    parser.add_argument('--input', required=True, help='输入h5ad文件路径')
    parser.add_argument('--output', default='cell_type_counts.csv', help='输出CSV文件路径')
    parser.add_argument('--plot', action='store_true', help='是否生成可视化图表')
    args = parser.parse_args()

    try:
        # 读取数据文件
        adata = sc.read_h5ad(args.input)
        
        # 检查cell_type列是否存在
        if 'cell_type' not in adata.obs.columns:
            raise ValueError("adata.obs中不存在'cell_type'列！请检查列名")
        
        # 统计细胞类型数量
        cell_type_counts = adata.obs['cell_type'].value_counts().reset_index()
        cell_type_counts.columns = ['CellType', 'Count']
        
        # 保存结果
        cell_type_counts.to_csv(args.output, index=False)
        print(f"✅ 细胞类型统计结果已保存至: {os.path.abspath(args.output)}")
        print("细胞类型分布统计:")
        print(cell_type_counts)
        
        # 可选：生成可视化
        if args.plot:
            plot_path = os.path.splitext(args.output)[0] + '.png'
            generate_plot(cell_type_counts, plot_path)
            print(f" 可视化图表已保存至: {os.path.abspath(plot_path)}")
            
    except FileNotFoundError:
        print(f"❌ 错误：文件不存在 - {args.input}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 发生错误: {str(e)}")
        sys.exit(1)

def generate_plot(counts_df, save_path):
    """生成两种类型的可视化图表"""
    plt.figure(figsize=(15, 6))
    
    # 水平条形图
    plt.subplot(1, 2, 1)
    counts_sorted = counts_df.sort_values('Count', ascending=True)
    plt.barh(counts_sorted['CellType'], counts_sorted['Count'], color='#4c72b0')
    plt.xlabel('细胞数量', fontsize=12)
    plt.ylabel('细胞类型', fontsize=12)
    plt.title('细胞类型分布', fontsize=14)
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    
    # 饼图
    plt.subplot(1, 2, 2)
    explode = [0.05] * len(counts_df)  # 轻微分离各扇形
    plt.pie(counts_df['Count'], 
            labels=counts_df['CellType'],
            autopct='%1.1f%%',
            startangle=90,
            explode=explode,
            wedgeprops={'edgecolor': 'w', 'linewidth': 1})
    plt.axis('equal')
    plt.title('细胞类型比例', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main()
