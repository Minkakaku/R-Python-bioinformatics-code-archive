import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

def calculate_and_plot_celltype_proportions(h5ad_path, output_prefix="MPC_analysis"):
    """
    统计MPC样本的细胞类型占比并绘制堆积柱状图
    
    参数:
    h5ad_path: h5ad文件路径
    output_prefix: 输出文件前缀
    """
    
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    print(f"obs列: {adata.obs.columns.tolist()}")
    
    # 检查必要的列是否存在
    required_cols = ['group', 'patient', 'annotation']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"在adata.obs中未找到'{col}'列")
    
    # 筛选MPC样本
    print("\n筛选MPC样本...")
    mpc_adata = adata[adata.obs['group'] == 'MPC', :].copy()
    
    if mpc_adata.n_obs == 0:
        raise ValueError("没有找到MPC样本!")
    
    print(f"MPC样本细胞数: {mpc_adata.n_obs}")
    
    # 获取MPC的patient列表
    mpc_patients = mpc_adata.obs['patient'].unique()
    print(f"MPC patient数量: {len(mpc_patients)}")
    print(f"MPC patient列表: {sorted(mpc_patients.tolist())}")
    
    # 计算每个patient的细胞数量
    patient_counts = mpc_adata.obs['patient'].value_counts()
    print(f"\n各Patient细胞数量统计:")
    for patient, count in patient_counts.items():
        print(f"  {patient}: {count} 个细胞")
    
    # 统计细胞类型
    celltypes = mpc_adata.obs['annotation'].unique()
    print(f"\n细胞类型数量: {len(celltypes)}")
    print(f"细胞类型列表: {sorted(celltypes.tolist())}")
    
    # 计算每个patient中每种细胞类型的占比
    print("\n计算每个patient的细胞类型占比...")
    
    # 创建交叉表
    patient_celltype_counts = pd.crosstab(
        mpc_adata.obs['patient'], 
        mpc_adata.obs['annotation']
    )
    
    # 计算每个patient的总细胞数
    patient_totals = patient_celltype_counts.sum(axis=1)
    
    # 计算占比（百分比）
    patient_celltype_proportions = (patient_celltype_counts.div(patient_totals, axis=0) * 100).round(2)
    
    # 保存统计结果
    proportions_file = f"{output_prefix}_celltype_proportions.csv"
    patient_celltype_proportions.to_csv(proportions_file, encoding='utf-8-sig')
    print(f"细胞类型占比数据已保存到: {proportions_file}")
    
    # 保存每个patient的总细胞数
    totals_file = f"{output_prefix}_patient_totals.csv"
    patient_totals.to_csv(totals_file, encoding='utf-8-sig')
    print(f"patient总细胞数已保存到: {totals_file}")
    
    # 绘制堆积柱状图
    print("\n绘制堆积柱状图...")
    plot_stacked_bar_chart(patient_celltype_proportions, patient_totals, output_prefix)
    
    # 创建汇总统计
    create_summary_statistics(mpc_adata, patient_celltype_proportions, output_prefix)
    
    return mpc_adata, patient_celltype_proportions, patient_totals

def plot_stacked_bar_chart(proportions, patient_totals, output_prefix):
    """
    绘制堆积柱状图
    """
    # 设置图形大小
    plt.figure(figsize=(max(12, len(proportions)*0.5), 8))
    
    # 准备数据
    patient_names = proportions.index.tolist()
    celltypes = proportions.columns.tolist()
    
    # 按patient的总细胞数排序
    sorted_patients = patient_totals.sort_values(ascending=False).index.tolist()
    proportions = proportions.loc[sorted_patients]
    
    # 设置颜色方案
    n_colors = len(celltypes)
    colors = plt.cm.Set3(np.linspace(0, 1, n_colors))
    
    # 创建底部值
    bottom = np.zeros(len(patient_names))
    
    # 绘制堆积柱状图
    for i, celltype in enumerate(celltypes):
        values = proportions[celltype].values
        
        # 只绘制有值的部分
        if values.sum() > 0:
            plt.bar(patient_names, values, bottom=bottom, 
                   color=colors[i], edgecolor='white', 
                   width=0.8, label=celltype)
            bottom += values
    
    # 设置图表属性
    plt.xlabel('Patient', fontsize=12)
    plt.ylabel('Proportion (%)', fontsize=12)
    plt.title('Cell Type Proportions by Patient (MPC Group)', fontsize=14, fontweight='bold')
    
    # 旋转x轴标签
    plt.xticks(rotation=45, ha='right')
    
    # 添加网格
    plt.grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # 添加图例
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, 
               title='Cell Types', title_fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图形
    output_file = f"{output_prefix}_stacked_bar_chart.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"堆积柱状图已保存到: {output_file}")
    
    # 同时保存为PDF格式
    output_file_pdf = f"{output_prefix}_stacked_bar_chart.pdf"
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"堆积柱状图(PDF)已保存到: {output_file_pdf}")
    
    plt.show()

def create_summary_statistics(adata, proportions, output_prefix):
    """
    创建汇总统计
    """
    print("\n创建汇总统计...")
    
    # 1. 计算总体细胞类型分布
    overall_celltype_counts = adata.obs['annotation'].value_counts()
    overall_total = overall_celltype_counts.sum()
    overall_proportions = (overall_celltype_counts / overall_total * 100).round(2)
    
    overall_stats = pd.DataFrame({
        'celltype': overall_celltype_counts.index,
        'count': overall_celltype_counts.values,
        'percentage': overall_proportions.values
    }).sort_values('percentage', ascending=False)
    
    overall_file = f"{output_prefix}_overall_celltype_stats.csv"
    overall_stats.to_csv(overall_file, index=False, encoding='utf-8-sig')
    print(f"总体细胞类型统计已保存到: {overall_file}")
    
    # 2. 计算每个细胞类型在不同patient中的平均值和标准差
    celltype_means = proportions.mean().round(2)
    celltype_stds = proportions.std().round(2)
    celltype_stats = pd.DataFrame({
        'celltype': celltype_means.index,
        'mean_percentage': celltype_means.values,
        'std_percentage': celltype_stds.values
    }).sort_values('mean_percentage', ascending=False)
    
    summary_file = f"{output_prefix}_celltype_summary_stats.csv"
    celltype_stats.to_csv(summary_file, index=False, encoding='utf-8-sig')
    print(f"细胞类型汇总统计已保存到: {summary_file}")
    
    # 3. 创建热图显示细胞类型分布
    create_heatmap(proportions, output_prefix)
    
    return overall_stats, celltype_stats

def create_heatmap(proportions, output_prefix):
    """
    创建热图显示细胞类型分布
    """
    print("创建热图...")
    
    # 设置图形大小
    plt.figure(figsize=(max(10, len(proportions.columns)*0.5), 
                        max(8, len(proportions.index)*0.4)))
    
    # 创建热图
    sns.heatmap(proportions.T, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Percentage (%)'})
    
    # 设置图表属性
    plt.title('Cell Type Proportions by Patient (Heatmap)', fontsize=14, fontweight='bold')
    plt.xlabel('Patient', fontsize=12)
    plt.ylabel('Cell Type', fontsize=12)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图形
    output_file = f"{output_prefix}_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"热图已保存到: {output_file}")
    
    plt.show()

def main():
    """
    主函数
    """
    # 文件路径
    h5ad_path = "MPCvsTNBC_all_sample_ct.h5ad"
    output_prefix = "MPC_analysis"
    
    try:
        # 执行分析
        mpc_adata, proportions, patient_totals = calculate_and_plot_celltype_proportions(
            h5ad_path, output_prefix
        )
        
        # 打印汇总信息
        print(f"\n=== 分析完成 ===")
        print(f"分析的样本组: MPC")
        print(f"patient数量: {len(proportions)}")
        print(f"细胞类型数量: {len(proportions.columns)}")
        print(f"总细胞数: {mpc_adata.n_obs}")
        
        # 打印前5个最常见的细胞类型
        overall_counts = mpc_adata.obs['annotation'].value_counts()
        print(f"\n前5个最常见的细胞类型:")
        for i, (celltype, count) in enumerate(overall_counts.head(5).items()):
            percentage = (count / mpc_adata.n_obs) * 100
            print(f"  {i+1}. {celltype}: {count} 个细胞 ({percentage:.2f}%)")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
