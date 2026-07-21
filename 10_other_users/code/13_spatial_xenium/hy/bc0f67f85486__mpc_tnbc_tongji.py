import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

def get_colors(n):
    """根据类别数量返回合适的颜色列表"""
    if n <= 10:
        cmap = plt.cm.tab10
    elif n <= 20:
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.gist_rainbow  # 多种颜色，但可能对比度下降
    return cmap(np.linspace(0, 1, n))

def save_proportions_csv(df, filename):
    """保存占比数据到 CSV"""
    df.to_csv(filename, encoding='utf-8-sig')
    print(f"数据已保存到: {filename}")

def plot_group_stacked_bar(proportions, group_totals, output_prefix="MPCvsTNBC"):
    """
    绘制分组堆积柱状图（MPC vs TNBC），柱子加宽，图例优化
    """
    groups = proportions.index.tolist()
    celltypes = proportions.columns.tolist()
    n_celltypes = len(celltypes)

    # 图形尺寸：宽度根据分组数微调，高度固定
    fig_width = max(8, 2.5 * len(groups))   # 柱子少时也至少有8英寸宽
    plt.figure(figsize=(fig_width, 6.5))

    colors = get_colors(n_celltypes)
    bottom = np.zeros(len(groups))
    bar_width = 0.7  # 柱子加宽（之前是0.6）

    for i, ct in enumerate(celltypes):
        values = proportions[ct].values
        if values.sum() > 0:
            plt.bar(groups, values, bottom=bottom,
                    color=colors[i], edgecolor='white',
                    width=bar_width, label=ct)
            bottom += values

    plt.xlabel('Group', fontsize=12)
    plt.ylabel('Proportion (%)', fontsize=12)
    plt.title('Cell Type Composition: MPC vs TNBC', fontsize=14, fontweight='bold')
    plt.grid(True, axis='y', alpha=0.3, linestyle='--')

    # 图例处理：如果类别很多，放在图下方并分多列
    if n_celltypes > 15:
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                   ncol=6, fontsize=7, title='Cell Types', title_fontsize=9,
                   frameon=False)
        plt.subplots_adjust(bottom=0.25)  # 为图例留空间
    else:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                   fontsize=8, title='Cell Types', title_fontsize=10)

    plt.tight_layout()

    # 保存
    png_file = f"{output_prefix}_stacked_bar_chart.png"
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    pdf_file = f"{output_prefix}_stacked_bar_chart.pdf"
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.show()
    print(f"分组柱状图已保存: {png_file}, {pdf_file}")

def plot_patient_stacked_bar(proportions, patient_totals, group_name, output_prefix):
    """
    绘制某个分组内各患者的堆积柱状图，宽度随患者数动态变化
    """
    # 按总细胞数降序排列患者
    sorted_patients = patient_totals.sort_values(ascending=False).index.tolist()
    proportions = proportions.loc[sorted_patients]  # 按顺序重排

    celltypes = proportions.columns.tolist()
    n_patients = len(sorted_patients)
    n_celltypes = len(celltypes)

    # 动态宽度：每个患者至少0.6英寸，最少10英寸
    fig_width = max(10, n_patients * 0.7)
    plt.figure(figsize=(fig_width, 7))

    colors = get_colors(n_celltypes)
    bottom = np.zeros(n_patients)

    for i, ct in enumerate(celltypes):
        values = proportions[ct].values
        if values.sum() > 0:
            plt.bar(sorted_patients, values, bottom=bottom,
                    color=colors[i], edgecolor='white',
                    width=0.75, label=ct)
            bottom += values

    plt.xlabel('Patient', fontsize=12)
    plt.ylabel('Proportion (%)', fontsize=12)
    plt.title(f'Cell Type Proportions by Patient ({group_name} Group)', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, axis='y', alpha=0.3, linestyle='--')

    # 图例处理
    if n_celltypes > 15:
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
                   ncol=6, fontsize=7, title='Cell Types', title_fontsize=9,
                   frameon=False)
        plt.subplots_adjust(bottom=0.3)
    else:
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                   fontsize=7, title='Cell Types', title_fontsize=9)

    plt.tight_layout()

    # 保存
    png_file = f"{output_prefix}_patient_stacked_bar_{group_name}.png"
    pdf_file = f"{output_prefix}_patient_stacked_bar_{group_name}.pdf"
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.show()
    print(f"{group_name}组患者柱状图已保存: {png_file}, {pdf_file}")

def calculate_and_plot_all(h5ad_path):
    """主分析流程：分组对比 + 各分组患者分析"""
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

    required_cols = ['group', 'patient', 'annotation']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"缺少必须列: '{col}'")

    # 1. 分组对比分析
    print("\n======= 1. MPC vs TNBC 分组分析 =======")
    groups = ['MPC', 'TNBC']
    filtered = adata[adata.obs['group'].isin(groups), :].copy()
    if filtered.n_obs == 0:
        raise ValueError("没有找到MPC或TNBC的细胞！")

    # 分组占比
    group_counts = pd.crosstab(filtered.obs['group'], filtered.obs['annotation'])
    group_totals = group_counts.sum(axis=1)
    group_props = (group_counts.div(group_totals, axis=0) * 100).round(2)

    save_proportions_csv(group_props, "MPCvsTNBC_group_celltype_proportions.csv")
    save_proportions_csv(group_totals.to_frame('total_cells'), "MPCvsTNBC_group_totals.csv")
    plot_group_stacked_bar(group_props, group_totals, output_prefix="MPCvsTNBC")

    # 2. 各分组患者分析
    for group in groups:
        print(f"\n======= 2. {group}组患者分析 =======")
        sub = filtered[filtered.obs['group'] == group, :].copy()
        if sub.n_obs == 0:
            print(f"{group}组无细胞，跳过")
            continue

        # 患者×细胞类型计数
        patient_counts = pd.crosstab(sub.obs['patient'], sub.obs['annotation'])
        patient_totals = patient_counts.sum(axis=1)
        patient_props = (patient_counts.div(patient_totals, axis=0) * 100).round(2)

        # 保存数据
        prefix = f"{group}_patient_analysis"
        save_proportions_csv(patient_props, f"{prefix}_celltype_proportions.csv")
        save_proportions_csv(patient_totals.to_frame('total_cells'), f"{prefix}_patient_totals.csv")

        # 绘制患者堆积柱状图
        plot_patient_stacked_bar(patient_props, patient_totals, group, output_prefix=prefix)

    print("\n======= 全部分析完成 =======")

if __name__ == "__main__":
    h5ad_path = "MPCvsTNBC_all_sample_ct.h5ad"  # 根据实际情况修改
    try:
        calculate_and_plot_all(h5ad_path)
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
