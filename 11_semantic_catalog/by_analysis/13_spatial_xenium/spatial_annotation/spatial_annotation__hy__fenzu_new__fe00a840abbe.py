import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ========== 顶刊配色方案 ==========
NPG_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#0072B5", "#BC3C29", "#FFDC91", "#EE4C97", "#631879",
    "#5F559B", "#A20056", "#808080", "#1B1919", "#F7F7F7"
]

def get_top_colors(n):
    """返回 n 种顶刊配色，循环使用色板"""
    colors = []
    for i in range(n):
        colors.append(NPG_COLORS[i % len(NPG_COLORS)])
    return colors

def remove_top_right_spines(ax):
    """去掉图表的上边框和右边框，只保留X轴和Y轴"""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # 可选：调整剩余边框的线宽
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

def save_csv(df, path):
    df.to_csv(path, encoding='utf-8-sig')
    print(f"数据已保存: {path}")

# ========== 绘图函数 ==========
def plot_grouped_patient_stacked(patient_props, patient_totals, patient_groups,
                                 output_prefix="MPCvsTNBC_patients"):
    """
    所有患者合并的堆积柱状图，顶部显示分组颜色条，图例在下方，无上右边框。
    """
    # 按组排序患者：先MPC后TNBC，组内按细胞总数降序
    mpc_order = patient_totals[patient_groups == 'MPC'].sort_values(ascending=False).index.tolist()
    tnbc_order = patient_totals[patient_groups == 'TNBC'].sort_values(ascending=False).index.tolist()
    patient_order = mpc_order + tnbc_order
    patient_props = patient_props.loc[patient_order]

    celltypes = patient_props.columns.tolist()
    n_patients = len(patient_order)
    n_celltypes = len(celltypes)

    # 动态尺寸（多留顶部空间给颜色条，底部空间给图例）
    fig_width = max(14, n_patients * 0.6)
    fig, ax = plt.subplots(figsize=(fig_width, 8.5))

    colors = get_top_colors(n_celltypes)

    # ---------- 绘制堆积柱状图 ----------
    bottom = np.zeros(n_patients)
    bar_width = 0.7
    x_positions = np.arange(n_patients)
    for i, ct in enumerate(celltypes):
        vals = patient_props[ct].values
        ax.bar(x_positions, vals, bottom=bottom,
               color=colors[i], edgecolor='white', width=bar_width, label=ct)
        bottom += vals

    # ---------- X轴设置（只显示患者ID，无标题）----------
    ax.set_xticks(x_positions)
    ax.set_xticklabels(patient_order, rotation=45, ha='right', fontsize=9)
    ax.set_xlabel('')
    ax.set_xlim(-0.5, n_patients - 0.5)

    # ---------- 在柱子顶部添加分组颜色条 ----------
    y_max = 100
    bar_height = y_max * 0.03                     # 颜色条高度
    y_bottom = y_max + y_max * 0.01              # 在柱子顶部稍上位置
    group_colors = ['#E64B35' if patient_groups[p] == 'MPC' else '#4DBBD5' for p in patient_order]
    ax.bar(x_positions, bar_height, bottom=y_bottom, color=group_colors,
           edgecolor='none', width=bar_width)
    # 调整Y轴上限，留出颜色条空间
    ax.set_ylim(0, y_max * 1.10)

    # 在颜色条上方添加分组文字标签
    mpc_x_center = (len(mpc_order) - 1) / 2 if mpc_order else 0
    tnbc_x_center = len(mpc_order) + (len(tnbc_order) - 1) / 2 if tnbc_order else 0
    text_y = y_bottom + bar_height * 2.0
    if mpc_order:
        ax.text(mpc_x_center, text_y, 'MPC', ha='center', fontsize=12,
                color='#E64B35', fontweight='bold')
    if tnbc_order:
        ax.text(tnbc_x_center, text_y, 'TNBC', ha='center', fontsize=12,
                color='#4DBBD5', fontweight='bold')

    # ---------- 图例移至下方 ----------
    # 先调整图表区域，为底部图例留空间
    plt.subplots_adjust(bottom=0.25)  # 调大bottom值
    ncol = min(6, n_celltypes)  # 每行最多6个
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              ncol=ncol, fontsize=8, title='Cell Types', title_fontsize=10,
              frameon=False)

    # ---------- 样式调整 ----------
    ax.set_ylabel('Proportion (%)', fontsize=12)
    ax.set_title('Cell Type Composition by Patient (MPC vs TNBC)', fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

    # 去掉上右边框
    remove_top_right_spines(ax)

    # 保存
    png_file = f"{output_prefix}_all_patients_stacked.png"
    pdf_file = f"{output_prefix}_all_patients_stacked.pdf"
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.show()
    print(f"合并患者堆积图已保存: {png_file}, {pdf_file}")

def plot_group_stacked_bar(proportions, group_totals, output_prefix="MPCvsTNBC"):
    """分组对比堆积柱状图，无上右边框，图例在下方"""
    groups = proportions.index.tolist()
    celltypes = proportions.columns.tolist()
    n_celltypes = len(celltypes)

    fig_width = max(8, 2.5 * len(groups))
    fig, ax = plt.subplots(figsize=(fig_width, 7.0))

    colors = get_top_colors(n_celltypes)
    bottom = np.zeros(len(groups))
    bar_width = 0.7

    for i, ct in enumerate(celltypes):
        values = proportions[ct].values
        ax.bar(groups, values, bottom=bottom,
               color=colors[i], edgecolor='white',
               width=bar_width, label=ct)
        bottom += values

    ax.set_xlabel('')
    ax.set_ylabel('Proportion (%)', fontsize=12)
    ax.set_title('Cell Type Composition: MPC vs TNBC', fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

    # 图例在下方
    plt.subplots_adjust(bottom=0.3)
    ncol = min(6, n_celltypes)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18),
              ncol=ncol, fontsize=8, title='Cell Types', title_fontsize=10,
              frameon=False)

    # 去边框
    remove_top_right_spines(ax)

    png_file = f"{output_prefix}_stacked_bar_chart.png"
    pdf_file = f"{output_prefix}_stacked_bar_chart.pdf"
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.show()
    print(f"分组柱状图已保存: {png_file}, {pdf_file}")

def plot_patient_stacked_bar(proportions, patient_totals, group_name, output_prefix):
    """各分组内部患者堆积柱状图，无上右边框，图例在下方"""
    sorted_patients = patient_totals.sort_values(ascending=False).index.tolist()
    proportions = proportions.loc[sorted_patients]

    celltypes = proportions.columns.tolist()
    n_patients = len(sorted_patients)
    n_celltypes = len(celltypes)

    fig_width = max(10, n_patients * 0.7)
    fig, ax = plt.subplots(figsize=(fig_width, 7.5))

    colors = get_top_colors(n_celltypes)
    bottom = np.zeros(n_patients)

    for i, ct in enumerate(celltypes):
        values = proportions[ct].values
        ax.bar(sorted_patients, values, bottom=bottom,
               color=colors[i], edgecolor='white',
               width=0.75, label=ct)
        bottom += values

    ax.set_xlabel('')
    ax.set_ylabel('Proportion (%)', fontsize=12)
    ax.set_title(f'Cell Type Proportions by Patient ({group_name} Group)', fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(range(n_patients))
    ax.set_xticklabels(sorted_patients, rotation=45, ha='right', fontsize=9)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)

    # 图例在下方
    plt.subplots_adjust(bottom=0.3)
    ncol = min(6, n_celltypes)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.18),
              ncol=ncol, fontsize=7, title='Cell Types', title_fontsize=9,
              frameon=False)

    # 去边框
    remove_top_right_spines(ax)

    png_file = f"{output_prefix}_patient_stacked_{group_name}.png"
    pdf_file = f"{output_prefix}_patient_stacked_{group_name}.pdf"
    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.show()
    print(f"{group_name}组患者柱状图已保存: {png_file}, {pdf_file}")

# ========== 主流程 ==========
def main():
    h5ad_path = "MPCvsTNBC_all_sample_ct.h5ad"  # 请根据实际情况修改
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

    required_cols = ['group', 'patient', 'annotation']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"缺少必须列: '{col}'")

    # 筛选两组数据
    groups = ['MPC', 'TNBC']
    filtered = adata[adata.obs['group'].isin(groups), :].copy()
    if filtered.n_obs == 0:
        raise ValueError("没有找到MPC或TNBC的细胞！")

    # ---------- 1. 分组对比图 ----------
    print("\n======= 1. MPC vs TNBC 分组对比图 =======")
    group_counts = pd.crosstab(filtered.obs['group'], filtered.obs['annotation'])
    group_totals = group_counts.sum(axis=1)
    group_props = (group_counts.div(group_totals, axis=0) * 100).round(2)
    save_csv(group_props, "MPCvsTNBC_group_celltype_proportions.csv")
    save_csv(group_totals.to_frame('total_cells'), "MPCvsTNBC_group_totals.csv")
    plot_group_stacked_bar(group_props, group_totals, output_prefix="MPCvsTNBC")

    # ---------- 2. 合并所有患者图（画在一起，顶部颜色条）----------
    print("\n======= 2. 合并所有患者堆积图 =======")
    all_patient_counts = pd.crosstab(filtered.obs['patient'], filtered.obs['annotation'])
    all_patient_totals = all_patient_counts.sum(axis=1)
    all_patient_props = (all_patient_counts.div(all_patient_totals, axis=0) * 100).round(2)

    # 患者分组信息
    patient_groups = filtered.obs[['patient', 'group']].drop_duplicates().set_index('patient')['group']

    save_csv(all_patient_props, "MPCvsTNBC_all_patients_proportions.csv")
    save_csv(all_patient_totals.to_frame('total_cells'), "MPCvsTNBC_all_patients_totals.csv")

    plot_grouped_patient_stacked(all_patient_props, all_patient_totals, patient_groups,
                                 output_prefix="MPCvsTNBC")

    # ---------- 3. 各组内部患者图 ----------
    for group in groups:
        print(f"\n======= 3. {group}组内部患者图 =======")
        sub = filtered[filtered.obs['group'] == group, :]
        if sub.n_obs == 0:
            continue
        pat_counts = pd.crosstab(sub.obs['patient'], sub.obs['annotation'])
        pat_totals = pat_counts.sum(axis=1)
        pat_props = (pat_counts.div(pat_totals, axis=0) * 100).round(2)
        prefix = f"{group}_patient_analysis"
        save_csv(pat_props, f"{prefix}_celltype_proportions.csv")
        save_csv(pat_totals.to_frame('total_cells'), f"{prefix}_patient_totals.csv")
        plot_patient_stacked_bar(pat_props, pat_totals, group, output_prefix=prefix)

    print("\n======= 全部分析完成！ =======")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
