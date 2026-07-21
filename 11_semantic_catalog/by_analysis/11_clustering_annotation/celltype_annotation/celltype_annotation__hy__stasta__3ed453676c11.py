# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools
from matplotlib import font_manager
from matplotlib.patches import Patch
from statannotations.Annotator import Annotator

# ==================== 0. 字体与样式设置 ====================
# 如果字体路径不对，请注释掉下面这行，改用系统默认
font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')

plt.rcParams['font.family'] = 'arial'
plt.rcParams['font.weight'] = 'heavy'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ==================== 1. 参数调整区 ====================
Y_AXIS_PADDING_FACTOR = 1.5  # ������ Y轴顶部预留空间倍数，不够可以改成 1.6 或 1.7
LINE_HEIGHT = 0.025          # ������ 显著性横线之间的间距 (0.025 的比例)

# ==================== 2. 读取数据 ====================
file_path = "All_filter.h5ad"
print(f"正在读取文件: {file_path}...")
adata = sc.read_h5ad(file_path)

# ==================== 3. 数据处理：计算细胞比例 ====================
required_cols = ['sample', 'celltype', 'group_final']
if not all(col in adata.obs.columns for col in required_cols):
    raise ValueError(f"adata.obs 中必须包含以下列: {required_cols}")

sample_to_group = adata.obs.groupby('sample', observed=False)['group_final'].first()
counts = adata.obs.groupby(['sample', 'celltype'], observed=False).size().unstack(fill_value=0)
totals = counts.sum(axis=1)
proportions = counts.div(totals, axis=0) * 100
proportions['group_final'] = sample_to_group

df_plot = proportions.reset_index().melt(
    id_vars=['sample', 'group_final'],
    var_name='celltype',
    value_name='proportion'
)

# ==================== 4. 准备绘图参数 ====================
celltypes = df_plot['celltype'].unique()

# 指定分组顺序 
order = ['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R']
# 注意：请根据您真实数据里的 group 名称修改上面这行。如果名称没变，可直接保留。
# 如果不想固定顺序，改成: order = sorted(df_plot['group_final'].unique())

pairs = list(itertools.combinations(order, 2))

output_dir = "stasta"
os.makedirs(output_dir, exist_ok=True)
print(f"图片将保存至文件夹: {output_dir}")

palette = sns.color_palette("Set2", n_colors=len(order))

# ==================== 5. 循环每个细胞类型单独出图 ====================
for ct in celltypes:
    subset = df_plot[df_plot['celltype'] == ct]
    
    fig, ax = plt.subplots(figsize=(10, 7)) 
    
    # 绘制箱线图
    sns.boxplot(
        x='group_final', y='proportion', data=subset,
        order=order, hue='group_final', palette=palette,
        width=0.5, fliersize=0, legend=False, ax=ax
    )
    
    # 绘制散点抖动图
    sns.stripplot(
        x='group_final', y='proportion', data=subset,
        order=order, color='black', size=4,
        alpha=0.5, jitter=True, ax=ax
    )

    # ==================== 5.1 添加显著性检验 ====================
    try:
        annotator = Annotator(
            ax=ax, pairs=pairs, data=subset,
            x='group_final', y='proportion', order=order
        )
        annotator.configure(
            test='Mann-Whitney',
            text_format='star',
            loc='inside',               
            hide_non_significant=True,  
            line_height=LINE_HEIGHT,    # ������ 应用自定义的线间距
            text_offset=2,              
            line_offset_to_group=0.05   
        )
        annotator.apply_and_annotate()
    except Exception as e:
        print(f"注意: 细胞类型 {ct} 的显著性检验因数据原因跳过: {e}")
        
    # ==================== 5.2 调整 Y 轴范围 ====================
    max_val = subset['proportion'].max()
    if max_val > 0:
        # ������ 大幅提高顶部余量至 1.5 倍 (确保线条能全部塞进图里)
        ax.set_ylim(0, max_val * Y_AXIS_PADDING_FACTOR)  
    else:
        ax.set_ylim(0, 5)

    # ==================== 5.3 设置标签和标题 ====================
    ax.set_title(ct, fontsize=18, fontweight='bold') 
    ax.set_ylabel('Cell Proportion (%)', fontsize=16) 
    ax.set_xlabel('') 
    
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    legend_elements = [Patch(facecolor=palette[i], label=group) for i, group in enumerate(order)]
    ax.legend(handles=legend_elements, loc='lower center', 
              bbox_to_anchor=(0.5, -0.15), ncol=len(order), frameon=False, fontsize=14)

    # ==================== 6. 保存当前图的 PDF ====================
    safe_ct_name = ct.replace("/", "_").replace(" ", "_")
    output_file = os.path.join(output_dir, f"{safe_ct_name}.pdf")
    
    plt.savefig(output_file, format='pdf', bbox_inches='tight')
    plt.close() 
    print(f"已保存: {output_file}")

print("\n全部单图生成完毕！")
