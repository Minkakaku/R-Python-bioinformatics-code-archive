import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.colors import ListedColormap
import pandas as pd

# 设置matplotlib以支持PDF输出
plt.rcParams['pdf.fonttype'] = 42  # 确保文本在PDF中可编辑
plt.rcParams['ps.fonttype'] = 42

# 1. 读取h5ad文件
adata = sc.read_h5ad("merged_TMA_data.h5ad")

# 2. 检查所需列是否存在
required_cols = ['mucus_zone', 'sample', 'CNs_10']
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"列 '{col}' 在adata.obs中不存在")

# 3. 创建输出文件夹
output_dir = "CN_2_7_8"
os.makedirs(output_dir, exist_ok=True)
print(f"创建输出文件夹: {output_dir}")

# 4. 获取所有样本
samples = adata.obs['sample'].unique()
print(f"找到 {len(samples)} 个样本: {samples}")

# 5. 定义配色方案
# 创建一个新的颜色列
adata.obs['plot_color'] = None

# 首先，为CN2, CN7, CN8指定颜色
cn2_mask = adata.obs['CNs_10'] == 'CN2'
cn7_mask = adata.obs['CNs_10'] == 'CN7'
cn8_mask = adata.obs['CNs_10'] == 'CN8'

# 应用CN颜色
adata.obs.loc[cn2_mask, 'plot_color'] = 'red'
adata.obs.loc[cn7_mask, 'plot_color'] = 'blue'
adata.obs.loc[cn8_mask, 'plot_color'] = 'yellow'

# 然后，为其余细胞应用mucus_zone颜色
remaining_mask = ~cn2_mask & ~cn7_mask & ~cn8_mask
inner_mask = (adata.obs['mucus_zone'] == 'inner') & remaining_mask
outer_mask = (adata.obs['mucus_zone'] == 'outer') & remaining_mask

# 深灰色和浅灰色
adata.obs.loc[inner_mask, 'plot_color'] = '#333333'  # 深灰色
adata.obs.loc[outer_mask, 'plot_color'] = '#CCCCCC'  # 浅灰色

# 6. 检查是否有未着色的细胞
uncolored = adata.obs['plot_color'].isna().sum()
if uncolored > 0:
    print(f"警告: 有 {uncolored} 个细胞未着色，将使用默认颜色")
    adata.obs.loc[adata.obs['plot_color'].isna(), 'plot_color'] = '#888888'  # 中性灰色

# 7. 检查空间坐标
if 'spatial' in adata.obsm:
    spatial_key = 'spatial'
elif 'X_spatial' in adata.obsm:
    spatial_key = 'X_spatial'
elif hasattr(adata, 'uns') and 'spatial' in adata.uns:
    print("空间坐标在adata.uns['spatial']中")
    if 'spatial' in adata.obsm:
        spatial_key = 'spatial'
    else:
        raise ValueError("无法从uns中提取空间坐标")
elif 'X_umap' in adata.obsm:
    print("使用UMAP坐标进行可视化")
    spatial_key = 'X_umap'
elif 'X_tsne' in adata.obsm:
    print("使用t-SNE坐标进行可视化")
    spatial_key = 'X_tsne'
else:
    raise ValueError("未找到空间坐标。请确保数据包含空间坐标。")

# 8. 为每个样本绘制空间图
for sample in samples:
    print(f"处理样本: {sample}")
    
    # 提取当前样本的数据
    sample_mask = adata.obs['sample'] == sample
    sample_adata = adata[sample_mask].copy()
    
    if sample_adata.n_obs == 0:
        print(f"  警告: 样本 {sample} 没有细胞，跳过")
        continue
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # 子图1: 按mucus_zone和CNs_10着色
    ax1 = axes[0]
    
    # 获取当前样本的颜色
    colors = sample_adata.obs['plot_color']
    
    # 绘制散点图
    if spatial_key in sample_adata.obsm:
        coords = sample_adata.obsm[spatial_key]
        if coords.shape[1] >= 2:
            x = coords[:, 0]
            y = coords[:, 1]
            
            # 绘制所有点
            scatter = ax1.scatter(x, y, c=colors, s=10, edgecolors='none')
            
            # 设置坐标轴 (使用英文标签)
            ax1.set_aspect('equal')
            ax1.invert_yaxis()  # 通常空间数据的y轴是反的
            ax1.set_xlabel("X coordinate")
            ax1.set_ylabel("Y coordinate")
            ax1.set_title(f"Sample {sample}\nColored by mucus_zone and CN2/7/8")
        else:
            print(f"  警告: 样本 {sample} 的空间坐标维度不足")
            ax1.text(0.5, 0.5, f"Insufficient spatial coordinates\n{coords.shape[1]}D", 
                     ha='center', va='center', transform=ax1.transAxes)
    else:
        print(f"  警告: 样本 {sample} 没有空间坐标")
        ax1.text(0.5, 0.5, "No spatial coordinates", ha='center', va='center', transform=ax1.transAxes)
    
    # 创建图例
    from matplotlib.patches import Patch
    
    # 定义图例元素 (使用英文标签)
    legend_elements = [
        Patch(facecolor='red', edgecolor='black', label='CN2'),
        Patch(facecolor='blue', edgecolor='black', label='CN7'),
        Patch(facecolor='yellow', edgecolor='black', label='CN8'),
        Patch(facecolor='#333333', edgecolor='black', label='mucus_zone: inner'),
        Patch(facecolor='#CCCCCC', edgecolor='black', label='mucus_zone: outer'),
    ]
    
    # 如果有未着色的细胞，添加到图例
    if (sample_adata.obs['plot_color'] == '#888888').any():
        legend_elements.append(Patch(facecolor='#888888', edgecolor='black', label='Other'))
    
    ax1.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    
    # 子图2: 细胞类型/CNs_10的统计信息
    ax2 = axes[1]
    
    # 统计当前样本中各种类型的细胞数量
    cn_counts = {}
    for cn in ['CN2', 'CN7', 'CN8']:
        cn_counts[cn] = (sample_adata.obs['CNs_10'] == cn).sum()
    
    # 统计mucus_zone
    inner_count = (sample_adata.obs['mucus_zone'] == 'inner').sum()
    outer_count = (sample_adata.obs['mucus_zone'] == 'outer').sum()
    
    # 创建统计文本 (使用英文)
    stats_text = f"Sample: {sample}\n"
    stats_text += f"Total cells: {sample_adata.n_obs}\n\n"
    
    stats_text += "CNs_10 distribution:\n"
    for cn, count in cn_counts.items():
        if count > 0:
            percentage = count / sample_adata.n_obs * 100
            stats_text += f"  {cn}: {count} ({percentage:.1f}%)\n"
    
    stats_text += "\nmucus_zone distribution:\n"
    stats_text += f"  inner: {inner_count} ({inner_count/sample_adata.n_obs*100:.1f}%)\n"
    stats_text += f"  outer: {outer_count} ({outer_count/sample_adata.n_obs*100:.1f}%)\n"
    
    # 其他类型的CNs_10
    other_cns = sample_adata.obs['CNs_10'][~sample_adata.obs['CNs_10'].isin(['CN2', 'CN7', 'CN8'])].unique()
    if len(other_cns) > 0:
        stats_text += f"\nOther CNs_10 types: {', '.join(other_cns)}"
    
    ax2.text(0.05, 0.95, stats_text, transform=ax2.transAxes, 
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis('off')
    ax2.set_title("Sample Statistics")  # 使用英文标题
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图形为PDF和PNG两种格式
    safe_sample_name = str(sample).replace('/', '_').replace('\\', '_').replace(':', '_')
    
    # 保存为PDF（矢量格式，适合出版）
    pdf_path = os.path.join(output_dir, f"sample_{safe_sample_name}.pdf")
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', format='pdf')
    print(f"  已保存PDF: {pdf_path}")
    
    # 保存为PNG（位图格式，适合网页展示）
    png_path = os.path.join(output_dir, f"sample_{safe_sample_name}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight', format='png')
    print(f"  已保存PNG: {png_path}")
    
    plt.close()
    
# 9. 创建汇总统计
print("\n创建汇总统计...")
summary_data = []

for sample in samples:
    sample_mask = adata.obs['sample'] == sample
    if sample_mask.sum() == 0:
        continue
    
    sample_adata = adata[sample_mask]
    
    # 计算统计
    total_cells = sample_adata.n_obs
    cn2_count = (sample_adata.obs['CNs_10'] == 'CN2').sum()
    cn7_count = (sample_adata.obs['CNs_10'] == 'CN7').sum()
    cn8_count = (sample_adata.obs['CNs_10'] == 'CN8').sum()
    inner_count = (sample_adata.obs['mucus_zone'] == 'inner').sum()
    outer_count = (sample_adata.obs['mucus_zone'] == 'outer').sum()
    
    summary_data.append({
        'sample': sample,
        'total_cells': total_cells,
        'CN2_count': cn2_count,
        'CN2_percent': cn2_count/total_cells*100 if total_cells>0 else 0,
        'CN7_count': cn7_count,
        'CN7_percent': cn7_count/total_cells*100 if total_cells>0 else 0,
        'CN8_count': cn8_count,
        'CN8_percent': cn8_count/total_cells*100 if total_cells>0 else 0,
        'inner_count': inner_count,
        'inner_percent': inner_count/total_cells*100 if total_cells>0 else 0,
        'outer_count': outer_count,
        'outer_percent': outer_count/total_cells*100 if total_cells>0 else 0,
    })

# 保存汇总统计为CSV
summary_df = pd.DataFrame(summary_data)
summary_path = os.path.join(output_dir, "sample_summary.csv")
summary_df.to_csv(summary_path, index=False)
print(f"汇总统计已保存: {summary_path}")

# 10. 创建所有样本的概览图
print("创建所有样本概览图...")
n_samples = len(samples)
n_cols = min(4, n_samples)
n_rows = (n_samples + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))

if n_rows == 1 and n_cols == 1:
    axes = np.array([axes])
elif n_rows == 1 or n_cols == 1:
    axes = axes.reshape(-1)

for idx, (sample, ax) in enumerate(zip(samples, axes.flatten())):
    sample_mask = adata.obs['sample'] == sample
    sample_adata = adata[sample_mask]
    
    if sample_adata.n_obs == 0 or spatial_key not in sample_adata.obsm:
        ax.text(0.5, 0.5, f"Sample {sample}\nNo data", 
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        continue
    
    coords = sample_adata.obsm[spatial_key]
    if coords.shape[1] >= 2:
        x = coords[:, 0]
        y = coords[:, 1]
        colors = sample_adata.obs['plot_color']
        
        ax.scatter(x, y, c=colors, s=5, edgecolors='none')
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_title(f"Sample {sample}\n({sample_adata.n_obs} cells)")  # 使用英文
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis='both', which='both', length=0)
    else:
        ax.text(0.5, 0.5, f"Sample {sample}\nNo spatial coordinates", 
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')

# 隐藏多余的子图
for idx in range(len(samples), len(axes.flatten())):
    axes.flatten()[idx].axis('off')

plt.suptitle("Spatial Distribution of All Samples", fontsize=16, y=1.02)  # 使用英文标题
plt.tight_layout()

# 保存概览图为PDF和PNG
overview_pdf_path = os.path.join(output_dir, "all_samples_overview.pdf")
plt.savefig(overview_pdf_path, dpi=300, bbox_inches='tight', format='pdf')
print(f"概览图PDF已保存: {overview_pdf_path}")

overview_png_path = os.path.join(output_dir, "all_samples_overview.png")
plt.savefig(overview_png_path, dpi=300, bbox_inches='tight', format='png')
print(f"概览图PNG已保存: {overview_png_path}")

plt.close()

print(f"\n完成! 所有结果已保存到 '{output_dir}' 文件夹中。")
print(f"共处理了 {len(samples)} 个样本。")
print(f"每个样本生成了PDF和PNG两种格式的图像。")
