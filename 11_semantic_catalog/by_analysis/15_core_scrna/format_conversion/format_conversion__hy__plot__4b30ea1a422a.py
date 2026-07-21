import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np

# 设置图形参数
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
plt.rcParams['font.family'] = ['Arial', 'sans-serif']  # 设置字体

# 1. 读取h5ad文件
adata = sc.read_h5ad('result.h5ad')  # 请替换为你的文件路径

# 检查数据中是否包含必要的注释信息
print("AnnData对象基本信息:")
print(adata)
print("\n观察数据（obs）中的列:", list(adata.obs.columns))
print("\n样本列名称:", adata.obs['sample'].unique()[:5])  # 显示前5个样本名称
print("细胞类型列名称:", adata.obs['cell_type'].unique()[:10])  # 显示前10个细胞类型

# 2. 获取所有唯一样本名称
samples = adata.obs['sample'].unique().tolist()
print(f"\n总共有 {len(samples)} 个样本需要处理")

# 3. 创建输出目录
output_dir = "umap_by_sample"
os.makedirs(output_dir, exist_ok=True)

# 4. 为每个样本生成UMAP图
for i, sample in enumerate(samples):
    print(f"正在处理样本 {i+1}/{len(samples)}: {sample}")
    
    # 提取当前样本的数据
    sample_mask = adata.obs['sample'] == sample
    sample_adata = adata[sample_mask].copy()
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # 绘制UMAP图
    sc.pl.umap(
        sample_adata,
        color='cell_type',           # 按细胞类型着色
        title=f'{sample} - Cell Types',  # 标题包含样本名称
        ax=ax,
        show=False,                 # 不立即显示，以便我们后续保存
        legend_loc='right margin',  # 图例位置
        size=20,                    # 点的大小
        frameon=False,              # 不显示边框
        alpha=0.7                   # 点的透明度
    )
    
    # 调整图形外观
    ax.set_xlabel('UMAP1', fontsize=12)
    ax.set_ylabel('UMAP2', fontsize=12)
    ax.title.set_size(14)
    
    # 保存图形
    plt.tight_layout()
    plt.savefig(f'{output_dir}/umap_{sample}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/umap_{sample}.pdf', bbox_inches='tight')
    
    # 关闭图形以释放内存
    plt.close(fig)

print(f"\n所有样本的UMAP图已保存到 {output_dir} 目录")

# 5. (可选) 生成一个汇总图，显示所有样本的缩略图
def create_summary_plot(adata, samples, output_dir):
    """创建包含所有样本缩略图的汇总图"""
    n_samples = len(samples)
    n_cols = 4  # 每行显示的图数量
    n_rows = (n_samples + n_cols - 1) // n_cols  # 计算需要的行数
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3))
    axes = axes.flatten()  # 将轴数组展平，便于迭代
    
    for i, sample in enumerate(samples):
        if i < len(axes):
            # 提取当前样本数据
            sample_mask = adata.obs['sample'] == sample
            sample_adata = adata[sample_mask].copy()
            
            # 绘制到对应的子图中
            sc.pl.umap(
                sample_adata,
                color='cell_type',
                title=sample,
                ax=axes[i],
                show=False,
                legend_loc=None,  # 缩略图中省略图例
                size=10,
                frameon=False
            )
            # 调整子图标题大小
            axes[i].title.set_size(10)
    
    # 隐藏多余的子图
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/all_samples_summary.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

# 执行汇总图创建
create_summary_plot(adata, samples, output_dir)
print("汇总图已生成")

# 6. 输出处理摘要
print("\n=== 处理完成摘要 ===")
print(f"处理的样本总数: {len(samples)}")
print(f"输出目录: {output_dir}")
print("每个样本生成的文件:")
print("  - umap_[样本名].png (高分辨率PNG)")
print("  - umap_[样本名].pdf (矢量PDF)")
print("  - all_samples_summary.png (所有样本汇总图)")
