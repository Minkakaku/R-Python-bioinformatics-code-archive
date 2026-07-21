import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from matplotlib import rcParams
from scipy.sparse import issparse
import warnings
from matplotlib.colors import LinearSegmentedColormap
warnings.filterwarnings('ignore')

# 创建输出目录
output_dir = "spa_gene"
os.makedirs(output_dir, exist_ok=True)

# 设置绘图风格
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 1. 读取h5ad文件
print("正在读取数据...")
adata = sc.read_h5ad("merged_TMA_data.h5ad")
print(f"数据形状: {adata.shape}")
print(f"样本信息列名: {adata.obs.columns.tolist()}")
print(f"空间坐标信息: {adata.obsm.keys()}")

# 2. 提取raw矩阵作为新数据集
if adata.raw is not None:
    print("使用adata.raw作为新数据集...")
    raw_adata = adata.raw.to_adata()
    raw_adata.obsm = adata.obsm.copy()
    raw_adata.obs = adata.obs.copy()
else:
    print("未找到raw矩阵，使用原始数据...")
    raw_adata = adata.copy()
    sc.pp.normalize_total(raw_adata, target_sum=1e4)
    sc.pp.log1p(raw_adata)

# 3. 基因列表
gene_list = [
    "VEGFA", "TGFB2", "SNAI2", "HAS2", "BNIP3", 
    "NTRK2", "CXCR4", "THY1", "FSCN1", "SOX8", 
    "LGMN", "CYBA", "ITGA6","MUC5AC","MUC19","FOXA2","ATOH1","IL4","IL13","STAT6","TGFB1","CD44","GATA3",'FGFR2','SOX9','MCAM',
    "JAG1","PDGFB","PDGFRB","CXCL12","CD46","PLAU","PLAUR","THY1","ITGAX","ITGB2","ADM","CALCRL","LGALS9","COL2A1","SDC4","CDH1",
    "TFPI","JAM2","F11R","EGFR","CD44","COL2A1","NECTIN2","TIGIT","CYBA","COL4A1","PDGFA","JAG1","CD46","PLAU","PLAUR","LGALS9"
]

print("\n检查基因在数据集中的存在情况...")
# 检查哪些基因在数据集中存在
available_genes = []
for gene in gene_list:
    # 尝试多种匹配方式
    if gene in raw_adata.var_names:
        available_genes.append(gene)
        print(f"✓ 找到 '{gene}'")
    else:
        # 尝试大小写不敏感匹配
        var_names_upper = [g.upper() for g in raw_adata.var_names]
        if gene.upper() in var_names_upper:
            idx = var_names_upper.index(gene.upper())
            matched_gene = raw_adata.var_names[idx]
            available_genes.append(matched_gene)
            print(f"✓ 找到不区分大小写的匹配 '{gene}' -> '{matched_gene}'")
        else:
            print(f"✗ 基因 '{gene}' 未找到")

if not available_genes:
    print("\n错误: 未找到任何目标基因!")
    print("前20个可用基因:")
    print(raw_adata.var_names[:20].tolist())
    exit()

print(f"\n成功找到 {len(available_genes)}/{len(gene_list)} 个基因")

# 4. 确定样本列名
possible_sample_cols = ['sample', 'Sample', 'patient', 'Patient', 'batch', 'Batch', 'orig.ident', 'slice', 'Slide']
sample_col = None

for col in possible_sample_cols:
    if col in raw_adata.obs.columns:
        sample_col = col
        print(f"使用 '{col}' 作为样本列")
        break

if sample_col is None:
    print("未找到样本列，将所有数据视为一个样本")
    raw_adata.obs['sample'] = 'All_Samples'
    sample_col = 'sample'

# 5. 确定空间坐标
spatial_keys = ['spatial', 'X_spatial', 'coordinates', 'obsm/spatial', 'X_umap']
spatial_key = None

for key in spatial_keys:
    if key in raw_adata.obsm:
        spatial_key = key
        print(f"使用 '{key}' 作为空间坐标")
        break

if spatial_key is None:
    # 尝试查找任何2D坐标
    for key in raw_adata.obsm.keys():
        if raw_adata.obsm[key].shape[1] >= 2:
            spatial_key = key
            print(f"使用 '{key}' 作为空间坐标（推断）")
            break

if spatial_key is None:
    print("警告: 未找到空间坐标信息")
    print("尝试使用PCA或UMAP坐标...")
    if 'X_pca' in raw_adata.obsm and raw_adata.obsm['X_pca'].shape[1] >= 2:
        spatial_key = 'X_pca'
        print(f"使用 '{spatial_key}' 作为坐标")
    elif 'X_umap' in raw_adata.obsm and raw_adata.obsm['X_umap'].shape[1] >= 2:
        spatial_key = 'X_umap'
        print(f"使用 '{spatial_key}' 作为坐标")
    else:
        print("错误: 无法找到可用的坐标信息")
        exit()

# 6. 获取所有样本
samples = sorted(raw_adata.obs[sample_col].unique())
print(f"\n找到 {len(samples)} 个样本: {samples}")

def get_gene_expression(adata_subset, gene):
    """获取基因表达数据的辅助函数"""
    if gene not in adata_subset.var_names:
        return None
    
    # 尝试从layers中获取
    if hasattr(adata_subset, 'layers'):
        if 'counts' in adata_subset.layers:
            expr_data = adata_subset[:, gene].layers['counts']
        elif 'spliced' in adata_subset.layers:
            expr_data = adata_subset[:, gene].layers['spliced']
        else:
            expr_data = adata_subset[:, gene].X
    else:
        expr_data = adata_subset[:, gene].X
    
    # 将稀疏矩阵转换为密集数组
    if issparse(expr_data):
        expr = expr_data.toarray().flatten()
    else:
        expr = expr_data.flatten()
    
    return expr

# 创建自定义颜色映射：从灰白色到红色
def create_white_to_red_cmap():
    """创建从灰白色到红色的渐变颜色映射"""
    colors = ['#F5F5F5', '#FFCCCC', '#FF9999', '#FF6666', '#FF3333', '#FF0000', '#CC0000']
    n_bins = 256
    return LinearSegmentedColormap.from_list('white_to_red', colors, N=n_bins)

# 辅助函数：绘制单一样本单基因的空间表达图
def plot_single_gene_sample_spatial(gene, sample, expr, spatial_coords, cell_count, 
                                   vmin, vmax, output_dir, white_to_red_cmap):
    """绘制单个基因在单个样本中的空间表达图"""
    # 创建输出目录
    gene_sample_dir = os.path.join(output_dir, "individual_plots", gene, sample)
    os.makedirs(gene_sample_dir, exist_ok=True)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # 动态调整点的大小
    if cell_count > 10000:
        point_size = 2
    elif cell_count > 5000:
        point_size = 3
    elif cell_count > 2000:
        point_size = 4
    else:
        point_size = 5
    
    # 绘制散点图
    scatter = ax.scatter(
        spatial_coords[:, 0],
        spatial_coords[:, 1],
        c=expr,
        cmap=white_to_red_cmap,
        s=point_size,
        edgecolors='none',
        alpha=0.7,
        vmin=vmin,
        vmax=vmax,
        linewidths=0
    )
    
    # 设置标题和标签
    ax.set_title(f'Sample: {sample} - Gene: {gene}', fontsize=16, fontweight='bold')
    ax.set_xlabel("X coordinate", fontsize=12)
    ax.set_ylabel("Y coordinate", fontsize=12)
    
    # 移除刻度
    ax.set_xticks([])
    ax.set_yticks([])
    
    # 设置坐标轴范围
    x_range = spatial_coords[:, 0].max() - spatial_coords[:, 0].min()
    y_range = spatial_coords[:, 1].max() - spatial_coords[:, 1].min()
    padding = 0.05
    
    if x_range > 0 and y_range > 0:
        ax.set_xlim([spatial_coords[:, 0].min() - padding * x_range, 
                   spatial_coords[:, 0].max() + padding * x_range])
        ax.set_ylim([spatial_coords[:, 1].min() - padding * y_range, 
                   spatial_coords[:, 1].max() + padding * y_range])
    
    # 添加颜色条
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Expression Level', fontsize=12)
    
    # 添加统计信息
    mean_expr = np.mean(expr)
    max_expr = np.max(expr)
    cells_positive = np.sum(expr > 0)
    percent_positive = cells_positive / len(expr) * 100
    std_expr = np.std(expr)
    
    stats_text = (f"Cells: {cell_count:,}\n"
                  f"Mean: {mean_expr:.3f}\n"
                  f"Max: {max_expr:.3f}\n"
                  f"Std: {std_expr:.3f}\n"
                  f"Positive: {percent_positive:.1f}%")
    
    ax.text(0.02, 0.98, stats_text,
           transform=ax.transAxes,
           fontsize=10,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    plt.tight_layout()
    
    # 保存为两种格式
    base_filename = f"{gene}_{sample}_spatial"
    
    # 保存为PNG
    png_path = os.path.join(gene_sample_dir, f"{base_filename}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    # 保存为PDF
    pdf_path = os.path.join(gene_sample_dir, f"{base_filename}.pdf")
    plt.savefig(pdf_path, bbox_inches='tight')
    
    plt.close(fig)
    
    return png_path, pdf_path, {
        'gene': gene,
        'sample': sample,
        'mean_expression': mean_expr,
        'max_expression': max_expr,
        'std_expression': std_expr,
        'cells_positive': cells_positive,
        'percent_positive': percent_positive,
        'total_cells': cell_count
    }

# 7. 为每个基因创建所有样本的分布图
print("\n开始为每个基因创建所有样本的分布图...")
print("="*60)

# 创建主输出目录
gene_wise_dir = os.path.join(output_dir, "gene_wise_plots")
individual_dir = os.path.join(output_dir, "individual_plots")
os.makedirs(gene_wise_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)

# 统计信息收集
all_gene_stats = []
individual_plot_stats = []  # 存储单个图的统计信息

# 创建自定义颜色映射
white_to_red_cmap = create_white_to_red_cmap()

# 对每个基因进行循环
for gene_idx, gene in enumerate(available_genes, 1):
    print(f"\n处理基因 {gene_idx}/{len(available_genes)}: {gene}")
    
    # 计算每个样本中该基因的表达统计
    gene_stats = []
    sample_expressions = {}  # 存储每个样本的表达数据
    
    for sample in samples:
        # 提取当前样本的数据
        sample_mask = raw_adata.obs[sample_col] == sample
        sample_adata = raw_adata[sample_mask].copy()
        
        if sample_adata.shape[0] == 0:
            continue
        
        # 获取基因表达
        expr = get_gene_expression(sample_adata, gene)
        
        if expr is not None:
            # 计算统计信息
            mean_expr = np.mean(expr)
            max_expr = np.max(expr)
            min_expr = np.min(expr)
            std_expr = np.std(expr)
            cells_detected = np.sum(expr > 0)
            total_cells = len(expr)
            detection_rate = cells_detected / total_cells * 100
            
            gene_stats.append({
                'Gene': gene,
                'Sample': sample,
                'Mean_Expression': mean_expr,
                'Max_Expression': max_expr,
                'Min_Expression': min_expr,
                'Std_Expression': std_expr,
                'Cells_Detected': cells_detected,
                'Total_Cells': total_cells,
                'Detection_Rate': detection_rate
            })
            
            # 保存表达数据用于后续可视化
            sample_expressions[sample] = {
                'expression': expr,
                'spatial_coords': sample_adata.obsm[spatial_key][:, :2] if spatial_key in sample_adata.obsm else None,
                'cell_count': total_cells
            }
    
    if not gene_stats:
        print(f"  警告: 基因 {gene} 在所有样本中均未检测到表达")
        continue
    
    # 创建基因统计DataFrame
    gene_stats_df = pd.DataFrame(gene_stats)
    
    # 保存基因的统计信息
    stats_path = os.path.join(gene_wise_dir, f"{gene}_statistics.csv")
    gene_stats_df.to_csv(stats_path, index=False)
    print(f"  统计信息已保存: {stats_path}")
    
    # 将统计信息添加到总表
    all_gene_stats.append(gene_stats_df)
    
    # 计算样本数量，确定图形布局
    n_samples = len(sample_expressions)
    
    if n_samples == 0:
        continue
    
    # 获取所有样本的表达范围，用于统一颜色映射
    all_expressions = []
    for sample_data in sample_expressions.values():
        all_expressions.extend(sample_data['expression'])
    
    if all_expressions:
        # 使用百分位数避免极端值影响颜色映射
        vmin = np.percentile(all_expressions, 1)  # 使用1st百分位数
        vmax = np.percentile(all_expressions, 99)  # 使用99th百分位数
        
        # 确保vmin和vmax不同
        if vmax - vmin < 0.001:
            vmax = vmin + 1
    else:
        vmin = 0
        vmax = 1
    
    # 第一部分：为每个样本创建单独的图
    print(f"  为基因 {gene} 创建 {n_samples} 个单独样本图...")
    for sample, sample_data in sample_expressions.items():
        expr = sample_data['expression']
        spatial_coords = sample_data['spatial_coords']
        cell_count = sample_data['cell_count']
        
        if spatial_coords is not None and len(expr) > 0:
            try:
                png_path, pdf_path, stats = plot_single_gene_sample_spatial(
                    gene, sample, expr, spatial_coords, cell_count,
                    vmin, vmax, output_dir, white_to_red_cmap
                )
                individual_plot_stats.append(stats)
                print(f"    ✓ 样本 {sample}: PNG和PDF已保存")
            except Exception as e:
                print(f"    ✗ 样本 {sample} 出错: {e}")
    
    # 第二部分：创建综合图（所有样本在一个图中）
    print(f"  创建基因 {gene} 的综合图...")
    
    # 计算布局：每行最多4个子图
    n_cols = min(4, n_samples)
    n_rows = int(np.ceil(n_samples / n_cols))
    
    # 创建图形
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5*n_cols, 4.5*n_rows))
    
    # 如果只有一个子图
    if n_samples == 1:
        axes = np.array([axes])
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    if n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # 设置主标题
    fig.suptitle(f'Gene: {gene} - Spatial Expression Across Samples', fontsize=16, y=1.02)
    
    # 为每个样本绘制子图
    for idx, (sample, sample_data) in enumerate(sample_expressions.items()):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        expr = sample_data['expression']
        spatial_coords = sample_data['spatial_coords']
        cell_count = sample_data['cell_count']
        
        if spatial_coords is not None and len(expr) > 0:
            # 动态调整点的大小
            if cell_count > 10000:
                point_size = 2
            elif cell_count > 5000:
                point_size = 3
            elif cell_count > 2000:
                point_size = 4
            else:
                point_size = 5
            
            # 绘制散点图
            scatter = ax.scatter(
                spatial_coords[:, 0],
                spatial_coords[:, 1],
                c=expr,
                cmap=white_to_red_cmap,
                s=point_size,
                edgecolors='none',
                alpha=0.7,
                vmin=vmin,
                vmax=vmax,
                linewidths=0
            )
            
            # 设置子图标题
            ax.set_title(f"Sample: {sample}", fontsize=12, fontweight='bold')
            ax.set_xlabel("X coordinate", fontsize=9)
            ax.set_ylabel("Y coordinate", fontsize=9)
            
            # 移除刻度
            ax.set_xticks([])
            ax.set_yticks([])
            
            # 设置坐标轴范围
            x_range = spatial_coords[:, 0].max() - spatial_coords[:, 0].min()
            y_range = spatial_coords[:, 1].max() - spatial_coords[:, 1].min()
            padding = 0.05
            
            if x_range > 0 and y_range > 0:
                ax.set_xlim([spatial_coords[:, 0].min() - padding * x_range, 
                           spatial_coords[:, 0].max() + padding * x_range])
                ax.set_ylim([spatial_coords[:, 1].min() - padding * y_range, 
                           spatial_coords[:, 1].max() + padding * y_range])
            
            # 添加表达统计
            mean_expr = np.mean(expr)
            cells_positive = np.sum(expr > 0)
            percent_positive = cells_positive / len(expr) * 100
            
            ax.text(0.02, 0.98, 
                   f"Cells: {cell_count:,}\nPos: {percent_positive:.1f}%", 
                   transform=ax.transAxes,
                   fontsize=8,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, f"No data\nfor {sample}", 
                   ha='center', va='center', transform=ax.transAxes, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
    
    # 隐藏多余的子图
    for idx in range(len(sample_expressions), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')
    
    # 添加全局颜色条
    if all_expressions:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        sm = plt.cm.ScalarMappable(cmap=white_to_red_cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Expression Level', fontsize=10)
    
    # 调整布局
    plt.tight_layout(rect=[0, 0, 0.9, 0.98])
    
    # 保存综合图为两种格式
    base_filename = f"{gene}_all_samples_spatial"
    
    # PNG格式
    png_path = os.path.join(gene_wise_dir, f"{base_filename}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    # PDF格式
    pdf_path = os.path.join(gene_wise_dir, f"{base_filename}.pdf")
    plt.savefig(pdf_path, bbox_inches='tight')
    
    plt.close(fig)
    print(f"  综合图已保存: {png_path} 和 {pdf_path}")
    
    # 第三部分：创建表达水平汇总图（箱线图）
    if len(gene_stats) > 1:
        print(f"  创建基因 {gene} 的统计图...")
        
        # 准备数据用于箱线图
        plot_data = []
        plot_labels = []
        
        for sample in samples:
            if sample in sample_expressions:
                expr = sample_expressions[sample]['expression']
                if len(expr) > 0:
                    plot_data.append(expr)
                    plot_labels.append(sample)
        
        if plot_data:
            # 创建箱线图
            fig2, ax2 = plt.subplots(figsize=(max(8, len(samples)*0.3), 6))
            
            boxprops = dict(facecolor='lightcoral', color='darkred', linewidth=1.5)
            whiskerprops = dict(color='darkred', linewidth=1.5)
            capprops = dict(color='darkred', linewidth=1.5)
            medianprops = dict(color='red', linewidth=2)
            
            bp = ax2.boxplot(plot_data, labels=plot_labels, showfliers=False,
                           patch_artist=True,
                           boxprops=boxprops, whiskerprops=whiskerprops,
                           capprops=capprops, medianprops=medianprops)
            
            ax2.set_facecolor('#f8f8f8')
            ax2.set_title(f'Gene: {gene} - Expression Distribution Across Samples', fontsize=14)
            ax2.set_ylabel('Expression Level', fontsize=12)
            ax2.set_xlabel('Sample', fontsize=12)
            plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')
            plt.tight_layout()
            
            # 保存为两种格式
            boxplot_base = f"{gene}_expression_boxplot"
            boxplot_png = os.path.join(gene_wise_dir, f"{boxplot_base}.png")
            boxplot_pdf = os.path.join(gene_wise_dir, f"{boxplot_base}.pdf")
            plt.savefig(boxplot_png, dpi=300, bbox_inches='tight')
            plt.savefig(boxplot_pdf, bbox_inches='tight')
            plt.close(fig2)
            print(f"  箱线图已保存: {boxplot_png} 和 {boxplot_pdf}")
        
        # 创建检测率条形图
        if gene_stats:
            fig3, ax3 = plt.subplots(figsize=(max(8, len(samples)*0.3), 6))
            
            detection_rates = []
            sample_labels = []
            
            for stat in gene_stats:
                detection_rates.append(stat['Detection_Rate'])
                sample_labels.append(stat['Sample'])
            
            if detection_rates:
                colors = plt.cm.Reds(np.linspace(0.3, 1, len(detection_rates)))
                bars = ax3.bar(range(len(detection_rates)), detection_rates, color=colors)
                
                ax3.set_facecolor('#f8f8f8')
                ax3.set_title(f'Gene: {gene} - Detection Rate Across Samples', fontsize=14)
                ax3.set_ylabel('Detection Rate (%)', fontsize=12)
                ax3.set_xlabel('Sample', fontsize=12)
                ax3.set_xticks(range(len(detection_rates)))
                ax3.set_xticklabels(sample_labels, rotation=45, ha='right')
                
                for i, bar in enumerate(bars):
                    height = bar.get_height()
                    ax3.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                           f'{detection_rates[i]:.1f}%', ha='center', va='bottom', fontsize=9)
                
                max_rate = max(detection_rates) if detection_rates else 100
                ax3.set_ylim([0, min(110, max_rate + 10)])
                
                plt.tight_layout()
                
                # 保存为两种格式
                detection_base = f"{gene}_detection_rate"
                detection_png = os.path.join(gene_wise_dir, f"{detection_base}.png")
                detection_pdf = os.path.join(gene_wise_dir, f"{detection_base}.pdf")
                plt.savefig(detection_png, dpi=300, bbox_inches='tight')
                plt.savefig(detection_pdf, bbox_inches='tight')
                plt.close(fig3)
                print(f"  检测率图已保存: {detection_png} 和 {detection_pdf}")

# 8. 保存单个图的统计信息
if individual_plot_stats:
    print("\n保存单个图的统计信息...")
    individual_stats_df = pd.DataFrame(individual_plot_stats)
    individual_stats_path = os.path.join(output_dir, "individual_plots_statistics.csv")
    individual_stats_df.to_csv(individual_stats_path, index=False)
    print(f"单个图统计信息已保存到: {individual_stats_path}")

# 9. 创建所有基因的综合汇总
print("\n" + "="*60)
print("创建所有基因的综合汇总...")

if all_gene_stats:
    # 合并所有基因的统计信息
    combined_stats = pd.concat(all_gene_stats, ignore_index=True)
    combined_path = os.path.join(output_dir, "all_genes_combined_statistics.csv")
    combined_stats.to_csv(combined_path, index=False)
    print(f"综合统计数据已保存到: {combined_path}")
    
    # 创建热图：每个基因在每个样本中的平均表达
    try:
        heatmap_data = combined_stats.pivot(index='Sample', columns='Gene', values='Mean_Expression')
        heatmap_data = heatmap_data.reindex(samples)
        
        # 绘制热图
        fig4, ax4 = plt.subplots(figsize=(max(10, len(available_genes)*0.8), max(6, len(samples)*0.5)))
        sns.heatmap(
            heatmap_data,
            cmap=white_to_red_cmap,
            annot=True,
            fmt='.2f',
            linewidths=0.5,
            cbar_kws={'label': 'Mean Expression'},
            square=True,
            ax=ax4
        )
        ax4.set_title('Mean Gene Expression Across Samples', fontsize=14, fontweight='bold')
        ax4.set_xlabel('Gene')
        ax4.set_ylabel('Sample')
        plt.tight_layout()
        
        # 保存为两种格式
        heatmap_png = os.path.join(output_dir, "all_genes_heatmap.png")
        heatmap_pdf = os.path.join(output_dir, "all_genes_heatmap.pdf")
        plt.savefig(heatmap_png, dpi=300, bbox_inches='tight')
        plt.savefig(heatmap_pdf, bbox_inches='tight')
        plt.close(fig4)
        print(f"热图已保存到: {heatmap_png} 和 {heatmap_pdf}")
        
        # 保存热图数据
        heatmap_data_path = os.path.join(output_dir, "all_genes_expression_matrix.csv")
        heatmap_data.to_csv(heatmap_data_path)
        print(f"表达矩阵已保存到: {heatmap_data_path}")
        
    except Exception as e:
        print(f"创建热图时出错: {e}")
    
    # 创建检测率热图
    try:
        detection_heatmap_data = combined_stats.pivot(index='Sample', columns='Gene', values='Detection_Rate')
        detection_heatmap_data = detection_heatmap_data.reindex(samples)
        
        fig5, ax5 = plt.subplots(figsize=(max(10, len(available_genes)*0.8), max(6, len(samples)*0.5)))
        sns.heatmap(
            detection_heatmap_data,
            cmap='Blues',
            annot=True,
            fmt='.1f',
            linewidths=0.5,
            cbar_kws={'label': 'Detection Rate (%)'},
            square=True,
            ax=ax5
        )
        ax5.set_title('Gene Detection Rate Across Samples', fontsize=14, fontweight='bold')
        ax5.set_xlabel('Gene')
        ax5.set_ylabel('Sample')
        plt.tight_layout()
        
        # 保存为两种格式
        detection_heatmap_png = os.path.join(output_dir, "all_genes_detection_rate_heatmap.png")
        detection_heatmap_pdf = os.path.join(output_dir, "all_genes_detection_rate_heatmap.pdf")
        plt.savefig(detection_heatmap_png, dpi=300, bbox_inches='tight')
        plt.savefig(detection_heatmap_pdf, bbox_inches='tight')
        plt.close(fig5)
        print(f"检测率热图已保存到: {detection_heatmap_png} 和 {detection_heatmap_pdf}")
        
    except Exception as e:
        print(f"创建检测率热图时出错: {e}")

# 10. 创建配置文件记录运行参数
config_info = {
    "input_file": "merged_TMA_data.h5ad",
    "output_dir": output_dir,
    "gene_wise_dir": gene_wise_dir,
    "individual_dir": individual_dir,
    "total_cells": raw_adata.shape[0],
    "total_genes": raw_adata.shape[1],
    "samples_found": len(samples),
    "genes_requested": len(gene_list),
    "genes_found": len(available_genes),
    "sample_column": sample_col,
    "spatial_coordinates": spatial_key,
    "used_raw_matrix": adata.raw is not None,
    "output_formats": "PNG and PDF",
    "visualization_settings": "点大小自适应(2-5), 颜色映射:灰白到红, 透明度:0.7",
    "individual_plots_generated": len(individual_plot_stats) if individual_plot_stats else 0
}

config_df = pd.DataFrame([config_info])
config_path = os.path.join(output_dir, "analysis_configuration.csv")
config_df.to_csv(config_path, index=False)
print(f"\n配置文件已保存到: {config_path}")

print("\n" + "="*60)
print("脚本执行完成!")
print(f"所有结果已保存到 '{output_dir}' 目录")
print(f"输出目录结构:")
print(f"  {output_dir}/")
print(f"  ├── gene_wise_plots/       # 按基因组织的综合图")
print(f"  │   ├── {available_genes[0]}_all_samples_spatial.png/pdf")
print(f"  │   ├── {available_genes[0]}_expression_boxplot.png/pdf")
print(f"  │   ├── {available_genes[0]}_detection_rate.png/pdf")
print(f"  │   ├── {available_genes[0]}_statistics.csv")
print(f"  │   └── ...")
print(f"  ├── individual_plots/      # 每个基因-样本的单独图")
print(f"  │   ├── {available_genes[0]}/")
print(f"  │   │   ├── {samples[0] if samples else 'sample'}/")
print(f"  │   │   │   ├── {available_genes[0]}_{samples[0] if samples else 'sample'}_spatial.png")
print(f"  │   │   │   └── {available_genes[0]}_{samples[0] if samples else 'sample'}_spatial.pdf")
print(f"  │   │   └── ...")
print(f"  │   └── ...")
print(f"  ├── all_genes_heatmap.png/pdf")
print(f"  ├── all_genes_detection_rate_heatmap.png/pdf")
print(f"  ├── all_genes_expression_matrix.csv")
print(f"  ├── all_genes_combined_statistics.csv")
print(f"  ├── individual_plots_statistics.csv")
print(f"  └── analysis_configuration.csv")
print(f"\n统计信息:")
print(f"  - 共处理了 {len(samples)} 个样本")
print(f"  - 成功分析了 {len(available_genes)} 个基因")
print(f"  - 生成了 {len(individual_plot_stats) if individual_plot_stats else 0} 个单独样本图")
print(f"  - 所有图都保存为PNG和PDF两种格式")
print("="*60)
