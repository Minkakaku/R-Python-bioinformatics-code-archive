import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from scipy.interpolate import griddata, Rbf
from scipy.ndimage import gaussian_filter
import matplotlib.patches as mpatches
import os
import warnings
warnings.filterwarnings('ignore')

def create_output_folder(folder_name='spatial_muc19_smooth_contour'):
    """
    创建输出文件夹
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

def check_required_columns(adata):
    """
    检查必要的列是否存在
    """
    required_cols = ['sample', 'mucus_zone']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        raise ValueError(f"adata.obs中未找到以下列: {missing_cols}")
    
    if 'spatial' not in adata.obsm:
        raise ValueError("adata.obsm中未找到'spatial'坐标信息")
    
    print(f"数据形状: {adata.shape}")
    print(f"样本数量: {adata.obs['sample'].nunique()}")
    print(f"空间坐标维度: {adata.obsm['spatial'].shape}")

def get_normalized_muc19_expression(adata):
    """
    从adata.raw矩阵中获取标准化后的MUC19表达量
    """
    # 检查是否有raw矩阵
    if adata.raw is None:
        print("警告: adata没有raw矩阵，尝试使用adata.X")
        exp_matrix = adata.X
        var_names = adata.var_names
    else:
        print("使用adata.raw矩阵（标准化数据）")
        exp_matrix = adata.raw.X
        var_names = adata.raw.var_names
    
    # 检查MUC19是否存在
    muc19_genes = [gene for gene in var_names if 'MUC19' in gene.upper()]
    
    if not muc19_genes:
        print("警告: 未找到MUC19基因")
        # 检查类似的基因名
        all_genes = list(var_names)
        muc_genes = [gene for gene in all_genes if 'MUC' in gene.upper()]
        print(f"找到的MUC基因: {muc_genes[:10]}")
        return None, None
    
    # 使用第一个匹配的MUC19基因
    muc19_gene = muc19_genes[0]
    print(f"使用基因: {muc19_gene}")
    
    # 获取表达量
    if hasattr(exp_matrix, 'toarray'):
        exp_matrix = exp_matrix.toarray()
    
    # 找到基因的索引
    gene_idx = list(var_names).index(muc19_gene)
    
    # 提取表达量
    muc19_expression = exp_matrix[:, gene_idx]
    
    # 如果是稀疏格式，转换为密集数组
    if hasattr(muc19_expression, 'toarray'):
        muc19_expression = muc19_expression.toarray().flatten()
    elif hasattr(muc19_expression, 'A1'):
        muc19_expression = muc19_expression.A1
    elif isinstance(muc19_expression, np.matrix):
        muc19_expression = muc19_expression.A1
    
    # 确保是numpy数组
    muc19_expression = np.array(muc19_expression)
    
    print(f"MUC19表达量统计 (标准化后):")
    print(f"  最小值: {muc19_expression.min():.4f}")
    print(f"  最大值: {muc19_expression.max():.4f}")
    print(f"  平均值: {muc19_expression.mean():.4f}")
    print(f"  中位数: {np.median(muc19_expression):.4f}")
    
    return muc19_expression, muc19_gene

def create_smooth_contour_grid(x, y, z, grid_resolution=500, method='linear', smooth_sigma=3.0):
    """
    创建非常平滑的等高线网格
    """
    # 创建高分辨率网格
    xi = np.linspace(x.min(), x.max(), grid_resolution)
    yi = np.linspace(y.min(), y.max(), grid_resolution)
    xi, yi = np.meshgrid(xi, yi)
    
    # 插值到网格
    try:
        # 使用径向基函数插值获得更平滑的结果
        if len(z) > 10:  # 确保有足够的数据点
            # 尝试RBF插值
            rbf = Rbf(x, y, z, function='multiquadric', smooth=0.1)
            zi = rbf(xi, yi)
        else:
            # 如果点太少，使用griddata
            zi = griddata((x, y), z, (xi, yi), method=method, fill_value=0)
        
        # 应用高斯平滑
        if smooth_sigma > 0:
            zi = gaussian_filter(zi, sigma=smooth_sigma)
        
        return xi, yi, zi, True
    except Exception as e:
        print(f"RBF插值失败: {e}")
        # 回退到简单方法
        try:
            zi = griddata((x, y), z, (xi, yi), method=method, fill_value=0)
            if smooth_sigma > 0:
                zi = gaussian_filter(zi, sigma=smooth_sigma)
            return xi, yi, zi, True
        except Exception as e2:
            print(f"网格插值失败: {e2}")
            return None, None, None, False

def plot_sample_spatial_with_smooth_contour(adata_sample, sample_name, global_muc19_expression, muc19_gene, 
                                           adata_indices, save_dir='spatial_muc19_smooth_contour', dpi=300):
    """
    为单个样本绘制空间图，使用非常平滑的等高线
    """
    # 获取坐标
    coords = adata_sample.obsm['spatial']
    x = coords[:, 0]
    y = coords[:, 1]
    
    # 获取mucus_zone信息
    mucus_zone = adata_sample.obs['mucus_zone'].values
    
    # 获取当前样本的MUC19表达量
    if global_muc19_expression is not None:
        sample_indices = adata_sample.obs.index
        sample_mask = adata_indices.isin(sample_indices)
        sample_muc19_exp = global_muc19_expression[sample_mask]
    else:
        sample_muc19_exp = None
    
    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax1, ax2, ax3, ax4 = axes.flatten()
    
    # 图1：基本空间图（mucus_zone着色）
    color_map = {'inner': 'darkgray', 'outer': 'lightgray'}
    
    # 绘制inner细胞
    inner_mask = mucus_zone == 'inner'
    if inner_mask.any():
        ax1.scatter(x[inner_mask], y[inner_mask], 
                   c=color_map['inner'], s=8, alpha=0.7, 
                   label='inner', edgecolors='none')
    
    # 绘制outer细胞
    outer_mask = mucus_zone == 'outer'
    if outer_mask.any():
        ax1.scatter(x[outer_mask], y[outer_mask], 
                   c=color_map['outer'], s=8, alpha=0.7, 
                   label='outer', edgecolors='none')
    
    ax1.set_title(f'Sample: {sample_name}\nMucus Zone Distribution', fontsize=12)
    ax1.set_xlabel('X coordinate', fontsize=10)
    ax1.set_ylabel('Y coordinate', fontsize=10)
    ax1.set_aspect('equal')
    
    # 添加图例
    inner_patch = mpatches.Patch(color=color_map['inner'], label='inner')
    outer_patch = mpatches.Patch(color=color_map['outer'], label='outer')
    ax1.legend(handles=[inner_patch, outer_patch], loc='upper right')
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # 图2：MUC19表达量散点图
    if sample_muc19_exp is not None and len(sample_muc19_exp) > 0:
        # 创建颜色映射
        norm = colors.Normalize(vmin=sample_muc19_exp.min(), vmax=sample_muc19_exp.max())
        cmap = plt.cm.viridis
        
        # 绘制散点图
        scatter = ax2.scatter(x, y, c=sample_muc19_exp, s=8, 
                             cmap=cmap, norm=norm, alpha=0.7, edgecolors='none')
        
        ax2.set_title(f'Sample: {sample_name}\n{muc19_gene} Expression', fontsize=12)
        ax2.set_xlabel('X coordinate', fontsize=10)
        ax2.set_ylabel('Y coordinate', fontsize=10)
        ax2.set_aspect('equal')
        
        # 添加颜色条
        cbar = plt.colorbar(scatter, ax=ax2)
        cbar.set_label(f'{muc19_gene} Expression (Normalized)', fontsize=10)
        ax2.grid(True, alpha=0.3, linestyle='--')
    else:
        ax2.text(0.5, 0.5, 'No MUC19 expression data', 
                transform=ax2.transAxes, ha='center', va='center')
        ax2.set_title(f'Sample: {sample_name}', fontsize=12)
    
    # 图3：非常平滑的MUC19表达量等高线图
    if sample_muc19_exp is not None and len(sample_muc19_exp) > 0:
        # 创建非常平滑的等高线网格
        xi, yi, zi, success = create_smooth_contour_grid(x, y, sample_muc19_exp, 
                                                          grid_resolution=500, method='cubic', smooth_sigma=3.0)
        
        if success:
            # 绘制散点图（根据mucus_zone着色）
            if inner_mask.any():
                ax3.scatter(x[inner_mask], y[inner_mask], 
                           c=color_map['inner'], s=5, alpha=0.7, 
                           edgecolors='none')
            
            if outer_mask.any():
                ax3.scatter(x[outer_mask], y[outer_mask], 
                           c=color_map['outer'], s=5, alpha=0.7, 
                           edgecolors='none')
            
            # 绘制非常平滑的等高线
            exp_min = sample_muc19_exp.min()
            exp_max = sample_muc19_exp.max()
            
            # 使用更多的等高线级别
            if exp_max - exp_min < 0.1:
                contour_levels = np.linspace(exp_min, exp_max, 20)
            else:
                # 使用更密集的等高线级别
                contour_levels = np.linspace(exp_min, exp_max, 30)
            
            # 绘制等高线 - 使用更细的线宽和更柔和的颜色
            contour_plot = ax3.contour(xi, yi, zi, levels=contour_levels, 
                                       colors='red', linewidths=1.0, alpha=0.7)
            
            # 添加等高线标签
            ax3.clabel(contour_plot, inline=True, fontsize=7, fmt='%.3f')
            
            ax3.set_title(f'Sample: {sample_name}\n{muc19_gene} Smooth Contour', fontsize=12)
            ax3.set_xlabel('X coordinate', fontsize=10)
            ax3.set_ylabel('Y coordinate', fontsize=10)
            ax3.set_aspect('equal')
            ax3.grid(True, alpha=0.3, linestyle='--')
        else:
            ax3.text(0.5, 0.5, 'Contour generation failed', 
                    transform=ax3.transAxes, ha='center', va='center')
            ax3.set_title(f'Sample: {sample_name}', fontsize=12)
    else:
        ax3.text(0.5, 0.5, 'No MUC19 data', 
                transform=ax3.transAxes, ha='center', va='center')
        ax3.set_title(f'Sample: {sample_name}', fontsize=12)
    
    # 图4：平滑的MUC19表达量热图叠加
    if sample_muc19_exp is not None and len(sample_muc19_exp) > 0 and success:
        # 绘制平滑的热图
        im = ax4.contourf(xi, yi, zi, levels=30, cmap='viridis', alpha=0.8)
        
        # 叠加散点图
        if inner_mask.any():
            ax4.scatter(x[inner_mask], y[inner_mask], 
                       c=color_map['inner'], s=5, alpha=0.7, 
                       edgecolors='black', linewidth=0.5)
        
        if outer_mask.any():
            ax4.scatter(x[outer_mask], y[outer_mask], 
                       c=color_map['outer'], s=5, alpha=0.7, 
                       edgecolors='black', linewidth=0.5)
        
        ax4.set_title(f'Sample: {sample_name}\n{muc19_gene} Expression Heatmap', fontsize=12)
        ax4.set_xlabel('X coordinate', fontsize=10)
        ax4.set_ylabel('Y coordinate', fontsize=10)
        ax4.set_aspect('equal')
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax4)
        cbar.set_label(f'{muc19_gene} Expression (Normalized)', fontsize=10)
        ax4.grid(True, alpha=0.3, linestyle='--')
    else:
        ax4.text(0.5, 0.5, 'Heatmap generation failed or no MUC19 data', 
                transform=ax4.transAxes, ha='center', va='center')
        ax4.set_title(f'Sample: {sample_name}', fontsize=12)
    
    plt.suptitle(f'Sample: {sample_name} - Spatial Analysis with Smooth MUC19 Expression', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存图形
    fig_path_png = os.path.join(save_dir, f'sample_{sample_name}_smooth_contour.png')
    fig_path_pdf = os.path.join(save_dir, f'sample_{sample_name}_smooth_contour.pdf')
    
    fig.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
    fig.savefig(fig_path_pdf, bbox_inches='tight')
    
    print(f"  保存图形: {fig_path_png}")
    
    plt.close(fig)
    
    return fig_path_png

def main():
    """
    主函数
    """
    # 数据路径
    h5ad_path = 'extracted_samples.h5ad'
    
    # 检查文件是否存在
    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"未找到h5ad文件: {h5ad_path}")
    
    # 创建输出文件夹
    save_dir = create_output_folder('spatial_muc19_smooth_contour')
    
    print("=" * 60)
    print("     Spatial Analysis with Smooth MUC19 Expression Contours")
    print("=" * 60)
    
    # 1. 加载数据
    print(f"\n加载数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    
    # 2. 检查必要的列
    check_required_columns(adata)
    
    # 3. 获取标准化后的MUC19表达量
    muc19_result = get_normalized_muc19_expression(adata)
    if muc19_result is None:
        muc19_expression = None
        muc19_gene = 'MUC19'
        print("警告: 无法获取MUC19表达量，将只绘制空间分布图")
    else:
        muc19_expression, muc19_gene = muc19_result
    
    # 4. 获取所有样本
    samples = adata.obs['sample'].unique()
    print(f"\n发现 {len(samples)} 个样本")
    print(f"样本列表: {list(samples)}")
    
    # 5. 为每个样本绘制图形
    print(f"\n开始为每个样本绘制空间图...")
    saved_files = []
    
    for sample in samples:
        print(f"\n处理样本: {sample}")
        
        # 提取当前样本的数据
        sample_mask = adata.obs['sample'] == sample
        adata_sample = adata[sample_mask]
        
        print(f"  细胞数量: {len(adata_sample)}")
        print(f"  inner细胞: {(adata_sample.obs['mucus_zone'] == 'inner').sum()}")
        print(f"  outer细胞: {(adata_sample.obs['mucus_zone'] == 'outer').sum()}")
        
        # 绘制图形
        fig_path = plot_sample_spatial_with_smooth_contour(
            adata_sample, sample, muc19_expression, muc19_gene, 
            adata.obs.index, save_dir
        )
        saved_files.append(fig_path)
    
    print("\n" + "=" * 60)
    print("分析完成！")
    print(f"所有结果已保存到: {save_dir}")
    print("=" * 60)
    print(f"生成的图形文件:")
    print(f"  - 每个样本的单独图: {len(samples)} 个")
    print(f"  - 文件格式: sample_样本名_smooth_contour.png/pdf")
    print("=" * 60)
    
    return adata, muc19_expression, muc19_gene

# 运行主函数
if __name__ == "__main__":
    adata, muc19_expression, muc19_gene = main()
