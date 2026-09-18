# -*- coding: utf-8 -*-
"""
角化鳞状细胞微环境分析脚本
功能：计算50μm半径内细胞组成比例，绘制分组箱线图
依赖：scanpy==1.10.0, scipy==1.13.0, pandas==2.2.2, seaborn==0.13.2
"""
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from tqdm import tqdm

# ----------------------------
# 1. 数据加载与预处理
# ----------------------------
def load_and_validate(h5ad_path: str):
    """加载h5ad数据并进行完整性校验"""
    adata = sc.read_h5ad(h5ad_path)
    
    # 关键数据校验 [1,2,3](@ref)
    assert 'spatial' in adata.obsm, "空间坐标缺失: 需存储在adata.obsm['spatial']"
    assert 'celltype' in adata.obs, "细胞类型列缺失: 需命名为'celltype'"
    assert 'group' in adata.obs, "分组信息列缺失: 需命名为'group'"
    assert "CCL19+ CAF" in adata.obs['celltype'].unique(), (
        "目标细胞类型未找到: 请检查'Keratinizing Squamous Cells'是否存在于celltype列")
    
    print(f"✔ 数据加载成功: {adata.shape[0]}个细胞, {adata.shape[1]}个基因")
    return adata

# ----------------------------
# 2. 邻域细胞比例计算
# ----------------------------
def calculate_neighborhood_ratios(
    adata, 
    target_cell: str = "Keratinizing Squamous Cells", 
    radius: float = 50.0
) -> pd.DataFrame:
    """
    计算目标细胞周围半径内的细胞组成
    返回包含比例和分组信息的DataFrame
    """
    # 提取目标细胞坐标 [2,3](@ref)
    target_mask = adata.obs['celltype'] == target_cell
    target_coords = adata.obsm['spatial'][target_mask]
    all_coords = adata.obsm['spatial']
    
    # 构建空间索引树 [3,4](@ref)
    kdtree = cKDTree(all_coords)
    celltypes = adata.obs['celltype'].cat.categories
    
    # 收集结果
    results = []
    group_labels = adata.obs.loc[target_mask, 'group'].values
    
    # 遍历每个目标细胞 [4](@ref)
    for i, coord in enumerate(tqdm(target_coords, desc="分析微环境")):
        # 查询半径内细胞索引 (排除自身)
        neighbor_idxs = kdtree.query_ball_point(coord, r=radius)
        neighbor_idxs = [idx for idx in neighbor_idxs if not target_mask[idx]]
        
        # 计算细胞类型比例
        neighbor_types = adata.obs.iloc[neighbor_idxs]['celltype']
        type_counts = neighbor_types.value_counts(normalize=True)
        
        # 构建结果记录
        record = {'cell_id': adata.obs.index[target_mask][i], 'group': group_labels[i]}
        for ct in celltypes:
            record[ct] = type_counts.get(ct, 0.0)  # 未出现类型设为0
        results.append(record)
    
    return pd.DataFrame(results)

# ----------------------------
# 3. 数据可视化
# ----------------------------
def plot_neighborhood_composition(df: pd.DataFrame, min_proportion: float = 0.05):
    """绘制分组箱线图展示微环境组成"""
    # 转换长格式并过滤稀有细胞 [5](@ref)
    plot_df = df.melt(
        id_vars=['group', 'cell_id'],
        var_name='celltype',
        value_name='proportion'
    )
    
    # 过滤低频细胞类型
    valid_types = plot_df.groupby('celltype')['proportion'].mean()
    valid_types = valid_types[valid_types > min_proportion].index
    plot_df = plot_df[plot_df['celltype'].isin(valid_types)]
    
    # 绘图设置
    plt.figure(figsize=(14, 8))
    sns.set_style("whitegrid", {'grid.linestyle': ':'})
    palette = sns.color_palette("Set2", n_colors=plot_df['group'].nunique())
    
    # 绘制分组箱线图 [5](@ref)
    ax = sns.boxplot(
        data=plot_df,
        x='celltype',
        y='proportion',
        hue='group',
        palette=palette,
        linewidth=2.5,
        width=0.7,
        showfliers=False,  # 不显示异常值
        dodge=True        # 分组并列显示
    )
    
    # 美化图形
    plt.title(f"Keratinizing Squamous Cells", 
              fontsize=16, pad=20)
    plt.xlabel("Cell Type", fontsize=14)
    plt.ylabel("Percentage", fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title='Group', frameon=False, fontsize=12)
    
    # 添加统计标识
    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(f"Keratin_microenvironment_{radius}um.pdf", dpi=300)
    plt.show()

# ----------------------------
# 主执行流程
# ----------------------------
if __name__ == "__main__":
    # 配置参数
    h5ad_path = "merged_single_cell_data.h5ad"  # 替换为实际路径
    target_cell = "CCL19+ CAF"
    radius = 200.0  # 微环境半径(μm)
    
    # 执行分析流程
    adata = load_and_validate(h5ad_path)
    neighbor_df = calculate_neighborhood_ratios(adata, target_cell, radius)
    plot_neighborhood_composition(neighbor_df)
    
    # 保存结果
    neighbor_df.to_csv(f"neighborhood_composition_{radius}um.csv", index=False)
    print("✅ 分析完成！结果已保存至CSV文件和PNG图片")
