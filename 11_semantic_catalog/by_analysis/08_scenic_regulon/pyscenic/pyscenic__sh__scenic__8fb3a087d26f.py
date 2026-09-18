import numpy as np
import scanpy as sc
import loompy as lp
import os
import argparse
import pandas as pd
from anndata import AnnData

def process_and_subset_adata(ad_all, output_dir, sample_fraction=0.3, min_cells=10):
    """
    处理 adata 数据：抽样、过滤基因、保存结果
    Args:
        ad_all: 原始 AnnData 对象
        output_dir: 输出目录路径
        sample_fraction: 抽样比例（默认0.3）
        min_cells: 基因最小表达细胞数（默认10）
    Returns:
        adata: 处理后的 AnnData 对象
    """
    # 检查细胞数量是否大于4万
    if ad_all.n_obs > 40000:
        print(f"细胞数 {ad_all.n_obs} > 40,000，按比例 {sample_fraction} 抽样...")
        n_samples = int(ad_all.n_obs * sample_fraction)
        random_indices = np.random.choice(ad_all.n_obs, n_samples, replace=False)
        adata = ad_all[random_indices, :].copy()
    else:
        print(f"细胞数 {ad_all.n_obs} ≤ 40,000，跳过抽样")
        adata = ad_all.copy()
    
    # 过滤低表达基因
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"过滤后基因数: {adata.n_vars}, 细胞数: {adata.n_obs}")
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存处理后的数据
    adata.write(os.path.join(output_dir, 'sub_anno.h5ad'), compression="gzip")
    return adata

def create_loom_for_scenic(adata, output_dir):
    """
    为 SCENIC 创建 loom 文件，检查并使用 adata.raw 属性
    Args:
        adata: AnnData 对象
        output_dir: 输出目录
    """
    loom_path = os.path.join(output_dir, "sce.loom")
    
    # 检查是否有 raw 属性
    if hasattr(adata, 'raw') and adata.raw is not None:
        print("使用 adata.raw 创建 loom 文件")
        row_attrs = {"Gene": np.array(adata.raw.var_names)}
        col_attrs = {
            "CellID": np.array(adata.raw.obs_names),
            "nGene": np.array(adata.raw.X.getnnz(axis=1))
        }
        lp.create(loom_path, adata.raw.X.transpose(), row_attrs, col_attrs)
    else:
        print("使用 adata 创建 loom 文件")
        row_attrs = {"Gene": np.array(adata.var_names)}
        col_attrs = {
            "CellID": np.array(adata.obs_names),
            "nGene": np.array(adata.X.getnnz(axis=1))
        }
        lp.create(loom_path, adata.X.transpose(), row_attrs, col_attrs)
    
    print(f"loom 文件已保存: {loom_path}")

def save_umap_and_metadata(adata, output_dir):
    """
    保存已计算的UMAP坐标和metadata到CSV文件
    Args:
        adata: AnnData 对象，需已包含UMAP坐标
        output_dir: 输出目录
    """
    # 检查是否存在UMAP坐标
    if 'X_umap' not in adata.obsm:
        raise ValueError("输入的AnnData对象不包含UMAP坐标，请确保已计算UMAP")
    
    # 保存UMAP坐标
    umap_path = os.path.join(output_dir, 'umap_coordinates.csv')
    umap_df = pd.DataFrame(adata.obsm['X_umap'], 
                          index=adata.obs_names, 
                          columns=['UMAP_1', 'UMAP_2'])
    umap_df.to_csv(umap_path)
    print(f"UMAP坐标已保存: {umap_path}")
    
    # 保存metadata
    metadata_path = os.path.join(output_dir, 'metadata.csv')
    metadata_df = adata.obs
    metadata_df.to_csv(metadata_path)
    print(f"metadata已保存: {metadata_path}")

def run_pyscenic(loom_path, output_dir, species="human"):
    """
    运行 pySCENIC 分析
    Args:
        loom_path: 输入 loom 文件路径
        output_dir: 输出目录
        species: 物种，可选值为 "human" 或 "mouse"
    """
    print('=====================>pyscenic will be running soon===========================')
    # 定义输出文件路径
    grn_output = os.path.join(output_dir, "grn_result.csv")
    ctx_output = os.path.join(output_dir, "ctx_result.csv")
    auc_output = os.path.join(output_dir, "aue_result.loom")
    
    # 根据物种选择对应的SCENIC脚本
    if species.lower() == "human":
        scenic_script = "/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/pyscenic.sh"
    elif species.lower() == "mouse":
        scenic_script = "/mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC/new_code/pyscenic_mouse.sh"
    else:
        raise ValueError(f"不支持的物种: {species}，请选择 'human' 或 'mouse'")
    
    # 检查 SCENIC 脚本是否存在
    if not os.path.exists(scenic_script):
        raise FileNotFoundError(f"SCENIC 脚本未找到: {scenic_script}")
    
    # 执行 SCENIC 命令
    cmd = f"sh {scenic_script} {loom_path} {grn_output} {ctx_output} {auc_output} {output_dir}"
    print("运行 SCENIC 命令:", cmd)
    os.system(cmd)

# =============== 主流程 ===============
if __name__ == "__main__":
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description='执行SCENIC分析流程')
    parser.add_argument('--adata_path', required=True, help='输入AnnData文件的绝对路径')
    parser.add_argument('--output_dir', required=True, help='输出目录的绝对路径')
    parser.add_argument('--species', required=True, choices=['human', 'mouse'], help='物种：human或mouse')
    parser.add_argument('--sample_fraction', type=float, default=0.3, help='抽样比例（默认0.3）')
    
    args = parser.parse_args()
    
    # 加载数据
    print(f"加载数据: {args.adata_path}")
    ad_all = sc.read(args.adata_path)
    # ad_all = ad_all[ad_all.obs['cell_type'].isin(['Non-classical Monocytes']), :]
    # 1. 处理数据（抽样+过滤）
    adata = process_and_subset_adata(ad_all, args.output_dir, sample_fraction=args.sample_fraction)
    
    # 2. 保存已计算的UMAP坐标和metadata
    save_umap_and_metadata(adata, args.output_dir)
    
    # 3. 创建 loom 文件（检查 raw 属性）
    create_loom_for_scenic(adata, args.output_dir)
    
    # 4. 运行 SCENIC（根据物种选择脚本）
    run_pyscenic(os.path.join(args.output_dir, "sce.loom"), args.output_dir, species=args.species)
