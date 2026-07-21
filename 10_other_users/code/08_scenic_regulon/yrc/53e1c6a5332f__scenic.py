import numpy as np
import scanpy as sc
import loompy as lp
import os
import argparse
import pandas as pd
from anndata import AnnData

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
    parser.add_argument('--adata_path', required=False, help='输入AnnData文件的绝对路径')
    parser.add_argument('--output_dir', required=True, help='输出目录的绝对路径')
    parser.add_argument('--species', required=True, choices=['human', 'mouse'], help='物种：human或mouse')
    parser.add_argument('--sample_fraction', type=float, default=0.3, help='抽样比例（默认0.3）')
    
    args = parser.parse_args()
    
    # 4. 运行 SCENIC（根据物种选择脚本）
    run_pyscenic(os.path.join(args.output_dir, "sce.loom"), args.output_dir, species=args.species)
