import sys
import os
import argparse
import warnings

# ================= 配置区域 (Hardcoded Paths) =================
# 1. CellPhoneDB 数据库路径 (写死)
DB_FILE_PATH = "/home/lyh/Script/cellphone_db/db/v5.0/cellphonedb.zip"

# 2. 绘图函数脚本所在路径
PLOT_FUNC_PATH = "/home/lyh/Script/cellphone_db/plot_fuction"
# ============================================================

def parse_args():
    """
    只处理参数解析，不依赖任何重型库，保证瞬间响应
    """
    parser = argparse.ArgumentParser(description="Run CellPhoneDB (Omicverse) with optional splitting.")
    
    parser.add_argument("-i", "--input", type=str, required=True, 
                        help="Path to input .h5ad file")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Path to output directory")
    parser.add_argument("-c", "--celltype", type=str, required=True, 
                        help="Column name for cell type annotation (e.g., sub2)")
    parser.add_argument("-s", "--split_col", type=str, default=None, required=False, 
                        help="[Optional] Column to split data by (e.g., Group). If not provided, runs on all data.")
    
    return parser.parse_args()

def main():
    # 1. 先解析参数 (瞬间完成)
    args = parse_args()
    
    # 2. 基础路径检查
    if not os.path.exists(DB_FILE_PATH):
        print(f"Error: Database file not found at: {DB_FILE_PATH}")
        sys.exit(1)

    # 3. 参数验证通过后，才开始加载重型库
    print("Parameters parsed. Loading libraries (Scanpy, Omicverse)...")
    
    # --- 延迟导入区域 START ---
    try:
        import scanpy as sc
        import omicverse as ov
        import pandas as pd
        import numpy as np
        import traceback # 用于打印详细错误
        
        # 忽略警告
        warnings.filterwarnings('ignore')
        
        # 导入自定义绘图函数
        process_cellchat_analysis = None
        if os.path.exists(PLOT_FUNC_PATH):
            sys.path.append(PLOT_FUNC_PATH)
            try:
                # 假设绘图库里有 import scanpy 等，所以也放在这里导入
                from cellphonedb_plot import process_cellchat_analysis
            except ImportError:
                print(f"Warning: Could not import 'process_cellchat_analysis' from {PLOT_FUNC_PATH}")
            except Exception as e:
                 print(f"Warning: Error importing plot function: {e}")
        else:
            print(f"Warning: Path {PLOT_FUNC_PATH} does not exist. Plotting functions might fail.")
            
    except ImportError as e:
        print(f"Critical Error: Missing required library. {e}")
        sys.exit(1)
    # --- 延迟导入区域 END ---

    # 初始化 Omicverse
    ov.plot_set()
    print(f"Omicverse version: {ov.__version__}")

    # 颜色定义
    color3 = [
        '#C04851', '#2775B6', '#61AC85', '#F1CA17', '#FC7930',
        '#22A2C3', '#E3B4B8', '#61649F', '#F7DE98', '#621624',
        '#E16C96', '#BDAEAD', '#8CC269', '#FBCD31', '#5CB3CC',
        '#EE4866', '#624941', '#F9CB8B', '#70887D', '#C08EAF',
        '#51C4D3', '#2D0C13', '#E5D3AA', '#FF9900', '#CAD3C3',
        '#863020', '#ea7293', '#80766e', '#f8bc31', '#2f90b9'
    ]

    # 参数赋值
    adata_path = args.input
    Output = args.output
    celltype_col = args.celltype
    split_col = args.split_col
    
    # 设置工作目录
    if not os.path.exists(Output):
        os.makedirs(Output)
    os.chdir(Output)
    print("Output dir is " + os.getcwd())

    # 读取数据
    print(f"Loading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    
    # 设置颜色字典
    unique_celltypes = adata.obs[celltype_col].unique()
    if len(unique_celltypes) > len(color3):
        print("Warning: More cell types than defined colors. Colors will recycle.")
        extended_colors = color3 * (len(unique_celltypes) // len(color3) + 1)
        color_dict = dict(zip(unique_celltypes, extended_colors[:len(unique_celltypes)]))
    else:
        color_dict = dict(zip(unique_celltypes, color3))

    # --- 核心逻辑：判断是否需要分割 ---
    tasks = [] 

    if split_col:
        print(f"Splitting data by column: '{split_col}'")
        if split_col not in adata.obs.columns:
            raise ValueError(f"Column '{split_col}' not found in adata.obs")
        
        for sample in adata.obs[split_col].unique():
            tmp = adata[adata.obs[split_col] == sample].copy()
            tasks.append((str(sample), tmp))
    else:
        print("No split column provided. Using all data as one group.")
        tasks.append(("All_Samples", adata.copy()))

    # --- 执行分析循环 ---
    for sample_name, tmp_adata in tasks:
        try:
            save_dir = os.path.join(Output, sample_name)
            os.makedirs(save_dir, exist_ok=True)
            
            if tmp_adata.n_obs < 100:
                print(f"Skipping {sample_name} - only {tmp_adata.n_obs} cells (less than 100)")
                continue
            
            print(f"Processing: {sample_name} (Cells: {tmp_adata.n_obs})")
            
            # 运行 CellPhoneDB
            cpdb_results, adata_cpdb = ov.single.run_cellphonedb_v5(
                tmp_adata,
                celltype_key=celltype_col,
                cpdb_file_path=DB_FILE_PATH, 
                output_dir=save_dir,
                temp_dir="./"
            )
            
            ov.utils.save(cpdb_results, os.path.join(save_dir, 'gex_cpdb_test.pkl'))
            adata_cpdb.write(os.path.join(save_dir, 'gex_cpdb_ad.h5ad'))
            
        except Exception as e:
            print(f"Error processing {sample_name}: {str(e)}")
            continue

    # --- 绘图部分 (已修复) ---
    if process_cellchat_analysis:
        print("Starting plotting process...")
        try:
            # === FIX START: 确保 group_col 始终有值 ===
            plot_group_col = split_col # 默认使用传入的参数
            
            if plot_group_col is None:
                print("Notice: No split column provided. Creating dummy column 'All_Analysis_Group' for plotting compatibility.")
                # 在 adata.obs 中新建一列，所有细胞的值都设为 'All_Samples'
                # 这样绘图函数就会把所有细胞当做同一个组来处理，而不会因为 None 报错
                adata.obs['All_Analysis_Group'] = 'All_Samples' 
                plot_group_col = 'All_Analysis_Group'
            # === FIX END ===

            process_cellchat_analysis(
                color_dict, 
                base_path=Output, 
            )
            print("Plotting finished.")
            
        except Exception as e:
            print(f"Error during plotting: {str(e)}")
            traceback.print_exc() # 打印完整报错堆栈，方便排查
    else:
        print("Skipping plotting because function could not be loaded.")

if __name__ == "__main__":
    main()