import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import sys
import os

def add_CNs10_with_options(h5ad_path, xls_path, output_path=None, fill_na='keep', sep='\t'):
    """
    添加CNs_10信息到h5ad文件，包含更多选项
    
    参数:
    h5ad_path: h5ad文件路径
    xls_path: XLS特征文件路径
    output_path: 输出文件路径
    fill_na: 如何处理未匹配的细胞
             'keep': 保持NaN
             'unknown': 填充为'unknown'
             'drop': 删除未匹配的细胞
    sep: XLS文件分隔符
    """
    
    # 读取h5ad文件
    print(f"读取h5ad文件: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列
    required_cols = ['cell_id', 'sample']
    for col in required_cols:
        if col not in adata.obs.columns:
            print(f"错误: h5ad文件中缺少'{col}'列")
            print(f"可用的列: {list(adata.obs.columns)}")
            return None
    
    # 读取XLS文件
    print(f"读取XLS文件: {xls_path}")
    try:
        # 尝试不同的读取方式
        if xls_path.endswith('.xls') or xls_path.endswith('.xlsx'):
            features_df = pd.read_excel(xls_path)
        else:
            # 尝试自动检测分隔符
            features_df = pd.read_csv(xls_path, sep=sep, engine='python')
    except Exception as e:
        print(f"读取XLS文件失败: {e}")
        return None
    
    print(f"XLS文件: {len(features_df)} 行, {len(features_df.columns)} 列")
    
    # 检查必要列
    if 'CNs_10' not in features_df.columns:
        print(f"错误: XLS文件中缺少'CNs_10'列")
        print(f"XLS文件列名: {list(features_df.columns)}")
        return None
    
    # 确保列名一致
    cell_id_col = 'cell_id' if 'cell_id' in features_df.columns else \
                 next((c for c in features_df.columns if 'cell_id' in c.lower()), None)
    sample_col = 'sample' if 'sample' in features_df.columns else \
                next((c for c in features_df.columns if 'sample' in c.lower()), None)
    
    if not cell_id_col or not sample_col:
        print("错误: 在XLS文件中找不到cell_id或sample列")
        return None
    
    # 重命名列以便统一
    features_df = features_df.rename(columns={cell_id_col: 'cell_id', sample_col: 'sample'})
    
    # 去重
    duplicate_mask = features_df.duplicated(subset=['cell_id', 'sample'], keep=False)
    if duplicate_mask.any():
        print(f"警告: 发现 {duplicate_mask.sum()} 个重复的cell_id-sample组合")
        features_df = features_df.drop_duplicates(subset=['cell_id', 'sample'], keep='first')
    
    # 创建合并键
    adata.obs['merge_key'] = adata.obs['cell_id'].astype(str) + '_' + adata.obs['sample'].astype(str)
    features_df['merge_key'] = features_df['cell_id'].astype(str) + '_' + features_df['sample'].astype(str)
    
    # 创建CNs_10值的映射字典
    cns_dict = dict(zip(features_df['merge_key'], features_df['CNs_10']))
    
    # 添加CNs_10列
    adata.obs['CNs_10'] = adata.obs['merge_key'].map(cns_dict)
    
    # 处理未匹配的细胞
    unmatched_count = adata.obs['CNs_10'].isna().sum()
    print(f"\n匹配统计:")
    print(f"  总细胞数: {adata.n_obs}")
    print(f"  成功匹配: {adata.n_obs - unmatched_count}")
    print(f"  未匹配: {unmatched_count}")
    
    if fill_na == 'unknown':
        adata.obs['CNs_10'] = adata.obs['CNs_10'].fillna('unknown')
        print("  未匹配的细胞填充为'unknown'")
    elif fill_na == 'drop':
        original_size = adata.n_obs
        adata = adata[adata.obs['CNs_10'].notna(), :].copy()
        print(f"  删除未匹配细胞: {original_size} -> {adata.n_obs}")
    
    # 删除临时列
    adata.obs.drop('merge_key', axis=1, inplace=True, errors='ignore')
    
    # CNs_10值统计
    print(f"\nCNs_10值统计:")
    print(adata.obs['CNs_10'].describe())
    
    # 保存结果
    if output_path:
        adata.write(output_path)
        print(f"\n结果保存到: {output_path}")
    
    return adata

def main_cli():
    """命令行界面"""
    parser = argparse.ArgumentParser(description='将CNs_10信息添加到h5ad文件')
    parser.add_argument('--h5ad', '-i', required=True, help='输入h5ad文件路径')
    parser.add_argument('--xls', '-x', required=True, help='XLS特征文件路径')
    parser.add_argument('--output', '-o', default=None, help='输出h5ad文件路径')
    parser.add_argument('--fill-na', choices=['keep', 'unknown', 'drop'], 
                       default='keep', help='如何处理未匹配的细胞')
    parser.add_argument('--sep', default='\t', help='XLS文件分隔符')
    
    args = parser.parse_args()
    
    # 检查文件是否存在
    for file_path, desc in [(args.h5ad, 'h5ad文件'), (args.xls, 'XLS文件')]:
        if not os.path.exists(file_path):
            print(f"错误: 找不到{desc}: {file_path}")
            sys.exit(1)
    
    print("=" * 60)
    print("CNs_10信息添加工具")
    print("=" * 60)
    
    adata = add_CNs10_with_options(
        h5ad_path=args.h5ad,
        xls_path=args.xls,
        output_path=args.output,
        fill_na=args.fill_na,
        sep=args.sep
    )
    
    if adata is not None:
        print("\n处理完成!")
        if args.output:
            print(f"结果已保存到: {args.output}")

if __name__ == "__main__":
    # 直接运行
    # main()
    
    # 或使用命令行界面
    main_cli()
