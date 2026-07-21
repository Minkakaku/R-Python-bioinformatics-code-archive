import scanpy as sc
import pandas as pd
import numpy as np
import sys
from pathlib import Path

def main():
    # 文件路径
    input_file = "all_raw.h5ad"  # 请替换为您的输入文件名
    output_file = "all_filter.h5ad"
    
    # 检查输入文件是否存在
    if not Path(input_file).exists():
        print(f"错误: 文件 {input_file} 不存在！")
        sys.exit(1)
    
    print("1. 读取h5ad数据...")
    adata = sc.read_h5ad(input_file)
    print(f"   原始数据形状: {adata.shape}")
    print(f"   细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
    
    # 检查是否有celltype列
    if 'celltype' not in adata.obs.columns:
        print("错误: 数据中缺少 'celltype' 列")
        print(f"可用列: {adata.obs.columns.tolist()}")
        sys.exit(1)
    
    # 2. 修改celltype列
    print("\n2. 修改celltype列...")
    print(f"   原始celltype类型: {adata.obs['celltype'].nunique()} 种")
    print(f"   原始celltype分布:")
    celltype_counts = adata.obs['celltype'].value_counts()
    for i, (ctype, count) in enumerate(celltype_counts.items()):
        print(f"     {i+1:2d}. {ctype:30s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    # 定义需要修改的映射关系
    celltype_mapping = {
        "M1-like Inflammatory": "CTSL+ PLAUR+ Macro",
        "M2-like": "STAB1+ Macro",
        "B_GC_Proliferating": "GCB",
        "B_Memory_like": "MBC",
        "B_Mast": "Mast",
        "Angiogenic":"Angiogenic EC",
        "Arteriole":"Arteriole EC",
        "Lymphatic":"Lymphatic EC",
        "Venule":"Venule EC"
    }
    
    # 备份原始celltype
    #adata.obs['celltype_original'] = adata.obs['celltype'].copy()
    
    # 应用映射
    adata.obs['celltype'] = adata.obs['celltype'].map(lambda x: celltype_mapping.get(x, x))
    
    print(f"   修改后celltype类型: {adata.obs['celltype'].nunique()} 种")
    print(f"   修改后的celltype分布:")
    celltype_counts_after = adata.obs['celltype'].value_counts()
    for i, (ctype, count) in enumerate(celltype_counts_after.items()):
        print(f"     {i+1:2d}. {ctype:30s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    # 3. 提取指定的细胞类型
    print("\n3. 提取指定的细胞类型...")
    
    # 定义要保留的细胞类型列表
    keep_celltypes = [
        "CXCL5+ Inf-stem-tumor",
        "Prol-Metab-Epi",
        "EMT-like Inv Ent",
        "Tumor Epi",
        "WNT+ ISC",
        "MUC5B+ CGC",
        "Hypo-Stem-tumor",
        "Hypo-Metab",
        "Inf-Hyp Epi",
        "Early stress",
        "Metab stress",
        "MUC5AC+ MGC",
        "Mat Ent",
        "PCSK1+ AGC",
        "IGF2+ prolif-stem-tumor",
        "CTSL+ PLAUR+ Macro",  # 从M1-like映射而来
        "STAB1+ Macro",       # 从M2-like映射而来
        "CD4+ T",
        "γδ T",
        "CD8+ T",
        "Treg",
        "MKI67+ T",
        "DC",
        "GCB",                # 从B_GC_proliferating映射而来
        "MBC",                # 从B_Memory_like映射而来
        "Mast",               # 从B_Mast映射而来
        "Plasmas",
        "ADAMDEC1+ Homeo-Immuno Fib",
        "Act Strom CAF",
        "matri-prod CAF",
        "RSPO3+ GREM1+ Crypt Fib",
        "WNT5B+ VEGFA+ Villus Fib",
        "Myofib",
        "apCAF",
        "Com Fib",
        "CXCL12+ Fib",
        "Pericyte-like Fib",
        "Neutrophils",
        "Angiogenic EC", 
        "Arteriole EC",
        "Lymphatic EC",
        "Venule EC"
    ]
    
    # 检查要保留的细胞类型是否存在
    existing_keep_types = [ctype for ctype in keep_celltypes if ctype in adata.obs['celltype'].unique()]
    non_existing_types = [ctype for ctype in keep_celltypes if ctype not in adata.obs['celltype'].unique()]
    
    if existing_keep_types:
        print(f"   将保留以下细胞类型: {existing_keep_types}")
        
        # 统计要保留的细胞数量
        keep_mask = adata.obs['celltype'].isin(existing_keep_types)
        cells_to_keep = keep_mask.sum()
        
        print(f"   将保留的细胞数: {cells_to_keep} ({cells_to_keep/adata.shape[0]:.2%})")
        print(f"   保留前细胞数: {adata.shape[0]}")
        
        # 提取这些细胞
        adata_filtered = adata[keep_mask, :].copy()
        
        print(f"   保留后细胞数: {adata_filtered.shape[0]}")
        print(f"   保留的细胞比例: {adata_filtered.shape[0]/adata.shape[0]:.2%}")
        
        # 统计每个被保留类型的细胞数
        print(f"\n   被保留的细胞类型分布:")
        for ctype in existing_keep_types:
            count = (adata.obs['celltype'] == ctype).sum()
            if count > 0:
                print(f"     {ctype:30s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    if non_existing_types:
        print(f"\n   注意: 以下指定的保留类型在数据中不存在: {non_existing_types}")
    
    # 如果没有要保留的类型，则使用原始数据
    if not existing_keep_types:
        print("   没有找到需要保留的细胞类型，使用原始数据")
        adata_filtered = adata.copy()
    
    # 4. 检查过滤后的数据
    print(f"\n4. 过滤后数据统计:")
    print(f"   细胞数: {adata_filtered.shape[0]}, 基因数: {adata_filtered.shape[1]}")
    print(f"   celltype类型: {adata_filtered.obs['celltype'].nunique()} 种")
    print(f"   过滤后celltype分布:")
    
    filtered_celltype_counts = adata_filtered.obs['celltype'].value_counts()
    for i, (ctype, count) in enumerate(filtered_celltype_counts.items()):
        print(f"     {i+1:2d}. {ctype:30s}: {count:6d} 个细胞 ({count/adata_filtered.shape[0]:.2%})")
    
    # 5. 保存过滤后的数据
    print(f"\n5. 保存数据到 {output_file}...")
    adata_filtered.write_h5ad(output_file, compression='gzip')
    print(f"   文件已保存")
    
    # 6. 保存处理信息
    info_file = "filtering_info.txt"
    with open(info_file, 'w') as f:
        f.write("细胞过滤信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"输入文件: {input_file}\n")
        f.write(f"输出文件: {output_file}\n")
        f.write(f"原始细胞数: {adata.shape[0]}\n")
        f.write(f"过滤后细胞数: {adata_filtered.shape[0]}\n")
        f.write(f"剔除细胞数: {adata.shape[0] - adata_filtered.shape[0]}\n")
        f.write(f"保留比例: {adata_filtered.shape[0]/adata.shape[0]:.2%}\n\n")
        
        f.write("celltype修改映射:\n")
        for old, new in celltype_mapping.items():
            f.write(f"  {old} -> {new}\n")
        
        f.write("\n保留的细胞类型:\n")
        for ctype in keep_celltypes:
            f.write(f"  {ctype}\n")
        
        f.write("\n原始celltype分布:\n")
        for ctype, count in celltype_counts.items():
            f.write(f"  {ctype}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n修改后celltype分布:\n")
        for ctype, count in celltype_counts_after.items():
            f.write(f"  {ctype}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n过滤后celltype分布:\n")
        for ctype, count in filtered_celltype_counts.items():
            f.write(f"  {ctype}: {count} 个细胞 ({count/adata_filtered.shape[0]:.2%})\n")
        
        f.write("\n被保留的细胞类型统计:\n")
        for ctype in existing_keep_types:
            count = (adata.obs['celltype'] == ctype).sum()
            if count > 0:
                f.write(f"  {ctype}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
    
    print(f"\n6. 处理信息已保存到: {info_file}")
    
    # 7. 创建细胞类型统计表
    stats_df = pd.DataFrame({
        'celltype': filtered_celltype_counts.index,
        'count': filtered_celltype_counts.values,
        'percentage': (filtered_celltype_counts.values / adata_filtered.shape[0] * 100).round(2)
    })
    stats_df = stats_df.sort_values('count', ascending=False)
    stats_df.to_csv("filtered_celltype_statistics.csv", index=False, encoding='utf-8-sig')
    print(f"   细胞类型统计已保存到: filtered_celltype_statistics.csv")
    
    print("\n处理完成!")

if __name__ == "__main__":
    main()
