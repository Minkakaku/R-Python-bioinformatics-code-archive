import scanpy as sc
import pandas as pd
import numpy as np
import random

# 设置随机种子以保证结果可重复
np.random.seed(42)
random.seed(42)

# 1. 读取原始h5ad文件
print("读取原始h5ad文件...")
adata = sc.read_h5ad("MPC_samples.h5ad")

# 2. 检查数据结构
print(f"原始数据形状: {adata.shape}")
print(f"细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
print("\n原始数据列名:")
print(adata.obs.columns)

# 3. 查找细胞类型注释列
annotation_cols = [col for col in adata.obs.columns if any(keyword in col.lower() for keyword in ["cell", "type", "anno", "cluster", "label"])]

if len(annotation_cols) == 0:
    print("警告: 未找到细胞类型注释列，请检查列名。")
    print("可用的列名:")
    print(adata.obs.columns.tolist())
    cell_type_col = input("请输入细胞类型列名: ")
else:
    print(f"找到可能的细胞类型列: {annotation_cols}")
    if len(annotation_cols) > 1:
        for i, col in enumerate(annotation_cols):
            print(f"{i}: {col} (唯一值数量: {adata.obs[col].nunique()})")
        idx = int(input("请选择细胞类型列的编号: "))
        cell_type_col = annotation_cols[idx]
    else:
        cell_type_col = annotation_cols[0]

print(f"使用的细胞类型列: {cell_type_col}")

# 4. 查看原始细胞类型分布
original_cell_type_counts = adata.obs[cell_type_col].value_counts()
print(f"\n原始细胞类型分布 (共{len(original_cell_type_counts)}种类型):")
for i, (cell_type, count) in enumerate(original_cell_type_counts.items()):
    print(f"{i+1:2d}. {cell_type:20s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")

# 5. 去除Renal_Epithelium细胞
print(f"\n去除Renal_Epithelium细胞...")
original_renal_count = original_cell_type_counts.get("Renal_Epithelium", 0)
adata_no_renal = adata[adata.obs[cell_type_col] != "Renal_Epithelium", :].copy()

print(f"去除前细胞总数: {adata.shape[0]}")
print(f"去除后细胞总数: {adata_no_renal.shape[0]}")
print(f"去除Renal_Epithelium细胞数: {original_renal_count}")

# 6. 查看去除Renal_Epithelium后的细胞类型分布
no_renal_cell_type_counts = adata_no_renal.obs[cell_type_col].value_counts()
print(f"\n去除Renal_Epithelium后的细胞类型分布 (共{len(no_renal_cell_type_counts)}种类型):")
for i, (cell_type, count) in enumerate(no_renal_cell_type_counts.items()):
    print(f"{i+1:2d}. {cell_type:20s}: {count:6d} 个细胞 ({count/adata_no_renal.shape[0]:.2%})")

# 7. 计算去除Renal_Epithelium后的细胞类型比例
no_renal_total = adata_no_renal.shape[0]
target_total = 30000

if no_renal_total <= target_total:
    print(f"\n去除Renal_Epithelium后，数据只有 {no_renal_total} 个细胞，少于目标 {target_total} 个，将使用全部细胞。")
    adata_subset = adata_no_renal.copy()
else:
    print(f"\n基于去除Renal_Epithelium后的比例，计算分层抽样策略，目标抽取 {target_total} 个细胞...")
    
    # 计算去除Renal_Epithelium后每个细胞类型的比例
    cell_type_proportions = no_renal_cell_type_counts / no_renal_total
    
    # 计算每个细胞类型应抽取的细胞数（四舍五入）
    target_counts = (cell_type_proportions * target_total).round().astype(int)
    
    # 调整总数使其等于target_total
    total_adjust = target_total - target_counts.sum()
    
    if total_adjust != 0:
        print(f"初步计算总数: {target_counts.sum()}，需要调整 {total_adjust} 个细胞")
        
        # 如果总数不符合，调整细胞数量最多的几个细胞类型
        if total_adjust > 0:
            # 如果总数太少，增加细胞数量最多的几个细胞类型
            largest_types = no_renal_cell_type_counts.nlargest(abs(total_adjust)).index
            for cell_type in largest_types:
                target_counts[cell_type] += 1
                total_adjust -= 1
                if total_adjust == 0:
                    break
        else:
            # 如果总数太多，减少细胞数量最多的几个细胞类型
            largest_types = no_renal_cell_type_counts.nlargest(abs(total_adjust)).index
            for cell_type in largest_types:
                if target_counts[cell_type] > 1:  # 确保至少保留一个细胞
                    target_counts[cell_type] -= 1
                    total_adjust += 1
                    if total_adjust == 0:
                        break
    
    if total_adjust != 0:
        print(f"警告: 无法精确调整到 {target_total} 个细胞，最终总数: {target_counts.sum()}")
    
    print("\n去除Renal_Epithelium后，每个细胞类型的抽样策略:")
    for i, (cell_type, target_count) in enumerate(target_counts.items()):
        original_count = no_renal_cell_type_counts[cell_type]
        proportion = original_count / no_renal_total
        print(f"{i+1:2d}. {cell_type:20s}: 原始 {original_count:4d} -> 抽取 {target_count:4d} ({proportion:.2%})")
    
    # 8. 执行分层抽样
    sampled_indices = []
    
    for cell_type, n_samples in target_counts.items():
        # 获取该细胞类型的所有细胞索引
        cell_type_mask = adata_no_renal.obs[cell_type_col] == cell_type
        cell_type_indices = adata_no_renal.obs.index[cell_type_mask].tolist()
        
        if len(cell_type_indices) <= n_samples:
            # 如果细胞类型数量少于目标抽样数，抽取全部
            sampled_indices.extend(cell_type_indices)
            print(f"  {cell_type}: 只有 {len(cell_type_indices)} 个细胞，抽取全部")
        else:
            # 随机抽取指定数量的细胞
            sampled_cells = random.sample(cell_type_indices, n_samples)
            sampled_indices.extend(sampled_cells)
    
    # 9. 创建子集
    adata_subset = adata_no_renal[sampled_indices, :].copy()
    print(f"\n成功抽取 {len(sampled_indices)} 个细胞，子集形状: {adata_subset.shape}")

# 10. 查看子集中的细胞类型分布
subset_cell_type_counts = adata_subset.obs[cell_type_col].value_counts()
print(f"\n子集中的细胞类型分布:")
for i, (cell_type, count) in enumerate(subset_cell_type_counts.items()):
    print(f"{i+1:2d}. {cell_type:20s}: {count:6d} 个细胞 ({count/adata_subset.shape[0]:.2%})")

# 11. 比较原始数据、去除Renal_Epithelium后和子集的细胞类型分布
print(f"\n细胞类型分布比较:")
print(f"{'细胞类型':<20} {'原始数据':<12} {'去Renal后':<12} {'子集':<12} {'变化':<10}")
print("-" * 60)
for cell_type in original_cell_type_counts.index:
    if cell_type == "Renal_Epithelium":
        original_count = original_cell_type_counts.get(cell_type, 0)
        no_renal_count = 0
        subset_count = 0
        original_pct = original_count / adata.shape[0] * 100
        no_renal_pct = 0
        subset_pct = 0
        pct_change = subset_pct - original_pct
        print(f"{cell_type:<20} {original_count:<6} ({original_pct:5.2f}%) {no_renal_count:<6} ({no_renal_pct:5.2f}%) {subset_count:<6} ({subset_pct:5.2f}%) {'N/A':>10}")
    else:
        original_count = original_cell_type_counts.get(cell_type, 0)
        no_renal_count = no_renal_cell_type_counts.get(cell_type, 0)
        subset_count = subset_cell_type_counts.get(cell_type, 0)
        original_pct = original_count / adata.shape[0] * 100
        no_renal_pct = no_renal_count / no_renal_total * 100 if no_renal_total > 0 else 0
        subset_pct = subset_count / adata_subset.shape[0] * 100 if adata_subset.shape[0] > 0 else 0
        pct_change = subset_pct - original_pct
        print(f"{cell_type:<20} {original_count:<6} ({original_pct:5.2f}%) {no_renal_count:<6} ({no_renal_pct:5.2f}%) {subset_count:<6} ({subset_pct:5.2f}%) {pct_change:>+6.2f}%")

# 12. 保存子集到新的h5ad文件
output_file = "Epi_subset_no_renal.h5ad"
adata_subset.write_h5ad(output_file, compression='gzip')
print(f"\n子集已保存到: {output_file}")

# 13. 保存抽样信息到文本文件
with open("sampling_info_no_renal.txt", "w") as f:
    f.write("细胞类型分层抽样信息 (去除Renal_Epithelium后)\n")
    f.write("=" * 60 + "\n")
    f.write(f"原始数据细胞数: {adata.shape[0]}\n")
    f.write(f"去除Renal_Epithelium后细胞数: {no_renal_total}\n")
    f.write(f"目标抽取细胞数: {target_total}\n")
    f.write(f"实际抽取细胞数: {adata_subset.shape[0]}\n")
    f.write(f"使用的细胞类型列: {cell_type_col}\n\n")
    
    f.write("细胞类型分布比较:\n")
    f.write("-" * 60 + "\n")
    f.write(f"{'细胞类型':<20} {'原始数量':<10} {'原始比例':<10} {'去Renal后数量':<12} {'去Renal后比例':<12} {'抽取数量':<10} {'抽取比例':<10}\n")
    f.write("-" * 60 + "\n")
    
    for cell_type in original_cell_type_counts.index:
        if cell_type == "Renal_Epithelium":
            original_count = original_cell_type_counts.get(cell_type, 0)
            original_pct = original_count / adata.shape[0] * 100
            f.write(f"{cell_type:<20} {original_count:<10} {original_pct:9.2f}% {'0':<12} {'0.00%':<12} {'0':<10} {'0.00%':<10}\n")
        else:
            original_count = original_cell_type_counts.get(cell_type, 0)
            no_renal_count = no_renal_cell_type_counts.get(cell_type, 0)
            subset_count = subset_cell_type_counts.get(cell_type, 0)
            original_pct = original_count / adata.shape[0] * 100
            no_renal_pct = no_renal_count / no_renal_total * 100 if no_renal_total > 0 else 0
            subset_pct = subset_count / adata_subset.shape[0] * 100 if adata_subset.shape[0] > 0 else 0
            f.write(f"{cell_type:<20} {original_count:<10} {original_pct:9.2f}% {no_renal_count:<12} {no_renal_pct:9.2f}% {subset_count:<10} {subset_pct:9.2f}%\n")

print(f"抽样信息已保存到: sampling_info_no_renal.txt")

# 14. 可选的: 创建简单的可视化
try:
    import matplotlib.pyplot as plt
    
    # 创建包含三个子图的比较图
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # 原始数据
    axes[0].pie(original_cell_type_counts.values, labels=original_cell_type_counts.index, autopct='%1.1f%%')
    axes[0].set_title(f'原始数据 ({adata.shape[0]} 个细胞)')
    
    # 去除Renal_Epithelium后
    axes[1].pie(no_renal_cell_type_counts.values, labels=no_renal_cell_type_counts.index, autopct='%1.1f%%')
    axes[1].set_title(f'去除Renal后 ({no_renal_total} 个细胞)')
    
    # 子集数据
    axes[2].pie(subset_cell_type_counts.values, labels=subset_cell_type_counts.index, autopct='%1.1f%%')
    axes[2].set_title(f'抽取子集 ({adata_subset.shape[0]} 个细胞)')
    
    plt.tight_layout()
    plt.savefig('cell_type_distribution_comparison_no_renal.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"细胞类型分布比较图已保存到: cell_type_distribution_comparison_no_renal.png")
    
except ImportError:
    print("未安装matplotlib，跳过可视化")

print("\n数据处理完成！")
