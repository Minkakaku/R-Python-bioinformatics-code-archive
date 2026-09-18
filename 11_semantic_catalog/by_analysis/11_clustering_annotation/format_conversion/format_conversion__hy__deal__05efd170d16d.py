import scanpy as sc
import numpy as np
import pandas as pd

# 1. 读取源数据
print("读取源数据 adata_celltype.h5ad...")
adata_source = sc.read_h5ad('adata_celltype.h5ad')
print(f"源数据原始形状: {adata_source.shape}")

# 2. 将raw矩阵赋值给X矩阵
if adata_source.raw is not None:
    print("将raw矩阵赋值给X矩阵...")
    adata_source.X = adata_source.raw.X.copy()
    print(f"X矩阵替换完成，数据类型: {type(adata_source.X)}")
    
    # 如果需要，也可以更新var
    if adata_source.raw.var is not None:
        adata_source.var = adata_source.raw.var.copy()
        print("var已更新为raw.var")
else:
    print("警告: 源数据没有raw属性，保持X矩阵不变")

# 3. 读取目标数据
print("\n读取目标数据 Epi_new_nodup.h5ad...")
adata_target = sc.read_h5ad('Epi_new_nodup.h5ad')
print(f"目标数据形状: {adata_target.shape}")

# 4. 检查基因顺序
print("\n检查基因顺序...")
source_genes = list(adata_source.var_names)
target_genes = list(adata_target.var_names)

# 4.1 检查基因名称是否完全相同
if set(source_genes) == set(target_genes):
    if source_genes == target_genes:
        print("✅ 基因顺序完全相同")
        gene_match_type = "exact"
    else:
        print("⚠️ 基因集合相同但顺序不同")
        print("  将重排源数据的基因顺序以匹配目标数据")
        
        # 创建基因顺序映射
        gene_order_map = {}
        for i, gene in enumerate(target_genes):
            gene_order_map[gene] = i
        
        # 检查源数据中是否包含所有目标基因
        missing_genes = []
        for gene in target_genes:
            if gene not in source_genes:
                missing_genes.append(gene)
        
        if missing_genes:
            print(f"错误: 源数据缺少 {len(missing_genes)} 个目标数据中的基因:")
            for gene in missing_genes[:10]:  # 只显示前10个
                print(f"  - {gene}")
            if len(missing_genes) > 10:
                print(f"  ... 还有 {len(missing_genes)-10} 个")
            exit()
        
        # 获取重排索引
        reorder_indices = []
        for gene in target_genes:
            idx = source_genes.index(gene)
            reorder_indices.append(idx)
        
        print(f"基因重排索引计算完成")
        print(f"  源数据基因数: {len(source_genes)}")
        print(f"  目标数据基因数: {len(target_genes)}")
        print(f"  重排索引数: {len(reorder_indices)}")
        
        # 重排源数据的基因
        print("开始重排源数据的基因...")
        
        # 重排X矩阵
        if hasattr(adata_source.X, 'toarray'):
            # 如果是稀疏矩阵
            adata_source.X = adata_source.X[:, reorder_indices].copy()
        else:
            # 如果是密集矩阵
            adata_source.X = adata_source.X[:, reorder_indices]
        print("  X矩阵重排完成")
        
        # 重排var
        adata_source.var = adata_source.var.iloc[reorder_indices].copy()
        adata_source.var.index = target_genes
        print("  var重排完成")
        
        # 重排varm
        for key in list(adata_source.varm.keys()):
            if adata_source.varm[key].shape[0] == len(source_genes):
                adata_source.varm[key] = adata_source.varm[key][reorder_indices]
                print(f"  varm['{key}']重排完成")
        
        gene_match_type = "reordered"
        print("✅ 基因重排完成")
        
else:
    # 基因集合不同
    common_genes = set(source_genes) & set(target_genes)
    source_only_genes = set(source_genes) - set(target_genes)
    target_only_genes = set(target_genes) - set(source_genes)
    
    print(f"⚠️ 基因集合不完全相同:")
    print(f"  共同基因: {len(common_genes)} 个")
    print(f"  源数据特有基因: {len(source_only_genes)} 个")
    print(f"  目标数据特有基因: {len(target_only_genes)} 个")
    
    if len(common_genes) == 0:
        print("❌ 错误: 没有共同基因，无法合并")
        exit()
    
    print("\n处理策略:")
    print("  1. 只保留共同基因")
    print("  2. 按照目标数据的基因顺序排列")
    
    # 获取共同基因，按目标数据的顺序
    common_gene_order = [gene for gene in target_genes if gene in common_genes]
    
    # 获取源数据中共同基因的索引
    source_common_indices = []
    for gene in common_gene_order:
        if gene in source_genes:
            source_common_indices.append(source_genes.index(gene))
        else:
            print(f"错误: 基因 {gene} 应该在共同基因中但未在源数据中找到")
            exit()
    
    # 获取目标数据中共同基因的索引
    target_common_indices = []
    for gene in common_gene_order:
        if gene in target_genes:
            target_common_indices.append(target_genes.index(gene))
    
    print(f"\n开始处理源数据，只保留共同基因...")
    print(f"  原始形状: {adata_source.shape}")
    
    # 重排源数据的X矩阵
    if hasattr(adata_source.X, 'toarray'):
        adata_source.X = adata_source.X[:, source_common_indices].copy()
    else:
        adata_source.X = adata_source.X[:, source_common_indices]
    
    # 重排源数据的var
    adata_source.var = adata_source.var.iloc[source_common_indices].copy()
    adata_source.var.index = common_gene_order
    
    # 重排源数据的varm
    for key in list(adata_source.varm.keys()):
        if adata_source.varm[key].shape[0] == len(source_genes):
            adata_source.varm[key] = adata_source.varm[key][source_common_indices]
    
    print(f"  处理后的形状: {adata_source.shape}")
    
    print(f"\n开始处理目标数据，只保留共同基因...")
    print(f"  原始形状: {adata_target.shape}")
    
    # 重排目标数据的X矩阵
    if hasattr(adata_target.X, 'toarray'):
        adata_target.X = adata_target.X[:, target_common_indices].copy()
    else:
        adata_target.X = adata_target.X[:, target_common_indices]
    
    # 重排目标数据的var
    adata_target.var = adata_target.var.iloc[target_common_indices].copy()
    adata_target.var.index = common_gene_order
    
    print(f"  处理后的形状: {adata_target.shape}")
    
    gene_match_type = "intersection"
    print("✅ 基因交集处理完成")

# 5. 确保cell_id和sample列都存在
if 'cell_id' not in adata_source.obs.columns or 'sample' not in adata_source.obs.columns:
    print("❌ 错误: 源数据缺少'cell_id'或'sample'列")
    exit()

if 'cell_id' not in adata_target.obs.columns or 'sample' not in adata_target.obs.columns:
    print("❌ 错误: 目标数据缺少'cell_id'或'sample'列")
    exit()

# 6. 创建复合键用于细胞匹配
print("\n创建复合键用于细胞匹配...")
adata_source.obs['composite_key'] = adata_source.obs['cell_id'].astype(str) + '_' + adata_source.obs['sample'].astype(str)
adata_target.obs['composite_key'] = adata_target.obs['cell_id'].astype(str) + '_' + adata_target.obs['sample'].astype(str)

# 7. 检查匹配情况
source_cells = len(adata_source)
target_cells = len(adata_target)
matched_keys = set(adata_target.obs['composite_key']).intersection(set(adata_source.obs['composite_key']))
match_count = len(matched_keys)

print(f"源数据细胞数: {source_cells}")
print(f"目标数据细胞数: {target_cells}")
print(f"匹配的细胞数: {match_count}")
print(f"匹配率: {match_count/target_cells*100:.1f}%")

if match_count == 0:
    print("❌ 错误: 没有找到匹配的细胞")
    exit()

# 8. 创建映射字典
source_key_to_index = dict(zip(adata_source.obs['composite_key'], adata_source.obs.index))

# 9. 按照目标数据的顺序获取源数据中对应细胞的索引
ordered_source_indices = []
missing_cells = []

for i, key in enumerate(adata_target.obs['composite_key']):
    if key in source_key_to_index:
        ordered_source_indices.append(source_key_to_index[key])
    else:
        ordered_source_indices.append(None)
        missing_cells.append(i)

if missing_cells:
    print(f"⚠️ 警告: 有 {len(missing_cells)} 个细胞在源数据中找不到匹配")

# 10. 开始合并信息
print("\n开始合并信息...")

# 10.1 合并obs中的列
print("合并obs中的列...")
obs_cols_to_add = [col for col in adata_source.obs.columns 
                  if col not in ['cell_id', 'sample', 'composite_key'] 
                  and col not in adata_target.obs.columns]

for col in obs_cols_to_add:
    print(f"  添加列: {col}")
    adata_target.obs[col] = np.nan
    
    for i, (target_idx, source_idx) in enumerate(zip(adata_target.obs.index, ordered_source_indices)):
        if source_idx is not None:
            adata_target.obs.at[target_idx, col] = adata_source.obs.at[source_idx, col]

# 10.2 合并obsm中的内容
print("\n合并obsm中的内容...")
for key in adata_source.obsm.keys():
    if key not in adata_target.obsm.keys():
        print(f"  添加obsm: {key}")
        source_shape = adata_source.obsm[key].shape
        target_shape = (adata_target.shape[0], source_shape[1])
        new_array = np.full(target_shape, np.nan)
        
        for i, (target_idx, source_idx) in enumerate(zip(adata_target.obs.index, ordered_source_indices)):
            if source_idx is not None:
                source_loc = adata_source.obs.index.get_loc(source_idx)
                new_array[i] = adata_source.obsm[key][source_loc]
        
        adata_target.obsm[key] = new_array

# 10.3 合并varm中的内容
print("\n合并varm中的内容...")
for key in adata_source.varm.keys():
    if key not in adata_target.varm.keys():
        print(f"  添加varm: {key}")
        adata_target.varm[key] = adata_source.varm[key].copy()

# 10.4 合并obsp中的内容
print("\n合并obsp中的内容...")
for key in adata_source.obsp.keys():
    if key not in adata_target.obsp.keys():
        print(f"  添加obsp: {key}")
        n_cells_target = adata_target.shape[0]
        new_matrix = np.full((n_cells_target, n_cells_target), np.nan)
        
        target_indices_list = list(adata_target.obs.index)
        
        for i, (target_idx_i, source_idx_i) in enumerate(zip(target_indices_list, ordered_source_indices)):
            if source_idx_i is not None:
                source_loc_i = adata_source.obs.index.get_loc(source_idx_i)
                for j, (target_idx_j, source_idx_j) in enumerate(zip(target_indices_list, ordered_source_indices)):
                    if source_idx_j is not None:
                        source_loc_j = adata_source.obs.index.get_loc(source_idx_j)
                        new_matrix[i, j] = adata_source.obsp[key][source_loc_i, source_loc_j]
        
        adata_target.obsp[key] = new_matrix

# 10.5 合并uns中的内容
print("\n合并uns中的内容...")
for key in adata_source.uns.keys():
    if key not in adata_target.uns.keys():
        print(f"  添加uns: {key}")
        adata_target.uns[key] = adata_source.uns[key]

# 11. 清理临时列
if 'composite_key' in adata_target.obs.columns:
    adata_target.obs = adata_target.obs.drop(columns=['composite_key'])

# 12. 保存结果
output_file = 'Epi_new_nodup_merged.h5ad'
print(f"\n保存合并结果到 {output_file}...")
adata_target.write(output_file)
print(f"✅ 合并完成，已保存为 {output_file}")

# 13. 验证结果
print("\n验证合并结果:")
print(f"最终数据形状: {adata_target.shape}")
print(f"obs列数: {adata_target.obs.shape[1]}")
print(f"obsm键: {list(adata_target.obsm.keys())}")
print(f"varm键: {list(adata_target.varm.keys())}")
print(f"obsp键: {list(adata_target.obsp.keys())}")
print(f"uns键: {list(adata_target.uns.keys())[:10]}...")  # 只显示前10个

# 14. 保存合并日志
log_file = 'merge_log.txt'
with open(log_file, 'w') as f:
    f.write(f"数据合并日志\n")
    f.write(f"============\n")
    f.write(f"合并时间: {pd.Timestamp.now()}\n")
    f.write(f"源数据文件: adata_celltype.h5ad\n")
    f.write(f"目标数据文件: Epi_new_nodup.h5ad\n")
    f.write(f"输出文件: {output_file}\n")
    f.write(f"\n基因处理情况:\n")
    f.write(f"  处理类型: {gene_match_type}\n")
    f.write(f"  最终基因数: {adata_target.shape[1]}\n")
    f.write(f"\n细胞匹配情况:\n")
    f.write(f"  源数据细胞数: {source_cells}\n")
    f.write(f"  目标数据细胞数: {target_cells}\n")
    f.write(f"  匹配细胞数: {match_count}\n")
    f.write(f"  匹配率: {match_count/target_cells*100:.1f}%\n")
    f.write(f"  未匹配细胞数: {len(missing_cells)}\n")

print(f"合并日志已保存到 {log_file}")
