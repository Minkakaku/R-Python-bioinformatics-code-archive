import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

def add_state_to_h5ad_optimized(h5ad_path, state_csv_path, output_path=None):
    """
    使用pandas向量化操作优化版本，快速将state.csv中的group信息匹配到h5ad数据
    
    参数:
    h5ad_path: h5ad文件路径
    state_csv_path: state.csv文件路径
    output_path: 输出文件路径，如果为None则覆盖原文件
    """
    
    print("=" * 60)
    print("使用优化方法将state.csv中的group信息匹配到h5ad数据")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    print(f"   观测变量(obs)列: {list(adata.obs.columns)}")
    
    # 检查必需的列
    required_cols = ['cell_id', 'sample']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"h5ad.obs中缺少必需的列: {missing_cols}")
    
    # 2. 读取state.csv数据
    print(f"\n2. 读取state.csv数据: {state_csv_path}")
    state_df = pd.read_csv(state_csv_path)
    print(f"   state.csv形状: {state_df.shape}")
    print(f"   state.csv列: {list(state_df.columns)}")
    
    # 检查state.csv的列
    if not {'cell_id', 'sample', 'group'}.issubset(state_df.columns):
        raise ValueError("state.csv必须包含'cell_id', 'sample', 'group'列")
    
    # 显示state.csv的前几行
    print(f"\n   state.csv前5行:")
    print(state_df.head())
    
    # 3. 检查state.csv中是否有重复的(cell_id, sample)组合
    print(f"\n3. 检查数据...")
    state_duplicates = state_df.duplicated(subset=['cell_id', 'sample'], keep=False)
    if state_duplicates.any():
        print(f"   警告: state.csv中有 {state_duplicates.sum()} 个重复的(cell_id, sample)组合")
        # 显示重复的行
        duplicates = state_df[state_duplicates]
        print(f"   重复行示例:")
        print(duplicates.head())
        # 去重，保留第一个
        original_len = len(state_df)
        state_df = state_df.drop_duplicates(subset=['cell_id', 'sample'], keep='first')
        print(f"   已去重: 从 {original_len} 行减少到 {len(state_df)} 行")
    
    # 4. 使用字典映射进行匹配，避免行数变化
    print(f"\n4. 使用字典映射进行匹配...")
    
    # 创建映射字典
    state_mapping = {}
    for _, row in state_df.iterrows():
        key = (str(row['cell_id']), str(row['sample']))
        state_mapping[key] = row['group']
    
    print(f"   创建了 {len(state_mapping)} 个映射条目")
    
    # 5. 批量匹配
    print(f"\n5. 批量匹配group信息到h5ad...")
    
    # 创建state列，初始值为NaN
    adata.obs['state'] = np.nan
    
    # 准备数据
    h5ad_cell_ids = adata.obs['cell_id'].astype(str).tolist()
    h5ad_samples = adata.obs['sample'].astype(str).tolist()
    h5ad_indices = adata.obs.index.tolist()
    
    matched_count = 0
    unmatched_count = 0
    
    # 批量查找
    for idx, cell_id, sample in zip(h5ad_indices, h5ad_cell_ids, h5ad_samples):
        key = (cell_id, sample)
        if key in state_mapping:
            adata.obs.at[idx, 'state'] = state_mapping[key]
            matched_count += 1
        else:
            unmatched_count += 1
    
    # 6. 统计匹配结果
    print(f"\n6. 匹配结果统计:")
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   匹配成功: {matched_count} ({matched_count/adata.n_obs*100:.2f}%)")
    print(f"   未匹配: {unmatched_count} ({unmatched_count/adata.n_obs*100:.2f}%)")
    
    # 统计每个group的细胞数
    state_counts = adata.obs['state'].value_counts(dropna=False)
    print(f"\n   state列分布:")
    for state, count in state_counts.items():
        percentage = count / adata.n_obs * 100
        state_display = 'NaN' if pd.isna(state) else state
        print(f"     {state_display}: {count} ({percentage:.2f}%)")
    
    # 7. 保存结果
    if output_path is None:
        output_path = h5ad_path
        print(f"\n7. 正在覆盖原始文件: {output_path}")
    else:
        print(f"\n7. 正在保存到新文件: {output_path}")
    
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    # 8. 创建匹配结果报告
    report_path = os.path.splitext(output_path)[0] + "_matching_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("State信息匹配报告（优化版本）\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("文件信息:\n")
        f.write(f"  h5ad文件: {h5ad_path}\n")
        f.write(f"  state.csv文件: {state_csv_path}\n")
        f.write(f"  输出文件: {output_path}\n\n")
        
        f.write("数据统计:\n")
        f.write(f"  h5ad总细胞数: {adata.n_obs}\n")
        f.write(f"  h5ad总基因数: {adata.n_vars}\n")
        f.write(f"  state.csv记录数: {len(state_df)}\n")
        f.write(f"  state.csv唯一cell_id数: {state_df['cell_id'].nunique()}\n")
        f.write(f"  state.csv唯一sample数: {state_df['sample'].nunique()}\n\n")
        
        f.write("匹配结果:\n")
        f.write(f"  匹配成功: {matched_count} ({matched_count/adata.n_obs*100:.2f}%)\n")
        f.write(f"  未匹配: {unmatched_count} ({unmatched_count/adata.n_obs*100:.2f}%)\n\n")
        
        f.write("state列分布:\n")
        for state, count in state_counts.items():
            percentage = count / adata.n_obs * 100
            state_display = 'NaN' if pd.isna(state) else state
            f.write(f"  {state_display}: {count} ({percentage:.2f}%)\n")
        
        f.write("\nstate.csv中group分布:\n")
        group_counts = state_df['group'].value_counts()
        for group, count in group_counts.items():
            percentage = count / len(state_df) * 100
            f.write(f"  {group}: {count} ({percentage:.2f}%)\n")
        
        # 添加未匹配的细胞信息
        f.write("\n未匹配的细胞样本分布（前20个）:\n")
        unmatched_samples = adata.obs[adata.obs['state'].isna()]['sample'].value_counts().head(20)
        for sample, count in unmatched_samples.items():
            f.write(f"  {sample}: {count} 个细胞\n")
    
    print(f"\n8. 匹配报告已保存到: {report_path}")
    
    return adata

def add_state_to_h5ad_ultra_fast(h5ad_path, state_csv_path, output_path=None):
    """
    使用pandas合并的极速版本，但更安全的实现
    """
    
    print("=" * 60)
    print("使用极速方法将state.csv中的group信息匹配到h5ad数据")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    
    # 检查必需的列
    required_cols = ['cell_id', 'sample']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"h5ad.obs中缺少必需的列: {missing_cols}")
    
    # 2. 读取state.csv数据
    print(f"\n2. 读取state.csv数据: {state_csv_path}")
    state_df = pd.read_csv(state_csv_path)
    print(f"   state.csv形状: {state_df.shape}")
    
    # 检查state.csv的列
    if not {'cell_id', 'sample', 'group'}.issubset(state_df.columns):
        raise ValueError("state.csv必须包含'cell_id', 'sample', 'group'列")
    
    # 3. 使用pandas合并数据 - 更安全的实现
    print(f"\n3. 使用pandas合并数据（安全实现）...")
    
    # 创建h5ad.obs的DataFrame副本
    h5ad_obs = adata.obs.copy()
    
    # 确保列名一致
    h5ad_obs['cell_id'] = h5ad_obs['cell_id'].astype(str)
    h5ad_obs['sample'] = h5ad_obs['sample'].astype(str)
    
    state_df['cell_id'] = state_df['cell_id'].astype(str)
    state_df['sample'] = state_df['sample'].astype(str)
    
    # 检查重复
    duplicates = state_df.duplicated(subset=['cell_id', 'sample'], keep=False)
    if duplicates.any():
        print(f"   警告: state.csv中有 {duplicates.sum()} 个重复的(cell_id, sample)组合，将去重")
        state_df = state_df.drop_duplicates(subset=['cell_id', 'sample'], keep='first')
    
    # 创建唯一的行标识符
    h5ad_obs['row_id'] = range(len(h5ad_obs))
    
    # 创建复合键
    h5ad_obs['composite_key'] = h5ad_obs['cell_id'] + '|||' + h5ad_obs['sample']
    state_df['composite_key'] = state_df['cell_id'] + '|||' + state_df['sample']
    
    # 创建映射字典
    state_mapping = state_df.set_index('composite_key')['group'].to_dict()
    
    # 使用map进行匹配
    h5ad_obs['state'] = h5ad_obs['composite_key'].map(state_mapping)
    
    # 恢复顺序
    h5ad_obs = h5ad_obs.sort_values('row_id')
    
    # 4. 统计匹配结果
    print(f"\n4. 匹配结果统计:")
    
    matched_count = h5ad_obs['state'].notna().sum()
    unmatched_count = h5ad_obs['state'].isna().sum()
    
    print(f"   总细胞数: {len(h5ad_obs)}")
    print(f"   匹配成功: {matched_count} ({matched_count/len(h5ad_obs)*100:.2f}%)")
    print(f"   未匹配: {unmatched_count} ({unmatched_count/len(h5ad_obs)*100:.2f}%)")
    
    # 5. 将group信息添加回h5ad
    print(f"\n5. 将group信息添加回h5ad...")
    
    adata.obs['state'] = h5ad_obs['state'].values
    
    # 统计每个group的细胞数
    state_counts = adata.obs['state'].value_counts(dropna=False)
    print(f"\n   state列分布:")
    for state, count in state_counts.items():
        percentage = count / len(h5ad_obs) * 100
        state_display = 'NaN' if pd.isna(state) else state
        print(f"     {state_display}: {count} ({percentage:.2f}%)")
    
    # 6. 保存结果
    if output_path is None:
        output_path = h5ad_path
        print(f"\n6. 正在覆盖原始文件: {output_path}")
    else:
        print(f"\n6. 正在保存到新文件: {output_path}")
    
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    return adata

def add_state_to_h5ad_chunked(h5ad_path, state_csv_path, output_path=None, chunk_size=100000):
    """
    分块处理方法，适合超大内存数据集
    """
    
    print("=" * 60)
    print("使用分块处理方法将state.csv中的group信息匹配到h5ad数据")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    
    # 检查必需的列
    required_cols = ['cell_id', 'sample']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"h5ad.obs中缺少必需的列: {missing_cols}")
    
    # 2. 读取state.csv数据
    print(f"\n2. 读取state.csv数据: {state_csv_path}")
    state_df = pd.read_csv(state_csv_path)
    print(f"   state.csv形状: {state_df.shape}")
    
    # 检查state.csv的列
    if not {'cell_id', 'sample', 'group'}.issubset(state_df.columns):
        raise ValueError("state.csv必须包含'cell_id', 'sample', 'group'列")
    
    # 创建映射字典
    state_mapping = {}
    for _, row in state_df.iterrows():
        key = (str(row['cell_id']), str(row['sample']))
        state_mapping[key] = row['group']
    
    print(f"   创建了 {len(state_mapping)} 个映射条目")
    
    # 3. 分块处理
    print(f"\n3. 分块处理数据 (块大小: {chunk_size})...")
    
    # 创建state列
    adata.obs['state'] = np.nan
    
    total_cells = adata.n_obs
    matched_count = 0
    unmatched_count = 0
    
    # 分块处理
    for start in range(0, total_cells, chunk_size):
        end = min(start + chunk_size, total_cells)
        chunk_size_current = end - start
        
        print(f"   处理块 {start//chunk_size + 1}/{(total_cells + chunk_size - 1)//chunk_size} "
              f"({start} - {end-1}, 共 {chunk_size_current} 个细胞)")
        
        # 获取当前块的数据
        chunk_indices = adata.obs.index[start:end]
        
        for idx in chunk_indices:
            row = adata.obs.loc[idx]
            cell_id = str(row['cell_id'])
            sample = str(row['sample'])
            key = (cell_id, sample)
            
            if key in state_mapping:
                adata.obs.at[idx, 'state'] = state_mapping[key]
                matched_count += 1
            else:
                unmatched_count += 1
    
    # 4. 统计匹配结果
    print(f"\n4. 匹配结果统计:")
    print(f"   总细胞数: {total_cells}")
    print(f"   匹配成功: {matched_count} ({matched_count/total_cells*100:.2f}%)")
    print(f"   未匹配: {unmatched_count} ({unmatched_count/total_cells*100:.2f}%)")
    
    # 5. 保存结果
    if output_path is None:
        output_path = h5ad_path
        print(f"\n5. 正在覆盖原始文件: {output_path}")
    else:
        print(f"\n5. 正在保存到新文件: {output_path}")
    
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    return adata

if __name__ == "__main__":
    # 设置文件路径
    h5ad_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/all_filter_modified.h5ad"
    state_csv_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/state.csv"
    
    # 输出路径
    output_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_with_state_new.h5ad"
    
    print("选择匹配方法:")
    print("1. 优化方法（推荐，使用字典映射）")
    print("2. 极速方法（使用pandas合并，适合数据干净的情况）")
    print("3. 分块处理方法（适合超大内存数据集）")
    
    choice = input("\n请选择 (1/2/3, 默认1): ").strip() or "1"
    
    if choice == "1":
        # 使用方法1：优化方法
        adata = add_state_to_h5ad_optimized(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path,
            output_path=output_path
        )
        
    elif choice == "2":
        # 使用方法2：极速方法
        adata = add_state_to_h5ad_ultra_fast(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path,
            output_path=output_path
        )
        
    elif choice == "3":
        # 使用方法3：分块处理
        chunk_size = int(input("请输入块大小 (默认100000): ").strip() or "100000")
        adata = add_state_to_h5ad_chunked(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path,
            output_path=output_path,
            chunk_size=chunk_size
        )
    
    print("\n" + "=" * 60)
    print("处理完成!")
    print("=" * 60)
