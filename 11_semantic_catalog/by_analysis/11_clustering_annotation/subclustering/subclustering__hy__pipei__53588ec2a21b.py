import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

def add_state_to_h5ad(h5ad_path, state_csv_path, output_path=None):
    """
    将state.csv中的group信息匹配到h5ad数据的obs中
    
    参数:
    h5ad_path: h5ad文件路径
    state_csv_path: state.csv文件路径
    output_path: 输出文件路径，如果为None则覆盖原文件
    """
    
    print("=" * 60)
    print("开始将state.csv中的group信息匹配到h5ad数据")
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
    
    # 3. 创建映射字典
    print(f"\n3. 创建映射字典...")
    # 使用(cell_id, sample)作为复合键
    state_mapping = {}
    for _, row in state_df.iterrows():
        key = (str(row['cell_id']), str(row['sample']))
        state_mapping[key] = row['group']
    
    print(f"   创建了 {len(state_mapping)} 个映射条目")
    
    # 4. 匹配group信息到h5ad
    print(f"\n4. 匹配group信息到h5ad...")
    
    # 创建state列，初始值为NaN
    adata.obs['state'] = np.nan
    
    matched_count = 0
    unmatched_count = 0
    
    # 遍历h5ad中的每个细胞
    for idx, row in adata.obs.iterrows():
        cell_id = str(row['cell_id'])
        sample = str(row['sample'])
        key = (cell_id, sample)
        
        if key in state_mapping:
            adata.obs.at[idx, 'state'] = state_mapping[key]
            matched_count += 1
        else:
            unmatched_count += 1
    
    # 5. 统计匹配结果
    print(f"\n5. 匹配结果统计:")
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   匹配成功: {matched_count} ({matched_count/adata.n_obs*100:.2f}%)")
    print(f"   未匹配: {unmatched_count} ({unmatched_count/adata.n_obs*100:.2f}%)")
    
    # 统计每个group的细胞数
    state_counts = adata.obs['state'].value_counts(dropna=False)
    print(f"\n   state列分布:")
    for state, count in state_counts.items():
        percentage = count / adata.n_obs * 100
        print(f"     {state}: {count} ({percentage:.2f}%)")
    
    # 6. 保存结果
    if output_path is None:
        # 如果没有指定输出路径，修改原始文件
        output_path = h5ad_path
        print(f"\n6. 正在覆盖原始文件: {output_path}")
    else:
        print(f"\n6. 正在保存到新文件: {output_path}")
    
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    # 7. 创建匹配结果报告
    report_path = os.path.splitext(output_path)[0] + "_matching_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("State信息匹配报告\n")
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
            f.write(f"  {state}: {count} ({percentage:.2f}%)\n")
        
        f.write("\nstate.csv中group分布:\n")
        group_counts = state_df['group'].value_counts()
        for group, count in group_counts.items():
            percentage = count / len(state_df) * 100
            f.write(f"  {group}: {count} ({percentage:.2f}%)\n")
    
    print(f"\n7. 匹配报告已保存到: {report_path}")
    
    return adata

# 使用pandas优化版本（处理大型数据集更快）
def add_state_to_h5ad_pandas(h5ad_path, state_csv_path, output_path=None):
    """
    使用pandas向量化操作，处理大型数据集更快
    """
    
    print("=" * 60)
    print("使用pandas优化版本匹配state信息")
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
    
    # 3. 使用pandas合并数据
    print(f"\n3. 使用pandas合并数据...")
    
    # 创建h5ad.obs的临时副本
    h5ad_obs = adata.obs.copy()
    
    # 重置索引，确保有唯一的标识符
    h5ad_obs['original_index'] = h5ad_obs.index
    
    # 在state_df中创建复合键
    state_df['composite_key'] = state_df['cell_id'].astype(str) + '_' + state_df['sample'].astype(str)
    state_mapping = state_df.set_index('composite_key')['group'].to_dict()
    
    # 在h5ad_obs中创建相同的复合键
    h5ad_obs['composite_key'] = h5ad_obs['cell_id'].astype(str) + '_' + h5ad_obs['sample'].astype(str)
    
    # 使用map函数进行匹配
    h5ad_obs['state'] = h5ad_obs['composite_key'].map(state_mapping)
    
    # 恢复原始索引
    h5ad_obs = h5ad_obs.set_index('original_index')
    
    # 更新adata.obs
    adata.obs['state'] = h5ad_obs['state']
    
    # 4. 统计匹配结果
    print(f"\n4. 匹配结果统计:")
    matched_count = adata.obs['state'].notna().sum()
    unmatched_count = adata.obs['state'].isna().sum()
    
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   匹配成功: {matched_count} ({matched_count/adata.n_obs*100:.2f}%)")
    print(f"   未匹配: {unmatched_count} ({unmatched_count/adata.n_obs*100:.2f}%)")
    
    # 5. 保存结果
    if output_path is None:
        output_path = h5ad_path
    
    adata.write_h5ad(output_path)
    print(f"\n5. 数据已保存到: {output_path}")
    
    return adata

# 批量处理多个h5ad文件
def add_state_to_multiple_h5ad(state_csv_path, h5ad_files, output_dir=None):
    """
    批量处理多个h5ad文件
    """
    
    print("=" * 60)
    print("批量处理多个h5ad文件")
    print("=" * 60)
    
    # 读取state.csv
    state_df = pd.read_csv(state_csv_path)
    
    # 按sample分组state.csv
    state_by_sample = {}
    for sample in state_df['sample'].unique():
        sample_state = state_df[state_df['sample'] == sample]
        state_by_sample[sample] = sample_state.set_index('cell_id')['group'].to_dict()
    
    print(f"state.csv包含 {len(state_by_sample)} 个样本的数据")
    
    for h5ad_path in h5ad_files:
        print(f"\n处理文件: {h5ad_path}")
        
        try:
            adata = sc.read_h5ad(h5ad_path)
            
            # 检查是否有sample列
            if 'sample' not in adata.obs.columns:
                print(f"  警告: 文件 {h5ad_path} 没有sample列，跳过")
                continue
            
            # 获取当前样本
            sample_name = adata.obs['sample'].iloc[0] if adata.obs['sample'].nunique() == 1 else None
            
            if sample_name is None:
                print(f"  警告: 文件 {h5ad_path} 包含多个样本，跳过")
                continue
            
            if sample_name not in state_by_sample:
                print(f"  警告: 样本 {sample_name} 不在state.csv中，跳过")
                continue
            
            # 获取该样本的映射
            sample_mapping = state_by_sample[sample_name]
            
            # 匹配state信息
            adata.obs['state'] = np.nan
            
            matched = 0
            for idx, row in adata.obs.iterrows():
                cell_id = str(row.get('cell_id', ''))
                if cell_id and cell_id in sample_mapping:
                    adata.obs.at[idx, 'state'] = sample_mapping[cell_id]
                    matched += 1
            
            # 保存结果
            if output_dir:
                filename = os.path.basename(h5ad_path)
                output_path = os.path.join(output_dir, f"with_state_{filename}")
            else:
                output_path = h5ad_path
            
            adata.write_h5ad(output_path)
            print(f"  处理完成: 匹配了 {matched}/{adata.n_obs} 个细胞")
            
        except Exception as e:
            print(f"  处理文件 {h5ad_path} 时出错: {e}")

# 验证匹配结果的函数
def verify_matching(h5ad_path, state_csv_path):
    """
    验证匹配结果的准确性
    """
    
    print("=" * 60)
    print("验证匹配结果")
    print("=" * 60)
    
    # 读取数据
    adata = sc.read_h5ad(h5ad_path)
    state_df = pd.read_csv(state_csv_path)
    
    # 检查是否有state列
    if 'state' not in adata.obs.columns:
        print("错误: h5ad文件中没有'state'列")
        return
    
    # 创建state.csv的映射
    state_mapping = {}
    for _, row in state_df.iterrows():
        key = (str(row['cell_id']), str(row['sample']))
        state_mapping[key] = row['group']
    
    # 验证匹配
    correct = 0
    incorrect = 0
    missing = 0
    
    for idx, row in adata.obs.iterrows():
        cell_id = str(row.get('cell_id', ''))
        sample = str(row.get('sample', ''))
        state = row.get('state', None)
        
        key = (cell_id, sample)
        
        if pd.isna(state):
            # 没有匹配到state
            if key in state_mapping:
                # 应该匹配但没有匹配
                incorrect += 1
            else:
                # 正确未匹配
                correct += 1
        else:
            # 有state值
            if key in state_mapping:
                if state_mapping[key] == state:
                    correct += 1
                else:
                    incorrect += 1
            else:
                incorrect += 1
    
    print(f"总细胞数: {adata.n_obs}")
    print(f"正确匹配: {correct} ({correct/adata.n_obs*100:.2f}%)")
    print(f"错误匹配: {incorrect} ({incorrect/adata.n_obs*100:.2f}%)")
    
    if incorrect > 0:
        print("警告: 存在匹配错误，请检查数据！")

if __name__ == "__main__":
    # 设置文件路径
    h5ad_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/all_filter_modified.h5ad"
    state_csv_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/state.csv"
    
    # 输出路径（如果不指定，会覆盖原文件）
    output_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_with_state.h5ad"
    
    print("选择匹配方法:")
    print("1. 基本方法（推荐）")
    print("2. pandas优化方法（适合大数据集）")
    print("3. 批量处理方法（多个h5ad文件）")
    
    choice = input("\n请选择 (1/2/3, 默认1): ").strip() or "1"
    
    if choice == "1":
        # 使用方法1：基本方法
        adata = add_state_to_h5ad(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path,
            output_path=output_path
        )
        
        # 验证匹配结果
        verify_input = input("\n是否验证匹配结果? (y/n, 默认n): ").strip().lower() or "n"
        if verify_input == "y":
            verify_matching(output_path, state_csv_path)
    
    elif choice == "2":
        # 使用方法2：pandas优化方法
        adata = add_state_to_h5ad_pandas(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path,
            output_path=output_path
        )
    
    elif choice == "3":
        # 使用方法3：批量处理
        h5ad_files = [
            "/path/to/file1.h5ad",
            "/path/to/file2.h5ad",
            # 添加更多文件路径
        ]
        output_dir = "./output_with_state"
        
        add_state_to_multiple_h5ad(
            state_csv_path=state_csv_path,
            h5ad_files=h5ad_files,
            output_dir=output_dir
        )
    
    print("\n" + "=" * 60)
    print("处理完成!")
    print("=" * 60)
