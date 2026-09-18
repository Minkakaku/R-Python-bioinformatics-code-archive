
import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

def add_state_to_h5ad_with_detailed_check(h5ad_path, state_csv_path, output_path=None):
    """
    详细检查并匹配state.csv到h5ad数据，找出不匹配的样本
    
    参数:
    h5ad_path: h5ad文件路径
    state_csv_path: state.csv文件路径
    output_path: 输出文件路径
    """
    
    print("=" * 60)
    print("详细检查并匹配state.csv到h5ad数据")
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
    
    # 3. 详细检查数据
    print(f"\n3. 详细检查数据...")
    
    # 检查h5ad中的样本
    h5ad_samples = set(adata.obs['sample'].astype(str).unique())
    print(f"   h5ad中共有 {len(h5ad_samples)} 个样本")
    print(f"   h5ad样本列表: {sorted(list(h5ad_samples))[:10]}")  # 只显示前10个
    
    # 检查state.csv中的样本
    state_samples = set(state_df['sample'].astype(str).unique())
    print(f"   state.csv中共有 {len(state_samples)} 个样本")
    print(f"   state.csv样本列表: {sorted(list(state_samples))[:10]}")
    
    # 找出不匹配的样本
    h5ad_only_samples = h5ad_samples - state_samples
    state_only_samples = state_samples - h5ad_samples
    
    print(f"\n   在h5ad中但不在state.csv中的样本 ({len(h5ad_only_samples)}个):")
    if h5ad_only_samples:
        for sample in sorted(list(h5ad_only_samples)):
            # 统计这个样本在h5ad中的细胞数
            count = (adata.obs['sample'].astype(str) == sample).sum()
            print(f"     - {sample}: {count} 个细胞")
    else:
        print("     无")
    
    print(f"\n   在state.csv中但不在h5ad中的样本 ({len(state_only_samples)}个):")
    if state_only_samples:
        for sample in sorted(list(state_only_samples)):
            # 统计这个样本在state.csv中的记录数
            count = (state_df['sample'].astype(str) == sample).sum()
            print(f"     - {sample}: {count} 条记录")
    else:
        print("     无")
    
    # 4. 检查重复
    print(f"\n4. 检查重复记录...")
    
    # 检查h5ad中的重复(cell_id, sample)
    h5ad_duplicates = adata.obs.duplicated(subset=['cell_id', 'sample'], keep=False)
    if h5ad_duplicates.any():
        print(f"   警告: h5ad中有 {h5ad_duplicates.sum()} 个重复的(cell_id, sample)组合")
    
    # 检查state.csv中的重复
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
    
    # 5. 使用字典映射进行匹配
    print(f"\n5. 使用字典映射进行匹配...")
    
    # 创建映射字典
    state_mapping = {}
    for _, row in state_df.iterrows():
        key = (str(row['cell_id']), str(row['sample']))
        state_mapping[key] = row['group']
    
    print(f"   创建了 {len(state_mapping)} 个映射条目")
    
    # 6. 批量匹配
    print(f"\n6. 批量匹配group信息到h5ad...")
    
    # 创建state列，初始值为NaN
    adata.obs['state'] = np.nan
    
    # 准备数据
    h5ad_cell_ids = adata.obs['cell_id'].astype(str).tolist()
    h5ad_samples_list = adata.obs['sample'].astype(str).tolist()
    h5ad_indices = adata.obs.index.tolist()
    
    matched_count = 0
    unmatched_count = 0
    
    # 按样本统计匹配情况
    sample_match_stats = {}
    
    # 批量查找
    for idx, cell_id, sample in zip(h5ad_indices, h5ad_cell_ids, h5ad_samples_list):
        key = (cell_id, sample)
        if key in state_mapping:
            adata.obs.at[idx, 'state'] = state_mapping[key]
            matched_count += 1
            
            # 统计样本匹配
            if sample not in sample_match_stats:
                sample_match_stats[sample] = {'matched': 0, 'unmatched': 0}
            sample_match_stats[sample]['matched'] += 1
        else:
            unmatched_count += 1
            
            # 统计样本不匹配
            if sample not in sample_match_stats:
                sample_match_stats[sample] = {'matched': 0, 'unmatched': 0}
            sample_match_stats[sample]['unmatched'] += 1
    
    # 7. 统计匹配结果
    print(f"\n7. 匹配结果统计:")
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   匹配成功: {matched_count} ({matched_count/adata.n_obs*100:.2f}%)")
    print(f"   未匹配: {unmatched_count} ({unmatched_count/adata.n_obs*100:.2f}%)")
    
    # 按样本显示匹配率
    print(f"\n   各样本匹配率:")
    for sample in sorted(sample_match_stats.keys()):
        stats = sample_match_stats[sample]
        total = stats['matched'] + stats['unmatched']
        if total > 0:
            match_rate = stats['matched'] / total * 100
            if match_rate < 100:
                print(f"     {sample}: {match_rate:.2f}% ({stats['matched']}/{total} 个细胞) - 不完整匹配!")
            else:
                print(f"     {sample}: {match_rate:.2f}% ({stats['matched']}/{total} 个细胞)")
    
    # 找出完全未匹配的样本
    completely_unmatched_samples = []
    for sample, stats in sample_match_stats.items():
        if stats['matched'] == 0 and stats['unmatched'] > 0:
            completely_unmatched_samples.append(sample)
    
    if completely_unmatched_samples:
        print(f"\n   完全未匹配的样本 ({len(completely_unmatched_samples)}个):")
        for sample in completely_unmatched_samples:
            count = sample_match_stats[sample]['unmatched']
            print(f"     - {sample}: {count} 个细胞")
            
            # 检查这个样本是否在state.csv中
            if sample in h5ad_only_samples:
                print(f"       原因: 样本不在state.csv中")
            else:
                # 检查是否有cell_id匹配问题
                sample_cells = adata[adata.obs['sample'] == sample, :]
                sample_cell_ids = set(sample_cells.obs['cell_id'].astype(str))
                
                # 获取state.csv中这个样本的cell_id
                state_sample_df = state_df[state_df['sample'].astype(str) == sample]
                state_cell_ids = set(state_sample_df['cell_id'].astype(str))
                
                # 计算交集
                common_cell_ids = sample_cell_ids.intersection(state_cell_ids)
                print(f"       样本在state.csv中有 {len(state_cell_ids)} 个cell_id记录")
                print(f"       h5ad中有 {len(sample_cell_ids)} 个cell_id")
                print(f"       共同cell_id: {len(common_cell_ids)} 个")
                
                if len(common_cell_ids) == 0:
                    print(f"       原因: cell_id完全不匹配")
                else:
                    print(f"       原因: 部分cell_id不匹配")
    
    # 8. 保存结果
    if output_path is None:
        output_path = h5ad_path
        print(f"\n8. 正在覆盖原始文件: {output_path}")
    else:
        print(f"\n8. 正在保存到新文件: {output_path}")
    
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    # 9. 创建详细的匹配结果报告
    report_path = os.path.splitext(output_path)[0] + "_detailed_matching_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("State信息匹配详细报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("文件信息:\n")
        f.write(f"  h5ad文件: {h5ad_path}\n")
        f.write(f"  state.csv文件: {state_csv_path}\n")
        f.write(f"  输出文件: {output_path}\n\n")
        
        f.write("数据统计:\n")
        f.write(f"  h5ad总细胞数: {adata.n_obs}\n")
        f.write(f"  h5ad总基因数: {adata.n_vars}\n")
        f.write(f"  h5ad样本数: {len(h5ad_samples)}\n")
        f.write(f"  state.csv记录数: {len(state_df)}\n")
        f.write(f"  state.csv样本数: {len(state_samples)}\n")
        f.write(f"  state.csv唯一cell_id数: {state_df['cell_id'].nunique()}\n")
        f.write(f"  state.csv唯一sample数: {state_df['sample'].nunique()}\n\n")
        
        f.write("样本匹配情况:\n")
        f.write(f"  h5ad中但不在state.csv中的样本 ({len(h5ad_only_samples)}个):\n")
        for sample in sorted(list(h5ad_only_samples)):
            count = (adata.obs['sample'].astype(str) == sample).sum()
            f.write(f"    {sample}: {count} 个细胞\n")
        f.write(f"\n  state.csv中但不在h5ad中的样本 ({len(state_only_samples)}个):\n")
        for sample in sorted(list(state_only_samples)):
            count = (state_df['sample'].astype(str) == sample).sum()
            f.write(f"    {sample}: {count} 条记录\n")
        f.write("\n")
        
        f.write("匹配结果:\n")
        f.write(f"  匹配成功: {matched_count} ({matched_count/adata.n_obs*100:.2f}%)\n")
        f.write(f"  未匹配: {unmatched_count} ({unmatched_count/adata.n_obs*100:.2f}%)\n\n")
        
        f.write("各样本匹配详情:\n")
        for sample in sorted(sample_match_stats.keys()):
            stats = sample_match_stats[sample]
            total = stats['matched'] + stats['unmatched']
            match_rate = stats['matched'] / total * 100 if total > 0 else 0
            f.write(f"  {sample}: {match_rate:.2f}% ({stats['matched']}/{total} 个细胞)\n")
        f.write("\n")
        
        f.write("完全未匹配的样本:\n")
        if completely_unmatched_samples:
            for sample in completely_unmatched_samples:
                count = sample_match_stats[sample]['unmatched']
                f.write(f"  {sample}: {count} 个细胞未匹配\n")
        else:
            f.write("  无\n")
    
    print(f"\n9. 详细匹配报告已保存到: {report_path}")
    
    return adata

def quick_check_unmatched_samples(h5ad_path, state_csv_path):
    """
    快速检查不匹配的样本
    """
    
    print("=" * 60)
    print("快速检查不匹配的样本")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    
    # 2. 读取state.csv数据
    print(f"\n2. 读取state.csv数据: {state_csv_path}")
    state_df = pd.read_csv(state_csv_path)
    
    # 3. 检查样本
    h5ad_samples = set(adata.obs['sample'].astype(str).unique())
    state_samples = set(state_df['sample'].astype(str).unique())
    
    print(f"\n3. 样本检查:")
    print(f"   h5ad样本数: {len(h5ad_samples)}")
    print(f"   state.csv样本数: {len(state_samples)}")
    
    # 找出不匹配的样本
    h5ad_only_samples = h5ad_samples - state_samples
    state_only_samples = state_samples - h5ad_samples
    
    print(f"\n4. 不匹配的样本:")
    print(f"   在h5ad中但不在state.csv中的样本 ({len(h5ad_only_samples)}个):")
    for sample in sorted(list(h5ad_only_samples)):
        count = (adata.obs['sample'].astype(str) == sample).sum()
        print(f"     - {sample}: {count} 个细胞")
    
    print(f"\n   在state.csv中但不在h5ad中的样本 ({len(state_only_samples)}个):")
    for sample in sorted(list(state_only_samples)):
        count = (state_df['sample'].astype(str) == sample).sum()
        print(f"     - {sample}: {count} 条记录")
    
    # 4. 如果h5ad已经有state列，检查不匹配的情况
    if 'state' in adata.obs.columns:
        print(f"\n5. 现有state列分析:")
        
        # 检查state列中NaN的比例
        nan_count = adata.obs['state'].isna().sum()
        nan_percentage = nan_count / adata.n_obs * 100
        print(f"   state列中NaN值: {nan_count} ({nan_percentage:.2f}%)")
        
        # 按样本检查
        print(f"\n   各样本未匹配情况:")
        for sample in sorted(list(h5ad_samples)):
            sample_data = adata[adata.obs['sample'].astype(str) == sample, :]
            nan_in_sample = sample_data.obs['state'].isna().sum()
            if nan_in_sample > 0:
                percentage = nan_in_sample / len(sample_data.obs) * 100
                print(f"     {sample}: {nan_in_sample}/{len(sample_data.obs)} ({percentage:.2f}%) 未匹配")

if __name__ == "__main__":
    # 设置文件路径
    h5ad_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/all_filter_modified.h5ad"
    state_csv_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/state.csv"
    
    # 输出路径
    output_path = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_with_state_checked.h5ad"
    
    print("选择检查方式:")
    print("1. 快速检查不匹配的样本")
    print("2. 详细匹配并分析原因")
    
    choice = input("\n请选择 (1/2, 默认1): ").strip() or "1"
    
    if choice == "1":
        # 快速检查
        quick_check_unmatched_samples(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path
        )
        
    elif choice == "2":
        # 详细匹配
        adata = add_state_to_h5ad_with_detailed_check(
            h5ad_path=h5ad_path,
            state_csv_path=state_csv_path,
            output_path=output_path
        )
    
    print("\n" + "=" * 60)
    print("检查完成!")
    print("=" * 60)
