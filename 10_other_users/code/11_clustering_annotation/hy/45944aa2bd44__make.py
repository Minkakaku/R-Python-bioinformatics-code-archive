import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

def filter_state_nan_and_save(h5ad_path, output_path):
    """
    过滤掉state列为NaN的细胞并保存
    
    参数:
    h5ad_path: 输入h5ad文件路径
    output_path: 输出h5ad文件路径
    """
    
    print("=" * 60)
    print("过滤state列为NaN的细胞")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   原始数据形状: {adata.shape}")
    
    # 检查是否有state列
    if 'state' not in adata.obs.columns:
        raise ValueError("adata.obs中缺少'state'列")
    
    # 2. 统计过滤前的state分布
    print(f"\n2. 过滤前state分布:")
    state_counts_before = adata.obs['state'].value_counts(dropna=False)
    for state, count in state_counts_before.items():
        percentage = count / adata.n_obs * 100
        state_display = 'NaN' if pd.isna(state) else state
        print(f"   {state_display}: {count} ({percentage:.2f}%)")
    
    # 3. 过滤掉state为NaN的细胞
    print(f"\n3. 过滤state为NaN的细胞...")
    
    # 记录原始细胞数
    original_cells = adata.n_obs
    
    # 过滤掉state为NaN的细胞
    adata_filtered = adata[~adata.obs['state'].isna()].copy()
    
    # 统计过滤结果
    filtered_cells = adata_filtered.n_obs
    removed_cells = original_cells - filtered_cells
    removed_percentage = removed_cells / original_cells * 100
    
    print(f"   原始细胞数: {original_cells}")
    print(f"   过滤后细胞数: {filtered_cells}")
    print(f"   移除细胞数: {removed_cells} ({removed_percentage:.2f}%)")
    
    # 4. 保存过滤后的数据
    print(f"\n4. 保存过滤后的数据: {output_path}")
    adata_filtered.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    # 5. 创建过滤报告
    report_path = os.path.splitext(output_path)[0] + "_filter_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("State过滤报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("文件信息:\n")
        f.write(f"  输入文件: {h5ad_path}\n")
        f.write(f"  输出文件: {output_path}\n\n")
        
        f.write("过滤统计:\n")
        f.write(f"  原始细胞数: {original_cells}\n")
        f.write(f"  过滤后细胞数: {filtered_cells}\n")
        f.write(f"  移除细胞数: {removed_cells} ({removed_percentage:.2f}%)\n\n")
        
        f.write("过滤前state分布:\n")
        for state, count in state_counts_before.items():
            percentage = count / original_cells * 100
            state_display = 'NaN' if pd.isna(state) else state
            f.write(f"  {state_display}: {count} ({percentage:.2f}%)\n")
        
        f.write("\n过滤后state分布:\n")
        state_counts_after = adata_filtered.obs['state'].value_counts()
        for state, count in state_counts_after.items():
            percentage = count / filtered_cells * 100
            f.write(f"  {state}: {count} ({percentage:.2f}%)\n")
    
    print(f"\n5. 过滤报告已保存到: {report_path}")
    
    return adata_filtered

def add_sample_new_from_releate(h5ad_path, releate_path, output_path):
    """
    根据releate文件为h5ad数据添加sample_new列 - 优化版
    
    参数:
    h5ad_path: 输入h5ad文件路径
    releate_path: releate文件路径
    output_path: 输出h5ad文件路径
    """
    
    print("=" * 60)
    print("根据releate文件添加sample_new列 - 优化版")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    print(f"   观测变量(obs)列: {list(adata.obs.columns)}")
    
    # 检查必要的列
    required_cols = ['sample', 'state']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"adata.obs中缺少必需的列: {missing_cols}")
    
    # 2. 读取releate文件
    print(f"\n2. 读取releate文件: {releate_path}")
    
    # 尝试不同的分隔符
    try:
        releate_df = pd.read_csv(releate_path, sep='\t')
        print(f"   使用制表符分隔符读取成功")
    except:
        try:
            releate_df = pd.read_csv(releate_path, sep=',')
            print(f"   使用逗号分隔符读取成功")
        except:
            try:
                releate_df = pd.read_csv(releate_path, sep=' ')
                print(f"   使用空格分隔符读取成功")
            except:
                raise ValueError(f"无法读取releate文件，请检查文件格式")
    
    print(f"   releate文件形状: {releate_df.shape}")
    print(f"   releate文件列: {list(releate_df.columns)}")
    
    # 检查releate文件的列
    required_releate_cols = ['sample', 'sample_new', 'state']
    missing_releate_cols = [col for col in required_releate_cols if col not in releate_df.columns]
    if missing_releate_cols:
        raise ValueError(f"releate文件中缺少必需的列: {missing_releate_cols}")
    
    # 显示releate文件的前几行
    print(f"\n   releate文件前10行:")
    print(releate_df.head(10))
    
    # 3. 检查数据
    print(f"\n3. 检查数据...")
    
    # 确保列是字符串类型
    for col in ['sample', 'state']:
        releate_df[col] = releate_df[col].astype(str)
        adata.obs[col] = adata.obs[col].astype(str)
    
    # 获取h5ad中的唯一(sample, state)组合
    h5ad_combinations = set(zip(adata.obs['sample'], adata.obs['state']))
    print(f"   h5ad中有 {len(h5ad_combinations)} 个唯一的(sample, state)组合")
    
    # 获取releate中的唯一(sample, state)组合
    releate_combinations = set(zip(releate_df['sample'], releate_df['state']))
    print(f"   releate中有 {len(releate_combinations)} 个唯一的(sample, state)组合")
    
    # 找出不匹配的组合
    matched_combinations = h5ad_combinations.intersection(releate_combinations)
    unmatched_h5ad = h5ad_combinations - releate_combinations
    unmatched_releate = releate_combinations - h5ad_combinations
    
    print(f"\n   匹配情况:")
    print(f"     匹配的组合: {len(matched_combinations)} 个")
    print(f"     h5ad中未匹配的组合: {len(unmatched_h5ad)} 个")
    print(f"     releate中未匹配的组合: {len(unmatched_releate)} 个")
    
    # 显示部分未匹配的组合
    if unmatched_h5ad:
        print(f"\n   h5ad中未匹配的组合示例（前10个）:")
        for i, (sample, state) in enumerate(sorted(list(unmatched_h5ad))[:10]):
            count = ((adata.obs['sample'] == sample) & (adata.obs['state'] == state)).sum()
            print(f"     ({sample}, {state}): {count} 个细胞")
    
    if unmatched_releate:
        print(f"\n   releate中未匹配的组合示例（前10个）:")
        for i, (sample, state) in enumerate(sorted(list(unmatched_releate))[:10]):
            sample_new = releate_df[
                (releate_df['sample'] == sample) & (releate_df['state'] == state)
            ]['sample_new'].iloc[0]
            print(f"     ({sample}, {state}): 应映射到 '{sample_new}'")
    
    # 4. 创建映射字典 - 优化为pandas Series用于向量化操作
    print(f"\n4. 创建映射...")
    
    # 创建从(sample, state)到sample_new的映射Series
    # 先确保releate_df中的(sample, state)组合是唯一的
    duplicate_mask = releate_df.duplicated(subset=['sample', 'state'], keep=False)
    if duplicate_mask.any():
        print(f"   警告: releate文件中有重复的(sample, state)组合")
        print(f"   重复组合数量: {duplicate_mask.sum()}")
        # 保留第一个出现的
        releate_df = releate_df.drop_duplicates(subset=['sample', 'state'], keep='first')
    
    # 创建映射Series
    mapping_series = releate_df.set_index(['sample', 'state'])['sample_new'].astype(str)
    print(f"   创建了 {len(mapping_series)} 个映射条目")
    
    # 5. 应用映射 - 使用向量化操作替代循环
    print(f"\n5. 应用映射添加sample_new列...")
    
    # 创建临时MultiIndex用于映射
    temp_index = pd.MultiIndex.from_arrays([adata.obs['sample'], adata.obs['state']], names=['sample', 'state'])
    
    # 使用map进行向量化映射
    adata.obs['sample_new'] = pd.Series(temp_index).map(mapping_series)
    
    # 统计匹配情况
    matched_cells = adata.obs['sample_new'].notna().sum()
    unmatched_cells = adata.obs['sample_new'].isna().sum()
    
    # 对于没有映射的，使用原sample值
    adata.obs['sample_new'] = adata.obs['sample_new'].fillna(adata.obs['sample'])
    
    # 6. 统计结果
    print(f"\n6. 匹配结果统计:")
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   成功匹配: {matched_cells} ({matched_cells/adata.n_obs*100:.2f}%)")
    print(f"   未匹配(保留原sample): {unmatched_cells} ({unmatched_cells/adata.n_obs*100:.2f}%)")
    
    # 查看sample_new的分布
    sample_new_counts = adata.obs['sample_new'].value_counts(dropna=False)
    print(f"\n   sample_new列分布（前20个）:")
    for sample_new, count in sample_new_counts.head(20).items():
        percentage = count / adata.n_obs * 100
        print(f"     {sample_new}: {count} ({percentage:.2f}%)")
    
    if len(sample_new_counts) > 20:
        print(f"     ... 还有 {len(sample_new_counts) - 20} 个值未显示")
    
    # 7. 保存结果
    print(f"\n7. 正在保存到文件: {output_path}")
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    # 8. 创建匹配报告
    report_path = os.path.splitext(output_path)[0] + "_matching_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("sample_new列匹配报告 - 优化版\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("文件信息:\n")
        f.write(f"  h5ad文件: {h5ad_path}\n")
        f.write(f"  releate文件: {releate_path}\n")
        f.write(f"  输出文件: {output_path}\n\n")
        
        f.write("数据统计:\n")
        f.write(f"  总细胞数: {adata.n_obs}\n")
        f.write(f"  总基因数: {adata.n_vars}\n")
        f.write(f"  h5ad中唯一的(sample, state)组合数: {len(h5ad_combinations)}\n")
        f.write(f"  releate中唯一的(sample, state)组合数: {len(releate_combinations)}\n\n")
        
        f.write("匹配统计:\n")
        f.write(f"  匹配的组合数: {len(matched_combinations)}\n")
        f.write(f"  h5ad中未匹配的组合数: {len(unmatched_h5ad)}\n")
        f.write(f"  releate中未匹配的组合数: {len(unmatched_releate)}\n\n")
        
        f.write("细胞匹配结果:\n")
        f.write(f"  成功匹配: {matched_cells} ({matched_cells/adata.n_obs*100:.2f}%)\n")
        f.write(f"  未匹配(保留原sample): {unmatched_cells} ({unmatched_cells/adata.n_obs*100:.2f}%)\n\n")
        
        f.write("映射关系:\n")
        for (sample, state), sample_new in mapping_series.items():
            f.write(f"  ({sample}, {state}) -> {sample_new}\n")
        f.write("\n")
        
        f.write("sample_new列分布:\n")
        for sample_new, count in sample_new_counts.items():
            percentage = count / adata.n_obs * 100
            f.write(f"  {sample_new}: {count} ({percentage:.2f}%)\n")
        
        if unmatched_h5ad:
            f.write(f"\nh5ad中未匹配的组合（前50个）:\n")
            for i, (sample, state) in enumerate(sorted(list(unmatched_h5ad))[:50]):
                count = ((adata.obs['sample'] == sample) & (adata.obs['state'] == state)).sum()
                f.write(f"  ({sample}, {state}): {count} 个细胞\n")
        
        if unmatched_releate:
            f.write(f"\nreleate中未匹配的组合:\n")
            for i, (sample, state) in enumerate(sorted(list(unmatched_releate))):
                sample_new = releate_df[
                    (releate_df['sample'] == sample) & (releate_df['state'] == state)
                ]['sample_new'].iloc[0]
                f.write(f"  ({sample}, {state}): 应映射到 '{sample_new}'\n")
    
    print(f"\n8. 匹配报告已保存到: {report_path}")
    
    return adata

def add_sample_new_from_releate_ultrafast(h5ad_path, releate_path, output_path):
    """
    根据releate文件为h5ad数据添加sample_new列 - 极速版
    
    这个版本使用了pandas的merge操作，速度最快
    """
    
    print("=" * 60)
    print("根据releate文件添加sample_new列 - 极速版")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    print(f"   观测变量(obs)列: {list(adata.obs.columns)}")
    
    # 检查必要的列
    required_cols = ['sample', 'state']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"adata.obs中缺少必需的列: {missing_cols}")
    
    # 2. 读取releate文件
    print(f"\n2. 读取releate文件: {releate_path}")
    
    # 尝试不同的分隔符
    try:
        releate_df = pd.read_csv(releate_path, sep='\t')
        print(f"   使用制表符分隔符读取成功")
    except:
        try:
            releate_df = pd.read_csv(releate_path, sep=',')
            print(f"   使用逗号分隔符读取成功")
        except:
            try:
                releate_df = pd.read_csv(releate_path, sep=' ')
                print(f"   使用空格分隔符读取成功")
            except:
                raise ValueError(f"无法读取releate文件，请检查文件格式")
    
    print(f"   releate文件形状: {releate_df.shape}")
    print(f"   releate文件列: {list(releate_df.columns)}")
    
    # 检查releate文件的列
    required_releate_cols = ['sample', 'sample_new', 'state']
    missing_releate_cols = [col for col in required_releate_cols if col not in releate_df.columns]
    if missing_releate_cols:
        raise ValueError(f"releate文件中缺少必需的列: {missing_releate_cols}")
    
    # 确保列是字符串类型
    for col in ['sample', 'state', 'sample_new']:
        releate_df[col] = releate_df[col].astype(str)
        if col in ['sample', 'state']:
            adata.obs[col] = adata.obs[col].astype(str)
    
    # 3. 使用merge进行快速匹配
    print(f"\n3. 使用pandas merge进行快速匹配...")
    
    # 保存原始索引
    original_index = adata.obs.index
    
    # 创建临时DataFrame用于合并
    obs_df = adata.obs.reset_index(drop=False)
    
    # 去重releate_df，避免合并时产生重复行
    releate_df_unique = releate_df.drop_duplicates(subset=['sample', 'state'])
    
    # 使用左连接合并
    merged_df = obs_df.merge(releate_df_unique[['sample', 'state', 'sample_new']], 
                            on=['sample', 'state'], 
                            how='left', 
                            suffixes=('', '_mapped'))
    
    # 确保顺序不变
    merged_df = merged_df.set_index(original_index.name or 'index')
    merged_df.index = original_index
    
    # 填充没有匹配到的sample_new
    merged_df['sample_new'] = merged_df['sample_new'].fillna(merged_df['sample'])
    
    # 更新adata.obs
    adata.obs['sample_new'] = merged_df['sample_new']
    
    # 统计匹配情况
    matched_cells = (adata.obs['sample_new'] != adata.obs['sample']).sum()
    unmatched_cells = (adata.obs['sample_new'] == adata.obs['sample']).sum()
    
    print(f"\n4. 匹配结果统计:")
    print(f"   总细胞数: {adata.n_obs}")
    print(f"   成功匹配: {matched_cells} ({matched_cells/adata.n_obs*100:.2f}%)")
    print(f"   未匹配(保留原sample): {unmatched_cells} ({unmatched_cells/adata.n_obs*100:.2f}%)")
    
    # 查看sample_new的分布
    sample_new_counts = adata.obs['sample_new'].value_counts(dropna=False)
    print(f"\n   sample_new列分布（前20个）:")
    for sample_new, count in sample_new_counts.head(20).items():
        percentage = count / adata.n_obs * 100
        print(f"     {sample_new}: {count} ({percentage:.2f}%)")
    
    if len(sample_new_counts) > 20:
        print(f"     ... 还有 {len(sample_new_counts) - 20} 个值未显示")
    
    # 5. 保存结果
    print(f"\n5. 正在保存到文件: {output_path}")
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    return adata

def process_all_steps():
    """
    完整处理流程
    """
    
    print("=" * 60)
    print("完整处理流程：过滤state为NaN的细胞并添加sample_new列")
    print("=" * 60)
    
    # 设置文件路径
    input_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_with_state_new.h5ad"
    releate_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/releate"  # 请确认文件扩展名
    intermediate_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state.h5ad"
    final_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state_samp_new.h5ad"
    
    # 步骤1: 过滤state为NaN的细胞
    print(f"\n步骤1: 过滤state为NaN的细胞")
    print("-" * 40)
    
    adata_filtered = filter_state_nan_and_save(
        h5ad_path=input_h5ad,
        output_path=intermediate_h5ad
    )
    
    # 步骤2: 根据releate文件添加sample_new列
    print(f"\n步骤2: 根据releate文件添加sample_new列")
    print("-" * 40)
    
    # 使用极速版
    adata_final = add_sample_new_from_releate_ultrafast(
        h5ad_path=intermediate_h5ad,
        releate_path=releate_file,
        output_path=final_h5ad
    )
    
    # 最终统计
    print(f"\n最终结果统计:")
    print(f"  最终数据形状: {adata_final.shape}")
    print(f"  最终文件: {final_h5ad}")
    
    # 显示最终样本分布
    print(f"\n最终样本分布:")
    sample_counts = adata_final.obs['sample_new'].value_counts()
    for sample, count in sample_counts.head(20).items():
        percentage = count / adata_final.n_obs * 100
        print(f"  {sample}: {count} ({percentage:.2f}%)")
    
    if len(sample_counts) > 20:
        print(f"  ... 还有 {len(sample_counts) - 20} 个样本未显示")
    
    print("\n" + "=" * 60)
    print("所有处理完成!")
    print("=" * 60)
    
    return adata_final

def main():
    """
    主函数，提供交互式选择
    """
    
    print("选择处理方式:")
    print("1. 完整流程：先过滤state为NaN，再添加sample_new列")
    print("2. 只过滤state为NaN的细胞")
    print("3. 只为已有h5ad添加sample_new列（优化版）")
    print("4. 只为已有h5ad添加sample_new列（极速版）")
    
    choice = input("\n请选择 (1/2/3/4, 默认1): ").strip() or "1"
    
    if choice == "1":
        # 完整流程
        adata_final = process_all_steps()
    
    elif choice == "2":
        # 只过滤state为NaN
        input_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_with_state_new.h5ad"
        output_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state.h5ad"
        
        adata_filtered = filter_state_nan_and_save(
            h5ad_path=input_h5ad,
            output_path=output_h5ad
        )
    
    elif choice == "3":
        # 只为已有h5ad添加sample_new列（优化版）
        input_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state.h5ad"
        releate_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/releate"
        output_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state_samp_new.h5ad"
        
        adata_final = add_sample_new_from_releate(
            h5ad_path=input_h5ad,
            releate_path=releate_file,
            output_path=output_h5ad
        )
    
    elif choice == "4":
        # 只为已有h5ad添加sample_new列（极速版）
        input_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state.h5ad"
        releate_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/releate"
        output_h5ad = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_state_samp_new.h5ad"
        
        adata_final = add_sample_new_from_releate_ultrafast(
            h5ad_path=input_h5ad,
            releate_path=releate_file,
            output_path=output_h5ad
        )

if __name__ == "__main__":
    main()
