import scanpy as sc
import pandas as pd
import numpy as np
import os

def add_state_new_column(h5ad_path, output_path):
    """
    根据sample_new列映射关系添加state_new列
    """
    
    print("=" * 60)
    print("为h5ad数据添加state_new列")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    print(f"   观测变量(obs)列: {list(adata.obs.columns)}")
    
    # 检查sample_new列是否存在
    if 'sample_new' not in adata.obs.columns:
        raise ValueError("adata.obs中缺少'sample_new'列")
    
    # 2. 定义映射关系
    print(f"\n2. 定义映射关系...")
    
    mapping_dict = {
        "N-P1-1": "N",
        "N-P2-1": "N",
        "N-P2-3": "N",
        "N-P3-3": "N",
        "N-P4-2": "N",
        "N-P4-3": "N",
        "N-P5": "N",
        "N-P6-1": "N",
        "N-P6-2": "N",
        "A-Ind-P1-1": "AI",
        "A-Ind-P1-2": "AI",
        "A-Ind-P2": "AI",
        "A-Ind-P4": "AI",
        "A-Ind-P5": "AI",
        "A-P1-1": "AA",
        "A-P2-1": "AA",
        "A-P2-2": "AA",
        "A-P3-1": "AA",
        "A-P3-2": "AA",
        "A-P3-3": "AA",
        "A-P6-1": "AA",
        "A-P6-2": "AA",
        "C-Ind-P5": "C",
        "C-P1-1": "C",
        "C-P2-1": "C",
        "C-P2-2": "C",
        "C-P3-1": "C",
        "C-P3-2": "C",
        "C-P3-3": "C",
        "C-P6-1": "C",
        "C-P6-2": "C"
    }
    
    print(f"   映射关系数量: {len(mapping_dict)}")
    print(f"   映射值类别: {set(mapping_dict.values())}")
    
    # 显示映射关系
    print(f"\n   映射关系详情:")
    for key, value in mapping_dict.items():
        print(f"     {key} -> {value}")
    
    # 3. 应用映射添加state_new列
    print(f"\n3. 应用映射添加state_new列...")
    
    # 确保sample_new列为字符串类型
    adata.obs['sample_new'] = adata.obs['sample_new'].astype(str)
    
    # 使用map函数进行向量化映射
    adata.obs['state_new'] = adata.obs['sample_new'].map(mapping_dict)
    
    # 统计映射结果
    total_cells = adata.n_obs
    mapped_cells = adata.obs['state_new'].notna().sum()
    unmapped_cells = adata.obs['state_new'].isna().sum()
    unmapped_percentage = unmapped_cells / total_cells * 100 if total_cells > 0 else 0
    
    print(f"\n   映射结果统计:")
    print(f"     总细胞数: {total_cells}")
    print(f"     成功映射: {mapped_cells} ({mapped_cells/total_cells*100:.2f}%)")
    print(f"     未映射: {unmapped_cells} ({unmapped_percentage:.2f}%)")
    
    # 显示未映射的sample_new值
    if unmapped_cells > 0:
        unmapped_samples = adata.obs[adata.obs['state_new'].isna()]['sample_new'].unique()
        print(f"\n   未映射的sample_new值 (前20个):")
        for i, sample in enumerate(unmapped_samples[:20]):
            count = (adata.obs['sample_new'] == sample).sum()
            print(f"     {sample}: {count} 个细胞")
        
        if len(unmapped_samples) > 20:
            print(f"     ... 还有 {len(unmapped_samples) - 20} 个值未显示")
    
    # 4. 查看state_new的分布
    print(f"\n4. state_new列分布:")
    state_new_counts = adata.obs['state_new'].value_counts(dropna=False)
    
    for state, count in state_new_counts.items():
        percentage = count / total_cells * 100
        state_display = '未映射' if pd.isna(state) else state
        print(f"     {state_display}: {count} ({percentage:.2f}%)")
    
    # 5. 保存结果
    print(f"\n5. 保存到文件: {output_path}")
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    # 6. 创建映射报告
    report_path = os.path.splitext(output_path)[0] + "_mapping_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("state_new列映射报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("文件信息:\n")
        f.write(f"  输入文件: {h5ad_path}\n")
        f.write(f"  输出文件: {output_path}\n")
        f.write(f"  报告文件: {report_path}\n\n")
        
        f.write("数据统计:\n")
        f.write(f"  总细胞数: {total_cells}\n")
        f.write(f"  总基因数: {adata.n_vars}\n")
        f.write(f"  sample_new唯一值数: {adata.obs['sample_new'].nunique()}\n\n")
        
        f.write("映射结果:\n")
        f.write(f"  成功映射: {mapped_cells} ({mapped_cells/total_cells*100:.2f}%)\n")
        f.write(f"  未映射: {unmapped_cells} ({unmapped_percentage:.2f}%)\n\n")
        
        f.write("映射关系:\n")
        for key, value in mapping_dict.items():
            f.write(f"  {key} -> {value}\n")
        f.write("\n")
        
        f.write("state_new列分布:\n")
        for state, count in state_new_counts.items():
            percentage = count / total_cells * 100
            state_display = '未映射' if pd.isna(state) else state
            f.write(f"  {state_display}: {count} ({percentage:.2f}%)\n")
        
        if unmapped_cells > 0:
            f.write(f"\n未映射的sample_new值:\n")
            for sample in sorted(unmapped_samples):
                count = (adata.obs['sample_new'] == sample).sum()
                f.write(f"  {sample}: {count} 个细胞\n")
    
    print(f"\n6. 映射报告已保存到: {report_path}")
    
    return adata

def add_state_new_column_advanced(h5ad_path, output_path, custom_mapping=None):
    """
    为h5ad数据添加state_new列 - 增强版
    支持自定义映射关系
    """
    
    print("=" * 60)
    print("为h5ad数据添加state_new列 - 增强版")
    print("=" * 60)
    
    # 1. 读取h5ad数据
    print(f"\n1. 读取h5ad数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"   数据形状: {adata.shape}")
    
    # 检查sample_new列是否存在
    if 'sample_new' not in adata.obs.columns:
        raise ValueError("adata.obs中缺少'sample_new'列")
    
    # 2. 使用提供的映射或默认映射
    if custom_mapping is not None:
        mapping_dict = custom_mapping
        print(f"   使用自定义映射关系")
    else:
        # 默认映射关系
        mapping_dict = {
            "N-P1-1": "N",
            "N-P2-1": "N",
            "N-P2-3": "N",
            "N-P3-3": "N",
            "N-P4-2": "N",
            "N-P4-3": "N",
            "N-P5": "N",
            "N-P6-1": "N",
            "N-P6-2": "N",
            "A-Ind-P1-1": "AI",
            "A-Ind-P1-2": "AI",
            "A-Ind-P2": "AI",
            "A-Ind-P4": "AI",
            "A-Ind-P5": "AI",
            "A-P1-1": "AA",
            "A-P2-1": "AA",
            "A-P2-2": "AA",
            "A-P3-1": "AA",
            "A-P3-2": "AA",
            "A-P3-3": "AA",
            "A-P6-1": "AA",
            "A-P6-2": "AA",
            "C-Ind-P5": "C",
            "C-P1-1": "C",
            "C-P2-1": "C",
            "C-P2-2": "C",
            "C-P3-1": "C",
            "C-P3-2": "C",
            "C-P3-3": "C",
            "C-P6-1": "C",
            "C-P6-2": "C"
        }
        print(f"   使用默认映射关系")
    
    print(f"   映射关系数量: {len(mapping_dict)}")
    
    # 3. 应用映射
    print(f"\n2. 应用映射...")
    
    # 确保sample_new列为字符串类型
    adata.obs['sample_new'] = adata.obs['sample_new'].astype(str)
    
    # 使用向量化映射
    adata.obs['state_new'] = adata.obs['sample_new'].map(mapping_dict)
    
    # 统计
    total_cells = adata.n_obs
    mapped_cells = adata.obs['state_new'].notna().sum()
    
    print(f"   成功映射: {mapped_cells}/{total_cells} ({mapped_cells/total_cells*100:.2f}%)")
    
    # 4. 显示state_new分布
    print(f"\n3. state_new分布:")
    for state, count in adata.obs['state_new'].value_counts(dropna=True).items():
        percentage = count / total_cells * 100
        print(f"   {state}: {count} ({percentage:.2f}%)")
    
    # 5. 保存结果
    print(f"\n4. 保存到文件: {output_path}")
    adata.write_h5ad(output_path)
    print(f"   数据已保存!")
    
    return adata

def main():
    """
    主函数
    """
    
    # 设置文件路径
    input_file = "All_filter_state_samp_new.h5ad"
    output_file = "All_filter_state_new_samp_new.h5ad"
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在")
        print(f"当前目录: {os.getcwd()}")
        print(f"目录内容: {os.listdir('.')}")
        return
    
    print(f"当前工作目录: {os.getcwd()}")
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    
    # 执行映射
    try:
        adata_processed = add_state_new_column(input_file, output_file)
        
        # 额外检查：验证结果
        print(f"\n" + "="*60)
        print("验证结果")
        print("="*60)
        
        print(f"\n1. 检查新添加的列:")
        print(f"   obs列: {list(adata_processed.obs.columns)}")
        
        print(f"\n2. 检查state_new列的值:")
        unique_states = adata_processed.obs['state_new'].dropna().unique()
        print(f"   state_new唯一值: {unique_states}")
        
        print(f"\n3. 样本到state_new的映射示例:")
        sample_state_mapping = adata_processed.obs[['sample_new', 'state_new']].drop_duplicates()
        for i, (sample, state) in enumerate(sample_state_mapping.head(10).values):
            print(f"   {sample} -> {state}")
        
        if len(sample_state_mapping) > 10:
            print(f"   ... 还有 {len(sample_state_mapping) - 10} 个映射未显示")
        
        print(f"\n" + "="*60)
        print("处理完成!")
        print("="*60)
        
    except Exception as e:
        print(f"处理过程中发生错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
