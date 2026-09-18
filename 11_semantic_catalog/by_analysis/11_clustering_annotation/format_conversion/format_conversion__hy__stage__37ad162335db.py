import pandas as pd
import scanpy as sc

def add_stage_column(adata_path, csv_files_dir, output_path=None):
    """
    为h5ad数据添加stage列，基于obs中的'cell_id'列与CSV文件进行匹配
    """
    # 1. 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(adata_path)
    print(f"数据形状: {adata.shape}")
    
    # 检查是否存在'cell_id'列
    if 'cell_id' not in adata.obs.columns:
        raise KeyError("在adata.obs中未找到'cell_id'列，请检查数据。")
    
    # 2. 初始化stage列为"others"
    adata.obs['stage'] = 'others'
    print("已初始化stage列为'others'")
    
    # 3. 直接指定的样本stage赋值
    # 009D-2为Normal
    adata.obs.loc[adata.obs['sample'] == '009D-2', 'stage'] = 'Normal'
    
    # SR202500002C-1和SR202500002C-2为Adenoma
    adata.obs.loc[adata.obs['sample'].isin(['SR202500002C-1', 'SR202500002C-2']), 'stage'] = 'Adenoma'
    
    # SR202500005F-2为Normal
    adata.obs.loc[adata.obs['sample'] == 'SR202500005F-2', 'stage'] = 'Normal'
    
    print("已完成直接指定的样本stage赋值")
    
    # 4. 从CSV文件读取并赋值的样本列表
    csv_samples = ['009D-1', 'SR202500003B', 'SR202500005F-1']
    
    for sample in csv_samples:
        csv_file_path = f"{csv_files_dir}/{sample}.csv"
        
        try:
            print(f"正在处理样本 {sample} 的CSV文件...")
            
            # 读取CSV文件，第一列是cell_id，第二列是stage
            stage_df = pd.read_csv(csv_file_path, header=None, names=['cell_id', 'stage'])
            
            # 创建cell_id到stage的映射字典
            cell_to_stage_map = dict(zip(stage_df['cell_id'], stage_df['stage']))
            
            # 获取当前样本的所有细胞行索引和对应的cell_id
            sample_cell_mask = adata.obs['sample'] == sample
            sample_cell_ids = adata.obs.loc[sample_cell_mask, 'cell_id']
            
            # 统计匹配情况
            matched_cells = 0
            total_cells = len(sample_cell_ids)
            
            # 为在CSV文件中存在的细胞设置stage
            for idx, cell_id in sample_cell_ids.items():
                if cell_id in cell_to_stage_map:
                    # 使用行索引idx来定位具体的行
                    adata.obs.loc[idx, 'stage'] = cell_to_stage_map[cell_id]
                    matched_cells += 1
            
            print(f"样本 {sample} 处理完成: {matched_cells}/{total_cells} 个细胞成功匹配")
            
        except FileNotFoundError:
            print(f"警告: 文件 {csv_file_path} 不存在，样本 {sample} 的stage将保持为'others'")
            continue
        except Exception as e:
            print(f"处理样本 {sample} 的CSV文件时出错: {e}")
            continue
    
    # 5. 验证结果
    print("\n" + "="*50)
    print("stage列赋值结果统计")
    print("="*50)
    
    stage_counts = adata.obs['stage'].value_counts()
    print("\n各stage的细胞数量:")
    for stage, count in stage_counts.items():
        print(f"  {stage}: {count} 个细胞")
    
    print("\n各样本的stage分布:")
    sample_stage_dist = pd.crosstab(adata.obs['sample'], adata.obs['stage'], margins=True)
    print(sample_stage_dist)
    
    # 6. 保存结果（可选）
    if output_path:
        print(f"\n正在保存结果到 {output_path}")
        adata.write_h5ad(output_path)
        print("保存完成")
    else:
        print("\n注意: 未指定输出路径，结果未保存到文件")
    
    return adata

# 使用示例
if __name__ == "__main__":
    # 设置文件路径
    h5ad_file = "adata_celltype.h5ad"  # 替换为您的h5ad文件路径
    csv_directory = "sample"  # 替换为CSV文件所在目录
    output_file = "adata_celltype_stage.h5ad"  # 输出文件路径（可选）
    
    # 执行函数
    try:
        result_adata = add_stage_column(h5ad_file, csv_directory, output_file)
        print("\n" + "="*50)
        print("处理完成!")
        print("="*50)
        
        # 打印前几行数据进行检查
        print("\n检查前几行数据:")
        print(result_adata.obs[['sample', 'cell_id', 'stage']].head())
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
