import pandas as pd

def merge_cnv_scores(a_file, cell_info_file, output_file=None):
    """
    将a文件中的cnv_score合并到cell_info_res0.4.csv文件中
    
    参数:
    a_file: 包含cell, cluster, cnv_score三列的文件路径
    cell_info_file: 包含cell_id和res_0.4两列的文件路径  
    output_file: 输出文件路径，默认为None会覆盖原文件
    """
    # 读取两个文件
    try:
        # 读取a文件，假设是逗号分隔
        df_a = pd.read_csv(a_file)
        print(f"成功读取{a_file}，包含{len(df_a)}行数据")
        
        # 读取cell_info文件
        df_cell_info = pd.read_csv(cell_info_file)
        print(f"成功读取{cell_info_file}，包含{len(df_cell_info)}行数据")
        
    except FileNotFoundError as e:
        print(f"文件读取错误: {e}")
        return
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return
    
    # 检查必要的列是否存在
    required_a_columns = ['cell', 'cnv_score']
    required_cell_info_columns = ['cell_id']
    
    missing_a_columns = [col for col in required_a_columns if col not in df_a.columns]
    missing_cell_info_columns = [col for col in required_cell_info_columns if col not in df_cell_info.columns]
    
    if missing_a_columns:
        print(f"错误: a文件中缺少必要的列: {missing_a_columns}")
        print(f"a文件中的列: {list(df_a.columns)}")
        return
        
    if missing_cell_info_columns:
        print(f"错误: cell_info文件中缺少必要的列: {missing_cell_info_columns}")
        print(f"cell_info文件中的列: {list(df_cell_info.columns)}")
        return
    
    # 显示数据预览
    print("\n数据预览:")
    print("a文件前5行:")
    print(df_a.head())
    print("\ncell_info文件前5行:")
    print(df_cell_info.head())
    
    # 合并数据 - 使用左连接，保留cell_info中的所有行
    merged_df = df_cell_info.merge(
        df_a[['cell', 'cnv_score']], 
        left_on='cell_id', 
        right_on='cell', 
        how='left'
    )
    
    # 删除重复的cell列（来自a文件）
    if 'cell' in merged_df.columns:
        merged_df = merged_df.drop('cell', axis=1)
    
    print(f"\n合并后数据形状: {merged_df.shape}")
    print(f"成功匹配的细胞数量: {merged_df['cnv_score'].notna().sum()}")
    print(f"未匹配的细胞数量: {merged_df['cnv_score'].isna().sum()}")
    
    # 检查是否有重复的cell_id
    duplicate_cells = merged_df[merged_df.duplicated('cell_id', keep=False)]
    if len(duplicate_cells) > 0:
        print(f"警告: 发现{len(duplicate_cells)}个重复的cell_id")
        print("重复的cell_id示例:")
        print(duplicate_cells['cell_id'].value_counts().head())
    
    # 确定输出文件路径
    if output_file is None:
        output_file = cell_info_file  # 覆盖原文件
        print(f"将覆盖原文件: {cell_info_file}")
    else:
        print(f"将保存为新文件: {output_file}")
    
    # 保存结果
    try:
        merged_df.to_csv(output_file, index=False)
        print(f"\n合并完成! 结果已保存到: {output_file}")
        
        # 显示合并结果的统计信息
        print("\n合并结果统计:")
        print(f"总细胞数: {len(merged_df)}")
        print(f"具有cnv_score的细胞数: {merged_df['cnv_score'].notna().sum()}")
        if 'res_0.4' in merged_df.columns:
            print(f"各cluster的细胞分布:")
            print(merged_df['res_0.4'].value_counts().sort_index())
        
    except Exception as e:
        print(f"保存文件时发生错误: {e}")
        return
    
    return merged_df

# 使用示例
if __name__ == "__main__":
    # 文件路径
    a_file_path = "a"  # 您的a文件路径
    cell_info_file_path = "cell_info_res0.4.csv"  # 您的cell_info文件路径
    
    # 执行合并
    result_df = merge_cnv_scores(
        a_file=a_file_path,
        cell_info_file=cell_info_file_path,
        output_file="cell_info_with_cnv_0.4.csv"  # 可选: 指定输出文件名，不指定则覆盖原文件
    )
    
    if result_df is not None:
        print("\n合并后的数据前5行:")
        print(result_df.head())
