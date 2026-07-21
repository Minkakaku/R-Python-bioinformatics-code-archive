import pandas as pd
import os
import glob
from tqdm import tqdm  # 进度条库，可选安装

def extract_cell_ids_from_csv_files():
    """
    从多个CSV文件中提取文件名和Cell ID
    """
    # 查找所有CSV文件
    csv_files = glob.glob("*.csv")
    
    if not csv_files:
        print("在当前目录下没有找到CSV文件")
        return None
    
    print(f"找到 {len(csv_files)} 个CSV文件")
    
    # 存储结果
    all_data = []
    
    # 处理每个CSV文件
    for csv_file in tqdm(csv_files, desc="处理CSV文件中"):
        try:
            # 读取CSV文件，跳过以#开头的注释行
            # 尝试不同的编码格式
            try:
                df = pd.read_csv(csv_file, comment='#')
            except UnicodeDecodeError:
                # 如果默认编码失败，尝试其他编码
                try:
                    df = pd.read_csv(csv_file, comment='#', encoding='latin1')
                except:
                    df = pd.read_csv(csv_file, comment='#', encoding='cp1252')
            
            # 检查是否有'Cell ID'列
            if 'Cell ID' in df.columns:
                cell_ids = df['Cell ID'].tolist()
            elif 'Cell_ID' in df.columns:
                cell_ids = df['Cell_ID'].tolist()
            elif 'cell_id' in df.columns:
                cell_ids = df['cell_id'].tolist()
            else:
                # 尝试第一列
                first_col = df.columns[0]
                cell_ids = df[first_col].tolist()
                
                # 检查第一列的值是否包含Cell ID的特征
                sample_ids = df[first_col].head(5).astype(str).tolist()
                if all(any(char in id_str for char in ['-', '_']) for id_str in sample_ids):
                    print(f"  警告: 文件 {csv_file} 使用第一列作为Cell ID")
                else:
                    print(f"  警告: 文件 {csv_file} 没有找到明确的Cell ID列")
                    continue
            
            # 为每个Cell ID添加文件名
            for cell_id in cell_ids:
                all_data.append({
                    'Filename': csv_file,
                    'Cell_ID': str(cell_id)
                })
            
        except Exception as e:
            print(f"  处理文件 {csv_file} 时出错: {e}")
    
    # 转换为DataFrame
    if all_data:
        result_df = pd.DataFrame(all_data)
        return result_df
    else:
        print("没有提取到任何数据")
        return None

# 主程序
if __name__ == "__main__":
    print("开始处理CSV文件...")
    
    # 提取数据
    result_df = extract_cell_ids_from_csv_files()
    
    if result_df is not None:
        # 统计信息
        print(f"\n处理完成!")
        print(f"总记录数: {len(result_df)}")
        print(f"处理的文件数: {result_df['Filename'].nunique()}")
        print(f"平均每个文件的Cell ID数: {len(result_df)/result_df['Filename'].nunique():.1f}")
        
        # 查看前几行
        print("\n前10行数据:")
        print(result_df.head(10))
        
        # 保存结果
        output_file = "filename_cellid_mapping.csv"
        result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\n结果已保存到: {output_file}")
        
        # 保存为TSV格式
        tsv_file = "filename_cellid_mapping.tsv"
        result_df.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8-sig')
        print(f"结果已保存到: {tsv_file}")
        
        # 保存为Excel格式
        excel_file = "filename_cellid_mapping.xlsx"
        result_df.to_excel(excel_file, index=False)
        print(f"结果已保存到: {excel_file}")
        
        # 生成统计报告
        print("\n统计报告:")
        print("="*50)
        
        # 按文件名统计
        file_stats = result_df.groupby('Filename').size()
        print(f"\n文件统计 (前10个文件):")
        for filename, count in file_stats.head(10).items():
            print(f"  {filename}: {count} 个Cell ID")
        
        if len(file_stats) > 10:
            print(f"  ... 还有 {len(file_stats)-10} 个文件")
        
        # 总统计
        print(f"\n总计:")
        print(f"  Cell ID总数: {len(result_df)}")
        print(f"  文件总数: {result_df['Filename'].nunique()}")
        print(f"  平均每个文件: {len(result_df)/result_df['Filename'].nunique():.1f} 个Cell ID")
        print(f"  最大文件: {file_stats.idxmax()} ({file_stats.max()} 个Cell ID)")
        print(f"  最小文件: {file_stats.idxmin()} ({file_stats.min()} 个Cell ID)")
        
        # 保存统计信息
        stats_df = pd.DataFrame({
            'Filename': file_stats.index,
            'Cell_Count': file_stats.values
        }).sort_values('Cell_Count', ascending=False)
        
        stats_df.to_csv("file_statistics.csv", index=False, encoding='utf-8-sig')
        print(f"\n文件统计已保存到: file_statistics.csv")
