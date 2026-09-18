import os
import pandas as pd
import glob
from pathlib import Path

def filter_csv_by_cell_ids(id_file_path, bak_files_pattern='*.csv.bak'):
    """
    从.csv.bak文件中过滤出出现在id文件中的cell_id行，并追加保存为.csv文件
    
    参数:
    ----------
    id_file_path : str
        包含cell_id的文本文件路径
    bak_files_pattern : str
        .csv.bak文件的匹配模式，默认为'*.csv.bak'
    """
    
    # 1. 读取id文件中的所有cell_id
    cell_ids_set = set()
    
    id_file = Path(id_file_path)
    if not id_file.exists():
        print(f"错误: id文件不存在: {id_file_path}")
        return
    
    try:
        with open(id_file, 'r') as f:
            for line in f:
                cell_id = line.strip()
                if cell_id:  # 跳过空行
                    cell_ids_set.add(cell_id)
    except Exception as e:
        print(f"读取id文件 {id_file} 时出错: {e}")
        return
    
    print(f"从id文件中收集了 {len(cell_ids_set)} 个唯一的cell_id")
    
    if len(cell_ids_set) == 0:
        print("警告: 没有读取到任何cell_id")
        return
    
    # 2. 查找所有的.csv.bak文件
    bak_files = list(Path('.').glob(bak_files_pattern))
    
    if len(bak_files) == 0:
        print(f"没有找到匹配 {bak_files_pattern} 的文件")
        return
    
    print(f"找到了 {len(bak_files)} 个.csv.bak文件")
    
    # 3. 处理每个.csv.bak文件
    processed_count = 0
    for bak_file in bak_files:
        try:
            # 生成输出文件名（去掉.bak后缀）
            output_file = bak_file.with_suffix('')  # 移除.bak后缀
            
            # 读取.bak文件的前两行注释
            with open(bak_file, 'r') as f:
                # 读取前两行注释
                comment_lines = []
                for _ in range(2):
                    line = f.readline()
                    if line.startswith('#'):
                        comment_lines.append(line)
                    else:
                        # 如果不是注释行，放回去
                        f.seek(0)
                        break
                
                # 读取数据
                df = pd.read_csv(f)
            
            print(f"处理文件: {bak_file.name}")
            print(f"  原始行数: {len(df)}")
            
            # 检查是否有'Cell ID'列
            if 'Cell ID' not in df.columns:
                # 尝试查找大小写不敏感的列名
                cell_id_col = None
                for col in df.columns:
                    if col.lower() == 'cell id':
                        cell_id_col = col
                        break
                
                if cell_id_col is None:
                    print(f"  警告: 文件 {bak_file.name} 中没有找到'Cell ID'列，跳过")
                    continue
                
                # 重命名为标准列名
                df = df.rename(columns={cell_id_col: 'Cell ID'})
            
            # 过滤出在cell_ids_set中的行
            filtered_df = df[df['Cell ID'].isin(cell_ids_set)]
            
            print(f"  过滤后行数: {len(filtered_df)}")
            
            if len(filtered_df) > 0:
                # 检查输出文件是否存在
                file_exists = output_file.exists()
                
                # 追加模式写入
                mode = 'a' if file_exists else 'w'
                with open(output_file, mode) as f:
                    # 如果文件不存在，写入注释行和表头
                    if not file_exists:
                        for comment in comment_lines:
                            f.write(comment)
                        # 写入表头
                        filtered_df.to_csv(f, index=False)
                    else:
                        # 如果文件已存在，只追加数据，不写表头
                        filtered_df.to_csv(f, index=False, header=False)
                
                if file_exists:
                    print(f"  已追加到现有文件: {output_file.name}")
                else:
                    print(f"  已创建新文件: {output_file.name}")
                
                processed_count += 1
            else:
                print(f"  警告: 没有匹配的cell_id，跳过")
        
        except Exception as e:
            print(f"处理文件 {bak_file.name} 时出错: {e}")
    
    print(f"\n处理完成! 共处理了 {processed_count} 个文件")

# 高级版本：支持更多功能
def filter_csv_by_cell_ids_advanced(id_file_path, bak_files_pattern='*.csv.bak', 
                                    output_dir=None, id_column='Cell ID', 
                                    skip_comments=2, verbose=True, 
                                    encoding='utf-8', append=True):
    """
    高级版本的过滤函数，支持更多选项
    
    参数:
    ----------
    id_file_path : str
        包含cell_id的文本文件路径
    bak_files_pattern : str
        .csv.bak文件的匹配模式
    output_dir : str, optional
        输出目录，如果为None则输出到当前目录
    id_column : str, optional
        cell_id列的列名，默认为'Cell ID'
    skip_comments : int, optional
        跳过的注释行数，默认为2
    verbose : bool, optional
        是否显示详细信息，默认为True
    encoding : str, optional
        文件编码，默认为'utf-8'
    append : bool, optional
        是否追加到现有文件，默认为True
    """
    
    # 创建输出目录
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = Path('.')
    
    # 1. 读取id文件中的所有cell_id
    cell_ids_set = set()
    id_file = Path(id_file_path)
    
    if not id_file.exists():
        print(f"错误: id文件不存在: {id_file_path}")
        return
    
    try:
        with open(id_file, 'r', encoding=encoding) as f:
            for line in f:
                cell_id = line.strip()
                if cell_id and not cell_id.startswith('#'):  # 跳过空行和注释
                    cell_ids_set.add(cell_id)
    except Exception as e:
        print(f"读取id文件 {id_file} 时出错: {e}")
        return
    
    if verbose:
        print(f"从id文件中收集了 {len(cell_ids_set)} 个唯一的cell_id")
    
    if len(cell_ids_set) == 0:
        print("错误: 没有读取到任何cell_id")
        return
    
    # 2. 查找所有的.csv.bak文件
    bak_files = list(Path('.').rglob(bak_files_pattern))
    
    if len(bak_files) == 0:
        print(f"没有找到匹配 {bak_files_pattern} 的文件")
        return
    
    if verbose:
        print(f"找到了 {len(bak_files)} 个.csv.bak文件")
    
    # 3. 处理每个.csv.bak文件
    stats = {
        'total_files': len(bak_files),
        'processed': 0,
        'skipped': 0,
        'total_rows_before': 0,
        'total_rows_after': 0,
        'appended_to': 0,
        'created_new': 0
    }
    
    for bak_file in bak_files:
        try:
            # 生成输出文件名
            if output_dir:
                output_file = output_path / bak_file.with_suffix('').name
            else:
                output_file = bak_file.with_suffix('')
            
            # 检查输出文件是否存在
            file_exists = output_file.exists()
            
            # 如果文件存在且不允许追加，则跳过
            if file_exists and not append:
                if verbose:
                    print(f"跳过已存在的文件: {output_file.name}")
                stats['skipped'] += 1
                continue
            
            # 读取文件
            with open(bak_file, 'r', encoding=encoding) as f:
                # 读取注释行
                comment_lines = []
                for _ in range(skip_comments):
                    line = f.readline()
                    if line.startswith('#'):
                        comment_lines.append(line)
                    else:
                        # 如果不是注释行，放回去
                        f.seek(0)
                        comment_lines = []
                        break
                
                # 读取数据
                try:
                    df = pd.read_csv(f)
                except pd.errors.EmptyDataError:
                    print(f"  警告: 文件 {bak_file.name} 为空，跳过")
                    stats['skipped'] += 1
                    continue
            
            stats['total_rows_before'] += len(df)
            
            if verbose:
                print(f"处理文件: {bak_file.name}")
                print(f"  原始行数: {len(df)}")
            
            # 查找cell_id列
            actual_id_column = None
            for col in df.columns:
                if col.lower() == id_column.lower():
                    actual_id_column = col
                    break
            
            if actual_id_column is None:
                # 尝试其他可能的列名
                possible_names = ['Cell ID', 'cell_id', 'Cell_ID', 'cell id', 'cellid']
                for col in df.columns:
                    if any(name in col.lower() for name in [p.lower() for p in possible_names]):
                        actual_id_column = col
                        break
            
            if actual_id_column is None:
                print(f"  警告: 文件 {bak_file.name} 中没有找到'{id_column}'列，可用列: {list(df.columns)}")
                stats['skipped'] += 1
                continue
            
            # 重命名为标准列名以便后续处理
            df_renamed = df.rename(columns={actual_id_column: 'cell_id_temp'})
            
            # 过滤出在cell_ids_set中的行
            filtered_df = df_renamed[df_renamed['cell_id_temp'].isin(cell_ids_set)]
            
            # 恢复原始列名
            filtered_df = filtered_df.rename(columns={'cell_id_temp': actual_id_column})
            
            stats['total_rows_after'] += len(filtered_df)
            
            if verbose:
                print(f"  过滤后行数: {len(filtered_df)}")
                print(f"  匹配率: {len(filtered_df)/len(df)*100:.1f}%")
            
            if len(filtered_df) > 0:
                # 写入新的.csv文件
                mode = 'a' if file_exists and append else 'w'
                
                with open(output_file, mode, encoding=encoding) as f:
                    # 如果是新文件，写入注释行和表头
                    if mode == 'w':
                        for comment in comment_lines:
                            f.write(comment)
                        filtered_df.to_csv(f, index=False)
                        stats['created_new'] += 1
                    else:
                        # 追加模式，不写表头
                        filtered_df.to_csv(f, index=False, header=False)
                        stats['appended_to'] += 1
                
                if verbose:
                    if mode == 'w':
                        print(f"  已创建新文件: {output_file}")
                    else:
                        print(f"  已追加到现有文件: {output_file}")
                
                stats['processed'] += 1
            else:
                if verbose:
                    print(f"  警告: 没有匹配的cell_id，跳过")
                stats['skipped'] += 1
        
        except Exception as e:
            print(f"处理文件 {bak_file.name} 时出错: {e}")
            stats['skipped'] += 1
    
    # 打印统计信息
    print("\n" + "="*50)
    print("处理完成!")
    print(f"总文件数: {stats['total_files']}")
    print(f"成功处理: {stats['processed']}")
    print(f"  新创建文件: {stats['created_new']}")
    print(f"  追加到现有文件: {stats['appended_to']}")
    print(f"跳过/失败: {stats['skipped']}")
    print(f"总行数 (处理前): {stats['total_rows_before']}")
    print(f"总行数 (处理后): {stats['total_rows_after']}")
    if stats['total_rows_before'] > 0:
        print(f"总体保留率: {stats['total_rows_after']/stats['total_rows_before']*100:.1f}%")
    print("="*50)

# 使用示例
if __name__ == "__main__":
    import sys
    
    # 简单用法示例
    if len(sys.argv) > 1:
        id_file = sys.argv[1]
        
        if len(sys.argv) > 2:
            pattern = sys.argv[2]
        else:
            pattern = '*.csv.bak'
        
        filter_csv_by_cell_ids(id_file, pattern)
    else:
        # 如果没有命令行参数，使用硬编码的示例
        print("使用示例:")
        print("1. 基本用法（追加模式）:")
        print("   filter_csv_by_cell_ids('cell_ids.txt')")
        print("2. 高级用法（指定输出目录）:")
        print("   filter_csv_by_cell_ids_advanced('cell_ids.txt', output_dir='./filtered')")
        print("3. 高级用法（覆盖模式）:")
        print("   filter_csv_by_cell_ids_advanced('cell_ids.txt', append=False)")
        
        # 示例调用
        id_file_path = "all_id"  # 修改为您的id文件路径
        
        if Path(id_file_path).exists():
            print(f"\n尝试处理文件: {id_file_path}")
            filter_csv_by_cell_ids(id_file_path)
        else:
            print(f"\n示例文件不存在: {id_file_path}")
            print("请创建一个包含cell_id的文件，例如:")
            print("khoemamh-1")
            print("khoemamh-2")
            print("khoemamh-3")
