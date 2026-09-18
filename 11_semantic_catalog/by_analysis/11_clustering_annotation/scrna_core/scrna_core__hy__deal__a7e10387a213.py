import os
import csv
import glob

def process_csv_files():
    """
    处理所有CSV文件，提取信息并整合到state.csv
    输出格式：
    cell_id, sample, group
    """
    
    # 获取当前目录下所有CSV文件
    csv_files = glob.glob("*.csv")
    
    if not csv_files:
        print("当前目录下没有找到CSV文件")
        return
    
    print(f"找到 {len(csv_files)} 个CSV文件")
    print("开始处理...")
    
    all_data = []
    total_cells = 0
    
    for csv_file in csv_files:
        # 提取样本名（完整文件名）
        sample_name = os.path.basename(csv_file)
        
        # 提取组名：用"_"分割的第一个字母
        group_name = sample_name.split('_')[0]
        
        print(f"处理文件: {csv_file}")
        print(f"  组名: {group_name}")
        
        cell_count = 0
        
        try:
            with open(csv_file, 'r', encoding='utf-8') as file:
                # 读取前3行（注释和标题）
                lines = file.readlines()
                
                if len(lines) < 4:
                    print(f"  警告: 文件 {csv_file} 行数不足4行，跳过")
                    continue
                
                # 验证文件格式
                if not lines[0].startswith("#Selection name:"):
                    print(f"  警告: 文件 {csv_file} 格式不符合预期，第一行不是注释")
                
                if not lines[2].startswith("Cell ID"):
                    print(f"  警告: 文件 {csv_file} 第三行不是标题行")
                
                # 从第4行开始处理（索引3）
                for i in range(3, len(lines)):
                    line = lines[i].strip()
                    if not line:  # 跳过空行
                        continue
                    
                    # 分割行，使用逗号分隔
                    parts = line.split(',')
                    if len(parts) >= 1:
                        cell_id = parts[0].strip()
                        if cell_id:  # 确保cell_id不为空
                            # 添加到数据列表
                            all_data.append([cell_id, sample_name, group_name])
                            cell_count += 1
                            total_cells += 1
        
        except Exception as e:
            print(f"  处理文件 {csv_file} 时出错: {e}")
            continue
        
        print(f"  提取了 {cell_count} 个细胞ID")
    
    # 将所有数据写入state.csv文件
    if all_data:
        output_file = "state.csv"
        print(f"\n正在写入 {len(all_data)} 条记录到 {output_file}...")
        
        try:
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # 写入标题行
                writer.writerow(["cell_id", "sample", "group"])
                # 写入数据
                writer.writerows(all_data)
            
            print(f"成功创建: {output_file}")
            print(f"总计处理: {len(csv_files)} 个文件")
            print(f"总计细胞: {total_cells} 个")
            
            # 显示统计信息
            print("\n统计信息:")
            group_stats = {}
            for _, _, group in all_data:
                group_stats[group] = group_stats.get(group, 0) + 1
            
            for group, count in group_stats.items():
                print(f"  组 {group}: {count} 个细胞")
            
        except Exception as e:
            print(f"写入文件 {output_file} 时出错: {e}")
    else:
        print("没有提取到任何数据")

# 另一个版本：使用pandas库（更简洁）
def process_csv_files_with_pandas():
    """
    使用pandas库处理CSV文件
    需要先安装pandas: pip install pandas
    """
    try:
        import pandas as pd
    except ImportError:
        print("pandas库未安装，使用基本方法处理")
        process_csv_files()
        return
    
    # 获取当前目录下所有CSV文件
    csv_files = glob.glob("*.csv")
    
    if not csv_files:
        print("当前目录下没有找到CSV文件")
        return
    
    print(f"找到 {len(csv_files)} 个CSV文件")
    print("开始处理...")
    
    all_dfs = []
    
    for csv_file in csv_files:
        # 提取样本名（完整文件名）
        sample_name = os.path.basename(csv_file)
        
        # 提取组名：用"_"分割的第一个字母
        group_name = sample_name.split('_')[0]
        
        print(f"处理文件: {csv_file}")
        
        try:
            # 使用pandas读取CSV，跳过前3行
            df = pd.read_csv(csv_file, skiprows=3)
            
            # 确保有Cell ID列
            if 'Cell ID' in df.columns:
                # 创建新DataFrame
                new_df = pd.DataFrame({
                    'cell_id': df['Cell ID'],
                    'sample': sample_name,
                    'group': group_name
                })
                
                all_dfs.append(new_df)
                print(f"  提取了 {len(new_df)} 个细胞ID")
            else:
                print(f"  警告: 文件 {csv_file} 没有'Cell ID'列")
                
        except Exception as e:
            print(f"  处理文件 {csv_file} 时出错: {e}")
    
    # 合并所有DataFrame
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        
        # 保存到CSV文件
        output_file = "state.csv"
        combined_df.to_csv(output_file, index=False)
        
        print(f"\n成功创建: {output_file}")
        print(f"总计处理: {len(csv_files)} 个文件")
        print(f"总计细胞: {len(combined_df)} 个")
        
        # 显示统计信息
        print("\n统计信息:")
        group_stats = combined_df['group'].value_counts()
        for group, count in group_stats.items():
            print(f"  组 {group}: {count} 个细胞")
        
        # 显示前几行数据
        print("\n输出文件前5行:")
        print(combined_df.head())
    else:
        print("没有提取到任何数据")

# 创建一个更详细的版本，包含文件验证
def process_csv_files_detailed():
    """
    详细处理CSV文件，包含文件格式验证
    """
    
    # 获取当前目录下所有CSV文件
    csv_files = glob.glob("*.csv")
    
    if not csv_files:
        print("当前目录下没有找到CSV文件")
        return
    
    print("=" * 60)
    print(f"找到 {len(csv_files)} 个CSV文件")
    print("=" * 60)
    
    all_data = []
    total_cells = 0
    processed_files = 0
    
    for csv_file in csv_files:
        # 提取样本名（完整文件名）
        sample_name = os.path.basename(csv_file)
        
        # 提取组名：用"_"分割的第一个字母
        group_name = sample_name.split('_')[0]
        
        print(f"\n[{processed_files + 1}/{len(csv_files)}] 处理文件: {sample_name}")
        print(f"  组名: {group_name}")
        
        cell_count = 0
        
        try:
            with open(csv_file, 'r', encoding='utf-8') as file:
                # 读取前3行
                lines = []
                for _ in range(3):
                    lines.append(file.readline().strip())
                
                # 验证文件格式
                if len(lines) < 3:
                    print(f"  错误: 文件 {csv_file} 行数不足3行，跳过")
                    continue
                
                # 检查注释行
                if not lines[0].startswith("#Selection name:"):
                    print(f"  警告: 文件 {csv_file} 第一行不是'#Selection name:'")
                
                if not lines[1].startswith("#Area"):
                    print(f"  警告: 文件 {csv_file} 第二行不是区域面积注释")
                
                # 检查标题行
                expected_headers = ["Cell ID", "Cluster", "Transcripts", "Area (µm^2)"]
                actual_headers = [h.strip() for h in lines[2].split(',')]
                
                if actual_headers[:len(expected_headers)] != expected_headers:
                    print(f"  警告: 文件 {csv_file} 标题行不符合预期")
                    print(f"    预期: {expected_headers}")
                    print(f"    实际: {actual_headers}")
                
                # 处理数据行
                for line_num, line in enumerate(file, start=4):  # 从第4行开始
                    line = line.strip()
                    if not line:  # 跳过空行
                        continue
                    
                    parts = line.split(',')
                    if len(parts) >= 1:
                        cell_id = parts[0].strip()
                        if cell_id:  # 确保cell_id不为空
                            all_data.append([cell_id, sample_name, group_name])
                            cell_count += 1
        
        except FileNotFoundError:
            print(f"  错误: 文件 {csv_file} 不存在")
            continue
        except Exception as e:
            print(f"  处理文件 {csv_file} 时出错: {e}")
            continue
        
        print(f"  提取了 {cell_count} 个细胞ID")
        total_cells += cell_count
        processed_files += 1
    
    # 写入输出文件
    if all_data:
        output_file = "state.csv"
        print(f"\n" + "=" * 60)
        print(f"正在写入 {len(all_data)} 条记录到 {output_file}")
        print("=" * 60)
        
        try:
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # 写入标题行
                writer.writerow(["cell_id", "sample", "group"])
                # 写入数据
                writer.writerows(all_data)
            
            print(f"✓ 成功创建: {output_file}")
            print(f"✓ 处理文件: {processed_files}/{len(csv_files)} 个")
            print(f"✓ 总计细胞: {total_cells} 个")
            
            # 显示详细统计信息
            print("\n" + "=" * 60)
            print("详细统计信息:")
            print("=" * 60)
            
            # 按组统计
            group_stats = {}
            for _, _, group in all_data:
                group_stats[group] = group_stats.get(group, 0) + 1
            
            for group, count in group_stats.items():
                percentage = (count / total_cells) * 100
                print(f"  组 {group}: {count} 个细胞 ({percentage:.1f}%)")
            
            # 按样本统计
            sample_stats = {}
            for _, sample, _ in all_data:
                sample_stats[sample] = sample_stats.get(sample, 0) + 1
            
            print(f"\n按样本统计 (前10个):")
            sorted_samples = sorted(sample_stats.items(), key=lambda x: x[1], reverse=True)
            for sample, count in sorted_samples[:10]:
                print(f"  {sample}: {count} 个细胞")
            
            if len(sorted_samples) > 10:
                print(f"  ... 还有 {len(sorted_samples) - 10} 个样本")
            
            # 显示输出文件的前几行
            print(f"\n输出文件 {output_file} 的前5行:")
            with open(output_file, 'r', encoding='utf-8') as f:
                for i, line in enumerate(f):
                    if i < 6:  # 包括标题行
                        print(f"  {line.strip()}")
                    else:
                        break
            
        except Exception as e:
            print(f"写入文件 {output_file} 时出错: {e}")
    else:
        print("没有提取到任何数据")

# 创建辅助函数：验证单个文件
def validate_csv_file(csv_file):
    """验证CSV文件格式"""
    try:
        with open(csv_file, 'r', encoding='utf-8') as file:
            lines = []
            for _ in range(3):
                lines.append(file.readline().strip())
            
            if len(lines) < 3:
                return False, "文件行数不足3行"
            
            if not lines[0].startswith("#Selection name:"):
                return False, "第一行不是'#Selection name:'"
            
            if not lines[1].startswith("#Area"):
                return False, "第二行不是区域面积注释"
            
            expected_headers = ["Cell ID", "Cluster", "Transcripts", "Area (µm^2)"]
            actual_headers = [h.strip() for h in lines[2].split(',')]
            
            if actual_headers[:len(expected_headers)] != expected_headers:
                return False, f"标题行不符合预期: {actual_headers}"
            
            return True, "格式正确"
    
    except Exception as e:
        return False, f"读取文件出错: {e}"

if __name__ == "__main__":
    print("=" * 60)
    print("CSV文件处理工具")
    print("功能：提取所有CSV文件的cell_id，生成state.csv")
    print("=" * 60)
    
    # 使用方法1：基本方法
    process_csv_files()
    
    # 使用方法2：使用pandas（如果已安装）
    # process_csv_files_with_pandas()
    
    # 使用方法3：详细版本
    # process_csv_files_detailed()
    
    print("\n" + "=" * 60)
    print("处理完成！")
    print("=" * 60)
