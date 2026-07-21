import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def add_region_to_h5ad_with_unmatched_samples(csv_file, h5ad_file, output_file=None):
    """
    读取CSV文件和h5ad数据，根据匹配情况添加region列
    并输出在CSV中但未匹配上的样本
    
    参数:
    - csv_file: CSV文件路径
    - h5ad_file: h5ad文件路径
    - output_file: 输出文件路径，如果不指定则覆盖原文件
    """
    
    # 1. 读取CSV文件
    print(f"正在读取CSV文件: {csv_file}")
    try:
        # 尝试不同的编码方式读取CSV文件
        try:
            csv_data = pd.read_csv(csv_file, encoding='utf-8-sig')
        except:
            try:
                csv_data = pd.read_csv(csv_file, encoding='utf-8')
            except:
                csv_data = pd.read_csv(csv_file, encoding='latin1')
        
        # 处理BOM字符
        if csv_data.columns[0].startswith('\ufeff'):
            csv_data.columns = [col.replace('\ufeff', '') for col in csv_data.columns]
        
        print(f"CSV数据形状: {csv_data.shape}")
        print(f"CSV列名: {csv_data.columns.tolist()}")
        
        # 重命名列，确保列名一致
        column_mapping = {}
        for col in csv_data.columns:
            if 'file' in col.lower() or 'sample' in col.lower() or 'filename' in col.lower():
                column_mapping[col] = 'sample'
            elif 'cell' in col.lower() and 'id' in col.lower():
                column_mapping[col] = 'cell_id'
        
        if column_mapping:
            csv_data = csv_data.rename(columns=column_mapping)
            print(f"重命名列: {column_mapping}")
        
        # 确保有sample和cell_id列
        if 'sample' not in csv_data.columns or 'cell_id' not in csv_data.columns:
            # 尝试使用前两列
            if len(csv_data.columns) >= 2:
                csv_data = csv_data.rename(columns={
                    csv_data.columns[0]: 'sample',
                    csv_data.columns[1]: 'cell_id'
                })
                print(f"使用前两列作为sample和cell_id")
            else:
                print("错误: CSV文件没有找到sample和cell_id列")
                return None, None
        
        # 清理数据
        csv_data['sample'] = csv_data['sample'].astype(str).str.strip()
        csv_data['cell_id'] = csv_data['cell_id'].astype(str).str.strip()
        
        # 删除可能的空值
        csv_data = csv_data.dropna(subset=['sample', 'cell_id'])
        
        print(f"CSV中tumor样本-细胞对数量: {len(csv_data)}")
        print(f"CSV中唯一的sample数量: {csv_data['sample'].nunique()}")
        
        # 显示前几行
        print("\nCSV文件前5行:")
        print(csv_data.head())
        
    except Exception as e:
        print(f"读取CSV文件失败: {e}")
        return None, None
    
    # 2. 读取h5ad文件
    print(f"\n正在读取h5ad文件: {h5ad_file}")
    try:
        adata = sc.read_h5ad(h5ad_file)
        print(f"h5ad数据形状: {adata.shape}")
        print(f"obs列: {adata.obs.columns.tolist()}")
        
        # 检查必需列
        if 'sample' not in adata.obs.columns or 'cell_id' not in adata.obs.columns:
            print("错误: h5ad数据缺少'sample'或'cell_id'列")
            print(f"可用列: {adata.obs.columns.tolist()}")
            return None, None
        
        # 清理h5ad中的sample和cell_id列
        adata.obs['sample'] = adata.obs['sample'].astype(str).str.strip()
        adata.obs['cell_id'] = adata.obs['cell_id'].astype(str).str.strip()
        
        print(f"h5ad中唯一的sample数量: {adata.obs['sample'].nunique()}")
        print(f"h5ad中sample列示例: {adata.obs['sample'].unique()[:10]}")
        print(f"h5ad中cell_id列示例: {adata.obs['cell_id'].unique()[:10]}")
        
    except Exception as e:
        print(f"读取h5ad文件失败: {e}")
        return None, None
    
    # 3. 创建匹配集合
    print(f"\n创建匹配集合...")
    
    # 创建CSV中的(sample, cell_id)对集合
    csv_pairs = set()
    for _, row in csv_data.iterrows():
        csv_sample = str(row['sample']).strip()
        csv_cell_id = str(row['cell_id']).strip()
        csv_pairs.add((csv_sample, csv_cell_id))
    
    print(f"CSV中唯一的(sample, cell_id)对数量: {len(csv_pairs)}")
    
    # 4. 检查样本名匹配情况
    print(f"\n检查样本名匹配情况...")
    
    h5ad_samples = set(adata.obs['sample'].unique())
    csv_samples = set(csv_data['sample'].unique())
    
    print(f"h5ad中sample数量: {len(h5ad_samples)}")
    print(f"CSV中sample数量: {len(csv_samples)}")
    
    # 查找完全匹配的样本
    exact_matches = h5ad_samples.intersection(csv_samples)
    print(f"完全匹配的样本数: {len(exact_matches)}")
    
    # 查找在CSV中但不在h5ad中的样本
    unmatched_csv_samples = csv_samples - h5ad_samples
    print(f"在CSV中但不在h5ad中的样本数: {len(unmatched_csv_samples)}")
    
    if unmatched_csv_samples:
        print("在CSV中但不在h5ad中的样本:")
        for i, sample in enumerate(sorted(unmatched_csv_samples)):
            print(f"  {i+1:3d}. {sample}")
        
        # 统计这些样本在CSV中的记录数
        unmatched_samples_data = []
        for sample in unmatched_csv_samples:
            count = len(csv_data[csv_data['sample'] == sample])
            unmatched_samples_data.append({
                'sample': sample,
                'cell_count_in_csv': count
            })
        
        unmatched_df = pd.DataFrame(unmatched_samples_data)
        print("\n在CSV中但不在h5ad中的样本详情:")
        print(unmatched_df)
    
    # 5. 添加region列
    print(f"\n添加region列...")
    
    # 初始化所有细胞为normal
    adata.obs['region'] = 'normal'
    
    # 统计匹配数量
    tumor_count = 0
    normal_count = 0
    
    # 用于记录每个样本的匹配情况
    sample_match_stats = {}
    
    # 初始化所有CSV样本的统计
    for sample in csv_samples:
        sample_match_stats[sample] = {
            'matched_cells': 0,
            'total_cells_in_csv': len(csv_data[csv_data['sample'] == sample])
        }
    
    # 方法1: 使用循环（适用于中等规模数据）
    print("使用循环进行匹配...")
    for idx, row in adata.obs.iterrows():
        h5ad_sample = str(row['sample']).strip()
        h5ad_cell_id = str(row['cell_id']).strip()
        
        # 检查是否在CSV中
        if (h5ad_sample, h5ad_cell_id) in csv_pairs:
            adata.obs.at[idx, 'region'] = 'tumor'
            tumor_count += 1
            
            # 记录匹配
            if h5ad_sample in sample_match_stats:
                sample_match_stats[h5ad_sample]['matched_cells'] += 1
        else:
            normal_count += 1
    
    print(f"匹配结果: tumor={tumor_count}, normal={normal_count}")
    
    # 6. 分析匹配结果
    print(f"\n分析匹配结果...")
    
    # 统计各样本的匹配情况
    sample_stats = adata.obs.groupby(['sample', 'region']).size().unstack(fill_value=0)
    
    # 找出在CSV中有匹配但匹配数为0的样本
    csv_samples_with_matches = []
    for sample in csv_samples:
        if sample in h5ad_samples:  # 样本在h5ad中存在
            tumor_cells = sample_stats.loc[sample, 'tumor'] if sample in sample_stats.index and 'tumor' in sample_stats.columns else 0
            csv_samples_with_matches.append({
                'sample': sample,
                'tumor_cells_in_h5ad': tumor_cells,
                'total_cells_in_csv': len(csv_data[csv_data['sample'] == sample]),
                'matched_ratio': tumor_cells / len(csv_data[csv_data['sample'] == sample]) if len(csv_data[csv_data['sample'] == sample]) > 0 else 0
            })
    
    csv_samples_df = pd.DataFrame(csv_samples_with_matches)
    
    # 找出在CSV中但在h5ad中没有匹配的样本
    zero_match_samples = csv_samples_df[csv_samples_df['tumor_cells_in_h5ad'] == 0]
    
    print(f"\n在CSV中有但在h5ad中没有匹配的样本 (匹配数为0): {len(zero_match_samples)}")
    if len(zero_match_samples) > 0:
        print("样本详情:")
        for i, row in zero_match_samples.iterrows():
            print(f"  样本: {row['sample']}, CSV中细胞数: {row['total_cells_in_csv']}")
    
    # 7. 保存结果
    if output_file is None:
        output_file = h5ad_file.replace('.h5ad', '_with_region.h5ad')
    
    print(f"\n保存结果到: {output_file}")
    adata.write_h5ad(output_file)
    
    # 8. 生成详细统计报告
    print(f"\n生成统计报告...")
    
    report_file = output_file.replace('.h5ad', '_report.txt')
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Region分类统计报告\n")
        f.write("="*60 + "\n\n")
        
        f.write("1. 输入文件信息\n")
        f.write("-"*40 + "\n")
        f.write(f"CSV文件: {csv_file}\n")
        f.write(f"  tumor样本-细胞对数量: {len(csv_data)}\n")
        f.write(f"  唯一sample数量: {csv_data['sample'].nunique()}\n")
        f.write(f"h5ad文件: {h5ad_file}\n")
        f.write(f"  细胞总数: {adata.shape[0]}\n")
        f.write(f"  唯一sample数量: {adata.obs['sample'].nunique()}\n\n")
        
        f.write("2. 样本名匹配情况\n")
        f.write("-"*40 + "\n")
        f.write(f"完全匹配的样本数: {len(exact_matches)}\n")
        f.write(f"在CSV中但不在h5ad中的样本数: {len(unmatched_csv_samples)}\n")
        
        if len(unmatched_csv_samples) > 0:
            f.write("\n在CSV中但不在h5ad中的样本:\n")
            for i, sample in enumerate(sorted(unmatched_csv_samples)):
                count = len(csv_data[csv_data['sample'] == sample])
                f.write(f"  {i+1:3d}. {sample} (在CSV中有{count}个细胞)\n")
        
        f.write("\n3. 匹配结果\n")
        f.write("-"*40 + "\n")
        f.write(f"tumor细胞数: {tumor_count} ({tumor_count/adata.shape[0]*100:.2f}%)\n")
        f.write(f"normal细胞数: {normal_count} ({normal_count/adata.shape[0]*100:.2f}%)\n\n")
        
        f.write("4. 在CSV中有但在h5ad中没有匹配的样本 (匹配数为0)\n")
        f.write("-"*40 + "\n")
        if len(zero_match_samples) > 0:
            f.write(f"共有 {len(zero_match_samples)} 个样本:\n")
            for i, row in zero_match_samples.iterrows():
                f.write(f"  {i+1:3d}. {row['sample']} (在CSV中有{row['total_cells_in_csv']}个细胞)\n")
        else:
            f.write("无\n")
        
        f.write("\n5. 各样本region分布\n")
        f.write("-"*40 + "\n")
        for sample in sorted(sample_stats.index):
            tumor_cells = sample_stats.loc[sample, 'tumor'] if 'tumor' in sample_stats.columns else 0
            normal_cells = sample_stats.loc[sample, 'normal'] if 'normal' in sample_stats.columns else 0
            total = tumor_cells + normal_cells
            
            if total > 0:
                tumor_pct = (tumor_cells / total) * 100
                f.write(f"{sample}: {tumor_cells}/{total} ({tumor_pct:.1f}%) tumor\n")
        
        f.write(f"\n6. 输出文件\n")
        f.write("-"*40 + "\n")
        f.write(f"h5ad文件: {output_file}\n")
        f.write(f"报告文件: {report_file}\n")
    
    print(f"报告已保存到: {report_file}")
    
    # 9. 保存未匹配样本到CSV文件
    if len(unmatched_csv_samples) > 0:
        unmatched_samples_file = output_file.replace('.h5ad', '_unmatched_samples.csv')
        unmatched_samples_data = []
        for sample in unmatched_csv_samples:
            count = len(csv_data[csv_data['sample'] == sample])
            unmatched_samples_data.append({
                'sample': sample,
                'cell_count_in_csv': count
            })
        
        unmatched_samples_df = pd.DataFrame(unmatched_samples_data)
        unmatched_samples_df.to_csv(unmatched_samples_file, index=False, encoding='utf-8-sig')
        print(f"未匹配样本列表已保存到: {unmatched_samples_file}")
    
    # 10. 保存零匹配样本到CSV文件
    if len(zero_match_samples) > 0:
        zero_match_file = output_file.replace('.h5ad', '_zero_match_samples.csv')
        zero_match_samples.to_csv(zero_match_file, index=False, encoding='utf-8-sig')
        print(f"零匹配样本列表已保存到: {zero_match_file}")
    
    return adata, zero_match_samples

# 主程序
if __name__ == "__main__":
    # 设置文件路径
    csv_file = "releate.csv"  # 你的CSV文件名
    h5ad_file = "All_filter.h5ad"  # 你的h5ad文件名
    output_file = "All_filter_with_region.h5ad"  # 输出文件名
    
    print("开始处理region分类...")
    print("="*60)
    
    try:
        # 使用完整版本
        adata_with_region, zero_match_samples = add_region_to_h5ad_with_unmatched_samples(
            csv_file=csv_file,
            h5ad_file=h5ad_file,
            output_file=output_file
        )
        
        if adata_with_region is not None:
            print(f"\n处理完成!")
            print(f"原始数据: {adata_with_region.shape}")
            print(f"region分布:")
            print(adata_with_region.obs['region'].value_counts())
            
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
