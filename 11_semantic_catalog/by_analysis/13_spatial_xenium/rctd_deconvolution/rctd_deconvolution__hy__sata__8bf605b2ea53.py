import scanpy as sc
import pandas as pd
import numpy as np

def analyze_celltype_distribution(h5ad_file, annotation_categories=None):
    """
    统计指定annotation类别的细胞在celltype列上的分布及占比
    
    参数:
    h5ad_file: h5ad文件路径
    annotation_categories: 需要统计的annotation类别列表，如果为None则使用默认列表
    """
    # 如果未提供annotation类别，使用默认列表
    if annotation_categories is None:
        annotation_categories = [
            'Unknown', 'Endo-CAF', 'B-CAF', 'Endo-Immu', 'B-T', 
            'T-Mast', 'Plasma-CAF', 'B-Epi', 'T-Epi', 'T-CAF'
        ]
    
    # 读取h5ad文件
    print("读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"数据形状: {adata.shape}")
    print(f"obs列名: {list(adata.obs.columns)}")
    
    # 检查必要的列是否存在
    if 'annotation' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中不存在'annotation'列")
    
    if 'celltype' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中不存在'celltype'列")
    
    # 获取总细胞数
    total_cells = adata.n_obs
    print(f"\n总细胞数: {total_cells}")
    
    # 检查数据中实际存在的annotation类别
    existing_categories = [cat for cat in annotation_categories if cat in adata.obs['annotation'].values]
    missing_categories = [cat for cat in annotation_categories if cat not in adata.obs['annotation'].values]
    
    if missing_categories:
        print(f"以下annotation类别在数据中不存在: {missing_categories}")
        print(f"将统计存在的类别: {existing_categories}")
    
    if not existing_categories:
        print("没有找到指定的annotation类别")
        return None
    
    # 创建结果存储字典
    results = {}
    
    print("\n" + "="*80)
    print("详细统计结果:")
    print("="*80)
    
    # 对每个存在的annotation类别进行统计
    for category in existing_categories:
        # 筛选当前annotation类别的细胞
        category_mask = adata.obs['annotation'] == category
        category_cells = adata.obs[category_mask]
        category_count = len(category_cells)
        
        # 计算该annotation类别在整体中的占比（使用round函数而不是方法）
        category_percent_total = round(category_count / total_cells * 100, 2)
        
        # 统计celltype分布
        celltype_dist = category_cells['celltype'].value_counts(dropna=False)
        
        # 计算每个celltype在当前annotation类别中的占比（使用round函数）
        if category_count > 0:
            # 使用列表推导式和round函数，避免对单个float使用.round()方法
            celltype_percentages = [round(count / category_count * 100, 2) for count in celltype_dist]
            celltype_percentages_series = pd.Series(celltype_percentages, index=celltype_dist.index)
            
            # 计算每个celltype在总细胞中的占比
            total_percentages = [round(count / total_cells * 100, 2) for count in celltype_dist]
        else:
            celltype_percentages_series = pd.Series([], dtype=float)
            total_percentages = []
        
        # 创建统计结果DataFrame
        stats_df = pd.DataFrame({
            '细胞数量': celltype_dist,
            f'在{category}中的占比(%)': celltype_percentages_series,
            '在总细胞中的占比(%)': total_percentages
        })
        
        # 存储结果
        results[category] = {
            'stats': stats_df,
            'total_count': category_count,
            'percent_total': category_percent_total
        }
        
        # 打印当前annotation类别的统计结果
        print(f"\n Annotation类别: {category}")
        print(f"   细胞数量: {category_count}个 ({category_percent_total:.2f}% 的总细胞)")
        print("-" * 60)
        
        if len(celltype_dist) > 0:
            # 打印详细的celltype分布
            for i, (celltype, count) in enumerate(celltype_dist.items(), 1):
                # 使用安全的四舍五入方法
                if category_count > 0:
                    pct_in_category = round(count / category_count * 100, 2)
                else:
                    pct_in_category = 0
                
                pct_total = round(count / total_cells * 100, 2)
                
                if pd.isna(celltype):
                    celltype_str = "NaN/空值"
                else:
                    celltype_str = str(celltype)
                    
                print(f"   {i}. {celltype_str}: {count}个细胞 "
                      f"({pct_in_category}% 的{category}细胞, {pct_total}% 的总细胞)")
        else:
            print("   该annotation类别中没有找到对应的celltype信息")
        
        print()
    
    # 生成汇总表格
    print("\n" + "="*80)
    print("汇总表格:")
    print("="*80)
    
    summary_data = []
    for category in existing_categories:
        if category in results:
            category_info = results[category]
            celltype_dist = category_info['stats']['细胞数量']
            
            for celltype, count in celltype_dist.items():
                if category_info['total_count'] > 0:
                    pct_in_category = round(count / category_info['total_count'] * 100, 2)
                else:
                    pct_in_category = 0
                
                pct_total = round(count / total_cells * 100, 2)
                
                summary_data.append({
                    'Annotation类别': category,
                    'CellType': celltype if not pd.isna(celltype) else 'NaN/空值',
                    '细胞数量': count,
                    f'在{category}中占比(%)': pct_in_category,
                    '在总细胞中占比(%)': pct_total,
                    '该Annotation总细胞数': category_info['total_count']
                })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        # 按Annotation类别和细胞数量排序
        summary_df = summary_df.sort_values(['Annotation类别', '细胞数量'], ascending=[True, False])
        print(summary_df.to_string(index=False, max_cols=6))
    else:
        print("没有找到匹配的数据")
        summary_df = pd.DataFrame()
    
    return results, summary_df

if __name__ == "__main__":
    # 文件路径
    h5ad_file = "adata_celltype_RCTD.h5ad"  # 替换为您的文件路径
    
    # 可选的annotation类别列表（如果与默认不同可以自定义）
    custom_categories = [
        'Unknown', 'Endo-CAF', 'B-CAF', 'Endo-Immu', 'B-T', 
        'T-Mast', 'Plasma-CAF', 'B-Epi', 'T-Epi', 'T-CAF'
    ]
    
    try:
        results, summary_df = analyze_celltype_distribution(h5ad_file, custom_categories)
        
        # 可选：保存详细结果到CSV文件
        if results is not None and not summary_df.empty:
            # 保存汇总结果
            summary_df.to_csv('annotation_celltype_summary.csv', index=False, encoding='utf-8-sig')
            print(f"\n汇总结果已保存至: annotation_celltype_summary.csv")
            
            # 保存每个annotation类别的详细统计
            for category, result_info in results.items():
                result_info['stats'].to_csv(f'{category}_celltype_stats.csv', encoding='utf-8-sig')
                print(f"{category}的详细统计已保存至: {category}_celltype_stats.csv")
                
        print("\n统计完成!")
        
    except Exception as e:
        print(f"处理过程中出现错误: {str(e)}")
        import traceback
        traceback.print_exc()
