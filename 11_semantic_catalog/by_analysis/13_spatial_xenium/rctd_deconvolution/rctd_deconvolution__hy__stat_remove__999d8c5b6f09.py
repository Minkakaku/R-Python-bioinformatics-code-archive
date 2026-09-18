import scanpy as sc
import pandas as pd

def count_cell_groups(h5ad_file, annotation_categories=None):
    """
    统计指定annotation类别的细胞总数及去除后的细胞总数
    
    参数:
    h5ad_file: h5ad文件路径
    annotation_categories: 需要统计的annotation类别列表
    """
    # 如果未提供annotation类别，使用默认列表
    if annotation_categories is None:
        annotation_categories = [
            'Unknown', 'Endo-CAF', 'B-CAF', 'Endo-Immu', 'B-T', 
            'T-Mast', 'Plasma-CAF', 'B-Epi', 'T-Epi', 'T-CAF'
        ]
    
    # 读取h5ad文件[7](@ref)
    print("读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"数据形状: {adata.shape}")
    
    # 检查annotation列是否存在[7](@ref)
    if 'annotation' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中不存在'annotation'列")
    
    # 获取总细胞数
    total_cells = adata.n_obs
    print(f"\n总细胞数: {total_cells}")
    
    # 检查数据中实际存在的annotation类别
    existing_categories = [cat for cat in annotation_categories if cat in adata.obs['annotation'].values]
    missing_categories = [cat for cat in annotation_categories if cat not in adata.obs['annotation'].values]
    
    if missing_categories:
        print(f"以下annotation类别在数据中不存在: {missing_categories}")
    
    if not existing_categories:
        print("没有找到指定的annotation类别")
        return None
    
    # 筛选指定annotation类别的细胞
    target_mask = adata.obs['annotation'].isin(existing_categories)
    target_cells = adata[target_mask]
    target_count = target_cells.n_obs
    
    # 计算去除这些细胞后的细胞数
    remaining_cells = adata[~target_mask]
    remaining_count = remaining_cells.n_obs
    
    # 打印统计结果
    print("\n" + "="*60)
    print("统计结果:")
    print("="*60)
    
    print(f"指定annotation类别的细胞总数: {target_count}")
    print(f"去除这些细胞后的细胞总数: {remaining_count}")
    print(f"指定类别细胞占比: {target_count/total_cells*100:.2f}%")
    print(f"剩余细胞占比: {remaining_count/total_cells*100:.2f}%")
    
    # 按类别详细统计
    print(f"\n各annotation类别的详细统计:")
    print("-" * 40)
    
    for category in existing_categories:
        category_count = (adata.obs['annotation'] == category).sum()
        category_percent = category_count / total_cells * 100
        print(f"{category}: {category_count}个细胞 ({category_percent:.2f}%)")
    
    # 验证计算是否正确
    print(f"\n验证:")
    print(f"指定类别细胞数 + 剩余细胞数 = {target_count} + {remaining_count} = {target_count + remaining_count}")
    print(f"总细胞数: {total_cells}")
    print(f"计算是否正确: {target_count + remaining_count == total_cells}")
    
    return {
        'total_cells': total_cells,
        'target_cells_count': target_count,
        'remaining_cells_count': remaining_count,
        'target_cells_percent': target_count/total_cells*100,
        'remaining_cells_percent': remaining_count/total_cells*100,
        'category_breakdown': {cat: (adata.obs['annotation'] == cat).sum() for cat in existing_categories}
    }

if __name__ == "__main__":
    # 文件路径
    h5ad_file = "adata_celltype_RCTD.h5ad"  # 替换为您的文件路径
    
    # 可选的annotation类别列表（如果与默认不同可以自定义）
    custom_categories = [
        'Unknown', 'Endo-CAF', 'B-CAF', 'Endo-Immu', 'B-T', 
        'T-Mast', 'Plasma-CAF', 'B-Epi', 'T-Epi', 'T-CAF'
    ]
    
    try:
        result = count_cell_groups(h5ad_file, custom_categories)
        
        if result is not None:
            # 可选：保存统计结果到文本文件
            with open('cell_count_statistics.txt', 'w', encoding='utf-8') as f:
                f.write("细胞统计结果\n")
                f.write("="*50 + "\n")
                f.write(f"总细胞数: {result['total_cells']}\n")
                f.write(f"指定annotation类别的细胞总数: {result['target_cells_count']}\n")
                f.write(f"去除这些细胞后的细胞总数: {result['remaining_cells_count']}\n")
                f.write(f"指定类别细胞占比: {result['target_cells_percent']:.2f}%\n")
                f.write(f"剩余细胞占比: {result['remaining_cells_percent']:.2f}%\n\n")
                
                f.write("各annotation类别的详细统计:\n")
                for category, count in result['category_breakdown'].items():
                    percent = count / result['total_cells'] * 100
                    f.write(f"{category}: {count}个细胞 ({percent:.2f}%)\n")
            
            print(f"\n详细统计结果已保存至: cell_count_statistics.txt")
                
        print("\n统计完成!")
        
    except Exception as e:
        print(f"处理过程中出现错误: {str(e)}")
