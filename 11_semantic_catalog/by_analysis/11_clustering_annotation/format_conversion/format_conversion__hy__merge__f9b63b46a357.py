import scanpy as sc
import anndata as ad
import os

def merge_h5ad_files(file_list, output_file="merged_data.h5ad", batch_key="batch"):
    """
    合并多个h5ad文件到一个AnnData对象中。

    参数
    ----------
    file_list : list of str
        需要合并的h5ad文件路径列表。
    output_file : str, optional
        合并后输出的文件名，默认为"merged_data.h5ad"。
    batch_key : str, optional
        在obs中添加的批次信息列名，默认为"batch"。
    """

    adata_list = []  # 用于存储读取的AnnData对象
    batch_names = [] # 用于存储批次名称（从文件名提取）

    print("开始读取h5ad文件...")
    
    # 1. 逐个读取h5ad文件
    for file_path in file_list:
        if not os.path.exists(file_path):
            print(f"警告: 文件 {file_path} 不存在，已跳过。")
            continue
        
        print(f"正在读取: {file_path}")
        try:
            # 读取单个h5ad文件
            adata = sc.read_h5ad(file_path)
            
            # 确保观测名（obs index，通常是细胞Barcode）唯一
            # 这是避免合并时出现"ValueError: cannot reindex from a duplicate axis"错误的关键步骤[6](@ref)
            adata.obs_names_make_unique()
            
            # 从文件名提取批次名（去掉路径和扩展名）
            batch_name = os.path.splitext(os.path.basename(file_path))[0]
            batch_names.append(batch_name)
            
            # 为当前数据集的观测数据（obs）添加批次来源信息[8](@ref)
            adata.obs[batch_key] = batch_name
            
            adata_list.append(adata)
            print(f"  成功读取，数据形状: {adata.shape} (细胞数 x 基因数)")
            
        except Exception as e:
            print(f"  读取文件 {file_path} 时出错: {e}")
            continue

    if len(adata_list) == 0:
        print("错误: 没有成功读取任何h5ad文件。")
        return None

    print(f"\n成功读取 {len(adata_list)} 个文件，开始合并...")

    # 2. 合并所有AnnData对象[7,8](@ref)
    # 使用join='inner'可以只保留所有数据集共有的基因（var），这是较安全的做法[8](@ref)
    # 如果希望保留所有基因（即使某些基因只存在于部分数据集），可使用join='outer'
    try:
        # 方法1: 使用anndata的concat函数[7](@ref)
        merged_adata = ad.concat(adata_list, 
                                  join='inner', 
                                  label=batch_key, 
                                  keys=batch_names,
                                  index_unique='-')  # 使合并后的细胞索引唯一，例如在原Barcode后加上批次名[6](@ref)
        
        # 方法2: 也可以使用scanpy的AnnData.concatenate，功能类似[8](@ref)
        # merged_adata = adata_list[0].concatenate(adata_list[1:], 
        #                                         batch_key=batch_key, 
        #                                         batch_categories=batch_names)
        
        print("合并成功！")
        print(f"合并后数据形状: {merged_adata.shape}")
        print(f"批次信息: {merged_adata.obs[batch_key].value_counts()}")
        
    except Exception as e:
        print(f"合并过程中出错: {e}")
        # 尝试使用更宽松的合并方式（保留所有基因）
        try:
            print("尝试使用outer方式合并（保留所有基因）...")
            merged_adata = ad.concat(adata_list, 
                                      join='outer', 
                                      label=batch_key, 
                                      keys=batch_names,
                                      index_unique='-')
            print("outer方式合并成功！")
            print(f"合并后数据形状: {merged_adata.shape}")
        except Exception as e2:
            print(f"outer方式合并也失败: {e2}")
            return None

    # 3. 保存合并后的数据
    try:
        merged_adata.write_h5ad(output_file)
        print(f"\n合并后的数据已保存至: {output_file}")
    except Exception as e:
        print(f"保存文件时出错: {e}")

    return merged_adata

if __name__ == "__main__":
    # 需要合并的文件列表（请根据实际情况修改路径）
    file_list = [
        "adata_celltype.h5ad",
        "celltype_B_cells.h5ad", 
        "celltype_Endothelials.h5ad",
        "celltype_Fibroblasts.h5ad",
        "celltype_Macrophages.h5ad",
        "celltype_Myeloids.h5ad",
        "celltype_Plasma.h5ad", 
        "celltype_Smooth_muscle_cells.h5ad",
        "celltype_T_cells.h5ad"
    ]
    
    # 检查文件是否存在
    existing_files = [f for f in file_list if os.path.exists(f)]
    missing_files = set(file_list) - set(existing_files)
    
    if missing_files:
        print("以下文件不存在，将不会被合并:")
        for f in missing_files:
            print(f"  {f}")
    
    if not existing_files:
        print("没有找到可合并的文件，请检查文件路径。")
    else:
        # 执行合并
        merged_data = merge_h5ad_files(existing_files, 
                                       output_file="merged_celltypes.h5ad",
                                       batch_key="source")
        
        if merged_data is not None:
            print("\n合并完成！可以开始后续分析。")
