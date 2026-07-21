# subcluster_mapping.py
import scanpy as sc
import pandas as pd
import numpy as np
import os
import logging
from typing import Union, List

def map_subcluster_annotations(
    main_adata_path: str,
    sub_adata_paths: Union[str, List[str]],
    output_path: str = None,
    main_annotation_col: str = "cell_type",
    sub_annotation_col: str = "cell_type",
    verbose: bool = True
) -> pd.DataFrame:
    """
    将亚群注释映射回原始大群数据集
    
    参数:
        main_adata_path: 主群AnnData文件路径(.h5ad)
        sub_adata_paths: 亚群AnnData文件路径或路径列表
        output_path: 输出CSV文件路径(可选)
        main_annotation_col: 主群中细胞注释列名(默认'cell_type')
        sub_annotation_col: 亚群中细胞注释列名(默认'cell_type')
        verbose: 是否打印进度信息
    
    返回:
        pd.DataFrame: 更新后的注释DataFrame
    """
    # 设置日志
    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING)
    logger = logging.getLogger(__name__)
    
    # 确保sub_adata_paths是列表
    if isinstance(sub_adata_paths, str):
        sub_adata_paths = [sub_adata_paths]
    
    logger.info(f"加载主群数据: {main_adata_path}")
    main_adata = sc.read_h5ad(main_adata_path)
    
    # 备份原始注释
    if f"{main_annotation_col}_old" in main_adata.obs.columns:
        logger.warning(f"列 '{main_annotation_col}_old' 已存在，将被覆盖")
    
    main_adata.obs[f"{main_annotation_col}_old"] = main_adata.obs[main_annotation_col].copy()
    
    # 初始化新注释列
    main_adata.obs[main_annotation_col] = main_adata.obs[f"{main_annotation_col}_old"].copy()
    
    # 存储映射统计信息
    mapping_stats = []
    total_cells = 0
    
    # 处理每个亚群
    for i, path in enumerate(sub_adata_paths):
        logger.info(f"处理亚群 {i+1}/{len(sub_adata_paths)}: {path}")
        
        if not os.path.exists(path):
            logger.error(f"文件不存在: {path}")
            continue
            
        try:
            sub_adata = sc.read_h5ad(path)
            
            # 检查注释列是否存在
            if sub_annotation_col not in sub_adata.obs.columns:
                logger.error(f"亚群数据中缺少注释列 '{sub_annotation_col}': {path}")
                continue
                
            # 检查细胞ID是否匹配
            common_cells = set(sub_adata.obs_names) & set(main_adata.obs_names)
            if len(common_cells) == 0:
                logger.error(f"在亚群 {path} 和主群之间找不到共同的细胞ID")
                continue
                
            # 创建映射字典
            mapping_dict = sub_adata.obs[sub_annotation_col].to_dict()

            # 找出新增注释类别
            new_categories = pd.Series(list(mapping_dict.values())).unique()

            # 如果是 Categorical 类型，则动态添加新类别
            if not pd.api.types.is_categorical_dtype(main_adata.obs[main_annotation_col]):
                main_adata.obs[main_annotation_col] = main_adata.obs[main_annotation_col].astype("category")

            existing_categories = main_adata.obs[main_annotation_col].cat.categories
            categories_to_add = [cat for cat in new_categories if cat not in existing_categories]
            if categories_to_add:
                main_adata.obs[main_annotation_col] = main_adata.obs[main_annotation_col].cat.add_categories(categories_to_add)


            # 正式更新注释
            for cell_id, new_annot in mapping_dict.items():
                if cell_id in main_adata.obs_names:
                    main_adata.obs.at[cell_id, main_annotation_col] = new_annot
                        
            # 记录统计信息
            num_mapped = len(common_cells)
            total_cells += num_mapped
            mapping_stats.append({
                "subcluster": os.path.basename(path),
                "cells_mapped": num_mapped,
                "original_annot": main_adata.obs.loc[list(common_cells), f"{main_annotation_col}_old"].value_counts().to_dict(),
                "new_annot": main_adata.obs.loc[list(common_cells), main_annotation_col].value_counts().to_dict()
            })
            
        except Exception as e:
            logger.error(f"处理亚群 {path} 时出错: {str(e)}")
    
    # 打印统计信息
    logger.info(f"总计映射细胞: {total_cells}/{len(main_adata.obs)}")
    for col in main_adata.obs.columns:
        if pd.api.types.is_categorical_dtype(main_adata.obs[col]):
            main_adata.obs[col] = main_adata.obs[col].cat.remove_unused_categories()
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        main_adata.obs.to_csv(output_path)
        logger.info(f"注释结果已保存至: {output_path}")
    
    # 返回映射统计
    return main_adata.obs, mapping_stats

if __name__ == "__main__":
    main_path = "/path/to/main_adata.h5ad"
    sub_paths = ["/path/to/sub1.h5ad", "/path/to/sub2.h5ad"]
    output_csv = "/path/to/updated_annotations.csv"

    updated_obs, stats = map_subcluster_annotations(
        main_adata_path=main_path,
        sub_adata_paths=sub_paths,
        output_path=output_csv
    )

    print("✅ 前5行更新后的注释：")
    print(updated_obs[[f"{main_annotation_col}_old", main_annotation_col]].head())

    print("\n📊 映射统计：")
    for stat in stats:
        print(f"\n亚群: {stat['subcluster']}")
        print(f"映射细胞数: {stat['cells_mapped']}")
        print(f"原始注释分布: {stat['original_annot']}")
        print(f"新注释分布: {stat['new_annot']}")