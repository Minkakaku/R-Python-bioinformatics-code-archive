import scanpy as sc
import pandas as pd

# 读取数据
adata = sc.read_h5ad('adata_celltype.h5ad')
print(f"原始数据: {adata.shape}")

# 检查cell_id_sample列是否存在
if 'cell_id_sample' in adata.obs.columns:
    # 将cell_id_sample列设置为索引
    adata.obs.index = adata.obs['cell_id_sample'].astype(str)
    adata.obs = adata.obs.drop(columns=['cell_id_sample'])
    
    # 去重
    if not adata.obs.index.is_unique:
        print("发现重复索引，开始去重...")
        adata = adata[~adata.obs.index.duplicated(keep='first'), :].copy()
        print(f"去重后: {adata.shape}")
    else:
        print("索引唯一，无需去重")
    
    # 保存结果
    adata.write('data_sp_ep.h5ad')
    print("已保存为 data_sp_ep.h5ad")
else:
    print("错误: 未找到'cell_id_sample'列")
    print(f"可用列: {list(adata.obs.columns)}")
