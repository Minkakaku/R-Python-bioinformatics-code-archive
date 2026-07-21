import scanpy as sc
import pandas as pd

# 读取h5ad数据
adata = sc.read_h5ad("result.h5ad")

# 提取sample和group列
if 'sample' in adata.obs.columns and 'group' in adata.obs.columns:
    result = adata.obs[['sample', 'group']].copy()
    
    # 保存到CSV
    result.to_csv("sample_group_info.csv", index=False)
    print("信息已保存到: sample_group_info.csv")
    
    # 显示基本信息
    print(f"样本数: {result['sample'].nunique()}")
    print(f"分组数: {result['group'].nunique()}")
    print(f"前10行数据:")
    print(result.head(10))
else:
    print("错误: 数据中缺少'sample'或'group'列")
