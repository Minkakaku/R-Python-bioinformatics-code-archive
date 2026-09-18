import scanpy as sc
import pandas as pd

# 1. 读取h5ad文件
print("正在读取h5ad文件...")
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10.h5ad")

# 2. 检查必需列
if 'sample_new' not in adata.obs.columns or 'state' not in adata.obs.columns:
    print(f"错误: 缺少必要列")
    print(f"可用列: {adata.obs.columns.tolist()}")
else:
    # 3. 提取样本-状态对应关系（去重）
    sample_state_df = adata.obs[['sample_new', 'state']].drop_duplicates().sort_values('sample_new')
    
    # 4. 保存为CSV
    sample_state_df.to_csv("sample_state_list.csv", index=False, encoding='utf-8-sig')
    print(f"已保存到: sample_state_list.csv")
    
    # 5. 打印列表
    print("\n样本-状态对应列表:")
    print("="*40)
    for i, row in sample_state_df.iterrows():
        print(f"{row['sample_new']}: {row['state']}")
