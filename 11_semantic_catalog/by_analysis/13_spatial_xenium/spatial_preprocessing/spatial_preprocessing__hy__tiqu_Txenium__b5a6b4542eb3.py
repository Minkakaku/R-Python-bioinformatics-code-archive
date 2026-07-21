import pandas as pd
import scanpy as sc

# 1. 读取 T 细胞类型列表（无表头，每行一个类型名）
t_types = pd.read_csv("bb", header=None)[0].str.strip().tolist()

# 2. 读取 AnnData
adata = sc.read_h5ad("output-XETG00305__0065812_CT.h5ad")

# 3. 确认必要列存在
assert "celltype" in adata.obs.columns, "obs 中缺少 'celltype' 列"
assert "sample" in adata.obs.columns, "obs 中缺少 'sample' 列，请根据实际数据修改列名"

# 4. 筛选出 T 细胞
t_mask = adata.obs["celltype"].isin(t_types)
adata_t = adata[t_mask].copy()

# 5. 提取 index（细胞观测名）与 sample
result = pd.DataFrame({
    "sample": adata_t.obs["cell_id_sample"].values
})

# 6. 保存为 CSV，仅包含 index 和 sample 两列
result.to_csv("xenium_T.csv", index=False)
