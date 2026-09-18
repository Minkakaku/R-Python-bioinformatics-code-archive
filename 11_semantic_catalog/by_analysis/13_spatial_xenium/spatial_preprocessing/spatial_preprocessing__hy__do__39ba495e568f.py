import scanpy as sc
import pandas as pd

# ========== 患者列表 ==========
mpc_patients = [
    "Z2038176-15","Q2324341-1","Z1530765-11","F1706967-11","F1805049-2",
    "Z1915172-7","Z1954693-9","ZF2207170-8","ZF2202916-8","ZF2208790-1",
    "Z2262137-5","Z2271577-12","ZF2304103-8","Z2323818-8","Z2336500-6",
    "Z2345459-7","Z2332034-5","Z2417954-8","Q2236500-1","Q2406275-2",
    "Z2101519-1","Z2123354-1","Z1642635-14","Z2224197-10","Z1736838-15",
    "Z1944219-14","Z1430967-9","Z2057500-35","Z2171766-3","Z2353385-10",
    "Z1635989-12","Z1814675-2","Z1826868-6","Z1952876-9","Z1802235-10"
]

tnbc_patients = [
    "Z2279427-2","H2510098-1","H2510098-2","Q2243285-1","Z2370717-2",
    "Q2337205-2","Q2233699-2","Q2239242-1","Q2319383-1","Q2322393-1",
    "Q2337350-1"
]

# ========== 读取数据 ==========
adata = sc.read("adata_with_band_by_sample_improved.h5ad")
print(f"原始细胞总数: {adata.n_obs}")

# 统一 patient 为字符串
adata.obs['patient'] = adata.obs['patient'].astype(str)

# 构建 cancer_type 映射字典
patient_to_cancer = {}
for p in mpc_patients:
    patient_to_cancer[p] = "MPC"
for p in tnbc_patients:
    patient_to_cancer[p] = "TNBC"

# 仅保留在列表中的患者
mask = adata.obs['patient'].isin(patient_to_cancer.keys())
adata = adata[mask].copy()
print(f"筛选后细胞数: {adata.n_obs}")

# 添加 cancer_type 列
adata.obs['cancer_type'] = adata.obs['patient'].map(patient_to_cancer)

# 添加 cancer_type_zone 列（cancer_type_mucus_zone）
if 'mucus_zone' in adata.obs.columns:
    adata.obs['cancer_type_zone'] = (
        adata.obs['cancer_type'].astype(str) + '_' + adata.obs['mucus_zone'].astype(str)
    )
    print("已创建 'cancer_type_zone' 列")
else:
    print("警告: 'mucus_zone' 列不存在，无法创建 'cancer_type_zone'")

# 打印分布信息
print("\ncancer_type 分布:")
print(adata.obs['cancer_type'].value_counts())
if 'cancer_type_zone' in adata.obs.columns:
    print("\ncancer_type_zone 分布:")
    print(adata.obs['cancer_type_zone'].value_counts())

# 保存结果
adata.write("final_all.h5ad")
print("\n已保存至: final_all.h5ad")
