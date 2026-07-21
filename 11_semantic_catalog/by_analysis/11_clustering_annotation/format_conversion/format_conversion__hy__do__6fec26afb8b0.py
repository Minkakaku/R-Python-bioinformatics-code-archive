import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

# 设置随机种子确保可重复性
np.random.seed(42)

# 读取数据
adata = sc.read('adata_cluster.h5ad')

# 1. 获取分类比例
class_counts = adata.obs['res_0.5'].value_counts(normalize=True)

# 2. 计算每类应抽取的细胞数（按比例分配）
total_cells = 50000
class_samples = (class_counts * total_cells).round().astype(int)

# 3. 调整总数确保精确为50000（由于四舍五入可能导致总数偏差）
diff = total_cells - class_samples.sum()
if diff != 0:
    # 将差额加到最大的类别中
    largest_class = class_counts.idxmax()
    class_samples[largest_class] += diff

# 4. 分层抽样
sampled_indices = []  # 存储抽样细胞的索引
sampled_cell_ids = []  # 存储抽样细胞的原始ID

for class_label, n_samples in class_samples.items():
    # 获取当前类别的所有细胞索引和ID
    class_mask = adata.obs['res_0.5'] == class_label
    class_indices = adata.obs.index[class_mask].tolist()
    
    # 随机抽取
    if len(class_indices) <= n_samples:
        # 如果类别细胞数少于需要样本数，取全部
        sampled_indices.extend(class_indices)
        sampled_cell_ids.extend(class_indices)  # 保留原始ID
    else:
        selected = np.random.choice(class_indices, n_samples, replace=False)
        sampled_indices.extend(selected)
        sampled_cell_ids.extend(selected)  # 保留原始ID

# 5. 创建子集
adata_sampled = adata[sampled_indices].copy()

# 6. 保存抽样细胞ID到CSV文件
pd.DataFrame({'cell_id': sampled_cell_ids}).to_csv('sampled_cell_ids.csv', index=False)

# 7. 保存抽样后的数据
adata_sampled.write('adata_sampled_cluster.h5ad')

# 8. 验证抽样结果
print("原始数据细胞数:", adata.n_obs)
print("抽样后细胞数:", adata_sampled.n_obs)
print("\n原始分类比例:")
print(adata.obs['res_0.5'].value_counts(normalize=True))
print("\n抽样后分类比例:")
print(adata_sampled.obs['res_0.5'].value_counts(normalize=True))
print("\n抽样细胞ID已保存至: sampled_cell_ids.csv")
