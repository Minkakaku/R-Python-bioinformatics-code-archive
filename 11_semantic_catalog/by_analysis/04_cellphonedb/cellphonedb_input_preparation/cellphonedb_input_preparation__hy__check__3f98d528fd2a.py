import pandas as pd
meta_cells = set(pd.read_csv("meta.csv")["cell_id"])
count_cells = set(pd.read_csv("count.csv", nrows=1).columns[1:])  # 跳过首列"Gene_id"
# 检查差异细胞
missing_in_meta = count_cells - meta_cells
missing_in_count = meta_cells - count_cells
print(f"Count中独有的细胞: {missing_in_meta}")
print(f"Meta中独有的细胞: {missing_in_count}")
