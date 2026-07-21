import pandas as pd
import anndata as ad
import geopandas as gpd
from shapely.geometry import Point
from spatialdata_io import xenium

# 1. 读数据
sdata = xenium(
    path="xenium_data_dir",
    nucleus_boundaries=True,
    transcripts=True,
    cells_boundaries=False,   # 我们只关心核，不读细胞轮廓
)

# 2. 转录本 -> GeoDataFrame
tx = sdata.points["transcripts"].compute()
tx["geometry"] = gpd.points_from_xy(tx["x"], tx["y"])
tx_gdf = gpd.GeoDataFrame(tx, geometry="geometry")

# 3. 核边界 -> GeoDataFrame（含 cell_id）
nuc = sdata.shapes["nucleus_boundaries"].reset_index()  # 通常叫 'index' 或 'cell_id'

# 4. 空间连接：转录本落在哪个核
tx_in_nuc = gpd.sjoin(tx_gdf, nuc, predicate="within")  # 得到列 'index'（即 cell_id）

# 5. 按 cell_id + gene 聚合计数
exp = (
    tx_in_nuc
    .groupby(["index", "feature_name"], observed=True)
    .size()
    .unstack(fill_value=0)
    .astype(int)
)

# 6. 让行名、列名符合 AnnData 规范
exp.index.name = "cell_id"
exp.columns.name = "gene"

# 7. 创建新的 AnnData
adata_nuc = ad.AnnData(
    X=exp.values,
    obs=pd.DataFrame(index=exp.index),
    var=pd.DataFrame(index=exp.columns)
)