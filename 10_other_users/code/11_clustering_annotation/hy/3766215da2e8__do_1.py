import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 准备数据 (这里创建DataFrame，实际使用时可以直接读取您的数据源)
data = {
    "Symbol": ["SOCS3", "SOCS3", "SOCS3", "SOCS3", "SOCS1", "SOCS3", "SOCS3", "SOCS3", "SOCS3",
               "STAT2", "STAT1", "JAK3", "STAT6", "TYK2", "JAK3", "IFNGR1", "IFNGR1"],
    "logfoldchanges": [-2.5291514, -2.0900233, -3.8919094, -2.829026, -0.835223, -5.2207017,
                       -2.5210917, -3.9351864, -4.1098604, 1.5638305, 0.97596735, 1.3301098,
                       0.5046044, 0.89681333, 1.0853907, 0.5783798, 0.86135674],
    "pct_nz_group": [0.20754290876242096, 0.31732840874115115, 0.48451829655861256,
                     0.26964636542239684, 0.13015717092337917, 0.4542728635682159,
                     0.3262987012987013, 0.42178217821782177, 0.3280828082808281,
                     0.16230489379180849, 0.1491608966083793, 0.2642073778664008,
                     0.235543369890329, 0.22504055449048813, 0.13641055891461437,
                     0.22872732635304527, 0.1823297692740574],
    "group": ["Kerationocytes", "Stratum basal cells", "Fibroblasts", "T cells", "T cells",
              "Mast cells", "Myeloid", "Macrophages", "Epithelials", "Fibroblasts",
              "Fibroblasts", "T cells", "T cells", "Myeloid", "Myeloid", "Myeloid", "Macrophages"]
}
df = pd.DataFrame(data)

# 设置图形和坐标轴
plt.figure(figsize=(10, 8))
ax = plt.subplot(1, 1, 1)

# 创建自定义的蓝色渐变颜色映射 (负值越深蓝色越深)
#colors = ["#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"]  # 浅蓝到深蓝
colors_neg = ["#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]
colors_pos = ["#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"]
colors = colors_neg + colors_pos
cmap = LinearSegmentedColormap.from_list("blue_cmap", colors)

# 获取logfoldchanges的范围用于颜色映射
vmin, vmax = df["logfoldchanges"].min(), df["logfoldchanges"].max()
norm = plt.Normalize(vmin, vmax)

# 获取唯一组和基因用于坐标轴
groups = df["group"].unique()
symbols = df["Symbol"].unique()

# 设置点的大小范围 (基于pct_nz_group值)
min_size = min(df["pct_nz_group"]) * 1000
max_size = max(df["pct_nz_group"]) * 1000

# 绘制点图
for group in groups:
    group_idx = list(groups).index(group)  # y轴位置
    for symbol in symbols:
        # 检查数据中是否有该组和基因的数据点
        condition = (df["group"] == group) & (df["Symbol"] == symbol)
        if not df[condition].empty:
            data_point = df[condition].iloc[0]
            x = np.where(symbols == symbol)[0][0]  # x轴位置
            color = cmap(norm(data_point["logfoldchanges"]))
            size = data_point["pct_nz_group"] * 1000  # 缩放点的大小
            plt.scatter(x, group_idx, s=size, c=[color], edgecolor="black", alpha=0.8)

# 添加轴标签
plt.xticks(range(len(symbols)), symbols, rotation=45, ha="right")
plt.yticks(range(len(groups)), groups)
plt.xlabel("Genes")
plt.ylabel("Cell Groups")

# 添加颜色条
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("logfoldchanges")

# 添加点的大小图例
for size in [0.1, 0.2, 0.4]:
    plt.scatter([], [], s=size * 1000, c="gray", edgecolor="black", label=f"{size:.1f}", alpha=0.6)
plt.legend(title="pct_nz_group", labelspacing=1.2, borderpad=1.2, frameon=True, loc="upper right")

# 设置布局并保存
plt.tight_layout()
plt.savefig("dotplot1.png", dpi=300, bbox_inches="tight")
plt.show()
