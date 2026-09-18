import scanpy as sc
import anndata
import rpy.robjects as ro
from rpy.robjects.packages import importr
from rpy.robjects import pandas2ri, numpy2ri
from rpy.robjects.conversion import localconverter
import os

# 1. 加载h5ad数据并提取必要信息
adata = sc.read_h5ad("/home/hy/Project/FWSP20251095/Redo/subcluster/sub/integrated_data.h5ad")  # 替换为你的h5ad文件路径
if 'celltype_final' not in adata.obs.columns:
    raise ValueError("adata.obs中缺少'celltype_final'列！")

# 2. 准备CellChat所需的输入数据
with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    # 提取表达矩阵（行为基因，列为细胞）
    exp_matrix = adata.to_df().T  # 转置为CellChat要求的格式
    ro.globalenv['data_input'] = exp_matrix
    
    # 提取细胞类型标签
    cell_types = adata.obs['celltype_final'].astype('category')
    ro.globalenv['meta'] = ro.DataFrame({'celltype_final': cell_types.cat.codes})
    ro.globalenv['celltype_names'] = ro.StrVector(cell_types.cat.categories)

# 3. 初始化R环境并加载CellChat
utils = importr('utils')
utils.chooseCRANmirror(ind=1)  # 选择CRAN镜像
if not ro.packages.isinstalled('CellChat'):
    utils.install_packages('CellChat')

cellchat = importr('CellChat')
patchwork = importr('patchwork')

# 4. 创建CellChat对象并设置数据库
ro.r('''
library(CellChat)
library(patchwork)

# 创建CellChat对象
cellchat <- createCellChat(
    object = data_input, 
    meta = meta, 
    group.by = "celltype_final"
)

# 设置细胞类型名称（避免数字标签）
cellchat@idents <- factor(
    cellchat@idents, 
    levels = seq_along(celltype_names) - 1, 
    labels = celltype_names
)

# 使用人类配体-受体数据库（若为小鼠改为CellChatDB.mouse）
cellchat@DB <- CellChatDB.human
''')

# 5. 计算细胞通讯概率
ro.r('''
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 计算通讯概率（使用Truncated Mean减少噪声）
cellchat <- computeCommunProb(
    cellchat, 
    type = "truncatedMean", 
    trim = 0.1  # 去除10%极端值
)

# 过滤低置信通讯（每组至少10个细胞）
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 聚合通路水平通讯
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
''')

# 6. 可视化关键结果
output_dir = "cellchat_results"
os.makedirs(output_dir, exist_ok=True)

# (1) 细胞互作热图
ro.r(f'''
pdf("{output_dir}/interaction_heatmap.pdf", width=10, height=8)
netVisual_heatmap(cellchat, measure = "count", font.size = 8)
dev.off()
''')

# (2) 通路水平气泡图
ro.r(f'''
pdf("{output_dir}/pathway_bubble.pdf", width=12, height=10)
netAnalysis_signalingRole_scatter(cellchat, font.size = 6)
dev.off()
''')

# (3) 所有细胞类型互作网络图
ro.r(f'''
pdf("{output_dir}/full_network.pdf", width=15, height=15)
netVisual_aggregate(
    cellchat, 
    layout = "circle", 
    edge.width.max = 8, 
    vertex.label.cex = 1.2
)
dev.off()
''')

# 7. 保存结果
ro.r(f'saveRDS(cellchat, file="{output_dir}/cellchat_object.rds")')
ro.r(f'write.csv(subsetCommunication(cellchat), file="{output_dir}/cell_communication_table.csv")')

print(f"分析完成！结果保存至: {output_dir}")
