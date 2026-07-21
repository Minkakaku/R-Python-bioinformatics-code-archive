import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# 1. 数据准备
adata = sc.read_h5ad("/home/hy/Project/FWSP20251095/Redo/subcluster/sub/integrated_data.h5ad")  # 替换为实际路径
adata.var_names_make_unique()

# 创建输出目录
os.makedirs("cellchat_results", exist_ok=True)
os.makedirs("cellphonedb_results", exist_ok=True)

# 准备CellPhoneDB输入文件
meta_data = adata.obs[['celltype_final']].reset_index()
meta_data.columns = ['Cell', 'cell_type']
meta_data.to_csv("cellphonedb_meta.txt", sep='\t', index=False)

counts_data = pd.DataFrame(adata.X.toarray(), 
                          index=adata.obs.index,
                          columns=adata.var.index)
counts_data.T.to_csv("cellphonedb_counts.txt", sep='\t')

# 2. CellPhoneDB分析
print("运行CellPhoneDB分析...")
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path='cellphonedb.zip',  # 需提前下载数据库
    meta_file_path='cellphonedb_meta.txt',
    counts_file_path='cellphonedb_counts.txt',
    counts_data='hgnc_symbol',
    output_path='cellphonedb_results',
    iterations=1000,
    threshold=0.1,
    threads=4
)

# 3. CellChat分析(R环境调用)
print("运行CellChat分析...")
pandas2ri.activate()
numpy2ri.activate()

# 初始化R环境
cellchat = importr('CellChat')
base = importr('base')
utils = importr('utils')
patchwork = importr('patchwork')

# 将数据转换为R对象
with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    r_data = ro.conversion.py2rpy(counts_data.T)
    r_meta = ro.conversion.py2rpy(meta_data.set_index('Cell'))

# 创建CellChat对象
ro.globalenv['data.input'] = r_data
ro.globalenv['meta'] = r_meta

ro.r('''
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
    cellchat <- setIdent(cellchat, ident.use = "cell_type")
    CellChatDB <- CellChatDB.human  # 根据物种选择
    cellchat@DB <- CellChatDB
    cellchat <- subsetData(cellchat)
    future::plan("multiprocess", workers = 4)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # 保存结果
    saveRDS(cellchat, "cellchat_results/cellchat_object.rds")
    
    # 基础可视化
    pdf("cellchat_results/cellchat_net_visualization.pdf")
    netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)))
    dev.off()
''')

# 4. 可视化函数
def visualize_cpdb_results():
    # CellPhoneDB点图
    means = pd.read_csv('cellphonedb_results/means.txt', sep='\t')
    pvals = pd.read_csv('cellphonedb_results/pvalues.txt', sep='\t')
    
    # 选择显著交互作用
    significant_means = means[pvals.iloc[:, -1] < 0.05]
    
    plt.figure(figsize=(12, 8))
    for i, row in significant_means.iterrows():
        plt.scatter(row['interacting_pair'], row.index,
                    s=-np.log10(pvals.iloc[i, -1]+1e-300)*10,
                    c=row.iloc[-1], cmap='viridis')
    plt.colorbar(label='Interaction Strength')
    plt.xlabel('Ligand-Receptor Pairs')
    plt.ylabel('Cell Type Pairs')
    plt.title('CellPhoneDB: Significant Interactions')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('cellphonedb_results/cpdb_dotplot.pdf')
    plt.close()

def visualize_cellchat_results():
    # CellChat网络图（使用R）
    ro.r('''
        cellchat <- readRDS("cellchat_results/cellchat_object.rds")
        pdf("cellchat_results/cellchat_interaction_heatmap.pdf")
        netVisual_heatmap(cellchat)
        dev.off()
        
        pdf("cellchat_results/cellchat_pathway_signaling.pdf")
        for(pathway in cellchat@netP$pathways){
            netVisual_aggregate(cellchat, signaling = pathway)
        }
        dev.off()
    ''')

# 5. 执行可视化
print("生成可视化结果...")
visualize_cpdb_results()
visualize_cellchat_results()

print("分析完成！结果保存在 cellchat_results/ 和 cellphonedb_results/ 目录")
