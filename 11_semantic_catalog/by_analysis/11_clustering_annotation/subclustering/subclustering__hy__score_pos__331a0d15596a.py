import os
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
from matplotlib.pyplot import rc_context

##读取内置数据集
adata=sc.read_h5ad("integrated_data.h5ad")
#选取这条GOBP_CELL_CYCLE通路上的基因集
#geneset =["AKT1","AKT2","AKT3","BCL2L1","CBL","CBLB","CBLC","CCND1","CCND2","CCND3","CISH","CLCF1","CNTF","CNTFR","CREBBP","CRLF2","CSF2","CSF2RA","CSF2RB","CSF3","CSF3R","CSH1","CTF1","EP300","EPO","EPOR","GH1","GH2","GHR","GRB2","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNE","IFNG","IFNGR1","IFNGR2","IFNK","IFNL1","IFNL2","IFNL3","IFNLR1","IFNW1","IL10","IL10RA","IL10RB","IL11","IL11RA","IL12A","IL12B","IL12RB1","IL12RB2","IL13","IL13RA1","IL13RA2","IL15","IL15RA","IL19","IL2","IL20","IL20RA","IL20RB","IL21","IL21R","IL22","IL22RA1","IL22RA2","IL23A","IL23R","IL24","IL26","IL2RA","IL2RB","IL2RG","IL3","IL3RA","IL4","IL4R","IL5","IL5RA","IL6","IL6R","IL6ST","IL7","IL7R","IL9","IL9R","IRF9","JAK1","JAK2","JAK3","LEP","LEPR","LIF","LIFR","MPL","MYC","OSM","OSMR","PIAS1","PIAS2","PIAS3","PIAS4","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","PIM1","PRL","PRLR","PTPN11","PTPN6","SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS7","SOS1","SOS2","SPRED1","SPRED2","SPRY1","SPRY2","SPRY3","SPRY4","STAM","STAM2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","TPO","TSLP","TYK2"]

geneset = [
    # 受体及配体
    "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", 
    "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNB1", "IFNE", 
    "IFNG", "IFNK", "IFNL1", "IFNL2", "IFNL3", "IFNW1", "IL2", "IL3", "IL4", 
    "IL5", "IL6", "IL7", "IL9", "IL10", "IL11", "IL12A", "IL12B", "IL13", 
    "IL15", "IL19", "IL20", "IL21", "IL22", "IL23A", "IL24", "IL26", "CSF2", 
    "CSF3", "CNTF", "EPO", "GH1", "GH2", "PRL", "TPO", "LEP", "LIF", "OSM",
    
    # JAK激酶家族
    "JAK1", "JAK2", "JAK3", "TYK2",
    
    # STAT转录因子
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    
    # 信号转导增强子
    "GRB2", "SOS1", "SOS2", "AKT1", "AKT2", "AKT3", 
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R5",
    "PIM1", "MYC", "CCND1", "CCND2", "CCND3",
    
    # 受体亚基
    "IL2RA", "IL2RB", "IL2RG", "IL3RA", "IL4R", "IL5RA", "IL6R", "IL6ST", 
    "IL7R", "IL9R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", 
    "IL13RA1", "IL13RA2", "IL15RA", "IL20RA", "IL20RB", "IL21R", 
    "IL22RA1", "IL22RA2", "IL23R", "CSF2RA", "CSF2RB", "CSF3R", "CNTFR", 
    "EPOR", "GHR", "LIFR", "OSMR", "PRLR", "MPL"
]

#进行基因集打分i
sc.tl.score_genes(adata,geneset,score_name='JAK_score')
print(adata.obs)
adata.write('adata_celltype_ad_score_pos.h5ad')
#结果可视化
#sc.pl.umap(adata,color='celltype_final',legend_loc='on data')
#sc.pl.umap(adata,color='JAK_score')
#plt.savefig(f"score.png", dpi=300, bbox_inches='tight')

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# 绘制疾病组的JAK_score
with plt.rc_context():
    sc.pl.umap(
        adata[adata.obs['group'] == 'disease'],
        color='JAK_score',
        #cmap='Reds',  # 使用红色调色板
        vmin=0, vmax=adata.obs['JAK_score'].max(),  # 设置统一的值域
        title='Disease Group - JAK/STAT positively regulated gene set',
        ax=axes[0],
        show=False
    )

# 绘制正常组的JAK_score
with plt.rc_context():
    sc.pl.umap(
        adata[adata.obs['group'] == 'normal'],
        color='JAK_score',
        #cmap='Blues',  # 使用蓝色调色板
        vmin=0, vmax=adata.obs['JAK_score'].max(),  # 设置统一的值域
        title='Normal Group - JAK/STAT positively regulated gene set',
        ax=axes[1],
        show=False
    )

# 调整布局并保存JAK_score图
plt.tight_layout()
plt.savefig("disease_vs_normal_JAK_pos_score.png", bbox_inches='tight')
plt.savefig("disease_vs_normal_JAK_pos_score.pdf", bbox_inches='tight')
plt.close('all')  # 关闭所有图形对象释放内存


library_names = adata.obs['sample'].drop_duplicates().to_list()
out_dir = '/home/hy/Project/FWSP20251095/Redo/subcluster/sub/score'
for library in library_names:
    with rc_context({'figure.figsize': (8, 8)}):
        sq.pl.spatial_scatter(
            adata[adata.obs['sample'] == library,:],
            library_id="spatial",
            color="JAK_score",
            title=f'{library}',
            shape=None,
            size=1,
            img=False,
            use_raw=True
        )
        plt.savefig(os.path.join(out_dir, f'spatial_{library}_pos_score.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, f'spatial_{library}_pos_score.pdf'), facecolor='white', bbox_inches="tight")

