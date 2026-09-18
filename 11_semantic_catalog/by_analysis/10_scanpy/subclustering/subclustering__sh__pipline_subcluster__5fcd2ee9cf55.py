import os, sys
import numpy as np
import scanpy as sc
import matplotlib
# from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
# from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from matplotlib import font_manager as fm

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    '''
    out_dir : 输出文件夹地址
    file_name : 文件名称
    fig_size : 输出图片的大小
    '''
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/Main_mapping/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512027077-01/06.Mission/Ccl3_Ccr5_gene_plot/"
sc.pl.umap(adata, color=["Ccl3","Ccr5"], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=["Ccl3","Ccr5"], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Ccl3","Ccr5"], groupby='sub_cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sub_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sub_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Ccl3","Ccr5"], groupby='group', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Ccl3","Ccr5"], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Ccl3","Ccr5"], groupby='sample', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Ccl3","Ccr5"], groupby=['cell_type','sample'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-sample.pdf'),bbox_inches='tight')

genes=['Ccl1','Ccl25','Ccl19','Ccl21b','Ccl21e','Ccl21f','Ccl21d','Ccl21a','Ccl3','Ccl4','Ccl17','Ccl22','Ccl5','Ccl8','Ccl9','Ccl6','Ccl7','Ccl12','Ccl2','Ccl11','Ccl24','Ccl26','Ccl27al','Ccl27b','Gm20878','Gm21586','Ccl27a','Gm21953','Ccl28','Ccl20','Cxcl1','Cxcl2','Cxcl3','Cxcl5','Ppbp','Pf4','Cxcl9','Cxcl10','Cxcl11','Cxcl13','Cxcl12','Cxcl16','Cxcl14','Cxcl15','Cxcl17','Xcl1','Cx3cl1','Il2','Il4','Il7','Il9','Il15','Il21','Tslp','Il3','Il5','Csf2','Il13','Epo','Gh','Prl','Prl6a1','Prl5a1','Thpo','Csf3','Lep','Il6','Il11','Il12a','Il12b','Il23a','Il27','Ebi3','Il31','Clcf1','Cntf','Ctf1','Ctf2','Lif','Osm','Il10','Il19','Il20','Il24','Il22b','Il22','Ifnl2','Ifnl3','Ifna1','Ifna11','Ifna2','Ifna4','Ifna5','Ifna6','Ifna7','Ifna9','Ifnab','Ifna13','Ifna16','Ifna15','Ifna12','Ifna14','Ifnb1','Ifnk','Ifne','Ifng','Il1a','Il1b','Il1rn','Il36rn','Il36a','Il36b','Il36g','Il1f10','Il18','Il33','Il17a','Il17f','Il17b','Il17c','Il17d','Il25','Il16','Il34','Csf1','Tnf','Lta','Ltb','Tnfsf14','Fasl','Tnfsf15','Tnfsf10','Eda','Klk1b4','Ngf','Tnfsf11','Tnfsf12','Cd70','Tnfsf8','Cd40lg','Tnfsf9','Tnfsf4','Tnfsf18','Tnfsf13','Tnfsf13b','Tgfb1','Tgfb2','Tgfb3','Gdf15','Gdf2','Bmp10','Inha','Bmp3','Gdf10','Gdf11','Mstn','Inhba','Inhbb','Gdf1','Gdf3','Nodal','Gdf9','Inhbc','Inhbe','Amh','Bmp2','Bmp4','Gdf5','Gdf6','Gdf7','Bmp15','Bmp5','Bmp6','Bmp7','Bmp8a','Bmp8b','Ccr8','Ccr9','Ackr4','Ccr7','Ccr4','Ccr5','Ccr3','Ccr2','Ccr1','Ccr1l1','Ccr10','Ccr6','Cxcr1','Cxcr2','Cxcr3','Cxcr5','Cxcr4','Ackr3','Cxcr6','Xcr1','Cx3cr1','Il2ra','Il2rb','Il2rg','Il4ra','Il7r','Il9r','Il15ra','Il21r','Crlf2','Il3ra','Csf2rb','Csf2rb2','Il5ra','Csf2ra','Il13ra1','Il13ra2','Epor','Ghr','Prlr','Mpl','Csf3r','Lepr','Il6ra','Il6st','Il11ra3','Il11ra1','Il11ra2','Il12rb1','Il12rb2','Il23r','Il27ra','Il31ra','Cntfr','Lifr','Osmr','Il10ra','Il10rb','Il20ra','Il20rb','Il22ra1','Ifnlr1','Ifnar1','Ifnar2','Ifngr1','Ifngr2','Il1r1','Il1rap','Il1r2','Il1rl2','Il18r1','Il18rap','Il1rl1','Il17ra','Il17rc','Il17rb','Il17re','Cd4','Csf1r','Tnfrsf1a','Tnfrsf1b','Ltbr','Tnfrsf14','Fas','Tnfrsf25','Tnfrsf10b','Tnfrsf21','Edar','Eda2r','Ngfr','Tnfrsf11b','Tnfrsf11a','Tnfrsf12a','Cd27','Tnfrsf8','Cd40','Tnfrsf9','Tnfrsf4','Tnfrsf18','Tnfrsf17','Tnfrsf13b','Tnfrsf13c','Tnfrsf19','Relt','Rell1','Rell2','Tgfbr1','Tgfbr2','Acvrl1','Acvr2a','Bmpr2','Acvr2b','Acvr1b','Acvr1c','Acvr1','Amhr2','Bmpr1a','Bmpr1b']
sc.tl.score_genes(adata, gene_list=genes, score_name='Cytokine-cytokine receptor interaction')
sc.pl.umap(adata, color='Cytokine-cytokine receptor interaction', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_Cytokine-cytokine receptor interaction', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='Cytokine-cytokine receptor interaction', groupby='group', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_Cytokine-cytokine receptor interaction', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['Cytokine-cytokine receptor interaction'], groupby=['sub_cell_type'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_Cytokine-cytokine receptor interaction', fig_size={'w': 6, 'h': 12})


out_dir="/mnt/NexData_file/home/sh/Projects/AE2512027077-01/06.Mission/App_Cd74_gene_plot/"
sc.pl.umap(adata, color=["App","Cd74"], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='sub_cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sub_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sub_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='group', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='sample', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby=['cell_type','sample'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-sample.pdf'),bbox_inches='tight')