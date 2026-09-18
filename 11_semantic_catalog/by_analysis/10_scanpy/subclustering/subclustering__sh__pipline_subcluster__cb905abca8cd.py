import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
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
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512027077-01/03.Mission/gene_plot/"
sc.pl.umap(adata, color=["Col1a1","Col1a2","Col2a1","Col6a1","Hbegf","Ibsp","Thbs1","Thbs2"], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=["Col1a1","Col1a2","Col2a1","Col6a1","Hbegf","Ibsp","Thbs1","Thbs2"], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Col1a1","Col1a2","Col2a1","Col6a1","Hbegf","Ibsp","Thbs1","Thbs2"], groupby='group', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["Col1a1","Col1a2","Col2a1","Col6a1","Hbegf","Ibsp","Thbs1","Thbs2"], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.pdf'),bbox_inches='tight')

adata_sub=adata[adata.obs.cell_type=="Osteoblast",]
sc.pl.dotplot(adata, var_names=["Col1a1","Col1a2","Col2a1","Col6a1","Hbegf","Ibsp","Thbs1","Thbs2"], groupby='sample', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.pdf'),bbox_inches='tight')

ECM_receptor_interaction=['Col1a1','Col1a2','Col2a1','Col4a1','Col4a2','Col4a3','Col4a4','Col4a5','Col4a6','Col6a1','Col6a2','Col6a3','Col6a6','Col6a5','Col6a4','Col9a1','Col9a2','Col9a3','Lama1','Lama2','Lama3','Lama5','Lama4','Lamb1','Lamb2','Lamb3','Lamc1','Lamc2','Lamc3','Chad','Reln','Thbs1','Comp','Thbs2','Thbs3','Thbs4','Fn1','Spp1','Vtn','Tnc','Tnr','Tnn','Tnxb','Npnt','Fras1','Frem2','Frem1','Dspp','Vwf','Ibsp','Dmp1','Mepe','Agrn','Hspg2','Itga1','Itga2','Itga2b','Itga3','Itga4','Itga5','Itga6','Itga7','Itga8','Itga9','Itga10','Itga11','Itgav','Itgb1','Itgb3','Itgb4','Itgb5','Itgb6','Itgb7','Itgb8','Cd44','Sdc1','Sdc4','Sv2a','Sv2b','Sv2c','Cd36','Gp5','Gp1ba','Gp1bb','Gp9','Gp6','Dag1','Cd47','Hmmr']
ErbB_signaling_pathway=['Egf','Tgfa','Areg','Egfr','Erbb2','Plcg1','Plcg2','Camk2d','Camk2a','Camk2b','Camk2g','Prkca','Prkcb','Prkcg','Cbl','Cblb','Stat5a','Stat5b','Src','Ptk2','Crk','Crkl','Abl1','Abl2','Nck1','Nck2','Pak1','Pak2','Pak3','Pak4','Pak5','Pak6','Map2k4','Map2k7','Mapk10','Mapk8','Mapk9','Jun','Elk1','Btc','Hbegf','Ereg','Erbb3','Nrg1','Nrg2','Erbb4','Shc1','Shc2','Shc3','Shc4','Grb2','Sos1','Sos2','Hras','Kras','Nras','Araf','Braf','Raf1','Map2k1','Map2k2','Mapk1','Mapk3','Myc','Gab1','Nrg3','Nrg4','Pik3ca','Pik3cd','Pik3cb','Pik3r1','Pik3r2','Pik3r3','Akt1','Akt2','Akt3','Mtor','Rps6kb2','Rps6kb1','Eif4ebp1','Bad','Gsk3b','Cdkn1b','Cdkn1a']
sc.tl.score_genes(adata, gene_list=ECM_receptor_interaction, score_name='ECM_receptor_interaction')
sc.tl.score_genes(adata, gene_list=ErbB_signaling_pathway, score_name='ErbB_signaling_pathway')

out_dir="/mnt/NexData_file/home/sh/Projects/AE2512027077-01/03.Mission/genesetscore_plot/"