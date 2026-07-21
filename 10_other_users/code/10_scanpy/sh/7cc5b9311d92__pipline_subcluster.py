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

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None):
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    if fig_size is None:
        fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'Fibroblasts', 'Pericytes', 'Endothelial cells', 'AT1', 'Cillated Epithelial cells', 'Mesothelial cell', 'SMC', 'Myeloid cells', 'Neutrophils', 'T cells', 'AT2', 'Plasma cells', 'Myofibroblasts', 'B cells', 'NK cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512127429-01/06.Mission/Myeloid/"
adata=sc.read("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/2026.03.09/Myeloid_cells/celltype/result.h5ad")

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

genes=['Abca1','Abi1','Acvr1b','Acvr2a','Adgre1','Adm','Adora2b','Adrm1','Ahr','Aplnr','Aqp9','Atp2a2','Atp2b1','Atp2c1','Axl','Bdkrb1','Best1','Bst2','Btg2','C3ar1','Calcrl','Ccl17','Ccl20','Ccl22','Ccl24','Ccl5','Ccl7','Ccr7','Ccrl2','Cd14','Cd40','Cd48','Cd69','Cd70','Cd82','Cdkn1a','Chst2','Clec5a','Cmklr1','Csf1','Csf3','Csf3r','Cx3cl1','Cxcl10','Cxcl11','Cxcl15','Cxcl5','Cxcl9','Cxcr6','Cybb','Dcbld2','Ebi3','Edn1','Eif2ak2','Emp3','Ereg','F3','Ffar2','Fpr1','Fzd5','Gabbr1','Gch1','Gna15','Gnai3','Gp1ba','Gpc3','Gpr132','Gpr183','Has2','Hbegf','Hif1a','Hpn','Hrh1','Icam1','Icam4','Icosl','Ifitm1','Ifnar1','Ifngr2','Il10','Il10ra','Il12b','Il15','Il15ra','Il18','Il18r1','Il18rap','Il1a','Il1b','Il1r1','Il2rb','Il4ra','Il6','Il7r','Inhba','Irak2','Irf1','Irf7','Itga5','Itgb3','Itgb8','Kcna3','Kcnj2','Kcnmb2','Kif1b','Klf6','Lamp3','Lck','Lcp2','Ldlr','Lif','Lpar1','Lta','Ly6e','Lyn','Marco','Mefv','Mep1a','Met','Mmp14','Msr1','Mxd1','Myc','Nampt','Ndp','Nfkb1','Nfkbia','Nlrp3','Nmi','Nmur1','Nod2','Npffr2','Olr1','Oprk1','Osm','Osmr','P2rx4','P2rx7','P2ry2','Pcdh7','Pde4b','Pdpn','Pik3r5','Plaur','Prok2','Psen1','Ptafr','Ptger2','Ptger4','Ptgir','Ptpre','Pvr','Raf1','Rasgrp1','Rela','Rgs1','Rgs16','Rhog','Ripk2','Rnf144b','Ros1','Rtp4','Scarf1','Scn1b','Sele','Selenos','Sell','Sema4d','Serpine1','Sgms2','Slamf1','Slc11a2','Slc1a2','Slc28a2','Slc31a1','Slc31a2','Slc4a4','Slc7a1','Slc7a2','Sphk1','Sri','Stab1','Tacr1','Tacr3','Tapbp','Timp1','Tlr1','Tlr2','Tlr3','Tnfaip6','Tnfrsf1b','Tnfrsf9','Tnfsf10','Tnfsf15','Tnfsf9','Tpbg','Vip']
sc.tl.score_genes(adata, gene_list=genes, score_name='HALLMARK_INFLAMMATORY_RESPONSE')
sc.pl.umap(adata, color='HALLMARK_INFLAMMATORY_RESPONSE', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_HALLMARK_INFLAMMATORY_RESPONSE', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='HALLMARK_INFLAMMATORY_RESPONSE', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_HALLMARK_INFLAMMATORY_RESPONSE', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['HALLMARK_INFLAMMATORY_RESPONSE'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_HALLMARK_INFLAMMATORY_RESPONSE', fig_size={'w': 6, 'h': 6})
genes=['Aadat','Acaa1a','Acaa2','Acadl','Acadm','Acads','Acadvl','Acat2','Aco2','Acot2','Acot8','Acox1','Acsl1','Acsl4','Acsl5','Acsm3','Acss1','Adh7','Adipor2','Adsl','Alad','Aldh1a1','Aldh3a1','Aldh3a2','Aldh9a1','Aldoa','Aoc3','Apex1','Aqp7','Auh','Bckdhb','Blvra','Bmpr1b','Bphl','Car2','Car4','Car6','Cbr1','Cbr3','Cd1d2','Cd36','Cel','Cidea','Cpox','Cpt1a','Cpt2','Crat','Cryz','Cyp1a1','Cyp4a10','D2hgdh','Decr1','Dhcr24','Dld','Dlst','Ech1','Echs1','Eci1','Eci2','Ehhadh','Elovl5','Eno2','Eno3','Ephx1','Erp29','Etfdh','Fabp1','Fabp2','Fasn','Fh1','Fmo1','G0s2','Gabarapl1','Gad2','Gapdhs','Gcdh','Glul','Gpd1','Gpd2','Grhpr','Gstz1','H2az1','Hadh','Hadhb','Hao2','Hccs','Hibch','Hmgcl','Hmgcs1','Hmgcs2','Hpgd','Hsd17b10','Hsd17b11','Hsd17b4','Hsd17b7','Hsdl2','Hsp90aa1','Hsph1','Idh1','Idh3b','Idh3g','Idi1','Inmt','Kmt5a','Ldha','Lgals1','Ltc4s','Maoa','Mcee','Mdh1','Mdh2','Me1','Metap1','Mgll','Mif','Mix23','Mlycd','Nbn','Ncaph2','Nsdhl','Nthl1','Odc1','Ostc','Pcbd1','Pdha1','Pdhb','Ppara','Prdx6','Psme1','Ptprg','Pts','Rap1gds1','Rdh1','Rdh11','Reep6','Retsat','S100a10','Sdha','Sdhc','Sdhd','Serinc1','Slc22a5','Sms-ps','Sucla2','Suclg1','Suclg2','Tdo2','Trp53inp2','Ube2l6','Ugdh','Urod','Uros','Vnn1','Xist','Ywhah']
sc.tl.score_genes(adata, gene_list=genes, score_name='HALLMARK_FATTY_ACID_METABOLISM')
sc.pl.umap(adata, color='HALLMARK_FATTY_ACID_METABOLISM', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_HALLMARK_FATTY_ACID_METABOLISM', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='HALLMARK_FATTY_ACID_METABOLISM', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_HALLMARK_FATTY_ACID_METABOLISM', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['HALLMARK_FATTY_ACID_METABOLISM'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_HALLMARK_FATTY_ACID_METABOLISM', fig_size={'w': 6, 'h': 6})



# 对Fibroblast注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'Fibroblasts', 'Pericytes', 'Endothelial cells', 'AT1', 'Cillated Epithelial cells', 'Mesothelial cell', 'SMC', 'Myeloid cells', 'Neutrophils', 'T cells', 'AT2', 'Plasma cells', 'Myofibroblasts', 'B cells', 'NK cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512127429-01/06.Mission/Fibroblasts/"
adata=sc.read("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/2026.01.21/Topic4/Fibroblasts/celltype/result.h5ad")

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

genes=['Abi3bp','Acta2','Adam12','Anpep','Aplp1','Areg','Basp1','Bdnf','Bgn','Bmp1','Cadm1','Cald1','Calu','Cap2','Capg','Ccn1','Ccn2','Cd44','Cdh11','Cdh2','Cdh6','Col11a1','Col12a1','Col16a1','Col1a1','Col1a2','Col3a1','Col4a1','Col4a2','Col5a1','Col5a2','Col5a3','Col6a2','Col6a3','Col7a1','Col8a2','Colgalt1','Comp','Copa','Crlf1','Cthrc1','Cxcl12','Cxcl15','Cxcl5','Dab2','Dcn','Dkk1','Dpysl3','Dst','Ecm1','Ecm2','Edil3','Efemp2','Eln','Emp3','Eno2','Fap','Fas','Fbln1','Fbln2','Fbln5','Fbn1','Fbn2','Fermt2','Fgf2','Flna','Fmod','Fn1','Foxc2','Fstl1','Fstl3','Fuca1','Fzd8','Gadd45a','Gadd45b','Gas1','Gem','Gja1','Glipr1','Gm21451','Gpc1','Gpx7','Grem1','Htra1','Id2','Igfbp2','Igfbp3','Igfbp4','Il15','Il6','Inhba','Itga2','Itga5','Itgav','Itgb1','Itgb3','Itgb5','Jun','Lama1','Lama2','Lama3','Lamc1','Lamc2','Lgals1','Lox','Loxl1','Lrp1','Lrrc15','Lum','Magee1','Matn2','Matn3','Mcm7','Mest','Mfap5','Mgp','Mmp14','Mmp2','Mmp3','Msx1','Myl9','Mylk','Nid2','Nnmt','Notch2','Nt5e','Ntm','Oxtr','P3h1','Pcolce','Pcolce2','Pdgfrb','Pdlim4','Pfn2','Plaur','Plod1','Plod2','Plod3','Pmepa1','Pmp22','Postn','Ppib','Prrx1','Pthlh','Ptx3','Pvr','Qsox1','Rgs4','Rhob','Sat1','Scg2','Sdc1','Sdc4','Serpine1','Serpine2','Serpinh1','Sfrp1','Sfrp4','Sgcb','Sgcd','Sgcg','Slc6a8','Slit2','Slit3','Snai2','Sntb1','Sparc','Spock1','Spp1','Tagln','Tfpi2','Tgfb1','Tgfbi','Tgfbr3','Tgm2','Thbs1','Thbs2','Thy1','Timp1','Timp3','Tnc','Tnfaip3','Tnfrsf11b','Tnfrsf12a','Tpm1','Tpm2','Tpm4','Vcam1','Vcan','Vegfa','Vegfc','Vim','Wipf1','Wnt5a']
sc.tl.score_genes(adata, gene_list=genes, score_name='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
sc.pl.umap(adata, color='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', fig_size={'w': 6, 'h': 6})
genes=['Acvr1','Apc','Arid4b','Bcar3','Bmp2','Bmpr1a','Bmpr2','Cdh1','Cdk9','Cdkn1c','Ctnnb1','Eng','Fkbp1a','Fnta','Furin','Hdac1','Hipk2','Id1','Id2','Id3','Ifngr2','Junb','Klf10','Ltbp2','Map3k7','Ncor2','Nog','Pmepa1','Ppm1a','Ppp1ca','Ppp1r15a','Rab31','Rhoa','Serpine1','Ski','Skil','Slc20a1','Smad1','Smad3','Smad6','Smad7','Smurf1','Smurf2','Sptbn1','Tgfb1','Tgfbr1','Tgif1','Thbs1','Tjp1','Trim33','Ube2d3','Wwtr1','Xiap']
sc.tl.score_genes(adata, gene_list=genes, score_name='HALLMARK_TGF_BETA_SIGNALING')
sc.pl.umap(adata, color='HALLMARK_TGF_BETA_SIGNALING', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_HALLMARK_TGF_BETA_SIGNALING', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='HALLMARK_TGF_BETA_SIGNALING', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_HALLMARK_TGF_BETA_SIGNALING', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['HALLMARK_TGF_BETA_SIGNALING'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_HALLMARK_TGF_BETA_SIGNALING', fig_size={'w': 6, 'h': 6})
genes=['Adam17','Axin1','Axin2','Ccnd2','Csnk1e','Ctnnb1','Cul1','Dkk1','Dkk4','Dll1','Dvl2','Frat1','Fzd1','Fzd8','Gnai1','Hdac11','Hdac2','Hdac5','Hey1','Hey2','Jag1','Jag2','Kat2a','Lef1','Maml1','Myc','Ncor2','Ncstn','Nkd1','Notch1','Notch4','Numb','Ppard','Psen2','Ptch1','Rbpj','Skp2','Tcf7','Trp53','Wnt1','Wnt5b','Wnt6']
sc.tl.score_genes(adata, gene_list=genes, score_name='HALLMARK_WNT_BETA_CATENIN_SIGNALING')
sc.pl.umap(adata, color='HALLMARK_WNT_BETA_CATENIN_SIGNALING', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_HALLMARK_WNT_BETA_CATENIN_SIGNALING', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='HALLMARK_WNT_BETA_CATENIN_SIGNALING', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_HALLMARK_WNT_BETA_CATENIN_SIGNALING', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['HALLMARK_WNT_BETA_CATENIN_SIGNALING'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_HALLMARK_WNT_BETA_CATENIN_SIGNALING', fig_size={'w': 6, 'h': 6})


# 对Endothelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'Fibroblasts', 'Pericytes', 'Endothelial cells', 'AT1', 'Cillated Epithelial cells', 'Mesothelial cell', 'SMC', 'Myeloid cells', 'Neutrophils', 'T cells', 'AT2', 'Plasma cells', 'Myofibroblasts', 'B cells', 'NK cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512127429-01/06.Mission/Endothelial/"
adata=sc.read("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/2026.01.05/Endothelial_cells/celltype/result.h5ad")

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
    

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512127429-01/06.Mission/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/2026.01.21/Topic4/Fibroblasts/celltype/result.h5ad")
bdata = sc.read("/mnt/NexData_file/home/wzb/backup/scRNA/SC2512127429/2026.03.09/Myeloid_cells/celltype/result.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')