import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import anndata as ad
from cycler import cycler
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt
import os
from itertools import combinations
sc.set_figure_params(figsize=(4, 4))
sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
import warnings

# 禁用所有警告
warnings.filterwarnings("ignore")
def save_fig(plt, out_dir, file_name, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()
class vdj_obj:
    def __init__(self,mdata = None,species= None,method = None,celltype= None,cutoff = 15,output = None):
        '''
        mdata : None
        species : mouse or human
        method : TCR or BCR
        cutoff :defalut 15
        output : None
        '''
        self.dat = mdata
        self.species = species
        self.output = output
        self.method = method
        self.cutoff = cutoff
        self.celltype = celltype
        self.tmp_clonotype = None

    def qc(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'1.Quality Control')
        method = self.method
        species = self.species
        ir.pp.index_chains(mdata)
        ir.tl.chain_qc(mdata)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"wspace": 0.5})
        if species == "human":
            if method == "TCR":
                geneuse = "CD3D"
                filename = "T"
            elif method == "BCR":
                geneuse = "CD79A"
                filename = "B"
        elif species == "mouse":
            if method == "TCR":
                geneuse = "Cd3d"
                filename = "T"
            elif method == "BCR":
                geneuse = "Cd79a"
                filename = "B"
        mu.pl.embedding(mdata, basis="gex:umap", color=geneuse, ax=ax0, show=False)
        mu.pl.embedding(mdata, basis="gex:umap", color=["airr:receptor_type"], ax=ax1)
        save_fig(fig,out_dir=outdir,file_name=filename,fig_size={'width':10,'height':5})
        # 1.1质控前结果
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        ir.pl.group_abundance(mdata, groupby="airr:receptor_subtype", target_col="gex:group",ax=axs[0])
        axs[0].set_title('Receptor Subtype')
        ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:group",ax=axs[1])
        axs[1].set_title('Chain Pairing barplot')
        mu.pl.embedding(mdata, basis="gex:umap", color="airr:chain_pairing", groups="multichain",ax=axs[2])
        axs[2].set_title('Chain Pairing Umap')
        fig.tight_layout()
        save_fig(fig,out_dir=outdir,file_name='1.1.QC_raw',fig_size={'width':15,'height':5})
        # 2.进行质量控制过滤
        mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: ~np.isin(x, ["multichain","orphan VDJ", "orphan VJ"]))
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        ir.pl.group_abundance(mdata, groupby="airr:receptor_subtype", target_col="gex:group",ax=axs[0])
        axs[0].set_title('Receptor Subtype')
        ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col="gex:group",ax=axs[1])
        axs[1].set_title('Chain Pairing barplot')
        mu.pl.embedding(mdata, basis="gex:umap", color="airr:chain_pairing",groups = "single pair",ax=axs[2])
        axs[2].set_title('Chain Pairing Umap')
        save_fig(fig,out_dir=outdir,file_name='1.2.QC_filter',fig_size={'width':15,'height':5})
        self.dat = mdata
        return self

    def def_cc(self):
        mdata = self.dat
        cutoff = self.cutoff
        celltype = self.celltype
        outdir = os.path.join(self.output,'2.Define clonotype clusters')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # 按照氨基酸的一致性BLOSUM矩阵用于氨基酸比对的评分（https://en.wikipedia.org/wiki/BLOSUM），一致性评分卡值cutoff
        ir.pp.ir_dist(mdata,metric="alignment",sequence="aa",cutoff=cutoff)
        ir.tl.define_clonotype_clusters(mdata,sequence="aa",metric="alignment",receptor_arms="all",dual_ir="any")
        ir.tl.clonotype_network(mdata,min_cells=3,sequence="aa",metric="alignment")
        if celltype is None:
            if "gex:cell_type" not in mdata.obs.columns.to_list():
                mdata.obs["gex:cell_type"] = mdata['gex'].obs["leiden"]
                mdata['gex'].obs["cell_type"] = mdata['gex'].obs["leiden"]
        else:
            mdata.obs["gex:cell_type"] = mdata['gex'].obs[celltype]
            mdata['gex'].obs["cell_type"] = mdata['gex'].obs[celltype]

        # 创建一个一行三列的图
        fig, axs = plt.subplots(1, 3,figsize = (30,10))
        ir.pl.clonotype_network(mdata, color="gex:sample", ax =axs[0])
        ir.pl.clonotype_network(mdata, color="gex:group",ax =axs[1])
        ir.pl.clonotype_network(mdata, color="gex:cell_type", ax =axs[2])
        save_fig(fig,out_dir=outdir,file_name='2.Network',fig_size={'width':20,'height':8})
        self.dat = mdata
        return self
    
    def clo_exp(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'3.Clonal expansion')
        # 3.clonotype分析
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # 3.1 Clonal expansion统计计算
        ir.tl.clonal_expansion(mdata,target_col='cc_aa_alignment',clip_at=6)
        mu.pl.embedding(mdata, basis="gex:umap", color=["airr:clonal_expansion", "airr:cc_aa_alignment_size",'gex:cell_type'],show=False)
        save_fig(plt,out_dir=outdir,file_name='3.1.Clonotype_analysis_umap',fig_size={'width':15,'height':5})
        fig,axs = plt.subplots(1,3,figsize = (15,5))
        ir.pl.clonal_expansion(mdata,target_col='cc_aa_alignment', groupby="gex:cell_type", clip_at=6, normalize=False,ax = axs[0])
        ir.pl.clonal_expansion(mdata,target_col='cc_aa_alignment', groupby="gex:cell_type", clip_at=6, normalize=True,ax = axs[1])
        ir.pl.alpha_diversity(mdata,target_col='cc_aa_alignment',metric="normalized_shannon_entropy", groupby="gex:cell_type",ax = axs[2])
        axs[2].set_title("")
        save_fig(fig,out_dir=outdir,file_name='3.2.Clonotype_expansion_barplot',fig_size={'width':15,'height':5})
        # 3.1.1 对不同细胞类型的数量进行统计
        clono_size = ir.tl.summarize_clonal_expansion(mdata,target_col='cc_aa_alignment',groupby='gex:group')
        clono_size.to_csv(os.path.join(outdir,'Clonotype_size_group_summary.tsv'),sep="\t")
        clono_size = ir.tl.summarize_clonal_expansion(mdata,target_col='cc_aa_alignment',groupby='gex:cell_type')
        clono_size.to_csv(os.path.join(outdir,'Clonotype_size_cell_type_summary.tsv'),sep="\t")
        # 3.2 Gini系数计算
        data = ir.tl.alpha_diversity(mdata, groupby="gex:cell_type", target_col='cc_aa_alignment', metric='gini_index',inplace=False)
        my_title = {
            "title": "Alpha diversity of cc_aa_alignment by gex:cell_type",
            "ylab": 'gini_index'
        }  
        ax =ir.pl.base.bar(data, style_kws=my_title,fig_kws = {'figsize': (6, 4), 'dpi': 300})
        ax.get_legend().remove()
        plt.gcf().canvas.draw()
        ax
        save_fig(plt,out_dir=outdir,file_name='3.3.Clonotype_expansion_barplot',fig_size={'width':5,'height':5})
        plt.clf() 
        self.dat = mdata
        return self
    def clo_abu(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'4.Clonotype abundance')
        # 4. Clonotype abundance克隆丰度分析
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # 4.1 cell_type
        fig,axs = plt.subplots(1,3,figsize = (15,5))
        ir.pl.group_abundance(mdata,groupby="airr:cc_aa_alignment",target_col="gex:cell_type",max_cols=10,ax = axs[0])
        ir.pl.group_abundance(mdata,groupby="airr:cc_aa_alignment",target_col="gex:cell_type",max_cols=10,normalize="gex:group",ax = axs[1])
        ir.pl.group_abundance(mdata,groupby="airr:cc_aa_alignment",target_col="gex:cell_type",max_cols=10,normalize=True,ax = axs[2])
        save_fig(fig,out_dir=outdir,file_name='4.1.Clonotype_abundance_celltype_barplot',fig_size={'width':15,'height':5})
        # 4.1.1 gene usage cell_type
        fig,axs = plt.subplots(1,2,figsize = (10,5))
        with ir.get.airr_context(mdata,"v_call"):
            ir.pl.group_abundance(mdata,groupby="VJ_1_v_call",target_col="gex:cell_type",max_cols=10,ax = axs[0])
            ir.pl.group_abundance(mdata,groupby="VJ_1_v_call",target_col="gex:cell_type",max_cols=10,normalize=True,ax = axs[1])
        save_fig(fig,out_dir=outdir,file_name='4.1.1.Clonotype_abundance_celltype_VJ_barplot',fig_size={'width':10,'height':5})
        fig,axs = plt.subplots(1,2,figsize = (10,5))
        with ir.get.airr_context(mdata,"v_call"):
            ir.pl.group_abundance(mdata, groupby="VDJ_1_v_call",target_col="gex:cell_type",max_cols=10,ax = axs[0])
            ir.pl.group_abundance(mdata,groupby="VDJ_1_v_call",target_col="gex:cell_type",max_cols=10,normalize=True,ax = axs[1])
        fig.tight_layout()
        save_fig(fig,out_dir=outdir,file_name='4.1.2.Clonotype_abundance_celltype_VDJ_barplot',fig_size={'width':10,'height':5})
        # 4.2 group
        fig,axs = plt.subplots(1,3,figsize = (15,5))
        ir.pl.group_abundance(mdata, groupby="airr:cc_aa_alignment",target_col="gex:group",max_cols=10,ax = axs[0])
        ir.pl.group_abundance(mdata,groupby="airr:cc_aa_alignment",target_col="gex:group",max_cols=10,normalize="gex:cell_type",ax = axs[1])
        ir.pl.group_abundance(mdata,groupby="airr:cc_aa_alignment",target_col="gex:group",max_cols=10,normalize=True,ax = axs[2])
        save_fig(fig,out_dir=outdir,file_name='4.2.Clonotype_abundance_group_barplot',fig_size={'width':15,'height':5})
        # 4.2.gene usage group
        fig,axs = plt.subplots(1,2,figsize = (10,5))
        with ir.get.airr_context(mdata, "v_call"):
            ir.pl.group_abundance(mdata,groupby="VJ_1_v_call",target_col="gex:group",max_cols=10,ax = axs[0])
            ir.pl.group_abundance(mdata,groupby="VJ_1_v_call",target_col="gex:group",max_cols=10,normalize=True,ax = axs[1])
        fig.tight_layout()
        save_fig(fig,out_dir=outdir,file_name='4.2.1.Clonotype_abundance_group_VJ_barplot',fig_size={'width':10,'height':5})
        fig,axs = plt.subplots(1,2,figsize = (10,5))
        with ir.get.airr_context(mdata, "v_call"):
            ir.pl.group_abundance(mdata, groupby="VDJ_1_v_call",target_col="gex:group",max_cols=10,ax = axs[0])
            ir.pl.group_abundance(mdata,groupby="VDJ_1_v_call",target_col="gex:group",max_cols=10,normalize=True,ax = axs[1])
        fig.tight_layout()
        save_fig(fig,out_dir=outdir,file_name='4.2.2.Clonotype_abundance_group_VDJ_barplot',fig_size={'width':10,'height':5})
        # 4.3 umap结果
        if(len(mdata['gex'].obs['sample'].drop_duplicates().to_list())>2):
            tmp_clonotype = mdata['airr'].obs['cc_aa_alignment'].value_counts().head(10).index.tolist()
            mdata.obs['airr:top_clonotype'] = np.where(mdata.obs['airr:cc_aa_alignment'].isin(tmp_clonotype),mdata.obs['airr:cc_aa_alignment'],np.nan)
            mu.pl.embedding(mdata,basis="gex:umap",color='airr:top_clonotype',show = False)
            plt.tight_layout()
            save_fig(plt,out_dir=outdir,file_name='4.3.Clonotype_abundance_umap',fig_size={'width':5,'height':5})
        else:
            tmp_clonotype = mdata['airr'].obs['cc_aa_alignment'].value_counts().head(10).index.tolist()
        # 4.4 vdj组成
        ir.pl.vdj_usage(mdata[mdata.obs['airr:cc_aa_alignment'].isin(tmp_clonotype),:],max_ribbons=10,max_labelled_segments=8   )
        save_fig(plt,out_dir=outdir,file_name='4.4.Top10_sankey_combination',fig_size={'width':8,'height':8})
        plt.clf()
        # 4.5 clonotype marker
        outdir = os.path.join(outdir,'4.5.Clonotype marker')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        with ir.get.obs_context(mdata["gex"], {"clone_id": mdata.obs["airr:cc_aa_alignment"]}) as tmp_ad:
            for item1, item2 in combinations(tmp_clonotype, 2):
                try:
                    sc.tl.rank_genes_groups(
                        tmp_ad, "clone_id", groups=[item1], reference=item2, method="wilcoxon"
                    )
                    sc.pl.rank_genes_groups_violin(tmp_ad, groups=item1, n_genes=10,show=False)
                    plt.tight_layout()
                    save_fig(plt,out_dir=outdir,file_name='_'.join([item1,item2]),fig_size={'width':5,'height':5})
                    plt.clf() 
                except:
                    continue
        self.tmp_clonotype = tmp_clonotype
        self.dat = mdata 
        return self
    
    def clo_siz(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'5.Clonotype size')
        # 5.spectral type diagram CDR3长度分析
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # 5.1 celltype
        fig,axs = plt.subplots(1,3,figsize = (15,5))
        ir.pl.spectratype(mdata,color="gex:cell_type",viztype="bar", ax = axs[0],fig_kws={"dpi": 300})
        axs[0].tick_params(axis='x', rotation=90)
        ir.pl.spectratype(mdata,color="gex:cell_type",viztype="curve",ax = axs[1],fig_kws={"dpi": 300},kde_kws={"kde_norm": False})
        ir.pl.spectratype(mdata,color="gex:cell_type",viztype="curve",curve_layout="shifted",ax = axs[2],fig_kws={"dpi": 300},kde_kws={"kde_norm": False})
        fig.tight_layout()
        save_fig(fig,out_dir=outdir,file_name='5.1.CDR3_length_celltype',fig_size={'width':15,'height':5})
        # 5.2 group
        fig,axs = plt.subplots(1,3,figsize = (15,5))
        ir.pl.spectratype(mdata,color="gex:group",viztype="bar",ax = axs[0],fig_kws={"dpi": 300})
        axs[0].tick_params(axis='x', rotation=90)
        ir.pl.spectratype(mdata,color="gex:group",viztype="curve",ax = axs[1],fig_kws={"dpi": 300},kde_kws={"kde_norm": False})
        ir.pl.spectratype(mdata,color="gex:group",viztype="curve",curve_layout="shifted",ax = axs[2],fig_kws={"dpi": 300},kde_kws={"kde_norm": False})
        fig.tight_layout()
        save_fig(fig,out_dir=outdir,file_name='5.2.CDR3_length_group',fig_size={'width':15,'height':5})
        self.dat = mdata 
        return self
    
    def com_rep(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'6.Comparing repertoires')
        # 6.Comparing repertoires表达模式比较
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if(len(mdata['gex'].obs['sample'].drop_duplicates().to_list())>1):
            df, dst, lk = ir.tl.repertoire_overlap(mdata, groupby="gex:sample", target_col = "cc_aa_alignment",inplace=False)
            df.head()
            ir.pl.repertoire_overlap(mdata,groupby = "gex:sample",target_col = "cc_aa_alignment",heatmap_cats=["gex:group"],yticklabels=True,xticklabels=True)
            save_fig(plt,out_dir=outdir,file_name='6.1.Repertoire_heatmap',fig_size={'width':10,'height':5})
            # 6.2 不同样本间的比较
            outdir1 = os.path.join(outdir,'6.2Comparing_repertoires_group')
            if not os.path.exists(outdir1):
                os.mkdir(outdir1)
            group = mdata['gex'].obs['group'].drop_duplicates().to_list()
            for item1, item2 in combinations(group, 2):
                tmp_group = [item1, item2]
                ir.pl.repertoire_overlap(mdata, "gex:group",pair_to_plot=tmp_group,target_col = "cc_aa_alignment")
                save_fig(plt,out_dir=outdir1,file_name='_'.join(tmp_group),fig_size={'width':5,'height':5})
            # 6.3 不同细胞类型间的比较
            outdir2 = os.path.join(os.path.join(outdir,'6.3Comparing_repertoires_celltype'))
            if not os.path.exists(outdir2):
                os.mkdir(outdir2)
            cell_type = mdata['gex'].obs['cell_type'].drop_duplicates().to_list()
            for item1, item2 in combinations(cell_type, 2):
                tmp_group = [item1, item2]
                ir.pl.repertoire_overlap(mdata, "gex:cell_type",pair_to_plot=tmp_group,target_col = "cc_aa_alignment")
                save_fig(plt,out_dir=outdir2,file_name='_'.join(tmp_group),fig_size={'width':5,'height':5})
        self.dat = mdata 
        return self
    
    def clo_modu(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'7.Clonotype modularity')
        # 7.Clonotype modularity计算ClonoCluster相似性
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        sc.pp.neighbors(mdata["gex"])
        if(len(mdata['gex'].obs['sample'].drop_duplicates().to_list())>2):
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            ir.tl.clonotype_modularity(mdata, target_col="cc_aa_alignment")
            ir.pl.clonotype_network(mdata,color="airr:clonotype_modularity",label_fontsize=9,panel_size=(6, 6),base_size=20,ax=axs[0])
            mu.pl.embedding(mdata, basis="gex:umap", color="airr:clonotype_modularity",ax = axs[1])
            save_fig(fig,out_dir=outdir,file_name='7.1.Clonotype_modularity_score',fig_size={'width':10,'height':5})
            ir.pl.clonotype_modularity(mdata, base_size=10,labels_quantile_cutoff = (0.85, 0.95))
            save_fig(plt,out_dir=outdir,file_name='7.2.Clonotype_modularity_volcano',fig_size={'width':10,'height':10})
            # 7.3绘制相似性前10
            clonotypes_top_modularity = list(mdata.obs.set_index("airr:cc_aa_alignment")["airr:clonotype_modularity"].sort_values(ascending=False).index.unique().values[:10])
            mdata.obs['airr:top10_modularity'] = np.where(mdata.obs['airr:cc_aa_alignment'].isin(clonotypes_top_modularity),mdata.obs['airr:cc_aa_alignment'],np.nan)
            mu.pl.embedding(mdata,basis="gex:umap",color="airr:top10_modularity",show=False)
            save_fig(plt,out_dir=outdir,file_name='7.3.Top10_Clonotype_modularity_umap',fig_size={'width':10,'height':10})
            # 7.4 violin_plot
            with ir.get.obs_context(mdata["gex"], {"cc_aa_alignment": mdata.obs["airr:cc_aa_alignment"]}) as tmp_ad:
                sc.tl.rank_genes_groups(tmp_ad,"cc_aa_alignment",groups=clonotypes_top_modularity,reference="rest",method="wilcoxon")
                fig, axs = plt.subplots(2, 5, figsize=(20, 13))
                for i, ct in enumerate(clonotypes_top_modularity):
                    ax = axs.flatten()[i]
                    sc.pl.rank_genes_groups_violin(tmp_ad, groups=[ct], n_genes=10, ax=ax, show=False, strip=False)
            fig.tight_layout()
            save_fig(fig,out_dir=outdir,file_name='7.4.Top10_Clonotype_modularity_violin',fig_size={'width':20,'height':13})
        self.dat = mdata
        return self
    
    def vdj_pre(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'8.VDJDB Predict')
        # 8.VDJDB Predict
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        vdjdb = sc.read_h5ad('/lvdata/wzb/pipline/scirpy/data/vdjdb.h5ad')
        ir.pp.ir_dist(mdata, vdjdb, metric="identity", sequence="aa")
        ir.tl.ir_query(mdata,vdjdb,metric="identity",sequence="aa",receptor_arms="any",dual_ir="any")
        df = ir.tl.ir_query_annotate_df(mdata,vdjdb,metric="identity",sequence="aa",include_ref_cols=["antigen.epitope", "antigen.gene", "antigen.species", "mhc.a", "mhc.b", "mhc.class", "reference.id", "vdjdb.score"])
        df.to_csv(os.path.join(outdir,'VDJDB_result.tsv'),sep="\t")
        ir.tl.ir_query_annotate(mdata,vdjdb,metric="identity",sequence="aa",include_ref_cols=["antigen.species"],strategy="most-frequent")
        
        mu.pl.embedding(mdata, "gex:umap", color="airr:antigen.species",show=False)
        save_fig(plt,out_dir=outdir,file_name='VDJDB_epitope_databases',fig_size={'width':15,'height':5})
        self.dat = mdata
        return self
    
    def mot_ana(self):
        mdata = self.dat
        tmp_clonotype = self.tmp_clonotype
        outdir = os.path.join(self.output,'9.Motif Analysis')
        # 9.Motif Analysis
        from palmotif import compute_pal_motif, mpl_logo, svg_logo,compute_motif
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        with ir.get.airr_context(mdata, "junction_aa", ["VJ_1", "VDJ_1", "VJ_2", "VDJ_2"]):
            for i in tmp_clonotype:
                cdr3_tmp = (mdata.obs.loc[lambda x: x["airr:cc_aa_alignment"] == i].astype(str).groupby(
                        ["VJ_1_junction_aa","VDJ_1_junction_aa","VJ_2_junction_aa","VDJ_2_junction_aa","airr:receptor_subtype"],
                        observed=True,
                        dropna=False).size().reset_index(name="n_cells"))
                tmp_VJ = cdr3_tmp['VJ_1_junction_aa']
                tmp_VDJ = cdr3_tmp['VDJ_1_junction_aa']
                centroid_VJ = cdr3_tmp["VJ_1_junction_aa"].value_counts().index[0]
                centroid_VDJ = cdr3_tmp["VDJ_1_junction_aa"].value_counts().index[0]
                vjplot,_=compute_pal_motif(seqs=tmp_VJ,centroid=centroid_VJ)
                vdjplot,_=compute_pal_motif(seqs=tmp_VDJ,centroid=centroid_VDJ)
                # 无序列算法
                svg_logo(vjplot,f'{outdir}/{i}_vj.svg')
                svg_logo(vdjplot,f'{outdir}/{i}_vdj.svg')
        self.dat = mdata
        return self

    def sha_clo(self):
        mdata = self.dat
        outdir = os.path.join(self.output,'10.Share Clonotype')
        # 10.Share Clonotype
        import seaborn as sns
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        mdata['airr'].obs[['cell_type','sample','group']] = mdata['gex'].obs[['cell_type','sample','group']]
        if(len(mdata['gex'].obs['sample'].drop_duplicates().to_list())>1):
            df_all = mdata['airr'].obs[['cell_type','sample','group','cc_aa_alignment']]
            heatlist = ['cell_type','sample','group']
            for p in heatlist:
                all_clonotype_size = df_all[p].value_counts()*2
                all_clonotype_size = all_clonotype_size.to_dict()
                vc_clone = df_all.cc_aa_alignment.value_counts()
                multiples = [k for k, a in vc_clone.items() if a > 1]
                df = df_all[df_all.cc_aa_alignment.isin(multiples)]
                clusters = df[p].unique().tolist() #cell_type信息
                size_clusters = {k: v for k, v in df[p].value_counts().items()}
                gpby_clonotype = df.groupby("cc_aa_alignment")[p].apply(list)
                gpby_clonotype = gpby_clonotype[gpby_clonotype.apply(len) > 1].copy()
                proba = {}
                for c1 in clusters:
                    tmp = clusters[:]
                    tmp.remove(c1)
                    for c2 in tmp:
                        proba[c1, c2] = np.sum(gpby_clonotype.apply(lambda x: x.count(c1)+x.count(c2)))
                df_shared_clonotype = pd.Series(proba)
                df_shared_clonotype
                L = clusters
                mux = pd.MultiIndex.from_product([L, L], names=df_shared_clonotype.index.names)
                df_shared_clonotype = df_shared_clonotype.reindex(mux)
                df_ustack = df_shared_clonotype.unstack()
                df_ustack = df_ustack[L]
                df_ustack = df_ustack.reindex(L)
                df_ustack = df_ustack.fillna(all_clonotype_size)
                sns.reset_defaults()
                pl = sns.clustermap(df_ustack,cmap="Blues",square=True,annot=True,fmt=".0f",col_cluster=False,row_cluster=False,linewidth=.5)
                save_fig(plt,out_dir=outdir,file_name=f'Heatmap_{p}_cell_withnumber',fig_size={'width':10,'height':5})
                pl = sns.clustermap(df_ustack,cmap="Blues",square=True,col_cluster=False,row_cluster=False,linewidth=.5)
                save_fig(plt,out_dir=outdir,file_name=f'Heatmap_{p}_cell_withoutnumber',fig_size={'width':10,'height':5})
        self.dat = mdata
        return self


    def save_h5mu(self):
        mdata = self.dat
        outdir = self.output
        import os
        os.environ["HDF5_OBJECT_HEADER_MAX_SIZE"] = "1000000"
        for i in ["v_call","d_call","j_call","c_call","cdr3","cdr3_aa"]:
            with ir.get.airr_context(mdata,i): 
                tmp = mdata.obs
            for j in ['VJ_1', 'VDJ_1', 'VJ_2', 'VDJ_2']:
                tmp_obs = j+"_"+i
                mdata.obs[tmp_obs] = tmp[tmp_obs]
        mdata.obs.to_csv(os.path.join(outdir,"Result.xls"),sep="\t")
        mdata.obs = mdata.obs.loc[:, mdata.obs.notna().any(axis=0)]
        mdata.write_h5mu(os.path.join(outdir,"Result.h5mu"))