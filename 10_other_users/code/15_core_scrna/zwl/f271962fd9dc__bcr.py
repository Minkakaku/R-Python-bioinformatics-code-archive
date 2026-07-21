import sys
sys.path.append("/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/script/")
from sciana import *
import scanpy as sc
import os
import anndata as ad
import scirpy as ir
import muon as mu
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from mudata import MuData
import os
import re

filedir = [
    os.path.join("/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/vdjdata/", i, 'filtered_contig_annotations.csv')
    for i in os.listdir("/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/vdjdata/")
]

# 先读入第一个文件
adata_tcr = ir.io.read_10x_vdj(filedir[0])

obs_name_list = adata_tcr.obs_names.to_list()
obs_name_list = [i+("_P2") for i in obs_name_list]
adata_tcr.obs_names = obs_name_list

adata_tcr.var_names_make_unique()

# 依次合并后续文件
for path in filedir[1:]:
    adata_tcr_tmp = ir.io.read_10x_vdj(path)
    sample = re.search(r'vdjdata/(\w+)/', path).group(1)
    obs_name_list = adata_tcr_tmp.obs_names.to_list()
    obs_name_list = [i+("_")+sample for i in obs_name_list]
    adata_tcr_tmp.obs_names = obs_name_list
    adata_tcr_tmp.var_names_make_unique()
    adata_tcr = ad.concat([adata_tcr, adata_tcr_tmp], join='outer', index_unique=None)

# 保证 obs 唯一
adata_tcr.obs_names_make_unique()

adata_tcr.shape

adata=sc.read_h5ad("/mnt/NexData_file/home/lst/project/FWAE20252674/subcluster/B_P/celltype/B_P/adata_anno.h5ad")
# meta=pd.read_csv('/mnt/NexData_file/home/lst/project/FWSC20240776/售后7/CD8T_meta.xls',sep='\t',index_col=0)
# adata.obs=meta
# adata.uns['log1p']['base'] = 2
# adata_cd8=adata[adata.obs['cell_type']!='Other']
obs_name_list = adata.obs_names.astype(str)+("_")+adata.obs['sample'].astype(str)
adata.obs_names = obs_name_list

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR_test/"
with rc_context({'figure.figsize': (10, 6), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight'}):
    sc.pl.umap(adata,color=['cell_type','group'])
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata1=adata
adata1 = adata1.raw.to_adata()
adata1.shape

# set(adata1.obs['cell_type'])

# umi_threshold=50
# adata =sc.read_h5ad('/mnt/NexData_file/home/lst/project/FWSC20240776/售后4/NPT/LN/整合/adata_anno.h5ad')
# adata=adata[~adata.obs['cell_type'].isin(['Myeloid','B cells'])]
# adata.obs["type"] = adata.obs["NPT_umi"].apply(lambda x: "NPT" if x > umi_threshold else "non-NPT")
# #adata.obs["type"] = adata.obs.apply(lambda x: x['type'] if x['cell_type']!='Myeloid' else "non-NPT",axis=1)
# #adata.obs["type"] = adata.obs.apply(lambda x: x['type'] if x['cell_type']!='B cells' else "non-NPT",axis=1)
# adata1=adata[adata.obs['type']=='NPT']
# adata1 = adata1.raw.to_adata()
# adata1.shape

# umi_threshold=50
# adata =sc.read_h5ad('/mnt/NexData_file/home/lst/project/FWSC20240776/售后5/recluster/L/celltype/adata_anno.h5ad')
# adata=adata[adata.obs['cell_type']!='Myeloid cells']
# adata.obs["type"] = adata.obs["NPT_umi"].apply(lambda x: "NPT" if x > umi_threshold else "non-NPT")
# #adata.obs["type"] = adata.obs.apply(lambda x: x['type'] if x['cell_type']!='Myeloid' else "non-NPT",axis=1)
# #adata.obs["type"] = adata.obs.apply(lambda x: x['type'] if x['cell_type']!='B cells' else "non-NPT",axis=1)
# adata1=adata[adata.obs['type']=='NPT']
# adata1 = adata1.raw.to_adata()
# adata1.shape

mdata = MuData({'gex':adata1,'airr':adata_tcr})

mu.pp.intersect_obs(mdata)

os.environ["HDF5_OBJECT_HEADER_MAX_SIZE"] = "1000000"
vdj = vdj_obj(mdata = mdata,species= "human",method = "BCR",cutoff = None,celltype = 'cell_type',output = "/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR/")
vdj = vdj.qc()
vdj = vdj.def_cc()
def save_fig(plt, out_dir, file_name, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()
def def_cc(vdjobj):
    mdata = vdjobj.dat
    cutoff = 15
    celltype = vdjobj.celltype
    outdir = os.path.join(vdjobj.output,'2.Define clonotype clusters')
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
    fig, axs = plt.subplots(1, 3,figsize = (36,12))
    ir.pl.clonotype_network(mdata, color="gex:sample", ax =axs[0],base_size=20,label_fontsize=9,panel_size=(7,7))
    ir.pl.clonotype_network(mdata, color="gex:group",ax =axs[1],base_size=20,label_fontsize=9,panel_size=(7,7))
    ir.pl.clonotype_network(mdata, color="gex:cell_type", ax =axs[2],base_size=20,label_fontsize=9,panel_size=(7,7))
    save_fig(fig,out_dir=outdir,file_name='2.Network',fig_size={'width':27,'height':10})
    vdjobj.dat = mdata
    return vdjobj
vdj = def_cc(vdj)
vdj = vdj.clo_exp()
vdj = vdj.clo_abu()
vdj = vdj.clo_siz()
vdj = vdj.com_rep()
vdj = vdj.clo_modu()
def clo_modu(vdjobj):
    mdata = vdjobj.dat
    outdir = os.path.join(vdjobj.output,'7.Clonotype modularity')
    # 7.Clonotype modularity计算ClonoCluster相似性
    sc.pp.neighbors(mdata["gex"])
    if(len(mdata['gex'].obs['sample'].drop_duplicates().to_list())>2):
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        ir.tl.clonotype_modularity(mdata, target_col="cc_aa_alignment")
        ir.pl.clonotype_network(mdata,color="airr:clonotype_modularity",label_fontsize=9,panel_size=(6, 6),base_size=5,ax=axs[0])
        mu.pl.embedding(mdata, basis="gex:umap", color="airr:clonotype_modularity",ax = axs[1])
        save_fig(fig,out_dir=outdir,file_name='7.1.Clonotype_modularity_score',fig_size={'width':10,'height':5})
        vdjobj.dat = mdata
        return vdjobj
vdj = clo_modu(vdj)
vdj = vdj.vdj_pre()
vdj = vdj.mot_ana()
vdj = vdj.sha_clo()
vdj.save_h5mu()


os.mkdir('/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR/0.Cellranger')
ir.io.write_airr(adata_tcr,os.path.join('/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR/0.Cellranger','filtered_contig_annotations.tsv'))

output = '/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR/'
pid = "FWAE20252674"
ptitle = "BCR"
paddress = "浙江大学医学院附属第一医院"
pcustmor = "姚秋妹"
sample_num = "13"
rmp = ''
with open("/lvdata/wzb/pipline/scirpy/tools/page/tools/vdj/report.rmd",'r',encoding='utf-8') as f:
    line = f.read()
    rmd = line.replace('%pid%', pid).replace('%ptitle%', ptitle).replace('%paddress%', paddress).replace('%pcustmor%', pcustmor).replace('%num_sample%', str(sample_num))
with open(os.path.join(output,'report.rmd'),'w',encoding = 'utf-8') as w:
    w.write(rmd)
f.close()
w.close()
#os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(output))

os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(output))

# adata1 =sc.read_h5ad('/mnt/NexData_file/home/lst/project/FWSC20240776/售后5/recluster/L/celltype/adata_anno.h5ad')
# adata1=adata[adata.obs['cell_type']!='Myeloid cells']
# adata1 = adata1.raw.to_adata()
# adata1.shape

# dirlist = os.listdir('/mnt/NexData_file/home/lst/project/FWSC20240776/售后6/VDJ_df/L')
# dirlist
# filedir = []
# for i in dirlist:
#     filedir.append(os.path.join('/mnt/NexData_file/home/lst/project/FWSC20240776/售后6/VDJ_df/L',i,'filtered_contig_annotations.csv'))
# filedir
# adata_tcr = ad.AnnData()
# for i in filedir:
#     adata_tcr_tmp = ir.io.read_10x_vdj(i)
#     adata_tcr_tmp.var_names_make_unique()
#     adata_tcr = ad.concat([adata_tcr,adata_tcr_tmp],join='outer'
#                         , index_unique=None
#                         )
#     adata_tcr.obs_names_make_unique()
# adata_tcr.shape

# adata_tcr

# for i in filedir:
#     try:
#         adata_tcr_tmp = ir.io.read_10x_vdj(i)
#         print(f"读取成功: {i}")
#         print(f"X 是否存在: {adata_tcr_tmp.X is not None}")
#     except Exception as e:
#         print(f"读取失败: {i}, 错误: {e}")
