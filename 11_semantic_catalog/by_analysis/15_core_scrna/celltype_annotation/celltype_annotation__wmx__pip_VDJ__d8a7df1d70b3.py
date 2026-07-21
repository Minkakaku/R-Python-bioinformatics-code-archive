import sys
# sys.path.append("/mnt/NexData_file/share/scRNA-tools/VDJ/")
from sciana import *
import scanpy as sc
import os
import anndata as ad
import scirpy as ir
import muon as mu
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from mudata import MuData
import re
import mudata as md
import argparse
import shutil

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


def clo_modu(vdjobj):
    mdata = vdjobj.dat
    outdir = os.path.join(vdjobj.output, '7.Clonotype modularity')
    
    # 1. 计算 GEX 的邻居关系
    sc.pp.neighbors(mdata["gex"])
    
    # --- 核心修复逻辑 ---
    # 统计每个克隆的细胞数量
    clonotype_counts = mdata['airr'].obs['cc_aa_alignment'].value_counts()
    
    # 找出细胞数 >= 2 的克隆名称
    valid_clonotypes = clonotype_counts[clonotype_counts >= 2].index.tolist()
    
    # 检查是否有符合条件的克隆
    if len(valid_clonotypes) < 1:
        print("Warning: No clonotypes with at least 2 cells found. Skipping modularity.")
        return vdjobj

    # 重要步骤：对 mdata 进行过滤，只保留细胞数 >= 2 的克隆对应的细胞
    # 否则 ir.tl.clonotype_modularity 会尝试对所有克隆（包括 1 个细胞的）做差异分析
    mask = mdata['airr'].obs['cc_aa_alignment'].isin(valid_clonotypes)
    mdata_sub = mdata[mask, :].copy() 
    # --- 修复逻辑结束 ---

    # 判断样本数量是否足够画图
    if(len(mdata_sub['gex'].obs['sample'].drop_duplicates().to_list()) >= 1):
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        
        try:
            # 在过滤后的子集上运行计算
            ir.tl.clonotype_modularity(mdata_sub, target_col="cc_aa_alignment")
            
            # 绘图逻辑
            ir.pl.clonotype_network(mdata_sub, color="airr:clonotype_modularity", 
                                    label_fontsize=9, panel_size=(6, 6), base_size=5, ax=axs[0])
            mu.pl.embedding(mdata_sub, basis="gex:umap", color="airr:clonotype_modularity", ax=axs[1])
            
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                
            save_fig(fig, out_dir=outdir, file_name='7.1.Clonotype_modularity_score', fig_size={'width':10, 'height':5})
            
            # 将计算结果写回原始对象（对应细胞的位置）
            # 这一步是为了保证下游分析能用到这个 Score
            mdata.obs.loc[mask, "airr:clonotype_modularity"] = mdata_sub.obs["airr:clonotype_modularity"]
            
        except Exception as e:
            print(f"Skipping modularity due to error: {e}")
            
        vdjobj.dat = mdata
        return vdjobj

def H5ADHandler(adata, groupby = 'cell_type', subsample = 'all', subcelltype = 'all', renamesample = None, renamecelltype = None):
    """
    :return: 处理后的AnnData对象
    """
    adata = sc.read_h5ad(adata)
    if adata.raw is not None:
        adata = adata.raw.to_adata()
    else:
        adata = adata
    adata.obs[groupby] = adata.obs[groupby].astype(str).astype(object).astype("category")
    if renamecelltype is not None:
        celltype_dict = {}
        with open(renamecelltype, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    key, celltypes = parts
                    celltype_list = [celltype.strip() for celltype in celltypes.split(',') if celltype.strip()]
                    celltype_dict[key] = celltype_list
        adata.obs[groupby] = adata.obs[groupby].replace(celltype_dict)
        print(f"修改细胞类型: {celltype_dict}")
    if renamesample is not None:
        sample_dict = {}
        with open(renamesample, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    key, samples = parts
                    sample_list = [sample.strip() for sample in samples.split(',') if sample.strip()]
                    sample_dict[key] = sample_list
        adata.obs['sample'] = adata.obs['sample'].replace(sample_dict)
        print(f"修改样本: {sample_dict}")
    if not subcelltype == 'all':
        subcelltype = subcelltype.split(',')
        print(f"细胞类型取子集: {subcelltype}")
        adata=adata[adata.obs[groupby].isin(subcelltype), :]
    if not subsample == 'all':
        subsample = subsample.split(',')
        print(f"样本取子集: {subsample}")
        adata=adata[adata.obs['sample'].isin(subsample), :]
    # adata.uns['log1p']['base'] = 2
		# 修改barcodes + 样本名便于一一对应匹配
    obs_name_list = adata.obs_names.astype(str).str[:18]+("_")+adata.obs['sample'].astype(str)
    adata.obs_names = obs_name_list
    return adata

def VDJHandler(file_path,method):
    """
    :return: 处理后的AnnData对象
    """
    filepath = file_path
    filedir = [
        os.path.join(filepath, i, 'filtered_contig_annotations.csv')
        for i in os.listdir(filepath)
    ]
    # 先读入第一个文件
    adata_tcr = ir.io.read_10x_vdj(filedir[0])
    obs_name_list = adata_tcr.obs_names.to_list()
    if(method=='BCR'):
        sample = re.search(r'vdj_b/(\w+)/', filedir[0]).group(1)
    else:
        sample = re.search(r'vdj_t/(\w+)/', filedir[0]).group(1)
    obs_name_list = [i+("_")+sample for i in obs_name_list]
    adata_tcr.obs_names = obs_name_list
    adata_tcr.var_names_make_unique()
    # 依次合并后续文件
    for path in filedir[1:]:
        adata_tcr_tmp = ir.io.read_10x_vdj(path)
        if(method=='BCR'):
            sample = re.search(r'vdj_b/(\w+)/', path).group(1)
        else:
            sample = re.search(r'vdj_t/(\w+)/', path).group(1)
        obs_name_list = adata_tcr_tmp.obs_names.to_list()
        obs_name_list = [i+("_")+sample for i in obs_name_list]
        adata_tcr_tmp.obs_names = obs_name_list
        adata_tcr_tmp.var_names_make_unique()
        adata_tcr = ad.concat([adata_tcr, adata_tcr_tmp], join='outer', index_unique=None)
    # 保证 obs 唯一
    adata_tcr.obs_names_make_unique()
    return adata_tcr

def main(data_path, file_path, out_dir, species, method, cutoff, groupby, sub_sample, sub_celltype, rename_sample, rename_celltype):
    """
    主函数：处理输入参数，调用类方法生成图表
    :param data_path: 输入目录，h5ad文件
    :param file_path: VDJ注释表文件路径
    :param out_dir: 输出目录
    :param species: 物种（"小鼠" 或 "人"）
    :param method: 分析类型（"BCR" 或 "TCR"）
    :param cutoff: 默认阈值15
    :param groupby: 展示分类对象，默认'cell_type'
    :param sub_sample: 样本取子集
    :param sub_subcelltype: 细胞类型取子集
    :param rename_sample: 样本重命名
    :param rename_celltype: 细胞类型重命名
    """
    # 创建输出目录
    os.makedirs(out_dir, exist_ok=True)
    # 输入目录中的h5ad文件
    print(f"Processing {data_path}")
    h5ad_handler = H5ADHandler(data_path, groupby=groupby, subsample=sub_sample, subcelltype=sub_celltype, renamesample=rename_sample, renamecelltype=rename_celltype)
    with rc_context({'figure.figsize': (12, 8), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight'}):
         sc.pl.umap(h5ad_handler,color=[groupby,'cell_type','sample','group'],show=False)
         plt.savefig(os.path.join(out_dir, 'dimplot_umap.png'))
         plt.savefig(os.path.join(out_dir, 'dimplot_umap.pdf'))
    h5ad_handler.shape
    VDJ_handler = VDJHandler(file_path,method)
    VDJ_handler.shape
    mdata = MuData({'gex':h5ad_handler,'airr':VDJ_handler})
    mu.pp.intersect_obs(mdata)
    os.environ["HDF5_OBJECT_HEADER_MAX_SIZE"] = "1000000"
    vdj = vdj_obj(mdata = mdata,species= species,method = method,cutoff = cutoff,celltype = groupby,output = out_dir)
    vdj = vdj.qc()
    vdj = vdj.def_cc()
    # vdj = def_cc(vdj)
    vdj = vdj.clo_exp()
    vdj = vdj.clo_abu()
    vdj = vdj.clo_siz()
    vdj = vdj.com_rep()
    vdj = vdj.clo_modu()
    # vdj = clo_modu(vdj)
    vdj = vdj.vdj_pre()
    vdj = vdj.mot_ana()
    vdj = vdj.sha_clo()
    
    # os.mkdir(f'{out_dir}0.Cellranger')
    outdir = os.path.join(out_dir,'0.Cellranger')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    ir.io.write_airr(VDJ_handler,os.path.join(f'{out_dir}/0.Cellranger/filtered_contig_annotations.tsv'))
    shutil.copy("/mnt/NexData_file/share/scRNA-tools/VDJ/report.rmd", out_dir)
    os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(out_dir))
    vdj.save_h5mu()
    # del mdata.obs['VJ_1_d_call']
    # del mdata.obs['VJ_2_d_call']
    # testdata = md.read("/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR/Result.h5mu")



if __name__ == "__main__":
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="Generate UMAP and Violin plots for single-cell data.")
    parser.add_argument('--data_path', type=str, required=True, help='Input h5ad Path')
    parser.add_argument('--file_path', type=str, required=True, help='Path to the VDJ annotations files (csv format).')
    parser.add_argument('--out_dir', type=str, required=True, help='Output directory containing results files.')
    parser.add_argument('--species', type=str, required=True, choices=['mouse', 'human'], help='Species: mouse or human.')
    parser.add_argument('--method', type=str, default='TCR', required=True, choices=['TCR', 'BCR'], help='TCR or BCR')
    parser.add_argument('--cutoff', type=str, default=15, required=False, help='一致性评分阈值')
    parser.add_argument('--groupby', type=str, default='cell_type', required=False, help='Output directory for saving plots.')
    parser.add_argument('--sub_sample', type=str, default='all', required=False, help='sub samples')
    parser.add_argument('--sub_celltype', type=str, default='all', required=False, help='sub cell types')
    parser.add_argument('--rename_sample', type=str, default=None, required=False, help='rename samples')
    parser.add_argument('--rename_celltype', type=str, default=None, required=False, help='rename cell types')
    args = parser.parse_args()

    # 调用主函数
    main(args.data_path, args.file_path, args.out_dir, args.species, args.method, args.cutoff, args.groupby, args.sub_sample, args.sub_celltype, args.rename_sample, args.rename_celltype)

# output = '/mnt/NexData_file/home/sh/Projects/FWAE20252674_BCR/BCR/'
# pid = "FWAE20252674"
# ptitle = "BCR"
# paddress = "浙江大学医学院附属第一医院"
# pcustmor = "姚秋妹"
# sample_num = "13"
# rmp = ''
# with open("/mnt/NexData_file/share/scRNA-tools/VDJ/report.rmd",'r',encoding='utf-8') as f:
#     line = f.read()
#     rmd = line.replace('%pid%', pid).replace('%ptitle%', ptitle).replace('%paddress%', paddress).replace('%pcustmor%', pcustmor).replace('%num_sample%', str(sample_num))
# with open(os.path.join(output,'report.rmd'),'w',encoding = 'utf-8') as w:
#     w.write(rmd)
# f.close()
# w.close()
#os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(output))
