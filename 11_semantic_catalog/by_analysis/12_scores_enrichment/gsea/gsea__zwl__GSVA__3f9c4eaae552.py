import gseapy as gp
import scanpy as sc
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def do_gsva(adata=None,group_list=['leiden'],outpath = os.getcwd(),
            gmtfile=None,heatmap_num = None,
            cmap = "viridis"):
    """
    -adata: input adata
    -group_list: GSVA group use
    -outpath: output file address
    -gmtfile: gmt file input address
    -heatmap_num: heatmap gsva num
    """
    os.chdir(outpath)
    df = sc.get.obs_df(adata,keys=group_list+list(adata.raw.var_names),use_raw=True)
    df_group = df.groupby(group_list).mean()
    if len(group_list) > 1:
        df_group.index = df_group.index.map(lambda x: '_'.join(map(str, x)))
    df_group.T.to_csv(os.path.join(outpath,'gene_exp.xls'),sep="\t")
    df = df_group.T
    es = gp.gsva(data=df,
             gene_sets=gmtfile,
             min_size = 0,
             max_size = 2000,
             output = None)
    res = es.res2d
    res = res.pivot(index ='Term',columns= 'Name',values = 'ES' )
    if not os.path.exists("GSVA"):
        os.mkdir('GSVA')
    outpath1 = os.path.join(outpath,'GSVA')
    filename = gmtfile.split('/')[-1].split('.')[0]
    res.to_csv(os.path.join(outpath1,f'{filename}_GSVA.xls'),sep="\t")
    print(res.index)
    if heatmap_num == None:
        heatmap_num = len(res.index)
    res = pd.read_csv(os.path.join(outpath1,f'{filename}_GSVA.xls'),sep="\t",header=0,index_col=0)
    g = sns.clustermap(res.head(heatmap_num), annot=False, cmap=cmap, row_cluster=False, col_cluster=False)
    # sns.clustermap(res.T, annot=False, cmap='coolwarm', row_cluster=False, col_cluster=False,
    #                 figsize=(20,10),
    #                 cbar_pos=(0, .2, .03, .4)
    # )
    # Rotate the x-axis tick labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    save_fig(plt,out_dir=outpath1,file_name=filename,fig_size={'width':len(res.columns),'height':heatmap_num/10*10})  