import os
import pandas as pd
import scanpy as sc
import streamlit as st
import numpy as np
import gseapy as gp
import random
from matplotlib.pyplot import rc_context
#from gseapy.plot import barplot
import matplotlib.pyplot as plt
from moudle.common import image_view
from zipfile import ZipFile
from gseapy.plot import barplot, dotplot

def callback():
    st.session_state['tool'] = None

def enrich(tmp, species, out):
    if species == 'Human':
        gene_sets = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human']
        #gene_sets = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'KEGG_2021_Human']
    else:
        gene_sets = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse']

    if list(tmp['names']) != []:
        enr = gp.enrichr(gene_list=list(tmp['names']),
                         gene_sets=gene_sets,
                         organism=species,
                         cutoff=0.05,
                         no_plot=True,
                         outdir=None,
                         verbose=False
                        )
        print(enr.results)
        for i in set(enr.results["Gene_set"]):
            try:
                enr_out = enr.results.loc[enr.results["Gene_set"] == i]
                enr_out[enr_out['P-value'] < 0.05].to_csv(os.path.join(out, i + '_Enrichment.xls'), sep='\t', index=None)
                
                barplot(enr.results.loc[enr.results["Gene_set"] == i], title=i, cutoff=0.05, figsize=(6, 8), top_term=20, ofname=os.path.join(out, i + '_barplot.png'))
                barplot(enr.results.loc[enr.results["Gene_set"] == i], title=i, cutoff=0.05, figsize=(6, 8), top_term=20, ofname=os.path.join(out, i + '_barplot.pdf'))
                dotplot(enr.results.loc[enr.results["Gene_set"] == i], title=i, cutoff=0.05, figsize=(6, 8), top_term=20, ofname=os.path.join(out, i + '_dotplot.png'))
                dotplot(enr.results.loc[enr.results["Gene_set"] == i], title=i, cutoff=0.05, figsize=(6, 8), top_term=20, ofname=os.path.join(out, i + '_dotplot.pdf'))
            except:
                continue
def find_first_png_file(directory):
    for root, dirs, files in os.walk(directory):
        if "All" in dirs:
            all_folder = os.path.join(root, "All")
            for file in os.listdir(all_folder):
                if file.endswith(".png"):
                    return str(os.path.join(all_folder, file))   
    return None
def function_enrichment():
    select_list=[]
    adata=None
    if 'function_enrichment' not in st.session_state:
        st.session_state['function_enrichment'] = {}
        st.session_state['function_enrichment']['out_dir'] = os.path.join('tmp/tools', 'function_enrichment', ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8)))
    out_dir = st.session_state['function_enrichment']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 差异富集分析')
        st.write('对单细胞数据进行亚群注释后，我们通常会对不同组别或不同细胞类型间差异基因功能富集分析。这里我们使用GO和KEGG数据库对差异基因进行功能富集。')
        uploaded_file = st.file_uploader(label='h5ad文件上传:', type='h5ad')
    if uploaded_file:
        adata = sc.read_h5ad(uploaded_file)
        context.write(adata)
        columns = adata.obs.columns
        select_list = [element for element in columns if 'louvain' in element or 'leiden' in element]
        select_list.append('cell_type')
        select_list.append('All')

    params, response, note = st.columns([2, 5, 3])
    with params.form('functionEnrichment'):
            select =  st.selectbox('select:', select_list)
            comparison = st.selectbox('comparison:', ('group','sample','leiden','cell_type','cluster'))
            logfc = st.number_input(label='logfc:', min_value=0.0, max_value=2.0, value=0.5, step=0.1)
            fdr = st.number_input(label='fdr:', min_value=0.01, max_value=1.0, value=0.05, step=0.01)
            pct = st.number_input(label='pct:', min_value=0.01, max_value=1.0, value=0.25, step=0.01)
            species = st.selectbox('species:', ('Human','Mouse'))
            up_down = st.checkbox('上下调分开', value=False, help='')
            submitted = st.form_submit_button('🚀Submit', type='primary')

    if submitted:
        adata.uns['log1p'] = {'base': 2}
        if select!='All':
            for i in set(adata.obs[select]):
                adata_tmp=adata[adata.obs[select]==i]
                sc.tl.rank_genes_groups(adata_tmp, comparison, method='wilcoxon', pts=True)
                result = adata_tmp.uns['rank_genes_groups']
                groups = result['names'].dtype.names
                out_dir1=os.path.join(out_dir,i)
                os.makedirs(out_dir1, exist_ok=True)
                all_diff = []
                for c in groups:
                    out = os.path.join(out_dir1, c)
                    os.makedirs(out, exist_ok=True)
                    tmp = sc.get.rank_genes_groups_df(adata_tmp, group=c)  
                    tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals']<fdr) & (tmp['pct_nz_group'] >= pct)]
                    tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
                    tmp['cluster'] = c
                    all_diff.append(tmp)

                    try:
                        if list(tmp['names']) != []:
                            os.makedirs(os.path.join(out,'All'), exist_ok=True)
                            #st.write('step1 done')
                            enrich(tmp,species,os.path.join(out,'All'))
                            
                    except:
                        try:
                            if list(tmp['names']) != []:
                                enrich(tmp,species,os.path.join(out,'All'))

                        except:
                            continue                      
                    if up_down:
                        tmp_up = tmp[(tmp['logfoldchanges']>logfc)]
                        tmp_down =tmp[(tmp['logfoldchanges']<-logfc)]
                        try:
                            if list(tmp_up['names']) != []:
                                os.makedirs(os.path.join(out,'Up'), exist_ok=True)
                                enrich(tmp_up,species,os.path.join(out,'Up'))
                        except:
                            try:
                                if list(tmp_up['names']) != []:
                                    enrich(tmp_up,species,os.path.join(out,'Up'))
                            except:
                                continue 
                        try:
                            if list(tmp_down['names']) != []:
                                os.makedirs(os.path.join(out,'Down'), exist_ok=True)
                                enrich(tmp_down,species,os.path.join(out,'Down'))
                        except:
                            try:
                                if list(tmp_down['names']) != []:
                                    enrich(tmp_down,species,os.path.join(out,'Down'))
                            except:
                                continue
        else:
            adata_tmp=adata
            sc.tl.rank_genes_groups(adata_tmp, comparison, method='wilcoxon', pts=True)
            result = adata_tmp.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            out_dir1=os.path.join(out_dir,'All')
            os.makedirs(out_dir1, exist_ok=True)
            all_diff = []
            for c in groups:
                out = os.path.join(out_dir1, c)
                os.makedirs(out, exist_ok=True)
                tmp = sc.get.rank_genes_groups_df(adata_tmp, group=c)  
                tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals']<fdr) & (tmp['pct_nz_group'] >= pct) ]
                tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
                tmp['cluster'] = c
                all_diff.append(tmp)

                try:
                    if list(tmp['names']) != []:
                        os.makedirs(os.path.join(out,'All'), exist_ok=True)
                        #st.write('step1 done')
                        enrich(tmp,species,os.path.join(out,'All'))
                        
                except:
                    try:
                        if list(tmp['names']) != []:
                            enrich(tmp,species,os.path.join(out,'All'))

                    except:
                        continue                      
                if up_down:
                    tmp_up = tmp[(tmp['logfoldchanges']>logfc)]
                    tmp_down =tmp[(tmp['logfoldchanges']<-logfc)]
                    try:
                        if list(tmp_up['names']) != []:
                            os.makedirs(os.path.join(out,'Up'), exist_ok=True)
                            enrich(tmp_up,species,os.path.join(out,'Up'))
                    except:
                        try:
                            if list(tmp_up['names']) != []:
                                enrich(tmp_up,species,os.path.join(out,'Up'))
                        except:
                            continue 
                    try:
                        if list(tmp_down['names']) != []:
                            os.makedirs(os.path.join(out,'Down'), exist_ok=True)
                            enrich(tmp_down,species,os.path.join(out,'Down'))
                    except:
                        try:
                            if list(tmp_down['names']) != []:
                                enrich(tmp_down,species,os.path.join(out,'Down'))
                        except:
                            continue
            for group in groups:
            # 提取 log2 fold change 和 p 值
                log2fc = result['logfoldchanges'][group]
                pvals = result['pvals'][group]
                pvals_adj = result['pvals_adj'][group]  # 如果有调整后的 p 值

                # 创建一个DataFrame
                deg_df = pd.DataFrame({
                    'gene': result['names'][group],
                    'log2FC': log2fc,
                    'pval': pvals,
                    'pval_adj': pvals_adj
                })

                # 根据 -log10(p-value) 排序
                deg_df['neg_log10_pval'] = -np.log10(deg_df['pval'])


                # 假设 deg_df 已经从差异表达结果中创建

                # 添加颜色映射的条件
                deg_df['color'] = '#bdbdbd'  # 默认灰色
                deg_df.loc[deg_df['log2FC'] > 1, 'color'] = '#da9287'  # log2FC > 1, 红色
                deg_df.loc[deg_df['log2FC'] < -1, 'color'] = '#55a1cd'  # log2FC < -1, 蓝色

                # 根据 log2FC 绝对值选择前5个高表达的基因
                top_genes = deg_df.loc[deg_df['log2FC'].abs().nlargest(5).index]

                # 设置绘图样式
                sns.set(style='white')  # 去掉网格

                # 创建火山图
                plt.figure(figsize=(8, 6))

                # 绘制散点图：去掉边框（edgecolor='none'），只保留填充色
                scatter = plt.scatter(deg_df['log2FC'], deg_df['neg_log10_pval'], 
                                    c=deg_df['color'], edgecolor='none', alpha=0.7, s=50)

                # 标注前5个高表达基因的名称
                # for i, row in top_genes.iterrows():
                #     plt.text(row['log2FC'], row['neg_log10_pval'], row['gene'], 
                #              color='black', fontsize=10, ha='center', va='center')

                # 设置标题和标签
                plt.title(group, fontsize=16)
                plt.xlabel('Log2 Fold Change', fontsize=12)
                plt.ylabel('-Log10 P-value', fontsize=12)

                # 绘制阈值线
                plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--')  # p-value < 0.05
                plt.axvline(x=1, color='gray', linestyle='--')  # Log2FC > 1
                plt.axvline(x=-1, color='gray', linestyle='--')  # Log2FC < -1

                # 去掉右侧和上方的框线，保留下框线和左框线
                sns.despine(left=False, bottom=False, right=True, top=True)

                # 创建图例
                legend_elements = [
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#da9287', markersize=10, label='Upregulated'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#55a1cd', markersize=10, label='Downregulated'),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#bdbdbd', markersize=10, label='No Change')
                ]

                # 添加图例
                plt.legend(handles=legend_elements, loc='upper left', fontsize=12)
                plt.savefig(os.path.join(out, f'volcano_{group}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
                plt.savefig(os.path.join(out, f'volcano_{group}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    png_file = find_first_png_file(out_dir)
    if png_file is not None:
        if os.path.exists(png_file):
                with response.expander("**:blue[结果解读📜]** ", expanded=True):
                    tabs = st.tabs(["🍓 Function enrichment"])        
                    image_view(tabs[0], png_file=png_file)
                with note.expander("**:blue[参数说明📜]** ", expanded=True):
                    st.markdown(f"""
                        - **select:** 需要分别进行富集分析的列；
                        - **comparison:** 比较组；
                        - **logfc:** 差异基因log Foldchange卡值，默认为0.5；
                        - **fdr:** 差异基因plavue 卡值，默认为0.05；
                        - **pct:** 基因表达细胞占比，默认0.25；
                        - **species:** 物种，可选小鼠和人；
                        - **上下调分开:** 是否需要上下调基因分别进行功能富集分析，默认合并上下调基因进行功能富集分析；       
                        - **结果文件位置：:rainbow[{st.session_state['function_enrichment']['out_dir']}]**
                        """)
                zip_filename = os.path.join(out_dir,"function_enrichment.zip")
                with ZipFile(zip_filename, "w") as zipf:
                    for root, dirs, files in os.walk(out_dir):
                        for file in files:
                            if not file.endswith(".zip"):
                                zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), out_dir))
                with open(zip_filename, "rb") as f:
                    bytes_data = f.read()    
                btn = note.download_button(
                                                label="Download",
                                                data=bytes_data,
                                                file_name='function_enrichment.zip',
                                                mime="application/zip",
                                                help='Download results.'
                                            )
    else:
        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            st.markdown("""👈**差异基因功能富集分析，选择分组(样本)信息，设置分析参数。**""") 
        with note.expander("**:blue[参数说明📜]** ", expanded=True):
            st.markdown(f"""
                - **select:** 需要分别进行富集分析的列；
                - **comparison:** 比较组；
                - **logfc:** 差异基因log Foldchange卡值，默认为0.5；
                - **fdr:** 差异基因plavue 卡值，默认为0.05；
                - **pct:** 基因表达细胞占比，默认0.25；
                - **species:** 物种，可选小鼠和人；
                - **上下调分开:** 是否需要上下调基因分别进行功能富集分析，默认合并上下调基因进行功能富集分析；      
                - **结果文件位置：:rainbow[{st.session_state['function_enrichment']['out_dir']}]**
                """)
