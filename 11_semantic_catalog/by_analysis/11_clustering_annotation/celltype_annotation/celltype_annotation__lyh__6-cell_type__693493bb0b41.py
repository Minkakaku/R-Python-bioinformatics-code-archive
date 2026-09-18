import os
import re
import json
import shutil
import pandas as pd
import scanpy as sc
import streamlit as st
import streamlit.components.v1 as components
from utils.utils import image_view
from utils.function import celltype_annotation, visua_gene_express
rmd_template = 'config/template/report.rmd'


st.subheader("🍊 :rainbow[Celltype Identification]")
st.info('''&emsp;&emsp;通过差异基因分析以及Marker 基因鉴定，我们可以确定每个细胞群中的真实细胞类型。目前，细胞类型鉴定主要有两种方法，**一种是基于已有细胞类型的Marker 基因**，通过软件或者人为鉴定这些特定Marker 所对应的细胞类型，即主要有监督和半监督这两种方式；
        **另一种是基于现有表达谱参考数据集**，基于相似性通过非监督的方式鉴定细胞类型。这两种方法各有优劣。前者人为干扰因素较多且可能对现有Marker 过于依赖。后者基于参考数据集进行鉴定，可有效避免人为干扰，但是也受制于参考数据集中的现有数据，只能鉴定到大类。
        目前的经验方法是，细胞注释要同时结合自动化注释和人工注释，这两种方法并不是独立割裂开来的，需要搭配使用，软件自动化细胞注释方法方便且系统化，完全依赖于参考数据集，注释的结果有时也并不是高置信度注释。结合两种方法进行鉴定细胞类型，使细胞类群鉴定结果更可靠。
        ''')

out_dir = os.path.join(st.session_state['pipe']['out_dir'], '06-Cell_type_annotation')      
params, response, note = st.columns([2, 5, 3])
with params.form('annotation'):
    index = 0
    if st.session_state['pipe']['organism'] == 'mouse':
        index = 1
    st.session_state['pipe']['organism'] = st.selectbox(label = 'Organism:', options = ('human', 'mouse'), index=index, help='选择样本物种信息，只支持人和小鼠。', key='opt_organism')
    tissue = st.selectbox(label = 'Tissue type:', options = ['All', 'bone', 'breaste', 'lung', 'pbmc'], help='主要用于后续细胞类型鉴定数据库比对，默认使用所有组织类型。', key='opt_tissue')        
    st.session_state['pipe']['database'] = st.selectbox(label = 'Database:', options = ('PanglaoDB', 'CellTypist'), help='选择细胞类型自动注释数据库', key='opt_database')               
    submitted = st.form_submit_button("🚀Annot", type='primary')
    if submitted:
        if st.session_state['pipe']['adata_rank_genes_groups']:
            with response:
                with st.spinner('Calculating...'):
                    st.session_state['pipe']['adata_cell_type'], st.session_state['pipe']['annotation_dict'] = celltype_annotation(st.session_state['pipe']['adata_rank_genes_groups'],
                                                                                                                    st.session_state['pipe']['database'],
                                                                                                                    st.session_state['pipe']['organism'],
                                                                                                                    out_dir)
                    st.session_state['pipe']['adata_cell_type'].write(os.path.join(st.session_state['pipe']['out_dir'], 'result.h5ad'))   
        else:
            response.error('请检查上一步是否正确执行！', icon="🚨")

        with open(os.path.join(st.session_state['pipe']['out_dir'], 'parameters.json'), 'w') as w:
            for k,v in st.session_state['pipe'].items():
                if 'adata' not in k and 'btn' not in k:
                    w.write(k+'\t'+str(v))
        #     w.write(json.dumps({k:v for k,v in st.session_state['pipe'].items() if 'adata' not in k or 'btn' not in k}))    

if st.session_state['pipe']['adata_cell_type']:
    with response.expander("**:blue[结果解读📜]** ", expanded=True):
        tabs = st.tabs(["🍊 UMAP-all", "🍊 UMAP-group", "🍊 TSNE-all", "🍊 Cell type proportion of samples"])
        image_view(tabs[0], png_file=os.path.join(out_dir, 'umap_cluster_celltype.png'), caption='样本整体UMAP降维聚类结果')

        if os.path.exists(os.path.join(out_dir, 'umap_split_group_celltype.png')):
            image_view(tabs[1], png_file=os.path.join(out_dir, 'umap_split_group_celltype.png'), caption='样本分组UMAP降维聚类结果')
        else:
            pass
        
        image_view(tabs[2], png_file=os.path.join(out_dir, 'tsne_cluster_celltype.png'), caption='样本整体TSNE降维聚类结果')
        image_view(tabs[3], png_file=os.path.join(out_dir, 'celltype_proportion_sample.png'), caption='各样本中细胞类型占比')

    with note.form('annotation_rename'):
        df = pd.DataFrame(st.session_state['pipe']['annotation_dict'].items(), columns=['Cluster', 'annotation'])
        df['cell_type'] = df['annotation']
        renames = st.data_editor(df, hide_index=True, disabled=('Cluster','annotation'), use_container_width=True).set_index('Cluster')['cell_type'].to_dict() 
        submitted = st.form_submit_button("🚀Rname", type='primary')
        if submitted:
            with response:
                with st.spinner('Calculating...'):
                    st.session_state['pipe']['adata_cell_type'], st.session_state['pipe']['annotation_dict'] = celltype_annotation(st.session_state['pipe']['adata_rank_genes_groups'],
                                                                                                                    st.session_state['pipe']['database'],
                                                                                                                    st.session_state['pipe']['organism'],
                                                                                                                    out_dir,
                                                                                                                    renames)
                    st.session_state['pipe']['adata_cell_type'].write(os.path.join(st.session_state['pipe']['out_dir'], 'result.h5ad'))
else:
    with response.expander("**:blue[结果解读📜]** ", expanded=True):
        st.markdown("""**分析结果文件位置：06-Cell_type_annotation/**""")   
    with note.expander("**:blue[参数说明📜]** ", expanded=True):
        st.markdown("""
                    - **Organism**: 样本物种信息;
                    - **Tissue type**: 样本组织类型;
                    - **Database**: 细胞类型自动注释参考数据库;
                    """)

if st.session_state['pipe']['adata_cell_type']:
    st.subheader('🍉 基因表达展示')
    st.info('&emsp;&emsp;整合PanglaoDB, CellMarker等数据库, 对单细胞测序细胞大类鉴定常用Marker基因进行汇总。细胞Cluster使用自动化注释后，可以将以下常用的细胞鉴定Marker的表达情况在UMAP/TSNE图上进行展示， 从而对自动注释结果进行校正。')
    params, response, note = st.columns([2, 5, 3])
    types = params.selectbox(label = '选择基因展示/预定义细胞类型Markers：', options = ('选择基因（支持多基因）', '使用预定义细胞类型Marker基因'), help='自定义基因查找，或使用预先定义的细胞类型Markers。')
    with params.form('gene_exp'):
        select_types = []
        if types == '选择基因（支持多基因）':
            genes = st.multiselect(label='基因:', options=st.session_state['pipe']['adata_filtered'].var_names)       
        else:
            genes = []
            celltyps = (st.session_state['pipe']['marker_dict'][st.session_state['pipe']['organism']]).keys()
            select_types = st.multiselect(label='细胞类型:', options=celltyps)       
            if select_types:
                for t in select_types:
                    genes.extend([i for i in st.session_state['pipe']['marker_dict'][st.session_state['pipe']['organism']][t] if i in st.session_state['pipe']['adata_filtered'].var_names])

        submitted = st.form_submit_button("🚀Plot", type='primary')
        
        if select_types:
            st.write({i:st.session_state['pipe']['marker_dict'][st.session_state['pipe']['organism']][i] for i in select_types})
        else:
            st.write(genes)

        with response.expander(f"**:red[基因表达]** ", expanded=True):
            if submitted:
                if not genes:
                    st.error('Please select genes to visualization of gene expression.', icon="🔥")
                else:
                    with st.spinner('Calculating...'):
                        st.session_state['pipe']['gene_expression'], fig = visua_gene_express(st.session_state['pipe']['adata_cell_type'], genes, out_dir)

            if st.session_state['pipe']['gene_expression']:
                st.image(st.session_state['pipe']['gene_expression'], caption=f"{genes}基因表达", use_column_width=True)

    with note.expander(f"**:blue[参数说明📜]** ", expanded=True):
        st.markdown("""- **选择基因(支持多基因)，在UMAP上展示选中基因的表达；或者使用预定义细胞类型Markers，在UMAP上展示这些Marker基因的表达。**""")


if st.session_state['pipe']['adata_cell_type']:
    st.subheader(f"🍑 生成报告  {os.path.abspath(st.session_state['pipe']['out_dir'])}")
    with st.form('report'): 
        cols = st.columns([3,5,4,3], gap='large')
        pid = cols[0].text_input('项目编号:', st.session_state['pipe']['pid'])
        ptitle = cols[1].text_input('项目名称:', '')
        paddress = cols[2].text_input('客户单位:', '')
        pcustmor = cols[3].text_input('客户姓名:', '')               
        submitted = st.form_submit_button("🚀Report", type='primary')

    out_dir = os.path.abspath(st.session_state['pipe']["out_dir"])
    
    if submitted:
        rmd = ''
        with open(rmd_template, 'r', encoding='utf-8') as f:
            line = f.read()
            if not os.path.exists(os.path.join(out_dir, '01-CellRanger', 'metrics_summary.csv')):
                line = re.sub(r'## 测序数据统计.*? 所有细胞的基因总数。', '', line, flags=re.DOTALL)
            sample_num = len(st.session_state['pipe']['adata_raw'].obs['sample'].drop_duplicates().tolist())
            cell_num_raw = st.session_state['pipe']['adata_raw'].X.shape[0]
            cell_num_after_filter = str(st.session_state['pipe']['adata_filtered'].X.shape[0])
            cell_types = ', '.join(st.session_state['pipe']['adata_cell_type'].obs['cell_type'].drop_duplicates().tolist())                   
            cluster_num = len(st.session_state['pipe']['adata_cell_type'].obs[f"{st.session_state['pipe']['cluster_method']}_{st.session_state['pipe']['resolution']}"].drop_duplicates().tolist())
            rmd = line.replace('%pid%', pid).replace('%ptitle%', ptitle).replace('%paddress%', paddress).replace('%pcustmor%', pcustmor).replace('%num_sample%', str(sample_num)).replace('%organism%', st.session_state['pipe']['organism'])
            rmd = rmd.replace('%cell_num_raw%', str(cell_num_raw)).replace('%cell_num_after_filter%', str(cell_num_after_filter)).replace('%cluster_num%', str(cluster_num)).replace('%cell_types%', cell_types) 
        with open(os.path.join(out_dir, 'report.rmd'), 'w', encoding='utf-8') as w:
            w.write(rmd)
        with st.spinner('Report generation ...'):   
            os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(out_dir))
        with st.spinner('File compression ...'):
            shutil.make_archive(f"tmp/{st.session_state['pipe']['pid']}", 'zip', out_dir)

    if os.path.exists(os.path.join(out_dir, 'report.html')): 
        # zip_report(out_dir, os.path.join(pid + '.zip'))
        cols = st.columns(8)
        if os.path.exists(f"tmp/{st.session_state['pipe']['pid']}.zip"):
            cols[7].download_button(
                                    label="📁 Download",
                                    data=open(f"tmp/{st.session_state['pipe']['pid']}.zip", "rb"),
                                    file_name=f"{st.session_state['pipe']['pid']}.zip",
                                    mime="application/octet-stream",
                                   )
        report = open(os.path.join(out_dir, 'report.html'))
        components.html(report.read(), height=1600, scrolling=True)


