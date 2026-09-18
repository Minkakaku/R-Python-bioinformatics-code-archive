import os
import random
import datetime
import pandas as pd
import scanpy as sc
import streamlit as st
from utils.utils import *
from streamlit_pills import pills
import plotly as py
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from utils.subcluster import make_report_dir, sub_cluster, image_view, visua_gene_express, celltype_rename


def form_recluster(adata):
    out_dir = st.session_state['subc']['out_dir']
    params, response, note = st.columns([2, 6, 2])
    with params.form('sub-cluster'):
        st.session_state['subc']['cell_type'] = st.multiselect(label = 'cell cype:', options = sorted(adata.obs['cell_type'].drop_duplicates().to_list()), placeholder='选择细胞类型')
        st.session_state['subc']['n_pcs'] = st.number_input(label='n_pcs:', min_value=5, max_value=50, value=st.session_state['subc']['n_pcs'], step=1)
        st.session_state['subc']['n_neighbors'] = st.number_input(label='n_neighbors:', min_value=10, max_value=50, value=st.session_state['subc']['n_neighbors'], step=1)
        st.session_state['subc']['resolution'] = st.number_input(label='resolution:', min_value=0.1, max_value=1.0, value=st.session_state['subc']['resolution'], step=0.1)
        st.session_state['subc']['dot_size'] = st.number_input(label='dot_size:', min_value=1, max_value=200, value=st.session_state['subc']['dot_size'], step=1, help='调整UMAP图中点的大小。')
        st.session_state['subc']['organism'] = st.radio("**functional enrichment:**", ["No", "human", "mouse"], captions = ["Do not Functional enrichment.", "Use human database.", "Use mouse database."], help='是否对各cluster Marker基因做功能富集分析。')
        submitted = st.form_submit_button("🚀提交", type='primary')
    if submitted:
        if not st.session_state['subc']['cell_type']:
            response.error('请至少选择1个细胞类型！', icon="🚨")
        else:
            sub_adata = adata[adata.obs['cell_type'].isin(st.session_state['subc']['cell_type']) , :]
            with st.spinner('数据加载中...'):
                st.session_state['subc']['sub_adata'] = sub_cluster(sub_adata, 
                                                                    st.session_state['subc']['n_pcs'], 
                                                                    st.session_state['subc']['n_neighbors'], 
                                                                    st.session_state['subc']['resolution'],
                                                                    st.session_state['subc']['organism'],
                                                                    st.session_state['subc']['dot_size'],
                                                                    st.session_state['subc']['out_dir'])
    if st.session_state['subc']['sub_adata']:
        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            tabs = st.tabs(["🍓 UMAP-all", "🍓 UMAP-group", "🍓 Markers", "🍓 Gene list", "🍓 Counts"])        
            tcols = tabs[0].columns([1,2])
            tcols[0].subheader(f"Cells Number: {st.session_state['subc']['sub_adata'].obs.shape[0]}")
            image_view(tcols[1], png_file=os.path.join(out_dir, 'dimplot_umap_cluster.png'), resize=(400,400))
            image_view(tabs[1], png_file=os.path.join(out_dir, 'umap_group_cluster.png'))
            image_view(tabs[2], png_file=os.path.join(out_dir, 'dotplot_marker.png'))           
            marker_gene = pd.read_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index_col='cluster')
            cluster = tabs[3].selectbox(label = 'Select cluster view marker gene:', options = sorted(marker_gene.index.drop_duplicates().to_list()))
            if not cluster:
                tabs[3].dataframe(marker_gene, use_container_width=True)
            else:    
                df_marker = marker_gene.loc[cluster]
                tabs[3].dataframe(df_marker, use_container_width=True)
            image_view(tabs[4], png_file=os.path.join(out_dir, 'umap_n_genes_by_counts.png'))
        with note.expander("**:blue[参数说明📜]** ", expanded=True):
            st.markdown(f"""
                - **Organism**: 选择一个或多个细胞类型，继续亚群细分;
                - **n_pcs**: 用于构建邻域图的主成分数量。默认使用前10个主成分；
                - **n_neighbors**: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved;
                - **resolution**: 亚群细分细胞聚类分辨率(默认值为0.3);
                - **orgnism**: 选择物种信息(用于功能富集分析);
                - **结果文件位置：:rainbow[{st.session_state['subc']['out_dir']}]**
                """)
    else:
        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            st.markdown("""👈**选择细胞类型，修改分析参数，点击提交。**""")
        with note.expander("**:blue[参数说明📜]** ", expanded=True):
            st.markdown("""
                        - **Organism**: 选择一个或多个细胞类型，继续亚群细分;
                        - **n_pcs**: 用于构建邻域图的主成分数量。默认使用前10个主成分；
                        - **n_neighbors**: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved;
                        - **resolution**: 亚群细分细胞聚类分辨率(默认值为0.3);
                        - **orgnism**: 选择物种信息(用于功能富集分析);
                        """)

def form_gene_expression(adata):
    st.subheader("② Gene expression")
    out_dir = st.session_state['subc']['out_dir']
    params, response, note = st.columns([2, 6, 2])
    with params.form('gene_exp'):
        genes = st.multiselect(label='visualization of gene expression:', options=adata.raw.to_adata().var_names)
        dot_size = st.number_input(label='dot_size:', min_value=1, max_value=200, value=25, step=1, help='调整UMAP图中点的大小。')
        submitted_view = st.form_submit_button("🚀Plot", type='primary')
        if submitted_view:
            if not genes:
                response.error('Please select genes to visualization of gene expression.', icon="🔥")
            else:
                st.session_state['subc']['gene_expression'] = visua_gene_express(st.session_state['subc']['sub_adata'], 
                                                                                 genes,
                                                                                 dot_size, 
                                                                                 out_dir)
    if st.session_state['subc']['gene_expression']:
        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            tabs = st.tabs(["🍓 Color map"])
            if os.path.exists(os.path.join(out_dir, '_'.join(genes)+'.png')):
                image_view(tabs[0], png_file=os.path.join(out_dir, '_'.join(genes)+'.png'))
                # image_view(tabs[1], png_file=os.path.join(out_dir, '_'.join(genes)+'_dotplot.png'))
        with note.expander("**:blue[参数说明📜]** ", expanded=True):
            st.markdown(f"""
                        - **选择基因(可多选)，展示基因表达；**
                        - **结果文件位置：:rainbow[{st.session_state['subc']['out_dir']}]**
                        """)
    else:
        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            st.markdown("""👈**选择基因(可多选)，展示基因表达**""")
        with note.expander("**:blue[参数说明📜]** ", expanded=True):
            st.markdown(f"""
                        - **选择基因(可多选)，展示基因表达；**
                        - **结果文件位置：:rainbow[{st.session_state['subc']['out_dir']}]**
                        """)

def form_rename():
    st.subheader("③ Cell type ident")
    out_dir = st.session_state['subc']['out_dir']
    params, response, note = st.columns([2, 6, 2])    
    with params.form('subcluster_rename'):
        tmp = sorted(st.session_state['subc']['sub_adata'].obs['leiden'].drop_duplicates().to_list())
        tmp_dict = dict(zip(tmp, tmp))
        df = pd.DataFrame(tmp_dict.items(), columns=['Cluster', 'cell_type'])
        renames = st.data_editor(df, hide_index=True, disabled=('Cluster'), use_container_width=True).set_index('Cluster')['cell_type'].to_dict() 
        submitted = st.form_submit_button("🚀Rname", type='primary')
        if submitted:
            with st.spinner('Calculating...'):
                st.session_state['subc']['sub_adata_celltype'] = celltype_rename(st.session_state['subc']['sub_adata'], renames, out_dir)
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                tabs = st.tabs(["🍓 UMAP-all", "🍓 UMAP-group", "🍓 Cell type proportion of samples"])       
                image_view(tabs[0], png_file=os.path.join(out_dir, 'dimplot_umap_celltype.png'))
                image_view(tabs[1], png_file=os.path.join(out_dir, 'umap_groups_celltype.png'))
                image_view(tabs[2], png_file=os.path.join(out_dir, 'celltype_proportion_sample.png'))
            with note.expander("**:blue[参数说明📜]** ", expanded=True):
                st.markdown(f"""
                            - **修改cell_type列信息，自定义细胞类型名称。**
                            - **结果文件位置：:rainbow[{st.session_state['subc']['out_dir']}]**
                            """)
        else:
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                st.markdown("""👈**修改cell_type列信息，自定义细胞类型名称。**""")   
            with note.expander("**:blue[参数说明📜]** ", expanded=True):
                st.markdown(f"""
                            - **修改cell_type列信息，自定义细胞类型名称。**
                            - **结果文件位置：:rainbow[{st.session_state['subc']['out_dir']}]**
                            """)

def subcluster():
    if 'subc' not in st.session_state:
        st.session_state['subc'] = {}
        st.session_state['subc']['out_dir'] = make_report_dir(os.path.join('tmp/tools/subcluster', ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))))
        st.session_state['subc']['sub_adata'] = None
        st.session_state['subc']['cell_type'] = None
        st.session_state['subc']['n_pcs'] = 10
        st.session_state['subc']['n_neighbors'] = 20
        st.session_state['subc']['resolution'] = 0.5
        st.session_state['subc']['dot_size'] = 25
        st.session_state['subc']['organism'] = 'No'
        st.session_state['subc']['gene_expression'] = None
        st.session_state['subc']['sub_adata_celltype'] = None
    out_dir = st.session_state['subc']['out_dir']

    context, _ = st.columns([8, 0.5])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader("🍓 Subclustering scRNA-seq datasets")
        st.write("After we have completed the scRNA-seq workflow and identified the various cell types present in our samples, we might decide that for a particular cell type, we would like to identify subtypes. For example, if we have a large cluster of CD4+ Helper T cells, we may want to identify subsets of Th1, Th2, Th17, Th9, and Tfh cells. To identify these cell subsets, we would subset the dataset to the cell type(s) of interest (e.g. CD4+ Helper T cells). ")
        uploaded_file = st.file_uploader(label="h5ad文件上传:", type='h5ad')
    if uploaded_file:
        adata = sc.read_h5ad(uploaded_file)
        st.write(adata)

        cols = st.columns([4, 8])
        with rc_context({'figure.figsize': (6, 6)}):
            sc.pl.umap(adata, 
                       color=["cell_type"],
                       frameon=True,
                       legend_loc='on data', 
                       legend_fontsize=9, 
                       legend_fontoutline=3, 
                       show=False)
            plt.savefig(os.path.join(out_dir, 'umap_all_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
        image_view(cols[0], png_file=os.path.join(out_dir, 'umap_all_celltype.png'),resize=(400,400))
            
        df_cell_num = adata.obs['cell_type'].value_counts().rename_axis('cell_type').reset_index(name='cell_counts')
        img_cell_num = px.bar(df_cell_num, 
                              x ='cell_type', 
                              y = 'cell_counts', 
                              color = 'cell_type', 
                              title='Number cells of every cell type.', 
                              text = 'cell_counts')
        cols[1].plotly_chart(img_cell_num, use_container_width=True)

        st.subheader("① Re-cluster")
        form_recluster(adata)
        if st.session_state['subc']['sub_adata']:
            form_gene_expression(adata)
            form_rename()
    else:
        context.info('上传.h5ad文件, 选择细胞类型cell_type继续亚群细分。', icon="ℹ️")
