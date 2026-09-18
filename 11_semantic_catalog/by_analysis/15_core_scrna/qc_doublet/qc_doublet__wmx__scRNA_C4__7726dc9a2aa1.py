import re
import os
import shutil
import pandas as pd
from PIL import Image
import streamlit as st
import streamlit_antd_components as sac
import streamlit.components.v1 as components
from config import out_dir, data_dir, marker_dict_file, rmd_template_C4, rmd_template_img_C4
from moudle.C4_scanpy_pipeline import make_report_dir, get_C4_mtx, data_qc, cell_filtering, dim_reduction_cluster, find_marker_gene, celltype_annotation, visua_gene_express

os.environ['OPENBLAS_NUM_THREADS'] = "64"
os.environ['MKL_NUM_THREADS'] = "64"
os.environ['OMP_NUM_THREADS'] = "64"

horizontal_bar = "<hr style='margin-top: 2; margin-bottom: 0; height: 1px; border: 1px solid #635985;'><br>"
hide_table_row_index = """
                        <style>
                            thead tr th:first-child {display:none}
                            tbody th {display:none}
                        </style>
                        """


def image_view(_container, png_file, caption='', resize=None):
    image = Image.open(png_file)
    if resize:
        new_image = image.resize(resize)
        _container.image(new_image, caption=caption)
    else:
        _container.image(image, caption=caption, use_container_width=True)


def data_server(pid, cols):
    tmp = []
    input_dir = os.path.join(data_dir, pid)
    sample_file = os.path.join(input_dir, 'sample.txt')
    st.session_state['out_dir'] = os.path.join(out_dir, pid)
    metrics_summary = os.path.join(st.session_state['out_dir'], '01-Dnbc4tools', 'metrics_summary.csv')
    
    if not os.path.exists(sample_file):
        tmp = []
        for sample in os.listdir(input_dir):
            if os.path.isdir(os.path.join(input_dir, sample)):
                tmp.append({"sample": sample, "rename": sample, "group": sample})
        sample_df = cols[0].data_editor(pd.DataFrame(tmp), hide_index=True, disabled=('sample',), use_container_width=False)             
        sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
        sample_df.to_csv(os.path.join(input_dir, 'sample.txt'), sep='\t', index=None)
    else:
        df = pd.read_csv(sample_file, sep='\t')
        sample_df = cols[0].data_editor(df, hide_index=True, disabled=('sample',), use_container_width=True)             
        sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
        sample_df.to_csv(os.path.join(input_dir, 'sample.txt'), sep='\t', index=None)

    submit = cols[0].button('🚀Submit', type='primary')
    
    if submit:
        with cols[1]:
            with st.spinner('数据加载中...'):
                make_report_dir(st.session_state['out_dir'])
                sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
                adata = get_C4_mtx(data_path=input_dir, sample_dict=sample_dict, out_dir=st.session_state['out_dir'])
                st.session_state['adata_qc'] = data_qc(adata, os.path.join(st.session_state['out_dir'], '02-Data_QC'))

    if 'adata_qc' in st.session_state:
        adata = st.session_state['adata_qc']
        qc_out_dir = os.path.join(st.session_state['out_dir'], '02-Data_QC')
        tab1, tab2, tab3 = cols[1].tabs(["🔔Summary", "📜violin plot", "📈scatter plot"])
        with tab1:
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("细胞总数", adata.X.shape[0])
            col2.metric("基因总数", adata.X.shape[1])
            col3.metric("样本个数", len(adata.obs['sample'].unique()))
            col4.metric("样本分组", len(adata.obs['group'].unique()))
            st.markdown(hide_table_row_index, unsafe_allow_html=True)
            # CellRanger Metrics
            metrics_summary = os.path.join(st.session_state['out_dir'], '01-Dnbc4tools', 'metrics_summary.csv')
            sample_list = list(adata.obs['sample'].drop_duplicates())
            try:
                metrics_sample = []
                for s in sample_list:
                    metrics = pd.read_csv(os.path.join(input_dir, s, 'output', 'metrics_summary.xls'),sep='\t')
                    metrics['Sample Name'] = s
                    cols = list(metrics)
                    cols = [cols[-1]] + cols[:-1]
                    metrics = metrics[cols]
                    metrics_sample.append(metrics)
                pd.concat(metrics_sample, ignore_index=True).to_csv(metrics_summary, index=None)
            except:
                pass

            if os.path.exists(metrics_summary):
                metrics_data = pd.read_csv(metrics_summary, sep=',')
                st.table(metrics_data)
            else:
                st.error(f' error: {metrics_summary}文件不存在!', icon="🚨")

        # QC
        with tab2:
            img_violin = Image.open(os.path.join(qc_out_dir, 'orig_nCount_nFeature_voilon.png'))
            st.image(img_violin, caption='注：nGene图：过滤前各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤前各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤前各个样本单个细胞中线粒体基因表达量的百分比分布情况；percent.ribo图：过滤前各个样本单个细胞中核糖体基因表达量的百分比分布情况。')

        with tab3:
            img_scatter = Image.open(os.path.join(qc_out_dir, 'orig_nCount_nFeature_scatter.png'))
            st.image(img_scatter, caption='注：左图表示细胞中检测到的UMI数目和线粒体基因含量的关系；右图表示细胞中检测到的UMI数目和表达基因数目的关系。')
    
    st.markdown(horizontal_bar, True)
    
      
def cell_filter():
    if 'adata_qc' in st.session_state: 
        st.subheader('② Cell filtering')
        out_dir = os.path.join(st.session_state['out_dir'], '03-Cell_filter')
        params, response = st.columns([2, 6])
        with params.form('cell_filter'):
            cols = st.columns(2)
            min_genes = cols[0].number_input(label='min_genes:', min_value=0, max_value=5000, value=200, step=1)
            st.session_state['min_genes'] = min_genes
            min_cells = cols[1].number_input(label='min_cells:', min_value=0, max_value=2000, value=10, step=1)
            pct_counts_mt = cols[0].number_input(label='pct_counts_mt:', min_value=0, max_value=100, value=20, step=1)
            st.session_state['pct_counts_mt'] = pct_counts_mt
            pct_counts_ribo = cols[1].number_input(label='pct_counts_ribo:', min_value=0, max_value=100, value=30, step=1)
            upper_lim = cols[0].number_input(label='max_genes:', min_value=500, max_value=20000, value=7000, step=100)
            st.session_state['upper_lim'] = upper_lim
            total_counts = cols[1].number_input(label='total_counts:', min_value=500, max_value=150000, value=25000, step=100)
            expected_doublet_rate = cols[0].number_input(label='doublet_rate(%):', min_value=0, max_value=100, value=5, step=1)               
            doublet_score = cols[1].number_input(label='doublet_score:', min_value=0.1, max_value=1.0, value=0.25, step=0.01)               
            submitted = st.form_submit_button("🚀Filter", type='primary')

        if submitted:
            with response:
                with st.spinner('Cell filtering calculating...'):
                    adata_qc = st.session_state['adata_qc'].copy()
                    st.session_state['adata_filtered'] = cell_filtering(adata_qc, 
                                                                        min_genes,
                                                                        min_cells,
                                                                        pct_counts_mt,
                                                                        pct_counts_ribo,
                                                                        upper_lim,
                                                                        total_counts,
                                                                        expected_doublet_rate,
                                                                        doublet_score,
                                                                        out_dir)

        if 'adata_filtered' in st.session_state: 
            adata_qc = st.session_state['adata_qc'].copy()
            adata_filtered = st.session_state['adata_filtered'].copy()
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("剩余细胞数", adata_filtered.X.shape[0], f'-{adata_qc.X.shape[0]-adata_filtered.X.shape[0]}')
                col2.metric("剩余基因数", adata_filtered.X.shape[1], f'-{adata_qc.X.shape[1]-adata_filtered.X.shape[1]}')
                tabs = st.tabs(["🚀violin plot (filtered)", "📈Scatter plot (filtered)", "🟠Scrublet score distribution"])
                img_filter_violin = Image.open(os.path.join(out_dir, 'filtered_nCount_nFeature_voilon.png'))
                tabs[0].image(img_filter_violin, caption='注：nGene图：过滤后各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤后各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤后各个样本单个细胞中线粒体基因表达量的百分比分布情况；percent.ribo图：过滤后各个样本单个细胞中核糖体基因表达量的百分比分布情况。', use_container_width=True)
                img_filter_scatter = Image.open(os.path.join(out_dir, 'filtered_nCount_nFeature_scatter.png'))
                tabs[1].image(img_filter_scatter, caption='注：nGene图：过滤后各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤后各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤后各个样本单个细胞中线粒体基因表达量的百分比分布情况；percent.ribo图：过滤后各个样本单个细胞中核糖体基因表达量的百分比分布情况。', use_container_width=True)  
                try:
                    img_scrublet = Image.open(os.path.join(out_dir, 'scrublet_score_distribution.png'))
                    tabs[2].image(img_scrublet, caption='注：双细胞打分分布图', use_container_width=True) 
                except:
                    pass                 
        else:
            with response.expander(f"**:blue[参数说明📜]** ", expanded=True):
                st.markdown("""
                            - **min_genes**: 去除表达基因数小于n的细胞(n默认值200);
                            - **min_cells**: 去除低于n个细胞表达的基因(n默认值10);
                            - **pct_counts_mt**: 过滤线粒体基因占比大于20%的细胞;
                            - **pct_counts_ribo**: 过滤核糖体基因占比大于30%的细胞;
                            - **max_genes**: 过滤基因总数大于n的细胞(n默认值为细胞总数的95%分位数);
                            - **doublet_rate(%)**: 期望双细胞或多细胞占总细胞数百分比(默认值以10x技术平台官方给出的1000个细胞多胞率进行计算);
                            """)
        
        st.markdown(horizontal_bar, True)
    
    
def dim_reduction_and_cluster():
    if 'adata_filtered' in st.session_state:
        st.subheader('③ Dimensionality Reduction and Clustering')
        out_dir = os.path.join(st.session_state['out_dir'], '04-Dimension_reduction_cluster')
        params, response = st.columns([2, 6])
        with params.form('clustering'):
            cols = st.columns(2)
            cluster_method = cols[0].selectbox(label='cluster_method:', options = ('leiden', 'louvain'), index=0, help = '选择聚类算法，默认使用leiden',)
            n_pcs = cols[1].number_input(label='n_pcs:', min_value=5, max_value=50, value=20, step=1)
            n_neighbors = cols[0].number_input(label='n_neighbors:', min_value=10, max_value=50, value=20, step=1)
            resolution = cols[1].number_input(label='resolution:', min_value=0.1, max_value=2.0, value=0.5, step=0.1)
            resolution = round(resolution, 1)
            remove_batch = st.radio("remove batch effect:",
                                    ["**No**", "**bbknn**", "**harmony**"], 
                                    captions = ["No batch effect.", "Use bbknn method.", "Use harmony method."], 
                                    help='选择去除批次效应方法。')
                       
            submitted = st.form_submit_button("🚀Clustering", type='primary')
            
        if submitted:
            with response:
                with st.spinner('Calculating...'):
                    adata_filtered = st.session_state['adata_filtered'].copy()
                    st.session_state['adata_cluster'] = dim_reduction_cluster(adata_filtered, 
                                                                              cluster_method,
                                                                              n_pcs,
                                                                              n_neighbors,
                                                                              resolution,
                                                                              remove_batch,
                                                                              out_dir)

        if 'adata_cluster' in st.session_state:
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                tabs = st.tabs(["🍊 Cluster", "🍊 Correlation and PCA", "🍊 Cell Number statistics", "🍊 PCA_loading", "🍊 Cluster Sample/Group"])
                cols = tabs[0].columns([2,2])
                image_view(cols[0], png_file=os.path.join(out_dir, 'umap_dimplot.png'), caption='样本整体UMAP降维聚类结果', resize=(500,500))
                image_view(cols[1], png_file=os.path.join(out_dir, 'cluster_total_counts_n_genes_by_counts.png'), caption='各细胞鉴定到的counts数及基因数', resize=(500,500))  
                cols = tabs[1].columns([2,2])
                image_view(cols[0], png_file=os.path.join(out_dir, 'cluster_correlation.png'), caption='各细胞Cluster相关性', resize=(600,600))
                image_view(cols[1], png_file=os.path.join(out_dir, 'pca_sample.png'), caption='PCA', resize=(600,600))            
                cell_statistics = pd.read_csv(os.path.join(out_dir, 'cluster_cells_statistics.xls'), sep='\t')
                tabs[2].data_editor(cell_statistics, hide_index=True, disabled=True, use_container_width=True)
                image_view(tabs[3], png_file=os.path.join(out_dir, 'pca_variance_ratio.png'), caption='注：PCA载荷图', resize=(600,600))
                if os.path.exists(os.path.join(out_dir, 'umap_split_sample_cluster.png')):
                    img_umap_sample = Image.open(os.path.join(out_dir, 'umap_split_sample_cluster.png'))
                    img_umap_group = Image.open(os.path.join(out_dir, 'umap_split_group_cluster.png'))
                    tabs[4].image(img_umap_sample, caption='注：按样本展示UMAP', use_container_width=True)
                    tabs[4].image(img_umap_group, caption='注：按分组展示UMAP', use_container_width=True)
        else:
            with response.expander(f"**:blue[参数说明📜]** ", expanded=True):
                st.markdown("""
                            - **cluster_method**: 单细胞图聚类方法，leiden or louvain(默认使用leiden);
                            - **n_pcs**: 用于构建邻域图的主成分数量。默认使用前20个主成分；
                            - **n_neighbors**: The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved;
                            - **resolution**: 细胞聚类分辨率(默认值为0.5);
                            - **remove batch effect**: 多样本合并分析时，选择是否消除批次效应，以及使用的算法;
                            """) 
        
        st.markdown(horizontal_bar, True)


def marker_gene():
    if 'adata_cluster' in st.session_state:
        st.subheader('④ Marker Gene')
        out_dir = os.path.join(st.session_state['out_dir'], '05-Marker_gene_analysis')
        params, response = st.columns([2, 6])
        with params.form('marker_gene'):
            cols = st.columns(2)
            rank_genes_groups_method = cols[0].selectbox(label='method:', options = ('wilcoxon', 't-test'), help = '选择计算Marker基因方法，默认使用wilcoxon')
            rank_genes_groups_fc = cols[1].number_input(label='log2fc:', min_value=0.1, max_value=5.0, value=0.5, step=0.1)
            rank_genes_groups_fdr = cols[0].number_input(label='fdr:', min_value=0.0, max_value=1.0, value=0.05, step=0.01)
            rank_genes_groups_pct = cols[1].number_input(label='pct:', min_value=0.0, max_value=1.0, value=0.25, step=0.01)
            submitted = st.form_submit_button("🚀Marker gene", type='primary')
        if submitted:
            with response:
                with st.spinner('Calculating...'):
                    adata_cluster = st.session_state['adata_cluster']
                    st.session_state['adata_rank_genes_groups'] = find_marker_gene(adata_cluster,
                                                                                   method=rank_genes_groups_method,
                                                                                   logfc=rank_genes_groups_fc,
                                                                                   fdr=rank_genes_groups_fdr,
                                                                                   pct=rank_genes_groups_pct,
                                                                                   out_dir=out_dir)   

        if 'adata_rank_genes_groups' in st.session_state:
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                tabs = st.tabs(["🍊 dotplot", "🍊 heatmap", "🍊 marker genes"])
                image_view(tabs[0], png_file=os.path.join(out_dir, 'dotplot_top3_marker_gene.png'), caption='Top3 marker gene dotplot')
                image_view(tabs[1], png_file=os.path.join(out_dir, 'heatmap_top10_marker_gene.png'), caption='Top10 marker gene heatmap')
                marker_gene = pd.read_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index_col='cluster')
                cluster = tabs[2].selectbox(label = 'Select cluster view marker gene:', options = sorted(marker_gene.index.drop_duplicates().to_list()))
                if not cluster:
                    tabs[2].dataframe(marker_gene, use_container_width=True)
                else:    
                    df_marker = marker_gene.loc[cluster]
                    tabs[2].dataframe(df_marker, use_container_width=True)
        else:
            with response.expander("**:blue[参数说明📜]** ", expanded=True):
                st.markdown("""
                            - **rank_genes_groups_method**: 选择计算基因差异的方法(默认使用wilcoxon);
                            - **log2fc**: 差异基因log2 FoldChange过滤阈值;
                            - **fdr**: 差异基因FDR过滤阈值;
                            - **pct**: 差异基因在处理组细胞中表达占比;
                            """)

        st.markdown(horizontal_bar, True)


def cell_type_annotation():
    if 'adata_rank_genes_groups' in st.session_state:
        st.subheader('⑤ Cell Type Annotation')
        out_dir = os.path.join(st.session_state['out_dir'], '06-Cell_type_annotation')
        params, response = st.columns([2, 6])
        with params.form('annotation'):
            cols = st.columns(2)
            organism = cols[0].selectbox(label = 'Organism:', options = ('human', 'mouse'), index=None, help='选择样本物种信息，只支持人和小鼠。', key='opt_organism')
            tissue = cols[1].selectbox(label = 'Tissue type:', options = ['All', 'bone', 'breaste', 'lung', 'pbmc'], help='主要用于后续细胞类型鉴定数据库比对，默认使用所有组织类型。', key='opt_tissue')        
            database = cols[0].selectbox(label = 'Database:', options = ('PanglaoDB', 'CellTypist'), help='选择细胞类型自动注释数据库', key='opt_database')               
            submitted = st.form_submit_button("🚀Annot", type='primary')
        if submitted:
            if not organism:
                params.exception(RuntimeError('🚨 You must choose sample organism to annotation!'))
            else:
                st.session_state['organism'] = organism
                with response:
                    with st.spinner('Annotation...'):
                        adata_rank_genes_groups = st.session_state['adata_rank_genes_groups']
                        try:
                            st.session_state['adata_cell_type'], st.session_state['annotation_dict'] = celltype_annotation(adata_rank_genes_groups,
                                                                                                                           database,
                                                                                                                           organism,
                                                                                                                           out_dir)
                        except:
                            st.error('物种信息是不是选错了！', icon="🚨")
                                   
    if 'adata_cell_type' in st.session_state:
        with params.form('annotation_rename'):
            df = pd.DataFrame(st.session_state['annotation_dict'].items(), columns=['Cluster', 'annotation'])          
            df['cell_type'] = df['annotation']
            edited_df = st.data_editor(df, hide_index=True, disabled=('Cluster','annotation'), use_container_width=True)   
            renames = edited_df.set_index('Cluster')['cell_type'].to_dict() 
            submitted = st.form_submit_button("🚀Rname", type='primary')
            if submitted:
                with response:
                    with st.spinner('Renaming...'):
                        adata_rank_genes_groups = st.session_state['adata_rank_genes_groups']
                        st.session_state['adata_cell_type'], annotation_dict = celltype_annotation(adata_rank_genes_groups,
                                                                                                   database,
                                                                                                   organism,
                                                                                                   out_dir,
                                                                                                   renames)

        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            #tabs = st.tabs(["🍊 UMAP-all", "🍊 UMAP-group", "🍊 Cell type proportion of samples","🍊 Cell type proportion of groups"])
            tabs = st.tabs(["🍊 UMAP-all", "🍊 UMAP-group", "🍊 Cell type proportion of samples"])
            image_view(tabs[0], png_file=os.path.join(out_dir, 'umap_cluster_celltype.png'), caption='样本整体UMAP降维聚类结果')
            if os.path.exists(os.path.join(out_dir, 'umap_split_group_celltype.png')):
                image_view(tabs[1], png_file=os.path.join(out_dir, 'umap_split_group_celltype.png'), caption='样本分组UMAP降维聚类结果')
            else:
                pass
            image_view(tabs[2], png_file=os.path.join(out_dir, 'celltype_proportion_sample.png'), caption='各样本中细胞类型占比', resize=(600,600))
            #image_view(tabs[3], png_file=os.path.join(out_dir, 'celltype_proportion_group.png'), caption='各组中细胞类型占比', resize=(600,600))

        st.markdown(horizontal_bar, True)


def gene_expression():
    data = pd.read_csv(marker_dict_file, sep='\t')
    marker_dict = {}
    if 'adata_filtered' in st.session_state:   
        for index, row in data.iterrows():
            cell_type = row['cell_type']
            organism = row['organism']
            genes = row['markers'].split(',')
            if organism not in marker_dict:
                marker_dict[organism] = {}
            if cell_type not in marker_dict[organism]:
                marker_dict[organism][cell_type] = []  
            marker_dict[organism][cell_type].extend(genes)

        st.subheader('⑥ Gene Expression')
        out_dir = os.path.join(st.session_state['out_dir'], '06-Cell_type_annotation')
        params, response = st.columns([2, 6])
        types = params.selectbox(label = '选择基因展示/预定义细胞类型Markers：', options = ('选择基因（支持多基因）', '使用预定义细胞类型Marker基因'), help='自定义基因查找，或使用预先定义的细胞类型Markers。')
        with params.form('gene_exp'):
            select_types = []
            if types == '选择基因（支持多基因）':
                genes = st.multiselect(label='基因:', options=st.session_state['adata_filtered'].var_names)       
            else:
                genes = []
                celltyps = marker_dict[st.session_state['organism']].keys()
                select_types = st.multiselect(label='细胞类型:', options=celltyps)       
                if select_types:
                    for t in select_types:
                        genes.extend([i for i in marker_dict[st.session_state['organism']][t] if i in st.session_state['adata_filtered'].var_names])

            submitted = st.form_submit_button("🚀Plot", type='primary')

            if select_types:
                st.write({i: marker_dict[st.session_state['organism']][i] for i in select_types})
            else:
                st.write(genes)

            with response.expander(f"**:red[基因表达]** ", expanded=True):
                #tabs = st.tabs(["🍊 UMAP", "🍊 Dotplot"])
                tabs = st.tabs(["🍊 UMAP"])
                if submitted:
                    if not genes:
                        st.error('Please select genes to visualization of gene expression.', icon="🔥")
                    else:
                        with st.spinner('Calculating...'):
                            st.session_state['gene_expression'], fig = visua_gene_express(st.session_state['adata_cell_type'], genes, out_dir)

                if 'gene_expression' in st.session_state:
                    if genes:
                        png_file = os.path.join(out_dir, '_'.join(genes) + '.png')
                        image_view(tabs[0], png_file=os.path.join(out_dir, '_'.join(genes)+'.png'), caption=f"{genes}基因表达")
                #if 'gene_expression1' in st.session_state:
                #    image_view(tabs[1], png_file=os.path.join(out_dir, 'markers_dotplot.png'), caption=f"{genes}基因表达及占比")


        st.markdown(horizontal_bar, True)


def generate_report():
    if 'adata_cell_type' in st.session_state:
        # sac.divider(align='center', size='lg', color='cyan')
        st.subheader(f"🍑 生成报告, 项目编号: {st.session_state['pid']}, 文件路径: {os.path.abspath(st.session_state['out_dir'])}")
        out_dir = os.path.abspath(st.session_state["out_dir"])
        zip_outdir = os.path.abspath(os.path.join(out_dir, '..'))
        with st.form('report'): 
            cols = st.columns([3,5,4,3], gap='large')
            pid = cols[0].text_input('项目编号:', st.session_state['pid'])
            ptitle = cols[1].text_input('项目名称:', '')
            paddress = cols[2].text_input('客户单位:', '')
            pcustmor = cols[3].text_input('客户姓名:', '')               
            submitted = st.form_submit_button("🚀Report", type='primary')
        
        if submitted:
            rmd = ''
            with open(rmd_template_C4, 'r', encoding='utf-8') as f:
                line = f.read()
                if not os.path.exists(os.path.join(out_dir, '01-Dnbc4tools', 'metrics_summary.csv')):
                    line = re.sub(r'## 测序数据统计.*? 所有细胞的基因总数。', '', line, flags=re.DOTALL)
                sample_num = len(st.session_state['adata_qc'].obs['sample'].drop_duplicates().tolist())
                cell_num_raw = st.session_state['adata_qc'].X.shape[0]
                cell_num_after_filter = str(st.session_state['adata_filtered'].X.shape[0])
                cell_types = ', '.join(st.session_state['adata_cell_type'].obs['cell_type'].drop_duplicates().tolist())                   
                cluster_num = len(st.session_state['adata_cell_type'].obs["cluster"].drop_duplicates().tolist())
                rmd = line.replace('%pid%', pid).replace('%ptitle%', ptitle).replace('%paddress%', paddress).replace('%pcustmor%', pcustmor).replace('%num_sample%', str(sample_num)).replace('%organism%', st.session_state['organism'])
                rmd = rmd.replace('%cell_num_raw%', str(cell_num_raw)).replace('%cell_num_after_filter%', str(cell_num_after_filter)).replace('%cluster_num%', str(cluster_num)).replace('%cell_types%', cell_types) 
                rmd = rmd.replace('%mt_pct%', str(st.session_state['pct_counts_mt']))
                rmd = rmd.replace('%min_genes%', str(st.session_state['min_genes']))
                rmd = rmd.replace('%upper_lim%', str(st.session_state['upper_lim']))
                rmd = rmd.replace('images_dir', rmd_template_img_C4)

            with open(os.path.join(out_dir, 'report.rmd'), 'w', encoding='utf-8') as w:
                w.write(rmd)
            with st.spinner('Report generation ...'):   
                os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(out_dir))
            with st.spinner('File compression ...'):
                st.toast('分析结果压缩中，请稍等!', icon='💗')
                st.session_state['adata_cell_type'].write(os.path.join(out_dir, 'result.h5ad'))
                shutil.make_archive(os.path.join(zip_outdir, st.session_state['pid']), 'zip', out_dir)        

        if os.path.exists(os.path.join(out_dir, 'report.html')): 
            cols = st.columns(8)
            if os.path.exists(os.path.join(zip_outdir, st.session_state['pid']+'.zip')):
                cols[7].download_button(
                                        label="📁 Download",
                                        data=open(os.path.join(zip_outdir, st.session_state['pid']+'.zip'), "rb"),
                                        file_name=f"{st.session_state['pid']}.zip",
                                        mime="application/octet-stream",
                                        )
            report = open(os.path.join(out_dir, 'report.html'))
            components.html(report.read(), height=1000, scrolling=True)   

 
def main():
    st.subheader('🍄 C4 scRNA标准分析')
    sac.alert(
        label='**使用说明**',
        description=f'''**选择项目编号，读取数据文件，按照分析步骤修改参数，最终自动生成标准分析报告。**''', 
        color='cyan',
        icon=True,
        size=15
    )
    tab1, tab2 = st.tabs(["📂 服务器数据加载", "📤 文件上传"])
    with tab1:
        st.subheader('① Data Loading')
        st.info('根据项目编号，加载/mnt/NexData_file/nas/outputdata/cellranger_out路径下CellRanger输出结果文件，选择项目编号后会根据样本名称生成样本信息表，用户可以在样本信息表中修改样本分组信息，用于后续分析', icon="📢")
        cols = st.columns([2, 6], gap='medium')
        pid = cols[0].selectbox('选择项目编号：',
                                os.listdir(data_dir),
                                index=None,
                                placeholder="Select project",
                                key='pid',
                                help='选择项目编号，加载CellRanger结果。')    
        if pid:
            cols[0].success(f"已选择项目编号: {pid}", icon="📂") 
            data_server(pid, cols)
            cell_filter()
            dim_reduction_and_cluster()
            marker_gene()
            cell_type_annotation()
            gene_expression()
            generate_report()


if __name__ == '__main__':
    main()
