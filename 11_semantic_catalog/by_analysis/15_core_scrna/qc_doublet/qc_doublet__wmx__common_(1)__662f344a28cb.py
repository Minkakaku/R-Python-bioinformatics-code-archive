import os
os.environ['OPENBLAS_NUM_THREADS'] = "64"
os.environ['MKL_NUM_THREADS'] = "64"
os.environ['OMP_NUM_THREADS'] = "64"
import re
import shutil
import random
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import streamlit as st
from PIL import Image
import celltypist
from celltypist import models
import scanpy.external as sce
import matplotlib.pyplot as plt
import streamlit_antd_components as sac
import streamlit.components.v1 as components
from matplotlib.pyplot import rc_context
from matplotlib import font_manager as fm
import warnings
from config import marker_dict_file
from moudle.scanpy_pipeline_10X import make_report_dir, get_10x_mtx
sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA_platform/scVision5/tmp/.scanpy_cached"
sc.settings.set_figure_params(figsize=(5,5), facecolor='white')
plt.tight_layout()


horizontal_bar = "<hr style='margin-top: 2; margin-bottom: 0; height: 1px; border: 1px solid #635985;'><br>"
hide_table_row_index = """
                        <style>
                            thead tr th:first-child {display:none}
                            tbody th {display:none}
                        </style>
                        """


def set_scanpy_color_palette(_adata, color_key, color_palette=None):
    """
    为 Scanpy 对象设置统一配色方案。
    
    参数:
    - adata: AnnData 对象
    - color_key: str，obs 中的列名，用于分组着色（如 'cell_type'）
    - color_palette: list，可选，自定义的颜色列表（最多 50 个）

    默认色板支持最多 50 种颜色，超出会报错。
    """
    default_palette = [
        "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
        "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
        "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
        "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
        "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
        "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
        "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
        "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
        "#f17d80", "#737495"
    ]
    palette = color_palette if color_palette else default_palette

    if color_key not in _adata.obs:
        raise KeyError(f"'{color_key}' 不存在于 adata.obs 中")

    categories = _adata.obs[color_key].astype(str).unique().tolist()
    n_categories = len(categories)

    if n_categories > len(palette):
        raise ValueError(f"类别数为 {n_categories}，超过了可用颜色数 {len(palette)}，请扩展色板。")

    # 设置为 Categorical 类型并指定顺序
    _adata.obs[color_key] = _adata.obs[color_key].astype("category")
    _adata.obs[color_key] = _adata.obs[color_key].cat.set_categories(categories)

    # 设置颜色
    _adata.uns[f"{color_key}_colors"] = palette[:n_categories]

def plot_celltype_proportion(_adata, out_dir,type='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = _adata.obs['cell_type'].cat.categories.to_list()
    values = list(_adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 创建堆积条形图
    ax = pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )

    # 设置X轴和Y轴标签
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')  # 设置X轴标签
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')  # 设置Y轴标签

    # 设置图例
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )

    # 设置坐标轴标签的字体大小和颜色
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    # 设置背景颜色
    ax.set_facecolor('white')

    # 设置图表边框加粗
    for spine in ax.spines.values():
        spine.set_edgecolor('black')  # 设置边框颜色
        spine.set_linewidth(1.5)

    # 设置图像大小，如果没有传入fig_size，则使用默认大小
    if fig_size is None:
        fig_size = {'width': len(set(_adata.obs[type])) * 2, 'height': 6}

    # 调整布局，确保图例不遮挡图形
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    # 保存图像
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)



def image_view(_container, png_file, caption='', resize=None):
    image = Image.open(png_file)
    if resize:
        new_image = image.resize(resize)
        _container.image(new_image, caption=caption)
    else:
        _container.image(image, caption=caption, use_container_width=True)

def cell_filter():
    if 'adata_qc' in st.session_state: 
        st.subheader('② Cell filtering')
        out_dir = os.path.join(st.session_state['out_dir'], '03-Cell_filter')
        params, response = st.columns([2, 6])
        with params.form('cell_filter'):
            cols = st.columns(2)
            min_genes = cols[0].number_input(label='min_genes:', min_value=0, max_value=5000, value=200, step=1)
            min_cells = cols[1].number_input(label='min_cells:', min_value=0, max_value=2000, value=10, step=1)
            pct_counts_mt = cols[0].number_input(label='pct_counts_mt:', min_value=0, max_value=100, value=20, step=1)
            pct_counts_ribo = cols[1].number_input(label='pct_counts_ribo:', min_value=0, max_value=100, value=30, step=1)
            upper_lim = cols[0].number_input(label='max_genes:', min_value=500, max_value=20000, value=7000, step=100)
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
            tabs = st.tabs(["🍊 UMAP-all", "🍊 UMAP-group", "🍊 Cell type proportion of samples","🍊 Cell type proportion of groups"])
            image_view(tabs[0], png_file=os.path.join(out_dir, 'umap_cluster_celltype.png'), caption='样本整体UMAP降维聚类结果')
            if os.path.exists(os.path.join(out_dir, 'umap_split_group_celltype.png')):
                image_view(tabs[1], png_file=os.path.join(out_dir, 'umap_split_group_celltype.png'), caption='样本分组UMAP降维聚类结果')
            else:
                pass
            image_view(tabs[2], png_file=os.path.join(out_dir, 'celltype_proportion_sample.png'), caption='各样本中细胞类型占比', resize=(600,600))
            image_view(tabs[3], png_file=os.path.join(out_dir, 'celltype_proportion_group.png'), caption='各组中细胞类型占比', resize=(600,600))

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
                tabs = st.tabs(["🍊 UMAP", "🍊 Dotplot"])
                #tabs = st.tabs(["🍊 UMAP"])
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
                        # dotplot 图展示
                        dotplot_file = os.path.join(out_dir, '_'.join(genes) + '_dotplot.png')
                        image_view(tabs[1], png_file=dotplot_file, caption=f"{genes} Dotplot")


        st.markdown(horizontal_bar, True)

def generate_report(rmd_template_path,rmd_template_img):
    if 'adata_cell_type' in st.session_state:
        st.subheader(f"🍑 生成报告, 项目编号: {st.session_state['pid']}, 文件路径: {os.path.abspath(st.session_state['out_dir'])}")
        out_dir = os.path.abspath(st.session_state["out_dir"])
        zip_outdir = os.path.abspath(os.path.join(out_dir, '..'))

        # 表单输入
        with st.form('report'): 
            cols = st.columns([3,5,4,3], gap='large')
            pid = cols[0].text_input('项目编号:', st.session_state['pid'])
            ptitle = cols[1].text_input('项目名称:', '')
            paddress = cols[2].text_input('客户单位:', '')
            pcustmor = cols[3].text_input('客户姓名:', '')               
            submitted = st.form_submit_button("🚀Report", type='primary')
        
        if submitted:
            # 保存字符串型 session_state 变量
            session_str_vars = {k: v for k, v in st.session_state.items() if isinstance(v, str)}
            csv_path = os.path.join(out_dir, 'session_vars.csv')
            pd.DataFrame(session_str_vars.items(), columns=['Key', 'Value']).to_csv(csv_path, index=False)

            # 读取并替换 RMD 模板
            with open(rmd_template_path, 'r', encoding='utf-8') as f:
                rmd_content = f.read()
                # 如果 metrics_summary.csv 缺失，移除相关部分
                if not os.path.exists(os.path.join(out_dir, '01-CellRanger', 'metrics_summary.csv')):
                    rmd_content = re.sub(r'## 测序数据统计.*? 所有细胞的基因总数。', '', rmd_content, flags=re.DOTALL)
                
                # 获取各类信息替换模板内容
                sample_num = len(st.session_state['adata_qc'].obs['sample'].drop_duplicates())
                cell_num_raw = st.session_state['adata_qc'].X.shape[0]
                cell_num_after_filter = st.session_state['adata_filtered'].X.shape[0]
                cell_types = ', '.join(st.session_state['adata_cell_type'].obs['cell_type'].drop_duplicates().tolist())
                cluster_num = len(st.session_state['adata_cell_type'].obs['cluster'].drop_duplicates())

                replacements = {
                    '%pid%': pid,
                    '%ptitle%': ptitle,
                    '%paddress%': paddress,
                    '%pcustmor%': pcustmor,
                    '%num_sample%': str(sample_num),
                    '%organism%': st.session_state['organism'],
                    '%cell_num_raw%': str(cell_num_raw),
                    '%cell_num_after_filter%': str(cell_num_after_filter),
                    '%cluster_num%': str(cluster_num),
                    '%cell_types%': cell_types,
                    'images_dir': rmd_template_img  # 确保这个变量是可访问的
                }

                for key, val in replacements.items():
                    rmd_content = rmd_content.replace(key, val)

            # 写入替换后的 RMD 文件
            rmd_output = os.path.join(out_dir, 'report.rmd')
            with open(rmd_output, 'w', encoding='utf-8') as w:
                w.write(rmd_content)

            # 渲染 RMarkdown
            with st.spinner('Report generation ...'):
                os.system(f"R -e \"rmarkdown::render('{rmd_output}', output_file='{out_dir}/report.html')\"")

            # 结果打包
            with st.spinner('File compression ...'):
                st.toast('分析结果压缩中，请稍等!', icon='💗')
                st.session_state['adata_cell_type'].write(os.path.join(out_dir, 'result.h5ad'))
                shutil.make_archive(os.path.join(zip_outdir, pid), 'zip', out_dir)

        # 报告展示与下载
        report_path = os.path.join(out_dir, 'report.html')
        zip_path = os.path.join(zip_outdir, f"{st.session_state['pid']}.zip")
        if os.path.exists(report_path):
            cols = st.columns(8)
            if os.path.exists(zip_path):
                cols[7].download_button(
                    label="📁 Download",
                    data=open(zip_path, "rb"),
                    file_name=f"{st.session_state['pid']}.zip",
                    mime="application/octet-stream"
                )
            report = open(report_path)
            components.html(report.read(), height=1000, scrolling=True)


def annotate_clusters_by_panglao(_adata, cluster_key='cluster', organism='human', obsm_key='ora_estimate', min_n=3):
    """
    使用 PanglaoDB marker 基于 cluster 平均表达执行 ORA 注释。

    参数：
    - _adata: AnnData 对象，包含已聚类的结果。
    - cluster_key: 用于聚类的 obs 列（默认为 'cluster'）。
    - organism: 'human' 或 'mouse'。
    - obsm_key: enrichment score 的 obsm 键名。
    - min_n: ORA 的最小 marker 数目。

    返回：
    - _adata: 增加了 'cell_type' 注释的 AnnData。
    - annotation_dict: cluster → cell_type 的映射字典。
    """
    # 保证 cluster 是字符串
    _adata.obs[cluster_key] = _adata.obs[cluster_key].astype(str)

    # 获取 PanglaoDB markers 并去重
    markers = dc.get_resource('PanglaoDB')
    markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]

    # 聚合为 cluster 平均表达矩阵
    cluster_expr = _adata.to_df().groupby(_adata.obs[cluster_key]).mean()
    _adata_cluster = sc.AnnData(cluster_expr)
    _adata_cluster.var_names = cluster_expr.columns.astype(str)
    _adata_cluster.obs_names = cluster_expr.index.astype(str)

    # 执行 ORA
    dc.run_ora(
        mat=_adata_cluster,
        net=markers,
        source='cell_type',
        target='genesymbol',
        min_n=min_n,
        use_raw=False,
        verbose=True
    )

    # enrichment 分数矩阵
    acts = dc.get_acts(_adata_cluster, obsm_key=obsm_key)
    ora_scores = acts.obsm['ora_estimate']

    # 构造成 DataFrame，行为 cluster，列为 cell_type
    acts_df = pd.DataFrame(
    ora_scores,
    index=acts.obs_names,
    columns=acts.var_names
    )

    # 自动注释（最大 enrichment cell type）
    annotation_dict = dc.assign_groups(acts_df)

    # 映射到每个细胞
    _adata.obs['cell_type'] = _adata.obs[cluster_key].map(annotation_dict)

    return _adata, annotation_dict


def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

@st.fragment
def data_qc(_adata, out_dir):   
    _adata.var['mt'] = _adata.var_names.str.upper().str.startswith('MT-') 
    print(_adata.var_names)
    _adata.var['ribo'] = _adata.var_names.str.upper().str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(_adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=90, multi_panel=True, show=False)
    axes = plt.gcf().axes
    for ax in axes:
        ax.set_xticklabels(ax.get_xticklabels(), horizontalalignment='right')
    save_fig(plt, out_dir, file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(_adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
    sc.pl.scatter(_adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
    save_fig(fig, out_dir, file_name='orig_nCount_nFeature_scatter')
    
    pd.DataFrame(_adata.obs).to_csv(os.path.join(out_dir, 'all_cell_features_raw.xls'), sep='\t')

    for s in _adata.obs['sample'].drop_duplicates().tolist():
        os.makedirs(os.path.join(out_dir, s), exist_ok=True)
        sc.pl.violin(_adata[_adata.obs['sample']==s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=90, multi_panel=True, show=False)
        save_fig(plt, os.path.join(out_dir, s), file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, os.path.join(out_dir, s), file_name='orig_nCount_nFeature_scatter')
        pd.DataFrame(_adata[_adata.obs['sample']==s].obs).to_csv(os.path.join(out_dir, s, 'cell_features_raw.xls'), sep='\t')
    
    return _adata

@st.fragment
def cell_filtering(_adata, min_genes, min_cells, pct_counts_mt, pct_counts_ribo, upper_lim, total_counts, expected_doublet_rate, threshold, out_dir):    
    sc.pp.filter_cells(_adata, min_genes=min_genes)
    sc.pp.filter_genes(_adata, min_cells=min_cells)
    _adata = _adata[(_adata.obs.n_genes_by_counts < upper_lim)]
    _adata = _adata[_adata.obs['pct_counts_mt'] < pct_counts_mt, :]
    _adata = _adata[_adata.obs['pct_counts_ribo'] < pct_counts_ribo, :]
    _adata = _adata[_adata.obs['total_counts'] < total_counts, :]
    
    sc.external.pp.scrublet(_adata, expected_doublet_rate=expected_doublet_rate / 100, threshold=threshold, batch_key='sample', random_state=1234)
    set_scanpy_color_palette(_adata, color_key="group",color_palette=["#feee7d", "#e8a0b8", "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", "#f17d80", "#737495"])
    samples = _adata.obs['sample'].drop_duplicates().tolist()
    if len(samples) > 1:
        sc.external.pl.scrublet_score_distribution(_adata)
        save_fig(plt, out_dir, file_name='scrublet_score_distribution', fig_size={'width':8, 'height':6})  
    _adata = _adata[_adata.obs['predicted_doublet']==False, :]

    pd.DataFrame(_adata.obs).to_csv(os.path.join(out_dir, 'all_cell_features_filtered.xls'), sep='\t')
    sc.pl.violin(_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby = 'sample', rotation=90, multi_panel=True, show=False)
    save_fig(plt, out_dir, file_name='filtered_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})  
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(_adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
    sc.pl.scatter(_adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
    save_fig(fig, out_dir, file_name='filtered_nCount_nFeature_scatter')      
    
    for s in _adata.obs['sample'].drop_duplicates().tolist():
        os.makedirs(os.path.join(out_dir, s), exist_ok=True)
        sc.pl.violin(_adata[_adata.obs['sample']==s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby = 'sample', rotation= 90, multi_panel=True, show=False)
        save_fig(plt, os.path.join(out_dir, s), file_name='filtered_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, os.path.join(out_dir, s), file_name='filtered_nCount_nFeature_scatter')       
        pd.DataFrame(_adata[_adata.obs['sample']==s].obs).to_csv(os.path.join(out_dir, s, 'cell_features_filtered.xls'), sep='\t')
    _adata.write(os.path.join(out_dir, 'adata_counts.h5ad'))

    return _adata


@st.fragment
def dim_reduction_cluster(_adata, method, n_pcs, n_neighbors, resolution, batch_effect, out_dir):
    _adata.layers["counts"] = _adata.X.copy()
    # data normalization and log10 transform
    sc.pp.normalize_total(_adata)
    sc.pp.log1p(_adata)
    _adata.raw = _adata
    # select highly variable gene and data scale 
    sc.pp.highly_variable_genes(_adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    sc.pl.highly_variable_genes(_adata, show=False)
    save_fig(plt, out_dir, file_name='highly_variable_genes', fig_size={'width':16, 'height':8})
    _adata = _adata[:, _adata.var['highly_variable']]
    sc.pp.regress_out(_adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(_adata, max_value=10)
    # PCA
    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pl.pca(_adata, color='sample', show=False)
    save_fig(plt, out_dir, file_name='pca_sample')
    sc.pl.pca(_adata, color='group', show=False)
    save_fig(plt, out_dir, file_name='pca_group')
    sc.pl.pca_variance_ratio(_adata, log=True, show=False)
    save_fig(plt, out_dir, file_name='pca_variance_ratio')

    if batch_effect == '**harmony**':
        sce.pp.harmony_integrate(_adata, 'sample', adjusted_basis='X_pca')
    # cluster
    sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)    

    # Remove batch effect
    if batch_effect == '**bbknn**':
        sce.pp.bbknn(_adata, batch_key='sample', n_pcs=n_pcs)  # running bbknn 1.3.6


    sc.tl.umap(_adata)
    if method == 'leiden':
        sc.tl.leiden(_adata, resolution=resolution, key_added='cluster')
    if method == 'louvain':
        sc.tl.louvain(_adata, resolution=resolution, key_added='cluster')



    # UMAP
    set_scanpy_color_palette(_adata, color_key="cluster")
    sc.pl.umap(_adata, 
               color='cluster',
               title='Cluster',
               legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=2, 
               show=False)
    save_fig(plt, out_dir, file_name='umap_dimplot')
    sc.pl.umap(_adata, color='sample', title='', show=False)
    save_fig(plt, out_dir, file_name='umap_sample_dimplot')
    sc.pl.umap(_adata, color='group', title='', show=False)
    save_fig(plt, out_dir, file_name='umap_group_dimplot')
      
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, 
                   color=['cluster', "total_counts", "n_genes_by_counts", "total_counts_mt"],
                   title=['Cluster', 'total_counts', 'n_genes_by_counts', "total_counts_mt"],
                   legend_loc='on data',
                   frameon=True,
                   legend_fontsize=9, 
                   legend_fontoutline=3, 
                   wspace=0.2, 
                   ncols=2, 
                   cmap='coolwarm', 
                   use_raw=True, 
                   show=False)
               
        plt.savefig(os.path.join(out_dir, 'cluster_total_counts_n_genes_by_counts.png'), facecolor='white', bbox_inches="tight", dpi=300)
    
    # write all cells features
    _adata.obs.to_csv(os.path.join(out_dir,'all_cells_features.xls'), sep='\t')
    
    # write cluster cell statistics result
    result = _adata.obs.groupby(['cluster', 'sample']).size().reset_index(name='number').pivot(index='cluster', columns='sample', values='number').reset_index()
    result.columns.name = None    
    result.to_csv(os.path.join(out_dir,'cluster_cells_statistics.xls'), sep='\t', index=None)
    
    # cluster correlation
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.correlation_matrix(_adata, 'cluster', show=False)
        plt.savefig(os.path.join(out_dir, 'cluster_correlation.png'), facecolor='white', bbox_inches="tight", dpi=300)

    groups = _adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["group"]==groups[i]], color="cluster", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_group_cluster', fig_size={'width':len(groups)*5, 'height':5})

    samples = _adata.obs['sample'].drop_duplicates().tolist()
    if len(samples) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
        fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
        for i in range(len(samples)):
            if i == len(samples) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["sample"]==samples[i]], color="cluster", title=samples[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_sample_cluster', fig_size={'width':len(samples)*5, 'height':5})



    return _adata

@st.fragment
def find_marker_gene(_adata, method, out_dir, logfc=1, fdr=0.05, pct=0.25):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
   
    sc.tl.rank_genes_groups(_adata, 'cluster', method=method, pts=True)
    sc.pl.rank_genes_groups(_adata, n_genes=25, sharey=False)
    save_fig(plt, out_dir, file_name='rank_genes_groups', fig_size={'width':20, 'height':12})
    
    result = _adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, f'cluster{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(_adata, group=c)
        tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)

        # enr_go = gp.enrichr(gene_list=list(tmp['names']),
        #         gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
        #         background='hsapiens_gene_ensembl',
        #         organism='Human', 
        #         cutoff=0.05,
        #         no_plot=True,
        #         outdir=out,
        #         verbose=False
        #         )

        # for f in os.listdir(out):
        #     if f.startswith('gseapy'):
        #         os.remove(os.path.join(out, f))

        sc.pl.rank_genes_groups_violin(_adata, groups=c, gene_names=list(tmp['names'].head(5)), jitter=0, size=0, show=False)
        save_fig(plt, out, file_name='top_5_genes_violin')

    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    # if st.session_state['pipe']['organism'] == 'human':
    #     gene_annot = os.path.join(annotDir, '9606_human_go_kegg_annotation.xls')
    # elif st.session_state['pipe']['organism'] == 'mouse':
    #     gene_annot = os.path.join(annotDir, '10090_mouse_go_kegg_annotation.xls')
    # gene_annot = pd.read_csv(gene_annot, sep='\t')
    # all_diff = pd.merge(all_diff, gene_annot, on='Symbol', how='left')
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'cluster_diff_statistics.xls'), sep='\t')

    sc.pl.rank_genes_groups_dotplot(
        _adata,
        n_genes=3,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=1,
        cmap='bwr',
        colorbar_title='log fold change',
        show=False
    )
    save_fig(plt, out_dir, file_name='dotplot_top3_marker_gene', fig_size={'width':20, 'height':8})

    sc.pl.rank_genes_groups_heatmap(
        _adata,
        n_genes=10,
        values_to_plot="logfoldchanges",
        min_logfoldchange=1,
        swap_axes=True,
        cmap='bwr',
        standard_scale='var',
        show_gene_labels=True,
        dendrogram=True,
        show=False
    )
    save_fig(plt, out_dir, file_name='heatmap_top10_marker_gene', fig_size={'width':26, 'height':20})

    #_adata.write(os.path.join(out_dir, 'adata_cluster.h5ad'))

    return _adata

@st.fragment
def celltype_annotation(_adata, database, organism, out_dir, cell_type=None):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)

    if organism == 'mouse':
        _adata.var['gene'] = _adata.var.index
        _adata.var.index = [i.upper() for i in _adata.var.index]
            
    if database == 'PanglaoDB':
        markers = dc.get_resource('PanglaoDB')        
        if cell_type:
            annotation_dict = cell_type
            _adata.obs['cell_type'] = [annotation_dict[clust] for clust in _adata.obs["cluster"]]
        else:
            _adata, annotation_dict = annotate_clusters_by_panglao(_adata, cluster_key='cluster', organism=organism)
    if database == 'CellTypist':
        # models.download_models(force_update = True)
        # models.models_path
        # models.models_description()
        if cell_type:
            annotation_dict = cell_type
            _adata.obs['cell_type'] = [annotation_dict[clust] for clust in _adata.obs["cluster"]]
        else:
            if organism == 'human':
                predictions = celltypist.annotate(_adata, model = 'Immune_All_High.pkl', majority_voting = True)
            if organism == 'mouse':   
                predictions = celltypist.annotate(_adata, model = 'Adult_Mouse_Gut.pkl', majority_voting = True)
            
            _adata = predictions.to_adata()
            _adata.obs['cell_type'] = _adata.obs['majority_voting']
            df = _adata.obs[["cluster", "cell_type"]]
            grouped = df.groupby('cluster')['cell_type'].apply(list).to_dict()
            annotation_dict = {k:'|'.join(set(v)) for k, v in grouped.items()}

    set_scanpy_color_palette(_adata, color_key="cell_type")
    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(_adata, color='cell_type', show=False, title='', legend_loc='on data', legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2)
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    groups = _adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["group"]==groups[i]], color="cell_type", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_group_celltype', fig_size={'width':len(groups)*5, 'height':5})

    result = _adata.obs.groupby(['cell_type', 'sample']).size().reset_index(name='number').pivot(index='cell_type', columns='sample', values='number').reset_index()
    result.columns.name = None
    # result['Total'] = result.sum(axis=1)
    # result['Percent'] = round((result['Total'] / sum(result['Total'])) * 100, 1)

    result.to_csv(os.path.join(out_dir, 'celltype_statistics.xls'), sep='\t', index=None)

    _adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')
    plot_celltype_proportion(_adata, out_dir=out_dir,type='group')
    plot_celltype_proportion(_adata, out_dir=out_dir,type='sample')
    # pd.DataFrame(_adata.obsm['X_tsne'], index=_adata.obs_names, columns=['tSNE_1', 'tSNE_2']).to_csv(os.path.join(out_dir, 'data_tsne.xls'), sep='\t')
    pd.DataFrame(_adata.obsm['X_umap'], index=_adata.obs_names, columns=['UMAP_1', 'UMAP_2']).to_csv(os.path.join(out_dir, 'data_umap.xls'), sep='\t')

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(_adata, color=["cluster","cell_type"], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    if organism == 'mouse':
        _adata.var.index = _adata.var['gene']

    _adata.write(os.path.join(out_dir, 'result.h5ad'))

    return _adata, annotation_dict


@st.fragment
def visua_gene_express(_adata, genes, out_dir):
    #st.set_option('deprecation.showPyplotGlobalUse', False)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color=genes, s=15, ncols=3, cmap='coolwarm', use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, '_'.join(genes)+'.png'), facecolor='white', bbox_inches="tight", dpi=200)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.dotplot(_adata, genes,groupby='cluster',cmap='coolwarm', use_raw=True, show=False,swap_axes=True)
        plt.savefig(os.path.join(out_dir, '_'.join(genes)+'_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=200)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.dotplot(_adata, genes,groupby='cell_type',cmap='coolwarm', use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, 'anno_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=200)
    gene_expression = Image.open(os.path.join(out_dir, '_'.join(genes)+'.png'))

    if len(genes) > 1:
        fig, axs = plt.subplots(len(genes), 1, figsize=(10, len(genes)*4), sharey=True)
        for i in range(len(genes)):
            sc.pl.violin(_adata, genes[i], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=axs[i], show=False)
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.violin(_adata, genes[0], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=ax, show=False)

    return gene_expression, fig

