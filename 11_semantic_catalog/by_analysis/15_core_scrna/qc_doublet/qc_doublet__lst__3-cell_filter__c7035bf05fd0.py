import os
import numpy as np
from PIL import Image
import streamlit as st
from utils.utils import image_view
from utils.function import cell_filtering
from streamlit_extras.metric_cards import style_metric_cards


st.subheader("🍊 :rainbow[Cell Filtering]")
st.info('''&emsp;10X Genomics 官方软件Cell Ranger对测序原始数据进行统计后，使用barcode 序列标记区分细胞，通过每个细胞内不同mRNA 的UMI 标签对单细胞转录组数据进行定量。
        虽然单细胞转录组测序采用的是barcode 以及UMI 序列信息去标记细胞和转录本信息，但也有可能会出现带有这些标记信息的凝胶珠中没有细胞或者细胞数目不止一个的情况；
        并且当细胞发生死亡或者裂解时，细胞中将含有大量的线粒体基因。  
        **因此需要根据细胞的UMI总数、单个细胞鉴定到的基因数目及线粒体基因比例3个指标进行高质量细胞筛选。**
        这一步目的对数据进行初步筛选和过滤，一些测序质量差的数据会让后续的分析产生噪音和干扰，因此我们需要将它们去除。\n''')

out_dir = os.path.join(st.session_state['pipe']['out_dir'], '03-Cell_filter')
params, response, note = st.columns([2, 5, 3])

with params.form('cell_filter'):

    if st.session_state['pipe']['adata_raw']:
        upper_lim = np.quantile(st.session_state['pipe']['adata_qc'].obs.n_genes_by_counts.values, 0.95)
        st.session_state['pipe']['upper_lim'] = int(upper_lim)
        expected_doublet_rate = (st.session_state['pipe']['adata_qc'].X.shape[0] / 1000) * 0.05 * 10
        st.session_state['pipe']['expected_doublet_rate'] = int(expected_doublet_rate)
    else:
        st.session_state['pipe']['upper_lim'] = 5000
        st.session_state['pipe']['expected_doublet_rate'] = 5

    st.session_state['pipe']['min_genes'] = st.number_input(label='min_genes:', min_value=0, max_value=5000, value=st.session_state['pipe']['min_genes'], step=1)
    st.session_state['pipe']['min_cells'] = st.number_input(label='min_cells:', min_value=0, max_value=2000, value=st.session_state['pipe']['min_cells'], step=1)
    st.session_state['pipe']['pct_counts_mt'] = st.number_input(label='pct_counts_mt:', min_value=0, max_value=100, value=st.session_state['pipe']['pct_counts_mt'], step=1)
    st.session_state['pipe']['pct_counts_ribo'] = st.number_input(label='pct_counts_ribo:', min_value=0, max_value=100, value=st.session_state['pipe']['pct_counts_ribo'], step=1)
    st.session_state['pipe']['upper_lim'] = st.number_input(label='max_genes:', min_value=500, max_value=10000, value=st.session_state['pipe']['upper_lim'], step=10)
    st.session_state['pipe']['expected_doublet_rate'] = st.number_input(label='expected_doublet_rate(%):', min_value=0, max_value=100, value=st.session_state['pipe']['expected_doublet_rate'], step=1)    
    submitted = st.form_submit_button("🚀继续", type='primary')
    if submitted:
        if st.session_state['pipe']['adata_qc']:
            with response:
                with st.spinner('Calculating...'):
                    st.session_state['pipe']['adata_filtered'] = cell_filtering(st.session_state['pipe']['adata_qc'], 
                                                                        st.session_state['pipe']['min_genes'],
                                                                        st.session_state['pipe']['min_cells'],
                                                                        st.session_state['pipe']['pct_counts_mt'],
                                                                        st.session_state['pipe']['pct_counts_ribo'],
                                                                        st.session_state['pipe']['upper_lim'],
                                                                        st.session_state['pipe']['expected_doublet_rate'],
                                                                        out_dir)
        else:
            response.error('请检查上一步Data Loading是否正确执行！', icon="🚨")

if st.session_state['pipe']['adata_filtered']:
    with response.expander("**:blue[结果解读📜]** ", expanded=True):
        col1, col2, col3, col4 = st.columns(4)
        style_metric_cards(border_left_color="#006400")
        col1.metric("原始细胞数", st.session_state['pipe']['adata_raw'].X.shape[0])
        col2.metric("原始基因数", st.session_state['pipe']['adata_raw'].X.shape[1])
        col3.metric("过滤后细胞数", st.session_state['pipe']['adata_filtered'].X.shape[0], st.session_state['pipe']['adata_filtered'].X.shape[0]-st.session_state['pipe']['adata_raw'].X.shape[0])
        col4.metric("过滤后基因数", st.session_state['pipe']['adata_filtered'].X.shape[1], st.session_state['pipe']['adata_filtered'].X.shape[1]-st.session_state['pipe']['adata_raw'].X.shape[1])
        tabs = st.tabs(["🚀violin plot (filtered)", "📈scatter plot (filtered)"])
        img_filter_violin = Image.open(os.path.join(out_dir, 'filtered_nCount_nFeature_voilon.png'))
        tabs[0].image(img_filter_violin, caption='注：nGene图：过滤后各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤后各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤后各个样本单个细胞中线粒体基因表达量的百分比分布情况；percent.ribo图：过滤后各个样本单个细胞中核糖体基因表达量的百分比分布情况。', use_column_width=True)
        img_filter_scatter = Image.open(os.path.join(out_dir, 'filtered_nCount_nFeature_scatter.png'))
        tabs[1].image(img_filter_scatter, caption='注：nGene图：过滤后各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤后各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤后各个样本单个细胞中线粒体基因表达量的百分比分布情况；percent.ribo图：过滤后各个样本单个细胞中核糖体基因表达量的百分比分布情况。', use_column_width=True)
    with note.expander(f"**:blue[参数说明📜]** ", expanded=True):
        st.markdown("""
                    - **min_genes**: 去除表达基因数小于n的细胞(n默认值200);
                    - **min_cells**: 去除低于n个细胞表达的基因(n默认值10);
                    - **pct_counts_mt**: 过滤线粒体基因占比大于20%的细胞;
                    - **pct_counts_ribo**: 过滤核糖体基因占比大于30%的细胞;
                    - **max_genes**: 过滤基因总数大于n的细胞(n默认值为细胞总数的95%分位数);
                    - **expected_doublet_rate(%)**: 期望双细胞或多细胞占总细胞数百分比(默认值以10x技术平台官方给出的1000个细胞多胞率进行计算);
                    """)
else:
    with response.expander("**:blue[结果解读📜]** ", expanded=True):
        st.markdown("""**分析结果文件位置：03-Cell_filter/**""")
    with note.expander("**:blue[参数说明📜]** ", expanded=True):
        st.markdown("""
                    - **min_genes**: 去除表达基因数小于n的细胞(n默认值200);
                    - **min_cells**: 去除低于n个细胞表达的基因(n默认值10);
                    - **pct_counts_mt**: 过滤线粒体基因占比大于20%的细胞;
                    - **pct_counts_ribo**: 过滤核糖体基因占比大于30%的细胞;
                    - **max_genes**: 过滤基因总数大于n的细胞(n默认值为细胞总数的95%分位数);
                    - **expected_doublet_rate(%)**: 期望双细胞或多细胞占总细胞数百分比(默认值以10x技术平台官方给出的1000个细胞多胞率进行计算);
                    """)
