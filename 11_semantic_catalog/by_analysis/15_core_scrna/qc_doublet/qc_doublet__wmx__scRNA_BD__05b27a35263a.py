import re
import os
os.environ['OPENBLAS_NUM_THREADS'] = "64"
os.environ['MKL_NUM_THREADS'] = "64"
os.environ['OMP_NUM_THREADS'] = "64"
import shutil
import pandas as pd
from PIL import Image
import streamlit as st
import streamlit_antd_components as sac
import streamlit.components.v1 as components
from moudle.common import *
from config import out_dir, data_dir, marker_dict_file, rmd_template_BD, rmd_template_img_BD
from moudle.BD_scanpy_pipeline import  get_BD_mtx,make_report_dir,extract_metrics,unzip_mex
import zipfile
import csv

horizontal_bar = "<hr style='margin-top: 2; margin-bottom: 0; height: 1px; border: 1px solid #635985;'><br>"
hide_table_row_index = """
                        <style>
                            thead tr th:first-child {display:none}
                            tbody th {display:none}
                        </style>
                        """
from collections import OrderedDict
def read_sample_dict(sample_file):
    df = pd.read_csv(sample_file, sep="\t")
    sample_dict = OrderedDict()
    for _, row in df.iterrows():
        sample_dict.setdefault(row["group"], []).append(row["sample"])
    return sample_dict

def data_server(pid, cols):
    tmp = []
    input_dir = os.path.join(data_dir, pid)
    sample_file = os.path.join(input_dir, 'sample.txt')
    st.session_state['out_dir'] = os.path.join(out_dir, pid)
    metrics_summary = os.path.join(st.session_state['out_dir'], '01-Rhapsody', 'metrics_summary.csv')
    
    if not os.path.exists(sample_file):
        tmp = []
        for sample in os.listdir(input_dir):
            if os.path.isdir(os.path.join(input_dir, sample)):
                tmp.append({"sample": sample, "rename": sample, "group": sample})
                unzip_mex(input_dir, sample)
        sample_df = cols[0].data_editor(pd.DataFrame(tmp), hide_index=True, disabled=('sample',), use_container_width=False)             
        #sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
        sample_dict=read_sample_dict(sample_file)
        sample_df.to_csv(os.path.join(input_dir, 'sample.txt'), sep='\t', index=None)
    else:
        df = pd.read_csv(sample_file, sep='\t')
        df = df.reset_index()
        sample_df = cols[0].data_editor(df, hide_index=True, disabled=('sample',), use_container_width=True)             
        sample_dict=read_sample_dict(sample_file)
        sample_df.to_csv(os.path.join(input_dir, 'sample.txt'), sep='\t', index=None)

    submit = cols[0].button('🚀Submit', type='primary')
    
    if submit:
        with cols[1]:
            with st.spinner('数据加载中...'):
                make_report_dir(st.session_state['out_dir'])
                sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
                adata = get_BD_mtx(data_path=input_dir, sample_dict=sample_dict, out_dir=st.session_state['out_dir'])
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
            metrics_summary = os.path.join(st.session_state['out_dir'], '01-Rhapsody', 'metrics_summary.csv')
            sample_list = list(adata.obs['sample'].drop_duplicates())
            try:
                metrics_sample = []
                for s in sample_list:
                    #在原样本文件夹下输出Sample_Metrics_Summary_need.csv
                    input_sample = os.path.join(input_dir, s, s+'_Metrics_Summary.csv')
                    output_sample = os.path.join(input_dir, s, s+'_Metrics_Summary_need.csv')
                    extract_metrics(
                        input_file=input_sample,
                        output_file=output_sample
                        )
                    #合并输出到01-Rhapsody/metrics_summary.csv
                    metrics = pd.read_csv(os.path.join(input_dir, s, s+'_Metrics_Summary_need.csv'),sep=',')
                    metrics['Sample Name'] = s
                    cols = list(metrics)
                    #cols = [cols[-1]] + cols[:-1]
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
    
      

 
def main():
    st.subheader('🍄 BD scRNA标准分析')
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
        st.info('根据项目编号，加载/mnt/NexData_file/nas/outputdata/cellranger_out路径下Rhapsody输出结果文件，选择项目编号后会根据样本名称生成样本信息表，用户可以在样本信息表中修改样本分组信息，用于后续分析', icon="📢")
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
            generate_report(rmd_template_BD, rmd_template_img_BD)


if __name__ == '__main__':
    main()
