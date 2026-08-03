import os
import pandas as pd
import streamlit as st
import random
import subprocess
import scanpy as sc
import matplotlib.pyplot as plt
from zipfile import ZipFile
import streamlit_antd_components as sac
import streamlit.components.v1 as components

def callback():
    st.session_state['tool'] = None

def run_SCENIC_analysis():
    params = st.session_state['SCENIC']
    h5ad = params['h5ad_path']
    outdir = os.path.join(params['out_dir'],'SCENIC')
    os.makedirs(outdir,exist_ok=True)
    species = params['species']
    sample_fraction = params['sample_fraction']

    cmd = [
    "python",
    "/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/scenic.py",
    "--adata_path", h5ad,
    "--output_dir", outdir,
    "--species",species,
    "--sample_fraction", str(sample_fraction),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"CellChat 分析失败:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

def SCENIC():
    if 'SCENIC' not in st.session_state:
        st.session_state['SCENIC'] = {}
        st.session_state['SCENIC']['out_dir'] = os.path.join(
            '/mnt/NexData_file/share/scRNA-tools/output', 'SCENIC',
            ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))
        )

    out_dir = st.session_state['SCENIC']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    #_, context, _ = st.columns([0.5, 8, 1])
    col1, context, col3=st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 转录因子分析')
        st.markdown(f"""
            - **species:** 物种（human 或 mouse）；
            - **sample_fraction:** 抽样比例，默认0.3；
            - **结果位置：** `[{out_dir}]`
            """)
        uploaded_file = st.file_uploader(label="h5ad文件上传：", type='h5ad')
        species=['human','mouse']
        sample_fraction = 0.3
        adata=None
        h5ad_path=''
        if uploaded_file:
            try:
                adata=sc.read_h5ad(uploaded_file)
                adata.write_h5ad(os.path.join(out_dir,'adata.h5ad'))
                h5ad_path=os.path.join(out_dir,'adata.h5ad')
            except Exception as e:
                st.error(f"读取 h5ad 文件失败: {str(e)}")
                st.stop()

    col_params, col_response = st.columns([2, 8])
    with col_params.form('SCENIC_run'):
        species = st.selectbox('物种:', species, index=0)       
        sample_fraction = st.number_input("抽样比例", min_value=0.1, value=0.3,step=0.1)
        submitted = st.form_submit_button('🚀 开始分析', type='primary', use_container_width=True)
            
    if uploaded_file and submitted:
        st.session_state['SCENIC'].update({
                'h5ad_path': h5ad_path,
                'species': species, 
                'sample_fraction': sample_fraction
            })

        with col_response.expander("**:blue[结果解读📜]** ", expanded=True):
            try:
                with st.spinner('SCENIC 分析中，请稍候...'):
                    run_SCENIC_analysis()

                report_path = os.path.join(out_dir, "report.html") 
                if os.path.exists(report_path):
                    st.success("✅ 分析完成！")
                    report = open(report_path)
                    components.html(report.read(), height=1000, scrolling=True)   

                else:
                    st.warning("分析完成，但未找到报告文件。")

                zip_path = os.path.join(out_dir, "SCENIC_results.zip")
                with ZipFile(zip_path, 'w') as zipf:
                    for root, _, files in os.walk(os.path.join(params['out_dir'],'SCENIC')):
                        for file in files:
                            if file != "uploaded.rds" and not file.endswith(".log"):
                                file_path = os.path.join(root, file)
                                arcname = os.path.relpath(file_path, out_dir)
                                zipf.write(file_path, arcname)

                with open(zip_path, "rb") as f:
                    st.download_button("📥 下载完整结果", f, file_name="SCENIC_results.zip", mime="application/zip", use_container_width=True)

            except Exception as e:
                st.error(f"分析失败: {str(e)}")
                st.stop()


