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

def run_scscore_analysis():
    params = st.session_state['scscore']
    h5ad = params['h5ad_path']
    outdir = os.path.join(params['out_dir'],'sc_score')
    os.makedirs(outdir,exist_ok=True)
    gmt = params['gmt_path']

    cmd = [
    "python",
    "/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/sc_score/sc_score.py",
    "--h5ad", h5ad,
    "--output", outdir,
    "--gmt",gmt
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Scscore 分析失败:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

def scscore():
    if 'scscore' not in st.session_state:
        st.session_state['scscore'] = {}
        st.session_state['scscore']['out_dir'] = os.path.join(
            '/mnt/NexData_file/share/scRNA-tools/output', 'scscore',
            ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))
        )

    out_dir = st.session_state['scscore']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 基因集评分分析')
        st.markdown(f"""
            - **结果位置：** `[{out_dir}]`
            """)
        uploaded_file1 = st.file_uploader(label="h5ad文件上传：", type='h5ad')
        uploaded_file2 = st.file_uploader(label="gmt文件上传：", type='gmt')
        gmt=None
        adata=None
        h5ad_path=''
        gmt_path=''
        if uploaded_file1:
            try:
                adata=sc.read_h5ad(uploaded_file1)
                adata.write_h5ad(os.path.join(out_dir,'adata.h5ad'))
                h5ad_path=os.path.join(out_dir,'adata.h5ad')

            except Exception as e:
                st.error(f"读取 h5ad 文件失败: {str(e)}")
                st.stop()
        if uploaded_file2:
            try:
                gmt_path=os.path.join(out_dir,'gene_set.gmt')
                gmt=pd.read_csv(uploaded_file2, sep="\t", header=None, index_col=0)
                gmt.to_csv(gmt_path,sep='\t',header=None,index=True)

            except Exception as e:
                st.error(f"读取 gmt 文件失败: {str(e)}")
                st.stop()

    params, response = st.columns([2, 8])
    with params.form('scscore_run'):
        submitted = st.form_submit_button('🚀 开始分析', type='primary', use_container_width=True)
            
    if uploaded_file1 and submitted:
        st.session_state['scscore'].update({
                'h5ad_path': h5ad_path,
                'gmt_path': gmt_path
            })

        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            try:
                with st.spinner('scscore 分析中，请稍候...'):
                    run_scscore_analysis()

                report_path = os.path.join(out_dir,"sc_score","report.html") 
                if os.path.exists(report_path):
                    st.success("✅ 分析完成！")
                    report = open(report_path)
                    components.html(report.read(), height=1000, scrolling=True)   

                else:
                    st.warning("分析完成，但未找到报告文件。")

                zip_path = os.path.join(out_dir, "scscore_results.zip")
                with ZipFile(zip_path, 'w') as zipf:
                    for root, _, files in os.walk(os.path.join(out_dir,'sc_score')):
                        for file in files:
                            if file != "uploaded.rds" and not file.endswith(".log"):
                                file_path = os.path.join(root, file)
                                arcname = os.path.relpath(file_path, out_dir)
                                zipf.write(file_path, arcname)

                with open(zip_path, "rb") as f:
                    st.download_button("📥 下载完整结果", f, file_name="scscore_results.zip", mime="application/zip", use_container_width=True)

            except Exception as e:
                st.error(f"分析失败: {str(e)}")
                st.stop()


