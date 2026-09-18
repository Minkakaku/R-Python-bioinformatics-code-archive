import os
import pandas as pd
import streamlit as st
import random
import subprocess
from zipfile import ZipFile
import streamlit_antd_components as sac
import streamlit.components.v1 as components

def callback():
    st.session_state['tool'] = None

def run_monocle2_analysis():
    params = st.session_state['monocle2']

    cmd = [
        "/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle2/bin/Rscript",
        "/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/pseudotime/monocle2.analysis.R",
        "--rds", params['rds_path'],
        "--output", os.path.join(params['out_dir'],'monocle2'),
    ]

    # 可选参数，存在则加入
    if params.get('type'):
        cmd += ['--type', params['type']]
    if params.get('downcluster'):
        cmd += ['--downcluster']
    if params.get('maxcells'):
        cmd += ['--maxcells', str(params['maxcells'])]
    if params.get('cluster'):
        cmd += ['--cluster', params['cluster']]
    if params.get('ptsize'):
        cmd += ['--ptsize', str(params['ptsize'])]
    if params.get('BEAM'):
        cmd += ['--BEAM']

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"monocle2 分析失败:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

def monocle2():
    if 'monocle2' not in st.session_state:
        st.session_state['monocle2'] = {}
        st.session_state['monocle2']['out_dir'] = os.path.join(
            '/mnt/NexData_file/share/scRNA-tools/output', 'monocle2',
            ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))
        )

    out_dir = st.session_state['monocle2']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 单细胞拟时序分析')
        st.write('使用monocle2对单细胞数据进行拟时序分析。')
        st.markdown(f"""
            - **type**：列名，该列用于可视化轨迹（默认：cell_type）；
            - **downcluster**：是否启用子聚类（down-cluster）模式，若启用将对数据进行进一步子聚类分析；
            - **maxcells**：若启用 downcluster，每个 cluster 最多抽取的细胞数量（默认 1000）；
            - **cluster**：聚类的列名，用于指定 trajectory 起始点，默认使用 leiden；
            - **ptsize**：图中点的大小（默认 1.0），可根据图形密集程度调节；
            - **BEAM**：是否启用 BEAM（Branch Expression Analysis Modeling）分析，用于识别轨迹分支上的基因表达变化（默认 FALSE）；
            - **结果位置：** `[{out_dir}]`
            """)
        uploaded_file = st.file_uploader(label="RDS文件上传：", type='rds')
        type_use_DE = None
        downcluster = False
        maxcells = 1000
        cluster_DE = None
        ptsize = 1
        beam_flag = False
        col_use=[]
        if uploaded_file:
            try:
                rds_path = os.path.join(out_dir, "uploaded.rds")
                with open(rds_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
                import rpy2
                from rpy2.robjects import r
                from rpy2.robjects import pandas2ri
                from rpy2.robjects.conversion import localconverter
                from rpy2.robjects import default_converter

                r(f'seurat_obj <- readRDS("{rds_path}")')

                with localconverter(default_converter + pandas2ri.converter):
                    meta_data = r("as.data.frame(seurat_obj@meta.data)")
                
                col_use = list(meta_data.columns.unique())

            except Exception as e:
                st.error(f"读取 RDS 文件失败: {str(e)}")
                st.stop()
    params, response = st.columns([2, 8])
    with params.form('monocle2_run'):
        species = st.selectbox('物种 (species):', ['human', 'mouse'], index=0)
        type_use= st.selectbox('可视化列名:', col_use, index=None)
        downcluster = st.checkbox("开启子聚类")
        maxcells = st.number_input("最大细胞数", min_value=100,value=1000, step=100)
        cluster = st.selectbox('聚类列名:', col_use, index=None)
        ptsize_val = st.number_input("点大小", min_value=0.1, value=1.0,step=0.1)
        beam_flag = st.checkbox("启用BEAM分析")
        submitted = st.form_submit_button('🚀 开始分析', type='primary', use_container_width=True)
            
    if uploaded_file and submitted:
        st.session_state['monocle2'].update({
                        'rds_path': rds_path,             # 必选
                        'type': type_use,                 # 可选
                        'downcluster': downcluster,  # 可选，布尔值
                        'maxcells': maxcells,        # 可选，数字
                        'cluster': cluster,          # 可选，列名
                        'ptsize': ptsize,            # 可选，数字
                        'BEAM': beam_flag                # 可选，布尔值
                    })

        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            try:
                with st.spinner('monocle2 分析中，请稍候...'):
                    run_monocle2_analysis()

                report_path = os.path.join(out_dir,'monocle2', "report.html") 

                if os.path.exists(report_path):
                    st.success("✅ 分析完成！")
                    report = open(report_path)
                    components.html(report.read(), height=1000, scrolling=True)   

                else:
                    st.warning("分析完成，但未找到报告文件。")

                zip_path = os.path.join(out_dir, "monocle2_results.zip")
                with ZipFile(zip_path, 'w') as zipf:
                    for root, _, files in os.walk(os.path.join(out_dir,'monocle2')):
                        for file in files:
                            if file != "uploaded.rds" and not file.endswith(".log"):
                                file_path = os.path.join(root, file)
                                arcname = os.path.relpath(file_path, out_dir)
                                zipf.write(file_path, arcname)

                with open(zip_path, "rb") as f:
                    st.download_button("📥 下载完整结果", f, file_name="monocle2_results.zip", mime="application/zip", use_container_width=True)

            except Exception as e:
                st.error(f"分析失败: {str(e)}")
                st.stop()


