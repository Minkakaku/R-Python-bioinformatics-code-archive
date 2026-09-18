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

def run_cellchat_analysis():
    params = st.session_state['cellchat_params']

    cmd = [
        "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/bin/Rscript",
        "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R",
        "--rds", params['rds_path'],
        "--outpath", os.path.join(params['out_dir'],'cellchat'),
        "--splitby", params['splitby'],
        "--species", params['species']
    ]

    if params['comparison'] != 'F':
        cmd += ["--comparison", params['comparison']]
    if params['subcelltype']:
        cmd += ["--subcelltype", ','.join(params['subcelltype'])]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"CellChat 分析失败:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

def cellchat():
    if 'cellchat' not in st.session_state:
        st.session_state['cellchat'] = {}
        st.session_state['cellchat']['out_dir'] = os.path.join(
            '/mnt/NexData_file/share/scRNA-tools/output', 'cellchat',
            ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))
        )

    out_dir = st.session_state['cellchat']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 细胞通讯分析')
        st.write('使用cellchat对单细胞数据进行细胞通讯分析。')
        st.markdown(f"""
            - **splitby:** 用于划分组别的列，如 `sample` 或 `group`；
            - **comparison:** 比较组名，格式为 `GroupA:GroupB`，默认 `F` 为单组分析；
            - **subcelltype:** 指定感兴趣的细胞类型子集进行分析；
            - **species:** 物种（human 或 mouse）；
            - **结果位置：** `[{out_dir}]`
            """)
        uploaded_file = st.file_uploader(label="RDS文件上传：", type='rds')
        comparison_options = ['F']
        group_options = []
        cell_type_options = []

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

                group_options = ['group','sample']

                for col in ['cell_type', 'CellType', 'celltype']:
                    if col in meta_data.columns:
                        cell_type_options = list(meta_data[col].unique())
                        break


            except Exception as e:
                st.error(f"读取 RDS 文件失败: {str(e)}")
                st.stop()
        splitby = st.selectbox('分组依据 (splitby):', group_options, index=None)
    params, response = st.columns([2, 8])
    with params.form('cellchat_run'):
        species = st.selectbox('物种 (species):', ['human', 'mouse'], index=0)
        if splitby:
            comparison_group = splitby
            group_values = list(meta_data[comparison_group].unique())
            from itertools import combinations
            comparison_options = ['F'] + [f"{a}:{b}" for a, b in combinations(group_values, 2)]        
            comparison = st.selectbox('比较组 (comparison):', comparison_options, index=0)
        subcelltype = st.multiselect('分析细胞类型 (subcelltype):', options=cell_type_options, default=cell_type_options)
        submitted = st.form_submit_button('🚀 开始分析', type='primary', use_container_width=True)
            
    if uploaded_file and submitted:
        st.session_state['cellchat_params'] = {
            'rds_path': rds_path,
            'out_dir': out_dir,
            'splitby': splitby,
            'comparison': comparison,
            'subcelltype': subcelltype,
            'species': species
        }

        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            try:
                with st.spinner('CellChat 分析中，请稍候...'):
                    run_cellchat_analysis()

                report_path = os.path.join(out_dir,'cellchat' ,"pip_report_cellchat_single.html") if comparison == 'F' else \
                    os.path.join(out_dir, 'cellchat' ,comparison.replace(":", "_"), "pip_report_cellchat.nb.html")
                if not os.path.exists(report_path):
                    report_path = os.path.join(out_dir,'cellchat',"pip_report_cellchat.html")

                if os.path.exists(report_path):
                    st.success("✅ 分析完成！")
                    report = open(report_path)
                    components.html(report.read(), height=1000, scrolling=True)   

                else:
                    st.warning("分析完成，但未找到报告文件。")

                zip_path = os.path.join(out_dir, "cellchat_results.zip")
                with ZipFile(zip_path, 'w') as zipf:
                    for root, _, files in os.walk(os.path.join(out_dir,'cellchat')):
                        for file in files:
                            if file != "uploaded.rds" and not file.endswith(".log"):
                                file_path = os.path.join(root, file)
                                arcname = os.path.relpath(file_path, out_dir)
                                zipf.write(file_path, arcname)

                with open(zip_path, "rb") as f:
                    st.download_button("📥 下载完整结果", f, file_name="cellchat_results.zip", mime="application/zip", use_container_width=True)

            except Exception as e:
                st.error(f"分析失败: {str(e)}")
                st.stop()


