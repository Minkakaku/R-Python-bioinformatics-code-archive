import os
import pandas as pd
import streamlit as st
import random
import subprocess
import scanpy as sc
import matplotlib.pyplot as plt
import sys
sys.path.append('/mnt/NexData_file/share/scRNA-tools/Proption/pyscript')
from draw_group import draw_group
from zipfile import ZipFile
import streamlit_antd_components as sac
import streamlit.components.v1 as components

WorkDir = '/mnt/NexData_file/share/scRNA-tools/Proption'
def callback():
    st.session_state['tool'] = None

def get_color(adata,method,outdir):
    keys = adata.obs[method].cat.categories.to_list()
    values = list(adata.uns[method+'_colors'])
    color_dict = dict(zip(keys,values))
    df = pd.DataFrame.from_dict(color_dict,orient='index',columns=[method])
    df.to_csv(os.path.join(outdir,method+'_color.xls'),sep="\t",index=True,index_label=method)
def run_Proption_analysis():
    params = st.session_state['Proption_params']
    adata = params['adata']
    outdir = params['out_dir']
    method = params['celltype']
    group = params['group']
    draw_group(adata,out_dir = outdir,method = method)
    get_color(adata,method=method,outdir=outdir)
    get_color(adata,method=group,outdir=outdir)
    adata.obs.to_csv(os.path.join(outdir,"meta.xls"),sep ="\t")
    # 2.进行分析

    cmd = [
    "Rscript",
    "/mnt/NexData_file/share/scRNA-tools/Proption/Rscript/cell_proportion_v1.R",
    "--MetaData", os.path.abspath(os.path.join(outdir, 'meta.xls')),
    "--Celltype", params['celltype'],
    "--Group", params['group'],
    "--Output", os.path.join(params['out_dir'],'Proption'),
    "--CelltypeColor", os.path.abspath(os.path.join(outdir,'Proption', method+'_color.xls')),
    "--GroupColor", os.path.abspath(os.path.join(outdir,'Proption', group+'_color.xls'))
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"CellChat 分析失败:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")


    # shell_script = '/mnt/NexData_file/share/scRNA-tools/Proption/2.pip.sh' # 替换为实际的shell脚本路径
    # os.system(f"""
    #     sh \
    #     {shell_script} \
    #     {os.path.abspath(os.path.join(outdir, 'meta.xls'))} \
    #     {method} \
    #     {group} \
    #     {os.path.abspath(outdir)} \
    #     {os.path.abspath(os.path.join(outdir, method+'_color.xls'))} \
    #     {os.path.abspath(os.path.join(outdir, group+'_color.xls'))}
    # """)

    # 生成结果报告
    report_rmd = '/mnt/NexData_file/share/scRNA-tools/Proption/report/report.rmd'
    os.system(f"""
        cp \
        {report_rmd} \
        {os.path.abspath(os.path.join(outdir,'Proption'))} \
    """)
    # shell_script = '/mnt/NexData_file/share/scRNA-tools/Proption/3.report.sh'
    # os.system(f"""
    #     Rscript -e rmarkdown::render('$OUTPUT_DIR/report.rmd',output_file='$OUTPUT_DIR/report.html')\
    #     {os.path.abspath(os.path.join(outdir))}
    # """)
    os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(os.path.join(outdir,'Proption')))
    # 删除中间文件
    os.system(f"""
              rm \
              {os.path.abspath(os.path.join(os.path.join(outdir,'Proption'),'report.knit.md'))}
              rm \
              {os.path.abspath(os.path.join(os.path.join(outdir,'Proption'),'report.rmd'))}
              """)

    # 去除¶
    file_path = os.path.join(outdir,'Proption','report.html')

    # 读取HTML文件内容
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()

    html_content = html_content.replace('¶', '')
    # 如果是HTML实体 &para; 或 &182;
    html_content = html_content.replace('&para;', '').replace('&182;', '')

    # 写回文件
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(html_content)

def Proption():
    if 'Proption' not in st.session_state:
        st.session_state['Proption'] = {}
        st.session_state['Proption']['out_dir'] = os.path.join(
            '/mnt/NexData_file/share/scRNA-tools/output', 'Proption',
            ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))
        )

    out_dir = st.session_state['Proption']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 细胞比例分析')
        st.markdown(f"""
            - **group_options:** 用于划分组别的列，如 `sample` 或 `group`；
            - **cell_type_options:** 细胞类型注释列名，如 cell_type；
            - **结果位置：** `[{out_dir}]`
            """)
        uploaded_file = st.file_uploader(label="h5ad文件上传：", type='h5ad')
        group=['sample','group']
        cell_type_options = []
        celltype=['cell_type']
        adata=None
        if uploaded_file:
            try:
                adata=sc.read_h5ad(uploaded_file)
                meta=adata.obs
                cell_type_options=list(meta.columns)
            except Exception as e:
                st.error(f"读取 h5ad 文件失败: {str(e)}")
                st.stop()

    params, response = st.columns([2, 8])
    with params.form('Proption_run'):
        group = st.selectbox('分组:', [ 'group','sample'], index=0)       
        celltype = st.selectbox('细胞类型:',cell_type_options, index=0)
        submitted = st.form_submit_button('🚀 开始分析', type='primary', use_container_width=True)
            
    if uploaded_file and submitted:
        st.session_state['Proption_params'] = {
            'out_dir': out_dir,
            'group': group,
            'celltype': celltype,
            'adata':adata
        }

        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            try:
                with st.spinner('Proption 分析中，请稍候...'):
                    run_Proption_analysis()

                report_path = os.path.join(out_dir,'Proption',"report.html") 
                if os.path.exists(report_path):
                    st.success("✅ 分析完成！")
                    report = open(report_path)
                    components.html(report.read(), height=1000, scrolling=True)   

                else:
                    st.warning("分析完成，但未找到报告文件。")

                zip_path = os.path.join(out_dir, "Proption_results.zip")
                with ZipFile(zip_path, 'w') as zipf:
                    for root, _, files in os.walk(os.path.join(out_dir,'Proption')):
                        for file in files:
                            if file != "uploaded.rds" and not file.endswith(".log"):
                                file_path = os.path.join(root, file)
                                arcname = os.path.relpath(file_path, out_dir)
                                zipf.write(file_path, arcname)

                with open(zip_path, "rb") as f:
                    st.download_button("📥 下载完整结果", f, file_name="Proption_results.zip", mime="application/zip", use_container_width=True)

            except Exception as e:
                st.error(f"分析失败: {str(e)}")
                st.stop()


