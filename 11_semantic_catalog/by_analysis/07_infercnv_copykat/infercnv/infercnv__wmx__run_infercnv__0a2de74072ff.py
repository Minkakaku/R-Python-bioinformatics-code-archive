import os
import pandas as pd
import streamlit as st
import random
import subprocess
import scanpy as sc
import matplotlib.pyplot as plt
from PIL import Image
from zipfile import ZipFile
import streamlit_antd_components as sac
import streamlit.components.v1 as components

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def callback():
    st.session_state['tool'] = None

def image_view(_container, png_file, caption='', resize=None):
    image = Image.open(png_file)
    if resize:
        new_image = image.resize(resize)
        _container.image(new_image, caption=caption)
    else:
        _container.image(image, caption=caption, use_container_width=True)

# 生成结果报告
def make_report(out_dir):
    report_rmd = '/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/inferCNV/report.rmd'
    out_dir=os.path.join(out_dir,'CNV')
    os.system(f"""
        cp \
        {report_rmd} \
        {os.path.abspath(out_dir)} \
    """)
    os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(out_dir))
    # 删除中间文件
    os.system(f"""
              rm \
              {os.path.abspath(os.path.join(out_dir,'report.knit.md'))}
              rm \
              {os.path.abspath(os.path.join(out_dir,'report.rmd'))}
              rm \
              {os.path.abspath(os.path.join(st.session_state['inferCNV']['out_dir'],'adata.h5ad'))}
              """)

    # 去除¶
    file_path = os.path.join(out_dir,'report.html')

    # 读取HTML文件内容
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()

    html_content = html_content.replace('¶', '')
    # 如果是HTML实体 &para; 或 &182;
    html_content = html_content.replace('&para;', '').replace('&182;', '')

    # 写回文件
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(html_content)


def run_inferCNV_analysis():
    params = st.session_state['inferCNV']
    out_dir = params['out_dir']
    h5ad = os.path.abspath(os.path.join(out_dir,'adata.h5ad'))
    species = params['species']
    reference_celltype = ','.join(params['reference_celltype'])
    cmd = [
    "/mnt/NexData_file/share/scRNA-tools/envs/pip_infercnv/bin/python3.10",
    "/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/inferCNV/pip_infercnv.py",
    "--h5ad",h5ad,
    "--output", out_dir,
    "--species", species,
    "--reference_celltype", reference_celltype
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"inferCNV 分析失败:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

def inferCNV():
    if 'inferCNV' not in st.session_state:
        st.session_state['inferCNV'] = {}
        st.session_state['inferCNV']['out_dir'] = os.path.join(
            '/mnt/NexData_file/share/scRNA-tools/output', 'inferCNV',
            ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))
        )

    out_dir = st.session_state['inferCNV']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 inferCNV分析')
        st.markdown(f"""
            - **species:** 物种（human 或 mouse）；
            - **reference_celltype:**作为参考的细胞类型，如 Endothelial cells,T cells；
            - **结果位置：** `[{out_dir}]`
            """) 
        uploaded_file = st.file_uploader(label="RDS文件上传：", type='h5ad')
        reference_celltype = ""
        species=['human','mouse']
        adata=None
        adata_anno=None
        if uploaded_file:
            try:
                adata=sc.read_h5ad(uploaded_file)
                adata.write_h5ad(os.path.join(out_dir,'adata.h5ad'))
                meta=adata.obs
                reference_celltype=list(meta['cell_type'].unique())
            except Exception as e:
                st.error(f"读取 h5ad 文件失败: {str(e)}")
                st.stop()

    params, response = st.columns([2, 8])
    with params.form('inferCNV_run'):
        species = st.selectbox('物种:',['human','mouse'], index=0)       
        reference_celltype = st.multiselect('参考细胞类型 (可多选):',options=reference_celltype,default=reference_celltype)
        submitted = st.form_submit_button('🚀 开始分析', type='primary', use_container_width=True)       
    if uploaded_file and submitted:
        # st.session_state['inferCNV'] = {
        #     'out_dir': out_dir,
        #     'species': species,
        #     'reference_celltype': reference_celltype,
        #     'adata':adata,
        #     'adata_anno':adata_anno
        # }
        st.session_state['inferCNV']['species']=species
        st.session_state['inferCNV']['reference_celltype']=reference_celltype
        st.session_state['inferCNV']['adata']=adata
        st.session_state['inferCNV']['adata_anno']=adata_anno

        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            try:
                with st.spinner('inferCNV 分析中，请稍候...'):
                    run_inferCNV_analysis()
                cnv_cnvleiden_bycell_path=os.path.join(out_dir,'CNV', 'cnv_cnvleiden_bycell.png')
                cnv_cnvleiden_bycluster_path=os.path.join(out_dir,'CNV','cnv_cnvleiden_bycluster.png')
                cnv_vlnplot_path=os.path.join(out_dir,'CNV', 'cnv_vlnplot.png')

                if os.path.exists(cnv_cnvleiden_bycell_path):
                    st.success("✅ inferCNV分析完成,请注释！")

                tabs = st.tabs(["🍊 cnv_cnvleiden_bycell", "🍊 cnv_cnvleiden_bycluster", "🍊 cnv_vlnplot"])
                image_view(tabs[0], png_file=cnv_cnvleiden_bycell_path, caption='')
                image_view(tabs[1], png_file=cnv_cnvleiden_bycluster_path, caption='')
                image_view(tabs[2], png_file=cnv_vlnplot_path, caption='')

            except Exception as e:
                st.error(f"分析失败: {str(e)}")
                st.stop()

def inferCNV_anno():
    
    if 'adata' in st.session_state['inferCNV']:
        st.subheader('🍓 inferCNV注释')
        out_dir = st.session_state['inferCNV']['out_dir']
        adata=st.session_state['inferCNV']['adata']
        meta1=pd.read_csv(os.path.join(out_dir,'CNV','metadata.csv'),sep=',',header=0,index_col=0)
        adata.obs=meta1
        meta=adata.obs
        method = 'cnv_leiden'
        params, response = st.columns([2, 8])
        with params.form('annotation_inferCNV'):
            meta['inferCNV']='normal'
            
            annotation_dict = (
                meta[[method, 'inferCNV']]
                .drop_duplicates()
                .set_index(method)['inferCNV']
                .to_dict()
                )
            st.session_state['inferCNV']['annotation_dict']=annotation_dict
            df = pd.DataFrame(annotation_dict.items(), columns=['cnv_leiden', 'annotation'])          
            df['inferCNV'] = df['annotation']
            edited_df = st.data_editor(df, hide_index=True, disabled=('cnv_leiden','annotation'), use_container_width=True)   
            renames = edited_df.set_index('cnv_leiden')['inferCNV'].to_dict() 
            submitted1 = st.form_submit_button("🚀Rname", type='primary')  
        if submitted1:
            with response:
                try:
                    with st.spinner('Renaming...'):
                        adata_anno=adata.copy()
                        adata_anno.obs['inferCNV'] = [renames[clust] for clust in adata_anno.obs[method]]
                        adata_anno.obs[['cell_type','cnv_leiden','cnv_score_cell','cnv_score_cluster','inferCNV']].to_csv(os.path.join(out_dir,'CNV','meta_inferCNV.xls'), sep='\t')
                        umap_cnv_path=os.path.join(out_dir,'CNV','umap_cnv.png')
                        sc.pl.umap(adata_anno, color='inferCNV', show=False,cmap='coolwarm',size=25)
                        save_fig(plt,file_name='umap_cnv',out_dir=os.path.join(out_dir,'CNV'),fig_size={'w':7,'h':5})
                        make_report(out_dir)
                    if os.path.exists(umap_cnv_path):
                        st.success("✅ 分析完成！")
                    #tabs = st.tabs(["🍊 umap_cnv"])
                    #image_view(tabs[0], png_file=umap_cnv_path, caption='')
                    #make_report(out_dir)
                    report_path = os.path.join(out_dir,'CNV' ,"report.html") 
                    report = open(report_path)
                    components.html(report.read(), height=1000, scrolling=True)
                    zip_path = os.path.join(out_dir, "inferCNV_results.zip")
                    with ZipFile(zip_path, 'w') as zipf:
                        for root, _, files in os.walk(os.path.join(out_dir,'CNV')):
                            for file in files:
                                if file != "uploaded.rds" and not file.endswith(".log"):
                                    file_path = os.path.join(root, file)
                                    arcname = os.path.relpath(file_path, out_dir)
                                    zipf.write(file_path, arcname)

                    with open(zip_path, "rb") as f:
                        st.download_button("📥 下载完整结果", f, file_name="inferCNV_results.zip", mime="application/zip", use_container_width=True)

                except Exception as e:
                    st.error(f"分析失败: {str(e)}")
                    st.stop()  
