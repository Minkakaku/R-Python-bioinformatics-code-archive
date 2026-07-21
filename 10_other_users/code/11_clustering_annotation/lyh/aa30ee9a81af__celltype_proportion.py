import os
import random
import scanpy as sc
import streamlit as st
from utils.utils import *
from utils.subcluster import image_view


def celltype_proportion():
    if 'cp' not in st.session_state:
        st.session_state['cp'] = {}
        st.session_state['cp']['out_dir'] = os.path.join('tmp/tools', 'cell_proprtion', ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8)))
    out_dir = st.session_state['cp']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader("🍓 细胞比例计算及可视化")
        st.write("细胞分群命名完成之后，大多数研究会关注各个样本细胞群的比例，尤其是免疫学的研究中，各群免疫细胞比例的变化是比较受关注的。 ")
        uploaded_file = st.file_uploader(label="h5ad文件上传:", type='h5ad')
    if uploaded_file:
        adata = sc.read_h5ad(uploaded_file)
        context.write(adata)
        columns = adata.obs.columns
        if 'cell_type' in columns:
            types = ['cell_type']
            adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')
        else:
            types = columns
        params, response, note = context.columns([2, 5, 3])
        with params.form('cell_proportion'):
            option = st.selectbox('Select cell_type or cluster:', types)
            x_axis = st.selectbox('x-axis:', ("sample", "group"))
            fig_width = st.number_input(label='figure width:', min_value=4.0, max_value=12.0, value=4.0, step=0.5)
            fig_heigh = st.number_input(label='figure heigh:', min_value=3.0, max_value=10.0, value=3.5, step=0.5)
            submitted = st.form_submit_button("🚀Submit", type='primary')
        if submitted:
            os.system(f"Rscript ./utils/Rscript/celltype_proportion.R {os.path.join(out_dir, 'all_cells_features.xls')} {x_axis} {out_dir} {fig_width} {fig_heigh}")    
        if os.path.exists(os.path.join(out_dir, 'celltype_proportion.png')):
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                tabs = st.tabs(["🍓 Cell Proportion"])        
                image_view(tabs[0], png_file=os.path.join(out_dir, 'celltype_proportion.png'))
            with note.expander("**:blue[参数说明📜]** ", expanded=True):
                st.markdown(f"""
                    - **Select cell_type or cluster:** 上传的h5ad文件中，adata.obs中必须包含cell_type列(默认按样本和cell_type进行细胞比例统计);
                    - **x-axis:** x轴展示数据；
                    - **figure width::** 设置图片宽度；
                    - **figure heigh:** 设置图片高度；
                    - **结果文件位置：:rainbow[{st.session_state['cp']['out_dir']}]**
                    """)
            with open(os.path.join(out_dir, 'celltype_proportion.png'), "rb") as file:
                btn = note.download_button(
                                            label="Download",
                                            data=file,
                                            file_name="celltype_proportion.png",
                                            mime="image/png", type='primary',
                                            help='Download cell proportion of samples.'
                                          )
        else:
            with response.expander("**:blue[结果解读📜]** ", expanded=True):
                st.markdown("""👈**细胞类型比例计算，选择分组信息，设置图片宽、高。**""") 
            with note.expander("**:blue[参数说明📜]** ", expanded=True):
                st.markdown(f"""
                    - **Select cell_type or cluster:** 上传的h5ad文件中，adata.obs中必须包含cell_type列(默认按样本和cell_type进行细胞比例统计);
                    - **x-axis:** x轴展示数据；
                    - **figure width::** 设置图片宽度；
                    - **figure heigh:** 设置图片高度；
                    """)
    else:
        context.info('上传.h5ad文件(adata.obs中必须包含cell_type列或cluster列), 计算细胞比例及柱状图进行可视化。', icon="ℹ️")
