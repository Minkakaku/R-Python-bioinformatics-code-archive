import os
import pandas as pd
import scanpy as sc
import streamlit as st
import numpy as np
import gseapy as gp
import random
from matplotlib.pyplot import rc_context
from gseapy.plot import barplot
import matplotlib.pyplot as plt
from moudle.common import image_view
from zipfile import ZipFile
from pathlib import Path

def callback():
    st.session_state['tool'] = None


def save_fig(plt, out_dir, fig_size={'width':5, 'height':5}):
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, 'dotplot.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, 'dotplot.pdf'))
def sc_dotplot():
    st.session_state['pipe'] = {}
    data = pd.read_csv('/home/lst/functionenrichment/config/annotation/marker_dict.txt', sep='\t')
    result = {}
    dot_fig=None
    for index, row in data.iterrows():
        cell_type = row['cell_type']
        organism = row['organism']
        genes = row['markers'].split(',')
        if organism not in result:
            result[organism] = {}
        if cell_type not in result[organism]:
            result[organism][cell_type] = []  
        result[organism][cell_type].extend(genes)
    st.session_state['pipe']['marker_dict'] = result
    st.session_state['pipe']['gene_expression'] = ''
    select_list=[]
    adata=None
    if 'sc_dotplot' not in st.session_state:
        st.session_state['sc_dotplot'] = {}
        st.session_state['sc_dotplot']['selected_option'] = None
        st.session_state['sc_dotplot']['selected_folde'] = None
        st.session_state['sc_dotplot']['selected_file'] = None
        st.session_state['sc_dotplot']['gene_list'] = [] 
        st.session_state['sc_dotplot']['out_dir'] = os.path.join('tmp/tools', 'sc_dotplot', ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8)))
    out_dir = st.session_state['sc_dotplot']['out_dir']
    os.makedirs(out_dir, exist_ok=True)
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader('🍓 单细胞基因表达dotplot')
        st.write('对单细数据，选择感兴趣的基因，使用dotplot展示不同分组基因表达均值及表达细胞占比。')
        uploaded_file = st.file_uploader(label='h5ad文件上传:', type='h5ad')
    if uploaded_file:
        adata = sc.read_h5ad(uploaded_file)
        context.write(adata)
        #columns = adata.obs.columns
        select_list = ['使用预定义的基因集','选择基因','31服务器项目导入']
    
    
    params, response, note = st.columns([2, 5, 3])
    with params.form('sc_dotplot'):
        genes=None
        selected_option =  st.selectbox('选择兴趣基因:', select_list)  
        if selected_option=='使用预定义的基因集':
                genes = {}
                st.session_state['pipe']['organism']=st.selectbox('物种:',['human','mouse'])
                celltyps = (st.session_state['pipe']['marker_dict'][st.session_state['pipe']['organism']]).keys()
                select_types = st.multiselect(label='细胞类型:', options=celltyps)       
                if select_types:
                    genes={k: v for k, v in st.session_state['pipe']['marker_dict'][st.session_state['pipe']['organism']].items() if k in select_types}

        elif  selected_option== '选择基因':
                genes=None
                genes = st.multiselect(label='基因:', options=adata.raw.var_names)
                # if submitted:
                #     dot_fig=sc.pl.dotplot(adata,genes,groupby=groupby1,cmap='coolwarm', use_raw=use_raw1,dot_max=dot_max1, show=False,return_fig=True)
                #     save_fig(dot_fig, out_dir, fig_size={'width':fig_width, 'height':fig_height})        
        groupby1 = st.multiselect(label='分组:', options= ['group','cluster','cell_type','sample','leiden','louvain_0.8'])
        use_raw1 = st.checkbox('use_raw', value=False, help='')
        dot_max1= st.number_input(label='dot_max:', min_value=0.1, max_value=1.0, value=0.6, step=0.1)
        fig_width=st.number_input(label='fig_width:', min_value=1.0, max_value=30.0, value=6.0, step=0.5)
        fig_height=st.number_input(label='fig_height:', min_value=1.0, max_value=30.0, value=5.0, step=0.5)
        #dot_min= st.number_input(label='dot_min:', min_value=0, max_value=0.5, value=0, step=0.1)
        submitted = st.form_submit_button()
    if submitted:
        dot_fig=sc.pl.dotplot(adata,genes,groupby=groupby1,cmap='coolwarm', use_raw=use_raw1,dot_max=dot_max1, show=False,return_fig=True)
        save_fig(dot_fig, out_dir, fig_size={'width':fig_width, 'height':fig_height})
# if selected_option=='31服务器项目导入':
    #     base_dir = '/home/ztt/scRNA/scVision/tmp/results/'  
    #     folders = [folder.name for folder in Path(base_dir).iterdir() if folder.is_dir()]  
    #     selected_folder = st.selectbox('选择文件夹', folders)
    #     if selected_folder:  
    #         st.session_state.selected_folder = selected_folder  
          
    #     # 指向子文件夹并获取文件  
    #         subfolder_path = Path(base_dir) / selected_folder / '06-Cell_type_annotation'  
    #         files = [file.name for file in subfolder_path.iterdir() if file.suffix == '.png' and not file.name.startswith(('tsne', 'umap', 'cell'))]  
          
    #     # 显示文件选择  
    #     if files:  
    #         selected_file = st.selectbox('选择文件', files)  
    #         if selected_file:  
    #             # 更新会话状态并处理文件名  
    #             st.session_state.selected_file = selected_file  
    #             gene_list = selected_file.split('_')  
    #             #st.session_state.gene_list = gene_list  
    #             if submitted:
    #                 dot_fig=sc.pl.dotplot(adata,gene_list,groupby=groupby1,cmap='coolwarm', use_raw=use_raw1,dot_max=dot_max1, show=False,return_fig=True)
    #                 save_fig(dot_fig, out_dir, fig_size={'width':fig_width, 'height':fig_height})

    #     else:  
    #         st.write('没有找到符合条件的文件。') 



    
    png_file = os.path.join(out_dir, 'dotplot.png')
    if png_file is not None:
        if os.path.exists(png_file):
                with response.expander("**:blue[结果解读📜]** ", expanded=True):
                    tabs = st.tabs(["🍓 dotplot"])        
                    image_view(tabs[0], png_file=png_file)
                with note.expander("**:blue[参数说明📜]** ", expanded=True):
                    st.markdown(f"""
                        - **选择兴趣基因:** 选择需要展示的基因，其中选择预定义的基因集可以展示标准分析报告中预定义marker，选择基因可自由选择表达基因；
                        - **物种:** 当选择预定义的基因集时，需选定对应物种；
                        - **细胞类型:** 当选择预定义的基因集时，需选定对应细胞类型；
                        - **分组:** dotplot展示所选分组；
                        - **use_raw:** 是否使用原始数据；
                        - **dot_max:** 点的最大值；
                        - **fig_width:** 图片宽度；
                        - **fig_height:** 图片高度；        
                        - **结果文件位置：:rainbow[{st.session_state['sc_dotplot']['out_dir']}]**
                        """)
                zip_filename = os.path.join(out_dir,"dotplot.zip")
                with ZipFile(zip_filename, "w") as zipf:
                    for root, dirs, files in os.walk(out_dir):
                        for file in files:
                            if not file.endswith(".zip"):
                                zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), out_dir))
                with open(zip_filename, "rb") as f:
                    bytes_data = f.read()    
                btn = note.download_button(
                                                label="Download",
                                                data=bytes_data,
                                                file_name='dotplot.zip',
                                                mime="application/zip",
                                                help='Download results.'
                                            )
    else:
        with response.expander("**:blue[结果解读📜]** ", expanded=True):
            st.markdown("""👈**单细胞基因表达dotplot，选择分组信息，设置绘图参数。**""") 
        with note.expander("**:blue[参数说明📜]** ", expanded=True):
            st.markdown(f"""
                        - **选择兴趣基因:** 选择需要展示的基因，其中选择预定义的基因集可以展示标准分析报告中预定义marker，选择基因可自由选择表达基因；
                        - **物种:** 当选择预定义的基因集时，需选定对应物种；
                        - **细胞类型:** 当选择预定义的基因集时，需选定对应细胞类型；
                        - **分组:** dotplot展示所选分组；
                        - **use_raw:** 是否使用原始数据；
                        - **dot_max:** 点的最大值；
                        - **fig_width:** 图片宽度；
                        - **fig_height:** 图片高度；        
                        - **结果文件位置：:rainbow[{st.session_state['sc_dotplot']['out_dir']}]**
                        """)
  