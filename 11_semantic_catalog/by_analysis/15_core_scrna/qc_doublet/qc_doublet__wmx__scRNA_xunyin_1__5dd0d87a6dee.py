import re
import os
os.environ['OPENBLAS_NUM_THREADS'] = "64"
os.environ['MKL_NUM_THREADS'] = "64"
os.environ['OMP_NUM_THREADS'] = "64"
import shutil
import pandas as pd
from PIL import Image
import streamlit as st
import scanpy as sc
import streamlit_antd_components as sac
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
from moudle.common import *
from config import out_dir, data_dir,adata_dir, marker_dict_file, rmd_template_xunyin, rmd_template_img_xunyin
from moudle.xunyin_scanpy_pipeline import get_xunyin_mtx,make_report_dir


horizontal_bar = "<hr style='margin-top: 2; margin-bottom: 0; height: 1px; border: 1px solid #635985;'><br>"
hide_table_row_index = """
                        <style>
                            thead tr th:first-child {display:none}
                            tbody th {display:none}
                        </style>
                        """



def data_server(pid, cols):
    tmp = []
    input_dir = os.path.join(data_dir, pid)
    sample_file = os.path.join(input_dir, 'sample.txt')
    st.session_state['out_dir'] = os.path.join(out_dir, pid)
    metrics_summary = os.path.join(st.session_state['out_dir'], '01-seeksoultools', 'metrics_summary.csv')
    
    if not os.path.exists(sample_file):
        tmp = []
        for sample in os.listdir(input_dir):
            if os.path.isdir(os.path.join(input_dir, sample)):
                tmp.append({"sample": sample, "rename": sample, "group": sample})
        sample_df = cols[0].data_editor(pd.DataFrame(tmp), hide_index=True, disabled=('sample',), use_container_width=False)             
        sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
        sample_df.to_csv(os.path.join(input_dir, 'sample.txt'), sep='\t', index=None)
    else:
        df = pd.read_csv(sample_file, sep='\t')
        sample_df = cols[0].data_editor(df, hide_index=True, disabled=('sample',), use_container_width=True)             
        sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
        sample_df.to_csv(os.path.join(input_dir, 'sample.txt'), sep='\t', index=None)

    submit = cols[0].button('🚀Submit', type='primary')
    
    if submit:
        with cols[1]:
            with st.spinner('数据加载中...'):
                make_report_dir(st.session_state['out_dir'])
                sample_dict = sample_df.groupby('group')['sample'].apply(list).to_dict()
                adata = get_xunyin_mtx(data_path=input_dir, sample_dict=sample_dict, out_dir=st.session_state['out_dir'])
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
            metrics_summary = os.path.join(st.session_state['out_dir'], '01-seeksoultools', 'metrics_summary.csv')
            sample_list = list(adata.obs['sample'].drop_duplicates())
            try:
                metrics_sample = []
                for s in sample_list:
                    metrics = pd.read_csv(os.path.join(input_dir, s, s, s+'_summary.csv'))
                    metrics['Sample Name'] = s
                    cols = list(metrics)
                    cols = [cols[-1]] + cols[:-1]
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
      

    # 质量控制和数据可视化
    if 'adata_qc' in st.session_state:
        adata = st.session_state['adata_qc']
        tab1, tab2, tab3 = cols[1].tabs(["🔔Summary", "📜Violin Plot", "📈Scatter Plot"])
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
        axes = plt.gcf().axes
        for ax in axes:
            ax.set_xticklabels(ax.get_xticklabels(), horizontalalignment='right')
        save_fig(plt, qc_out_dir, file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, qc_out_dir, file_name='orig_nCount_nFeature_scatter')
        
        pd.DataFrame(adata.obs).to_csv(os.path.join(qc_out_dir, 'all_cell_features_raw.xls'), sep='\t')
        with tab1:
            # 显示QC数据的基本信息
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("细胞总数", adata.X.shape[0])
            col2.metric("基因总数", adata.X.shape[1])
            col3.metric("样本个数", len(adata.obs['sample'].unique()))
            col4.metric("样本分组", len(adata.obs['group'].unique()))

            st.markdown(hide_table_row_index, unsafe_allow_html=True)

        # QC小提琴图
        with tab2:
            img_violin = Image.open(os.path.join(qc_out_dir, 'orig_nCount_nFeature_voilon.png'))
            st.image(img_violin, caption='注：nGene图：过滤前各个样本单个细胞中检测到的基因数量分布情况；nUMI图：过滤前各个样本单个细胞中检测到的nUMI分布情况；percent.mt图：过滤前各个样本单个细胞中线粒体基因表达量的百分比分布情况；percent.ribo图：过滤前各个样本单个细胞中核糖体基因表达量的百分比分布情况。')

        # QC散点图
        with tab3:
            img_scatter = Image.open(os.path.join(qc_out_dir, 'orig_nCount_nFeature_scatter.png'))
            st.image(img_scatter, caption='注：左图表示细胞中检测到的UMI数目和线粒体基因含量的关系；右图表示细胞中检测到的UMI数目和表达基因数目的关系。')

    st.markdown(horizontal_bar, True)
      

def generate_report(pid):
    if 'adata_cell_type' in st.session_state:
        # sac.divider(align='center', size='lg', color='cyan')
        st.subheader(f"🍑 生成报告, 项目编号: {pid}, 文件路径: {os.path.abspath(st.session_state['out_dir'])}")
        out_dir = os.path.abspath(st.session_state["out_dir"])
        zip_outdir = os.path.abspath(os.path.join(out_dir, '..'))
        with st.form('report'): 
            cols = st.columns([3,5,4,3], gap='large')
            pid = cols[0].text_input('项目编号:', pid)
            ptitle = cols[1].text_input('项目名称:', '')
            paddress = cols[2].text_input('客户单位:', '')
            pcustmor = cols[3].text_input('客户姓名:', '')               
            submitted = st.form_submit_button("🚀Report", type='primary')
        
        if submitted:
            rmd = ''
            with open(rmd_template_xunyin, 'r', encoding='utf-8') as f:
                line = f.read()
                if not os.path.exists(os.path.join(out_dir, '01-Dnbc4tools', 'metrics_summary.csv')):
                    line = re.sub(r'## 测序数据统计.*? 所有细胞的基因总数。', '', line, flags=re.DOTALL)
                sample_num = len(st.session_state['adata_qc'].obs['sample'].drop_duplicates().tolist())
                cell_num_raw = st.session_state['adata_qc'].X.shape[0]
                cell_num_after_filter = str(st.session_state['adata_filtered'].X.shape[0])
                cell_types = ', '.join(st.session_state['adata_cell_type'].obs['cell_type'].drop_duplicates().tolist())                   
                cluster_num = len(st.session_state['adata_cell_type'].obs["cluster"].drop_duplicates().tolist())
                rmd = line.replace('%pid%', pid).replace('%ptitle%', ptitle).replace('%paddress%', paddress).replace('%pcustmor%', pcustmor).replace('%num_sample%', str(sample_num)).replace('%organism%', st.session_state['organism'])
                rmd = rmd.replace('%cell_num_raw%', str(cell_num_raw)).replace('%cell_num_after_filter%', str(cell_num_after_filter)).replace('%cluster_num%', str(cluster_num)).replace('%cell_types%', cell_types) 
                rmd = rmd.replace('%mt_pct%', str(st.session_state['pct_counts_mt']))
                rmd = rmd.replace('%min_genes%', str(st.session_state['min_genes']))
                rmd = rmd.replace('%upper_lim%', str(st.session_state['upper_lim']))
                rmd = rmd.replace('images_dir', rmd_template_img_xunyin)

            with open(os.path.join(out_dir, 'report.rmd'), 'w', encoding='utf-8') as w:
                w.write(rmd)
            with st.spinner('Report generation ...'):   
                os.system("R -e \"rmarkdown::render('{0}/report.rmd',output_file='{0}/report.html')\"".format(out_dir))
            with st.spinner('File compression ...'):
                st.toast('分析结果压缩中，请稍等!', icon='💗')
                st.session_state['adata_cell_type'].write(os.path.join(out_dir, 'result.h5ad'))
                shutil.make_archive(os.path.join(zip_outdir, pid), 'zip', out_dir)        

        if os.path.exists(os.path.join(out_dir, 'report.html')): 
            cols = st.columns(8)
            if os.path.exists(os.path.join(zip_outdir, pid+'.zip')):
                cols[7].download_button(
                                        label="📁 Download",
                                        data=open(os.path.join(zip_outdir, pid+'.zip'), "rb"),
                                        file_name=f"{pid}.zip",
                                        mime="application/octet-stream",
                                        )
            report = open(os.path.join(out_dir, 'report.html'))
            components.html(report.read(), height=1000, scrolling=True)   

 
def main():
    st.subheader('🍄 寻因 scRNA标准分析')
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
        st.info('根据项目编号，加载/mnt/NexData_file/nas/outputdata/cellranger_out路径下CellRanger输出结果文件，选择项目编号后会根据样本名称生成样本信息表，用户可以在样本信息表中修改样本分组信息，用于后续分析', icon="📢")
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
            generate_report(pid)


if __name__ == '__main__':
    main()
