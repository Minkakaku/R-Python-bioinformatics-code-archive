import streamlit as st
from streamlit_pills import pills
from utils.widget.upload_file import save_upload_files, load_data


st.subheader("🍊 :rainbow[Data Loading]")
st.info('''&emsp;&emsp;上传单细胞测序数据用于后续分析，支持10x CellRanger输出结果文件上传。每个样本文件夹内包含CellRanger输出结果filtered_feature_bc_matrix文件夹内的barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz三个文件（必须有且名称不可修改），以及metrics_summary.csv和web_summary.html（允许缺失不上传）。
        文件上传后会根据样本名称生成样本信息表，用户可以在样本信息表中修改样本分组信息，便于后续分析。
        ''')
types = pills("", ["文件上传", "本地数据加载"], ["🔍", "📂"], label_visibility="collapsed")
if types == '文件上传':
    expander = st.expander(f"**:red[{types}]** ", expanded=True)
    cols = expander.columns([4, 4], gap='medium')
    files = cols[0].file_uploader("Upload CellRanger output directory", type=["zip", "tar"], help='将CellRanger输出结果压缩后上传。', accept_multiple_files=False)
    if files is not None:
        st.session_state['pipe']['pid'] = files.name.replace('.zip','').replace('.tar','')
        st.session_state['pipe']['data_dir'] = save_upload_files(files) 
    else:
        cols[0].info(f'上传CellRanger输出结果，首先将数据文件进行压缩，支持.zip和.tar文件格式。', icon="📢")    
    if st.session_state['pipe']['data_dir']:
        load_data(st.session_state['pipe']['data_dir'], cols[1])    
else:
    st.write('本地数据加载功能暂时不可用！')

