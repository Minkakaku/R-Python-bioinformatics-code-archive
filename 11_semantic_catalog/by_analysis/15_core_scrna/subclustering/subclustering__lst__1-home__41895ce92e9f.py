import os
from PIL import Image  
import streamlit as st
import streamlit_antd_components as sac


with st.container():
    st.markdown(f'''
        <style>
        .stApp .main .block-container{{
            padding-top:30px
        }}
        iframe{{
            display:block;
        }}
        </style>
        ''', unsafe_allow_html=True)
    sac.alert('**SubCluster tools added complete.**', banner=True, closable=False, color='success')

st.subheader("🧰 scVision: A web-drived platform for rapid analysis of single cell data using Scanpy.")
st.write('''<p style="text-indent: 2em;"><strong>
            Large-scale single-cell RNA sequencing (scRNA-seq) has emerged as a robust method for dissecting cellular heterogeneity at single-cell resolution. 
            However, to meet the increasingly high computational demands of non-programming experts, a user-friendly, scalable, and accessible online platform for analyzing scRNA-seq data is urgently needed. 
            Here, we have developed a web-based platform scVision use Scanpy (a scalable toolkit for analyzing single-cell gene expression data) that enables online massive single-cell transcriptome analysis, improving interactivity and reproducibility using high-quality visualization frameworks.
            provides easy access to interactive visualization, customized parameters, and publication-quality graphs. Furthermore, it comprehensively integrates preprocessing, clustering, developmental trajectory inference, cell-cell communication, cell-type annotation, subcluster analysis, and pathway enrichment.</strong></p>
         ''', unsafe_allow_html=True)

st.subheader("🍁 Workflow")
cols = st.columns([2, 6, 4], gap='medium')
workflow = Image.open(os.path.join('config', 'images', 'sc_workflow_integration.png'))
cols[1].image(workflow, caption='Scanpy Workflow Integration')

st.subheader("🧲 Development environment")
st.write('''
            - <h5>python > v3.8</h5>
            - <h5>R v4.3.1</h5>
            - <h5>streamlit v1.25.0</h5>
            - <h5>scanpy v1.9.3</h5>
            - <h5>pandas v1.5.3</h5>
            - <h5>scipy v1.10.1</h5>
         ''', unsafe_allow_html=True)
