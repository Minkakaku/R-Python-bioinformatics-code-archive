import os
import pandas as pd
import scanpy as sc
import streamlit as st
import numpy as np
import random
import matplotlib.pyplot as plt
from zipfile import ZipFile

# 创建临时目录
out_dir = os.path.join('/home/lst/functionenrichment/tmp/tools/sc_paga', ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8)))
os.makedirs(out_dir, exist_ok=True)

# sc_score函数实现所有流程：文件上传、基因集打分和 UMAP 绘制
def sc_paga():
    st.title("🧬 Batch Gene Set Scoring")

    # 创建三列布局
    col1, col2, col3 = st.columns([1, 2, 1])

    # 左侧栏：参数输入
    with col1:
        st.header("🛠️ Parameters")

        # Step 1: 上传 H5AD 文件
        h5ad_file = st.file_uploader("Upload your H5AD file", type=["h5ad"])
        # 输入参数
        groupby = st.text_input("输入分组参数 (groupby):")  # 假设使用Streamlit的文本输入
        threshold = st.number_input("输入阈值 (threshold):", min_value=0.0, value=0.1)
        # 提交按钮
        submit_button = st.button("Submit")

    # 中间栏：结果显示
    with col2:
        st.header("📊 Results")

        if submit_button and h5ad_file is not None and gmt_file is not None:
            # 将文件存入 out_dir 目录
            h5ad_path = os.path.join("out_dir", h5ad_file.name)

            with open(h5ad_path, "wb") as f:
                f.write(h5ad_file.getbuffer())

            # 使用 scanpy 读取 H5AD 文件
            adata = sc.read_h5ad(h5ad_path)
            st.success(f"Successfully loaded {h5ad_file.name}")
            st.write(f"Data shape: {adata.shape}")

            # paga轨迹分析
            out_dir = 'out_dir/paga/'
            os.makedirs(out_dir, exist_ok=True)


            # PAGA分析
            sc.tl.paga(adata, groups=groupby)
    
            # 生成并保存PAGA比较图
            sc.pl.paga_compare(adata, color=groupby, show=False, threshold=threshold,
                       node_size_scale=2, edge_width_scale=0.4, save=os.path.join(out_dir, 'paga'))

    # 清理操作：删除输入文件和结果文件夹
    #try:
        #os.remove(h5ad_file)  # 删除h5ad文件

    # 右侧栏：操作及结果说明
    with col3:
         st.subheader("操作说明")
         st.markdown("""
        - **h5ad文件上传:** 输入单细胞数据文件；
        - **输入分组参数:** 输入paga分组；
        - **输入阈值:** 输入paga阈值；
        - **结果展示:** 基因集打分后的UMAP图展示；
        - **结果文件位置:** 结果存放在服务器的临时目录中。
        """)

# Streamlit 主函数
def main():
    # 调用 sc_score 函数执行所有功能
    sc_paga()

if __name__ == "__main__":
    main()
