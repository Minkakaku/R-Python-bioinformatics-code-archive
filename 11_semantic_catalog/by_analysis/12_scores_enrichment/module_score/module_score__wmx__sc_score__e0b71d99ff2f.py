import os
import pandas as pd
import scanpy as sc
import streamlit as st
import random
import matplotlib.pyplot as plt
from PIL import Image

# 创建临时目录
out_dir = os.path.join('/mnt/NexData_file/share/scRNA_platform/scVision5/tmp/tools/sc_score', ''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8)))
os.makedirs(out_dir, exist_ok=True)

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

# sc_score函数实现所有流程：文件上传、基因集打分和 UMAP 绘制
@st.fragment
def sc_score():
    st.title("🧬 Batch Gene Set Scoring")

    # 创建三列布局
    col1, col2, col3 = st.columns([1, 2, 1])

    # 左侧栏：参数输入
    with col1:
        st.header("🛠️ Parameters")

        # Step 1: 上传 H5AD 文件
        h5ad_file = st.file_uploader("Upload your H5AD file", type=["h5ad"])
        
        # Step 2: 上传 GMT 文件
        gmt_file = st.file_uploader("Upload your GMT file", type=["gmt"])

        # 提交按钮
        submit_button = st.button("Submit")

    # 中间栏：结果显示
    with col2:
        st.header("Results")

        if submit_button and h5ad_file is not None and gmt_file is not None:
            # 将文件存入 out_dir 目录
            h5ad_path = os.path.join(out_dir, h5ad_file.name)
            gmt_path = os.path.join(out_dir, gmt_file.name)

            with open(h5ad_path, "wb") as f:
                f.write(h5ad_file.getbuffer())
            with open(gmt_path, "wb") as f:
                f.write(gmt_file.getbuffer())

            # 使用 scanpy 读取 H5AD 文件
            adata = sc.read_h5ad(h5ad_path)
            adata = adata.raw.to_adata()
            adata.raw = adata
            adata.obs_names_make_unique()
            st.success(f"Successfully loaded {h5ad_file.name}")
            st.write(f"Data shape: {adata.shape}")

            # 读取 GMT 文件
            gene_sets = pd.read_csv(gmt_path, sep="\t", header=None, index_col=0)
            st.success(f"Successfully loaded {gmt_file.name}")
            st.write(gene_sets.head())

            # 基因集打分并生成 UMAP 图
            out_dir_umap = os.path.join(out_dir, 'umap')
            os.makedirs(out_dir_umap, exist_ok=True)

            num_rows = len(gene_sets)
            
            # List of terms (gene sets) for the dropdown menu
            terms = gene_sets.index.tolist()
            term_plot_paths = {}  # Store the paths of the UMAP plots for all terms

            # Generate UMAP plots for all terms and store them
            for i in range(num_rows):
                # 获取当前基因集的名称和基因列表
                term = gene_sets.index[i]
                genes = gene_sets.iloc[i, 2:].dropna().to_list()  # 获取基因集

                # 使用 scanpy 进行基因集打分
                sc.tl.score_genes(adata, gene_list=genes, score_name=term)

                # Generate UMAP plot for the selected term
                sc.pl.umap(adata, color=term, size=25, show=False)
                file_name = f'umap_{term}'
                save_fig(plt, out_dir=out_dir_umap, file_name=file_name, fig_size={'w':8, 'h':6})

                # Store the path to the plot for later display
                term_plot_paths[term] = os.path.join(out_dir_umap, f'{file_name}.png')

            # Store the paths of UMAP plots in session_state to persist them across interactions
            st.session_state['term_plot_paths'] = term_plot_paths

            # Provide a dropdown for selecting the term
            selected_term = st.selectbox("Select a gene set term to visualize:", terms)

            if selected_term:
                # Display the UMAP plot for the selected term
                img_umap = Image.open(st.session_state['term_plot_paths'][selected_term])
                st.image(img_umap, caption=f'UMAP plot for gene set: {selected_term}')

            # Remove files after processing
            os.remove(h5ad_path)
            os.remove(gmt_path)

        elif submit_button:
            st.warning("Please upload both H5AD and GMT files.")

    # 右侧栏：操作及结果说明
    with col3:
         st.subheader("操作说明")
         st.markdown("""
        - **h5ad文件上传:** 输入单细胞数据文件；
        - **gmt文件上传:** 输入基因集文件；
        - **结果展示:** 基因集打分后的UMAP图展示；
        - **结果文件位置:** 结果存放在服务器的临时目录中。
        """)

# Streamlit 主函数
def main():
    # 调用 sc_score 函数执行所有功能
    sc_score()

if __name__ == "__main__":
    main()
