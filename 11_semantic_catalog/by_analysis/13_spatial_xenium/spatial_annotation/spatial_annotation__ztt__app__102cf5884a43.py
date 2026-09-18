import gradio as gr
import os
from pages.data_load import data_load_tab
from pages.preprocess import data_preprocess_tab
from pages.cell_annot import cell_annot_tab
from pages.annotation import create_annotation_tab
from pages.report import create_report_tab


css_path = os.path.join(os.path.dirname(__file__), "style.css")
try:
    with open(css_path, "r", encoding="utf-8") as f:
        css = f.read()
except FileNotFoundError:
    css = ""

def main():
    with gr.Blocks(css=css, title="Spatial Transcriptomics Data Analysis Platform", theme=gr.themes.Soft()) as app:
        # 全局状态存储
        adata_raw = gr.State(value=None)
        adata_filter = gr.State(value=None)
        adata_cluster = gr.State(value=None)
        adata_celltype = gr.State(value=None)
        adata_state = gr.State(value=None)
        output_dir_state = gr.State(value=None)
        
        # 页面头部
        gr.HTML("""
        <div class="page-header">
            <div class="header-content">
                <h1 class="page-title">🔬 空间转录组数据分析</h1>
                <p class="page-description">
                    Visium HD、Xenium、Stereo空间转录组技术平台，从数据预处理、质量控制到细胞社区等完整分析流程；调整参数，生成标准分析报告。
                 </p> 
            </div>
        </div>
        """)
        
        # 主内容区域
        with gr.Column(elem_classes=["main-content"]):
            with gr.Tabs() as tabs:
                with gr.TabItem("📁 数据选择", id=0):
                    data_load_tab(adata_raw, output_dir_state)
                
                # 数据质控标签页
                with gr.TabItem("🔍 细胞过滤", id=1):
                    # 传入原始数据状态、聚类后的全局状态与输出目录状态
                    data_preprocess_tab(adata_raw, adata_cluster, output_dir_state)
                
                # 降维聚类标签页
                with gr.TabItem("📊 细胞类型注释", id=2):
                    cell_annot_tab(adata_cluster, output_dir_state)
                
                # # 细胞类型注释标签页
                # with gr.TabItem("🏷️ 细胞类型注释", id=3):
                #     create_annotation_tab(adata_state)
                
                # # 生成报告标签页
                # with gr.TabItem("📄 生成报告", id=4):
                #     create_report_tab(adata_state)
        
        # 页面底部
        gr.HTML("""
        <div class="page-footer">
            <div class="footer-content">
                <div class="footer-links">
                    <a href="#" class="footer-link">📚 使用文档</a>
                    <a href="#" class="footer-link">🔧 API 接口</a>
                    <a href="#" class="footer-link">💡 功能建议</a>
                    <a href="#" class="footer-link">🐛 问题反馈</a>
                    <a href="#" class="footer-link">📧 联系我们</a>
                </div>
                <div class="footer-copyright">
                    © 2024 单细胞数据分析平台. 基于 Gradio 构建 | 版本 v1.0.0
                </div>
            </div>
        </div>
        """)
    
    return app

demo = main()

if __name__ == "__main__":
    demo.launch(
        server_name="0.0.0.0",
        server_port=7860
    )