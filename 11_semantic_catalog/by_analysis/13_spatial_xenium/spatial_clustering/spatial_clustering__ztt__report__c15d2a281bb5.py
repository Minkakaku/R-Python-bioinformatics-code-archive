import gradio as gr
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from utils.data_handler import DataHandler
from utils.plotter import Plotter
import io
import base64

def create_report_tab(adata_state):
    data_handler = DataHandler()
    plotter = Plotter()
    
    with gr.Row(equal_height=True):
        # 左侧参数面板
        with gr.Column(scale=1):
            gr.HTML("<h3>📋 报告设置</h3>")
            
            with gr.Group():
                gr.HTML("<h4>报告内容</h4>")
                
                include_sections = gr.CheckboxGroup(
                    choices=["数据概览", "质控结果", "降维聚类", "细胞类型注释", "标记基因", "方法描述"],
                    value=["数据概览", "质控结果", "降维聚类", "细胞类型注释"],
                    label="选择包含的章节",
                    interactive=True
                )
                
                report_format = gr.Radio(
                    choices=["HTML", "PDF", "Word", "Markdown"],
                    value="HTML",
                    label="报告格式",
                    interactive=True
                )
            
            with gr.Group():
                gr.HTML("<h4>图表设置</h4>")
                
                plot_quality = gr.Dropdown(
                    choices=["标准 (150 DPI)", "高 (300 DPI)", "超高 (600 DPI)"],
                    value="高 (300 DPI)",
                    label="图片质量",
                    interactive=True
                )
                
                plot_width = gr.Slider(
                    minimum=800,
                    maximum=2000,
                    value=1200,
                    step=100,
                    label="图片宽度 (像素)",
                    interactive=True
                )
                
                plot_height = gr.Slider(
                    minimum=600,
                    maximum=1500,
                    value=800,
                    step=100,
                    label="图片高度 (像素)",
                    interactive=True
                )
            
            with gr.Group():
                gr.HTML("<h4>报告信息</h4>")
                
                report_title = gr.Textbox(
                    label="报告标题",
                    value="单细胞数据分析报告",
                    interactive=True
                )
                
                author_name = gr.Textbox(
                    label="作者",
                    placeholder="请输入您的姓名",
                    interactive=True
                )
                
                project_description = gr.TextArea(
                    label="项目描述",
                    placeholder="请简要描述该项目...",
                    lines=3,
                    interactive=True
                )
            
            with gr.Group():
                gr.HTML("<h4>高级选项</h4>")
                
                include_code = gr.Checkbox(
                    value=False,
                    label="包含分析代码",
                    interactive=True
                )
                
                include_parameters = gr.Checkbox(
                    value=True,
                    label="包含分析参数",
                    interactive=True
                )
                
                create_summary_table = gr.Checkbox(
                    value=True,
                    label="创建结果汇总表",
                    interactive=True
                )
            
            generate_report_btn = gr.Button("📊 生成报告", elem_classes="custom-button")
            preview_report_btn = gr.Button("👁️ 预览报告", variant="secondary")
        
        # 右侧结果展示面板
        with gr.Column(scale=2):
            gr.HTML("<h3>📄 报告预览</h3>")
            
            with gr.Tabs():
                with gr.TabItem("报告预览"):
                    report_preview = gr.HTML()
                
                with gr.TabItem("报告摘要"):
                    report_summary = gr.Dataframe(
                        headers=["章节", "状态", "图表数量", "表格数量"],
                        label="报告摘要",
                        elem_classes="data-table"
                    )
                
                with gr.TabItem("图表预览"):
                    with gr.Row():
                        summary_plot1 = gr.Plot(label="数据概览", scale=1)
                        summary_plot2 = gr.Plot(label="质控结果", scale=1)
                    with gr.Row():
                        summary_plot3 = gr.Plot(label="聚类结果", scale=1)
                        summary_plot4 = gr.Plot(label="细胞类型", scale=1)
                
                with gr.TabItem("下载报告"):
                    with gr.Row():
                        download_html = gr.File(label="HTML报告", scale=1)
                        download_pdf = gr.File(label="PDF报告", scale=1)
                    with gr.Row():
                        download_summary = gr.File(label="汇总表格", scale=1)
                        download_plots = gr.File(label="图表压缩包", scale=1)
            
            # 生成状态
            report_status = gr.Textbox(
                label="报告状态",
                value="等待生成报告...",
                interactive=False,
                lines=3
            )
    
    # 事件处理函数
    def generate_report(include_sections, report_format, plot_quality, plot_width, plot_height,
                       report_title, author_name, project_description, include_code,
                       include_parameters, create_summary_table):
        try:
            status_msg = "🔄 开始生成报告...\n"
            
            # 获取当前时间
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # 模拟报告生成过程
            status_msg += "📊 收集分析结果...\n"
            
            # 这里应该收集之前的分析结果
            # 现在使用模拟数据
            n_cells = 3000
            n_genes = 1800
            n_clusters = 8
            n_celltypes = 6
            
            status_msg += "📈 生成图表...\n"
            
            # 生成模拟图表
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
            
            # 模拟数据概览图
            ax1.bar(['细胞数', '基因数'], [n_cells, n_genes])
            ax1.set_title('数据概览')
            
            # 模拟质控图
            ax2.hist(np.random.normal(1000, 200, 1000))
            ax2.set_title('质控分布')
            
            # 模拟聚类图
            ax3.scatter(np.random.randn(100), np.random.randn(100), c=np.random.randint(0, 8, 100))
            ax3.set_title('聚类结果')
            
            # 模拟细胞类型图
            celltype_counts = [500, 400, 350, 300, 250, 200]
            ax4.pie(celltype_counts, labels=['T cells', 'B cells', 'Monocytes', 'NK cells', 'Megakaryocytes', 'Erythrocytes'])
            ax4.set_title('细胞类型分布')
            
            plt.tight_layout()
            
            # 保存图表
            plots_path = "/tmp/summary_plots.png"
            plt.savefig(plots_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            status_msg += "📝 编写报告内容...\n"
            
            # 生成报告预览
            report_content = f"""
            <div style='padding: 20px; font-family: Arial, sans-serif;'>
                <h1 style='color: #2d3748; text-align: center;'>{report_title}</h1>
                <div style='text-align: center; color: #4a5568; margin-bottom: 30px;'>
                    <p><strong>作者:</strong> {author_name or '未指定'}</p>
                    <p><strong>生成时间:</strong> {current_time}</p>
                    <p><strong>项目描述:</strong> {project_description or '未提供'}</p>
                </div>
                
                <div style='background: #f7fafc; padding: 20px; border-radius: 10px; margin-bottom: 20px;'>
                    <h2>📊 分析摘要</h2>
                    <table style='width: 100%; border-collapse: collapse;'>
                        <tr><td><strong>总细胞数:</strong></td><td>{n_cells:,}</td></tr>
                        <tr><td><strong>检测到的基因数:</strong></td><td>{n_genes:,}</td></tr>
                        <tr><td><strong>识别的聚类数:</strong></td><td>{n_clusters}</td></tr>
                        <tr><td><strong>注释的细胞类型数:</strong></td><td>{n_celltypes}</td></tr>
                    </table>
                </div>
                
                <div style='background: #edf2f7; padding: 20px; border-radius: 10px;'>
                    <h2>📋 包含章节</h2>
                    <ul>
                        {''.join([f'<li>{section}</li>' for section in include_sections])}
                    </ul>
                </div>
                
                <div style='text-align: center; margin-top: 30px;'>
                    <p style='color: #48bb78; font-size: 18px;'>✅ 报告生成成功！</p>
                </div>
            </div>
            """
            
            # 生成报告摘要表格
            summary_data = []
            for section in include_sections:
                summary_data.append([
                    section,
                    "✅ 完成",
                    str(np.random.randint(3, 8)),
                    str(np.random.randint(1, 4))
                ])
            
            summary_df = pd.DataFrame(summary_data, columns=['章节', '状态', '图表数量', '表格数量'])
            
            # 模拟文件生成
            html_file = "/tmp/report.html"
            pdf_file = "/tmp/report.pdf"
            summary_file = "/tmp/summary.csv"
            plots_zip = "/tmp/plots.zip"
            
            # 创建模拟文件
            with open(html_file, 'w') as f:
                f.write(report_content)
            
            # 创建汇总表格
            summary_df.to_csv(summary_file, index=False)
            
            status_msg += "✅ 报告生成完成！"
            
            return (
                gr.update(value=report_content),
                gr.update(value=summary_df),
                gr.update(value=fig),
                gr.update(value=fig),
                gr.update(value=fig),
                gr.update(value=fig),
                gr.update(value=html_file),
                gr.update(value=pdf_file),
                gr.update(value=summary_file),
                gr.update(value=plots_zip),
                status_msg
            )
            
        except Exception as e:
            error_msg = f"❌ 报告生成失败: {str(e)}"
            return (
                gr.update(value=f"<p style='color: red;'>{error_msg}</p>"),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                error_msg
            )
    
    def preview_report():
        try:
            preview_msg = "👁️ 预览功能开发中..."
            return preview_msg
        except Exception as e:
            return f"❌ 预览失败: {str(e)}"
    
    # 事件绑定
    generate_report_btn.click(
        fn=generate_report,
        inputs=[include_sections, report_format, plot_quality, plot_width, plot_height,
                report_title, author_name, project_description, include_code,
                include_parameters, create_summary_table],
        outputs=[report_preview, report_summary, summary_plot1, summary_plot2, summary_plot3, summary_plot4,
                download_html, download_pdf, download_summary, download_plots, report_status]
    )
    
    preview_report_btn.click(
        fn=preview_report,
        outputs=[report_status]
    )