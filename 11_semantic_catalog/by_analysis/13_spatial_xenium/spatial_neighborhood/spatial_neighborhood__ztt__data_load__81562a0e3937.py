import os
import pandas as pd
import gradio as gr
from functions.plotter import Plotter
from functions.xenium_analysis import load_xenium_data
from datetime import datetime

data_dir_xenium = './data/xenium'
data_dir_hd = '/mnt/NexData_file/nas/outputdata/spaceranger_out'
data_dir_stereo = './data/stereo'
output_dir = './output'
output_dir = './output'

# 简化的状态卡片 HTML 生成工具函数（保持原有样式与内容）
def status_card_html(kind: str, inner_html: str) -> str:
    cls = {
        'wait': 'status-wait',
        'success': 'status-success',
        'error': 'status-error',
        'info': 'status-info',
    }.get(kind, 'status-info')
    # 将内容包裹在 status-content，以适配全局CSS的flex布局并实现垂直换行
    return f"<div class='status-card {cls}'><div class='status-content'>{inner_html}</div></div>"

def get_data_dir(data_type: str):
    paths = {"Xenium": data_dir_xenium, "Stereo": data_dir_stereo, "VisiumHD": data_dir_hd}
    path = paths.get(data_type)
    if not path or not os.path.exists(path):
        return ["暂无可用数据文件夹"]
    folders = [e.name for e in os.scandir(path) if e.is_dir()]
    return ["请选择"] + folders if folders else ["暂无可用数据文件夹"]

def update_folder_dropdown(data_type: str):
    folders = get_data_dir(data_type)
    default = folders[0] if folders and not folders[0].startswith("暂无") else None
    return gr.update(choices=folders, value=default)

def process_uploaded_file(file):
    if file is None:
        return gr.update(visible=False)
    ext = os.path.splitext(file.name)[1].lower()
    readers = {
        '.csv': pd.read_csv,
        '.tsv': lambda p: pd.read_csv(p, sep='\t'),
        '.txt': lambda p: pd.read_csv(p, sep='\t'),
        '.xlsx': pd.read_excel,
    }
    reader = readers.get(ext)
    if not reader:
        return gr.update(visible=False)
    try:
        df = reader(file.name)
        return gr.update(value=df, visible=True)
    except Exception as e:
        return gr.update(value=f"读取文件时出错: {e}", visible=True)

def data_load_tab(adata_raw, output_dir_state):
    with gr.Row():
        with gr.Column(scale=1):
            with gr.Accordion("🔔 选择数据类型和项目编号", open=True, elem_classes=["data-selection-accordion"]):
                data_type_dropdown = gr.Dropdown(
                    choices=["请选择", "VisiumHD", "Xenium", "Stereo"],
                    label="数据类型",
                    value="请选择",
                    interactive=True,
                    scale=1
                )
                available_folders = get_data_dir(data_type_dropdown.value)
                folder_dropdown = gr.Dropdown(
                    choices=available_folders,
                    label="项目编号",
                    value=available_folders[0] if available_folders and available_folders[0] != "暂无可用数据文件夹" else None,
                    interactive=True,
                    scale=1
                )
                sample_file_upload = gr.File(
                    label="上传Sample文件",
                    file_types=[".csv", ".txt", ".tsv", ".xlsx"],
                    file_count="single",
                    interactive=True
                )
                sample_dataframe = gr.Dataframe(
                    visible=False,
                )
            
            submit_btn = gr.Button("🚀 提交", elem_classes="custom-button")
        
        # 右侧结果展示面板
        with gr.Column(scale=2, elem_classes=["result-panel-compact"]):
            with gr.Tabs():
                with gr.TabItem("📄 数据概览"):
                    # 将“数据概览”和“QC质控图”合并到同一个卡片容器中
                    with gr.Column(scale=1, elem_classes=["card-section"]):
                        # 顶部两列：基本信息/文件夹列表 与 质量指标
                        with gr.Row(equal_height=True):
                            with gr.Column(scale=1):
                                basic_info = gr.HTML(
                                    value=status_card_html(
                                        'wait',
                                        "<h4>📋 等待数据选择</h4><br/><p>选择数据文件夹后，这里将显示数据的基本信息</p>"
                                    )
                                )
                                folder_list_dataframe = gr.Dataframe(
                                    headers=["文件夹名称"],
                                    value=pd.DataFrame({"文件夹名称": ["请先选择数据地址"]}),
                                    elem_classes=["data-table-compact"],
                                    interactive=False,
                                    label="📁 数据地址文件夹内容",
                                    visible=True
                                )
                            with gr.Column(scale=1):
                                quality_metrics = gr.Dataframe(
                                    headers=["Label", "Value"],
                                    elem_classes=["data-table-compact"],
                                    column_widths=["40%", "60%"],
                                    interactive=False,
                                    visible=False
                                )
                        # 底部：QC质控图，放入同一卡片容器中
                        with gr.Row():
                            qc_plot = gr.Plot(
                                label="QC 质控图",
                                elem_classes=["plot-compact"],
                                visible=False
                            )

    # 事件处理函数
    def process_selected_data(data_type, folder_name, uploaded_file, output_dir_state):
        # 提交后创建输出目录结构
        try:
            if not folder_name or folder_name == "请选择" or str(folder_name).startswith("暂无可用"):
                raise ValueError("未选择有效的项目编号")
            os.makedirs(output_dir, exist_ok=True)
            ts = datetime.now().strftime('%Y%m%d_%H%M%S')
            project_output = os.path.join(output_dir, f"{folder_name}_{ts}")
            os.makedirs(project_output, exist_ok=True)
            subdirs = [
                '01-Data_QC',
                '02-Cell_filter',
                '03-Dimension_reduction_cluster',
                '04-Marker_gene_analysis',
                '05-Cell_type_annotation',
                '06-Computation_spatial_statistics',
                '07-Cell_Neighborhood_analysis',
                '08-Proximity_analysis',
                '09-Nearest_neighbor_analysis',
            ]
            for sd in subdirs:
                os.makedirs(os.path.join(project_output, sd), exist_ok=True)
        except Exception as e:
            error_msg = f"📁 输出目录创建失败: {str(e)}"
            return (
                gr.update(value=status_card_html('error', f"<p>{error_msg}</p>")),
                gr.update(value=None, visible=False),
                gr.update(value=None, visible=False),
                None,
                None
            )

        if data_type == 'Xenium':
            try:
                adata = load_xenium_data(
                    data_dir=os.path.join(data_dir_xenium, folder_name),
                    sample_info=uploaded_file,
                    out_dir=project_output
                )
                qc_metrics = pd.DataFrame({
                    '名称': ['鉴定细胞数', '鉴定基因数', '样本名称', '分组名称', '细胞鉴定基因中位数'],
                    '数值': [
                        f"{adata.n_obs:,}",
                        f"{adata.n_vars:,}",
                        f"{', '.join(adata.obs['sample'].unique())}",
                        f"{', '.join(adata.obs['group'].unique())}",
                        f"{adata.obs['n_genes_by_counts'].median():.0f}",
                    ]
                })

                plot_handler = Plotter()
                qc_plot = plot_handler.plot_qc_violin(adata)

                basic_info = status_card_html(
                    'success',
                    f"<p>✅ 数据加载成功: {adata.n_obs:,} 细胞 × {adata.n_vars:,} 基因</p><//b><p>🗂️ 输出目录：{project_output}</p>"
                )
                return (
                    basic_info,
                    gr.update(value=qc_metrics, visible=True),
                    gr.update(value=qc_plot, visible=True),
                    adata,
                    project_output  # 更新输出目录状态
                )

            except Exception as e:
                error_msg = f"❌ 处理失败: {str(e)}"
                return (
                    gr.update(value=status_card_html('error', f"<p>{error_msg}</p><p>🗂️ 输出目录：{project_output}</p>")),
                    gr.update(value=None, visible=False),
                    gr.update(value=None, visible=False),
                    None,
                    project_output  # 即使数据处理失败，仍保留已创建的输出目录
                )

    def get_subfolder_list(data_type, folder_name):
        paths = {"Xenium": data_dir_xenium, "Stereo": data_dir_stereo, "VisiumHD": data_dir_hd}
        data_path = paths.get(data_type)
        if not data_path:
            return pd.DataFrame({"文件夹名称": ["无效的数据类型"]})
        full_path = os.path.join(data_path, folder_name)
        if not os.path.exists(full_path):
            return pd.DataFrame({"文件夹名称": ["路径不存在"]})
        try:
            subfolders = sorted([e.name for e in os.scandir(full_path) if e.is_dir()])
            return pd.DataFrame({"文件夹名称": subfolders or ["该文件夹为空"]})
        except Exception as e:
            return pd.DataFrame({"文件夹名称": [f"读取失败: {str(e)}"]})

    def update_summary_display(data_type, folder_name):
        valid = folder_name and folder_name != "请选择" and not str(folder_name).startswith("暂无可用")
        if not valid:
            html = status_card_html('wait', "<h4>📋 等待数据选择</h4><br/><p>选择数据文件夹后，这里将显示数据的基本信息</p>")
            return html, gr.update(visible=False)
        folder_list = get_subfolder_list(data_type, folder_name)
        first_cell = folder_list.iloc[0, 0] if not folder_list.empty else ""
        has_items = first_cell not in ["路径不存在", "该文件夹为空"] and not str(first_cell).startswith("读取失败")
        count = len(folder_list) if has_items else 0
        html = status_card_html('success', f"<h4>📁 数据地址已选择</h4><p><strong>数据类型:</strong> {data_type} | <strong>文件夹:</strong> {folder_name} | <strong>子文件夹:</strong> {count} 个</p>")
        return html, gr.update(value=folder_list, visible=True)

    # 绑定数据类型下拉框change事件，动态更新文件夹下拉框
    data_type_dropdown.change(
        fn=update_folder_dropdown,
        inputs=[data_type_dropdown],
        outputs=[folder_dropdown]
    )

    # 绑定数据地址选择事件，更新Summary显示内容
    folder_dropdown.change(
        fn=update_summary_display,
        inputs=[data_type_dropdown, folder_dropdown],
        outputs=[basic_info, folder_list_dataframe]
    )

    # 绑定文件上传事件，自动更新dataframe和可见性
    sample_file_upload.change(
        fn=process_uploaded_file,
        inputs=sample_file_upload,
        outputs=[sample_dataframe]
    )

    submit_btn.click(
        fn=process_selected_data,
        inputs=[data_type_dropdown, folder_dropdown, sample_file_upload, output_dir_state],
        outputs=[basic_info, quality_metrics, qc_plot, adata_raw, output_dir_state]
    )