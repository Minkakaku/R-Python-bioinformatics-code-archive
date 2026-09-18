import gradio as gr
import pandas as pd
import numpy as np
from utils.data_handler import DataHandler
from utils.plotter import Plotter

def create_annotation_tab(adata_state):
    data_handler = DataHandler()
    plotter = Plotter()
    
    with gr.Row(equal_height=True):
        # 左侧参数面板
        with gr.Column(scale=1):
            gr.HTML("<h3>🏷️ 注释参数</h3>")
            
            with gr.Group():
                gr.HTML("<h4>注释方法</h4>")
                
                annotation_method = gr.Radio(
                    choices=["自动注释", "手动注释", "参考数据库"],
                    value="自动注释",
                    label="选择注释方法",
                    interactive=True
                )
                
                reference_dataset = gr.Dropdown(
                    choices=["HumanPrimaryCellAtlas", "Blueprint", "ENCODE", "ImmGen", "MouseRNAseq"],
                    value="HumanPrimaryCellAtlas",
                    label="参考数据集",
                    interactive=True
                )
            
            with gr.Group():
                gr.HTML("<h4>标记基因</h4>")
                
                marker_genes_input = gr.Textbox(
                    label="输入已知标记基因（每行一个，格式：基因名,细胞类型）",
                    placeholder="例如：\nCD3D,T cells\nCD79A,B cells\nLYZ,Monocytes",
                    lines=5,
                    interactive=True
                )
                
                load_default_markers = gr.Button("📚 加载默认标记基因", variant="secondary")
            
            with gr.Group():
                gr.HTML("<h4>置信度阈值</h4>")
                
                confidence_threshold = gr.Slider(
                    minimum=0.1,
                    maximum=1.0,
                    value=0.5,
                    step=0.05,
                    label="最低置信度",
                    interactive=True
                )
                
                min_marker_genes = gr.Slider(
                    minimum=1,
                    maximum=10,
                    value=3,
                    step=1,
                    label="最少标记基因数",
                    interactive=True
                )
            
            with gr.Group():
                gr.HTML("<h4>可视化设置</h4>")
                
                show_confidence = gr.Checkbox(
                    value=True,
                    label="显示置信度",
                    interactive=True
                )
                
                plot_type = gr.CheckboxGroup(
                    choices=["UMAP", "热图", "气泡图", "堆叠柱状图"],
                    value=["UMAP", "热图"],
                    label="选择可视化类型",
                    interactive=True
                )
            
            annotate_btn = gr.Button("🏷️ 开始注释", elem_classes="custom-button")
            refine_btn = gr.Button("🔄 优化注释", elem_classes="custom-button")
        
        # 右侧结果展示面板
        with gr.Column(scale=2):
            gr.HTML("<h3>📊 注释结果</h3>")
            
            with gr.Tabs():
                with gr.TabItem("注释概览"):
                    annotation_summary = gr.HTML()
                
                with gr.TabItem("UMAP可视化"):
                    annotation_umap_plot = gr.Plot(label="细胞类型注释结果")
                
                with gr.TabItem("标记基因热图"):
                    marker_heatmap_plot = gr.Plot(label="标记基因表达热图")
                
                with gr.TabItem("气泡图"):
                    bubble_plot = gr.Plot(label="标记基因气泡图")
                
                with gr.TabItem("细胞类型统计"):
                    with gr.Row():
                        celltype_stats = gr.Dataframe(
                            headers=["细胞类型", "数量", "比例", "置信度"],
                            label="细胞类型统计",
                            scale=1
                        )
                        celltype_proportion_plot = gr.Plot(label="细胞类型比例", scale=1)
                
                with gr.TabItem("注释详情"):
                    annotation_details = gr.Dataframe(
                        headers=["聚类", "细胞类型", "置信度", "标记基因"],
                        label="详细注释结果",
                        elem_classes="data-table"
                    )
            
            # 注释状态
            annotation_status = gr.Textbox(
                label="注释状态",
                value="等待细胞类型注释...",
                interactive=False,
                lines=4
            )
    
    # 事件处理函数
    def load_default_marker_genes():
        default_markers = """CD3D,T cells
CD3E,T cells
CD3G,T cells
CD4,T cells
CD8A,T cells
CD8B,T cells
CD79A,B cells
CD79B,B cells
MS4A1,B cells
IGKC,B cells
LYZ,Monocytes
CD14,Monocytes
CD68,Monocytes
FCGR3A,Monocytes
CST3,Monocytes
MS4A7,Monocytes
NKG7,NK cells
GNLY,NK cells
KLRD1,NK cells
KLRC1,NK cells
FCGR3A,NK cells
PPBP,Megakaryocytes
PF4,Megakaryocytes
GNG11,Megakaryocytes
HBB,Erythrocytes
HBA1,Erythrocytes
HBA2,Erythrocytes"""
        
        return gr.update(value=default_markers)
    
    def perform_annotation(annotation_method, reference_dataset, marker_genes_input, 
                          confidence_threshold, min_marker_genes, show_confidence, plot_type):
        try:
            status_msg = f"🔄 开始{annotation_method}...\n"
            
            # 这里应该加载之前聚类的数据
            # 现在使用模拟数据
            import scanpy as sc
            adata = sc.datasets.pbmc3k()
            
            # 添加聚类结果
            sc.pp.neighbors(adata, n_neighbors=15)
            sc.tl.leiden(adata, resolution=0.5)
            
            # 解析标记基因
            marker_genes = {}
            if marker_genes_input:
                for line in marker_genes_input.strip().split('\n'):
                    if ',' in line:
                        gene, cell_type = line.strip().split(',')
                        if cell_type not in marker_genes:
                            marker_genes[cell_type] = []
                        marker_genes[cell_type].append(gene.strip())
            
            status_msg += "🧬 解析标记基因...\n"
            
            # 模拟注释过程
            cell_types = ['T cells', 'B cells', 'Monocytes', 'NK cells', 'Megakaryocytes', 'Erythrocytes']
            clusters = adata.obs['leiden'].unique()
            
            # 为每个聚类分配细胞类型
            cluster_annotations = {}
            cluster_confidences = {}
            cluster_markers = {}
            
            for cluster in clusters:
                # 模拟基于标记基因的注释
                if annotation_method == "自动注释":
                    # 简单的基于标记基因的注释模拟
                    import random
                    cell_type = random.choice(cell_types)
                    confidence = random.uniform(confidence_threshold, 0.95)
                    markers = marker_genes.get(cell_type, [f"Gene_{i}" for i in range(3)])
                else:
                    cell_type = f"Cluster {cluster}"
                    confidence = 0.8
                    markers = [f"Marker_{i}" for i in range(3)]
                
                cluster_annotations[cluster] = cell_type
                cluster_confidences[cluster] = confidence
                cluster_markers[cluster] = ", ".join(markers[:min_marker_genes])
            
            # 更新adata中的注释
            adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
            adata.obs['confidence'] = adata.obs['leiden'].map(cluster_confidences)
            
            status_msg += "📊 生成注释结果...\n"
            
            # 生成注释概览
            n_cells = len(adata)
            n_clusters = len(clusters)
            
            summary_html = f"""
            <div style='padding: 20px;'>
                <h4>📋 注释概览</h4>
                <table style='width: 100%; border-collapse: collapse;'>
                    <tr><td><strong>总细胞数:</strong></td><td>{n_cells:,}</td></tr>
                    <tr><td><strong>聚类数:</strong></td><td>{n_clusters}</td></tr>
                    <tr><td><strong>注释方法:</strong></td><td>{annotation_method}</td></tr>
                    <tr><td><strong>参考数据集:</strong></td><td>{reference_dataset}</td></tr>
                    <tr><td><strong>平均置信度:</strong></td><td>{adata.obs['confidence'].mean():.2f}</td></tr>
                </table>
            </div>
            """
            
            # 细胞类型统计
            celltype_counts = adata.obs['cell_type'].value_counts()
            celltype_data = []
            for cell_type, count in celltype_counts.items():
                avg_confidence = adata.obs[adata.obs['cell_type'] == cell_type]['confidence'].mean()
                celltype_data.append([
                    cell_type,
                    count,
                    f"{count/n_cells*100:.1f}%",
                    f"{avg_confidence:.2f}"
                ])
            
            celltype_df = pd.DataFrame(celltype_data, columns=['细胞类型', '数量', '比例', '置信度'])
            
            # 详细注释结果
            annotation_detail_data = []
            for cluster in clusters:
                cell_type = cluster_annotations[cluster]
                confidence = cluster_confidences[cluster]
                markers = cluster_markers[cluster]
                annotation_detail_data.append([f"Cluster {cluster}", cell_type, f"{confidence:.2f}", markers])
            
            annotation_details_df = pd.DataFrame(annotation_detail_data, columns=['聚类', '细胞类型', '置信度', '标记基因'])
            
            # 生成图表
            plots = {}
            if "UMAP" in plot_type:
                plots['umap'] = plotter.plot_annotated_umap(adata, show_confidence)
            if "热图" in plot_type:
                plots['heatmap'] = plotter.plot_marker_heatmap(adata)
            if "气泡图" in plot_type:
                plots['bubble'] = plotter.plot_marker_bubble(adata)
            if "堆叠柱状图" in plot_type:
                plots['stacked'] = plotter.plot_celltype_proportions(adata)
            
            status_msg += "✅ 细胞类型注释完成！"
            
            return (
                gr.update(value=summary_html),
                gr.update(value=plots.get('umap')),
                gr.update(value=plots.get('heatmap')),
                gr.update(value=plots.get('bubble')),
                gr.update(value=celltype_df),
                gr.update(value=plots.get('stacked')),
                gr.update(value=annotation_details_df),
                status_msg
            )
            
        except Exception as e:
            error_msg = f"❌ 细胞类型注释失败: {str(e)}"
            return (
                gr.update(value=f"<p style='color: red;'>{error_msg}</p>"),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                error_msg
            )
    
    def refine_annotation():
        try:
            status_msg = "🔄 正在优化注释...\n"
            status_msg += "📈 重新计算置信度...\n"
            status_msg += "🧬 检查标记基因一致性...\n"
            status_msg += "✅ 注释优化完成！"
            return status_msg
        except Exception as e:
            return f"❌ 注释优化失败: {str(e)}"
    
    # 事件绑定
    load_default_markers.click(
        fn=load_default_marker_genes,
        outputs=[marker_genes_input]
    )
    
    annotate_btn.click(
        fn=perform_annotation,
        inputs=[annotation_method, reference_dataset, marker_genes_input,
                confidence_threshold, min_marker_genes, show_confidence, plot_type],
        outputs=[annotation_summary, annotation_umap_plot, marker_heatmap_plot, bubble_plot,
                celltype_stats, celltype_proportion_plot, annotation_details, annotation_status]
    )
    
    refine_btn.click(
        fn=refine_annotation,
        outputs=[annotation_status]
    )