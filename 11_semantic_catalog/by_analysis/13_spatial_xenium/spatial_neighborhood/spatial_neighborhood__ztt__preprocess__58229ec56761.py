import gradio as gr
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from functions.plotter import Plotter
from functions.xenium_analysis import data_preprocess

def create_status_html(status_type, title, content=""):
    """创建状态HTML的通用函数（使用全局CSS的卡片样式，含图标）"""
    status_class_map = {
        'success': 'status-success',
        'error': 'status-error',
        'waiting': 'status-wait',
        'info': 'status-wait',
        'clustering': 'status-success'
    }
    icon_map = {
        'success': '✅',
        'error': '❌',
        'waiting': '⏳',
        'info': 'ℹ️',
        'clustering': '✅'
    }
    cls = status_class_map.get(status_type, 'status-wait')
    icon = icon_map.get(status_type, '⏳')
    return f"""
    <div class='status-card {cls}'>
        <span class='status-icon'>{icon}</span>
        <div class='status-content'>
            <h4>{title}</h4>
            {f"<p>{content}</p>" if content else ''}
        </div>
    </div>
    """

# 前端包装函数：调用纯计算 data_preprocess，并生成HTML与图表
def apply_cell_filter_ui(
    adata,
    filter_params_df,
    analysis_options,
    hvgs,
    npcs,
    n_neighbors,
    resolution,
    feature_species,
    feature_method,
    feature_top,
    feature_lfc
):
    result = data_preprocess(
        adata,
        filter_params_df,
        analysis_options,
        hvgs,
        npcs,
        n_neighbors,
        resolution,
        feature_species,
        feature_method,
        feature_top,
        feature_lfc
    )
    # 错误处理
    if result is None or ('error' in result and result['error']):
        err = result['error'] if isinstance(result, dict) else '未知错误'
        error_html = create_status_html('error', '❌ 分析失败', err)
        return None, error_html, None, error_html, None, None, None, pd.DataFrame(columns=["group", "gene", "score", "logfoldchanges", "pvals_adj"]) 

    # 过滤HTML与QC图
    filter_html = ""
    filter_qc_plot = None
    adata_processed = None
    if result.get('filter'):
        f = result['filter']
        adata_processed = f.get('adata_processed')
        sample_details = ""
        for stat in f.get('sample_stats', []):
            sample_details += f"""
            <tr>
                <td style='padding: 5px; border: 1px solid #ddd;'>{stat['sample']}</td>
                <td style='padding: 5px; border: 1px solid #ddd;'>{stat['original']:,}</td>
                <td style='padding: 5px; border: 1px solid #ddd;'>{stat['filtered']:,}</td>
                <td style='padding: 5px; border: 1px solid #ddd;'>{stat['removed']:,}</td>
                <td style='padding: 5px; border: 1px solid #ddd;'>{stat['removal_rate']:.1f}%</td>
            </tr>
            """
        filter_html = f"""
        <div class='result-card result-success'>
            <h5>✅ 细胞过滤完成</h5>
            <div class='result-section-lg'>
                <strong>总体统计：</strong>
                原始细胞数: {f['original_cells']:,} → 过滤后: {f['filtered_cells']:,} (移除 {f['removed_cells']:,}, {f['removal_rate']:.1f}%)
            </div>
            <div class='result-section-sm'>
                <strong>各样本过滤详情：</strong>
            </div>
            <table class='result-table'>
                <thead>
                    <tr>
                        <th>样本</th>
                        <th class='text-center'>原始细胞</th>
                        <th class='text-center'>过滤后</th>
                        <th class='text-center'>移除</th>
                        <th class='text-center'>移除率</th>
                    </tr>
                </thead>
                <tbody>
                    {sample_details}
                </tbody>
            </table>
        </div>
        """
        # QC 图（Plotter 小提琴图）
        try:
            filter_qc_plot = Plotter().plot_qc_violin(adata_processed)
        except Exception:
            filter_qc_plot = None

    # 聚类HTML与UMAP图片
    clustering_html = ""
    fig_umap = None
    fig_umap_group = None
    adata_cluster = result.get('cluster')
    if adata_cluster is not None and ("降维聚类" in analysis_options):
        params = result.get('params', {})
        clustering_info = f"高变基因: {params.get('hvgs')}, 主成分: {params.get('npcs')}, 邻居数: {params.get('n_neighbors')}, 分辨率: {params.get('resolution')}<br>批次效应校正: harmony"
        clustering_html = f"""
        <div class='result-card result-clustering'>
            <h5>✅ 降维聚类完成</h5>
            <p>{clustering_info}</p>
        </div>
        """
        try:
            # 若UMAP未计算则确保可用
            has_umap = hasattr(adata_cluster, "obsm") and ("X_umap" in adata_cluster.obsm.keys())
            if not has_umap:
                sc.pp.neighbors(adata_cluster)
                sc.tl.umap(adata_cluster)
            # 总体UMAP
            sc.pl.umap(adata_cluster, color=None, show=False)
            fig_umap = plt.gcf(); fig_umap.set_dpi(200); fig_umap.tight_layout()
            # 分组UMAP
            group_col = None
            for col in adata_cluster.obs.columns:
                if col.startswith('res_'):
                    group_col = col; break
            if group_col is None:
                group_col = 'leiden' if 'leiden' in adata_cluster.obs.columns else None
            if group_col:
                sc.pl.umap(adata_cluster, color=group_col, legend_loc='on data', show=False)
                fig_umap_group = plt.gcf(); fig_umap_group.set_dpi(200); fig_umap_group.tight_layout()
        except Exception:
            fig_umap = None
            fig_umap_group = None

    # 特征基因 DotPlot 与表
    feature_dotplot = None
    feature_table_df = pd.DataFrame(columns=["group", "gene", "score", "logfoldchanges", "pvals_adj"]) 
    marker_res = result.get('marker')
    if marker_res is not None and ("特征基因筛选" in analysis_options):
        try:
            sc.pl.rank_genes_groups_dotplot(
                marker_res['adata_for_marker'],
                groupby=marker_res['groupby_key'],
                standard_scale="var",
                n_genes=int(feature_top or 5),
                show=False,
                figsize=(20, 5),
            )
            feature_dotplot = plt.gcf(); feature_dotplot.set_dpi(200); feature_dotplot.tight_layout()
        except Exception:
            feature_dotplot = None
        try:
            feature_table_df = marker_res.get('feature_table_df', feature_table_df)
        except Exception:
            pass

    # 输出顺序需与页面绑定一致
    return adata_cluster, filter_html, filter_qc_plot, clustering_html, fig_umap, fig_umap_group, feature_dotplot, feature_table_df

def data_preprocess_tab(adata_raw, adata_cluster, output_dir_state):
    # 从adata_raw生成过滤参数DataFrame
    def generate_filter_parameters(adata):
        if adata is not None:
            try:
                dic = {}
                for index, sample in enumerate(adata.obs['sample'].drop_duplicates()):
                    tmp = adata[adata.obs['sample'] == sample]
                    min_genes = int(np.quantile(tmp.obs['n_genes_by_counts'], 0.05))
                    max_genes = int(np.quantile(tmp.obs['n_genes_by_counts'], 1))
                    dic.setdefault('sample', []).append(sample)
                    dic.setdefault('min_genes', []).append(min_genes)
                    dic.setdefault('max_genes', []).append(max_genes)
                    dic.setdefault('min_GenesPerUMI', []).append(0.7)
                    dic.setdefault('max_GenesPerUMI', []).append(1)
                
                df = pd.DataFrame(dic)
                status_html = create_status_html('success', '参数已自动生成')
                return gr.update(value=df, visible=True), gr.update(value=status_html), gr.update(interactive=True)
                
            except Exception as e:
                error_html = create_status_html('error', '❌ 参数生成失败', str(e))
                return gr.update(value=pd.DataFrame(), visible=False), gr.update(value=error_html), gr.update(interactive=False)
        else:
            waiting_html = create_status_html('waiting', '📋 等待数据加载', '请先在"数据选择"页面加载数据，系统将自动生成过滤参数')
            return gr.update(value=pd.DataFrame(), visible=False), gr.update(value=waiting_html), gr.update(interactive=False)
    
    # 细胞过滤函数
    # 细胞过滤/降维聚类/特征基因筛选逻辑已抽离至 functions.xenium_analysis.data_preprocess
    
    with gr.Row(equal_height=True):
        # 左侧参数面板（优化为分组卡片 + 标题样式）
        with gr.Column(scale=1):

            # 过滤参数设置卡片
            with gr.Group(elem_classes=["card-section"]):
                gr.HTML("<h4>📊 过滤参数设置</h4>")

                # 状态显示
                status_display = gr.HTML(
                    value=create_status_html('waiting', '📋 等待数据加载', '请先在"数据选择"页面加载数据')
                )

                # 可编辑的参数DataFrame（切换为紧凑表格样式）
                filter_params = gr.Dataframe(
                    headers=["sample", "min_genes", "max_genes", "min_GenesPerUMI", "max_GenesPerUMI"],
                    datatype=["str", "number", "number", "number", "number"],
                    value=pd.DataFrame(),
                    interactive=True,
                    visible=False,
                    label="过滤参数 (可编辑)",
                    elem_classes=["data-table-compact"]
                )

            # 分析功能选择卡片
            with gr.Group(elem_classes=["card-section"]):
                gr.HTML("<h4>🔧 分析功能选择</h4>")
                analysis_options = gr.CheckboxGroup(
                    choices=["细胞过滤", "降维聚类", "特征基因筛选"],
                    value=["细胞过滤"],
                    label="选择要执行的分析步骤",
                    interactive=True
                )

                # 将“降维聚类参数”与“特征基因筛选参数”并入该卡片
                with gr.Group(visible=False) as dim_params_group:
                    gr.HTML("<h5>📊 降维聚类参数</h5>")

                    with gr.Row():
                        hvgs_param = gr.Number(
                            value=2000,
                            label="高变基因数量 (hvgs)",
                            minimum=500,
                            maximum=5000,
                            step=100,
                            interactive=True
                        )
                        npcs_param = gr.Number(
                            value=30,
                            label="主成分数量 (npcs)",
                            minimum=10,
                            maximum=100,
                            step=5,
                            interactive=True
                        )

                    with gr.Row():
                        n_neighbors_param = gr.Number(
                            value=20,
                            label="邻居数量 (n_neighbors)",
                            minimum=5,
                            maximum=100,
                            step=5,
                            interactive=True
                        )
                        resolution_param = gr.Number(
                            value=1.0,
                            label="聚类分辨率 (resolution)",
                            minimum=0.1,
                            maximum=3.0,
                            step=0.1,
                            interactive=True
                        )

                    # 并入：特征基因筛选参数
                    with gr.Group(visible=False) as feature_params_box:
                        gr.HTML("<h5>🧬 特征基因筛选参数</h5>")
                        with gr.Row():
                            feature_species_dd = gr.Dropdown(
                                choices=["Human", "Mouse"],
                                value="Human",
                                label="物种 (Species)",
                                interactive=True
                            )
                            feature_method_dd = gr.Dropdown(
                                choices=["t-test", "wilcoxon"],
                                value="t-test",
                                label="统计方法",
                                interactive=True
                            )
                            top_genes_param = gr.Number(
                                value=5,
                                label="Top基因数",
                                minimum=1,
                                maximum=50,
                                step=1,
                                interactive=True
                            )
                            feature_lfc_param = gr.Number(
                                value=0.25,
                                label="差异倍数阈值 |logFC|",
                                minimum=0.0,
                                maximum=5.0,
                                step=0.05,
                                interactive=True
                            )


            # 操作按钮卡片（与上方保持一致的卡片风格）
            with gr.Group(elem_classes=["card-section"]):
                apply_filter_btn = gr.Button(
                    "🚀 执行分析", 
                    elem_classes="custom-button",
                    interactive=False
                )


        # 右侧结果展示面板
        with gr.Column(scale=2, elem_classes=["result-panel-compact"]):
            with gr.Tabs() as result_tabs:
                with gr.TabItem("过滤结果", visible=True) as filter_tab:
                    filter_result = gr.HTML(
                        value=create_status_html('waiting', '📋 等待分析操作', '设置参数后，点击"执行分析"查看结果')
                    )
                    
                    # 过滤后的质控图表
                    filter_qc_plots = gr.Plot(
                        label="过滤后质控图表",
                        elem_classes=["plot-compact"],
                        visible=True
                    )
                
                with gr.TabItem("聚类结果", visible=False) as clustering_tab:
                    clustering_result = gr.HTML(
                        value=create_status_html('waiting', '📊 等待聚类分析', '选择"降维聚类"并执行分析后查看结果')
                    )
                    
                    # 聚类可视化图表（多图显示）
                    with gr.Row():
                        with gr.Column(scale=3):
                            cluster_umap_plot = gr.Plot(
                                label="UMAP总体",
                                elem_classes=["plot-compact"],
                                visible=True
                            )
                        # 右侧控制面板：点大小与图例位置
                        with gr.Column(scale=1):
                            gr.HTML("<div style='text-align:right'><strong>绘图控制</strong></div>")
                            umap_point_size = gr.Slider(
                                minimum=1.0,
                                maximum=10.0,
                                value=5.0,
                                step=0.5,
                                label="UMAP点大小",
                                interactive=True
                            )
                            umap_legend_loc = gr.Dropdown(
                                choices=["on data", "right margin", "none"],
                                value="on data",
                                label="图例位置",
                                interactive=True
                            )
                    with gr.Row():
                        cluster_group_umap_plot = gr.Plot(
                            label="UMAP按group拆分",
                            elem_classes=["plot-compact"],
                            visible=True
                        )

                # 新增：特征基因结果 Tab
                with gr.TabItem("特征基因", visible=False) as feature_tab:
                    with gr.Row():
                        feature_dotplot_plot = gr.Plot(
                            label="Rank Genes Groups DotPlot",
                            elem_classes=["plot-compact"],
                            visible=True
                        )
                    feature_genes_table = gr.Dataframe(
                        headers=["group", "gene", "score", "logfoldchanges", "pvals_adj"],
                        datatype=["str", "str", "number", "number", "number"],
                        row_count=(10, "dynamic"),
                        col_count=5,
                        interactive=False,
                        label="每簇 Top 基因列表"
                    )
    
    # 条件显示降维聚类参数的函数（现已嵌入“分析功能选择”卡片）
    def toggle_dim_params(analysis_options):
        """根据分析选项显示/隐藏降维聚类与其内部参数"""
        show_dim_params = "降维聚类" in analysis_options
        return gr.Group(visible=show_dim_params)

    # 条件显示特征基因筛选参数的函数（嵌入同卡片）
    def toggle_feature_params(analysis_options):
        show_feature_params = "特征基因筛选" in analysis_options
        return gr.Group(visible=show_feature_params)
    
    # 动态显示结果tab的函数
    def toggle_result_tabs(analysis_options):
        """根据分析选项显示/隐藏结果tab"""
        show_filter_tab = "细胞过滤" in analysis_options
        show_clustering_tab = "降维聚类" in analysis_options
        show_feature_tab = "特征基因筛选" in analysis_options
        return (
            gr.TabItem(visible=show_filter_tab),
            gr.TabItem(visible=show_clustering_tab),
            gr.TabItem(visible=show_feature_tab)
        )

    # UMAP重绘函数：根据点大小与图例位置更新UMAP总体与分组图
    def redraw_umap(adata_cluster_state, point_size, legend_loc_str):
        try:
            if adata_cluster_state is None:
                fig_err = plt.figure(figsize=(6, 2))
                plt.text(0.5, 0.5, "未检测到聚类结果，请先执行降维聚类", ha="center", va="center")
                plt.axis("off")
                return gr.update(value=fig_err), gr.update(value=fig_err)

            # 将legend位置字符串转换为参数（scanpy支持'on data'/'right margin'/'none'）
            legend_loc = 'none' if str(legend_loc_str).lower() == 'none' else str(legend_loc_str)

            # 检测聚类列
            def find_cluster_key(adata_obj):
                try:
                    obs_df = getattr(adata_obj, 'obs')
                    cols = list(obs_df.columns)
                    candidates = [
                        "leiden", "louvain", "res", "cluster", "clusters",
                        "seurat_clusters", "kmeans", "phenograph"
                    ]
                    for c in candidates:
                        if c in cols:
                            return c
                    lowered = [str(c).lower() for c in cols]
                    for idx, col in enumerate(lowered):
                        for kw in candidates:
                            if kw in col:
                                return cols[idx]
                    viable = []
                    for col in cols:
                        try:
                            dtype_str = str(obs_df[col].dtype).lower()
                            n_groups = int(obs_df[col].nunique())
                            if n_groups >= 2 and ("category" in dtype_str or "object" in dtype_str):
                                viable.append((col, n_groups))
                        except Exception:
                            continue
                    if viable:
                        viable.sort(key=lambda x: x[1], reverse=True)
                        return viable[0][0]
                except Exception:
                    pass
                return None

            groupby_key = find_cluster_key(adata_cluster_state)
            if groupby_key is None:
                fig_err = plt.figure(figsize=(6, 2))
                plt.text(0.5, 0.5, "未找到聚类或注释列，无法绘制UMAP", ha="center", va="center")
                plt.axis("off")
                return gr.update(value=fig_err), gr.update(value=fig_err)

            # 若UMAP未计算则尝试计算
            try:
                has_umap = hasattr(adata_cluster_state, "obsm") and ("X_umap" in adata_cluster_state.obsm.keys())
                if not has_umap:
                    sc.pp.neighbors(adata_cluster_state)
                    sc.tl.umap(adata_cluster_state)
            except Exception:
                pass

            # 绘制UMAP总体
            fig_umap, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
            sc.pl.umap(
                adata_cluster_state,
                color=[groupby_key],
                title=['Cluster'],
                size=float(point_size),
                legend_loc=legend_loc,
                legend_fontsize=10,
                legend_fontoutline=2,
                frameon=False,
                show=False,
                ax=ax
            )

            # 按group拆分绘制
            fig_group = None
            try:
                groups = list(adata_cluster_state.obs.get('group', pd.Series()).dropna().drop_duplicates())
            except Exception:
                groups = []
            if groups and len(groups) > 1:
                fig_group, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
                for i, g in enumerate(groups):
                    sc.pl.umap(
                        adata_cluster_state[adata_cluster_state.obs["group"]==g],
                        color=groupby_key,
                        title=str(g),
                        ax=axs[i],
                        size=float(point_size),
                        legend_loc=legend_loc or 'on data',
                        legend_fontsize=12,
                        legend_fontoutline=2,
                        frameon=False,
                        show=False
                    )
            else:
                fig_group = plt.figure(figsize=(6, 2))
                plt.text(0.5, 0.5, "按group拆分不可用（缺少'group'或只有一个组）", ha="center", va="center")
                plt.axis("off")

            return gr.update(value=fig_umap), gr.update(value=fig_group)
        except Exception as e:
            fig_err = plt.figure(figsize=(6, 2))
            plt.text(0.5, 0.5, f"重绘失败: {e}", ha="center", va="center")
            plt.axis("off")
            return gr.update(value=fig_err), gr.update(value=fig_err)

    # 计算并更新UMAP点大小：基于聚类后的adata细胞数
    def calc_umap_point_size(adata_cluster_state):
        try:
            n_cells = int(getattr(adata_cluster_state, 'n_obs', adata_cluster_state.obs.shape[0]))
            base = 120000.0 / float(n_cells)
            size = base * 5.0 if base < 1.0 else base
            size = max(1.0, min(10.0, float(size)))
            return gr.update(value=size)
        except Exception:
            return gr.update(value=5.0)
    
    # 事件绑定
    
    # 当adata_raw变化时，自动生成过滤参数
    adata_raw.change(
        fn=generate_filter_parameters,
        inputs=[adata_raw],
        outputs=[filter_params, status_display, apply_filter_btn]
    )
    
    # 当分析选项变化时，显示/隐藏降维聚类参数和结果tab
    analysis_options.change(
        fn=toggle_dim_params,
        inputs=[analysis_options],
        outputs=[dim_params_group]
    )

    # 当分析选项变化时，显示/隐藏特征基因筛选参数
    analysis_options.change(
        fn=toggle_feature_params,
        inputs=[analysis_options],
        outputs=[feature_params_box]
    )
    
    # 当分析选项变化时，动态显示/隐藏结果tab
    analysis_options.change(
        fn=toggle_result_tabs,
        inputs=[analysis_options],
        outputs=[filter_tab, clustering_tab, feature_tab]
    )
    
    # 应用分析按钮点击事件
    apply_filter_btn.click(
        fn=apply_cell_filter_ui,
        inputs=[adata_raw, filter_params, analysis_options, hvgs_param, npcs_param, n_neighbors_param, resolution_param, feature_species_dd, feature_method_dd, top_genes_param, feature_lfc_param],
        # 输出：保存聚类后的adata到全局状态、过滤结果HTML、过滤QC图、聚类结果HTML、UMAP总体图、UMAP分组图、特征基因DotPlot、Top表
        outputs=[adata_cluster, filter_result, filter_qc_plots, clustering_result, cluster_umap_plot, cluster_group_umap_plot, feature_dotplot_plot, feature_genes_table]
    )

    # 聚类完成后，根据细胞数自动设定UMAP点大小滑块的值
    apply_filter_btn.click(
        fn=calc_umap_point_size,
        inputs=[adata_cluster],
        outputs=[umap_point_size]
    )

    # 绑定控件变化即重绘UMAP
    umap_point_size.change(
        fn=redraw_umap,
        inputs=[adata_cluster, umap_point_size, umap_legend_loc],
        outputs=[cluster_umap_plot, cluster_group_umap_plot]
    )
    umap_legend_loc.change(
        fn=redraw_umap,
        inputs=[adata_cluster, umap_point_size, umap_legend_loc],
        outputs=[cluster_umap_plot, cluster_group_umap_plot]
    )