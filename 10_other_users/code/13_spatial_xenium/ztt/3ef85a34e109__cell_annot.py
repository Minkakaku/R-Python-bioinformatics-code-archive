import gradio as gr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


def cell_annot_tab(adata_cluster, output_dir_state):
    """
    手动细胞类型注释页面（优化版）：
    - 增加物种下拉（人/小鼠），展示预设的可编辑 marker 表；点击按钮基于表生成表达 UMAP 与 DotPlot。
    - 右侧显示表达图与聚类命名预览；聚类命名逻辑保持不变。
    """

    def init_cluster_map(adata):
        # 根据传入的 adata_cluster 尝试自动填充聚类ID列表
        try:
            if adata is None:
                return pd.DataFrame(columns=["Cluster", "名称"])
            obs = getattr(adata, "obs", None)
            if obs is None or len(getattr(obs, "columns", [])) == 0:
                return pd.DataFrame(columns=["Cluster", "名称"])
            candidates = ["leiden", "louvain", "res", "cluster", "clusters"]
            col = next((c for c in candidates if c in obs.columns), None)
            if col is None:
                return pd.DataFrame(columns=["Cluster", "名称"])
            counts = obs[col].astype(str).value_counts().sort_index()
            return pd.DataFrame({"Cluster": counts.index, "名称": [""] * len(counts)})
        except Exception:
            return pd.DataFrame(columns=["Cluster", "名称"])

    def find_cluster_key(adata):
        try:
            if adata is None:
                return None
            obs = getattr(adata, "obs", None)
            if obs is None:
                return None
            candidates = ["leiden", "louvain", "res", "cluster", "clusters"]
            return next((c for c in candidates if c in obs.columns), None)
        except Exception:
            return None

    # 物种预设 marker 表
    PRESET_MARKERS = {
        "人": [
            ["T细胞", "CD3D, CD3E, TRAC"],
            ["B细胞", "MS4A1, CD79A"],
            ["巨噬细胞", "LYZ"],
            ["成纤维细胞", "COL1A1, COL1A2"],
            ["内皮细胞", "KDR, PECAM1"],
        ],
        "小鼠": [
            ["T细胞", "Cd3d, Cd3e, Trac"],
            ["B细胞", "Ms4a1, Cd79a"],
            ["巨噬细胞", "Lyz2"],
            ["成纤维细胞", "Col1a1, Col1a2"],
            ["内皮细胞", "Kdr, Pecam1"],
        ],
    }

    def preset_table(species):
        return PRESET_MARKERS.get(species or "人", PRESET_MARKERS["人"])[:]

    initial_cluster_df = init_cluster_map(adata_cluster)

    # 工具函数：获取可选的 groupby 列
    def get_groupby_choices(adata):
        try:
            obs = getattr(adata, "obs", None)
            if obs is None:
                return []
            candidates = ["leiden", "louvain", "res", "cluster", "clusters"]
            return [c for c in candidates if c in obs.columns]
        except Exception:
            return []

    groupby_default = find_cluster_key(adata_cluster)
    groupby_choices = get_groupby_choices(adata_cluster)

    # 两个 Accordion，每个内部采用左右两列：左参数设置，右结果展示

    # 2) 基因表达
    with gr.Accordion("🧪 基因表达", open=True):
        with gr.Row(equal_height=True):
            # 左列：参数
            with gr.Column(scale=1):
                species_dd = gr.Dropdown(choices=["人", "小鼠"], value="人", label="物种")
                marker_table = gr.Dataframe(
                    headers=["细胞类型", "基因列表"],
                    datatype=["str", "str"],
                    row_count=(5, "dynamic"),
                    col_count=2,
                    value=preset_table("人"),
                    label="细胞类型与标记基因（逗号分隔，可编辑）",
                    interactive=True,
                )
                gen_expr_btn = gr.Button("🎯 生成表达 UMAP 与 DotPlot", elem_classes="custom-button")
            # 右列：结果
            with gr.Column(scale=2):
                with gr.Tabs():
                    with gr.TabItem("UMAP"):
                        expr_umap_fig = gr.Plot(label="UMAP 基因表达")
                    with gr.TabItem("DotPlot"):
                        expr_dotplot_fig = gr.Plot(label="DotPlot 基因表达")
                # 状态文本框已移除，结果直接在图中展示

    # 3) 细胞注释
    with gr.Accordion("📝 细胞注释", open=True):
        with gr.Row(equal_height=True):
            # 左列：参数
            with gr.Column(scale=1):
                cluster_map_table = gr.Dataframe(
                    headers=["Cluster", "名称"],
                    datatype=["str", "str"],
                    row_count=(8, "dynamic"),
                    col_count=2,
                    value=initial_cluster_df if not initial_cluster_df.empty else [["0", ""], ["1", ""], ["2", ""], ["3", ""]],
                    label="为各聚类指定细胞类型名称",
                    interactive=True,
                )
                submit_btn = gr.Button("✅ 提交注释并预览", elem_classes="custom-button")
            # 右列：结果
            with gr.Column(scale=2):
                with gr.Tabs():
                    with gr.TabItem("UMAP"):
                        annot_umap_fig = gr.Plot(label="UMAP（按注释着色）")
                    with gr.TabItem("空间图"):
                        annot_spatial_fig = gr.Plot(label="空间图（如有坐标）")
                    with gr.TabItem("映射表"):
                        mapping_preview = gr.Dataframe(
                            headers=["Cluster", "名称"],
                            datatype=["str", "str"],
                            row_count=(8, "dynamic"),
                            col_count=2,
                            label="聚类命名预览",
                            interactive=False,
                        )
                # 状态文本框已移除，结果直接在图与表中展示

    # --- 事件处理（接入 Scanpy，不写入全局状态） ---

    def generate_expression_plots(marker_rows):
        try:
            if adata_cluster is None:
                raise ValueError("全局变量 adata_cluster 为空，请先完成聚类步骤")
            groupby = find_cluster_key(adata_cluster)
            if groupby is None:
                raise ValueError("未找到聚类列（如 'leiden'/'louvain'），请先完成聚类步骤以生成聚类列")

            df = pd.DataFrame(marker_rows, columns=["细胞类型", "基因列表"]) if marker_rows is not None else pd.DataFrame()
            df = df.dropna(how="all")
            if df.empty:
                raise ValueError("请至少填写一行细胞类型及其基因列表")

            # genes 用于 UMAP 表达上色；var_names 用于 DotPlot
            var_names = {}
            genes_all = []
            for _, row in df.iterrows():
                ct = str(row[0]).strip()
                genes = [g.strip() for g in str(row[1]).split(",") if str(row[1]) and g.strip()]
                if ct and len(genes) > 0:
                    var_names[ct] = genes
                    genes_all.extend(genes)
            if len(var_names) == 0:
                raise ValueError("请填写有效的细胞类型和至少一个基因")

            # 过滤不存在的基因
            genes_all = [g for g in genes_all if g in list(getattr(adata_cluster, "var_names", []))]
            if len(genes_all) == 0:
                raise ValueError("所填基因均不在对象中，请检查基因名与物种大小写")

            # UMAP 若未计算则尝试计算
            try:
                has_umap = hasattr(adata_cluster, "obsm") and ("X_umap" in adata_cluster.obsm.keys())
                if not has_umap:
                    sc.pp.neighbors(adata_cluster)
                    sc.tl.umap(adata_cluster)
            except Exception:
                pass

            # 生成 UMAP（多基因颜色）
            sc.pl.umap(adata_cluster, color=genes_all, show=False)
            fig_umap = plt.gcf()

            # 生成 DotPlot（按 groupby 分组）
            dp = sc.pl.dotplot(
                adata_cluster,
                var_names=var_names,
                groupby=groupby,
                standard_scale="var",
                show=False,
            )
            fig_dp = getattr(dp, "figure", plt.gcf())

            return gr.update(value=fig_umap), gr.update(value=fig_dp)
        except Exception as e:
            fig = plt.figure(figsize=(6, 2))
            plt.text(0.5, 0.5, f"表达绘制失败: {e}", ha="center", va="center")
            plt.axis("off")
            return gr.update(value=fig), gr.update(value=fig)

    def submit_annotation(mapping_rows):
        try:
            df = pd.DataFrame(mapping_rows, columns=["Cluster", "名称"]) if mapping_rows is not None else pd.DataFrame(columns=["Cluster", "名称"])
            df = df.dropna(how="all")
            if df.empty:
                raise ValueError("请填写至少一个聚类名称")
            groupby = find_cluster_key(adata_cluster)
            if groupby is None:
                raise ValueError("未找到聚类列（如 'leiden'/'louvain'），请先完成聚类步骤以生成聚类列")

            # 生成预览：复制对象并添加注释列
            try:
                ad = adata_cluster.copy()
            except Exception:
                ad = adata_cluster
            mapping = {str(r["Cluster"]): str(r["名称"]).strip() for _, r in df.iterrows() if str(r["Cluster"]) != ""}
            col_name = "manual_celltype"
            ad.obs[col_name] = ad.obs[groupby].astype(str).map(mapping).fillna("")

            # UMAP 注释预览
            sc.pl.umap(ad, color=[groupby, col_name], show=False)
            fig_umap = plt.gcf()

            # 空间图（如有坐标）
            fig_spa = plt.figure(figsize=(4, 3))
            try:
                if hasattr(ad, "obsm") and ("spatial" in ad.obsm.keys()):
                    plt.close(fig_spa)
                    sc.pl.spatial(ad, color=col_name, show=False)
                    fig_spa = plt.gcf()
                else:
                    ax = fig_spa.add_subplot(111)
                    ax.text(0.5, 0.5, "无空间坐标，仅展示UMAP注释", ha="center", va="center")
                    ax.axis("off")
            except Exception:
                ax = fig_spa.add_subplot(111)
                ax.text(0.5, 0.5, "空间图绘制失败", ha="center", va="center")
                ax.axis("off")

            return (
                gr.update(value=fig_umap),
                gr.update(value=fig_spa),
                gr.update(value=df),
            )
        except Exception as e:
            fig = plt.figure(figsize=(6, 2))
            plt.text(0.5, 0.5, f"注释失败: {e}", ha="center", va="center")
            plt.axis("off")
            return gr.update(value=fig), gr.update(value=fig), gr.update(value=pd.DataFrame(columns=["Cluster", "名称"]))

    # 事件绑定

    gen_expr_btn.click(
        fn=generate_expression_plots,
        inputs=[marker_table],
        outputs=[expr_umap_fig, expr_dotplot_fig],
    )

    submit_btn.click(
        fn=submit_annotation,
        inputs=[cluster_map_table],
        outputs=[annot_umap_fig, annot_spatial_fig, mapping_preview],
    )

    def on_species_change(sp):
        return gr.update(value=preset_table(sp))

    species_dd.change(
        fn=on_species_change,
        inputs=[species_dd],
        outputs=[marker_table],
    )