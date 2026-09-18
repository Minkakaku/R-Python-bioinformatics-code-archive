import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import math
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

# 基础设定，输出信息设置为 3 显示详细 ;设置图片的基本参数
def setup_scanpy(verbosity: int = 3, dpi: int = 300) -> None:
    """Set basic scanpy/scvelo plotting options."""
    sc.settings.verbosity = verbosity
    sc.settings.set_figure_params(dpi=dpi, facecolor="white")

# 路径创建函数,保存路径的设定及创建
def make_dir(outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    return outdir

# 图片保存函数，保存pdf和png
def save_current_fig(outdir: str, file_name: str, dpi: int = 300, rasterize: bool = True) -> None:
    """Save current matplotlib figure as png and pdf."""
    make_dir(outdir)
    fig = plt.gcf()

    for ax in fig.axes:
        ax.grid(False) # 去除网格，如有需求可以修改
        if rasterize:
            for artist in ax.collections:
                try:
                    artist.set_rasterized(True)
                except Exception:
                    pass

    fig.savefig(os.path.join(outdir, f"{file_name}.png"), bbox_inches="tight", dpi=dpi)
    fig.savefig(os.path.join(outdir, f"{file_name}.pdf"), bbox_inches="tight", dpi=dpi)
    plt.close(fig)

# 函数，防止名称中有奇怪的地方
def safe_name(x: object) -> str:
    x = str(x)
    x = re.sub(r"[\\/:*?\"<>|\s]+", "_", x)
    x = x.replace("+", "pos")
    return x.strip("_")

# 查询文件夹中的loom文件
def find_loom_files(root_dir: str) -> List[Path]:
    """Recursively find all .loom files under root_dir."""
    root = Path(root_dir).expanduser().resolve() 
    loom_files = sorted(root.rglob("*.loom"))
    return loom_files

# 如果有需要重命名，需要传入重命名的文件，这里的文件是csv格式
def read_rename_file(rename_file: str) -> Dict[str, str]:
    """
    read the file for rename the sample names of loom data
    The first column is old sample name in loom obs_names, the second column is new sample name.[第一列是旧名称，第二列是新名称]
    """
    rename_df = pd.read_csv(rename_file, header=0)
    rename_df.iloc[:, 0] = rename_df.iloc[:, 0].astype(str).str.strip() 
    rename_df.iloc[:, 1] = rename_df.iloc[:, 1].astype(str).str.strip()

    dup_old = rename_df[rename_df.iloc[:, 0].duplicated(keep=False)] # 如果旧名有重复，会提示
    if dup_old.shape[0] > 0:
        print("Warning: old sample names are duplicated; later values will overwrite earlier values:") 
        print(dup_old)

    sample_map = dict(zip(rename_df.iloc[:, 0], rename_df.iloc[:, 1]))
    return sample_map 

# 根据上面csv生成的字典重新命名loom的barcode
def fix_loom_obsname_once(x: str, sample_map: Optional[Dict[str, str]] = None) -> str:
    """
    Convert loom obs_names from sample:barcodex to new_sample:barcode-1.
    Only map sample name once. Do not recursively map.
    """
    if ":" not in x: 
        sample, bc = "", x
    else:
        sample, bc = x.split(":", 1)

    sample_map = sample_map or {}
    new_sample = sample_map.get(sample, sample) # 从loom原始barcode获得名字，在dict里寻找新名字

    if bc.endswith("x"):
        bc = bc[:-1] + "-1"

    if new_sample == "":
        return bc
    return f"{new_sample}:{bc}"

# 执行上面的重命名,如果不传入重命名文件就会只处理末尾的barcode变化
def rename_loom_obs_names(adata_l, rename_file: Optional[str] = None):
    """Rename loom obs_names using optional rename_file."""
    sample_map = read_rename_file(rename_file) if rename_file else {}
    old_obs_names = adata_l.obs_names.copy()

    adata_l.obs_names = [fix_loom_obsname_once(x, sample_map) for x in old_obs_names]
    
    print("===== sample rename map =====")
    if sample_map:
        for old, new in sample_map.items():
            print(f"{old} -> {new}")
    else:
        print("No rename_file provided. Only barcode suffix x -> -1 was applied.")

    # 打印信息确认一下
    print("\n===== obs_names before/after first 10 =====")
    for old, new in zip(old_obs_names[:10], adata_l.obs_names[:10]):
        print(f"{old} -> {new}")

    print("\nunique:", adata_l.obs_names.is_unique)
    return adata_l

# 读入数据及整合函数
def read_and_concat_loom(loom_dir: str, outdir: str, loom_h5ad_name: str = "loom.h5ad"):
    """Read all loom files, concatenate them, and save/read a cached loom h5ad."""
    make_dir(outdir)
    loom_h5ad = os.path.join(outdir, loom_h5ad_name)

    files = find_loom_files(loom_dir)
    print("===== loom files =====")
    for f in files:
        print(f)
    print(f"Total .loom files found: {len(files)}")

    if len(files) == 0:
        raise FileNotFoundError(f"No .loom files found in {loom_dir}")

    adatas_list = []
    for f in files:
        print(f"Processing {f}...")
        data = scv.read_loom(str(f))
        data.var_names_make_unique()
        adatas_list.append(data)

    print(f"Total .loom files processed: {len(adatas_list)}")

    if len(adatas_list) == 1:
        adata_l = adatas_list[0]
    else:
        adata_l = adatas_list[0].concatenate(adatas_list[1:], index_unique=None)

    adata_l.write_h5ad(loom_h5ad) # 保存整合完的loom文件
    print(f"Saved loom h5ad: {loom_h5ad}")
    return adata_l

# 读入h5ad
def read_input_h5ad(input_h5ad: str):
    """Read h5ad. Optionally replace with raw.to_adata()."""
    print(f"[Read h5ad] {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    print(adata)
    print("obs head:")
    print(adata.obs.head())
    return adata

# 给h5ad文件进行修改barcode
def prefix_obs_names_by_sample(adata, sample_col: str = "sample", sep: str = ":"):
    """Change obs_names to sample:barcode for matching loom obs_names."""
    if sample_col not in adata.obs.columns:
        raise KeyError(f"sample_col '{sample_col}' not found in adata.obs")

    old = adata.obs_names.copy()
    adata.obs_names = [f"{s}{sep}{bc}" for s, bc in zip(adata.obs[sample_col].astype(str), old)]
    adata.obs_names_make_unique()

    print("===== adata obs_names before/after first 10 =====")
    for o, n in zip(old[:10], adata.obs_names[:10]):
        print(f"{o} -> {n}")
    print("\nunique:", adata.obs_names.is_unique)
    return adata


def merge_velocity_data(adata, adata_l, outdir: str, sample_col: str = "sample", groupby: str = "cell_type"):
    print("[Merge] scv.utils.merge(adata, adata_l)")
    adata_f = scv.utils.merge(adata, adata_l)
    print(adata_f)
    return adata_f

# 对读入的h5ad文件进行umap绘制
def plot_input_umap(adata, outdir: str, groupby: str = "cell_type", prefix: str = "input") -> None:
    if "X_umap" not in adata.obsm.keys():
        print("Warning: X_umap not found. Skip input UMAP plot.")
        return
    if groupby not in adata.obs.columns:
        print(f"Warning: {groupby} not found in adata.obs. Skip input UMAP plot.")
        return
    sc.pl.umap(adata, color=[groupby], show=False, legend_loc="on data",legend_fontsize = 8)
    save_current_fig(outdir, f"{prefix}_umap_{groupby}")

# 进行确认，loom修改后的名称的sample是否包含有全部的读入的h5ad的sample，如果没有就中断，然后根据打印信息确认

def check_loom_sample_contains_adata_sample(
    adata,
    adata_l,
    sample_col: str = "sample",
    sep: str = ":"
):
    if sample_col not in adata.obs.columns:
        raise ValueError(f"sample_col '{sample_col}' not found in adata.obs")
    
    # adata 里的 sample
    adata_samples = set(
        adata.obs[sample_col]
        .astype(str)
        .str.strip()
        .dropna()
        .unique()
    )

    # loom obs_names 里 ":" 前面的 sample
    loom_obs_names = adata_l.obs_names.astype(str)

    loom_samples = set(
        x.split(sep, 1)[0].strip()
        for x in loom_obs_names
    )

    missing_samples = sorted(adata_samples - loom_samples)
    extra_loom_samples = sorted(loom_samples - adata_samples)

    print("===== sample check before merge =====")
    print("\nSamples in adata:")
    print(sorted(adata_samples))

    print("\nSamples in loom:")
    print(sorted(loom_samples))

    if len(missing_samples) > 0:
        print("\nERROR: samples in adata but not found in loom:")
        print(missing_samples)

        raise ValueError(
            "Sample check failed before merge. "
            "Some samples in adata are not found in loom obs_names prefix. "
            "Please check rename_file or sample_col."
        )

    print("\nSample check passed before merge.")
    return True



# 预处理函数 ，都是常规默认参数 https://scvelo.readthedocs.io/en/stable/VelocityBasics.html
def run_scvelo_preprocess(
    adata,
    min_shared_counts: int = 20,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 30,
):
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=min_shared_counts,
        n_top_genes=n_top_genes,
    )
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors) 
    scv.pp.neighbors(adata)
    return adata

# 执行分析函数,这里mode默认是"stochastic"，还有'deterministic',或 'dynamical' 可选
def run_velocity(adata, mode: str = "stochastic"): 
    if mode == "dynamical": 
        scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata)
    return adata


# 绘图函数
def plot_velocity_results( 
    adata,
    outdir: str,
    groupby: str = "cell_type",
    basis: str = "umap",
    prefix: str = "scvelo",
    rank_velocity_genes: bool = True,
    top_gene_n: int = 5,
) -> None:
    make_dir(outdir)

    # 1. stream plot
    scv.pl.velocity_embedding_stream(
        adata,
        basis=basis,
        color=groupby,
        legend_loc="on data",
        legend_fontsize=5,
        dpi=300,
        show=False,
    )
    save_current_fig(outdir, f"{prefix}_stream")

    # 2. arrow/grid plot
    scv.pl.velocity_embedding_grid(
        adata,
        basis=basis,
        color=groupby,
        title="",
        scale=0.25,
        show=False,
    )
    save_current_fig(outdir, f"{prefix}_arrow")

    # 3. pseudotime
    try:
        scv.tl.velocity_pseudotime(adata)
        scv.pl.scatter(adata, basis=basis, color="velocity_pseudotime", cmap="gnuplot", show=False)
        save_current_fig(outdir, f"{prefix}_pseudotime")
    except Exception as e:
        print(f"Warning: velocity_pseudotime failed: {e}")

    # 4. spliced/unspliced proportions
    try:
        scv.pl.proportions(adata, groupby=groupby, show=False)
        save_current_fig(outdir, f"spliced_unspliced_proportions_{groupby}")
    except Exception as e:
        print(f"Warning: proportions plot failed: {e}")

    # 5. confidence
    try:
        scv.tl.velocity_confidence(adata)
        keys = ["velocity_length", "velocity_confidence"]
        df = adata.obs.groupby(groupby, observed=False)[keys].mean().T
        df.to_csv(os.path.join(outdir, "velocity_confidence_by_group.csv"), encoding="utf-8-sig")
        scv.pl.scatter(adata, c=keys, cmap="coolwarm", show=False, perc=[5, 95])
        save_current_fig(outdir, f"{prefix}_velocity_confidence")
    except Exception as e:
        print(f"Warning: velocity_confidence failed: {e}")

    # 6. PAGA
    try:
        try:
            scv.tl.paga(adata, groups=groupby)
        except TypeError:
            scv.tl.paga(adata, groupby=groupby)
        scv.pl.paga(
            adata,
            basis=basis,
            size=50, # 大小
            alpha=0.1, # 颜色透明度
            min_edge_width=2,
            node_size_scale=1.5,
            show=False,
        )
        fig = plt.gcf()
        for ax in fig.axes:
            ax.grid(False)
            for txt in ax.texts:
                txt.set_visible(False)
        fig.savefig(os.path.join(outdir, f"{prefix}_PAGA.png"), bbox_inches="tight", dpi=300)
        fig.savefig(os.path.join(outdir, f"{prefix}_PAGA.pdf"), bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"Warning: PAGA failed: {e}")

    # 7. velocity genes
    if rank_velocity_genes:
        try:
            scv.tl.rank_velocity_genes(adata, groupby=groupby, min_corr=0.3)
            df_rank = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])
            df_rank.to_csv(os.path.join(outdir, f"rank_velocity_genes_by_{groupby}.csv"), encoding="utf-8-sig")

            gene_dir = os.path.join(outdir, "velocity_gene_plots")
            make_dir(gene_dir)

            for ct in df_rank.columns:
                genes = df_rank[ct].dropna().head(top_gene_n).tolist()
                genes = [g for g in genes if g in adata.var_names]
                if len(genes) == 0:
                    print(f"{ct}: no valid genes, skip")
                    continue
                cols = math.ceil(math.sqrt(top_gene_n))
                wid = 9*cols + 1
                hei = math.ceil(top_gene_n / cols) * 4.5
                print(ct, genes)
                scv.pl.velocity(
                    adata,
                    var_names=genes,
                    color=groupby,
                    basis=basis,
                    ncols= cols,
                    show=False,
                )
                fig = plt.gcf()
                fig.set_size_inches(wid, hei) # 根据cols，计算图片尺寸
                fig.subplots_adjust(left=0.05, right=0.98, bottom=0.06, top=0.95, wspace=0.35, hspace=0.35)
                for ax in fig.axes:
                    ax.grid(False)
                fig.savefig(os.path.join(gene_dir, f"{safe_name(ct)}_top{top_gene_n}_velocity.png"), bbox_inches="tight", dpi=300)
                fig.savefig(os.path.join(gene_dir, f"{safe_name(ct)}_top{top_gene_n}_velocity.pdf"), bbox_inches="tight", dpi=300)
                plt.close(fig)
        except Exception as e:
            print(f"Warning: rank_velocity_genes failed: {e}")


def run_velocity_pipeline(
    input_h5ad: str,
    loom_dir: str,
    output_dir: str,
    rename_file: Optional[str] = None,
    sample_col: str = "sample",
    groupby: str = "cell_type",
    add_sample_prefix: bool = True,
    mode: str = "stochastic", # 默认参数，未设置传参
    basis: str = "umap",
    n_top_genes: int = 2000, # 默认参数，未设置传参
    n_pcs: int = 30,# 默认参数，未设置传参
    n_neighbors: int = 30,# 默认参数，未设置传参
    min_shared_counts: int = 20,# 默认参数，未设置传参
    prefix: str = "scvelo", 
    rank_velocity_genes: bool = True,
    top_gene_n: int = 5, # 默认绘制top5基因的
):
    make_dir(output_dir) # 创建输出目录
    setup_scanpy()  # 基础设定，输出信息设置为 3 显示详细 ;设置图片的基本参数

    adata = read_input_h5ad(input_h5ad) # 读入h5ad
    plot_input_umap(adata, output_dir, groupby=groupby, prefix="input") # umap检查

    adata_l = read_and_concat_loom(loom_dir, output_dir) # 读入loom文件夹内的文件
    adata_l = rename_loom_obs_names(adata_l, rename_file=rename_file) # loom barcode的重命名

    if add_sample_prefix:
        adata = prefix_obs_names_by_sample(adata, sample_col=sample_col, sep=":") # 修改h5ad的barcode命名
    
    check_loom_sample_contains_adata_sample(
        adata=adata,
        adata_l=adata_l,
        sample_col=sample_col,
        sep=":"
    )
    adata_f = merge_velocity_data(adata, adata_l, outdir=output_dir, sample_col=sample_col, groupby=groupby) # 合并loom和h5ad

    adata_f = run_scvelo_preprocess(
        adata_f,
        min_shared_counts=min_shared_counts,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
    ) # 预处理
    adata_f = run_velocity(adata_f, mode=mode) # 分析
    h5ad_out = os.path.join(output_dir, f"{prefix}_{mode}.h5ad")
    adata_f.write_h5ad(h5ad_out) # 保存
    print(f"Saved velocity result: {h5ad_out}")
    plot_velocity_results(
        adata_f,
        outdir=output_dir,
        groupby=groupby,
        basis=basis,
        prefix=f"{prefix}_{mode}",
        rank_velocity_genes=rank_velocity_genes,
        top_gene_n=top_gene_n,
    ) # 绘图

    return adata_f
