#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Robust visualization script for endothelial -> CAF communication (inner vs outer)
Server-safe, no GUI dependency, outputs PNG + PDF
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path

# Force non-GUI backend for server environments
plt.switch_backend('Agg')


def clean_column_names(df):
    """Strip whitespace, lowercase all column names to avoid case/spacing errors"""
    df.columns = df.columns.str.strip().str.lower()
    return df


def check_required_columns(df, required_cols, file_path):
    """Validate required columns exist, exit gracefully if missing"""
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        print(f"\nERROR: Missing required columns in {file_path}: {missing}")
        print(f"INFO: Actual columns in file: {list(df.columns)}")
        print(f"HINT: Columns were auto-stripped and lowercased. Check spelling.")
        sys.exit(1)


def save_figure(fig_name, dpi=300):
    """Save both PNG (preview) and PDF (publication vector format)"""
    plt.savefig(f"{fig_name}.png", dpi=dpi, bbox_inches='tight')
    plt.savefig(f"{fig_name}.pdf", dpi=dpi, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Visualize endothelial-CAF communication differences")
    parser.add_argument("--inner", required=True, help="Path to inner mucinous zone CSV")
    parser.add_argument("--outer", required=True, help="Path to outer mucinous zone CSV")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    args = parser.parse_args()

    # Create output dir if not exists
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Load and clean data
    try:
        inner = pd.read_csv(args.inner,sep="\t")
        outer = pd.read_csv(args.outer,sep="\t")
    except FileNotFoundError as e:
        print(f"\nERROR: File not found: {e.filename}")
        sys.exit(1)

    inner = clean_column_names(inner)
    outer = clean_column_names(outer)

    # Validate required columns (match CLEANED column names, e.g. "ligand" not "Ligand")
    required_cols = ["ligand", "receptor", "target", "prob"]
    check_required_columns(inner, required_cols, args.inner)
    check_required_columns(outer, required_cols, args.outer)

    # Add region labels
    inner["region"] = "Inner"
    outer["region"] = "Outer"
    df = pd.concat([inner, outer], ignore_index=True)

    # Build plotting dataframe, aggregate duplicates to avoid pivot errors
    df_plot = df.copy()
    df_plot["interaction"] = df_plot["ligand"] + "-" + df_plot["receptor"]
    df_plot = df_plot.groupby(
        ["interaction", "target", "region"], as_index=False
    )["prob"].mean()

    # Calculate log2 fold change (Inner / Outer)
    prob_mat = df_plot.pivot_table(
        index=["interaction", "target"],
        columns="region",
        values="prob",
        fill_value=0
    ).reset_index()

    # Ensure both regions exist even if one has no data
    for col in ["Inner", "Outer"]:
        if col not in prob_mat.columns:
            prob_mat[col] = 0.0

    prob_mat["log2FC"] = np.log2(
        (prob_mat["Inner"] + 1e-5) / (prob_mat["Outer"] + 1e-5)
    )

    # Plot 1: Full dot plot (Inner vs Outer)
    plt.figure(figsize=(12, 10))
    plot_data = df_plot[df_plot["prob"] > 0].copy()
    if not plot_data.empty:
        sns.scatterplot(
            data=plot_data,
            x="target",
            y="interaction",
            size="prob",
            hue="region",
            sizes=(20, 400),
            alpha=0.85,
            palette={"Inner": "#D73027", "Outer": "#4575B4"}
        )
        plt.title("Endothelium -> CAF Communication: Inner vs Outer Mucinous Zone")
        plt.xlabel("Target CAF Subtype")
        plt.ylabel("Ligand-Receptor Interaction")
        plt.xticks(rotation=45, ha="right")
        plt.legend(title="Region", bbox_to_anchor=(1.05, 1), loc="upper left")
        save_figure(Path(args.outdir) / "endothelium_CAF_inner_vs_outer", dpi=300)

    # Plot 2: Differential heatmap (log2FC)
    plt.figure(figsize=(12, 10))
    heat_data = prob_mat.pivot(index="interaction", columns="target", values="log2FC")
    if not heat_data.empty:
        sns.heatmap(
            heat_data,
            cmap="RdBu_r",
            center=0,
            linewidths=0.5,
            linecolor="gray",
            cbar_kws={"label": "log2FC (Inner / Outer)"}
        )
        plt.title("Differential Endothelial Signaling to CAFs (Red = Inner Enriched)")
        plt.xticks(rotation=45, ha="right")
        save_figure(Path(args.outdir) / "endothelium_CAF_log2FC_heatmap", dpi=300)

    # Plot 3: Key pathway zoom (PDGF/NOTCH, suitable for main figure)
    highlight = [
        "PDGFB-PDGFRB",
        "DLL4-NOTCH2",
        "DLL4-NOTCH3",
        "JAG1-NOTCH2",
        "JAG1-NOTCH3"
    ]
    plt.figure(figsize=(10, 4))
    highlight_data = df_plot[df_plot["interaction"].isin(highlight)]
    if not highlight_data.empty:
        sns.scatterplot(
            data=highlight_data,
            x="target",
            y="interaction",
            size="prob",
            hue="region",
            sizes=(50, 500),
            alpha=0.9,
            palette={"Inner": "#D73027", "Outer": "#4575B4"}
        )
        plt.title("Key Tumor-Endothelial Signals (Inner Enriched)")
        plt.xlabel("Target CAF Subtype")
        plt.ylabel("Interaction")
        plt.xticks(rotation=45, ha="right")
        plt.legend(title="Region", bbox_to_anchor=(1.05, 1), loc="upper left")
        save_figure(Path(args.outdir) / "endothelium_CAF_key_interactions", dpi=300)

    print(f"\nSUCCESS: Run complete. Outputs saved to: {Path(args.outdir).absolute()}")
    print("Generated files:")
    print("  1. endothelium_CAF_inner_vs_outer.png/.pdf")
    print("  2. endothelium_CAF_log2FC_heatmap.png/.pdf")
    print("  3. endothelium_CAF_key_interactions.png/.pdf")


if __name__ == "__main__":
    main()
