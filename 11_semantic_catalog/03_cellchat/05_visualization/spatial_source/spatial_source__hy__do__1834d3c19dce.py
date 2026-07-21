#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare CellChat results: Inner vs Outer (Target: EMT-Epi)
Includes: Bubble Plot + Heatmap (log2FC) + Key Pathway Plot
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path

plt.switch_backend('Agg')

def clean_columns(df):
    df.columns = df.columns.str.strip().str.lower()
    return df

def check_cols(df, req, path):
    missing = [c for c in req if c not in df.columns]
    if missing:
        print(f"ERROR: Missing columns in {path}: {missing}")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)

def save_fig(name, dpi=300):
    plt.savefig(f"{name}.png", dpi=dpi, bbox_inches='tight')
    plt.savefig(f"{name}.pdf", dpi=dpi, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inner", required=True)
    parser.add_argument("--outer", required=True)
    parser.add_argument("--outdir", default="./plots")
    args = parser.parse_args()

    Path(args.outdir).mkdir(exist_ok=True)

    # Read TSV
    try:
        inner = pd.read_csv(args.inner, sep='\t')
        outer = pd.read_csv(args.outer, sep='\t')
    except Exception as e:
        print(f"Failed to read files: {e}")
        sys.exit(1)

    inner = clean_columns(inner)
    outer = clean_columns(outer)

    req = ["source", "target", "ligand", "receptor", "prob"]
    check_cols(inner, req, args.inner)
    check_cols(outer, req, args.outer)

    # Add region
    inner["region"] = "Inner"
    outer["region"] = "Outer"
    df = pd.concat([inner, outer], ignore_index=True)
    
    # Optional: Filter for EMT-Epi targets (uncomment if needed)
    # df = df[df["target"].str.contains("EMT|Epi", case=False, na=False)]
    
    df["interaction"] = df["ligand"] + "-" + df["receptor"]
    
    # Aggregate
    df_agg = df.groupby(["source", "target", "interaction", "region"], as_index=False)["prob"].mean()
    plot_data = df_agg[df_agg["prob"] > 0]

    # ==================== Plot 1: Bubble Plot ====================
    plt.figure(figsize=(14, 10))
    sns.scatterplot(
        data=plot_data,
        x="source",
        y="interaction",
        size="prob",
        hue="region",
        sizes=(20, 500),
        alpha=0.8,
        palette={"Inner": "#D73027", "Outer": "#4575B4"}
    )
    plt.title("Signaling to EMT-Epi: Inner vs Outer", fontsize=16)
    plt.xlabel("Source Cell Type", fontsize=12)
    plt.ylabel("Ligand-Receptor Interaction", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Region", bbox_to_anchor=(1.05, 1), loc='upper left')
    save_fig(Path(args.outdir) / "EMT_Epi_BubblePlot")

    # ==================== Plot 2: Heatmap (log2FC) ====================
    # Pivot table for probability
    prob_pivot = df_agg.pivot_table(
        index=["interaction", "source"],  # Group by interaction and source
        columns="region",
        values="prob",
        fill_value=0
    ).reset_index()

    # Ensure both regions exist
    for reg in ["Inner", "Outer"]:
        if reg not in prob_pivot.columns:
            prob_pivot[reg] = 0.0

    # Calculate log2FC (Inner / Outer)
    prob_pivot["log2FC"] = np.log2(
        (prob_pivot["Inner"] + 1e-5) / (prob_pivot["Outer"] + 1e-5)
    )

    # Reshape for heatmap: rows=interaction, cols=source, values=log2FC
    heatmap_data = prob_pivot.pivot_table(
        index="interaction",
        columns="source",
        values="log2FC",
        fill_value=0
    )

    plt.figure(figsize=(14, 10))
    sns.heatmap(
        heatmap_data,
        cmap="RdBu_r",  # Red-Blue diverging: Red=Inner enriched, Blue=Outer enriched
        center=0,
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "log2FC (Inner / Outer)"}
    )
    plt.title("Differential Signaling to EMT-Epi (Inner vs Outer)", fontsize=16)
    plt.xlabel("Source Cell Type", fontsize=12)
    plt.ylabel("Ligand-Receptor Interaction", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    save_fig(Path(args.outdir) / "EMT_Epi_DiffHeatmap")

    # ==================== Plot 3: Key Pathways ====================
    highlight = ["PDGFB-PDGFRB", "DLL4-NOTCH3", "JAG1-NOTCH3", "COL4A1-CD44", "THBS1-CD47"]
    highlight_data = plot_data[plot_data["interaction"].isin(highlight)]
    
    if not highlight_data.empty:
        plt.figure(figsize=(12, 5))
        sns.scatterplot(
            data=highlight_data,
            x="source",
            y="interaction",
            size="prob",
            hue="region",
            sizes=(50, 600),
            alpha=0.9,
            palette={"Inner": "#D73027", "Outer": "#4575B4"}
        )
        plt.title("Key Signaling Pathways to EMT-Epi", fontsize=14)
        plt.xlabel("Source Cell Type", fontsize=12)
        plt.ylabel("Interaction", fontsize=12)
        plt.xticks(rotation=45, ha="right")
        plt.legend(title="Region", bbox_to_anchor=(1.05, 1), loc='upper left')
        save_fig(Path(args.outdir) / "EMT_Epi_KeySignaling")

    print(f"\n✅ Done! Generated 3 figures in {args.outdir}:")
    print("1. EMT_Epi_BubblePlot.png/pdf")
    print("2. EMT_Epi_DiffHeatmap.png/pdf  <-- Your Heatmap")
    print("3. EMT_Epi_KeySignaling.png/pdf")

if __name__ == "__main__":
    main()
