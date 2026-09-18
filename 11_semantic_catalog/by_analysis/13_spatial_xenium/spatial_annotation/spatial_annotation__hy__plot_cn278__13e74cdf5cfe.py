#!/usr/bin/env python3
"""
细胞类型富集分析火山图 - CNS主刊发表级别
分析目标克隆(CN2, CN7, CN8)与其他细胞在细胞类型组成上的差异
增加功能：字体放大、所有PNG同时输出PDF
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import warnings
import os
warnings.filterwarnings('ignore')
import matplotlib.font_manager as font_manager   # ★ 必须导入

# 字体路径
font_path = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'

# 添加字体并全局启用
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'Arial'          # 首字母大写！
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42



# ============================================================================
# 0. 全局设置 - CNS级别样式 (字体进一步放大)
# ============================================================================
#matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
#matplotlib.rcParams['mathtext.fontset'] = 'stix'

# 高质量图形设置 + 字体放大
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['font.size'] = 16          # 从14 → 16
plt.rcParams['axes.titlesize'] = 20     # 从18 → 20
plt.rcParams['axes.labelsize'] = 18     # 从16 → 18
plt.rcParams['xtick.labelsize'] = 16    # 从14 → 16
plt.rcParams['ytick.labelsize'] = 16    # 从14 → 16
plt.rcParams['legend.fontsize'] = 14    # 从12 → 14
plt.rcParams['figure.titlesize'] = 24   # 从20 → 24
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['figure.max_open_warning'] = 50

output_dir = 'volcano'
os.makedirs(output_dir, exist_ok=True)

# ============================================================================
# 1. 数据加载和预处理
# ============================================================================
print("=" * 80)
print("步骤1: 加载数据并去除LYVE1+_TAM")
print("=" * 80)

adata = sc.read_h5ad('adata_with_CNs10_with_patient.h5ad')
print(f"原始数据维度: {adata.shape}")

if 'annotation' in adata.obs.columns:
    n_before = adata.shape[0]
    adata = adata[adata.obs['annotation'] != 'LYVE1+_TAM'].copy()
    print(f"去除LYVE1+_TAM: {n_before - adata.shape[0]} 个细胞")
    print(f"过滤后数据维度: {adata.shape}")

cn_col = None
for col in adata.obs.columns:
    if 'CN' in col.upper():
        sample_vals = adata.obs[col].dropna().astype(str).head(20)
        if any('CN' in str(v) for v in sample_vals):
            cn_col = col
            break

annot_col = 'annotation'
if annot_col not in adata.obs.columns:
    for col in adata.obs.columns:
        if col.lower() in ['annotation', 'cell_type', 'celltype', 'label']:
            annot_col = col
            break

print(f"CN列: '{cn_col}', 标注列: '{annot_col}'")

# ============================================================================
# 2. 分析函数
# ============================================================================
def run_fisher_test(adata, cn_col, annot_col, target_cns, reference_label):
    """Fisher精确检验"""
    adata_target = adata[adata.obs[cn_col].astype(str).isin(target_cns)].copy()
    adata_other = adata[~adata.obs[cn_col].astype(str).isin(target_cns)].copy()
    
    print(f"  目标 ({','.join(target_cns)}): {adata_target.shape[0]} 个细胞")
    print(f"  参考 ({reference_label}): {adata_other.shape[0]} 个细胞")
    
    all_cell_types = sorted(list(set(
        list(adata_target.obs[annot_col].dropna().unique()) + 
        list(adata_other.obs[annot_col].dropna().unique())
    )))
    
    results = []
    for cell_type in all_cell_types:
        target_count = (adata_target.obs[annot_col] == cell_type).sum()
        other_count = (adata_other.obs[annot_col] == cell_type).sum()
        
        table = [[target_count, len(adata_target) - target_count],
                 [other_count, len(adata_other) - other_count]]
        
        odds_ratio, p_val = fisher_exact(table)
        
        target_prop = target_count / len(adata_target) if len(adata_target) > 0 else 0
        other_prop = other_count / len(adata_other) if len(adata_other) > 0 else 0
        fc = target_prop / other_prop if other_prop > 0 else np.inf
        
        results.append({
            'Cell_Type': cell_type,
            'Target_Count': target_count,
            'Other_Count': other_count,
            'Target_Prop': target_prop,
            'Other_Prop': other_prop,
            'Fold_Enrichment': fc,
            'Odds_Ratio': odds_ratio,
            'Log2_Odds_Ratio': np.log2(odds_ratio) if odds_ratio > 0 else -np.inf,
            'P_Value': p_val
        })
    
    df = pd.DataFrame(results)
    df['FDR'] = multipletests(df['P_Value'], method='fdr_bh')[1]
    
    df = df.replace([np.inf, -np.inf], np.nan)
    max_fc = df['Log2_Odds_Ratio'].abs().max()
    if pd.isna(max_fc) or max_fc == 0:
        max_fc = 5
    df['Log2_Odds_Ratio'] = df['Log2_Odds_Ratio'].fillna(max_fc * 1.1)
    
    df['Significance'] = 'Not Significant'
    df.loc[(df['FDR'] < 0.05) & (df['Log2_Odds_Ratio'] > 0), 'Significance'] = 'Enriched (FDR<0.05)'
    df.loc[(df['FDR'] < 0.05) & (df['Log2_Odds_Ratio'] < 0), 'Significance'] = 'Depleted (FDR<0.05)'
    df.loc[(df['FDR'] < 0.01) & (df['Log2_Odds_Ratio'] > 0), 'Significance'] = 'Enriched (FDR<0.01)'
    df.loc[(df['FDR'] < 0.01) & (df['Log2_Odds_Ratio'] < 0), 'Significance'] = 'Depleted (FDR<0.01)'
    
    return df


def plot_volcano(fisher_df, target_cns, reference_label, output_file, only_positive=False):
    """绘制CNS级别火山图，同时保存PNG和PDF"""
    fig = plt.figure(figsize=(26, 18))  # 略微加大图形
    
    colors = {
        'Not Significant': '#BDBDBD',
        'Depleted (FDR<0.05)': '#67A9CF', 
        'Enriched (FDR<0.05)': '#EF8A62',
        'Depleted (FDR<0.01)': '#2166AC',
        'Enriched (FDR<0.01)': '#B2182B'
    }
    
    plot_df = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] > 0)].copy() if only_positive else fisher_df.copy()
    title_suffix = '\n(Only Significantly Enriched, FDR < 0.05)' if only_positive else ''
    
    # ===== A: 火山图 =====
    ax1 = plt.subplot(2, 2, (1, 2))
    
    for sig in ['Not Significant', 'Depleted (FDR<0.05)', 'Enriched (FDR<0.05)',
                'Depleted (FDR<0.01)', 'Enriched (FDR<0.01)']:
        mask = plot_df['Significance'] == sig
        if mask.sum() > 0:
            subset = plot_df[mask]
            sizes = 200 if '0.01' in sig else 120
            edge_width = 1.5 if '0.01' in sig else 0.5
            edge_color = '#333333' if '0.01' in sig else 'none'
            ax1.scatter(subset['Log2_Odds_Ratio'], -np.log10(subset['P_Value']),
                       c=colors[sig], s=sizes, alpha=0.85, 
                       edgecolors=edge_color, linewidth=edge_width,
                       label=f'{sig} (n={mask.sum()})',
                       zorder=3 if '0.01' in sig else 2)
    
    # 阈值线
    ax1.axhline(-np.log10(0.05), color='#333333', linestyle='--', alpha=0.6, linewidth=2.0)
    ax1.axhline(-np.log10(0.01), color='#999999', linestyle='--', alpha=0.4, linewidth=1.5)
    ax1.axvline(0, color='#333333', linestyle='-', alpha=0.3, linewidth=1.5)
    
    # 阈值文字 (加大字体)
    ax1.text(0.98, -np.log10(0.05) + 0.1, 'FDR = 0.05', transform=ax1.get_yaxis_transform(),
            fontsize=12, color='#333333', ha='right', va='bottom', style='italic')
    ax1.text(0.98, -np.log10(0.01) + 0.1, 'FDR = 0.01', transform=ax1.get_yaxis_transform(),
            fontsize=12, color='#999999', ha='right', va='bottom', style='italic')
    
    ax1.set_xlabel('log₂(Odds Ratio)', fontweight='bold', fontsize=20, labelpad=10)
    ax1.set_ylabel('−log₁₀(P-value)', fontweight='bold', fontsize=20, labelpad=10)
    
    target_name = ' + '.join(target_cns)
    ax1.set_title(f'Cell Type Enrichment Analysis\n{target_name} vs. {reference_label}{title_suffix}', 
                  fontweight='bold', fontsize=22, pad=15)
    
    legend = ax1.legend(fontsize=14, loc='upper left', frameon=True, 
                       fancybox=True, shadow=False, framealpha=0.95,
                       edgecolor='#CCCCCC', ncol=1,
                       title='Significance Level', title_fontsize=15)
    legend.get_title().set_fontweight('bold')
    
    ax1.grid(True, alpha=0.15, linestyle='-', linewidth=0.5)
    ax1.set_axisbelow(True)
    
    for spine in ax1.spines.values():
        spine.set_linewidth(1.5)
    
    if len(plot_df) > 0:
        max_x = max(abs(plot_df['Log2_Odds_Ratio'].min()), abs(plot_df['Log2_Odds_Ratio'].max())) * 1.2
        ax1.set_xlim(-max_x, max_x)
    
    # 标注显著细胞类型
    sig_data = plot_df[plot_df['FDR'] < 0.05].copy()
    if len(sig_data) > 0:
        sig_data = sig_data.sort_values('P_Value')
        annotations = []
        for idx, row in sig_data.iterrows():
            x, y = row['Log2_Odds_Ratio'], -np.log10(row['P_Value'])
            too_close = False
            for ann_x, ann_y in annotations:
                if abs(x - ann_x) < max_x * 0.08 and abs(y - ann_y) < 1.5:
                    too_close = True
                    break
            if not too_close or row['FDR'] < 0.001:
                annotations.append((x, y))
                fontweight = 'bold' if row['FDR'] < 0.01 else 'normal'
                fontsize = 10 if row['FDR'] < 0.01 else 9
                color = '#B2182B' if row['FDR'] < 0.01 else '#333333'
                ax1.annotate(row['Cell_Type'], (x, y),
                            xytext=(10, 10), textcoords='offset points', 
                            fontsize=fontsize, fontweight=fontweight, color=color,
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                     edgecolor='#CCCCCC', alpha=0.9),
                            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', 
                                          color='#666666', alpha=0.6, lw=0.8))
    
    # ===== B: 显著富集细胞类型 =====
    ax2 = plt.subplot(2, 2, 3)
    
    if only_positive:
        enriched_df = plot_df.nlargest(15, 'Log2_Odds_Ratio')
    else:
        enriched_df = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] > 0)].nlargest(15, 'Log2_Odds_Ratio')
    
    if len(enriched_df) > 0:
        colors_bar = plt.cm.Reds(np.linspace(0.4, 0.9, len(enriched_df)))[::-1]
        y_pos = range(len(enriched_df))
        ax2.barh(y_pos, enriched_df['Log2_Odds_Ratio'], 
                color=colors_bar, edgecolor='#333333', linewidth=0.8, height=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(enriched_df['Cell_Type'], fontsize=14)
        for i, (_, row) in enumerate(enriched_df.iterrows()):
            ax2.text(row['Log2_Odds_Ratio'] + max(enriched_df['Log2_Odds_Ratio']) * 0.02, i, 
                    f'{row["Log2_Odds_Ratio"]:.2f}', 
                    va='center', fontsize=13, fontweight='bold', color='#333333')
    else:
        ax2.text(0.5, 0.5, 'No significantly enriched cell types', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=16, color='#999999')
    
    ax2.set_xlabel('log₂(Odds Ratio)', fontweight='bold', fontsize=18, labelpad=10)
    ax2.set_title(f'Significantly Enriched in Target', fontweight='bold', fontsize=20, pad=12)
    ax2.axvline(0, color='#333333', linestyle='-', alpha=0.5, linewidth=1.5)
    ax2.grid(True, alpha=0.15, axis='x', linestyle='-', linewidth=0.5)
    ax2.set_axisbelow(True)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.5)
    
    # ===== C: 显著缺失/景观图 =====
    ax3 = plt.subplot(2, 2, 4)
    
    if not only_positive:
        depleted_df = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] < 0)].nsmallest(15, 'Log2_Odds_Ratio')
        if len(depleted_df) > 0:
            colors_bar = plt.cm.Blues(np.linspace(0.4, 0.9, len(depleted_df)))
            y_pos = range(len(depleted_df))
            ax3.barh(y_pos, depleted_df['Log2_Odds_Ratio'], 
                    color=colors_bar, edgecolor='#333333', linewidth=0.8, height=0.7)
            ax3.set_yticks(y_pos)
            ax3.set_yticklabels(depleted_df['Cell_Type'], fontsize=14)
            for i, (_, row) in enumerate(depleted_df.iterrows()):
                ax3.text(row['Log2_Odds_Ratio'] - max(abs(depleted_df['Log2_Odds_Ratio'])) * 0.02, i, 
                        f'{row["Log2_Odds_Ratio"]:.2f}', 
                        va='center', ha='right', fontsize=13, fontweight='bold', color='#333333')
        else:
            ax3.text(0.5, 0.5, 'No significantly depleted cell types', 
                    ha='center', va='center', transform=ax3.transAxes, fontsize=16, color='#999999')
        ax3.set_title(f'Significantly Depleted in Target', fontweight='bold', fontsize=20, pad=12)
    else:
        if len(plot_df) > 0:
            sorted_df = plot_df.sort_values('Log2_Odds_Ratio', ascending=False)
            colors_scatter = plt.cm.Reds(np.linspace(0.4, 0.9, len(sorted_df)))
            ax3.scatter(range(len(sorted_df)), sorted_df['Log2_Odds_Ratio'],
                       c=colors_scatter, s=120, edgecolors='#333333', linewidth=0.8, zorder=3)
            ax3.axhline(0, color='#333333', linestyle='-', alpha=0.3, linewidth=1.5)
            ax3.set_xlabel('Rank', fontweight='bold', fontsize=18, labelpad=10)
            ax3.set_ylabel('log₂(Odds Ratio)', fontweight='bold', fontsize=18, labelpad=10)
        ax3.set_title(f'Enrichment Landscape (n = {len(plot_df)})', fontweight='bold', fontsize=20, pad=12)
    
    ax3.set_xlabel('log₂(Odds Ratio)', fontweight='bold', fontsize=18, labelpad=10)
    ax3.axvline(0, color='#333333', linestyle='-', alpha=0.5, linewidth=1.5)
    ax3.grid(True, alpha=0.15, axis='x', linestyle='-', linewidth=0.5)
    ax3.set_axisbelow(True)
    for spine in ax3.spines.values():
        spine.set_linewidth(1.5)
    
    # 总标题
    fig.suptitle(f'Cell Type Composition Analysis\nTarget Clones: {target_name} vs. {reference_label}', 
                 fontweight='bold', fontsize=24, y=0.99)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存 PNG 和 PDF
    png_path = output_file
    pdf_path = output_file.replace('.png', '.pdf')
    plt.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"  ✓ 已保存: {png_path}")
    print(f"  ✓ 已保存: {pdf_path}")


def save_results(fisher_df, output_prefix):
    """保存Excel结果"""
    excel_file = f'{output_prefix}_results.xlsx'
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        fisher_df.to_excel(writer, sheet_name='Full_Results', index=False)
        enriched = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] > 0)]
        if len(enriched) > 0:
            enriched.to_excel(writer, sheet_name='Enriched', index=False)
        depleted = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] < 0)]
        if len(depleted) > 0:
            depleted.to_excel(writer, sheet_name='Depleted', index=False)
    print(f"  ✓ Excel: {excel_file}")


def print_summary(fisher_df, target_name):
    """打印摘要"""
    sig_enriched = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] > 0)]
    sig_depleted = fisher_df[(fisher_df['FDR'] < 0.05) & (fisher_df['Log2_Odds_Ratio'] < 0)]
    print(f"\n  {'='*60}")
    print(f"  {target_name}")
    print(f"  {'='*60}")
    print(f"  Significantly enriched (FDR<0.05): {len(sig_enriched)} cell types")
    print(f"  Significantly depleted (FDR<0.05): {len(sig_depleted)} cell types")
    if len(sig_enriched) > 0:
        print(f"\n  Top 10 Enriched:")
        for i, (_, row) in enumerate(sig_enriched.nlargest(10, 'Log2_Odds_Ratio').iterrows()):
            stars = '***' if row['FDR'] < 0.001 else ('**' if row['FDR'] < 0.01 else '*')
            print(f"    {i+1:2d}. {row['Cell_Type']:45s} log2FC = {row['Log2_Odds_Ratio']:6.2f}  FDR = {row['FDR']:.2e} {stars}")
    if len(sig_depleted) > 0:
        print(f"\n  Top 10 Depleted:")
        for i, (_, row) in enumerate(sig_depleted.nsmallest(10, 'Log2_Odds_Ratio').iterrows()):
            stars = '***' if row['FDR'] < 0.001 else ('**' if row['FDR'] < 0.01 else '*')
            print(f"    {i+1:2d}. {row['Cell_Type']:45s} log2FC = {row['Log2_Odds_Ratio']:6.2f}  FDR = {row['FDR']:.2e} {stars}")


# ============================================================================
# 3. 分析每个CN克隆
# ============================================================================
print("\n" + "=" * 80)
print("步骤2: 分别分析CN2, CN7, CN8 vs 其他克隆")
print("=" * 80)

all_target_cns = ['CN2', 'CN7', 'CN8']
reference_label = "Other CNs (excluding CN2,7,8)"
all_fisher_results = {}

for target_cn in all_target_cns:
    print(f"\n{'='*60}")
    print(f"Analyzing: {target_cn} vs {reference_label}")
    print(f"{'='*60}")
    fisher_df = run_fisher_test(adata, cn_col, annot_col, [target_cn], reference_label)
    all_fisher_results[target_cn] = fisher_df
    output_prefix = os.path.join(output_dir, f'CN278_{target_cn}')
    save_results(fisher_df, output_prefix)
    plot_volcano(fisher_df, [target_cn], reference_label, f'{output_prefix}_volcano.png', only_positive=False)
    print_summary(fisher_df, target_cn)


# ============================================================================
# 4. 综合分析 CN2+CN7+CN8 vs Others
# ============================================================================
print("\n" + "=" * 80)
print("步骤3: 综合分析 CN2+CN7+CN8 vs Others")
print("=" * 80)

fisher_combined = run_fisher_test(adata, cn_col, annot_col, all_target_cns, reference_label)
output_prefix_combined = os.path.join(output_dir, 'CN278_combined')
save_results(fisher_combined, output_prefix_combined)
plot_volcano(fisher_combined, all_target_cns, reference_label, 
             f'{output_prefix_combined}_enriched_only.png', only_positive=True)
plot_volcano(fisher_combined, all_target_cns, reference_label, 
             f'{output_prefix_combined}_full.png', only_positive=False)
print_summary(fisher_combined, 'CN2 + CN7 + CN8 (Combined)')


# ============================================================================
# 5. 一致性分析
# ============================================================================
print("\n" + "=" * 80)
print("步骤4: 一致性分析")
print("=" * 80)

enriched_sets = {}
for cn in all_target_cns:
    df = all_fisher_results[cn]
    enriched_sets[cn] = set(df[(df['FDR'] < 0.05) & (df['Log2_Odds_Ratio'] > 0)]['Cell_Type'].tolist())

for cn, types in enriched_sets.items():
    print(f"  {cn}: {len(types)} significantly enriched cell types")

common_enriched = set.intersection(*enriched_sets.values())
print(f"\n  Common enriched across all 3 clones: {len(common_enriched)}")
if common_enriched:
    for ct in sorted(common_enriched):
        print(f"    • {ct}")

for cn1, cn2 in combinations(all_target_cns, 2):
    common = enriched_sets[cn1] & enriched_sets[cn2]
    print(f"  {cn1} ∩ {cn2}: {len(common)}")

for cn in all_target_cns:
    others = set(all_target_cns) - {cn}
    unique = enriched_sets[cn] - set.union(*[enriched_sets[o] for o in others])
    print(f"  {cn} unique: {len(unique)}")

# 一致性热图
all_enriched_types = sorted(list(set.union(*enriched_sets.values())))
if len(all_enriched_types) > 0:
    consistency_matrix = pd.DataFrame(0.0, index=all_enriched_types, columns=all_target_cns)
    for cn in all_target_cns:
        df = all_fisher_results[cn]
        for ct in all_enriched_types:
            if ct in enriched_sets[cn]:
                row = df[df['Cell_Type'] == ct]
                if len(row) > 0:
                    consistency_matrix.loc[ct, cn] = row['Log2_Odds_Ratio'].values[0]
    
    fig, ax = plt.subplots(figsize=(12, max(10, len(all_enriched_types) * 0.45)))
    sns.heatmap(consistency_matrix, annot=True, fmt='.2f', cmap='YlOrRd', 
                cbar_kws={'label': 'log₂(Odds Ratio)', 'shrink': 0.8},
                ax=ax, linewidths=0.8, linecolor='white',
                vmin=0, annot_kws={'fontsize': 14, 'fontweight': 'bold'},
                cbar_ax=None)
    
    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel('log₂(Odds Ratio)', fontsize=16, fontweight='bold', labelpad=15)
    cbar.ax.tick_params(labelsize=14)
    
    ax.set_title('Enrichment Consistency Across Clones\n(Significantly Enriched Cell Types, FDR < 0.05)', 
                 fontweight='bold', fontsize=20, pad=15)
    ax.set_ylabel('Cell Type', fontweight='bold', fontsize=18, labelpad=10)
    ax.set_xlabel('CN Clone', fontweight='bold', fontsize=18, labelpad=10)
    ax.tick_params(axis='both', labelsize=14)
    
    plt.tight_layout()
    
    # 保存热图 PNG 和 PDF
    heatmap_png = os.path.join(output_dir, 'CN278_consistency_heatmap.png')
    heatmap_pdf = os.path.join(output_dir, 'CN278_consistency_heatmap.pdf')
    plt.savefig(heatmap_png, dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(heatmap_pdf, format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\n  ✓ Consistency heatmap: {heatmap_png}")
    print(f"  ✓ Consistency heatmap: {heatmap_pdf}")
    
    consistency_matrix.to_csv(os.path.join(output_dir, 'CN278_consistency_matrix.csv'))
    print(f"  ✓ Consistency matrix: CN278_consistency_matrix.csv")


# ============================================================================
# 6. 最终报告
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nOutput directory: {output_dir}/")
print(f"\nIndividual clone analysis:")
for cn in all_target_cns:
    print(f"  • CN278_{cn}_volcano.png / .pdf       - {cn} volcano plot")
    print(f"  • CN278_{cn}_results.xlsx              - {cn} detailed results")
print(f"\nCombined analysis:")
print(f"  • CN278_combined_enriched_only.png / .pdf  - Enriched only")
print(f"  • CN278_combined_full.png / .pdf           - Full volcano")
print(f"  • CN278_combined_results.xlsx              - Combined results")
print(f"\nConsistency analysis:")
print(f"  • CN278_consistency_heatmap.png / .pdf     - Heatmap")
print(f"  • CN278_consistency_matrix.csv             - Matrix")
print(f"\nKey finding:")
print(f"  {len(common_enriched)} cell types commonly enriched across CN2, CN7, CN8")
if common_enriched:
    for ct in sorted(common_enriched):
        print(f"    • {ct}")
print("=" * 80)
