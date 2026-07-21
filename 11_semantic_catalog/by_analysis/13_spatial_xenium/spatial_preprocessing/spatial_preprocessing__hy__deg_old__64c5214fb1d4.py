import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# 尝试导入gseapy用于富集分析
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    print("gseapy未安装，富集分析功能将不可用。请运行: pip install gseapy")
    GSEAPY_AVAILABLE = False

def save_raw_counts(adata):
    """将原始计数矩阵保存到raw属性中"""
    # 复制当前数据到raw属性
    adata.raw = adata.copy()
    print(f"已将原始计数矩阵保存到adata.raw")
    return adata

def normalize_data(adata, target_sum=1e4):
    """标准化数据"""
    print("开始数据标准化...")
    
    # 复制原始adata用于标准化
    adata_norm = adata.copy()
    
    # 标准化：使每个细胞的总计数相同
    sc.pp.normalize_total(adata_norm, target_sum=target_sum)
    
    # 对数变换
    sc.pp.log1p(adata_norm)
    
    # 保存标准化后的数据到原始adata对象的新层
    adata.layers['log1p_normalized'] = adata_norm.X.copy()
    
    print(f"数据标准化完成，标准化后的数据保存在adata.layers['log1p_normalized']")
    print(f"每个细胞的总计数标准化到{target_sum}")
    
    return adata

def get_comparison_groups(adata, group_col='group_final'):
    """获取所有两两比较组"""
    # 获取唯一的分组
    groups = adata.obs[group_col].unique()
    groups = [g for g in groups if str(g) != 'nan']  # 移除NaN值
    
    print(f"找到的分组: {groups}")
    
    # 生成所有两两组合
    comparisons = list(combinations(groups, 2))
    print(f"生成的比较组 ({len(comparisons)}个): {comparisons}")
    
    return groups, comparisons

def perform_dea_for_celltype(adata, celltype, group1, group2, 
                             celltype_col='celltype', group_col='group_final',
                             method='wilcoxon', n_genes=50):
    """对特定细胞类型执行差异表达分析"""
    
    # 提取特定细胞类型的数据
    celltype_data = adata[adata.obs[celltype_col] == celltype].copy()
    
    # 提取两个比较组的数据
    group1_data = celltype_data[celltype_data.obs[group_col] == group1].copy()
    group2_data = celltype_data[celltype_data.obs[group_col] == group2].copy()
    
    # 检查是否有足够的数据
    if group1_data.n_obs < 3 or group2_data.n_obs < 3:
        print(f"  警告: {celltype}中{group1}(n={group1_data.n_obs}) vs {group2}(n={group2_data.n_obs}) 细胞数不足，跳过")
        return None
    
    # 合并两个组的数据
    dea_data = sc.concat([group1_data, group2_data])
    
    # 创建新的分组标签
    dea_data.obs['comparison_group'] = pd.Categorical(
        [group1] * group1_data.n_obs + [group2] * group2_data.n_obs
    )
    
    # 使用原始计数进行差异分析
    dea_data.X = dea_data.raw.X
    
    # 对比较数据进行标准化
    sc.pp.normalize_total(dea_data, target_sum=1e4)
    sc.pp.log1p(dea_data)
    
    # 执行差异表达分析
    sc.tl.rank_genes_groups(
        dea_data, 
        groupby='comparison_group', 
        method=method,
        n_genes=n_genes,
        use_raw=False
    )
    
    # 获取差异表达结果
    de_results = sc.get.rank_genes_groups_df(dea_data, group=group1)
    
    # 添加元数据
    de_results['celltype'] = celltype
    de_results['comparison'] = f"{group1}_vs_{group2}"
    de_results['group1'] = group1
    de_results['group2'] = group2
    
    return de_results

def perform_dea_all_celltypes(adata, comparisons, output_dir, 
                              celltype_col='celltype', group_col='group_final'):
    """对所有细胞类型执行差异表达分析"""
    
    # 获取所有细胞类型
    celltypes = adata.obs[celltype_col].unique()
    celltypes = [ct for ct in celltypes if str(ct) != 'nan']  # 移除NaN值
    
    print(f"找到的细胞类型: {celltypes}")
    
    all_de_results = []
    
    for group1, group2 in comparisons:
        comparison_name = f"{group1}_vs_{group2}"
        comp_output_dir = os.path.join(output_dir, comparison_name)
        os.makedirs(comp_output_dir, exist_ok=True)
        
        print(f"\n处理比较组: {comparison_name}")
        print(f"输出目录: {comp_output_dir}")
        
        comp_results = []
        
        for celltype in celltypes:
            print(f"  分析细胞类型: {celltype}")
            
            de_results = perform_dea_for_celltype(
                adata, celltype, group1, group2, 
                celltype_col=celltype_col, group_col=group_col
            )
            
            if de_results is not None and len(de_results) > 0:
                comp_results.append(de_results)
                
                # 保存单个细胞类型的结果
                ct_output_dir = os.path.join(comp_output_dir, f"{celltype}")
                os.makedirs(ct_output_dir, exist_ok=True)
                
                de_file = os.path.join(ct_output_dir, f"de_results_{celltype}.csv")
                de_results.to_csv(de_file, index=False)
                print(f"    已保存到: {de_file}")
        
        if comp_results:
            # 合并所有细胞类型的结果
            comp_all_results = pd.concat(comp_results, ignore_index=True)
            
            # 保存完整结果
            comp_file = os.path.join(comp_output_dir, f"de_results_all_celltypes.csv")
            comp_all_results.to_csv(comp_file, index=False)
            
            all_de_results.append(comp_all_results)
            
            # 创建汇总统计
            summary_stats = comp_all_results.groupby('celltype').agg({
                'names': 'count',
                'logfoldchanges': ['mean', 'std'],
                'pvals_adj': lambda x: (x < 0.05).sum()
            }).round(4)
            
            summary_file = os.path.join(comp_output_dir, f"de_summary_statistics.csv")
            summary_stats.to_csv(summary_file)
            
            print(f"  已保存汇总统计: {summary_file}")
    
    if all_de_results:
        # 合并所有比较组的结果
        all_results = pd.concat(all_de_results, ignore_index=True)
        
        # 保存所有结果
        all_file = os.path.join(output_dir, "all_de_results.csv")
        all_results.to_csv(all_file, index=False)
        
        print(f"\n所有差异表达分析结果已保存到: {all_file}")
    
    return all_de_results

def perform_enrichment_analysis(de_results, output_dir, 
                                organism='Human', 
                                gene_sets=['KEGG_2019_Human', 'GO_Biological_Process_2021']):
    """执行富集分析"""
    
    if not GSEAPY_AVAILABLE:
        print("gseapy未安装，跳过富集分析")
        return None
    
    # 过滤显著差异表达基因
    sig_genes = de_results[
        (de_results['pvals_adj'] < 0.05) & 
        (abs(de_results['logfoldchanges']) > 1)
    ]
    
    if len(sig_genes) < 5:
        print(f"  显著差异基因数量不足({len(sig_genes)})，跳过富集分析")
        return None
    
    # 获取基因列表
    gene_list = sig_genes['names'].tolist()
    
    # 执行富集分析
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=None,
            cutoff=0.05
        )
        
        # 保存结果
        for gene_set, result_df in enr.results.items():
            if result_df is not None and len(result_df) > 0:
                # 创建文件名
                safe_gene_set = gene_set.replace(' ', '_').replace('/', '_')
                output_file = os.path.join(output_dir, f"enrichment_{safe_gene_set}.csv")
                result_df.to_csv(output_file, index=False)
                
                print(f"    已保存富集分析结果: {output_file}")
                
                # 创建富集分析图
                if len(result_df) > 0:
                    try:
                        # 选择前10个最显著的条目
                        top_terms = result_df.head(10)
                        
                        # 创建条形图
                        plt.figure(figsize=(10, 8))
                        plt.barh(range(len(top_terms)), 
                                -np.log10(top_terms['Adjusted P-value'].values))
                        plt.yticks(range(len(top_terms)), top_terms['Term'].values)
                        plt.xlabel('-log10(Adjusted P-value)')
                        plt.title(f'{gene_set} Enrichment - Top 10 Terms')
                        plt.tight_layout()
                        
                        # 保存图形
                        plot_file = os.path.join(output_dir, f"enrichment_{safe_gene_set}.png")
                        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                        plt.close()
                        
                        print(f"    已保存富集分析图: {plot_file}")
                    except Exception as e:
                        print(f"    创建富集分析图失败: {e}")
        
        return enr
        
    except Exception as e:
        print(f"  富集分析失败: {e}")
        return None

def run_enrichment_for_all_comparisons(de_results_list, comparisons, output_dir):
    """对所有比较组执行富集分析"""
    
    if not de_results_list:
        print("没有差异表达结果可用于富集分析")
        return
    
    for i, (comp_results, (group1, group2)) in enumerate(zip(de_results_list, comparisons)):
        comparison_name = f"{group1}_vs_{group2}"
        comp_output_dir = os.path.join(output_dir, comparison_name)
        
        print(f"\n对比较组 {comparison_name} 执行富集分析")
        
        # 对每个细胞类型执行富集分析
        celltypes = comp_results['celltype'].unique()
        
        for celltype in celltypes:
            print(f"  细胞类型: {celltype}")
            
            # 提取该细胞类型的结果
            ct_results = comp_results[comp_results['celltype'] == celltype]
            
            # 创建输出目录
            ct_output_dir = os.path.join(comp_output_dir, f"{celltype}")
            os.makedirs(ct_output_dir, exist_ok=True)
            
            # 执行富集分析
            perform_enrichment_analysis(ct_results, ct_output_dir)
        
        # 对所有细胞类型的合并结果执行富集分析
        print(f"  所有细胞类型合并:")
        comp_all_output_dir = os.path.join(comp_output_dir, "all_celltypes_combined")
        os.makedirs(comp_all_output_dir, exist_ok=True)
        
        # 合并所有显著基因
        sig_genes_all = []
        for celltype in celltypes:
            ct_results = comp_results[comp_results['celltype'] == celltype]
            sig_genes = ct_results[
                (ct_results['pvals_adj'] < 0.05) & 
                (abs(ct_results['logfoldchanges']) > 1)
            ]['names'].tolist()
            sig_genes_all.extend(sig_genes)
        
        # 去重
        sig_genes_all = list(set(sig_genes_all))
        
        if len(sig_genes_all) >= 5:
            # 创建临时DataFrame用于富集分析
            temp_df = pd.DataFrame({
                'names': sig_genes_all,
                'logfoldchanges': [1] * len(sig_genes_all),  # 占位符
                'pvals_adj': [0.001] * len(sig_genes_all)   # 占位符
            })
            
            perform_enrichment_analysis(temp_df, comp_all_output_dir)
        else:
            print(f"    合并的显著基因数量不足({len(sig_genes_all)})")

def create_summary_report(adata, comparisons, output_dir):
    """创建分析摘要报告"""
    
    report_file = os.path.join(output_dir, "analysis_summary.txt")
    
    with open(report_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("差异表达分析摘要报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("1. 数据基本信息\n")
        f.write(f"   总细胞数: {adata.n_obs}\n")
        f.write(f"   总基因数: {adata.n_vars}\n")
        f.write(f"   观测变量(obs): {list(adata.obs.columns)}\n\n")
        
        f.write("2. 分组信息\n")
        groups = adata.obs['group_final'].unique()
        f.write(f"   分组: {list(groups)}\n")
        
        for group in groups:
            n_cells = (adata.obs['group_final'] == group).sum()
            f.write(f"     {group}: {n_cells} 个细胞\n")
        f.write("\n")
        
        f.write("3. 细胞类型分布\n")
        celltypes = adata.obs['celltype'].unique()
        for celltype in celltypes:
            n_cells = (adata.obs['celltype'] == celltype).sum()
            f.write(f"     {celltype}: {n_cells} 个细胞\n")
        f.write("\n")
        
        f.write("4. 比较组信息\n")
        for i, (group1, group2) in enumerate(comparisons, 1):
            f.write(f"    {i}. {group1} vs {group2}\n")
        f.write("\n")
        
        f.write("5. 输出目录结构\n")
        f.write("   主输出目录/\n")
        f.write("   ├── group1_vs_group2/\n")
        f.write("   │   ├── celltype1/\n")
        f.write("   │   │   ├── de_results_celltype1.csv\n")
        f.write("   │   │   ├── enrichment_KEGG_2019_Human.csv\n")
        f.write("   │   │   ├── enrichment_GO_Biological_Process_2021.csv\n")
        f.write("   │   │   └── ...\n")
        f.write("   │   ├── celltype2/\n")
        f.write("   │   └── ...\n")
        f.write("   ├── group1_vs_group3/\n")
        f.write("   └── ...\n")
        f.write("\n")
        
        f.write("6. 分析参数\n")
        f.write("   标准化目标总和: 10,000 counts/cell\n")
        f.write("   差异表达分析方法: Wilcoxon rank-sum test\n")
        f.write("   显著差异阈值: |log2FC| > 1 且 adj.p < 0.05\n")
        f.write("   富集分析数据库: KEGG, GO Biological Process\n")
    
    print(f"分析摘要报告已保存到: {report_file}")

def main():
    # 设置路径
    input_file = "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/All_filter.h5ad"
    output_base_dir = "./DEA_Analysis"
    
    # 创建输出目录
    os.makedirs(output_base_dir, exist_ok=True)
    
    print("=" * 60)
    print("差异表达分析流程开始")
    print("=" * 60)
    
    # 1. 读取数据
    print("\n1. 读取数据...")
    adata = sc.read_h5ad(input_file)
    print(f"   读取完成: {adata.n_obs} 细胞, {adata.n_vars} 基因")
    
    # 检查必要的列是否存在
    required_cols = ['group_final', 'celltype']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"adata.obs中缺少必需的列: {col}")
    
    # 2. 保存原始计数矩阵
    print("\n2. 保存原始计数矩阵...")
    adata = save_raw_counts(adata)
    
    # 3. 标准化数据
    print("\n3. 标准化数据...")
    adata = normalize_data(adata, target_sum=1e4)
    
    # 4. 获取比较组
    print("\n4. 确定比较组...")
    groups, comparisons = get_comparison_groups(adata, group_col='group_final')
    
    # 5. 执行差异表达分析
    print("\n5. 执行差异表达分析...")
    de_results_list = perform_dea_all_celltypes(
        adata, comparisons, output_base_dir,
        celltype_col='celltype', group_col='group_final'
    )
    
    # 6. 执行富集分析
    if de_results_list and GSEAPY_AVAILABLE:
        print("\n6. 执行富集分析...")
        run_enrichment_for_all_comparisons(de_results_list, comparisons, output_base_dir)
    elif not GSEAPY_AVAILABLE:
        print("\n6. 跳过富集分析 (gseapy未安装)")
    
    # 7. 创建摘要报告
    print("\n7. 创建分析摘要...")
    create_summary_report(adata, comparisons, output_base_dir)
    
    # 8. 保存处理后的数据
    print("\n8. 保存处理后的数据...")
    processed_file = os.path.join(output_base_dir, "processed_data.h5ad")
    adata.write_h5ad(processed_file)
    print(f"   处理后的数据已保存到: {processed_file}")
    
    print("\n" + "=" * 60)
    print("分析完成!")
    print(f"所有结果已保存到: {output_base_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()
