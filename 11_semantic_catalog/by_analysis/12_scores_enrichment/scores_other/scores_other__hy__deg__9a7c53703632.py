import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse


target_genes = ["AKT1","AKT2","AKT3","BCL2L1","CBL","CBLB","CBLC","CCND1","CCND2","CCND3","CISH","CLCF1","CNTF","CNTFR","CREBBP","CRLF2","CSF2","CSF2RA","CSF2RB","CSF3","CSF3R","CSH1","CTF1","EP300","EPO","EPOR","GH1","GH2","GHR","GRB2","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNE","IFNG","IFNGR1","IFNGR2","IFNK","IFNL1","IFNL2","IFNL3","IFNLR1","IFNW1","IL10","IL10RA","IL10RB","IL11","IL11RA","IL12A","IL12B","IL12RB1","IL12RB2","IL13","IL13RA1","IL13RA2","IL15","IL15RA","IL19","IL2","IL20","IL20RA","IL20RB","IL21","IL21R","IL22","IL22RA1","IL22RA2","IL23A","IL23R","IL24","IL26","IL2RA","IL2RB","IL2RG","IL3","IL3RA","IL4","IL4R","IL5","IL5RA","IL6","IL6R","IL6ST","IL7","IL7R","IL9","IL9R","IRF9","JAK1","JAK2","JAK3","LEP","LEPR","LIF","LIFR","MPL","MYC","OSM","OSMR","PIAS1","PIAS2","PIAS3","PIAS4","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","PIM1","PRL","PRLR","PTPN11","PTPN6","SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS7","SOS1","SOS2","SPRED1","SPRED2","SPRY1","SPRY2","SPRY3","SPRY4","STAM","STAM2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","TPO","TSLP","TYK2"]  # 用户提供的完整列表
# 预处理与数据校验
def preprocess_data(adata_path):
    old_adata = sc.read_h5ad(adata_path)
    # 筛选单个样本
    # 筛选多个样本
    selected_samples = ["d2921","dwwl1_2","dwwl2_1","dxlx3_1","dxlx3_2","ndxz","nfxl","nwjy"]
    adata = old_adata[old_adata.obs["sample"].isin(selected_samples), :].copy()
    # 强制唯一化索引[6](@ref)
    adata.obs_names_make_unique(join="-dup")
    adata.var_names_make_unique()
    
    # 校验必须列存在[1,4](@ref)
    required_cols = ['celltype', 'group']
    missing = [c for c in required_cols if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing columns: {missing}")
    
    # 筛选目标基因
    target_genes = ["AKT1","AKT2","AKT3","BCL2L1","CBL","CBLB","CBLC","CCND1","CCND2","CCND3","CISH","CLCF1","CNTF","CNTFR","CREBBP","CRLF2","CSF2","CSF2RA","CSF2RB","CSF3","CSF3R","CSH1","CTF1","EP300","EPO","EPOR","GH1","GH2","GHR","GRB2","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNE","IFNG","IFNGR1","IFNGR2","IFNK","IFNL1","IFNL2","IFNL3","IFNLR1","IFNW1","IL10","IL10RA","IL10RB","IL11","IL11RA","IL12A","IL12B","IL12RB1","IL12RB2","IL13","IL13RA1","IL13RA2","IL15","IL15RA","IL19","IL2","IL20","IL20RA","IL20RB","IL21","IL21R","IL22","IL22RA1","IL22RA2","IL23A","IL23R","IL24","IL26","IL2RA","IL2RB","IL2RG","IL3","IL3RA","IL4","IL4R","IL5","IL5RA","IL6","IL6R","IL6ST","IL7","IL7R","IL9","IL9R","IRF9","JAK1","JAK2","JAK3","LEP","LEPR","LIF","LIFR","MPL","MYC","OSM","OSMR","PIAS1","PIAS2","PIAS3","PIAS4","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","PIM1","PRL","PRLR","PTPN11","PTPN6","SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS7","SOS1","SOS2","SPRED1","SPRED2","SPRY1","SPRY2","SPRY3","SPRY4","STAM","STAM2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","TPO","TSLP","TYK2"]  # 用户提供的完整列表
    gene_mask = adata.var_names.isin(target_genes)
    adata = adata[:, gene_mask].copy()
    return adata

# 分细胞类型执行差异分析
def celltype_specific_de(adata):
    results = []
    celltypes = adata.obs['celltype'].unique()
    
    for ct in celltypes:
        # 提取当前细胞类型子集
        ct_adata = adata[adata.obs['celltype'] == ct].copy()
        if ct_adata.n_obs < 10:  # 跳过细胞数过少的类型
            print(f"Skipping {ct} (cells < 10)")
            continue
            
        try:
            # 执行组间差异分析[8](@ref)
            sc.tl.rank_genes_groups(
                ct_adata,
                groupby='group',
                groups=['disease'],
                reference='normal',
                method='wilcoxon',
                pts=True,
                key_added=f'deg_{ct}'
            )
            
            # 提取结果并添加元数据
            deg_df = sc.get.rank_genes_groups_df(ct_adata, group='disease', key=f'deg_{ct}')
            deg_df.insert(0, 'celltype', ct)
            
            # 合并pts值并筛选[3,8](@ref)
            pts_vals = ct_adata.uns[f'deg_{ct}']['pts']['disease']
            deg_df['pts'] = deg_df['names'].map(pts_vals)
            
            # 过滤条件：adj.p<0.05, |logFC|>1, pts>0.1
            deg_df = deg_df[
                (deg_df.pvals_adj < 0.05) &
                (abs(deg_df.logfoldchanges) > 1) &
                (deg_df.pts > 0.1) &
                (deg_df['names'].isin(target_genes))
            ]
            
            results.append(deg_df)
            
        except Exception as e:
            print(f"Failed in {ct}: {str(e)}")
            continue
    
    return pd.concat(results)

# 保存结果文件
def save_results(combined_df):
    # 按细胞类型分割
    for (ct, ct_df) in combined_df.groupby('celltype'):
        # 处理特殊字符
        ct_clean = ct.replace(" ", "_")
        
        # 上下调分离
        up = ct_df[ct_df.logfoldchanges > 0]
        down = ct_df[ct_df.logfoldchanges < 0]
        
        # 保存文件[7](@ref)
        up.to_csv(f"upregulated_{ct_clean}.csv", index=False)
        down.to_csv(f"downregulated_{ct_clean}.csv", index=False)

if __name__ == "__main__":
    adata = preprocess_data("adata_celltype.h5ad")
    de_results = celltype_specific_de(adata)
    
    if not de_results.empty:
        save_results(de_results)
        print(f"Analyzed {len(de_results.celltype.unique())} cell types")
    else:
        print("No significant genes found")
