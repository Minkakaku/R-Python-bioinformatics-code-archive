import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
import multiprocessing
from spatialdata_io import xenium
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans
from scipy.spatial import KDTree
from scipy.stats import mannwhitneyu, ttest_ind

import warnings
warnings.filterwarnings("ignore")


# 多线程样本数据读入
def load_xenium_data(data_dir, sample_info, out_dir):
    adata = sc.read_h5ad('./adata_raw.h5ad')
    sc.pp.subsample(adata, 0.01)
    out_dir = os.path.join(out_dir, '01-Data_QC')
    os.makedirs(out_dir, exist_ok=True)
    print(out_dir, '-----')

    # 仅执行 QC 指标计算与派生变量，无绘图与文件写入
    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100), inplace=True)
    adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes_by_counts']) / np.log10(adata.obs['total_counts'])

    return adata

    # def sd_read_xenium(raw_dir, fov_name, sample_name, sample_group, result_list):
    #     sdata = xenium(raw_dir, cells_as_circles=True, n_jobs=6)
    #     if hasattr(sdata, 'tables'):
    #         for name, table in sdata.tables.items():
    #             table.obs['fov'] = fov_name
    #             table.obs['sample'] = sample_name
    #             table.obs['group'] = sample_group
    #             result_list.append(table)                        

    # with multiprocessing.Manager() as manager:
    #     result_list = manager.list()  # 使用Manager创建共享列表
    #     processes = []
    #     with open(sample_info, 'r') as f:
    #         line = f.readline()
    #         for line in f.readlines():
    #             line = line.strip().split('\t')
    #             raw_name = line[0]
    #             fov_name = line[1]
    #             sample_name = line[2]
    #             sample_group = line[3]
    #             p = multiprocessing.Process(
    #                     target=sd_read_xenium, 
    #                     args=(os.path.join(data_dir, raw_name), fov_name, sample_name, sample_group, result_list,)
    #                 )
    #             processes.append(p)
    #             p.start()

    #     for p in processes:
    #         p.join()


# ###################
#     result_list = []
#     with open(sample_info, 'r') as f:
#         line = f.readline()
#         for line in f.readlines():
#             line = line.strip().split('\t')
#             raw_name = line[0]
#             fov_name = line[1]
#             sample_name = line[2]
#             sample_group = line[3]
#             sdata = xenium(os.path.join(data_dir, raw_name), cells_as_circles=True, n_jobs=6)
#             if hasattr(sdata, 'tables'):
#                 for name, table in sdata.tables.items():
#                     table.obs['fov'] = fov_name
#                     table.obs['sample'] = sample_name
#                     table.obs['group'] = sample_group
#                     result_list.append(table)
# ######################

#     adata = result_list[0].concatenate(result_list[1:], index_unique=None)
#     adata.obs_names_make_unique()
#     adata.var_names_make_unique()
#     del adata.obs['batch']

#     if adata[adata.obs['sample'].str.contains('.csv_'), :].shape[0] > 0:
#         adata_single_sample = adata[~adata.obs['sample'].str.contains('.csv_')]
#         adata_multi_sample = adata[adata.obs['sample'].str.contains('.csv_')]
#         del adata_multi_sample.obs['sample']

#         tmp = []
#         sample2group = {}
#         with open(sample_info, 'r') as f:
#             line = f.readline()
#             for line in f.readlines():
#                 line = line.strip().split('\t')

#                 if ',' in line[2] and ',' in line[3]:
#                     sample2group.update(zip([i.split('.csv_')[1] for i in line[2].split(',')], line[3].split(',')))

#                 if '.csv_' in line[2] and line[4]:
#                     roi = line[2].split(',')
#                     for i in roi:
#                         i = i.split('.csv_')
#                         df = pd.read_csv(os.path.join(line[4], f'{i[0]}.csv'), skiprows=2)
#                         df.rename(columns={'Cell ID': 'cell_id'}, inplace=True)
#                         df['sample'] = i[1]
#                         df['fov'] = line[1]
#                         tmp.append(df[['cell_id', 'sample', 'fov']])

#         all_roi = pd.concat(tmp, axis=0)

#         merge_obs = pd.merge(adata_multi_sample.obs, all_roi, on=['cell_id', 'fov'], how='inner')

#         ### 为了防止cell_id重复
#         merge_obs['cell_id_fov'] = merge_obs['cell_id'].astype(str) + '-' + merge_obs['fov'].astype(str)
#         adata_multi_sample.obs['cell_id_fov'] = adata_multi_sample.obs['cell_id'].astype(str) + '-' + adata_multi_sample.obs['fov'].astype(str)
#         adata_multi_sample = adata_multi_sample[adata_multi_sample.obs['cell_id_fov'].isin(merge_obs['cell_id_fov'])]
#         del merge_obs['cell_id_fov']

#         adata_multi_sample.obs = merge_obs.copy()
#         adata_all = adata_single_sample.concatenate(adata_multi_sample, index_unique=None)
#         adata_all.obs['group'] = adata_all.obs.apply(lambda row: sample2group[row['sample']] if row['sample'] in sample2group else row['group'], axis=1)
#         adata_all.obs['cell_id_sample'] = adata_all.obs['cell_id'].astype(str) + '-' + adata_all.obs['sample'].astype(str)
#         del adata_all.obs['batch']
#         adata_all.obs.reset_index(inplace=True, drop=True)
#         adata_all.obs.index = adata_all.obs.index.astype(str)

#         # Data QC
#         sc.pp.calculate_qc_metrics(adata_all, percent_top=(10, 50, 100), inplace=True)
#         adata_all.obs['log10GenesPerUMI'] = np.log10(adata_all.obs['n_genes_by_counts']) / np.log10(adata_all.obs['total_counts'])
        
#         adata_all.write('./adata_raw.h5ad')

#         return adata_all
#     else:
#         adata.obs.reset_index(inplace=True, drop=True)


# 移除前端HTML生成函数，保持本文件仅包含计算逻辑


def data_preprocess(
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
    """纯计算逻辑：执行过滤、聚类与特征基因计算，不进行任何绘图或HTML生成。
    返回一个字典，供前端包装函数使用。
    """
    try:
        adata_processed = adata.copy()
        filter_info = None
        adata_cluster_result = None
        marker_info = None

        # 步骤1: 细胞过滤（计算统计与返回过滤后的 AnnData）
        if "细胞过滤" in analysis_options:
            if filter_params_df is None or len(filter_params_df) == 0:
                return {
                    'error': "过滤参数为空"
                }

            original_cells = adata_processed.n_obs
            adata_list = []
            sample_stats = []

            for _, row in filter_params_df.iterrows():
                sample = row['sample']
                min_genes = float(row['min_genes'])
                max_genes = float(row['max_genes'])
                min_genes_per_umi = float(row['min_GenesPerUMI'])
                max_genes_per_umi = float(row['max_GenesPerUMI'])

                tmp = adata_processed[adata_processed.obs['sample'] == sample]
                sample_original = tmp.n_obs

                tmp = tmp[
                    (tmp.obs['n_genes_by_counts'] > min_genes) &
                    (tmp.obs['n_genes_by_counts'] < max_genes) &
                    (tmp.obs['log10GenesPerUMI'] > min_genes_per_umi) &
                    (tmp.obs['log10GenesPerUMI'] < max_genes_per_umi)
                ]

                sample_filtered = tmp.n_obs
                sample_removed = sample_original - sample_filtered
                sample_removal_rate = (sample_removed / sample_original) * 100 if sample_original > 0 else 0

                sample_stats.append({
                    'sample': sample,
                    'original': sample_original,
                    'filtered': sample_filtered,
                    'removed': sample_removed,
                    'removal_rate': sample_removal_rate
                })

                if tmp.n_obs > 0:
                    adata_list.append(tmp)

            if adata_list:
                adata_processed = anndata.concat(adata_list, join='inner', index_unique=None)
                filtered_cells = adata_processed.n_obs
                removed_cells = original_cells - filtered_cells
                removal_rate = (removed_cells / original_cells) * 100
            else:
                return {
                    'error': "过滤后没有剩余细胞"
                }

            filter_info = {
                'original_cells': int(original_cells),
                'filtered_cells': int(filtered_cells),
                'removed_cells': int(removed_cells),
                'removal_rate': float(removal_rate),
                'sample_stats': sample_stats,
                'adata_processed': adata_processed
            }

        # 步骤2: 降维聚类（仅计算，不绘图）
        if "降维聚类" in analysis_options:
            # 使用过滤后的对象或原始对象
            ad_for_cluster = filter_info['adata_processed'] if filter_info else adata_processed

            sc.pp.normalize_total(ad_for_cluster)
            sc.pp.log1p(ad_for_cluster)
            sc.pp.highly_variable_genes(ad_for_cluster, flavor='seurat', batch_key='sample', n_top_genes=int(hvgs))
            ad_for_cluster = ad_for_cluster[:, ad_for_cluster.var['highly_variable']]
            sc.pp.regress_out(ad_for_cluster, ['total_counts'])
            sc.pp.scale(ad_for_cluster, max_value=10)
            sc.tl.pca(ad_for_cluster, svd_solver='arpack')
            # harmony / bbknn
            sc.external.pp.harmony_integrate(ad_for_cluster, 'sample', adjusted_basis='X_pca')
            sc.pp.neighbors(ad_for_cluster, metric='cosine', n_neighbors=int(n_neighbors), n_pcs=int(npcs))
            sc.tl.umap(ad_for_cluster, min_dist=0.3, spread=1.5)
            sc.tl.leiden(ad_for_cluster, resolution=float(resolution), key_added=f'res_{resolution}')

            adata_cluster_result = ad_for_cluster

        # 步骤3: 特征基因（仅计算排名与表，不绘图）
        if "特征基因筛选" in analysis_options:
            adata_for_marker = adata_cluster_result if adata_cluster_result is not None else (filter_info['adata_processed'] if filter_info else adata_processed)
            if adata_for_marker is None:
                return {
                    'error': "无法进行特征基因筛选：缺少有效的 AnnData"
                }

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

            groupby_key = find_cluster_key(adata_for_marker)
            if groupby_key is None:
                try:
                    obs_df = getattr(adata_for_marker, 'obs')
                    obs_cols = list(obs_df.columns)
                    annot_candidates = [
                        "celltype", "cell_type", "cell_type_label",
                        "annotation", "annot", "labels", "Label", "CellType"
                    ]
                    for c in annot_candidates:
                        if c in obs_cols and int(obs_df[c].nunique()) >= 2:
                            groupby_key = c
                            break
                except Exception:
                    groupby_key = None

            if groupby_key is None:
                sc.pp.neighbors(adata_for_marker)
                sc.tl.leiden(adata_for_marker, resolution=0.5)
                groupby_key = "leiden"

            try:
                getattr(adata_for_marker, 'uns', {}).pop('rank_genes_groups', None)
            except Exception:
                pass
            sc.tl.rank_genes_groups(adata_for_marker, groupby=groupby_key, use_raw=False, method=str(feature_method), pts=True)

            # 计算Top表
            try:
                df_all = sc.get.rank_genes_groups_df(adata_for_marker, None)
                lfc_thr = float(feature_lfc or 0.25)
                feature_table_df = (
                    df_all[df_all["logfoldchanges"].abs() >= lfc_thr]
                    .groupby("group", as_index=False)
                    .apply(lambda g: g.head(int(feature_top or 5)))
                    .reset_index(drop=True)
                )
                feature_table_df = feature_table_df[["group", "names", "scores", "logfoldchanges", "pvals_adj"]]
                feature_table_df.columns = ["group", "gene", "score", "logfoldchanges", "pvals_adj"]
            except Exception:
                feature_table_df = pd.DataFrame(columns=["group", "gene", "score", "logfoldchanges", "pvals_adj"]) 

            marker_info = {
                'adata_for_marker': adata_for_marker,
                'groupby_key': groupby_key,
                'feature_table_df': feature_table_df
            }

        return {
            'filter': filter_info,
            'cluster': adata_cluster_result,
            'marker': marker_info,
            'params': {
                'hvgs': int(hvgs),
                'npcs': int(npcs),
                'n_neighbors': int(n_neighbors),
                'resolution': float(resolution)
            }
        }

    except Exception as e:
        return {
            'error': str(e)
        }
#         adata.obs.index = adata.obs.index.astype(str)
#         # Data QC
#         sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100), inplace=True)
#         adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes_by_counts']) / np.log10(adata.obs['total_counts'])
#         return adata


