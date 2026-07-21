import os
import pandas as pd
import matplotlib.pyplot as plt
import omicverse as ov
import scanpy as sc

import os
import pandas as pd
import matplotlib.pyplot as plt
import omicverse as ov
import scanpy as sc


def collect_statistical_files(directory='.', exclude_second=True):
    """
    收集指定目录下特定的统计文件
    
    参数:
        directory (str): 要搜索的目录路径，默认为当前目录
        exclude_second (bool): 是否排除第二个匹配的文件
    
    返回:
        list: 匹配的文件列表
    """
    target_prefixes = [
        "statistical_analysis_deconvoluted",
        "statistical_analysis_means",
        "statistical_analysis_pvalues"
    ]
    
    matched_files = []
    
    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            for prefix in target_prefixes:
                if filename.startswith(prefix):
                    filepath = os.path.join(directory, filename)
                    matched_files.append(filepath)
                    break
    
    if exclude_second and len(matched_files) > 1:
        return matched_files[:1] + matched_files[2:]
    
    return matched_files


def process_cellchat_analysis(color_dict, base_path='.'):
    """
    处理CellChat分析流程，包括数据保存、可视化图形生成和通路分析
    
    参数:
        color_dict: 颜色字典，用于可视化
        base_path (str): 基础目录路径，默认为当前目录
    返回:
        dict: 处理结果的统计信息
    """
    
    def save_fig(fig_func, sample_path, *args, filename, **kwargs):
        """
        内部函数：保存图形文件
        """
        try:
            fig, ax = fig_func(*args, **kwargs)
            fig.savefig(os.path.join(sample_path, f"{filename}.png"), bbox_inches='tight')
            fig.savefig(os.path.join(sample_path, f"{filename}.pdf"), bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            print(f"警告: 保存图形 {filename} 失败: {str(e)}")
    

    
    # 获取样本列表
    samples = [item for item in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, item))]
    results = {}
    
    for sample in samples:
        try:
            sample_path = os.path.join(base_path, sample)
            print(f"正在处理样本: {sample}")
            
            # 复制当前样本的数据 - 使用参数指定的列名
            #tmp = adata[adata.obs[group_col] == sample].copy()
            
            # 加载CellChat结果
            cpdb_results = ov.utils.load(os.path.join(sample_path, 'gex_cpdb_test.pkl'))
            adata_cpdb = ov.read(os.path.join(sample_path, 'gex_cpdb_ad.h5ad'))
            
            # 保存TXT信息文件 - 直接使用to_csv方法[1,4](@ref)
            means_df = pd.DataFrame(adata_cpdb.layers['means'], 
                                   index=adata_cpdb.obs_names, 
                                   columns=adata_cpdb.var_names)
            pvals_df = pd.DataFrame(adata_cpdb.layers['pvalues'], 
                                   index=adata_cpdb.obs_names, 
                                   columns=adata_cpdb.var_names)
            
            long_df = means_df.reset_index().melt(id_vars='index', 
                                                 var_name='interaction_pair', 
                                                 value_name='mean_score')
            long_df = pd.merge(long_df, 
                             pvals_df.reset_index().melt(id_vars='index', 
                                                       var_name='interaction_pair', 
                                                       value_name='p_value'), 
                             on=['index', 'interaction_pair'])
            long_df = pd.merge(long_df, adata_cpdb.obs.reset_index(), on='index')
            long_df = pd.merge(long_df, adata_cpdb.var.reset_index(), 
                             left_on='interaction_pair', right_on='index')
            
            final_df = long_df.rename(columns={'index_x': 'cell_id'}).drop('index_y', axis=1, errors='ignore')
            
            # 直接使用to_csv保存为TXT文件，使用制表符分隔[2,4](@ref)
            txt_filepath = os.path.join(sample_path, 'infor.txt')
            try:
                final_df.to_csv(txt_filepath, sep='\t', index=False)
                print(f"数据已保存为TXT文件: {txt_filepath}")
                
                # 同时保存一个精简版本（前1000行）用于快速查看
                #if len(final_df) > 1000:
                    #sample_filepath = os.path.join(sample_path, 'infor_sample1000.txt')
                    #final_df.head(1000).to_csv(sample_filepath, sep='\t', index=False)
                    #print(f"同时保存了精简版本(前1000行): {sample_filepath}")
            except Exception as e:
                print(f"警告: 保存TXT文件失败: {str(e)}")
            
            
            # omicverse可视化流程
            viz = ov.pl.CellChatViz(adata_cpdb, palette=color_dict)
            count_matrix, weight_matrix = viz.compute_aggregated_network(
                pvalue_threshold=0.05, use_means=True)
            
            # 保存各种圆形图
            save_fig(viz.netVisual_circle, sample_path, weight_matrix, 
                    title="Interaction weights/strength", cmap='Reds', 
                    vertex_size_max=10, figsize=(12, 10), filename='weight_circle')
            
            save_fig(viz.netVisual_circle, sample_path, count_matrix, 
                    title="Number of interactions", cmap='Reds', 
                    vertex_size_max=10, figsize=(12, 10), filename='count_circle')
            
            # 通路分析
            pathway_comm = viz.compute_pathway_communication(method='mean', 
                                                          min_lr_pairs=2, 
                                                          min_expression=0.1)
            sig_pathways = viz.get_significant_pathways_v2(pathway_comm, 
                                                         strength_threshold=0.5, 
                                                         pvalue_threshold=0.05, 
                                                         min_significant_pairs=1)[0]
            
            pathway_output_dir = os.path.join(sample_path, 'pathway')
            os.makedirs(pathway_output_dir, exist_ok=True)
            
            for pathway in sig_pathways:
                try:
                    pathway_clean = pathway.replace("/", "_")
                    save_fig(viz.netVisual_aggregate, pathway_output_dir, 
                            signaling=pathway, layout="circle", figsize=(6, 5), 
                            vertex_size_max=10, filename=f"{pathway_clean}_circle")
                    
                    save_fig(viz.netVisual_chord_cell, pathway_output_dir, 
                            signaling=pathway, group_celltype=None, count_min=1, 
                            figsize=(6, 5), normalize_to_sender=True, 
                            filename=f"{pathway_clean}_chord")
                    
                    try:
                        h = viz.netVisual_heatmap_marsilea(signaling=pathway, 
                                                         color_heatmap="Reds", 
                                                         add_dendrogram=True)
                        h.render()
                        plt.savefig(f"{pathway_output_dir}/{pathway_clean}_heatmap.png", 
                                  bbox_inches="tight")
                        plt.savefig(f"{pathway_output_dir}/{pathway_clean}_heatmap.pdf", 
                                  bbox_inches="tight")
                        plt.close()
                    except Exception as e:
                        print(f"警告: 生成 {pathway} 热图失败: {str(e)}")
                        
                except Exception as e:
                    print(f"警告: 处理通路 {pathway} 失败: {str(e)}")
                    continue
            
            results[sample] = {
                'status': 'success',
                'significant_pathways': len(sig_pathways),
                'sender_receiver_pairs': len(sender_receiver)
            }
            print(f"✅ {sample} 处理完成")
            
        except Exception as e:
            print(f"❌ {sample} 处理失败: {str(e)}")
            results[sample] = {'status': 'failed', 'error': str(e)}
            continue
    
    return results

# 使用示例
# results = process_cellchat_analysis(adata, color_dict, base_path='.', celltype_col="spatial_sub", group_col="group")