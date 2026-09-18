import os
import anndata
import argparse
import warnings
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(prog='pipeline_cellphonedb.py', description="cellphonedb flow")
parser.add_argument('-i', action='store', dest='h5ad', required=True, help='adata.h5ad')
parser.add_argument('-o', action='store', dest='out_dir', help='Result output directory. default:cellphonedb_result')
parse = parser.parse_args()

if not parse.out_dir:
    parse.out_dir = 'cellphonedb_result'
if os.path.exists(parse.out_dir):
    os.system(f'rm -r {parse.out_dir}')
os.makedirs(parse.out_dir, exist_ok=True)

cpdb_file_path = '/home/hf/project/shen/mfcellphonedb/v5.0.0/cellphonedb.zip'
counts_file_path = parse.h5ad
out_path = parse.out_dir

adata = anndata.read_h5ad(counts_file_path)

assert adata.obs.index.name in ['barcode', 'cell_id', None], "adata obs index needs cell barcode"
meta = adata.obs[['cell_type']].copy()

meta.columns = ['cell_type']
meta.to_csv(os.path.join(out_path, 'meta.tsv'), index=True, sep='\t')
meta_file_path = os.path.join(out_path, 'meta.tsv')


try:

    deconvoluted, deconvoluted_percents, means, pvalues, significant_means, relevant_interactions = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,                 #    裺CellPhoneDB   ݿ ѹ    
        meta_file_path = meta_file_path,                 #    裺cell barcode  cell type  ӳ   ļ 
        counts_file_path = counts_file_path,             #    裺  һ   ı      h5ad  ʽ  
        counts_data = 'hgnc_symbol',                     #        Ļ   ע     ͣ     adata.var  ƥ 䣩
        iterations = 1000,                               #    ϴ ƴ   
        threshold = 0.01,                                #         ϸ        ֵ
        threads = 4,                                       #  滻threadsΪcores  v5.0.0 ²     
        debug_seed = 42,                                 #            
        result_precision = 3,                            #      Ծ ֵ  С  λ  
        pvalue = 0.05,                                   #       Pֵ  ֵ
        subsampling = False,                             #      Ӳ   
        subsampling_log = False,                         #  Ӳ     logת  
        subsampling_num_pc = 100,                        #  Ӳ        ɷ   
        subsampling_num_cells = 15000,                   #  Ӳ     ϸ    
        separator = '|',                                 # ϸ   Եķָ   
        debug = False,                                   #    õ   ģʽ
        output_path = out_path,                          #       ·  
        output_suffix = ""                               #     ļ   ׺
    )
except ValueError as e:

    all_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,
        meta_file_path = meta_file_path,
        counts_file_path = counts_file_path,
        counts_data = 'hgnc_symbol',
        iterations = 1000,
        threshold = 0.01,
        threads = 4,
        debug_seed = 42,
        result_precision = 3,
        pvalue = 0.05,
        subsampling = False,
        subsampling_log = False,
        subsampling_num_pc = 100,
        subsampling_num_cells = 15000,
        separator = '|',
        debug = False,
        output_path = out_path,
        output_suffix = ""
    )
    print(f"        ֵ      {len(all_results)}")

    deconvoluted = all_results[0]
    means = all_results[2]
    pvalues = all_results[3]
    significant_means = all_results[4]


kpy.plot_cpdb_heatmap(
    adata = adata,
    pvals = pvalues,
    celltype_key = "cell_type",
    figsize = (15,15),
    title = "Sum of significant interactions",
    symmetrical = False,
)
plt.savefig(os.path.join(out_path, 'heatmap.png'), facecolor='white', bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_path, 'heatmap.pdf'), facecolor='white')
plt.close()


cell_type_list = adata.obs['cell_type'].drop_duplicates().tolist()


from cellphonedb.utils import search_utils

for cell_type in cell_type_list:
    try:

        search_result = search_utils.search_analysis_results(
            query_cell_types_1 = [cell_type],
            query_cell_types_2 = cell_type_list,
            significant_means = significant_means,
            deconvoluted = deconvoluted,
            separator = '|',
            long_format = False
        )

        search_result['Total'] = search_result.iloc[:, 5:].sum(axis=1)
        df_sorted = search_result.sort_values(by='Total', ascending=False)
        top_10_rows = df_sorted.head(10)

        gene_columns = ['gene_a', 'gene_b']
        top_10_genes = top_10_rows[gene_columns].values.flatten().tolist()
        top_10_genes_unique = list(set([g for g in top_10_genes if pd.notna(g)]))
        

        p = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = cell_type,
            cell_type2 = ".",
            means = means,
            pvals = pvalues,
            genes = top_10_genes_unique,
            celltype_key = "cell_type",
            figsize = (20,10),
            title = f"{cell_type} Top10 Interactions",
            keep_significant_only = True,
            max_size = 6,
            highlight_size = 0.75,
            standard_scale = True
        )
        p.save(os.path.join(out_path, f'{cell_type}_Top10_bubble.png'), facecolor='white', bbox_inches='tight', dpi=300)
        p.save(os.path.join(out_path, f'{cell_type}_Top10_bubble.pdf'), facecolor='white')
        plt.close()
        
    except Exception as e:
        print(f"    {cell_type}ʱ      {e}")

        p = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = cell_type,
            cell_type2 = ".",
            means = means,
            pvals = pvalues,
            celltype_key = "cell_type",
            figsize = (20,10),
            title = f"{cell_type} Interactions",
            keep_significant_only = True,
            max_size = 6,
            highlight_size = 0.75,
            standard_scale = True
        )
        p.save(os.path.join(out_path, f'{cell_type}_bubble.png'), facecolor='white', bbox_inches='tight', dpi=300)
        p.save(os.path.join(out_path, f'{cell_type}_bubble.pdf'), facecolor='white')
        plt.close()


os.remove(meta_file_path)