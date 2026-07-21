import os
import pandas as pd
from typing import Dict, List
import scanpy as sc
# import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
# from Class_module import MarkerTable, H5ADHandler

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def get_return_value(marker_list):
    length = len(marker_list)
    if length > 4:
        return 20
    elif length > 0:
        return length*5

def preprocess_dict(marker_dict, adata) -> dict:
    """
    根据AnnData对象中的var_names与基因集-基因字典中各key的value取交集，
    重新赋给字典中的各key，并输出更新后的基因集-基因字典。
    :param adata: AnnData对象
    :return: 更新后的基因集-基因字典
    """
    var_names_set = set(adata.var_names)
    updated_marker_dict = {}
    for key, gene_list in marker_dict.items():
        # 将基因列表转换为集合并与 var_names_set 取交集
        intersected_genes = set(gene_list).intersection(var_names_set)
        # 如果交集中有基因，则更新到新的字典中
        if intersected_genes:
            updated_marker_dict[key] = list(intersected_genes)
    # 更新内部字典
    marker_dict = updated_marker_dict
    return updated_marker_dict

class MarkerTable:
    def __init__(self):
        pass
    
    @staticmethod
    def _read_single_gene_file(file_path: str) -> Dict[str, List[str]]:
        """
        处理单列gene的txt文件，生成基因集-基因字典。
        
        :param file_path: 输入文件路径
        :return: 基因集-基因字典，genelist为key的名字
        """
        with open(file_path, 'r') as file:
            genes = [line.strip() for line in file if line.strip()]
        gene_dict = {'genelist': genes}
        return gene_dict
    
    @staticmethod
    def _read_two_column_file(file_path: str, sep='\t') -> Dict[str, List[str]]:
        """
        处理两列的txt文件，第一列为key，第二列为逗号分割的gene，为value。
        :param file_path: 输入文件路径
        :param sep: 分隔符，默认为制表符
        :return: 基因集-基因字典
        """
        gene_dict = {}
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split(sep)
                if len(parts) == 2:
                    key, genes = parts
                    gene_list = [gene.strip() for gene in genes.split(',') if gene.strip()]
                    gene_dict[key] = gene_list
        return gene_dict
    
    @staticmethod
    def generate_gene_dict(input_file: str, species: str = None, tissue_types: List[str] = None, cellsub_types: List[str] = None) -> Dict[str, List[str]]:
        """
        根据输入文件类型生成基因集-基因字典。
        :param input_file: 输入文件路径
        :param species: 物种类型（mouse 或 human），仅适用于第一种格式
        :param tissue_types: 用户指定的一种或多种组织类型列表，仅适用于第一种格式
        :param cellsub_types: 用户指定的一种或多种细胞亚型列表，用于筛选 sub 列
        :return: 基因集-基因字典
        """
        try:
            # 尝试读取CSV/TSV文件，以确定是否为第一种格式
            df = pd.read_csv(input_file, sep='\t', nrows=0)
            if len(df.columns) > 2 and (species or tissue_types or cellsub_types):
                if tissue_types and cellsub_types:
                    raise ValueError("tissue_types 和 cellsub_types 不能同时使用。")
                else:
                    return MarkerTable._process_first_format(input_file, species, tissue_types, cellsub_types)
            # 如果不是第一种格式，则检查是单列还是两列文件
            with open(input_file, 'r') as file:
                first_line = file.readline().strip()
                if ',' in first_line:
                    return MarkerTable._read_two_column_file(input_file)
                else:
                    return MarkerTable._read_single_gene_file(input_file)
        except Exception as e:
            print(f"Error processing file: {e}")
            return {}
        
    @staticmethod
    def _process_first_format(input_file: str, species: str, tissue_types: List[str], cellsub_types: List[str]) -> Dict[str, List[str]]:
        """
        处理第一种输入表格格式，筛选指定的组织类型列（全称匹配）。
        :param input_file: 输入表格文件路径
        :param species: 物种类型（mouse 或 human）
        :param tissue_types: 用户指定的一种或多种组织类型列表
        :param cellsub_types: 用户指定的一种或多种细胞亚型列表，用于筛选 sub 列
        :return: 基因集-基因字典
        """
        df = pd.read_csv(input_file, sep='\t', header=0)
        if species == 'mouse':
            df.rename(columns={df.columns[0]: 'gene'}, inplace=True)
        elif species == 'human':
            df.rename(columns={df.columns[1]: 'gene'}, inplace=True)
        else:
            raise ValueError("Unsupported species type. Please choose 'mouse' or 'human'.")
        if tissue_types is not None:
            df['genelist'] = df['main'].astype(str) + '-' + df['sub'].astype(str)
            all_columns = df.columns.tolist()
            # 找到与提供的组织类型完全匹配的列
            tissue_columns = [col for col in all_columns[4:] if col in tissue_types]
            if not tissue_columns:
                print("警告：没有找到与提供的组织类型完全匹配的列。")
                return {}
            df_filtered = df[['gene', 'genelist'] + tissue_columns]
            df_filtered = df_filtered[df_filtered[tissue_columns].any(axis=1)]
        else:
            # 筛选 sub 列
            cellsub_types_list = cellsub_types.split(',')
            df_filtered = df[df['sub'].isin(cellsub_types_list)]
            if df_filtered.empty:
                print("警告：没有找到与提供的细胞亚型完全匹配的行。")
                return {}
            df_filtered['genelist'] = df_filtered['main'].astype(str) + '-' + df_filtered['sub'].astype(str)
        gene_dict = {}
        for _, row in df_filtered.iterrows():
            genelist = row['genelist']
            gene = row['gene']
            if genelist not in gene_dict:
                gene_dict[genelist] = []
            gene_dict[genelist].append(gene)
        return gene_dict


class H5ADHandler:
    def __init__(self, h5ad_path):
        # 读取h5ad文件
        self.adata = sc.read_h5ad(h5ad_path).raw.to_adata()
        
    def preprocess_data(self, marker_genes):
        """
        :param marker_genes: marker基因列表
        :return: 处理后的AnnData对象
        """
        adata = self.adata.copy()
        # 选择marker基因
        missing_genes = set(marker_genes) - set(adata.var_names)
        if missing_genes:
             print(f"警告：以下基因未找到: {missing_genes}")
        intersected_genes = list(adata.var_names.intersection(marker_genes))
        adata = adata[:, adata.var_names.isin(intersected_genes)]
        return adata
    
    def plot_umap(self, adata, marker_genes, output_path, prefix):
        """
        绘制UMAP图
        :param adata: AnnData对象
        :param output_path: 输出路径
        :param prefix: 文件前缀
        """
        sc.pl.umap(adata, color=marker_genes, show=False, cmap='coolwarm')
        save_fig(plt,file_name= f'{prefix}_umap',out_dir=output_path,fig_size={'w':get_return_value(marker_genes),'h':max(len(marker_genes)//4*5,5)})
        
    def plot_violin(self, adata, marker_genes, output_path, prefix):
        """
        绘制提琴图
        :param adata: AnnData对象
        :param marker_genes: marker基因列表
        :param output_path: 输出路径
        :param prefix: 文件前缀
        """
        sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=False, dendrogram=False)
        save_fig(plt,file_name= f'{prefix}_violin',out_dir=output_path,fig_size={'w':5,'h':5}) #这个地方需要算一下长宽
        
        
    def plot_bubble(self, adata, marker_genes, output_path, prefix):
        """
        绘制气泡图
        :param adata: AnnData对象
        :param marker_genes: marker基因列表
        :param output_path: 输出路径
        :param prefix: 文件前缀
        """
        sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', color_map='Purples', show=False)
        save_fig(plt,file_name= f'{prefix}_bubble',out_dir=output_path,fig_size={'w':5,'h':5}) #这个地方需要算一下长宽

def main(input_dir, marker_file, species, tissue_types, cellsub_types, output_dir):
    """
    主函数：处理输入参数，调用类方法生成图表
    :param input_dir: 输入目录，包含h5ad文件
    :param marker_file: marker表文件路径
    :param species: 物种（例如："小鼠" 或 "人"）
    :param tissue_types: 组织类型（例如："肠"、"胃" 等）
    :param output_dir: 输出目录
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    # 初始化MarkerTable类实例
    marker_dict = MarkerTable.generate_gene_dict(marker_file, species, tissue_types, cellsub_types)
    # 检查 input_dir 是目录还是单个文件
    input_path = Path(input_dir)
    if input_path.is_file() and input_path.suffix == '.h5ad':
        h5ad_files = input_path
    elif input_path.is_dir():
        h5ad_files = list(input_path.rglob('*.h5ad'))
        h5ad_files=h5ad_files[0]
    else:
        raise ValueError(f"Invalid input path: {input_dir}")
    # 输入目录中的h5ad文件选第一个
    print(f"Processing {h5ad_files}")
    h5ad_handler = H5ADHandler(h5ad_files)
    marker_sub = preprocess_dict(marker_dict,h5ad_handler.adata)

    for celltype, genes in marker_sub.items():
        processed_adata = h5ad_handler.preprocess_data(genes)
        # 绘制UMAP图
        h5ad_handler.plot_umap(processed_adata, genes, output_dir, celltype)
        # 绘制提琴图
        h5ad_handler.plot_violin(processed_adata, genes, output_dir, celltype)
        # 绘制气泡图
        h5ad_handler.plot_bubble(processed_adata, genes, output_dir, celltype)


if __name__ == "__main__":
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="Generate UMAP and Violin plots for single-cell data.")
    parser.add_argument('--input_dir', type=str, required=True, help='Input directory containing h5ad files.')
    parser.add_argument('--marker_file', type=str, required=True, help='Path to the marker table file (txt format).')
    parser.add_argument('--species', type=str, required=True, choices=['mouse', 'human'], help='Species: mouse or human.')
    parser.add_argument('--tissue_types', type=str, default=None, required=False, help='Tissue type (e.g., 肠, 胃).')
    parser.add_argument('--cellsub_types', type=str, default=None, required=False, help='cluster type (e.g., Tcell, NK).')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for saving plots.')

    args = parser.parse_args()

    # 调用主函数
    main(args.input_dir, args.marker_file, args.species, args.tissue_types, args.cellsub_types, args.output_dir)

# python generate_marker_plots.py --input_dir /path/to/h5ad/files --marker_file /path/to/marker_table.txt --species 小鼠 --tissue_type 肠 --output_dir /path/to/output/directory