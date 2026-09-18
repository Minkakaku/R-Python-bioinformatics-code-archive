import diopy
# 保存图片
import os
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import font_manager as fm
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size=None):
    plt.tight_layout()
    if fig_size is not None:
        width = fig_size.get('w', 5)  # 默认宽度5英寸
        height = fig_size.get('h', 5)  # 默认高度5英寸
        plt.gcf().set_size_inches(width, height)
    
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

# 加载亚群细分的函数
import sys
sys.path.append('/mnt/NexData_file/home/yrc/scrna-code/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion
import pandas as pd