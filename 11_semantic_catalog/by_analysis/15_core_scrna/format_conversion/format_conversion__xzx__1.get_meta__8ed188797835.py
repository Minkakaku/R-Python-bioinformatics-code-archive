import os
import matplotlib.pyplot as plt
import scanpy as sc
import sys
sys.path.append('/mnt/NexData_file/share/scRNA-tools/Proption/pyscript')
from draw_group import draw_group
import pandas as pd



# 需要传入参数
adata_file = '/mnt/NexData_file/share/scRNA-tools/Proption/test/anno.h5ad'
outdir = '/mnt/NexData_file/share/scRNA-tools/Proption/test/output'
method = "cell_type"
group = 'group'


WorkDir = os.getcwd()
adata = sc.read_h5ad(adata_file)

# 创建输出目录
os.makedirs(outdir,exist_ok=True)
# 依据method选择绘制
draw_group(adata,out_dir = outdir,method = method)

def get_color(adata,method):
    keys = adata.obs[method].cat.categories.to_list()
    values = list(adata.uns[method+'_colors'])
    color_dict = dict(zip(keys,values))
    df = pd.DataFrame.from_dict(color_dict,orient='index',columns=[method])
    df.to_csv(os.path.join(outdir,method+'_color.xls'),sep="\t",index=True,index_label=method)

get_color(adata,method=method)
get_color(adata,method=group)
adata.obs.to_csv(os.path.join(outdir,"meta.xls"),sep ="\t")

# 2.进行分析
shell_script = os.path.join(WorkDir,"2.pip.sh")  # 替换为实际的shell脚本路径
os.system(f"""
    sh \
    {shell_script} \
    {os.path.abspath(os.path.join(outdir, 'meta.xls'))} \
    {method} \
    {group} \
    {os.path.abspath(outdir)} \
    {os.path.abspath(os.path.join(outdir, method+'_color.xls'))} \
    {os.path.abspath(os.path.join(outdir, group+'_color.xls'))}
""")

# 生成结果报告
report_rmd = os.path.join(WorkDir,'/report/report.rmd')
os.system(f"""
    cp \
    {report_rmd} \
    {os.path.abspath(os.path.join(outdir))} \
""")
shell_script = os.path.join(WorkDir,"3.report.sh")
os.system(f"""
    sh \
    {shell_script} \
    {os.path.abspath(os.path.join(outdir))}
""")
# 删除中间文件
os.system(f"""
          rm \
          {os.path.abspath(os.path.join(outdir,'report.knit.md'))}
          rm \
          {os.path.abspath(os.path.join(outdir,'report.rmd'))}
          """)

# 去除¶
file_path = os.path.join(outdir,'report.html')

# 读取HTML文件内容
with open(file_path, 'r', encoding='utf-8') as file:
    html_content = file.read()

html_content = html_content.replace('¶', '')
# 如果是HTML实体 &para; 或 &182;
html_content = html_content.replace('&para;', '').replace('&182;', '')

# 写回文件
with open(file_path, 'w', encoding='utf-8') as file:
    file.write(html_content)
