import scanpy as sc

# 读取 h5ad 文件
adata = sc.read_h5ad('adata_celltype_RCTD_removed_with_patient.h5ad')  # 替换为您的文件路径

# 检查 annotation 列是否存在
if 'annotation' not in adata.obs.columns:
    print("错误：未找到 annotation 列")
else:
    # 查看替换前的值分布（可选）
    print("替换前的值分布：")
    print(adata.obs['annotation'].value_counts())
    
    # 执行字符串替换
    adata.obs['annotation'] = adata.obs['annotation'].str.replace('<200b>LYVE1+_TAM', 'LYVE1+_TAM')
    
    # 查看替换后的值分布（可选）
    print("\n替换后的值分布：")
    print(adata.obs['annotation'].value_counts())
    
    # 保存修改后的数据
    adata.write('modify.h5ad')  # 保存为新文件，避免覆盖原文件
    print("修改完成，已保存为新的 h5ad 文件")
