import anndata as ad
import os
from tqdm import tqdm

# 1. 设置输入和输出目录
input_dir = "/home/hy/Project/FWSP20251095/Redo/subcluster/sub"  # 替换为实际目录
output_dir = "/home/hy/Project/FWSP20251095/Redo/subcluster/sub/new"  # 替换为实际目录
os.makedirs(output_dir, exist_ok=True)  # 确保输出目录存在

# 2. 遍历目录中的所有h5ad文件
h5ad_files = [f for f in os.listdir(input_dir) if f.endswith('.h5ad')]
print(f"找到 {len(h5ad_files)} 个h5ad文件")

# 3. 批量处理函数
def add_celltype_sub(adata):
    """核心操作：复制celltype列到新列celltype_sub"""
    if 'celltype' not in adata.obs.columns:
        raise KeyError("当前文件缺失 'celltype' 列，请检查数据")
    adata.obs['celltype_sub'] = adata.obs['celltype']  # 直接复制列数据
    return adata

# 4. 处理文件并保存（带进度条和错误处理）
processed_files = []
failed_files = []

for filename in tqdm(h5ad_files, desc="处理文件中"):
    try:
        # 读取文件
        adata = ad.read_h5ad(os.path.join(input_dir, filename))
        
        # 添加新列
        adata = add_celltype_sub(adata)
        
        # 保存到新路径（避免覆盖原始文件）
        output_path = os.path.join(output_dir, f"modified_{filename}")
        adata.write_h5ad(output_path)
        processed_files.append(filename)
        
    except Exception as e:
        failed_files.append((filename, str(e)))
        print(f"\n处理失败: {filename} - 错误: {e}")

# 5. 输出结果报告
print("\n处理完成！")
print(f"成功处理文件数: {len(processed_files)}")
print(f"失败文件数: {len(failed_files)}")
if failed_files:
    print("\n失败详情:")
    for file, error in failed_files:
        print(f"  - {file}: {error}")
