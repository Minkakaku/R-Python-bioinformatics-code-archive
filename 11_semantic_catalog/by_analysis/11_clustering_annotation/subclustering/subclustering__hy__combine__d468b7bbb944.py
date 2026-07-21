import h5py
import anndata as ad
import os
import numpy as np

def merge_h5ad_files(input_files, output_file, merge_mode='union'):
    """
    合并多个h5ad文件到一个新文件
    :param input_files: 输入h5ad文件路径列表
    :param output_file: 输出文件路径
    :param merge_mode: 合并模式 ('union'保留所有特征, 'intersection'仅保留共同特征)
    """
    # 初始化输出文件
    with h5py.File(output_file, 'w') as h5_out:
        # 遍历每个输入文件
        for file_idx, file_path in enumerate(input_files): 
            with h5py.File(file_path, 'r') as h5_in:
                # 递归复制组和数据集
                def copy_datasets(name, obj):
                    # 重命名根组避免冲突（例如 /obs 变为 /file0_obs）
                    new_name = f"/file{file_idx}_{name.lstrip('/')}" if name != '/' else name
                    
                    if isinstance(obj, h5py.Dataset):
                        # 合并关键数据结构（obs/var/X）
                        if name.endswith('/obs'):
                            merge_obs(h5_out, obj, file_idx)
                        elif name.endswith('/var'):
                            merge_var(h5_out, obj, file_idx, merge_mode)
                        elif name.endswith('/X'):
                            merge_X(h5_out, obj, file_idx)
                        else:
                            # 直接复制非关键数据集
                            h5_in.copy(obj, h5_out, name=new_name)
                    elif isinstance(obj, h5py.Group):
                        # 创建组并递归处理内容
                        if new_name not in h5_out:
                            h5_out.create_group(new_name)
                        for sub_name in obj:
                            copy_datasets(f"{name}/{sub_name}", obj[sub_name])
                
                # 从根目录开始遍历
                copy_datasets('/', h5_in['/'])
    
    print(f"合并完成！输出文件: {output_file}")

def merge_obs(h5_out, obs_dataset, file_idx):
    """合并obs（细胞注释）数据"""
    obs_data = ad._io.h5ad.read_elem(obs_dataset)  # 解析h5ad的obs结构
    if 'file_idx' not in obs_data.columns:
        obs_data['batch'] = file_idx  # 添加批次标识
    # 将obs追加到输出文件
    if '/obs' not in h5_out:
        ad._io.h5ad.write_elem(h5_out, 'obs', obs_data)
    else:
        existing_obs = ad._io.h5ad.read_elem(h5_out['/obs'])
        merged_obs = existing_obs.concat([obs_data])
        ad._io.h5ad.write_elem(h5_out, 'obs', merged_obs)

def merge_var(h5_out, var_dataset, file_idx, merge_mode):
    """合并var（基因特征）数据"""
    var_data = ad._io.h5ad.read_elem(var_dataset)
    if '/var' not in h5_out:
        ad._io.h5ad.write_elem(h5_out, 'var', var_data)
        h5_out.attrs['feature_names'] = var_data.index.values  # 保存特征名
    else:
        existing_var = ad._io.h5ad.read_elem(h5_out['/var'])
        if merge_mode == 'union':
            # 取并集（自动去重）
            merged_var = existing_var.index.union(var_data.index)
        else:  # intersection
            merged_var = existing_var.index.intersection(var_data.index)
        # 更新特征矩阵索引
        h5_out.attrs['feature_names'] = list(merged_var)

def merge_X(h5_out, X_dataset, file_idx):
    """合并X（表达矩阵）数据"""
    X_data = ad._io.h5ad.read_elem(X_dataset)
    if '/X' not in h5_out:
        ad._io.h5ad.write_elem(h5_out, 'X', X_data)
    else:
        existing_X = ad._io.h5ad.read_elem(h5_out['/X'])
        # 沿行（细胞轴）拼接
        merged_X = np.vstack([existing_X, X_data])
        ad._io.h5ad.write_elem(h5_out, 'X', merged_X)

if __name__ == "__main__":
    # ===== 配置参数 =====
    input_files = [
    "Adipocytes.h5ad",
    "Endothelials.h5ad",
    "Epithelials.h5ad",
    "Fibroblasts.h5ad",
    "Kerationocytes_subcluster.h5ad",
    "Macro_subcluster.h5ad",
    "Mast_subcluster.h5ad",
    "Mixed_inflammatory_cells.h5ad",
    "Myeloid_subcluster.h5ad",
    "Nerve_cells.h5ad",
    "Pericytes.h5ad",
    "SMC.h5ad",
    "Stratum_basal_cells_subcluster.h5ad",
    "T_cells_subcluster.h5ad"
    ]
    output_file = "merged_data.h5ad"
    merge_mode = 'union'  # 可选'union'或'intersection'
    # ====================
    
    merge_h5ad_files(input_files, output_file, merge_mode)
