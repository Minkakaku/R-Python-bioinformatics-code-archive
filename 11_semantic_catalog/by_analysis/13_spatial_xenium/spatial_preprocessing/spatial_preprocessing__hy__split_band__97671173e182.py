import scanpy as sc
import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def split_h5ad_preserve_raw_structure(input_path, site_band_col='site_band', output_dir='./split_data'):
    """
    分割h5ad文件，正确保留原始raw矩阵结构
    
    原始数据：
    - adata.X: 所有基因的表达矩阵
    - adata.raw: 5001个高变基因的原始表达矩阵
    
    分割后每个文件：
    - 主数据X: 所有基因
    - raw矩阵: 5001个高变基因的原始表达矩阵
    
    Args:
        input_path: 输入h5ad文件路径
        site_band_col: 包含site_band信息的列名
        output_dir: 输出目录
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取原始数据
    print(f"读取数据: {input_path}")
    adata = sc.read_h5ad(input_path)
    
    # 检查site_band列
    if site_band_col not in adata.obs.columns:
        print(f"错误: 未找到列 '{site_band_col}'")
        print(f"可用列: {list(adata.obs.columns)}")
        return
    
    # 检查数据结构
    print(f"\n原始数据信息:")
    print(f"  主数据形状 (adata.X): {adata.shape}")
    print(f"  主数据基因数: {adata.n_vars}")
    
    if adata.raw is not None:
        if hasattr(adata.raw, 'to_adata'):
            raw_adata = adata.raw.to_adata()
            print(f"  raw矩阵形状: {raw_adata.shape}")
            print(f"  raw矩阵基因数: {raw_adata.n_vars}")
            
            # 确认是5001个基因
            if raw_adata.n_vars == 5001:
                print(f"  ✓ raw矩阵包含5001个高变基因")
            else:
                print(f"  注意: raw矩阵包含{raw_adata.n_vars}个基因，不是5001个")
        else:
            print(f"  raw矩阵类型: {type(adata.raw)}")
    else:
        print("  警告: 数据不包含raw矩阵")
        return
    
    # 获取原始raw的AnnData对象
    raw_adata = adata.raw.to_adata()
    
    # 验证细胞索引
    print(f"\n验证细胞索引:")
    print(f"  主数据细胞数: {len(adata.obs_names)}")
    print(f"  raw数据细胞数: {len(raw_adata.obs_names)}")
    
    # 检查是否有不匹配的细胞
    cells_in_main_not_raw = [cell for cell in adata.obs_names if cell not in raw_adata.obs_names]
    cells_in_raw_not_main = [cell for cell in raw_adata.obs_names if cell not in adata.obs_names]
    
    if cells_in_main_not_raw:
        print(f"  警告: {len(cells_in_main_not_raw)}个细胞在主数据中但不在raw中")
    if cells_in_raw_not_main:
        print(f"  警告: {len(cells_in_raw_not_main)}个细胞在raw中但不在主数据中")
    
    # 获取共同的细胞
    common_cells = set(adata.obs_names) & set(raw_adata.obs_names)
    print(f"  共同细胞数: {len(common_cells)}")
    
    if len(common_cells) == 0:
        print("错误: 主数据和raw数据没有共同细胞")
        return
    
    # 处理每个site_band分组
    for band in ['in', 'on', 'out']:
        if band not in adata.obs[site_band_col].values:
            print(f"\n跳过: 未找到 '{band}' 组")
            continue
        
        # 获取当前band的细胞索引
        band_mask = adata.obs[site_band_col] == band
        band_cells = adata.obs_names[band_mask]
        
        if len(band_cells) == 0:
            print(f"\n跳过: '{band}' 组没有细胞")
            continue
        
        print(f"\n处理 '{band}' 组: {len(band_cells)} 个细胞")
        
        # 1. 创建新的AnnData对象（主数据，所有基因）
        new_adata = adata[band_cells, :].copy()
        print(f"  主数据形状: {new_adata.shape}")
        
        # 2. 从原始raw中提取相同的细胞（5001个高变基因）
        # 检查哪些细胞在raw中
        available_cells = [cell for cell in band_cells if cell in raw_adata.obs_names]
        
        if len(available_cells) != len(band_cells):
            print(f"  警告: raw中只有 {len(available_cells)}/{len(band_cells)} 个细胞")
            if len(available_cells) == 0:
                print(f"  错误: 在raw中找不到任何细胞，跳过此组")
                continue
            
            # 如果只有部分细胞在raw中，也相应调整主数据
            new_adata = new_adata[available_cells, :].copy()
            print(f"  调整后主数据形状: {new_adata.shape}")
        
        # 从raw中提取相应的细胞
        new_raw_adata = raw_adata[available_cells, :].copy()
        print(f"  raw数据形状: {new_raw_adata.shape}")
        
        # 3. 设置raw属性
        new_adata.raw = new_raw_adata
        
        # 4. 保存文件
        output_path = os.path.join(output_dir, f'{band}.h5ad')
        new_adata.write(output_path)
        print(f"  保存到: {output_path}")
        
        # 5. 验证保存的文件
        print(f"  验证保存的文件...")
        try:
            loaded = sc.read_h5ad(output_path)
            print(f"    加载成功: 主数据形状={loaded.shape}")
            
            if loaded.raw is not None:
                loaded_raw = loaded.raw.to_adata() if hasattr(loaded.raw, 'to_adata') else loaded.raw
                print(f"    raw矩阵形状={loaded_raw.shape}")
                
                # 验证基因数
                if loaded_raw.shape[1] == 5001:
                    print(f"    ✓ raw矩阵包含5001个基因")
                else:
                    print(f"    raw矩阵包含{loaded_raw.shape[1]}个基因")
                
                # 验证是否是相同的矩阵
                original_raw_shape = new_raw_adata.shape
                loaded_raw_shape = loaded_raw.shape
                
                if original_raw_shape == loaded_raw_shape:
                    print(f"    ✓ raw矩阵形状匹配")
                else:
                    print(f"    ✗ raw矩阵形状不匹配")
                    
            else:
                print(f"    ✗ 加载的文件没有raw矩阵")
                
        except Exception as e:
            print(f"    验证失败: {e}")
    
    print(f"\n{'='*60}")
    print(f"分割完成! 输出目录: {output_dir}")
    
    # 打印汇总统计
    print(f"\n汇总:")
    for band in ['in', 'on', 'out']:
        file_path = os.path.join(output_dir, f'{band}.h5ad')
        if os.path.exists(file_path):
            try:
                data = sc.read_h5ad(file_path)
                raw_shape = data.raw.shape if data.raw is not None else None
                print(f"  {band}: 主数据={data.shape}, raw矩阵={raw_shape}")
            except Exception as e:
                print(f"  {band}: 加载失败 - {e}")
        else:
            print(f"  {band}: 文件不存在")

def verify_split_results(output_dir='./split_data'):
    """
    验证分割结果的详细信息
    """
    print(f"\n{'='*60}")
    print(f"详细验证分割结果")
    print(f"{'='*60}")
    
    for band in ['in', 'on', 'out']:
        file_path = os.path.join(output_dir, f'{band}.h5ad')
        
        if not os.path.exists(file_path):
            print(f"\n{band}: 文件不存在")
            continue
            
        print(f"\n{band}:")
        
        try:
            # 加载数据
            data = sc.read_h5ad(file_path)
            
            # 基本信息
            print(f"  文件名: {file_path}")
            print(f"  细胞数: {data.n_obs}")
            print(f"  基因数: {data.n_vars}")
            
            # 检查矩阵类型
            print(f"  X矩阵类型: {type(data.X)}")
            
            # 检查raw矩阵
            if data.raw is not None:
                raw_data = data.raw.to_adata() if hasattr(data.raw, 'to_adata') else data.raw
                print(f"  raw矩阵形状: {raw_data.shape}")
                print(f"  raw矩阵类型: {type(raw_data.X)}")
                
                # 比较基因数
                if data.n_vars != raw_data.n_vars:
                    print(f"  ✓ 主数据和raw数据基因数不同: {data.n_vars} vs {raw_data.n_vars}")
                else:
                    print(f"  ⚠ 主数据和raw数据基因数相同: 都是{data.n_vars}")
                
                # 检查是否是高变基因
                if raw_data.n_vars == 5001:
                    print(f"  ✓ raw矩阵包含5001个高变基因")
                elif 4800 <= raw_data.n_vars <= 5200:
                    print(f"  ✓ raw矩阵包含约{raw_data.n_vars}个高变基因")
                else:
                    print(f"  ⚠ raw矩阵包含{raw_data.n_vars}个基因")
            else:
                print(f"  ✗ 没有raw矩阵")
                
        except Exception as e:
            print(f"  验证失败: {e}")

# 主函数
if __name__ == "__main__":
    # 设置输入文件路径
    input_file = "adata_with_band_by_sample.h5ad"  # 修改为您的文件路径
    
    # 运行分割函数
    print("="*60)
    print("开始分割h5ad数据")
    print("目标: 按site_band分割，保留原始raw矩阵结构")
    print("="*60)
    
    split_h5ad_preserve_raw_structure(
        input_path=input_file,
        site_band_col='site_band',
        output_dir='./site_band_split'
    )
    
    # 运行验证
    verify_split_results(output_dir='./site_band_split')
    
    # 额外的验证：检查矩阵内容
    print(f"\n{'='*60}")
    print(f"额外的矩阵内容验证")
    print(f"{'='*60}")
    
    for band in ['in', 'on', 'out']:
        file_path = f'./site_band_split/{band}.h5ad'
        if os.path.exists(file_path):
            try:
                data = sc.read_h5ad(file_path)
                print(f"\n{band}:")
                
                # 检查X矩阵的稀疏性
                if hasattr(data.X, 'nnz'):
                    print(f"  X矩阵非零元素: {data.X.nnz}")
                    print(f"  X矩阵稀疏度: {data.X.nnz/(data.n_obs*data.n_vars)*100:.2f}%")
                else:
                    print(f"  X矩阵是密集矩阵")
                
                # 检查raw矩阵
                if data.raw is not None:
                    raw_data = data.raw.to_adata()
                    if hasattr(raw_data.X, 'nnz'):
                        print(f"  raw矩阵非零元素: {raw_data.X.nnz}")
                        print(f"  raw矩阵稀疏度: {raw_data.X.nnz/(raw_data.n_obs*raw_data.n_vars)*100:.2f}%")
                    else:
                        print(f"  raw矩阵是密集矩阵")
                        
            except Exception as e:
                print(f"\n{band}: 验证失败 - {e}")
