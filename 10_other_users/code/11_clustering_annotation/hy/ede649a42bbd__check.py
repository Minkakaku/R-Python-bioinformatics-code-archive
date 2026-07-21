import pandas as pd
import scanpy as sc

def diagnose_stage_issue(adata_path, csv_files_dir):
    """
    诊断h5ad文件与CSV文件匹配问题
    """
    print("=== 开始诊断h5ad文件与CSV文件的匹配问题 ===\n")
    
    # 1. 读取h5ad文件
    print("1. 读取h5ad文件...")
    try:
        adata = sc.read_h5ad(adata_path)
        print(f"   ✅ 成功读取h5ad文件")
        print(f"   数据形状: {adata.shape}")
        print(f"   obs列名: {list(adata.obs.columns)}")
    except Exception as e:
        print(f"   ❌ 读取h5ad文件失败: {e}")
        return
    
    # 2. 检查必要的列是否存在
    print("\n2. 检查必要的列...")
    necessary_cols = ['sample', 'cell_id']
    missing_cols = [col for col in necessary_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"   ❌ 缺失的列: {missing_cols}")
        print(f"   现有的列: {list(adata.obs.columns)}")
        return
    else:
        print(f"   ✅ 必要的列都存在: {necessary_cols}")
    
    # 3. 检查sample列的取值
    print("\n3. 检查sample列的取值...")
    unique_samples = adata.obs['sample'].unique()
    print(f"   所有样本: {list(unique_samples)}")
    
    # 需要从CSV处理的样本
    csv_samples = ['009D-1', 'SR202500003B', 'SR202500005F-1']
    print(f"   需要从CSV处理的样本: {csv_samples}")
    
    # 4. 检查cell_id列的基本信息
    print("\n4. 检查cell_id列...")
    print(f"   cell_id列数据类型: {adata.obs['cell_id'].dtype}")
    print(f"   cell_id示例值: {adata.obs['cell_id'].head().tolist()}")
    print(f"   cell_id唯一值数量: {adata.obs['cell_id'].nunique()}")
    print(f"   cell_id是否有空值: {adata.obs['cell_id'].isnull().sum()}")
    
    # 5. 逐个检查CSV文件
    print("\n5. 检查CSV文件...")
    for sample in csv_samples:
        csv_file_path = f"{csv_files_dir}/{sample}.csv"
        print(f"\n   检查样本 {sample}:")
        
        try:
            # 尝试读取CSV文件
            stage_df = pd.read_csv(csv_file_path, header=None)
            print(f"      ✅ CSV文件存在，形状: {stage_df.shape}")
            
            # 检查列数
            if len(stage_df.columns) >= 2:
                print(f"      第一列示例: {stage_df[0].head().tolist()}")
                print(f"      第二列示例: {stage_df[1].head().tolist()}")
                print(f"      CSV中细胞ID数据类型: {stage_df[0].dtype}")
                print(f"      CSV中stage数据类型: {stage_df[1].dtype}")
            else:
                print(f"      ❌ CSV文件列数不足: {len(stage_df.columns)}列，需要至少2列")
                
        except FileNotFoundError:
            print(f"      ❌ CSV文件不存在: {csv_file_path}")
        except Exception as e:
            print(f"      ❌ 读取CSV文件出错: {e}")
    
    # 6. 详细对比一个样本的匹配情况
    print("\n6. 详细匹配测试（以009D-1为例）...")
    test_sample = '009D-1'
    test_csv_path = f"{csv_files_dir}/{test_sample}.csv"
    
    try:
        # 读取测试CSV
        test_stage_df = pd.read_csv(test_csv_path, header=None, names=['cell_id', 'stage'])
        print(f"   CSV中的细胞ID数量: {len(test_stage_df)}")
        
        # 获取h5ad中该样本的细胞
        sample_cells = adata.obs[adata.obs['sample'] == test_sample]
        print(f"   h5ad中{test_sample}样本的细胞数量: {len(sample_cells)}")
        
        # 检查匹配情况
        h5ad_cell_ids = set(sample_cells['cell_id'].astype(str))
        csv_cell_ids = set(test_stage_df['cell_id'].astype(str))
        
        print(f"   CSV中有而h5ad中没有的细胞ID数量: {len(csv_cell_ids - h5ad_cell_ids)}")
        print(f"   h5ad中有而CSV中没有的细胞ID数量: {len(h5ad_cell_ids - csv_cell_ids)}")
        print(f"   共同存在的细胞ID数量: {len(h5ad_cell_ids & csv_cell_ids)}")
        
        if len(h5ad_cell_ids & csv_cell_ids) > 0:
            print("   ✅ 存在可匹配的细胞ID")
            print(f"   示例匹配的细胞ID: {list(h5ad_cell_ids & csv_cell_ids)[:3]}")
        else:
            print("   ❌ 没有共同的细胞ID，完全无法匹配")
            print(f"   h5ad中的示例ID: {list(h5ad_cell_ids)[:3]}")
            print(f"   CSV中的示例ID: {list(csv_cell_ids)[:3]}")
            
    except Exception as e:
        print(f"   测试匹配时出错: {e}")
    
    print("\n=== 诊断完成 ===")

# 使用示例
if __name__ == "__main__":
    # 请修改为您的实际路径
    h5ad_file = "adata_celltype.h5ad"  # 替换为您的h5ad文件路径
    csv_directory = "sample"  # 替换为CSV文件所在目录
    
    diagnose_stage_issue(h5ad_file, csv_directory)
