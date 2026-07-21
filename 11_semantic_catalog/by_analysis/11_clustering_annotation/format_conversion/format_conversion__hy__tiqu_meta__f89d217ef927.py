import scanpy as sc
import pandas as pd
import os
import glob

def batch_extract_features(h5ad_pattern="*.h5ad", output_dir="extracted_features"):
    """
    批量处理多个h5ad文件
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有h5ad文件
    h5ad_files = glob.glob(h5ad_pattern)
    
    if not h5ad_files:
        print(f"没有找到匹配 {h5ad_pattern} 的文件")
        return
    
    print(f"找到 {len(h5ad_files)} 个h5ad文件")
    
    all_features = []
    
    for h5ad_file in h5ad_files:
        print(f"\n处理: {h5ad_file}")
        
        try:
            # 读取数据
            adata = sc.read_h5ad(h5ad_file)
            
            # 目标列
            target_columns = ['cell_id', 'sample', 'sample_new', 'state', 'celltype', 'CNs_10', 'CNs_15']
            available_columns = [col for col in target_columns if col in adata.obs.columns]
            
            if not available_columns:
                print(f"  跳过: 没有目标列")
                continue
            
            # 提取数据
            feature_df = adata.obs[available_columns].copy()
            
            # 添加文件名列
            filename = os.path.basename(h5ad_file)
            feature_df['source_file'] = filename
            
            all_features.append(feature_df)
            
            # 单独保存每个文件
            output_name = filename.replace('.h5ad', '_features.csv')
            output_path = os.path.join(output_dir, output_name)
            feature_df.to_csv(output_path, index=False, encoding='utf-8-sig')
            
            print(f"  已保存: {output_path} ({len(feature_df)} 行)")
            
        except Exception as e:
            print(f"  处理失败: {e}")
    
    # 合并所有数据
    if all_features:
        combined_df = pd.concat(all_features, ignore_index=True)
        
        # 保存合并文件
        combined_path = os.path.join(output_dir, "all_features_combined.csv")
        combined_df.to_csv(combined_path, index=False, encoding='utf-8-sig')
        
        print(f"\n合并完成!")
        print(f"总行数: {len(combined_df)}")
        print(f"已保存到: {combined_path}")
        
        return combined_df
    else:
        print("\n没有成功提取任何数据")
        return None

# 运行批量处理
if __name__ == "__main__":
    combined_features = batch_extract_features(
        h5ad_pattern="All_filter_state_samp_new_with_CNs10-15.h5ad",
        output_dir="extracted_features"
    )
