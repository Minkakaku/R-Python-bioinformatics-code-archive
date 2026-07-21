import scanpy as sc
import pandas as pd

# 读取h5ad文件
adata = sc.read_h5ad('modify.h5ad')  # 替换为您的文件路径

# 定义映射关系字典
mapping_dict = {"Stress-responsive_glandular_Epithelium":"Tumor","ECM_CAF":"ECM_CAF","Proliferative_Epithelium":"Tumor","Luminal_Epithelium":"Luminal_Epithelium","Endo-Epi":"Endo-Epi","FGFR2+_Basal_Epithelium":"Tumor","MCAM+_CAF":"MCAM+_CAF","Basal_Epithelium":"Basal_Epithelium","CD8+_Tem":"CD8+_Tem","EMT_Epithelium":"Tumor","Epi-CAF":"Epi-CAF","myCAF":"myCAF","Arterial_endothelium":"Arterial_endothelium","cDC":"cDC","M2":"M2","MKI67+_T":"MKI67+_T","Neutrophil":"Neutrophil","B_follicular":"B_follicular","CD4+_Tex":"CD4+_Tex","EMT_myCAF":"EMT_myCAF","<200b>LYVE1+_TAM":"LYVE1+_TAM","CD8+_Tcm":"CD8+_Tcm","NK":"NK","Renal_Epithelium":"Renal_Epithelium","IFN-T":"IFN-T","MKI67+_Plasma":"MKI67+_Plasma","Plasma":"Plasma","B_germinal_center":"B_germinal_center","pDC":"pDC","EMT_M2":"EMT_M2","CD8+_Tex":"CD8+_Tex","Treg":"Treg","CD4+_Tcm":"CD4+_Tcm","TAM":"TAM","iCAF":"iCAF","Lymphatic_endothelium":"Lymphatic_endothelium","B_memory":"B_memory","Tumor_Vascular_endothelium":"Tumor_Vascular_endothelium","Vascular_endothelium":"Vascular_endothelium","LYVE1+_TAM":"LYVE1+_TAM"}

# 检查annotation列是否存在
if 'annotation' not in adata.obs.columns:
    print("错误：adata.obs中没有找到'annotation'列")
else:
    # 创建anno_new列，根据映射关系赋值
    adata.obs['anno_new'] = adata.obs['annotation'].map(mapping_dict)
    
    # 检查是否有未映射的值
    unmapped_values = adata.obs.loc[adata.obs['anno_new'].isna(), 'annotation'].unique()
    if len(unmapped_values) > 0:
        print(f"警告：以下annotation值没有映射关系，将保持为NaN：{unmapped_values}")
    
    # 查看新列的值分布
    print("anno_new列的值分布：")
    print(adata.obs['anno_new'].value_counts())
    
    # 可选：保存修改后的数据
    adata.write('modify_tumor.h5ad')
    print("修改完成，已保存为新的h5ad文件")

# 可选：导出到CSV查看
adata.obs[['annotation', 'anno_new']].to_csv('annotation_mapping_result.csv')
print("映射结果已保存为annotation_mapping_result.csv")
