import pandas as pd
import subprocess
import os
import sys
import tempfile

# 从命令行参数读取数据库路径
# Usage: python merge_multiple_annotations_step1.py <adduct_mod_path> <hmdb_metabolites_path> <inchikey_classyfire_path>
if len(sys.argv) < 4:
    print("Usage: python merge_multiple_annotations_step1.py <adduct_mod_tsv> <hmdb_metabolites_csv> <inchikey_classyfire_txt>")
    sys.exit(1)
adduct_mod_path = sys.argv[1]
hmdb_metabolites_path = sys.argv[2]
inchikey_classyfire_path = sys.argv[3]

# 读取5个 TSV 文件
df1 = pd.read_csv('result_CAMERA_annotated_inhouse_updated.tsv', sep='\t')
df2 = pd.read_csv('result_CAMERA_annotated_GNPS.tsv', sep='\t')
df3 = pd.read_csv('result_CAMERA_annotated_MoNA.tsv', sep='\t')
df4 = pd.read_csv('result_CAMERA_annotated_HMDB.tsv', sep='\t')
df5 = pd.read_csv('result_CAMERA_annotated_lipid.tsv', sep='\t')
df6 = pd.read_csv(adduct_mod_path, sep='\t')
df7 = pd.read_csv(hmdb_metabolites_path)

columns_to_add = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight', 'cas_registry_number', 'smiles', 'inchi', 'inchikey', 'kingdom', 'super_class', 'class', 'sub_class']
# 添加前缀以避免列名冲突
hmdb_subset = df7[columns_to_add].add_prefix('hmdb_')
hmdb_subset = hmdb_subset[hmdb_subset['hmdb_inchikey'].notna()]

# 转换为字典（raw 列为键，new 列为值）
mapping = dict(zip(df6['raw'], df6['new']))

# 动态找出 pcgroup 列所在的位置（包括它本身）
pcgroup_idx = df1.columns.get_loc('pcgroup') + 1  # 加1表示包含这一列
common = df1.iloc[:, :pcgroup_idx] 

# 分别提取每个文件中公共列以后的数据
extra1 = df1.iloc[:, pcgroup_idx:]
extra2 = df2.iloc[:, pcgroup_idx:]
extra3 = df3.iloc[:, pcgroup_idx:]
extra4 = df4.iloc[:, pcgroup_idx:]
extra5 = df5.iloc[:, pcgroup_idx:]

# 水平方向合并，公共列不变，后面的数据依次追加
merged_df = pd.concat([common, extra1, extra2, extra3, extra4, extra5], axis=1)

# 对公共部分中的CAMERA输出的isotopes和adduct数据内容进行清空（保留表头不变）
merged_df['isotopes'] = ''
merged_df['adduct'] = ''


# 定义一个函数来精简信息
def filter_columns(row):
    # 获取第39列的值并拆分成列表
    scores = str(row[pcgroup_idx+24]).split(';') if pd.notna(row[pcgroup_idx+24]) else []
    # 获取第38列的值并拆分成列表
    names = str(row[pcgroup_idx+23]).split(';') if pd.notna(row[pcgroup_idx+23]) else []
    # 获取第41列的值并拆分成列表
    inchikeys = str(row[pcgroup_idx+26]).split(';') if pd.notna(row[pcgroup_idx+26]) else []
    
    # 如果第39列有值
    if scores and any(scores):  # 确保scores列表不为空且至少有一个非空字符串
        # 将分数转换为浮点数
        scores = [float(score) for score in scores if score.strip()]  # 过滤掉空字符串
        # 找到最大分数的索引
        max_index = scores.index(max(scores))
        # 保留最大分数对应的名称和InChIKey
        row[pcgroup_idx+24] = str(scores[max_index])
        row[pcgroup_idx+23] = names[max_index] if max_index < len(names) and names[max_index].strip() else ''
        row[pcgroup_idx+26] = inchikeys[max_index] if max_index < len(inchikeys) and inchikeys[max_index].strip() else ''
    
    return row

# 应用filter_columns函数到每一行
merged_df = merged_df.apply(filter_columns, axis=1)

# 对合并后的结果进行筛选
def add_columns(row):
    #print(row[pcgroup_idx+2], row[pcgroup_idx+12], row[pcgroup_idx+19], row[pcgroup_idx+24], row[46])
    if pd.notna(row[pcgroup_idx+2]) and float(row[pcgroup_idx+2]) >= 0.5:  # 第18列（索引17）有值且 >= 0.5
        # name, adduct, inchikey, score, source, annotation
        new_col1, new_col2, new_col3, new_col4, new_col5, new_col6 = row[pcgroup_idx], row[pcgroup_idx+1], row[pcgroup_idx+4], row[pcgroup_idx+2], 'in-house database', 'Level1'  # 自建库注释
    else:
        # 选出最大值及其索引
        #selected_col = max([27, pcgroup_idx+19, 39, 46], key=lambda col: row[col] if pd.notna(row[col]) else float('-inf'))
        valid_values = {col: float(row[col]) for col in [pcgroup_idx+12, pcgroup_idx+19, pcgroup_idx+24, pcgroup_idx+31] if pd.notna(row[col])}
        selected_col = max(valid_values, key=valid_values.get) if valid_values else None
        
        if selected_col == pcgroup_idx+12:
            new_col1, new_col2, new_col3, new_col4, new_col5, new_col6 = row[pcgroup_idx+10], row[pcgroup_idx+11], row[pcgroup_idx+13], row[pcgroup_idx+12], 'GNPS', 'Level2'
        elif selected_col == pcgroup_idx+19:
            new_col1, new_col2, new_col3, new_col4, new_col5, new_col6= row[pcgroup_idx+16], row[pcgroup_idx+17], row[pcgroup_idx+20], row[pcgroup_idx+19], 'MoNA', 'Level2'
        elif selected_col == pcgroup_idx+24:
            p1, p2 = row[pcgroup_idx+23].rsplit(' ', 1) if isinstance(row[pcgroup_idx+23], str) and ' ' in row[pcgroup_idx+23] else (row[pcgroup_idx+23], None)
            new_col1, new_col2, new_col3, new_col4, new_col5, new_col6 = p1, p2, row[pcgroup_idx+26], row[pcgroup_idx+24], 'HMDB', 'Level2'
        elif selected_col == pcgroup_idx+31:
            new_col1, new_col2, new_col3, new_col4, new_col5, new_col6 = row[pcgroup_idx+29], row[pcgroup_idx+30], row[pcgroup_idx+32], row[pcgroup_idx+31], 'lipid', 'Level2'
        else:
            new_col1, new_col2, new_col3, new_col4, new_col5, new_col6 = None, None, None, None, None, None   
    
    return pd.Series([new_col1, new_col2, new_col3, new_col4, new_col5, new_col6])

new_columns = merged_df.apply(add_columns, axis=1)
# 检查哪些行的列数不为 6
invalid_rows = new_columns[new_columns.apply(lambda x: len(x) != 6, axis=1)]
if not invalid_rows.empty:
    print("以下行的返回值不是 6 列：")
    print(invalid_rows)
    
merged_df[['name_filtered', 'adduct_filtered', 'inchikey_filtered', 'score_filtered', 'source_filtered', 'Levels_filtered']] = merged_df.apply(add_columns, axis=1)

# 基于 inchikey_filtered 和 hmdb_inchikey 合并数据框
merged_df = pd.merge(merged_df, hmdb_subset, left_on='inchikey_filtered', right_on='hmdb_inchikey', how='left')

# 统一adduct格式
merged_df['adduct_filtered'] = merged_df['adduct_filtered'].replace(mapping)

# 找出未匹配的行（hmdb_inchikey 为 NaN），且 inchikey_filtered 不为空
unmatched = merged_df['hmdb_inchikey'].isna() & merged_df['inchikey_filtered'].notna()

####这里新增调用classyfireR，通过inchikey查询classification
# 示例：需要查询分类的 InChIKey 列表
inchikeys_to_query = merged_df.loc[unmatched, 'inchikey_filtered'].dropna().unique().tolist()

# 把 InChIKey 写入临时 csv（R 脚本读这个）
inchikey_file = 'inchikeys.csv'
pd.DataFrame({'inchikey': inchikeys_to_query}).to_csv(inchikey_file, index=False)

# R 输出结果保存路径
output_csv = 'classification_result.csv'
r_script_file = 'query_classyfire.R'
inchikey_classyfire = inchikey_classyfire_path

# 创建 .R 脚本
r_script = f"""
library(classyfireR)
library(dplyr)
inchikeys <- read.csv("{inchikey_file}", stringsAsFactors = FALSE)$inchikey



# 读取本地的分类信息
if (file.exists("{inchikey_classyfire}")) {{
  cache <- read.delim("{inchikey_classyfire}", sep = "\t", stringsAsFactors = FALSE)
}} else {{
  cache <- data.frame()
}}

# 定义查询函数
get_info <- function(key) {{
  # 如果缓存中已存在，就直接返回
  if (nrow(cache) > 0 && key %in% cache$inchikey) {{
    row <- cache[cache$inchikey == key, ]
    return(data.frame(inchikey = key,
                      hmdb_kingdom = row$hmdb_kingdom,
                      hmdb_super_class = row$hmdb_super_class,
                      hmdb_class = row$hmdb_class,
                      hmdb_sub_class = row$hmdb_sub_class,
                      stringsAsFactors = FALSE))
  }}

  # 否则调用 classfire 查询
  res <- tryCatch(get_classification(key), error=function(e) NULL)

  if (!is.null(res)) {{
    return(data.frame(inchikey = key,
                      hmdb_kingdom = res@classification$Classification[1],
                      hmdb_super_class = res@classification$Classification[2],
                      hmdb_class = res@classification$Classification[3],
                      hmdb_sub_class = res@classification$Classification[4],
                      stringsAsFactors = FALSE))
  }} else {{
    return(data.frame(inchikey = key,
                      hmdb_kingdom = NA,
                      hmdb_super_class = NA,
                      hmdb_class = NA,
                      hmdb_sub_class = NA,
                      stringsAsFactors = FALSE))
  }}
}}

results <- do.call(rbind, lapply(inchikeys, get_info))
write.csv(results, "{output_csv}", row.names = FALSE)

# 提取 hmdb_kingdom 为 NA 的 InChIKey
na_inchikeys <- results[is.na(results$hmdb_kingdom), "inchikey", drop = FALSE]

# 检查是否有 NA 值
if (nrow(na_inchikeys) > 0) {{
  # 如果有 NA，构建数据框
  output_df <- data.frame(
    InChIKey = na_inchikeys$inchikey,
    Kingdom = "",
    Superclass = "",
    Class = "",
    Subclass = "",
    stringsAsFactors = FALSE
  )
}} else {{
  # 如果没有 NA，生成空数据框（只有列名）
  output_df <- data.frame(
    InChIKey = character(0),
    Kingdom = character(0),
    Superclass = character(0),
    Class = character(0),
    Subclass = character(0),
    stringsAsFactors = FALSE
  )
}}

# 保存为制表符分隔的 txt 文件，NA 替换为空字符串
write.table(
  output_df,
  file = "classyfire_batch.out.txt",
  sep = "\\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = ""
)
"""

# 写入 query_classyfire.R 文件
with open(r_script_file, 'w') as f:
    f.write(r_script)

# 运行 R 脚本
subprocess.run(['Rscript', r_script_file], check=True)

# 删除 R 脚本
os.remove(r_script_file)

# 对未匹配的行，提取 inchikey_filtered 的前 14 个字符，提取hmdb_inchikey的前14个字符进行注释。
merged_df.loc[unmatched, 'inchikey_filtered_14'] = merged_df.loc[unmatched, 'inchikey_filtered'].str[:14]

# 对 hmdb_subset 提取 hmdb_inchikey 的前 14 个字符，并去重取第一个
hmdb_subset['hmdb_inchikey_14'] = hmdb_subset['hmdb_inchikey'].str[:14]
hmdb_subset_unique = hmdb_subset.drop_duplicates(subset='hmdb_inchikey_14', keep='first')

# 第二次合并，仅针对未匹配的行，基于前 14 个字符
unmatched_df = merged_df[unmatched].drop(columns=[col for col in hmdb_subset.columns if col in merged_df.columns])  # 避免列名冲突
unmatched_df = pd.merge(unmatched_df, hmdb_subset_unique, 
                        left_on='inchikey_filtered_14', 
                        right_on='hmdb_inchikey_14', 
                        how='left')

# 将第二次合并的结果更新到 merged_df 中
for col in hmdb_subset_unique.columns:
    if col in unmatched_df.columns:
        merged_df.loc[unmatched, col] = unmatched_df[col].values

# 新增逻辑：如果 hmdb_name 有值，则替换 name_filtered
mask = merged_df['hmdb_name'].notna()  # 筛选 hmdb_name 非空的行
merged_df.loc[mask, 'name_filtered'] = merged_df.loc[mask, 'hmdb_name']

# 清理临时列（可选）
merged_df = merged_df.drop(columns=['inchikey_filtered_14'], errors='ignore')

output_table = "merged_output_no_classyfire.tsv"
merged_df.to_csv(output_table, sep='\t', index=False, encoding='utf-8-sig')
print(f"合并结果已保存至 {output_table}")

