import pandas as pd

# 读取5个 TSV 文件
df1 = pd.read_csv('result_CAMERA_annotated_updated.tsv', sep='\t')
df2 = pd.read_csv('result_CAMERA_annotated_GNPS.tsv', sep='\t')
df3 = pd.read_csv('result_CAMERA_annotated_MoNA.tsv', sep='\t')
df4 = pd.read_csv('result_CAMERA_annotated_HMDB.tsv', sep='\t')
df5 = pd.read_csv('result_CAMERA_annotated_lipid.tsv', sep='\t')
df6 = pd.read_csv('/mnt/NexData_file/home/zjh/develop/metabolomics/database/adduct.mod.tsv', sep='\t')
df7 = pd.read_csv('/mnt/NexData_file/home/zjh/develop/metabolomics/database/hmdb_metabolites.csv')

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
            new_col1, new_col2, new_col3, new_col4, new_col5, new_col6 = row[pcgroup_idx+10], row[pcgroup_idx+11], row[28], row[pcgroup_idx+12], 'GNPS', 'Level2'
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

# 对未匹配的行，提取 inchikey_filtered 的前 14 个字符
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

# 输出到 Excel，并利用 xlsxwriter 设置合并单元格、颜色和字体为 Times New Roman
output_file = 'merged_output_raw.xlsx'
with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    # 将数据写入 Excel，从第二行开始写（第一行预留用于自定义表头）
    merged_df.to_excel(writer, index=False, startrow=1)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    
    # 定义统一的字体格式：Times New Roman
    default_font = 'Times New Roman'
    
    # 设置不同背景色的格式（同时指定字体）
    common_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#D9E1F2',  # 浅蓝色
        'border': 1,
        'font_name': default_font
    })
    table1_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#FCE4D6',  # 浅橙色
        'border': 1,
        'font_name': default_font
    })
    table2_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#E2EFDA',  # 浅绿色
        'border': 1,
        'font_name': default_font
    })
    table3_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#FFF2CC',  # 浅黄色
        'border': 1,
        'font_name': default_font
    })
    table4_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#DDEBF7',  # 浅蓝绿色
        'border': 1,
        'font_name': default_font
    })
    table5_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#F8CBAD',  # 浅橙红色
        'border': 1,
        'font_name': default_font
    })
    
    table6_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#C5A6E2',  # 淡紫色
        'border': 1,
        'font_name': default_font
    })

    table7_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'bg_color': '#FADADD',  # 淡紫色
        'border': 1,
        'font_name': default_font
    })
    
    # 设置第一行的合并单元格表头
    # 公共的前15列
    worksheet.merge_range(0, 0, 0, pcgroup_idx-1, '基础信息', common_fmt)
    
    # in-house database
    col_start = pcgroup_idx
    ncols1 = extra1.shape[1]
    col_end = col_start + ncols1 - 1
    worksheet.merge_range(0, col_start, 0, col_end, '自建库', table1_fmt)
    
    # GNPS database
    col_start = col_end + 1
    ncols2 = extra2.shape[1]
    col_end = col_start + ncols2 - 1
    worksheet.merge_range(0, col_start, 0, col_end, 'GNPS', table2_fmt)
    
    # MoNA database
    col_start = col_end + 1
    ncols3 = extra3.shape[1]
    col_end = col_start + ncols3 - 1
    worksheet.merge_range(0, col_start, 0, col_end, 'MoNA', table3_fmt)
    
    # HMDB database
    col_start = col_end + 1
    ncols4 = extra4.shape[1]
    col_end = col_start + ncols4 - 1
    worksheet.merge_range(0, col_start, 0, col_end, 'HMDB', table4_fmt)
    
    # lipid database
    col_start = col_end + 1
    ncols5 = extra5.shape[1]
    col_end = col_start + ncols5 - 1
    worksheet.merge_range(0, col_start, 0, col_end, 'lipid', table5_fmt)
    
    # filtered information
    col_start = col_end + 1
    col_end = col_start + 5
    worksheet.merge_range(0, col_start, 0, col_end, 'Filtered', table6_fmt)

    # annotation information
    col_start = col_end + 1
    col_end = col_start + 12
    worksheet.merge_range(0, col_start, 0, col_end, 'Annotation', table7_fmt)
    
    
    # 设置所有列默认使用 Times New Roman 字体（这里循环遍历所有列）
    total_cols = merged_df.shape[1]
    default_fmt = workbook.add_format({'font_name': default_font})
    for col_num in range(total_cols):
        worksheet.set_column(col_num, col_num, None, default_fmt)

    header_fmt = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'border': 1, 
        'font_name': default_font
    }) 
    ncols = merged_df.shape[1]
    for col in range(ncols):
        header_val = merged_df.columns[col]
        worksheet.write(1, col, header_val, header_fmt)

    
print(f"合并结果已保存至 {output_file}")

output_table = "merged_output_test.tsv"
merged_df.to_csv(output_table, sep='\t', index=False, encoding='utf-8-sig')
print(f"合并结果已保存至 {output_table}")

