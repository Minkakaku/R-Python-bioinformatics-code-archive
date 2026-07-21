import pandas as pd

#获得每个注释的区间范围
df1 = pd.read_csv('result_CAMERA_annotated_inhouse_updated.tsv', sep='\t')
df2 = pd.read_csv('result_CAMERA_annotated_GNPS.tsv', sep='\t')
df3 = pd.read_csv('result_CAMERA_annotated_MoNA.tsv', sep='\t')
df4 = pd.read_csv('result_CAMERA_annotated_HMDB.tsv', sep='\t')
df5 = pd.read_csv('result_CAMERA_annotated_lipid.tsv', sep='\t')

# 动态找出 pcgroup 列所在的位置（包括它本身）
pcgroup_idx = df1.columns.get_loc('pcgroup') + 1  # 加1表示包含这一列
common = df1.iloc[:, :pcgroup_idx] 

# 分别提取每个文件中公共列以后的数据
extra1 = df1.iloc[:, pcgroup_idx:]
extra2 = df2.iloc[:, pcgroup_idx:]
extra3 = df3.iloc[:, pcgroup_idx:]
extra4 = df4.iloc[:, pcgroup_idx:]
extra5 = df5.iloc[:, pcgroup_idx:]


# 读取 R 脚本的分类结果
classification_df = pd.read_csv("/mnt/NexData_file/home/zjh/project/Metabolomics/workflow/database/inchikey_classyfire.txt",sep="\t")

# 构建映射：以 inchikey 为键
classification_map = classification_df.set_index('inchikey')

# 首先读取列名
with open("merged_output_no_classyfire.tsv", 'r') as f:
    header = f.readline().strip('\n').split('\t')
    data = []
    for line in f:
        data.append(line.strip('\n').split('\t'))

# 创建DataFrame，允许重复列名
merged_df = pd.DataFrame(data, columns=header)
unmatched = merged_df['hmdb_inchikey'].isna() & merged_df['inchikey_filtered'].notna()

# 对于 unmatched 且 inchikey_filtered 不为空的行，按行补全字段
for col in ['hmdb_kingdom', 'hmdb_super_class', 'hmdb_class', 'hmdb_sub_class']:
    merged_df.loc[unmatched, col] = merged_df.loc[unmatched, 'inchikey_filtered'].map(classification_map[col])

# 输出到 Excel，并利用 xlsxwriter 设置合并单元格、颜色和字体为 Times New Roman
output_file = 'merged_output_added_classyfire.xlsx'
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

output_table = "merged_output_added_classyfire.tsv"
merged_df.to_csv(output_table, sep='\t', index=False, encoding='utf-8-sig')
print(f"合并结果已保存至 {output_table}")
