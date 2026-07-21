import pandas as pd
import glob
import os
from functools import reduce

# 设置文件所在的目录（这里假设脚本和 .tmp 文件在同一目录）
file_dir = "."  # 当前目录
# 如果要处理特定文件，可以直接给出列表：
# file_list = ["cn_prop.tmp", "cn_spa.tmp", "ct_prop.tmp",
#              "ct_region_den.tmp", "ct_region.tmp", "ct_spa.tmp"]
# 否则自动获取所有 .tmp 文件
file_list = glob.glob(os.path.join(file_dir, "*.tmp"))

# 如果文件列表为空，给出提示
if not file_list:
    raise FileNotFoundError("未找到任何 .tmp 文件，请检查路径。")

data_frames = []  # 保存处理后的每个表格

for fpath in file_list:
    # 读取文件，无表头，所有列作为字符串读入（避免数字类型冲突）
    df = pd.read_csv(fpath, sep='\t', header=None, dtype=str, keep_default_na=False)
    
    # 最后一行是真实的表头，其余为数据
    col_names = df.iloc[-1].tolist()      # 最后一行的所有字段
    data = df.iloc[:-1]                   # 上面所有行是数据
    
    # 将表头赋值给数据部分
    data.columns = col_names
    
    # 提取文件名前缀（去掉 .tmp 扩展名）
    prefix = os.path.splitext(os.path.basename(fpath))[0]
    
    # 构建新列名：第一列保持为 "sample"，其他列加上前缀
    # 假设第一列的表头名称确实是 "sample"；如果不是也可以保留原样，其他列则加前缀
    first_col = col_names[0]
    new_columns = [first_col] + [f"{prefix}_{c}" for c in col_names[1:]]
    data.columns = new_columns
    
    # 去除可能的重复 sample 行（如果有，保留第一个）
    data = data.drop_duplicates(subset=first_col, keep='first')
    
    data_frames.append(data)

# 按 sample 列合并所有数据（外连接，保留所有样本）
if len(data_frames) == 1:
    merged = data_frames[0]
else:
    # 使用 reduce 依次 merge
    merged = reduce(lambda left, right: pd.merge(left, right, on=first_col, how='outer', sort=False),
                    data_frames)

# 保存整合后的文件
output_file = "integrated.tsv"
merged.to_csv(output_file, sep='\t', index=False, na_rep='')

print(f"整合完成！共处理 {len(file_list)} 个文件，输出 → {output_file}")
print(f"总样本数：{len(merged)}，总列数（含sample）：{len(merged.columns)}")
