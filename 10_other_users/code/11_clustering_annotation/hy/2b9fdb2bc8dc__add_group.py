#!/usr/bin/env python3
"""
根据 cell_id 和 sample 匹配，为 all.csv 添加 group 列
"""

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='为all.csv添加group列')
    parser.add_argument('--all', default='all.csv', help='all.csv路径')
    parser.add_argument('--rrr', default='rrr.csv', help='rrr文件路径（包含cell_id,sample,group）')
    parser.add_argument('--output', default='all_with_group.csv', help='输出文件路径')
    parser.add_argument('--sep-all', default=',', help='all.csv分隔符')
    parser.add_argument('--sep-rrr', default='\t', help='rrr文件分隔符（默认制表符）')
    args = parser.parse_args()

    # 读取两个文件
    df_all = pd.read_csv(args.all, sep=args.sep_all)
    df_rrr = pd.read_csv(args.rrr, sep=args.sep_rrr)

    # 确保列名一致（去除可能存在的空格）
    df_all.columns = df_all.columns.str.strip()
    df_rrr.columns = df_rrr.columns.str.strip()

    # 检查必要列
    if not {'cell_id', 'sample', 'group'}.issubset(df_rrr.columns):
        raise ValueError("rrr文件必须包含 cell_id, sample, group 三列")
    if not {'cell_id', 'sample'}.issubset(df_all.columns):
        raise ValueError("all.csv必须包含 cell_id, sample 列")

    # 左连接匹配
    merged = df_all.merge(df_rrr[['cell_id', 'sample', 'group']],
                          on=['cell_id', 'sample'],
                          how='left')

    # 保存结果
    merged.to_csv(args.output, index=False)
    print(f"处理完成，结果已保存至: {args.output}")
    print(f"总行数: {len(merged)}，成功匹配到group的行数: {merged['group'].notna().sum()}")

if __name__ == '__main__':
    main()
