import pandas as pd
import sys
import argparse


def process_dataframe(file_path, output_file, output_file_filt, num_samples, enable_rt_filter=False):
    # 读取数据
    df = pd.read_csv(file_path, sep="\t", encoding='utf-8-sig')

    # 计算 intensity
    intensity_cols = df.columns[9:9+num_samples]
    df["intensity"] = df[intensity_cols].mean(axis=1)

    # 生成 Key
    def generate_key(row):
        mz_formatted = f"{row['mz']:.1f}"
        if pd.notna(row["inchikey_filtered"]):
            return f"{row['inchikey_filtered']}_{row['adduct_filtered']}_{mz_formatted}"
        return ""

    df["Key"] = df.apply(generate_key, axis=1)

    # 仅当 Key 不为空时计算 score_max 和 intensity_max
    df["score_max"] = ""
    df["intensity_max"] = ""
    df["Count"] = ""

    valid_keys = df[df["Key"] != ""]

    if not valid_keys.empty:
        valid_scores = valid_keys.dropna(subset=["score_filtered"])
        if not valid_scores.empty:
            df.loc[valid_scores.groupby("Key")["score_filtered"].idxmax(), "score_max"] = "Yes"

        valid_intensity = valid_keys.dropna(subset=["intensity"])
        if not valid_intensity.empty:
            df.loc[valid_intensity.groupby("Key")["intensity"].idxmax(), "intensity_max"] = "Yes"

        df.loc[df["Key"] != "", "Count"] = (df[["score_max", "intensity_max"]] == "Yes").sum(axis=1)

    # 删除指定范围的列
    drop_cols = df.columns[9+num_samples+3:9+num_samples+37]
    df = df.drop(columns=drop_cols, errors='ignore')

    df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
    df.loc[df["Count"] > 0, "Key2"] = df.loc[df["Count"] > 0].apply(
        lambda row: row["inchikey_filtered"] + f'_{row["rt"]:.1f}', axis=1
    )

    df["Count_rt"] = df["Key2"].map(df["Key2"].value_counts()).fillna(0).astype(int)

    # 保存结果
    df.to_csv(output_file, sep="\t", index=False)

    # 第一步：保留 Count >= 1 的行
    df = df[df['Count'] >= 1]

    # 第二步：根据条件过滤 rt 列（仅 RP 启用）
    if enable_rt_filter:
        df = df[
            (
                (df['hmdb_super_class'] == "Lipids and lipid-like molecules") | (df['source_filtered'] == "lipid")
            ) & (df['rt'] <= 660)
            |
            (
                ((df['hmdb_super_class'].isna()) & (df['source_filtered'] != "lipid")) | (df['hmdb_super_class'] != "Lipids and lipid-like molecules")
            ) & (df['rt'] <= 540)
        ]

    df.to_csv(output_file_filt, sep="\t", index=False)

    return output_file, output_file_filt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter after merge annotation")
    parser.add_argument("input_file", help="Input TSV file")
    parser.add_argument("output_clean", help="Output clean TSV file")
    parser.add_argument("output_filt", help="Output filtered TSV file")
    parser.add_argument("num_samples", type=int, help="Number of samples")
    parser.add_argument("--enable_rt_filter", action="store_true", default=False,
                        help="Enable RT-based filtering (used for RP, disabled for HILIC)")

    args = parser.parse_args()
    process_dataframe(args.input_file, args.output_clean, args.output_filt,
                      args.num_samples, args.enable_rt_filter)
