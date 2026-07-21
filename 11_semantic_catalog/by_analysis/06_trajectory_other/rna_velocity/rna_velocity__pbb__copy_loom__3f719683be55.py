import os
import shutil
import argparse

def copy_loom_files(project_dir, output_dir):
    if not os.path.exists(project_dir):
        raise ValueError(f"项目目录不存在: {project_dir}")

    os.makedirs(output_dir, exist_ok=True)

    copied = 0
    skipped = 0

    for sample in os.listdir(project_dir):
        sample_path = os.path.join(project_dir, sample)

        # 只处理目录
        if not os.path.isdir(sample_path):
            continue

        loom_path = os.path.join(sample_path, "velocyto")

        if not os.path.exists(loom_path):
            print(f"[跳过] 未找到velocyto目录: {sample}")
            skipped += 1
            continue

        # 查找 loom 文件
        loom_files = [f for f in os.listdir(loom_path) if f.endswith(".loom")]

        if not loom_files:
            print(f"[跳过] 未找到loom文件: {sample}")
            skipped += 1
            continue

        # 为每个样本创建单独输出目录
        sample_output_dir = os.path.join(output_dir, sample)
        os.makedirs(sample_output_dir, exist_ok=True)

        # 默认取第一个 loom 文件
        src_file = os.path.join(loom_path, loom_files[0])

        # 保留原文件名
        dest_file = os.path.join(sample_output_dir, loom_files[0])

        try:
            shutil.copy2(src_file, dest_file)
            print(f"[完成] {sample} -> {dest_file}")
            copied += 1
        except Exception as e:
            print(f"[失败] {sample}: {e}")
            skipped += 1

    print("\n=== 统计 ===")
    print(f"成功: {copied}")
    print(f"跳过/失败: {skipped}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="批量复制velocyto loom文件")
    parser.add_argument("-i", "--input", required=True, help="项目目录")
    parser.add_argument("-o", "--output", required=True, help="输出目录")

    args = parser.parse_args()

    copy_loom_files(args.input, args.output)