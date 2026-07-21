#!/usr/bin/env python3
import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_one_sample(sample_name, sample_path, input_dir, rmsk_gtf, ref_gtf):
    log_file = os.path.join(input_dir, f"{sample_name}.log")

    cmd = (
        f"/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/velocyto run10x "
        f"-m {rmsk_gtf} {sample_path} "
        f"{ref_gtf} -@ 32 "
        f"> {log_file} 2>&1"
    )

    print(f"[运行] {sample_name}", flush=True)
    print(cmd, flush=True)

    result = subprocess.run(cmd, shell=True)
    return sample_name, result.returncode


def run_velocyto(input_dir, species):
    # 物种判断
    if species.lower() == "human":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

    elif species.lower() == "mouse":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm10_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/refdata-gex-mm10-2020-A/genes/genes.gtf"

    else:
        print("❌ 物种参数错误，只能 human 或 mouse")
        sys.exit(1)

    # 获取样本列表
    samples = []
    for sample_name in os.listdir(input_dir):
        sample_path = os.path.join(input_dir, sample_name)
        if os.path.isdir(sample_path):
            samples.append((sample_name, sample_path))

    if not samples:
        print("❌ 没有找到样本目录")
        sys.exit(1)

    print(f"检测到 {len(samples)} 个样本，全部同时运行...\n", flush=True)

    # 不限制并行数 = 样本数
    with ThreadPoolExecutor(max_workers=len(samples)) as executor:
        futures = [
            executor.submit(
                run_one_sample,
                sample_name,
                sample_path,
                input_dir,
                rmsk_gtf,
                ref_gtf
            )
            for sample_name, sample_path in samples
        ]

        for future in as_completed(futures):
            sample_name, code = future.result()
            if code == 0:
                print(f"✅ {sample_name} 完成", flush=True)
            else:
                print(f"❌ {sample_name} 失败，返回码 {code}", flush=True)

    print("\n🎉 所有样本运行结束")


def main():
    if len(sys.argv) != 3:
        print("用法:")
        print("python run_velocyto.py <输入目录> <human/mouse>")
        sys.exit(1)

    input_dir = sys.argv[1]
    species = sys.argv[2]

    run_velocyto(input_dir, species)


if __name__ == "__main__":
    main()