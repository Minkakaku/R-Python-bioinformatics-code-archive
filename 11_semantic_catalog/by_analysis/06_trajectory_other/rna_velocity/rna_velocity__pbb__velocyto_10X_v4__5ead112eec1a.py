#!/usr/bin/env python3

import argparse
import subprocess
import sys


def get_gtf(species):
    """
    根据物种返回 gtf
    """

    if species.lower() == "human":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/refdata-gex-GRCh38-2024-A_genes.gtf"

    elif species.lower() == "mouse":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm39_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/refdata-gex-GRCm39-2024-A_genes.gtf"

    else:
        print(f"ERROR: 不支持的物种 {species}")
        sys.exit(1)

    return rmsk_gtf, ref_gtf


def prepare_bam(sample_name, sample_path):
    """
    复制 bam 到 velocyto 需要的位置
    """
    source_bam = os.path.join(
        sample_path,
        "per_sample_outs",
        sample_name,
        "count",
        "sample_alignments.bam"
    )

    target_bam = os.path.join(
        sample_path,
        "possorted_genome_bam.bam"
    )

    if not os.path.exists(source_bam):
        print(f"ERROR: 找不到 bam 文件:\n{source_bam}")
        sys.exit(1)

    # 已存在则跳过
    if os.path.exists(target_bam):
        print(f"已存在: {target_bam}")
        return

    print("移动 bam 文件...")
    print(f"{source_bam}")
    print("->")
    print(f"{target_bam}")

    shutil.move(source_bam, target_bam)

    print("bam 文件移动完成")


def run_velocyto(sample_name, species, sample_path):

    prepare_bam(sample_name, sample_path)
    rmsk_gtf, ref_gtf = get_gtf(species)

    cmd = f"""
/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/velocyto run10x \
-m {rmsk_gtf} \
{sample_path} \
{ref_gtf} \
-@ 32 \
> {sample_name}_loom.log 2>&1
"""

    print("运行命令：")
    print(cmd)

    subprocess.run(cmd, shell=True, check=True)

    #print(f"{sample_name} velocyto 完成")


def main():

    parser = argparse.ArgumentParser(description="10X v3 velocyto流程")
    parser.add_argument("--project_id", required=True)
    parser.add_argument("--sample_name", required=True)
    parser.add_argument(
        "--species",
        required=True,
        choices=["human", "mouse"]
    )
    parser.add_argument("--sample_path", required=True)
    args = parser.parse_args()

    run_velocyto(
        args.sample_name,
        args.species,
        args.sample_path
    )


if __name__ == "__main__":
    main()