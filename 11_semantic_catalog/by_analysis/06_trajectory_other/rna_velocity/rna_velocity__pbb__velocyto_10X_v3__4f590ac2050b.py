#!/usr/bin/env python3

import argparse
import subprocess
import sys


def get_gtf(species):
    """
    根据物种返回 gtf
    """

    if species == "human":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

    elif species == "mouse":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm10_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/refdata-gex-mm10-2020-A/genes/genes.gtf"

    else:
        print(f"ERROR: 不支持的物种 {species}")
        sys.exit(1)

    return rmsk_gtf, ref_gtf


def run_velocyto(sample_name, species, sample_path):
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