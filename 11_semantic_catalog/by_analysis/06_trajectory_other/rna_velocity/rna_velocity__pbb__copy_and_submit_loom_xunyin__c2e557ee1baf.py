#!/usr/bin/env python3

import os
import subprocess
import shutil
import sys


def copy_and_submit_sbatch(input_dir):

    # 获取所有样本目录（排除 loom）
    sample_dirs = [
        d for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d))
        and d not in ["loom", "M11_1"]
    ]

    if len(sample_dirs) == 0:
        print(f"No sample directories found in {input_dir}")
        return

    print(f"Found {len(sample_dirs)} samples")

    # loom总目录
    loom_base = os.path.join(input_dir, "loom")
    os.makedirs(loom_base, exist_ok=True)

    for sample in sample_dirs:

        print("=" * 80)
        print(f"Processing sample: {sample}")

        sample_dir = os.path.join(input_dir, sample)

        # barcode文件
        source_barcode = os.path.join(
            sample_dir,
            "Outs",
            "filtered_feature_bc_matrix",
            "barcodes.tsv.gz"
        )

        # loom/sample目录
        target_dir = os.path.join(loom_base, sample)
        os.makedirs(target_dir, exist_ok=True)

        # 目标barcode
        target_gz = os.path.join(target_dir, "barcodes.tsv.gz")
        target_txt = os.path.join(target_dir, "barcodes.tsv")

        # 拷贝并解压barcode
        if not os.path.exists(target_txt):

            if not os.path.exists(source_barcode):
                print(f"ERROR: barcode file not found")
                print(source_barcode)
                continue

            print(f"Copying barcode file")

            shutil.copy2(source_barcode, target_gz)

            print(f"Decompressing barcode file")

            subprocess.run(
                ["gunzip", "-f", target_gz],
                check=True
            )

        else:
            print("Barcode already exists")

        # bam文件
        bam_file = os.path.join(
            target_dir,
            f"{sample}_addtag.bam"
        )

        if not os.path.exists(bam_file):
            print(f"ERROR: bam file not found")
            print(bam_file)
            continue

        # sbatch命令
        sbatch_command = (
            f'sbatch '
            f'--job-name={sample}_loom '
            f'--partition=all '
            f'--output={target_dir}/{sample}_loom.log '
            f'--time=72:00:00 '
            f'--cpus-per-task=32 '
            f'--mem=256G '
            f'--wrap="'
            f'velocyto run '
            f'-b {target_txt} '
            f'-o {target_dir} '
            f'-m /mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf '
            f'{bam_file} '
            f'/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/GRCh38/genes/genes.gtf'
            f'"'
        )

        print("Submitting sbatch job...")
        print(sbatch_command)

        try:

            subprocess.run(
                sbatch_command,
                shell=True,
                check=True
            )

            print(f"Successfully submitted job for {sample}")

        except subprocess.CalledProcessError as e:

            print(f"Failed to submit job for {sample}")
            print(e)

    print("=" * 80)
    print("All jobs submitted")


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("Usage:")
        print("python copy_and_submit_loom.py <input_directory>")
        sys.exit(1)

    input_dir = sys.argv[1]

    if not os.path.exists(input_dir):
        print(f"ERROR: input directory not found")
        print(input_dir)
        sys.exit(1)

    copy_and_submit_sbatch(input_dir)