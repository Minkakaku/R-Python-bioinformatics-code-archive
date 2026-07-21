#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess


def get_species_files(species):
    """
    根据物种返回对应的 gtf 文件
    """
    species = species.lower()

    if species == "human":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/GRCh38/genes/genes.gtf"

    elif species == "mouse":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm10_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/mm10_ensemble_102/genes/genes.gtf"

    else:
        print("❌ 物种只能是 human 或 mouse")
        sys.exit(1)

    return rmsk_gtf, ref_gtf


def submit_jobs(input_dir, species):

    rmsk_gtf, ref_gtf = get_species_files(species)

    sample_dirs = [
        d for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d))
    ]

    for sample_name in sample_dirs:

        sample_path = os.path.join(input_dir, sample_name)

        # 输入文件
        in_bam = os.path.join(
            sample_path,
            "Analysis",
            "step2",
            "STAR",
            f"{sample_name}_SortedByCoordinate.bam"
        )

        umi_file = os.path.join(
            sample_path,
            "Analysis",
            "step3",
            "umi.xls"
        )

        barcode_gz = os.path.join(
            sample_path,
            "Outs",
            "filtered_feature_bc_matrix",
            "barcodes.tsv.gz"
        )

        # 输出目录
        velocyto_dir = os.path.join(sample_path, "velocyto")
        os.makedirs(velocyto_dir, exist_ok=True)

        # 输出文件
        out_bam = os.path.join(
            velocyto_dir,
            f"{sample_name}_addtag.bam"
        )

        barcode_tsv = os.path.join(
            velocyto_dir,
            "barcodes.tsv"
        )

        log_file = os.path.join(
            velocyto_dir,
            f"{sample_name}.log"
        )

        # 检查文件
        missing_files = []

        for f in [in_bam, umi_file, barcode_gz]:
            if not os.path.exists(f):
                missing_files.append(f)

        if missing_files:
            print(f"\n⚠️ 跳过 {sample_name}，缺少文件：")
            for f in missing_files:
                print(f"   {f}")
            continue

        # sbatch 内执行命令
        cmd = f"""
            mkdir -p {velocyto_dir} && \
            /mnt/NexData_file/share/scRNA-rawdata/software/seeksoultools.1.3.0a/seeksoultools utils addtag \
            --inbam {in_bam} \
            --umifile {umi_file} \
            --outbam {out_bam} && \
            cp {barcode_gz} {velocyto_dir}/ && \
            gunzip -f {velocyto_dir}/barcodes.tsv.gz && \
            /mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/velocyto run \
            -b {barcode_tsv} \
            -o {velocyto_dir} \
            -m {rmsk_gtf} \
            {out_bam} \
            {ref_gtf} \
            -@ 32 \
            > {log_file} 2>&1
            """

        sbatch_cmd = [
            "sbatch",
            f"--job-name={sample_name}_velocyto",
            "--partition=all",
            "--cpus-per-task=32",
            "--mem=256G",
            "--time=100:00:00",
            f"--output={log_file}",
            f'--wrap={cmd}'
        ]

        print(f"\n🚀 提交样本: {sample_name}")

        result = subprocess.run(
            sbatch_cmd,
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            print(result.stdout.strip())
        else:
            print(f"❌ 提交失败: {sample_name}")
            print(result.stderr)

    print("\n✅ 所有任务提交完成")


def main():

    parser = argparse.ArgumentParser(
        description="批量运行 addtag + velocyto"
    )

    parser.add_argument(
        "input_dir",
        help="输入目录"
    )

    parser.add_argument(
        "species",
        choices=["human", "mouse"],
        help="物种 human/mouse"
    )

    args = parser.parse_args()

    submit_jobs(
        args.input_dir,
        args.species
    )


if __name__ == "__main__":
    main()