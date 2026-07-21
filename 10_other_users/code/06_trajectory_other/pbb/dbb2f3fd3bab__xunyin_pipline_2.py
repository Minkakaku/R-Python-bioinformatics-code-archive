#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
import pandas as pd
from datetime import datetime

def generate_date_dir(base_dir, pid):
    current_date = datetime.now().strftime("%Y.%m.%d")

    out_dir = os.path.join(base_dir, pid, current_date)

    counter = 1
    while os.path.exists(out_dir):
        out_dir = os.path.join(base_dir, pid, f"{current_date}.{counter}")
        counter += 1

    return out_dir

def get_base_dir(project_id):
    if project_id.startswith("SC26"):
        return "/mnt/NexData_file/nas/rawdata/rawdata-sc-2026"
    elif project_id.startswith("SC25"):
        return "/mnt/NexData_file/nas/rawdata/rawdata-sc-2025"
    elif project_id.startswith("SC24"):
        return "/mnt/NexData_file/nas/rawdata/rawdata-sc-2024"
    else:
        raise ValueError(f"不支持的项目号: {project_id}")


def get_reference(species):
    if species == "human":
        genome_dir = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/GRCh38/star"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/GRCh38/genes/genes.gtf"
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        star_path_arg = ""

    elif species == "mouse":
        genome_dir = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/mm10_ensemble_102/star"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/mm10_ensemble_102/genes/genes.gtf"
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm10_rmsk.gtf"
        star_path_arg = ""

    else:
        raise ValueError(f"不支持的物种: {species}")

    return genome_dir, ref_gtf, rmsk_gtf, star_path_arg


def submit_job(pid, sample_name, species, customer_name, customer_company, out_dir):

    base_dir = get_base_dir(pid)
    genome_dir, ref_gtf, rmsk_gtf, star_path_arg = get_reference(species)

    #out_base = "/mnt/NexData_file/nas/outputdata/cellranger_out"
    #out_dir = os.path.join(out_base, pid)
    #os.makedirs(out_dir, exist_ok=True)

    sample_raw_dir = os.path.join(base_dir, pid, sample_name)

    fq1 = os.path.join(sample_raw_dir, "*R1_001.fastq.gz")
    fq2 = os.path.join(sample_raw_dir, "*R2_001.fastq.gz")

    r1_files = sorted(glob.glob(fq1))
    r2_files = sorted(glob.glob(fq2))

    if not (r1_files and r2_files):
        fq1 = os.path.join(sample_raw_dir, "*1.fq.gz")
        fq2 = os.path.join(sample_raw_dir, "*2.fq.gz")

        r1_files = sorted(glob.glob(fq1))
        r2_files = sorted(glob.glob(fq2))

    if not (r1_files and r2_files):
        print(f"❌ 缺少fastq文件: {sample_raw_dir}")
        return

    sample_out_dir = os.path.join(out_dir, sample_name)
    report_dir = os.path.join(sample_out_dir, "simple_report")
    velocyto_dir = os.path.join(sample_out_dir, "velocyto")

    xunyin_cmd = (
        f"/mnt/NexData_file/share/scRNA-rawdata/software/seeksoultools.1.3.0a/seeksoultools rna run "
        f"--fq1 {fq1} "
        f"--fq2 {fq2} "
        f"--samplename {sample_name} "
        f"--genomeDir {genome_dir} "
        f"--gtf {ref_gtf} "
        f"--outdir {out_dir} "
        f"--chemistry DDV2 "
        f"--core 32 "
        f"--include-introns "
        f"{star_path_arg}"
    )

    simple_report = "/mnt/NexData_file/share/scRNA-rawdata/script/abridged_report/simple_report.py"

    report_cmd = (
        f"/mnt/NexData_file/share/scRNA-rawdata/envs/scRNA_p/bin/python "
        f"{simple_report} "
        f"--platform xunyin "
        f"--pid {pid} "
        f"--psample {sample_name} "
        f"--paddress '{customer_company}' "
        f"--pcustmor '{customer_name}' "
        f"--input_dir {sample_out_dir} "
        f"--sample_name {sample_name} "
        f"--output_dir {report_dir} "
        f"--organism {species}"
    )

    in_bam = os.path.join(sample_out_dir, "Analysis", "step2", "STAR", f"{sample_name}_SortedByCoordinate.bam")
    umi_file = os.path.join(sample_out_dir, "Analysis", "step3", "umi.xls")
    barcode_gz = os.path.join(sample_out_dir, "Outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")

    out_bam = os.path.join(velocyto_dir, f"{sample_name}_addtag.bam")
    barcode_tsv = os.path.join(velocyto_dir, "barcodes.tsv")

    addtag_cmd = (
        f"/mnt/NexData_file/share/scRNA-rawdata/software/seeksoultools.1.3.0a/seeksoultools utils addtag "
        f"--inbam {in_bam} "
        f"--umifile {umi_file} "
        f"--outbam {out_bam}"
    )

    velocyto_cmd = (
        f"/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/velocyto run "
        f"-b {barcode_tsv} "
        f"-o {velocyto_dir} "
        f"-m {rmsk_gtf} "
        f"{out_bam} "
        f"{ref_gtf} "
        f"-@ 32"
    )

    pipeline_cmd = (
        f"mkdir -p {report_dir} {velocyto_dir} && "
        f"{xunyin_cmd} && "
        f"{report_cmd} && "
        f"{addtag_cmd} && "
        f"cp {barcode_gz} {velocyto_dir}/ && "
        f"gunzip -f {velocyto_dir}/barcodes.tsv.gz && "
        f"{velocyto_cmd}"
    )

    log_file = os.path.join(out_dir, f"{sample_name}.log")

    sbatch_cmd = [
        "sbatch",
        f"--job-name={sample_name}",
        "--partition=all",
        "--cpus-per-task=32",
        "--mem=256G",
        "--time=500:00:00",
        f"--output={log_file}",
        f"--wrap={pipeline_cmd}"
    ]

    print(f"\n===== 提交样本: {sample_name} =====")

    result = subprocess.run(sbatch_cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(result.stdout.strip())
    else:
        print(result.stderr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="SeekSoul RNA + Report + AddTag + Velocyto")
    parser.add_argument("xlsx_file", help="样本信息表")

    args = parser.parse_args()

    df = pd.read_excel(args.xlsx_file)

    species_map = {
        "人": "human",
        "小鼠": "mouse"
    }

    required_columns = ["项目编号", "物种中文名", "样本分析名", "项目联系人", "客户单位"]

    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Excel缺少列: {col}")

    out_base = "/mnt/NexData_file/nas/outputdata/cellranger_out"

    project_outdir = {}

    for _, row in df.iterrows():

        pid = str(row["项目编号"]).strip().split("-")[0]

        if pid not in project_outdir:
            project_outdir[pid] = generate_date_dir(out_base, pid)
            os.makedirs(project_outdir[pid], exist_ok=True)

        species_cn = str(row["物种中文名"]).strip()

        if species_cn not in species_map:
            print(f"❌ 跳过，不支持物种: {species_cn}")
            continue

        species = species_map[species_cn]

        sample_name = str(row["样本分析名"]).strip()
        customer_name = str(row["项目联系人"]).strip()
        customer_company = str(row["客户单位"]).strip()

        submit_job(
            pid,
            sample_name,
            species,
            customer_name,
            customer_company,
            project_outdir[pid]
        )