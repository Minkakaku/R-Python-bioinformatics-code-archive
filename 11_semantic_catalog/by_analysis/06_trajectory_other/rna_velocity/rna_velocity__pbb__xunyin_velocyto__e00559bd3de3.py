#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess


def get_species_files(species):

    species = species.lower()

    if species == "human":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/GRCh38/genes/genes.gtf"

    elif species == "mouse":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm10_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/seeksoul/mm10_ensemble_102/genes/genes.gtf"

    else:
        print("❌ species must be human or mouse")
        sys.exit(1)

    return rmsk_gtf, ref_gtf


def run_velocyto(sample_path, species, capture_type):

    sample_name = os.path.basename(sample_path.rstrip("/"))

    rmsk_gtf, ref_gtf = get_species_files(species)

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

    velocyto_dir = os.path.join(sample_path, "velocyto")
    os.makedirs(velocyto_dir, exist_ok=True)

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

    missing_files = []

    for f in [in_bam, umi_file, barcode_gz]:
        if not os.path.exists(f):
            missing_files.append(f)

    if missing_files:

        print(f"\n❌ 缺少文件:")

        for f in missing_files:
            print(f)

        sys.exit(1)

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
    -@ 32
    """

    print(cmd)

    subprocess.run(
        cmd,
        shell=True,
        executable="/bin/bash",
        check=True
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample_path",
        required=True
    )

    parser.add_argument(
        "--species",
        required=True,
        choices=["human", "mouse"]
    )

    parser.add_argument(
        "--capture_type",
        required=True
    )

    args = parser.parse_args()

    run_velocyto(
        args.sample_path,
        args.species,
        args.capture_type
    )


if __name__ == "__main__":
    main()