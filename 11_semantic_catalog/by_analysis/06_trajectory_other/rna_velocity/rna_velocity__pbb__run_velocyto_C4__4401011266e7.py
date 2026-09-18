#!/usr/bin/env python3
import os
import sys
import gzip
import shutil
import subprocess
import pysam

############################################################
# usage:
# python script.py input_dir output_dir human/mouse
############################################################


def fix_bam(input_bam, output_bam, threads, tmp_bam):

    #tmp_bam = "tmp.unsorted.bam"

    print(f"\n[FIX BAM] {input_bam}")

    infile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(tmp_bam, "wb", template=infile)

    total = 0
    kept = 0

    for read in infile:
        total += 1
        tags = dict(read.get_tags())

        cell = tags.get("DB", None)
        umi = tags.get("UR", None)

        if cell is None or umi is None:
            continue

        read.set_tag("CB", cell, value_type="Z", replace=True)
        read.set_tag("UB", umi, value_type="Z", replace=True)

        outfile.write(read)
        kept += 1

    infile.close()
    outfile.close()

    print(f"Total={total}, Kept={kept}")

    subprocess.run([
        "/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/samtools", "sort",
        "-@", str(threads),
        "-o", output_bam,
        tmp_bam
    ], check=True)

    subprocess.run(["/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/samtools", "index", output_bam], check=True)

    os.remove(tmp_bam)

    return output_bam


def run_velocyto(sample, bam, barcode, outdir, rmsk_gtf, ref_gtf):

    os.makedirs(outdir, exist_ok=True)

    log = os.path.join(outdir, f"{sample}.log")

    cmd = [
        "/mnt/NexData_file/share/scRNA-rawdata/envs/scRawdata/bin/velocyto",
        "run",
        "-b", barcode,
        "-o", outdir,
        "-m", rmsk_gtf,
        bam,
        ref_gtf
    ]

    with open(log, "w") as f:
        subprocess.Popen(cmd, stdout=f, stderr=f)


############################################################
# MAIN
############################################################
def main():

    if len(sys.argv) != 4:
        print("USAGE: python script.py <input_dir> <output_dir> <human|mouse>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    species = sys.argv[3].lower()

    threads = 16

    # reference
    if species == "human":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/GRCh38_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/C4/Human_GRCh38/genes.filter.gtf"

    elif species == "mouse":
        rmsk_gtf = "/mnt/NexData_file/share/scRNA-rawdata/velocyto_gtf/mm39_rmsk.gtf"
        ref_gtf = "/mnt/NexData_file/share/scRNA-rawdata/database/C4/Mouse_GRCm38/genes.filter.gtf"

    else:
        print("species must be human or mouse")
        sys.exit(1)

    for sample in os.listdir(input_dir):

        sample_path = os.path.join(input_dir, sample)
        if not os.path.isdir(sample_path):
            continue

        bam = os.path.join(sample_path, "output", "anno_decon_sorted.bam")
        barcode_gz = os.path.join(sample_path, "output", "filter_matrix", "barcodes.tsv.gz")

        if not os.path.exists(bam) or not os.path.exists(barcode_gz):
            print("[SKIP]", sample)
            continue

        print("\n======================")
        print("SAMPLE:", sample)
        print("======================")

        sample_out_dir = os.path.join(sample_path, "velocyto")
        os.makedirs(sample_out_dir, exist_ok=True)

        velocyto_out = sample_out_dir

        # fixed bam 也放这里
        fixed_bam = os.path.join(sample_out_dir, "fixed.sorted.bam")
        tmp_bam = os.path.join(sample_out_dir, "tmp.unsorted.bam")

        fix_bam(bam, fixed_bam, threads, tmp_bam)

        # barcode
        barcode_tsv = os.path.join(sample_out_dir, "barcodes.tsv")
        with gzip.open(barcode_gz, "rb") as fin, open(barcode_tsv, "wb") as fout:
            shutil.copyfileobj(fin, fout)

        # step3
        run_velocyto(sample, fixed_bam, barcode_tsv, velocyto_out, rmsk_gtf, ref_gtf)


if __name__ == "__main__":
    main()