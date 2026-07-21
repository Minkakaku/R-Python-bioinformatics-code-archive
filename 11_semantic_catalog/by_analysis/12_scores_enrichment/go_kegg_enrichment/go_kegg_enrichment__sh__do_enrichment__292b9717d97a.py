import os
import subprocess

base_dir = "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/Neutrophils/"
species = "mouse"

rscript = "/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/bin/Rscript"
r_code = "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/script/do_enrichment.R"

for level1 in os.listdir(base_dir):
    out = os.path.join(base_dir, level1)
    if not os.path.isdir(out):
        continue

    # for level2 in os.listdir(level1_path):
    #     out = os.path.join(level1_path, level2)
    #     if not os.path.isdir(out):
    #         continue

    deg_file = os.path.join(out, "differentially_expressed_genes.xls")
    if not os.path.exists(deg_file):
        print(f"[SKIP] no DEG file: {out}")
        continue

        kegg_out = os.path.join(out, "KEGG")
        os.makedirs(kegg_out, exist_ok=True)

    cmd = [
        rscript,
        r_code,
        "-f", deg_file,
        "-s", species,
        "-o", out
      ]

    print(f"[RUN] {out}")
    subprocess.run(cmd)