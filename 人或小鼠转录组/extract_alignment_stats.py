# Script to extract alignment statistics from HISAT2 summary results

import sys
import re
import os

hisat2_summary_file = sys.argv[1]
mapping_stats_output = sys.argv[2]

sample = os.path.basename(hisat2_summary_file).replace("_hisat2_alignment.summary", "")
total_reads_n = 0
total_mapped_reads_n = 0
total_unique_mapped_reads_n = 0
total_multiple_mapped_reads_n = 0

with open(hisat2_summary_file, 'r') as fin:
    # skip first line
    buf = fin.readline()
    buf = fin.readline()
    while buf:
        items = buf.strip().split(": ")
        print(items)
        description = items[0]
        count = items[1].split(" (")[0]
        if description == "Total pairs":
            total_reads_n = int(count) * 2
        elif description == "Aligned concordantly 1 time":
            total_unique_mapped_reads_n += int(count) * 2
        elif description == "Aligned 1 time":
            total_unique_mapped_reads_n += int(count)
        elif description == "Aligned concordantly >1 times":
            total_multiple_mapped_reads_n += int(count) * 2
        elif description == "Aligned >1 times":
            total_multiple_mapped_reads_n += int(count)
        buf = fin.readline()

total_mapped_reads_n = total_unique_mapped_reads_n + total_multiple_mapped_reads_n
total_mapped_reads_pct = round(total_mapped_reads_n / total_reads_n * 100, 2)
total_unique_mapped_reads_pct = round(total_unique_mapped_reads_n / total_reads_n * 100, 2)
total_multiple_mapped_reads_pct = round(total_multiple_mapped_reads_n / total_reads_n * 100, 2)

total_mapped_reads = "{}({}%)".format(total_mapped_reads_n, total_mapped_reads_pct)
total_unique_mapped_reads = "{}({}%)".format(total_unique_mapped_reads_n, total_unique_mapped_reads_pct)
total_multiple_mapped_reads = "{}({}%)".format(total_multiple_mapped_reads_n, total_multiple_mapped_reads_pct)

with open(mapping_stats_output, 'w') as fout:
    fout.write("\t".join(["Sample", "Total_Reads", "Total_Map", "Unique_Map", "Multiple_Map"]) + "\n")
    fout.write("\t".join([sample, str(total_reads_n), total_mapped_reads,
                          total_unique_mapped_reads, total_multiple_mapped_reads]) + "\n")
