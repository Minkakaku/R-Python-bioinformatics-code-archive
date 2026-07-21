# Script to extract statistics from fastp json output

import sys
import json
import os
# import pandas as pd
import glob

# Input
input_dir = sys.argv[1]

# Output
output_file = sys.argv[2]

# Read sample list
sample_list = [os.path.basename(f).replace("_fastp.json", "")
               for f in glob.glob(os.path.join(input_dir, "*_fastp.json"))]
sample_list.sort()

# Read fastp json output
with open(output_file, 'w') as fout:
    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("sample", "clean_reads_length_bp", "raw_data_read",
                                                                 "raw_data_bp", "filtered_reads_percentage",
                                                                 "clean_data_read", "clean_data_bp", "clean_data_gc",
                                                                 "clean_data_q20", "clean_data_q30"))
    for s in sample_list:
        json_file = os.path.join(input_dir, "{}_fastp.json".format(s, s))
        with open(json_file, 'r') as j:
            data = json.load(j)
            mean_length = "({}:{})".format(data["summary"]["after_filtering"]["read1_mean_length"],
                                           data["summary"]["after_filtering"]["read2_mean_length"])
            raw_data_read = int(data["summary"]["before_filtering"]["total_reads"])
            raw_data_bp = str(round(int(data["summary"]["before_filtering"]["total_bases"])/1000000000, 2)) + "G"
            clean_data_read = int(data["summary"]["after_filtering"]["total_reads"])
            clean_data_bp = str(round(int(data["summary"]["after_filtering"]["total_bases"])/1000000000, 2)) + "G"
            clean_data_gc = round(float(data["summary"]["after_filtering"]["gc_content"])*100, 2)
            clean_data_q20 = round(float(data["summary"]["after_filtering"]["q20_rate"])*100, 2)
            clean_data_q30 = round(float(data["summary"]["after_filtering"]["q30_rate"])*100, 2)
            raw_data_bp_int = int(data["summary"]["before_filtering"]["total_bases"])
            clean_data_bp_int = int(data["summary"]["after_filtering"]["total_bases"])
            filter_bp_percent = round((raw_data_bp_int - clean_data_bp_int) / raw_data_bp_int * 100, 2)
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(s, mean_length, raw_data_read, raw_data_bp,
                                                                         filter_bp_percent, clean_data_read, clean_data_bp,
                                                                         clean_data_gc, clean_data_q20, clean_data_q30))




