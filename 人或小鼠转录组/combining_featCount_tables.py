"""
This script creates a single TSV feature count table from the featureCounts
output tables in the target directory.

This combined feature count table can be used for differential expression
analysis (e.g. using DESeq2 or edgeR in R).
"""

# Loading required libraries
from time import gmtime, strftime
import pandas as pd
import argparse
import sys
import os
import glob

# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="Combines the featureCounts output tables in the target \
    directory.")

# -d PATH -o CUSTOM_FILENAME
parser.add_argument("-d", "--directory", dest="path",
                    help="path to target directory.")
parser.add_argument("-o", "--output", dest="custom_filename",
                    help="output filename.")
parser.add_argument("-m", "--metadata", dest="metadata",
                    help="path to metadata file.")

args = parser.parse_args()

# Changing to the target directory
path = args.path
# Loading sample order
df_header = ['Gene']
# samples = [os.path.basename(f).replace("feature_counts_", "")
#            for f in glob.glob(os.path.join(path, "feature_counts_*"))]
metadata = args.metadata
samples = list(pd.read_csv(metadata, header=0, sep="\t")["sample"])
df_header.extend(samples)
# Creating variables
fixed_headers = ["Gene", "Chromosome", "Start", "End", "Strand", "Length"]
target_file_prefix = "feature_counts_"
counts_table = pd.DataFrame()
output_filename = str()
target_file_count = 0
# Iterating through files in target directory, combining feature counts
# into one DataFrame object ("counts_table")
for filename in os.listdir(path):
    if filename.startswith(target_file_prefix):
        target_file_count = target_file_count + 1
        srr = filename.replace(target_file_prefix, "")
        each_file = os.path.join(path, filename)
        featCounts_df = pd.read_csv(each_file, sep="\t",
                                    lineterminator="\n",
                                    skiprows=1,
                                    header=0)
        featCounts_headers = fixed_headers.copy()
        featCounts_headers.append(srr)
        featCounts_df.columns = featCounts_headers
        gene_ids = featCounts_df["Gene"]
        counts = featCounts_df[srr]
        # Add the gene IDs and counts to the counts_table DataFrame as columns
        # if it's empty; otherwise add the counts only
        if counts_table.empty:
            counts_table = pd.concat([gene_ids, counts], axis=1, sort=False)
        else:
            counts_table = pd.concat([counts_table, counts], axis=1, sort=False)
        del featCounts_headers

counts_table = counts_table[df_header]
if target_file_count == 0:
    # Exiting script if there are no target files in the target directory
    print("ERROR: There are no featureCount files in the target directory. \n")
    parser.print_help(sys.stderr)
    sys.exit()
else:
    # Exporting counts_table DataFrame as a CSV file
    output_filename = args.custom_filename
    counts_table.to_csv(output_filename, index=False, sep='\t')
