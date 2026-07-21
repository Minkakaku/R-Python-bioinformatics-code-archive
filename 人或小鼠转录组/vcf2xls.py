import sys
import os

vcf = snakemake.input[0]
vcf_out = snakemake.output[0]

with open(vcf_out, 'w') as fout:
    fout.write(("{}"+"\t{}"*15+"\n").format("CHROM", "POS", "REF", "ALT", "QUAL",
                                            "DP", "AD", "GT", "GeneID", "GeneName",
                                            "FeatureID", "Biotype", "HGVS.c", "HGVS.p",
                                            "Effect", "Impact"))
    with open(vcf, 'r') as fin:
        buf = fin.readline()
        while buf:
            if not buf.startswith('#'):
                items = buf.strip().split("\t")
                CHROM = items[0]
                POS = items[1]
                REF = items[3]
                ALT = items[4]
                QUAL = items[5]
                COUNT = items[9]
                DP = COUNT.split(':')[2]
                AD = COUNT.split(':')[3]
                GT = COUNT.split(':')[0]
                INFO = items[7]
                INFO_items = INFO.split(';')

                for x in INFO_items:
                    if x.startswith("ANN="):
                        ANN = x.replace("ANN=", "")
                        ANN_items = ANN.split(',')
                        for y in ANN_items:
                            y_items = y.split('|')
                            GeneID = y_items[4]
                            GeneName = y_items[3]
                            FeatureID = y_items[6]
                            Biotype = y_items[7]
                            HGVS_c = y_items[9]
                            HGVS_p = y_items[10]
                            Effect = y_items[1]
                            Impact = y_items[2]
                            fout.write(("{}" + "\t{}" * 15 + "\n").format(
                                CHROM, POS, REF, ALT, QUAL,
                                DP, AD, GT, GeneID, GeneName,
                                FeatureID, Biotype, HGVS_c, HGVS_p,
                                Effect, Impact))
                            
            buf = fin.readline()
