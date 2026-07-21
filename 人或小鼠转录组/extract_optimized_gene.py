# Pyhthon script to extract genes with optimized structure from cuffcompare results

import gffutils


def fetch_gene_info(gene_name, db, name2id_dict):
    gid = name2id_dict[gene_name]
    gene_query = db[gid]
    pos = "{chrom}[{strand}]{start}-{end}".format(chrom=gene_query.chrom,
                                                  strand=gene_query.strand,
                                                  start=gene_query.start,
                                                  end=gene_query.end)
    return pos, gid


if __name__ == "__main__":
    cuffcmp_loci = snakemake.input[0]
    gffutils_db = snakemake.params["gffutils_db"]
    output = snakemake.output[0]

    # load gffutils database
    db = gffutils.FeatureDB(gffutils_db)

    # create a gene_name: gene_id dictionary
    gene_name2id_dict = {}
    for gene in db.all_features(featuretype="gene"):
        gene_id = gene.attributes["gene_id"][0]
        if "gene_name" in gene.attributes:
            gene_name = gene.attributes["gene_name"][0]
        else:
            gene_name = gene_id
        if gene_id not in gene_name2id_dict:
            gene_name2id_dict[gene_name] = gene_id

    # read cuffcompare .loci result file, convert format and output
    with open(output, 'w') as fout:
        fout.write("{}\t{}\t{}\t{}\t{}\n".format("Gene_ID", "Chromosome", "Strand", "Original_span", "Assembled_span"))
        with open(cuffcmp_loci, 'r') as f:
            buf = f.readline().strip()
            while buf:
                items = buf.split("\t")
                new_gene_id = items[0]
                new_gene_pos = items[1]
                new_gene_chr = new_gene_pos.split('[')[0]
                new_gene_strand = new_gene_pos.split('[')[1].split(']')[0]
                new_gene_start_end = new_gene_pos.split(']')[1]

                old_gene_list = [*set([x.split('|')[0] for x in items[2].split(',')])]
                if old_gene_list == ['-']:
                    buf = f.readline().strip()
                    continue
                elif len(old_gene_list) == 1:
                    old_gene_name = old_gene_list[0]
                    old_gene_pos, old_gene_id = fetch_gene_info(old_gene_name, db, gene_name2id_dict)
                    if old_gene_pos == new_gene_pos:
                        buf = f.readline().strip()
                        continue
                    else:
                        fout.write("{}\t{}\t{}\t{}\t{}\n".format(new_gene_id,
                                                                 new_gene_chr,
                                                                 new_gene_strand,
                                                                 "{}:{}".format(old_gene_id,
                                                                                old_gene_pos.split(']')[1]),
                                                                 new_gene_start_end))
                else:
                    for g in old_gene_list:
                        old_gene_name = g
                        old_gene_pos, old_gene_id = fetch_gene_info(old_gene_name, db, gene_name2id_dict)
                        fout.write("{}\t{}\t{}\t{}\t{}\n".format(new_gene_id,
                                                                 new_gene_chr,
                                                                 new_gene_strand,
                                                                 "{}:{}".format(old_gene_id,
                                                                                old_gene_pos.split(']')[1]),
                                                                 new_gene_start_end))
                buf = f.readline().strip()
