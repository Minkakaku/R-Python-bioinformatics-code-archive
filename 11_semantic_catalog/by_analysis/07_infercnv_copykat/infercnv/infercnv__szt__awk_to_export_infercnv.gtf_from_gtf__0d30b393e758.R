# Gene ordering file
# The gene ordering file provides the chromosomal location for each gene. 
# The format is tab-delimited and has no column header,
# simply providing the gene name, chromosome, and gene span:
#   
# WASH7P  chr1    14363   29806
# LINC00115       chr1    761586  762902
# NOC2L   chr1    879584  894689
# MIR200A chr1    1103243 1103332
# SDF4    chr1    1152288 1167411
# UBE2J2  chr1    1189289 1209265

cat GRCh38.gtf | grep -v '#' | \
  awk -F '\t' '{if($3 == "gene"){ split($9, a, ";\\s+");  for(i in a) { if(a[i] ~ /^gene_name/) {split(a[i], b, "\\s+"); print  b[2] "\t" $1 "\t" $4 "\t" $5 }} }}' | \
  sed -e 's/\"//g' > GRCh38.infercnv.gtf


cat mm10.gtf | grep -v '#' | \
  awk -F '\t' '{if($3 == "gene"){ split($9, a, ";\\s+");  for(i in a) { if(a[i] ~ /^gene_name/) {split(a[i], b, "\\s+"); print  b[2] "\t" $1 "\t" $4 "\t" $5 }} }}' | \
  sed -e 's/\"//g' > mm10.infercnv.gtf


# drop duplcated
read_tsv("gtf/GRCh38.infercnv.gtf", col_names = FALSE) %>% 
  dplyr::distinct(X1, .keep_all = TRUE) %>%
  write_tsv("gtf/GRCh38.infercnv.dropdup.gtf",col_names = FALSE)

read_tsv("gtf/mm10.infercnv.gtf", col_names = FALSE) %>% 
  dplyr::distinct(X1, .keep_all = TRUE) %>%
  write_tsv("gtf/mm10.infercnv.dropdup.gtf",col_names = FALSE)



