#!/usr/bin/bash
set -e

# Bash script to prepare AmpliSeqPipe to_customer files

# Command example:
#  1. Complete version output: $bash to_customer.sh output_demo_16s to_customer_demo_16s data/demo_16s config/metadata_demo_16s.tsv resources/readme.zip
#  2. Slim version output: $bash to_customer.sh output_demo_16s to_customer_demo_16s data/demo_16s config/metadata_demo_16s.tsv resources/readme.zip slim

# e.g. /home/lyw/pipeline/AmpliSeqPipe/output_demo_16S
old_dir=$(realpath ${1})
# e.g. /home/lyw/pipeline/AmpliSeqPipe/to_customer_demo_16S
new_dir=$(realpath ${2})
# e.g. metadata/metadata_demo_16S.tsv
metadata=$(realpath ${3})
# e.g. /home/lyw/pipeline/AmpliSeqPipe/resources/readme.zip
readme_zip=$(realpath ${4})
# (OPTIONAL) add one more parameter "slim" if output is slim version
slim=${5}

# rawdata (need to modify if the format is not correct)
mkdir -p "${new_dir}"/results/rawdata
{ cat "${metadata}" | sed 1d; printf '\n'; } | while IFS=$'\t' read -r sample group fq1 fq2
do
  [ -z "$sample" ] && continue
  fq1_path=$(realpath "$fq1")
  fq2_path=$(realpath "$fq2")
  cp -uv "${fq1_path}" "${new_dir}/results/rawdata/${sample}_1.raw.fastq.gz"
  cp -uv "${fq2_path}" "${new_dir}/results/rawdata/${sample}_2.raw.fastq.gz"
done

# cleandata
cp -ruv "${old_dir}"/01.data_clean/cleandata "${new_dir}"/results/

if [ "${slim}" == "slim" ]
then
  # 01.data_clean
  mkdir -p "${new_dir}"/results/01.data_clean
  cp -ruv "${old_dir}"/01.data_clean/multiqc "${new_dir}"/results/01.data_clean/
  cp -ruv "${old_dir}"/01.data_clean/qc_plot "${new_dir}"/results/01.data_clean/
  rm -rf "${new_dir}"/results/01.data_clean/qc_plot/*/*_GC.*
  cp -uv "${old_dir}"/02.ASVs_analysis/all_reads_qc_stats.xls "${new_dir}"/results/01.data_clean/
  cp -ruv "${old_dir}"/02.ASVs_analysis/reads_length_barplot "${new_dir}"/results/01.data_clean/

  # 02.ASVs_analysis
  mkdir -p "${new_dir}"/results/02.ASVs_analysis
  cp -uv "${old_dir}"/02.ASVs_analysis/ASVs_sequence.fasta "${new_dir}"/results/02.ASVs_analysis/
  cp -uv "${old_dir}"/02.ASVs_analysis/all_ASVs_count.xls "${new_dir}"/results/02.ASVs_analysis/
  cp -uv "${old_dir}"/02.ASVs_analysis/all_ASVs_count_rarefy.xls "${new_dir}"/results/02.ASVs_analysis/
  cp -uv "${old_dir}"/02.ASVs_analysis/all_ASVs_count_stats.xls "${new_dir}"/results/02.ASVs_analysis/
  cp -ruv "${old_dir}"/02.ASVs_analysis/ASVs_count_barplot "${new_dir}"/results/02.ASVs_analysis/
  cp -ruv "${old_dir}"/02.ASVs_analysis/ASVs_count_venn "${new_dir}"/results/02.ASVs_analysis/

  # 03.taxonomy
  mkdir -p "${new_dir}"/results/03.taxonomy
  cp -ruv "${old_dir}"/03.taxonomy/qiime2_classify_sklearn/ASVs_taxonomy.xls "${new_dir}"/results/03.taxonomy/
  cp -ruv "${old_dir}"/03.taxonomy/abundance "${new_dir}"/results/03.taxonomy/
  rm -rf "${new_dir}"/results/03.taxonomy/abundance/absolute/*.tmp.xls
  rm -rf "${new_dir}"/results/03.taxonomy/abundance/absolute/*.qza
  cp -ruv "${old_dir}"/03.taxonomy/barplot "${new_dir}"/results/03.taxonomy/
  cp -ruv "${old_dir}"/03.taxonomy/krona "${new_dir}"/results/03.taxonomy/

  # 04.alpha_diversity
  mkdir -p "${new_dir}"/results/04.alpha_diversity
  cp -uv "${old_dir}"/04.alpha_diversity/all.alpha_index.xls "${new_dir}"/results/04.alpha_diversity/
  cp -ruv "${old_dir}"/04.alpha_diversity/alpha_rarefaction_curve "${new_dir}"/results/04.alpha_diversity/

  # zip all results folders
  cd "${new_dir}"/results || exit 1
  zip -r 01.data_clean.zip 01.data_clean && rm -rf 01.data_clean
  zip -r 02.ASVs_analysis.zip 02.ASVs_analysis && rm -rf 02.ASVs_analysis
  zip -r 03.taxonomy.zip 03.taxonomy && rm -rf 03.taxonomy
  zip -r 04.alpha_diversity.zip 04.alpha_diversity && rm -rf 04.alpha_diversity
else
  # 01.data_clean
  mkdir -p "${new_dir}"/results/01.data_clean
  cp -ruv "${old_dir}"/01.data_clean/multiqc "${new_dir}"/results/01.data_clean/
  cp -ruv "${old_dir}"/01.data_clean/qc_plot "${new_dir}"/results/01.data_clean/
  rm -rf "${new_dir}"/results/01.data_clean/qc_plot/*/*_GC.*
  cp -uv "${old_dir}"/02.ASVs_analysis/all_reads_qc_stats.xls "${new_dir}"/results/01.data_clean/
  cp -ruv "${old_dir}"/02.ASVs_analysis/reads_length_barplot "${new_dir}"/results/01.data_clean/

  # 02.ASVs_analysis
  mkdir -p "${new_dir}"/results/02.ASVs_analysis
  cp -uv "${old_dir}"/02.ASVs_analysis/ASVs_sequence.fasta "${new_dir}"/results/02.ASVs_analysis/
  cp -uv "${old_dir}"/02.ASVs_analysis/all_ASVs_count.xls "${new_dir}"/results/02.ASVs_analysis/
  cp -uv "${old_dir}"/02.ASVs_analysis/all_ASVs_count_rarefy.xls "${new_dir}"/results/02.ASVs_analysis/
  cp -uv "${old_dir}"/02.ASVs_analysis/all_ASVs_count_stats.xls "${new_dir}"/results/02.ASVs_analysis/
  cp -ruv "${old_dir}"/02.ASVs_analysis/ASVs_count_barplot "${new_dir}"/results/02.ASVs_analysis/
  cp -ruv "${old_dir}"/02.ASVs_analysis/ASVs_count_venn "${new_dir}"/results/02.ASVs_analysis/

  # 03.taxonomy
  mkdir -p "${new_dir}"/results/03.taxonomy
  cp -ruv "${old_dir}"/03.taxonomy/qiime2_classify_sklearn/ASVs_taxonomy.xls "${new_dir}"/results/03.taxonomy/
  cp -ruv "${old_dir}"/03.taxonomy/abundance "${new_dir}"/results/03.taxonomy/
  rm -rf "${new_dir}"/results/03.taxonomy/abundance/absolute/*.tmp.xls
  rm -rf "${new_dir}"/results/03.taxonomy/abundance/absolute/*.qza
  cp -ruv "${old_dir}"/03.taxonomy/barplot "${new_dir}"/results/03.taxonomy/
  cp -ruv "${old_dir}"/03.taxonomy/heatmap "${new_dir}"/results/03.taxonomy/
  cp -ruv "${old_dir}"/03.taxonomy/krona "${new_dir}"/results/03.taxonomy/

  # 04.alpha_diversity
  mkdir -p "${new_dir}"/results/04.alpha_diversity
  cp -uv "${old_dir}"/04.alpha_diversity/all.alpha_index.xls "${new_dir}"/results/04.alpha_diversity/
  cp -ruv "${old_dir}"/04.alpha_diversity/alpha_stats "${new_dir}"/results/04.alpha_diversity/
  cp -ruv "${old_dir}"/04.alpha_diversity/alpha_rarefaction_curve "${new_dir}"/results/04.alpha_diversity/

  # 05.beta_diversity
  mkdir -p "${new_dir}"/results/05.beta_diversity
  cp -ruv "${old_dir}"/05.beta_diversity/beta_matrix "${new_dir}"/results/05.beta_diversity/
  cp -ruv "${old_dir}"/05.beta_diversity/beta_PCA "${new_dir}"/results/05.beta_diversity/
  cp -ruv "${old_dir}"/05.beta_diversity/beta_PCoA "${new_dir}"/results/05.beta_diversity/
  cp -ruv "${old_dir}"/05.beta_diversity/beta_NMDS "${new_dir}"/results/05.beta_diversity/
  cp -ruv "${old_dir}"/05.beta_diversity/beta_UPGMA/UPGMA_tree_barplot "${new_dir}"/results/05.beta_diversity/beta_UPGMA
  cp -ruv "${old_dir}"/05.beta_diversity/beta_heatmap "${new_dir}"/results/05.beta_diversity/
  cp -ruv "${old_dir}"/05.beta_diversity/beta_anosim "${new_dir}"/results/05.beta_diversity/

  # 06.taxonomic_biomarker
  cp -ruv "${old_dir}"/06.taxonomic_biomarker "${new_dir}"/results/

  # 07.association_analysis
  cp -ruv "${old_dir}"/07.association_analysis "${new_dir}"/results/

  # 08.function_analysis
  cp -ruv "${old_dir}"/08.function_analysis "${new_dir}"/results/
  if [ -d "${new_dir}"/results/08.function_analysis/qiime2_picrust2 ]
  then
    rm -rf "${new_dir}"/results/08.function_analysis/qiime2_picrust2
  elif [ -d "${new_dir}"/results/08.function_analysis/FUNGuild ]
  then
    rm -rf "${new_dir}"/results/08.function_analysis/FUNGuild/FUNGuild_input*
  fi

  # zip all results folders
  cd "${new_dir}"/results || exit 1
  zip -r 01.data_clean.zip 01.data_clean && rm -rf 01.data_clean
  zip -r 02.ASVs_analysis.zip 02.ASVs_analysis && rm -rf 02.ASVs_analysis
  zip -r 03.taxonomy.zip 03.taxonomy && rm -rf 03.taxonomy
  zip -r 04.alpha_diversity.zip 04.alpha_diversity && rm -rf 04.alpha_diversity
  zip -r 05.beta_diversity.zip 05.beta_diversity && rm -rf 05.beta_diversity
  zip -r 06.taxonomic_biomarker.zip 06.taxonomic_biomarker && rm -rf 06.taxonomic_biomarker
  zip -r 07.association_analysis.zip 07.association_analysis && rm -rf 07.association_analysis
  zip -r 08.function_analysis.zip 08.function_analysis && rm -rf 08.function_analysis
fi

# report
cp -uv "${old_dir}"/09.summary_report/summary_report.html "${new_dir}"/report.html

# other files
unzip -u -d "${new_dir}"/ ${readme_zip}
mv "${new_dir}"/readme/MD5Summer.exe "${new_dir}"/results/
mv "${new_dir}"/readme/广科安德数据校验方法.pdf "${new_dir}"/results/
mv "${new_dir}"/readme/结果文件说明文档.pdf "${new_dir}"/results/
mv "${new_dir}"/readme/软件和数据库版本.pdf "${new_dir}"/results/
rm -rf "${new_dir}"/readme

# generate md5
cd "${new_dir}"/results || exit 1
find . -type f -print0 | xargs -0 md5sum | sed -e '/\(md5.txt\|MD5Summer.exe\|report.html\|广科安德数据校验方法.pdf\|结果文件说明文档.pdf\|软件和数据库版本.pdf\)$/d' > md5.txt
