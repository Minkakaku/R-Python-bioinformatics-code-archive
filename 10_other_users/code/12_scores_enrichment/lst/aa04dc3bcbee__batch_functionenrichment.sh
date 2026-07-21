#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <base_dir>"
  echo "Example:"
  echo "  $0 /mnt/NexData_file/home/lst/project/FWSC20240744-02/售后10/neutrophils/DEG"
  exit 1
fi

base_dir="$1"
rscript_path="/mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R"

# 遍历DEG目录下所有一级子目录
for comp_dir in "$base_dir"/*/ ; do
  comp_name=$(basename "$comp_dir")  # 比如 PBS-vs-FCP-D14
  
  # 删除comp_dir下所有xls文件
  find "$comp_dir" -maxdepth 1 -type f -name "*.xls" -exec rm -f {} \;

  # 遍历comp_dir中所有以cluster开头的子目录
  for cluster_dir in "$comp_dir"/cluster*/ ; do
    [ -d "$cluster_dir" ] || continue  # 如果不存在，跳过

    cluster_name=$(basename "$cluster_dir")  # 例如 clusterXxx
    # 去除cluster前缀，加上_up后缀
    new_name="${cluster_name#cluster}_up"
    new_dir="$comp_dir/$new_name"

    # 重命名文件夹
    mv "$cluster_dir" "$new_dir"

    # 删除该文件夹下除differentially_expressed_genes.xls外的所有文件
    find "$new_dir" -maxdepth 1 -type f ! -name "differentially_expressed_genes.xls" -exec rm -f {} \;

    # 运行Rscript命令
    Rscript "$rscript_path" -f "$new_dir/differentially_expressed_genes.xls" -s mouse -o "$new_dir"
  done
done

