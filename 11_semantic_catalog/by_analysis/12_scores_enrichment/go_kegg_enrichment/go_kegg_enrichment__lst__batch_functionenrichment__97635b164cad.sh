#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <base_dir>"
  echo "Example:"
  echo "  $0 \"/mnt/NexData_file/home/lst/project/FWSC20240744-02/售后10/neutrophils/DEG\""
  exit 1
fi

base_dir="$1"
rscript_path="/mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R"

echo "开始处理目录：$base_dir"

for comp_dir in "$base_dir"/*/ ; do
  comp_name=$(basename "$comp_dir")
  echo "处理比较组目录: $comp_name"

  # 删除一级目录下所有xls文件
  echo "删除目录 $comp_dir 下所有一级xls文件..."
  find "$comp_dir" -maxdepth 1 -type f -name "*.xls" -exec rm -f {} \;

  # 遍历以cluster开头的子目录
  found_cluster=0
  for cluster_dir in "$comp_dir"/cluster*/ ; do
    if [ ! -d "$cluster_dir" ]; then
      echo "未找到以cluster开头的子目录，跳过"
      break
    fi

    found_cluster=1
    echo "检测到cluster目录: $cluster_dir"

    cluster_name=$(basename "$cluster_dir")
    new_name="${cluster_name#cluster}_up"
    new_dir="$comp_dir/$new_name"

    echo "准备重命名: $cluster_name -> $new_name"
    mv "$cluster_dir" "$new_dir"

    echo "删除 $new_dir 下除differentially_expressed_genes.xls外所有文件..."
    find "$new_dir" -maxdepth 1 -type f ! -name "differentially_expressed_genes.xls" -exec rm -f {} \;

    if [ -f "$new_dir/differentially_expressed_genes.xls" ]; then
      echo "找到差异表达文件: $new_dir/differentially_expressed_genes.xls"
      echo "开始调用Rscript进行富集分析..."
      Rscript "$rscript_path" -f "$new_dir/differentially_expressed_genes.xls" -s mouse -o "$new_dir"
      if [ $? -eq 0 ]; then
        echo "Rscript运行成功"
      else
        echo "Rscript运行失败，请检查R脚本及参数"
      fi
    else
      echo "警告: $new_dir/differentially_expressed_genes.xls 文件不存在，跳过Rscript调用"
    fi

  done

  if [ "$found_cluster" -eq 0 ]; then
    echo "该比较组目录无以cluster开头的子目录"
  fi

done

echo "全部处理完成"
