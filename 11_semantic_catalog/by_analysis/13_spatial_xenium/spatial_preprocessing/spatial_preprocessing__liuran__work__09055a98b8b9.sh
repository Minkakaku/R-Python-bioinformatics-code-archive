#!/bin/bash

# 激活conda环境
source /home/liuran/software/miniconda3/etc/profile.d/conda.sh
conda activate rpy

# 运行R脚本
R --vanilla < monocle_2.R > monocle_2.R.log 2>&1
