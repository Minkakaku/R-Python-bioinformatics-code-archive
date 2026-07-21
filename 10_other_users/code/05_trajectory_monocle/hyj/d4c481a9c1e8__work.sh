#!/bin/bash

# 激活conda环境
source /home/liuran/software/miniconda3/etc/profile.d/conda.sh
conda activate rpy

# 运行R脚本
R --vanilla < monocle.R > monocle.R.log 2>&1
