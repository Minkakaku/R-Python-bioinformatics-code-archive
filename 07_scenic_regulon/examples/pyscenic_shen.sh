## 接下来的运行转到bash脚本里运行
#!/bin/bash

# 设置数据库目录路径
dir=/home/hf/project/shen/pyscenic/data
tfs=$dir/allTFs_hg38.txt
feather=$dir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# 确保数据库文件完整无误
ls $tfs $feather $tbl

######################
# 运行pySCENIC for CD4
######################
input_loom=$dir/mf.loom

# 2.1 GRN推断（推测基因调控网络）
pyscenic grn \
--num_workers 12 \
--output adj.mf.tsv \
--method grnboost2 \
$input_loom \
$tfs


# 2.2 cisTarget（使用Motif和Rankings进行基因调控网络修饰）
pyscenic ctx \
adj.mf.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output reg.mf.csv \
--num_workers 12
# pyscenic ctx \


# 2.3 AUCell（计算每个细胞中的转录因子活性）
pyscenic aucell \
$input_loom \
reg.mf.csv \
--output out_mf_SCENIC.loom