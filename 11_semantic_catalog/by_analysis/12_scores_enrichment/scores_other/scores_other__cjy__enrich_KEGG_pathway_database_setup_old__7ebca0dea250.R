log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type=c("message"))

# 在 R 脚本开头添加
user_lib <- "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Pipeline/rnaseq_pipeline_hsa_mmu"  # 替换为容器内可写路径，如 /home/user/Rlibs
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

if (snakemake@params[["species"]] == "homo_sapiens") {
  # install.packages(snakemake@params[["kegg_hsa_pathway_database"]], repos=NULL, type = "source", quiet=TRUE)
  library(KEGG.hsa.db)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  install.packages(snakemake@params[["kegg_mmu_pathway_database"]], repos=NULL, type = "source", quiet=TRUE)
  library(KEGG.mmu.db)
}  

system(paste("touch", snakemake@output[[1]]), intern=TRUE)
