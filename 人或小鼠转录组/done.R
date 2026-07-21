# log <- file(snakemake@log[[1]], open="wt")
# sink(log)
# sink(log, type="message")

system(paste0("rm -rf ", snakemake@params[["output_path"]], "/flag"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/logs"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/tmp"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/01.data_clean/cleandata"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/01.data_clean/fastp/*.clean.fastq.gz"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/02.alignment_count/hisat2/*/*.sort.bam*"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/04.diffexp_analysis/Recommend.params.txt"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/04.diffexp_analysis/DESeq2_norm_counts_pca*"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/04.diffexp_analysis/*/DESeq2_norm_counts_pca*"), intern=TRUE)
system(paste0("rm -rf ", snakemake@params[["output_path"]], "/complete.done"), intern=TRUE)

system(paste0("echo 'Workflow completed successfully.' > ",  snakemake@output[["donefile"]]), intern=TRUE)

