log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(edgeR)
library(dplyr)
library(AnnotationHub)
library(ensembldb)

##############################
# Samples from each contrast #
##############################
params_lists = list(
  level1 = c("padj",2),
  level2 = c("padj",1.5),
  level3 = c("pvalue",2),
  level4 = c("pvalue",1.5),
  level5 = c("pvalue",1.2)
)

select_best_params <- function(res) {
  names = c()
  combs = c()
  ups = c()
  downs = c()
  for (params in c(1:length(params_lists))) {
    param_name = names(params_lists[params])
    param_comb = paste0(params_lists[[params]][1],", FC=", params_lists[[params]][2])
    res_new <- na.omit(data.frame( 
                  res[c(params_lists[[params]][1],"log2FoldChange")],
                  check.names=FALSE))
    up_num = dim(res_new[res_new[1] <= 0.05 & res_new[2] >= log2(as.numeric(params_lists[[params]][2])),])[1]
    down_num = dim(res_new[res_new[1] <= 0.05 & res_new[2] <= -log2(as.numeric(params_lists[[params]][2])),])[1]
    names = c(names,param_name)
    combs = c(combs,param_comb)
    ups = c(ups,up_num)
    downs = c(downs,down_num)
  }
  all_params = data.frame("param" = names, "comb" = combs, "up" = ups, "down" = downs)
  if (dim(all_params[all_params$up>=80 & all_params$down>=80,])[1] > 0) {
    final_level = all_params[all_params$up>=80 & all_params$down>=80,]$param[1]
  } else {
    final_level = "level5"
  }
  return(final_level)
}

all_level = c()

for (f in snakemake@input[['edger_tmps']]) {
  res <- read.table(file=f, header=TRUE, check.names=FALSE)
  # sort by p-value
  res <- res[order(res$pvalue),]

  level_select <- select_best_params(res)
  all_level <- c(all_level, level_select)
  # # Write dds object as RDS
  # print(snakemake@output[["rds_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))])
  # saveRDS(dds, file=snakemake@output[["rds_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))])
  # # Write normalized counts
  # norm_counts <- counts(dds, normalized=T)
  # print(snakemake@output[["tsv_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))])
  # write.table(data.frame("gene"=rownames(norm_counts), norm_counts, check.names=FALSE), 
  #             file=snakemake@output[["tsv_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))], 
  #             sep='\t', 
  #             row.names=F)
}

level_use <- rev(sort(unique(all_level)))[1]
data_params <- data.frame("p_type" = params_lists[[level_use]][1], "fc" = params_lists[[level_use]][2])

write.table(data_params, file=snakemake@output[["params_table"]], row.names=FALSE, sep='\t')

