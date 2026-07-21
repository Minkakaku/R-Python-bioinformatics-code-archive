log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

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

for (contrast in names(snakemake@params[["contrasts"]])) {
  
  print(contrast)
  contrast_groups <- unlist(strsplit(contrast, "-vs-"))
  
  # colData and countData must have the same sample order, but this is ensured
  # by the way we create the count matrix
  coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)
  coldata <- coldata[coldata$group %in% contrast_groups, , drop=FALSE]
  coldata <- coldata[order(row.names(coldata)), , drop=F]
  
  cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="Gene", check.names=FALSE)
  cts <- cts[ , rownames(coldata)]
  cts <- cts[ , order(names(cts))]
  
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData=coldata,
                                design=as.formula(snakemake@params[["model"]]))
  
  # remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  # normalization and preprocessing
  dds <- DESeq(dds, parallel=parallel)
  
  contrast <- c("group", contrast_groups)
  res <- results(dds, contrast=contrast, parallel=parallel)

  if (snakemake@params[["ashr"]] == "yes") {
    # shrink fold changes for lowly expressed genes
    # use ashr so we can use `contrast` as conversion to coef is not trivial
    # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
    res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")
  }

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



