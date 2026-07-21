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

###########################
# Samples from all groups #
###########################

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="Gene", check.names=FALSE)
cts <- cts[ , order(names(cts))]

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

# Write dds object as RDS
saveRDS(dds, file=snakemake@output[["rds_all"]])
# Write normalized counts
norm_counts <- counts(dds, normalized=T)
write.table(data.frame("gene"=rownames(norm_counts), norm_counts, check.names=FALSE), 
            file=snakemake@output[["tsv_all"]], 
            sep='\t', 
            row.names=F)

##############################
# Samples from each contrast #
##############################

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
  
  # Write dds object as RDS
  print(snakemake@output[["rds_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))])
  saveRDS(dds, file=snakemake@output[["rds_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))])
  # Write normalized counts
  norm_counts <- counts(dds, normalized=T)
  print(snakemake@output[["tsv_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))])
  write.table(data.frame("gene"=rownames(norm_counts), norm_counts, check.names=FALSE), 
              file=snakemake@output[["tsv_contrast"]][match(contrast, names(snakemake@params[["contrasts"]]))], 
              sep='\t', 
              row.names=F)
  
}


