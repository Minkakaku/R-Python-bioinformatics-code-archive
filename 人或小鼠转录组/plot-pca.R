log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
png(snakemake@output[["png"]], 
    width = 10, 
    height = 8, 
    units = "in",
    res = 300)
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
dev.off()
pdf(snakemake@output[["pdf"]], 
    width = 10, 
    height = 8)
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
dev.off()
