log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(AnnotationHub)
library(ensembldb)
library(ashr)
library(dplyr)

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# local Ensembl database used for converting gene ID to gene name and gene annotation
ensdbfile <- snakemake@params[["ensdb"]]

# data from deseq2-init
dds <- readRDS(snakemake@input[[1]])

contrast <- c("group", snakemake@params[["contrast"]])
res <- results(dds, contrast=contrast, parallel=parallel)

if (snakemake@params[["ashr"]] == "yes") {
  # shrink fold changes for lowly expressed genes
  # use ashr so we can use `contrast` as conversion to coef is not trivial
  # see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
  res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")
}

# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage

# generate MA plot
png(snakemake@output[["ma_png"]])
plotMA(res, ylim=c(-2,2))
dev.off()
pdf(snakemake@output[["ma_pdf"]])
plotMA(res, ylim=c(-2,2))
dev.off()

# convert gene id to gene name

case_samples <- rownames(colData(dds)[colData(dds)$group %in% contrast[2],])
control_samples <- rownames(colData(dds)[colData(dds)$group %in% contrast[3],])
all_samples <- c(case_samples, control_samples)
counts <- counts(dds)[rownames(res), all_samples]
counts <- cbind(counts, apply(counts[,case_samples], 1, mean))
counts <- cbind(counts, apply(counts[,control_samples], 1, mean))

# if any sample name is the same as group name, add "_name" to group name in the results table
case_group_new <- contrast[2]
for (s in case_samples) {
  if (s == contrast[2]) {
    case_group_new <- paste0(contrast[2], "_mean")
  }
}
control_group_new <- contrast[3]
for (s in control_samples) {
  if (s == contrast[3]) {
    control_group_new <- paste0(contrast[3], "_mean")
  }
}
colnames(counts)[-(1:(length(colnames(counts))-2))] <- c(case_group_new, control_group_new)

#res <- data.frame("gene"=rownames(res),res)
res <- data.frame("gene"=rownames(res), 
                  counts[rownames(res),], 
                  res[,-1],
                  check.names=FALSE)

edb <- EnsDb(ensdbfile)
g2g <- AnnotationDbi::select(
  edb,
  keys = res$gene,
  keytype = "GENEID",
  columns = c(
    "GENEID",
    "GENENAME",
    "SEQNAME",
    "GENESEQSTART",
    "GENESEQEND",
    "SEQSTRAND",
    "GENEBIOTYPE",
    "DESCRIPTION"
  ),
  multiVals = "first"
)

# Since 'multiVals = "first"' parameter doesn't work in AnnotationDbi::select (ONLY works for AnnotationDbi::mapIds)
# Manully make sure 1:1 annotation mapping for each gene
g2g <- g2g[!duplicated(g2g$GENEID),]

colnames(g2g) <- c("ensembl_gene_id", "external_gene_name",
                   "gene_chr", "gene_start", "gene_end", "gene_strand",
                   "gene_biotype", "gene_description")
g2g$gene_strand <- ifelse(g2g$gene_strand == 1, "+", ifelse(g2g$gene_strand == -1, "-", ""))

annotated <- merge(res, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$external_gene_name <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated <- annotated %>%
  dplyr::rename(gene_name = external_gene_name) %>%
  dplyr::select(-lfcSE)

annotated <- annotated[order(annotated$padj, annotated$pvalue),]

write.table(annotated, file=snakemake@output[["table"]], row.names=FALSE, sep='\t')

# results table for DEG only
select_params <- read.table(snakemake@input[[2]], header=TRUE, check.names=FALSE)

p_type <- snakemake@params[["p_type"]]
FC_cutoff <- snakemake@params[["FC_cutoff"]]
if (p_type == "auto") {
  p_type <- select_params[1,"p_type"]
} else {
  p_type <- snakemake@params[["p_type"]]
}

if (FC_cutoff == "auto") {
  FC_cutoff <- as.numeric(select_params[1,"fc"])
} else {
  FC_cutoff <- as.numeric(snakemake@params[["FC_cutoff"]])
}

log2FC_cutoff <- log2(FC_cutoff)
p_cutoff <- snakemake@params[["p_cutoff"]]
if (p_type == "padj") {
  annotated_deg <- annotated
  annotated_deg$regulate  <- as.factor(ifelse(is.na(annotated$padj), 'Normal', 
                                              ifelse(annotated$padj < p_cutoff & annotated$log2FoldChange > log2FC_cutoff, 'Up',
                                                     ifelse(annotated$padj < p_cutoff & annotated$log2FoldChange < -log2FC_cutoff, 'Down', 'Normal'))))
  annotated_deg <- annotated_deg[annotated_deg$regulate %in% c("Up", "Down"),]
} else if (p_type == "pvalue") {
  annotated_deg <- annotated
  annotated_deg$regulate  <- as.factor(ifelse(is.na(annotated$pvalue), 'Normal', 
                                              ifelse(annotated$pvalue < p_cutoff & annotated$log2FoldChange > log2FC_cutoff, 'Up',
                                                     ifelse(annotated$pvalue < p_cutoff & annotated$log2FoldChange < -log2FC_cutoff, 'Down', 'Normal'))))
  annotated_deg <- annotated_deg[annotated_deg$regulate %in% c("Up", "Down"),]
}

write.table(annotated_deg, file=snakemake@output[["table_deg"]], row.names=FALSE, sep='\t')
