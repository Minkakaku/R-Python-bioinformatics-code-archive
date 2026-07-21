log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(edgeR)
library(dplyr)
library(AnnotationHub)
library(ensembldb)

# local Ensembl database used for converting gene ID to gene name and gene annotation
ensdbfile <- snakemake@params[["ensdb"]]

# counts table
counts <- snakemake@input[[1]]
rawdata <- read.table(file = counts, sep = "\t", header = T, stringsAsFactors = F, check.names=FALSE)

metadata <- read.table(file = snakemake@params[["samples"]], sep = "\t", header = T)

case_group <- snakemake@params[["contrast"]][1]
control_group <- snakemake@params[["contrast"]][2]
case_samples <- metadata[metadata$group == case_group,]$sample
control_samples <- metadata[metadata$group == control_group,]$sample

samples <- c(case_samples, control_samples)
group <- c(rep(2, length(case_samples)), rep(1, length(control_samples)))

subdata <- rawdata[, c("Gene", samples)]

y <- DGEList(counts = subdata[, -1], genes = subdata[, 1], group = group)

## 2. Filtering
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]

## 3. Normalization
y <- calcNormFactors(y)

## 4. DEA
y_bcv <- y
# Note: hsa - 0.4; mmu - 0.1 
if (snakemake@params[["species"]] == "hsapiens") {
  bcv <- 0.4
} else if (snakemake@params[["species"]] == "mmusculus") {
  bcv <- 0.1
}  
et <- exactTest(y_bcv, dispersion = bcv ^ 2)

## 5. Output results
et_table <- et$table %>% 
  dplyr::select(-logCPM) %>%
  dplyr::rename(log2FoldChange = logFC, pvalue = PValue)

et_table$padj <- p.adjust(et_table$pvalue, "fdr")

genes <- y$genes
counts <- y$counts
counts <- cbind(counts, apply(counts[, case_samples, drop=FALSE], 1, mean))
counts <- cbind(counts, apply(counts[, control_samples, drop=FALSE], 1, mean))

# if any sample name is the same as group name, add "_name" to group name in the results table
case_group_new <- case_group
for (s in case_samples) {
  if (s == case_group) {
    case_group_new <- paste0(case_group, "_mean")
  }
}
control_group_new <- control_group
for (s in control_samples) {
  if (s == control_group) {
    control_group_new <- paste0(control_group, "_mean")
  }
}
colnames(counts)[-(1:(length(colnames(counts))-2))] <- c(case_group_new, control_group_new)

results <- cbind(genes, counts, et_table)
colnames(results)[1] <- "gene"

# convert gene id to gene name

edb <- EnsDb(ensdbfile)
g2g <- AnnotationDbi::select(
  edb,
  keys = results$gene,
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

annotated <- merge(results, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$external_gene_name <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated <- annotated %>%
  dplyr::rename(gene_name = external_gene_name)

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

system(paste("rm -f ", snakemake@input[["table_tmp"]]), intern=TRUE)

