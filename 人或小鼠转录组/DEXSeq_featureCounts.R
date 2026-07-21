log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# import library
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(DEXSeq)
library(stringr)
library(BiocParallel)
library(glue)

# input
bamDir <- dirname(dirname(snakemake@input[[1]]))
gtf <- snakemake@params[["gtf"]]
contrast <- snakemake@params[["contrast"]]
meta <- snakemake@params[["meta"]]
padj_cutoff <- as.double(snakemake@params[["padj_cutoff"]])

# output
dexseq_results <- snakemake@output[[1]]
dexseq_plots_zip <- snakemake@output[[2]]

# count reads
sample_df <- read.table(meta, header=TRUE, sep="\t")
sample_df <- sample_df[sample_df$group %in% contrast,]
samples <- sample_df$sample
groups <- sample_df$group


if (snakemake@params[["replicate"]] == "yes") {
  # set up parallelization
  BPPARAM <- MulticoreParam(snakemake@threads)
  
  # prepare the annotation
  txdb <- makeTxDbFromGFF(gtf)
  flattenedAnnotation <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
  names(flattenedAnnotation) <- sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
  
  bamFiles <- c()
  for (s in samples) {
    bam <- file.path(bamDir, s, paste0(s, ".sort.bam"))
    bamFiles <- c(bamFiles, bam)
  }
  
  bamFiles <- BamFileList(bamFiles)
  se <- summarizeOverlaps(
    flattenedAnnotation, 
    BamFileList(bamFiles), 
    singleEnd=FALSE,
    fragments=TRUE, 
    ignore.strand=TRUE, 
    BPPARAM=BPPARAM
  )
  
  # build a DEXSeqDataSet
  rownames(colData(se)) <- samples
  colData(se)$condition <- factor(groups, levels=replace(contrast,c(1,2),contrast[c(2,1)]))
  
  dxd <- DEXSeqDataSetFromSE(se, 
                             design= ~ sample + exon + condition:exon )
  
  # analyze DEU using DEXSeq
  ## normalisation
  dxd <- estimateSizeFactors(dxd)
  
  ## dispersion estimation
  dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)
  
  ## testing for differential exon usage
  dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
  dxr1 <- DEXSeqResults(dxd)
  
  ## output
  dxr1.sig <- as.data.frame(dxr1[dxr1$padj < padj_cutoff & !is.na(dxr1$padj),])
  dxr1.sig <- dxr1.sig %>%
    dplyr::rename(geneID = groupID) %>%
    dplyr::arrange(padj)
  dxr1.sig <- data.frame("groupID" = rownames(dxr1.sig),
                         dxr1.sig)
  
  trans <- dxr1.sig$transcripts
  transcripts <- c()
  for (lst in trans) {
    lst2str <- paste(unlist(lst, use.names=FALSE), collapse=',')
    transcripts <- c(transcripts, lst2str)
  }
  dxr1.sig$transcripts <- transcripts
  
  write.table(dxr1.sig, file=dexseq_results, row.names=FALSE, sep='\t')
  
  # plotDEXSeq
  dexseq_plots_dir <- gsub(".zip", "", dexseq_plots_zip)
  system(paste("mkdir", dexseq_plots_dir))
  
  work_dir <- getwd()
  setwd(dexseq_plots_dir)
  for (gene in unique(dxr1.sig$geneID)[1:50]) {
    png(paste0(basename(dirname(dexseq_plots_zip)),
               "_",
               gene,
               ".png"),
        width = 10,
        height = 8,
        units = "in",
        res = 300)
    plotDEXSeq(dxr1, gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, FDR=padj_cutoff)
    dev.off()
    pdf(paste0(basename(dirname(dexseq_plots_zip)),
               "_",
               gene,
               ".pdf"),
        width = 10,
        height = 8)
    plotDEXSeq(dxr1, gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, FDR=padj_cutoff)
    dev.off()
  }
  zip_file <-  paste0("../", basename(dexseq_plots_zip))
  system(paste("zip", zip_file, "./*"), intern=TRUE)
  setwd(work_dir)
  system(paste("rm -rf", dexseq_plots_dir), intern=TRUE)
} else if (snakemake@params[["replicate"]] == "no") {
  # if any group in this contrast does not have sample replicates, just generate empty result files
  
  # generate empty files
  control_group <- contrast[2]
  case_group <- contrast[1]
  log2fold <- paste("log2fold", case_group, control_group, sep="_")
  countdata <- paste(paste0("countData.", samples), collapse="\t")
  #headers <- glue("groupID\\tgeneID\\tfeatureID\\texonBaseMean\\tdispersion\\tstat\\tpvalue\\tpadj\\t{control_group}\\t{case_group}\\t{log2fold}\\tgenomicData.seqnames\\tgenomicData.start\\tgenomicData.end\\tgenomicData.width\\tgenomicData.strand\\t{countdata}\\ttranscripts\\n")
  headers <- glue("groupID\tgeneID\tfeatureID\texonBaseMean\tdispersion\tstat\tpvalue\tpadj\t{control_group}\t{case_group}\t{log2fold}\tgenomicData.seqnames\tgenomicData.start\tgenomicData.end\tgenomicData.width\tgenomicData.strand\t{countdata}\ttranscripts\n")
  #cmd1 <- paste0("echo -e '", headers, "' > ", dexseq_results)
  #system(cmd1, intern=TRUE)
  fileConn<-file(dexseq_results)
  writeLines(headers, fileConn)
  close(fileConn)
  
  cmd2 <- glue("echo | zip -q > {dexseq_plots_zip} && zip -dq {dexseq_plots_zip} -")
  system(cmd2, intern=TRUE)
}


