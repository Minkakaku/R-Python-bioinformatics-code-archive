library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
# library(ReactomePA)
library(DOSE)
library(ComplexUpset)
library(UpSetR)
library(tidyr)
library(ggsci)
library(dplyr)

GSEA_enrichment <- function(data,orgdb="org.Mm.eg.db",keggdb = "mmu",outdir){
  mydata <- data.frame(Gene = data$names, logFC = data$logfoldchanges) %>% arrange(desc(logFC))
  genelist <- mydata
    # bitr(mydata$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb) %>% 
    # inner_join(mydata, by = c("SYMBOL" = "Gene")) %>% arrange(desc(logFC))
  
  gsea_input <- genelist$logFC
  names(gsea_input) <- genelist$Gene
  head(gsea_input)
  
  geneset <- read.gmt("/mnt/NexData_file/home/zwl/Project/FWSC20253010/mmu.pathway.gmt")
  print(paste0('DEGgene in gmt number:',length(names(gsea_input)[names(gsea_input)%in% geneset$gene])))
  gg <- GSEA(gsea_input, TERM2GENE = geneset, verbose = F, pvalueCutoff = 0.8, pAdjustMethod = "BH")
  write.csv(gg@result, file = paste0(outdir, "/gsea_pathway.csv"), row.names = FALSE, quote = TRUE)
  mycol <- pal_nejm()(8)
  if (dim(gg@result)[1] != 0) {
    p = gseaplot2(gg, geneSetID = 1, color = mycol[1],
              title = "Specific KEGG",
              pvalue_table = F)
    ggsave(filename = paste0(outdir, "gsea_pathway.pdf"), width = 10, height = 7,p)
    ggsave(filename = paste0(outdir, "gsea_pathway.pdf"), width = 10, height = 7,p)
  }else{
    print('No gene enrichment')
  }
  
}


# setwd('/mnt/NexData_file/home/zwl/Project/FWSC20253010/require3/')
# input = '/mnt/NexData_file/home/zwl/Project/FWSC20253010/require3/cluster'
# output = '/mnt/NexData_file/home/zwl/Project/FWSC20253010/require3/cluster_GSEA/'
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20253010/require4/')
input = '/mnt/NexData_file/home/zwl/Project/FWSC20253010/require4/cell_type'
output = '/mnt/NexData_file/home/zwl/Project/FWSC20253010/require4/cell_type_GSEA/'
files <- list.files(input, pattern = "differentially_expressed_genes.xls", recursive = TRUE, full.names = TRUE)

for(n in files){
  cluster_name <- basename(dirname(n))
  print(cluster_name)
  setwd(output)
  dir.create(cluster_name)
  outdir = paste0(output,cluster_name,'/')
  
  data = read.table(file = n,sep = '\t',header = T)

  GSEA_enrichment(data,orgdb="org.Mm.eg.db",keggdb = "mmu",outdir)
}