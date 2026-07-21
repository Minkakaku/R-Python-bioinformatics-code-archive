setwd('/mnt/NexData_file/home/zwl/Project/SC2512157505/require2')
rm(list=ls())
gc()

library(clusterProfiler)
library(enrichplot)
library(myenrichplot)

gmt <- read.gmt("/mnt/NexData_file/home/zwl/Project/SC2512157505/require2/pathway.gmt")
gmt$gene <- stringr::str_to_title(gmt$gene)
head(gmt)

DEG_path = '/mnt/NexData_file/home/zwl/Project/SC2512157505/require1/Neutrophils/celltype/groupDEG/cluster'
folders = list.dirs(DEG_path, full.names = F, recursive = FALSE)

DEG_path = '/mnt/NexData_file/home/zwl/Project/SC2512157505/require1/Neutrophils/celltype/groupDEG/'
folders = 'ALL'
for (celltype in folders) {
  folders_g = list.dirs(paste0(DEG_path,'/',celltype), full.names = F, recursive = FALSE)
  for (group in folders_g) {
    pathway = paste0(DEG_path,'/',celltype,'/',group,'/0.',group,'_All_expressed_genes.xls')
    data <- read.table(pathway, header = TRUE, sep = "\t", fileEncoding = "UTF-8")

    data_all_sort <- data %>% arrange(desc(logfoldchanges))
    geneList <- data_all_sort$logfoldchanges
    names(geneList) <- data_all_sort$names
    # names(geneList) %in% gmt$gene
    
    res <- GSEA(geneList,TERM2GENE = gmt,  minGSSize = 1,
                pvalueCutoff = 0.9,pAdjustMethod = 'BH',
                maxGSSize = 500)
    write.csv(file = paste0(DEG_path,'/',celltype,'/',group,'/0.ALLpathway_GSEA_res.csv'),res@result)
    num = length(res@result$ID)
    print(num)
    if (num ==0) {
      print('No plot')
    }else{
      for (i in 1:num) {
        p = myenrichplot::gseaplot2(res, geneSetID = i,
                      pvalue_table = T,color = 'green2',
                      title = res$Description[i],
                      # pvalue_table_rownames = NULL,
                      # pvalue_table_columns = c("NES", "p.adjust")
                      )
        p
        ggsave(filename = paste0(DEG_path,'/',celltype,'/',group,'/0.',res@result[["ID"]][i],'_GSEA_plot.png'),
               p,w=8,h=6,dpi = 200)
        ggsave(filename = paste0(DEG_path,'/',celltype,'/',group,'/0.',res@result[["ID"]][i],'_GSEA_plot.pdf'),
               p,w=8,h=6,dpi = 200)
      }
    }
  }
}
