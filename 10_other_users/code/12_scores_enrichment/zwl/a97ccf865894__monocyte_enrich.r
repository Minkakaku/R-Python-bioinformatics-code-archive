
setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require5')
path = c('/mnt/NexData_file/home/zwl/Project/SC2511076348/require5/Myeloids2/celltype_enrich/clusterMono1/differentially_expressed_genes.xls',
         '/mnt/NexData_file/home/zwl/Project/SC2511076348/require5/Myeloids2/celltype_enrich/clusterMono2/differentially_expressed_genes.xls')
for (file_path in path) {
  output_dir <- dirname(normalizePath(file_path, winslash = "/", mustWork = TRUE))
  output_prefix <- file.path(output_dir, tools::file_path_sans_ext(basename(file_path)))
  
  data = read.table(file_path,sep = "\t",head=TRUE)
  print(dim(data))
  colnames(data) <- c("gene_name",'scores',"log2FC","pValue", "PValue","pct_nz_group" , "pct_nz_reference")
  #自行选择筛选自己感兴趣的基因集
  # data <- data[which(data$PValue <= 0.05 & data$log2FC >= 0.1 & data$pct_nz_group > 0.1),]
  # data <- data[which(data$PValue <= 0.05 & data$log2FC >= 0.1 & data$pct_nz_reference > 0.1),]
  print(dim(data))
  gene_set <- data$gene_name
  if(length(gene_set) > 0){
    GO_ALL <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
                       keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                       OrgDb=org.Hs.eg.db,  #物种对应的org包
                       ont = "ALL",   #CC细胞组件，MF分子功能，BP生物学过程
                       pvalueCutoff = 0.05,  #p值阈值
                       pAdjustMethod = "BH",  #多重假设检验校正方式
                       minGSSize = 3,   #注释的最小基因集，默认为10
                       maxGSSize = 500,  #注释的最大基因集，默认为500
                       qvalueCutoff = 1,  #q值阈值
                       readable = TRUE)  #基因ID转换为基因名
    write.table(GO_ALL, paste0(output_prefix, "_enrich.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)
    
    if(!is.null(GO_ALL) && nrow(GO_ALL) > 0){
      dotplot(GO_ALL,  #GO富集分析结果
              x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
              color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
              showCategory = 10,  #展示前20个点，默认为10个
              size = NULL,  #点的大小
              title = "GO_dotplot",  #设置图片的标题
              split='ONTOLOGY',
              label_format = 80)+
        facet_grid(ONTOLOGY~.,scale="free")
      
      ggsave(paste0(output_prefix, ".GO_dotplot.pdf"),width = 12,height = 9)
      ggsave(paste0(output_prefix, ".GO_dotplot.png"),width = 12,height = 9)
      
      barplot(GO_ALL,  #GO富集分析结果
              x = "Count",  #横坐标,默认Count,也可以为GeneRation
              color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
              showCategory = 10,  #展示前20个，默认为10个
              size = NULL,  
              title = "GO_barplot",  #设置图片的标题
              split='ONTOLOGY',
              label_format = 80)+
        facet_grid(ONTOLOGY~.,scale="free")
      
      ggsave(paste0(output_prefix, ".GO_barplot.pdf"),width = 12,height = 9)
      ggsave(paste0(output_prefix, ".GO_barplot.png"),width = 12,height = 9)
    }
    
    gene_symbol <- bitr(geneID = gene_set,  #感兴趣的基因集
                        fromType="SYMBOL",   #输入ID的类型
                        toType=c("ENTREZID"),   #输出ID的类型，可为多个
                        OrgDb="org.Hs.eg.db")  #物种注释数据库
    gene_KEGG <- gene_symbol[,2]
    
    KEGG<- enrichKEGG(gene = gene_KEGG,   #基因列表(同GO) 
                      organism = "hsa",  #物种
                      keyType = "kegg",  #指定的基因ID类型，默认为kegg
                      minGSSize = 3, 
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,  
                      pAdjustMethod = "BH",
                      qvalueCutoff = 1)
    
    id_map <- setNames(gene_symbol$SYMBOL, gene_symbol$ENTREZID)
    KEGG@result$geneID <- sapply(KEGG@result$geneID,
                                 function(x) paste(id_map[strsplit(x, "/")[[1]]], collapse = "/"))
    write.table(KEGG, paste0(output_prefix, "_enrich.KEGG.csv"), sep = ',', row.names = FALSE, quote = TRUE)
    
    KEGG@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", KEGG@result$Description)
    if(!is.null(KEGG) && nrow(KEGG) > 0){
      p <- dotplot(KEGG,  #GO富集分析结果
                   x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
                   color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
                   showCategory = 10,  #展示前20个点，默认为10个
                   size = NULL,  #点的大小
                   title = "KEGG_dotplot",  #设置图片的标题
                   label_format = 100)
      p + guides(colour = "Count")
      ggsave(paste0(output_prefix, ".KEGG_dotplot.pdf"),width = 12,height = 5)
      ggsave(paste0(output_prefix, ".KEGG_dotplot.png"),width = 12,height = 5)
      
      barplot(KEGG,  #GO富集分析结果
              x = "Count",  #横坐标,默认Count,也可以为GeneRation
              color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
              showCategory = 10,  #展示前20个，默认为10个
              size = NULL,  
              title = "KEGG_barplot",  #设置图片的标题
              label_format = 80)
      
      ggsave(paste0(output_prefix, ".KEGG_barplot.pdf"),width = 12,height = 5)
      ggsave(paste0(output_prefix, ".KEGG_barplot.png"),width = 12,height = 5)
    }
  }
}


