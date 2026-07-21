rm(list=ls())
gc()
gopathway = function(file_path,orgdb="org.Mm.eg.db",keggdb = "mmu",trend='ALL'){
  if (file.exists(file_path)) {
    output_dir <- dirname(normalizePath(file_path, winslash = "/", mustWork = TRUE))
    output_prefix <- file.path(output_dir, tools::file_path_sans_ext(basename(file_path)))
    
    data1 = read.table(file_path,sep = "\t",head=TRUE)
    print(dim(data1))
    colnames(data1) <- c("gene_name", "scores", "log2FC","pvals", "PValue","pct_nz_group" , "pct_nz_reference", "cluster")
    #自行选择筛选自己感兴趣的基因集
    if (trend =='up') {
      data <- data1[which(data1$PValue <= 0.05 & data1$log2FC >= 0.5 & data1$pct_nz_reference > 0.1),]
    }else if(trend =='down'){
      data <- data1[which(data1$PValue <= 0.05 & data1$log2FC <= -0.5 & data1$pct_nz_reference > 0.1),]
    }else{
      print('不分上下调')
      data = data1[which(data1$PValue <= 0.05 & (abs(data1$log2FC) > 0.5) & data1$pct_nz_reference > 0.1),]
    }
    
    print(dim(data))
    gene_set <- data$gene_name
    if(length(gene_set) > 0){
      GO_BP <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
                         keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                         OrgDb=orgdb,  #物种对应的org包
                         ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
                         pvalueCutoff = 0.05,  #p值阈值
                         pAdjustMethod = "BH",  #多重假设检验校正方式
                         minGSSize = 3,   #注释的最小基因集，默认为10
                         maxGSSize = 500,  #注释的最大基因集，默认为500
                         qvalueCutoff = 1,  #q值阈值
                         readable = TRUE)  #基因ID转换为基因名
      write.table(GO_BP, paste0(output_prefix,'_',trend, "_enrichBP.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)
      
      if(!is.null(GO_BP) && nrow(GO_BP) > 0){

        barplot(GO_BP,  #GO富集分析结果
                x = "Count",  #横坐标,默认Count,也可以为GeneRation
                color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
                showCategory = 10,  #展示前20个，默认为10个
                size = NULL,  
                title = "Barplot for Biological Process",  #设置图片的标题
                label_format = 80)
        
        ggsave( paste0(output_prefix,'_',trend, "_enrichBP.go_barplot.pdf"),width = 10,height = 6)
        ggsave( paste0(output_prefix,'_',trend, "_enrichBP.go_barplot.png"),width = 10,height = 6)
      }
      
      GO_CC <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
                        keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                        OrgDb=orgdb,  #物种对应的org包
                        ont = "CC",   #CC细胞组件，MF分子功能，BP生物学过程
                        pvalueCutoff = 0.05,  #p值阈值
                        pAdjustMethod = "BH",  #多重假设检验校正方式
                        minGSSize = 3,   #注释的最小基因集，默认为10
                        maxGSSize = 500,  #注释的最大基因集，默认为500
                        qvalueCutoff = 1,  #q值阈值
                        readable = TRUE)  #基因ID转换为基因名
      write.table(GO_CC, paste0(output_prefix,'_',trend, "_enrichCC.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)
      
      if(!is.null(GO_CC) && nrow(GO_CC) > 0){
        
        barplot(GO_CC,  #GO富集分析结果
                x = "Count",  #横坐标,默认Count,也可以为GeneRation
                color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
                showCategory = 10,  #展示前20个，默认为10个
                size = NULL,  
                title = "Barplot for Cell Component",   #设置图片的标题
                label_format = 80)
        
        ggsave( paste0(output_prefix,'_',trend, "_enrichCC.go_barplot.pdf"),width = 10,height = 6)
        ggsave( paste0(output_prefix,'_',trend, "_enrichCC.go_barplot.png"),width = 10,height = 6)
      }
      
      GO_MF <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
                        keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
                        OrgDb=orgdb,  #物种对应的org包
                        ont = "MF",   #CC细胞组件，MF分子功能，BP生物学过程
                        pvalueCutoff = 0.05,  #p值阈值
                        pAdjustMethod = "BH",  #多重假设检验校正方式
                        minGSSize = 3,   #注释的最小基因集，默认为10
                        maxGSSize = 500,  #注释的最大基因集，默认为500
                        qvalueCutoff = 1,  #q值阈值
                        readable = TRUE)  #基因ID转换为基因名
      write.table(GO_MF, paste0(output_prefix,'_',trend, "_enrichMF.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)
      
      if(!is.null(GO_MF) && nrow(GO_MF) > 0){
        
        barplot(GO_MF,  #GO富集分析结果
                x = "Count",  #横坐标,默认Count,也可以为GeneRation
                color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
                showCategory = 10,  #展示前20个，默认为10个
                size = NULL,  
                title = "Barplot for Molecular Function",  #设置图片的标题
                label_format = 80)
        
        ggsave( paste0(output_prefix,'_',trend, "_enrichMF.go_barplot.pdf"),width = 10,height = 6)
        ggsave( paste0(output_prefix,'_',trend, "_enrichMF.go_barplot.png"),width = 10,height = 6)
      }
      
      gene_symbol <- bitr(geneID = gene_set,  #感兴趣的基因集
                          fromType="SYMBOL",   #输入ID的类型
                          toType=c("ENTREZID"),   #输出ID的类型，可为多个
                          OrgDb="org.Mm.eg.db")  #物种注释数据库
      gene_KEGG <- gene_symbol[,2]
      
      KEGG<- enrichKEGG(gene = gene_KEGG,   #基因列表(同GO) 
                        organism = keggdb,  #物种
                        keyType = "kegg",  #指定的基因ID类型，默认为kegg
                        minGSSize = 3, 
                        maxGSSize = 500,
                        pvalueCutoff = 0.05,  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 1)
      
      id_map <- setNames(gene_symbol$SYMBOL, gene_symbol$ENTREZID)
      KEGG@result$geneID <- sapply(KEGG@result$geneID,
                                   function(x) paste(id_map[strsplit(x, "/")[[1]]], collapse = "/"))
      write.table(KEGG, paste0(output_prefix,'_',trend, "_enrich.KEGG.csv"), sep = ',', row.names = FALSE, quote = TRUE)
      
      KEGG@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", KEGG@result$Description)
      if(!is.null(KEGG) && nrow(KEGG) > 0){
        
        barplot(KEGG,  #GO富集分析结果
                x = "Count",  #横坐标,默认Count,也可以为GeneRation
                color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
                showCategory = 10,  #展示前20个，默认为10个
                size = NULL,  
                title = "KEGG_barplot",  #设置图片的标题
                label_format = 80)
        
        ggsave(paste0(output_prefix,'_',trend,  ".KEGG_barplot.pdf"),width = 10,height = 6)
        ggsave(paste0(output_prefix,'_',trend,  ".KEGG_barplot.png"),width = 10,height = 6)
        
        data_ALL = gene_symbol %>% left_join(data,by=c('SYMBOL'='gene_name'))
        geneList = data_ALL$log2FC
        names(geneList) = data_ALL$SYMBOL
        edox <- setReadable(KEGG, 'org.Mm.eg.db', 'ENTREZID')
        p1 <- cnetplot(KEGG,colorEdge = TRUE,foldChange=geneList,showCategory = 4)
        ggsave(paste0(output_prefix,'_',trend,  ".KEGG_cnetplot.pdf"),width = 10,height = 8,p1)
        ggsave(paste0(output_prefix,'_',trend,  ".KEGG_cnetplot.png"),width = 10,height = 8,p1)
      }
    }
  }else {
    # 文件不存在，跳过
    message("File not found: ", file_path, " - skipping")
    next  # 如果在循环里用 next，否则用 return(NULL) 或什么都不做
  }
}

setwd('/mnt/NexData_file/home/zwl/Project/SC2510316157/require4/')
path = c('/mnt/NexData_file/home/zwl/Project/SC2510316157/require4/ALL/Dgal_vs_control/0.Dgal_vs_control_All_expressed_genes.xls',
         '/mnt/NexData_file/home/zwl/Project/SC2510316157/require4/ALL/EGT_vs_control/0.EGT_vs_control_All_expressed_genes.xls',
         '/mnt/NexData_file/home/zwl/Project/SC2510316157/require4/ALL/EGT_vs_Dgal/0.EGT_vs_Dgal_All_expressed_genes.xls')
for(n in path){
  
  data = read.table(file = n,sep = "\t",head=TRUE)
  
  colnames(data) <- c("Genes","scores","avg_log2FC","pvals","p_val_adj","pct_nz_group","pct_nz_reference","cluster")
  # data <- data %>%
  #   mutate(avg_log2FC = ifelse(cluster != "G4", -logfoldchanges, logfoldchanges))%>%
  #   dplyr::filter(pct_nz_reference>0.1)
  data$log10p_val_adj = -log10(data$p_val_adj)
  
  #data$log10p_val_adj[is.infinite(data$log10p_val_adj)] = 300
  cut_off_PValue =0.05 #设置FDR的阈值
  cut_off_log2FC =1 #设置log2FC的阈值
  data$Sig = ifelse(data$p_val_adj < cut_off_PValue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                      abs(data$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                    ifelse(data$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
  data = data.frame(data)
  table(data$Sig)
  # data = data %>%
  #   filter(pct_nz_group >= 0.2) %>%
  #   filter(pct_nz_reference <= 0.6) %>%
  #   arrange(desc(abs(log2FC)),FDR)
  xlim_min<-min(data$avg_log2FC)-0.3
  xlim_max<-max(data$avg_log2FC)+0.3
  ylim_max <- min(ceiling(max(data$log10p_val_adj) / 50) * 50,300)
  ###绘图——基础火山图###
  colmap = data.frame(type = c("Down","no","Up"),color=c("#546de5","#d2dae2","#F08080"))
  coloruse = colmap[colmap$type %in% data$Sig,]
  
  p1 <- ggplot(data, aes(x =avg_log2FC, y=log10p_val_adj, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
    geom_point(alpha=0.8, size=1) +  #点的透明度、大小
    scale_color_manual(values=coloruse$color) + 
    # xlim(-7,7) +  #调整点的颜色和x轴的取值范围
    geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
    geom_hline(yintercept = -log10(cut_off_PValue), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
    labs(x="log2FC", y="-log10p_val_adj") +  #x、y轴标签
    ggtitle("Volcano Plot") + #标题
    scale_x_continuous(limits = c(xlim_min, xlim_max)) +
    scale_y_continuous(limits = c(0, ylim_max), breaks = seq(0, ylim_max, 50)) +
    theme_bw() + # 主题，help(theme)查找其他个性化设置
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position="right", 
          legend.title = element_blank()
    ) 
  p1
  
  Up_top_10 =(     #筛选差异显著上调的前10个Gene
    data %>%
      filter(Sig == 'Up') %>%
      arrange(desc(abs(avg_log2FC)),p_val_adj) %>%
      head(5)
  )
  Up_top_10 %>% gt() #数据制成表
  Down_top_10 = (       #筛选差异显著下调的前10个Gene
    data %>%
      filter(Sig == 'Down') %>%
      arrange( desc(abs(avg_log2FC)),p_val_adj) %>%
      head(5)
  )
  Down_top_10 %>% gt() #数据制成表
  
  
  p2 <- p1 +
    geom_label_repel(data = Up_top_10,
                     aes(avg_log2FC, log10p_val_adj, label = Genes),
                     size = 3, fill="#FFCCCC",
                     alpha = 0.65, color = "black")+
    geom_label_repel(data = Down_top_10,
                     aes(avg_log2FC, log10p_val_adj, label = Genes),
                     size = 3, fill="#CCFFFF",
                     alpha = 0.65, color = "black",max.overlaps=10)
  
  p2 #出图
  output_dir <- dirname(normalizePath(n, winslash = "/", mustWork = TRUE))
  output_prefix <- file.path(output_dir, tools::file_path_sans_ext(basename(n)))
  ggsave(paste0(output_prefix, "_volcano.pdf"), plot = p2, width = 8, height = 6)
  ggsave(paste0(output_prefix, "_volcano.png"), plot = p2, width = 8, height = 6)
  
  
  
}

for (file_path in path) {
  gopathway(file_path,orgdb="org.Mm.eg.db",keggdb = "mmu",trend='ALL')
}
