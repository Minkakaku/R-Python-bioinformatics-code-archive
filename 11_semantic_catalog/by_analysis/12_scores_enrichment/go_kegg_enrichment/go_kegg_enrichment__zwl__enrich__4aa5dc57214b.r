rm(list=ls())
gc()
gopathway = function(file_path,orgdb="org.Mm.eg.db",keggdb = "mmu",trend='ALL'){
  if (file.exists(file_path)) {
    output_dir <- dirname(normalizePath(file_path, winslash = "/", mustWork = TRUE))
    output_prefix <- file.path(output_dir, tools::file_path_sans_ext(basename(file_path)))
    
    data1 = read.table(file_path,sep = "\t",head=TRUE)
    colnames(data1)
    colnames(data1) <- c("gene_name", "scores", "log2FC","pvals", "PValue","pct_nz_group" , "pct_nz_reference")
    #自行选择筛选自己感兴趣的基因集
    if (trend =='up') {
      data <- data1[which(data1$PValue <= 0.05 & data1$log2FC >= 0.5 & data1$pct_nz_reference > 0.1),]
    }else if(trend =='down'){
      data <- data1[which(data1$PValue <= 0.05 & data1$log2FC <= -0.5 & data1$pct_nz_reference > 0.1),]
    }else{
      print('不分上下调')
      data = data1#[which(data1$PValue <= 0.05 & (abs(data1$log2FC) > 0.5) & data1$pct_nz_reference > 0.1),]
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
        
        # data_ALL = gene_symbol %>% left_join(data,by=c('SYMBOL'='gene_name'))
        # geneList = data_ALL$log2FC
        # names(geneList) = data_ALL$SYMBOL
        # edox <- setReadable(KEGG, 'org.Mm.eg.db', 'ENTREZID')
        # p1 <- cnetplot(KEGG,colorEdge = TRUE,foldChange=geneList,showCategory = 4)
        # ggsave(paste0(output_prefix,'_',trend,  ".KEGG_cnetplot.pdf"),width = 10,height = 8,p1)
        # ggsave(paste0(output_prefix,'_',trend,  ".KEGG_cnetplot.png"),width = 10,height = 8,p1)
      }
    }
  }else {
    # 文件不存在，跳过
    message("File not found: ", file_path, " - skipping")
    next  # 如果在循环里用 next，否则用 return(NULL) 或什么都不做
  }
}

setwd('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Myeloids_plot/enrich')
path = list.files(
  path = "/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Myeloids_plot/enrich",
  pattern = "^differentially_expressed_genes.xls$",  # 严格匹配 0xxxAll_expressed_genes.xls
  recursive = TRUE,  # 进入所有子文件夹
  full.names = TRUE   # 返回完整路径（方便直接读取）
)


for (file_path in path[3:10]) {
  gopathway(file_path,orgdb="org.Mm.eg.db",keggdb = "mmu",trend='ALL')
}


###################################monocle###########################################
# 沿时间轴的细胞密度图
setwd('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Macro_monocle')
HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Macro_monocle/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/Macro_monocle/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)

df = HSMM_myo@phenoData@data %>% dplyr::select(cell_type,Pseudotime,group) 
library(ggridges)
p_density <- ggplot(df, aes(x = Pseudotime, y = cell_type,fill =  cell_type)) +
  geom_density_ridges(rel_min_height = 0.005,alpha = 0.6)+
  theme_ridges()+
  scale_y_discrete(position = "right") + 
  theme(
    plot.background = element_rect('white'),
    legend.position = "none",
    panel.grid = element_blank() ,
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) 
p_density
ggsave(filename = "3.monocle.ridges_cortime.png",p_density,h=8,w=6,dpi = 150)
ggsave(filename = "3.monocle.ridges_cortime.pdf",p_density,h=8,w=6,dpi = 150)
diff_test_res <- differentialGeneTest(HSMM_myo[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(100) %>% row.names() 
write.csv(file = '3.obj_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE,
                            hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))

ggsave(filename = "3.monocle.heatmap_cortime.png",p,h=12,w=6,dpi = 150)
ggsave(filename = "3.monocle.heatmap_cortime.pdf",p,h=12,w=6,dpi = 150)
library(cowplot)
library(grid)

# 1️⃣ 把 pheatmap 转成 grob
heatmap_grob <- p$gtable

# 2️⃣ 拼图（上下结构）
final_plot <- plot_grid(
  p_density,
  heatmap_grob,
  ncol = 1,
  align = "v",
  rel_heights = c(0.3, 2)  # 上小下大，可以自己调
)

ggsave(filename = "3.monocle.heatmap_ridges_cortime.png",final_plot,h=14,w=6,dpi = 150)
ggsave(filename = "3.monocle.heatmap_ridges_cortime.pdf",final_plot,h=14,w=6,dpi = 150)

# 提取cluste对应的基因集
clusters <- cutree(p_pseudotime$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering, "clustering_B_cell.csv", row.names = F)

p_density
p_pseudotime

