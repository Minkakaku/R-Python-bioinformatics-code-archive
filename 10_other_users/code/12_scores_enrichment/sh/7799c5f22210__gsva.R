# conda activate "/mnt/NexData_file/share/scRNA-tools/envs/pip_GSVA/"
####分组GSVA热图####
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
# library(clusterProfiler)
library(patchwork)
library(limma)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20251753/02.Mission/script/KEGG.csv",sep = "\t",head=TRUE,strip.white = TRUE)
#转换成list
gene<-as.list(unique(geneset))
genesets<-lapply(gene, function(x) x[x != ""])

S <-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20251753/01.Mission/mapping_raw/Main_mapping.rds")
S <- subset(S,idents="Th17")
Idents(S) <- "group"
expr <- AverageExpression(S, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20251753/02.Mission/GSVA/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
pdf(paste0(outdir,pathway,"-heatmap.pdf"),width = 6,height = 8)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,pathway,"-heatmap.png"),width = 6,height = 8)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

library(stringr)
gsva.df <- gsva.df[,-1] 
rownames(gsva.df) <- gsub("_"," ",rownames(gsva.df))
rownames(gsva.df) <- str_to_sentence(str_to_lower(rownames(gsva.df))) 
# 设置细胞亚群名称
#colnames(gsva.df) <- c("Cluster 0", "Cluster 1", "Cluster 2","Cluster 3", "Cluster 4")

# 使用apply函数获取每列中最大的三个值的行索引
rows_to_keep <- unique(unlist(lapply(gsva.df, function(col) {
  order(col, decreasing = TRUE)[1:3] # 保留前3位
})))
# 根据索引保留这些行
filtered_gsva.df <- gsva.df[rows_to_keep, ]

library(pheatmap) 
# 使用pheatmap函数绘制热图
pheatmap(filtered_gsva.df,  # 输入的数据框，其中包含要展示的GSVA分析结果
         scale = "none",
         show_colnames = TRUE,  # 显示列名，即每个聚类的名称
         angle_col = "315",  # 列名的显示角度，设为315度使得列名从右上斜向左下显示
         cluster_rows = FALSE,  # 不进行行聚类，即不根据行数据对行进行聚类排序
         cluster_cols = FALSE,  # 不进行列聚类，即不根据列数据对列进行聚类排序
         border_color = "white",  # 设置单元格边框颜色为白色
         fontsize_col = 12,  # 设置列名的字体大小为12
         fontsize_row = 10)  # 设置行名的字体大小为10

pheatmap(gsva.df, show_colnames = T,cluster_cols=F,cluster_rows=F,
         scale = "none",angle_col = "315")



####多样本组间GSVA差异####
library(msigdbr) #提供MSigdb数据库基因集
library(GSVA)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(limma)

#KEGG
genesets <- msigdbr(species = "Mus musculus", category = "C2")
genesets <- genesets[genesets$gs_subcat%in%c("CP:KEGG"),]
genesets$gs_name <- gsub("KEGG_", "", genesets$gs_name)
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
#HALLMARK
genesets <- msigdbr(species = "Mus musculus", category = "H")
genesets$gs_name <- gsub("HALLMARK_", "", genesets$gs_name)
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)


S <- S[,S@meta.data$group %in% c("LPS","CON")]
Idents(S) <- S$sample
expr <- AverageExpression(S,assays = "RNA",slot = "data")[[1]]
expr <- as.matrix(expr)
gsva.res <- gsva(expr, genesets,method = "gsva")

#之前定义过分组信息
group_list <- c("LPS","LPS","LPS","CON","CON","CON")
#标准limma分析
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(expr)

contrast.matrix<-makeContrasts(contrasts = "LPS-CON",levels = design)
fit <- lmFit(gsva.res,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

#设置分组
diff$group <- ifelse( diff$logFC > 0 & diff$adj.P.Val < 0.05 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$adj.P.Val < 0.05 ,"down","noSig"))
diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>%
  arrange(t) %>% 
  rownames_to_column("ID")
diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

#diff2 <- diff2[abs(diff2$t)>4,]
diff2 <- top_n(diff2, n = 15, wt = abs(t))

#可视化
ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x=0,nudge_y=0,hjust =diff2$hjust,
            size=3)+
  labs(x="",
       y=paste0("t value of GSVA score\n","LPS vs CON"),
       title = "GSVA")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none" #去掉legend
  )
write.csv(diff,"Endo_gsva.csv")

####单样本组间GSVA差异####
library(msigdbr)
library(GSVA)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(limma)
genesets <- msigdbr(species = "Homo sapiens", category = "H")
genesets$gs_name <- gsub("HALLMARK_", "", genesets$gs_name)
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

expr <- as.matrix(A@assays$RNA@data)
gsva.res <- gsva(expr, genesets, method="gsva")

group_list <- A@meta.data[,c("group")]

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(expr)

contrast.matrix<-makeContrasts(contrasts = "P-C",levels = design)
fit <- lmFit(gsva.res,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

#设置分组
diff$group <- ifelse( diff$logFC > 0 & diff$adj.P.Val < 0.05 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$adj.P.Val < 0.05 ,"down","noSig"))
diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  filter(group != "noSig") %>%
  arrange(t) %>% 
  rownames_to_column("ID")
diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

#diff2 <- diff2[abs(diff2$t)>20,]
diff2 <- top_n(diff2, n = 15, wt = abs(t))

#可视化
ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x=0,nudge_y=0,hjust =diff2$hjust,
            size=3)+
  labs(x="",
       y=paste0("t value of GSVA score\n","P vs C"),
       title = "GSVA")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none" #去掉legend
  )
write.csv(diff,"KERA_COL7A1-_gsva.csv")





