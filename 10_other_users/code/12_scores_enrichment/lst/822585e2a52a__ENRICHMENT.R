library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(stringr)
library(forcats)
library(dplyr)
library(forcats)
library(org.Hs.eg.db)
library(clusterProfiler)

library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
setwd('F:/161_楼/')
diff<-read.csv('./table1.txt',sep = '\t')


diff = bitr(diff$gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)

go <- enrichGO(gene = diff$ENTREZID, 
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "ALL", 
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = T)

go.res <- data.frame(go)
goBP <- head(subset(go.res,subset = (ONTOLOGY == "BP")), 20)
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))
go.df <- rbind(goBP,goCC,goMF)
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
gr1.1 <- as.numeric(str_split(go.df$GeneRatio,"/",simplify = T)[,1])
gr2.1 <- as.numeric(str_split(go.df$GeneRatio,"/",simplify = T)[,2])
bg1.1 <- as.numeric(str_split(go.df$BgRatio,"/",simplify = T)[,1])
bg2.1 <- as.numeric(str_split(go.df$BgRatio,"/",simplify = T)[,2])
go.df$Gene.Ratio <- (gr1.1/gr2.1)
go.df <- arrange(go.df,Gene.Ratio)


ggplot(data = go.df,
       aes(x = Description, y = Count, fill = ONTOLOGY)) + # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip() +
  theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+ # 设置term名称过长时换行
  labs(x = "GO Terms",y = "Gene number")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11,lineheight = 0.8), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距

ggsave(filename = './GO_enrichment_bar.png', width=10, height = 10)
ggsave(filename = './GO_enrichment_bar.pdf', width=10, height = 10)
write.table(go.res,'./GO_enrichment.xls', sep='\t', row.names = F, quote =FALSE)


p<-ggplot(go.df, aes(x = Gene.Ratio,y =  reorder(Description,Gene.Ratio))) +
  geom_point(aes(color = p.adjust,size = Count)) +
  scale_color_gradient(low = "red", high = "blue") +
  xlab("Gene Ratio") +
  ylab("GO") +
  theme_bw() +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 50))+
  guides(color = guide_colorbar(reverse = TRUE))+
  theme(axis.text = element_text(lineheight = 0.7))
p<-p+facet_grid(ONTOLOGY~., scale = "free")
ggsave(p,file='./GO_enrichment_bubble.png',width = 9,height = 14)
ggsave(p,file='./GO_enrichment_bubble.pdf',width = 9,height = 14)

remotes::install_github("YuLab-SMU/createKEGGdb")
createKEGGdb::create_kegg_db('hsa')
install.packages("./KEGG.db_1.0.tar.gz",repos=NULL,type="source")




library('KEGG.db')

kegg <- enrichKEGG(gene = diff$ENTREZID,
                   keyType = "kegg",
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   use_internal_data = TRUE)

kegg_SYMBOL <- DOSE::setReadable(kegg,
                                 OrgDb = "org.Hs.eg.db",
                                 keyType = "ENTREZID")

data = as.data.frame(kegg_SYMBOL)
gr1 <- as.numeric(str_split(data$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(data$GeneRatio,"/",simplify = T)[,2])
bg1 <- as.numeric(str_split(data$BgRatio,"/",simplify = T)[,1])
bg2 <- as.numeric(str_split(data$BgRatio,"/",simplify = T)[,2])
data$Gene.Ratio <- (gr1/gr2)
data <- arrange(data,Gene.Ratio)
data$Description = factor(data$Description,levels = data$Description,ordered = T)

ggplot(data, aes(x = Gene.Ratio,y = reorder(Description,Gene.Ratio))) +
  geom_point(aes(color = p.adjust,size = Count)) +
  scale_color_gradient(low = "red", high = "blue") +
  xlab("Gene Ratio") +
  ylab("KEGG Pathway") +
  theme_bw() +
  guides(color = guide_colorbar(reverse = TRUE))
ggsave('./KEGG_enrichment.png')
ggsave('./KEGG_enrichment.pdf')
write.table(data.frame(Pathway=rownames(data),data),'./KEGG_enrichment.xls', sep='\t', row.names = F, quote =FALSE)










gr1.1 <- as.numeric(str_split(df$Ratio_in_study,"/",simplify = T)[,1])
gr2.1 <- as.numeric(str_split(df$Ratio_in_study,"/",simplify = T)[,2])
bg1.1 <- as.numeric(str_split(df$Ratio_in_pop,"/",simplify = T)[,1])
bg2.1 <- as.numeric(str_split(df$Ratio_in_pop,"/",simplify = T)[,2])
df$Richfactor <- (gr1.1/bg1.1)
df <- arrange(df,Richfactor)


p<-ggplot(df,aes(x=Richfactor,y=reorder(Description,-Padjust),size=Number,color=Padjust) )+
  geom_point()+
  theme_bw() +
  scale_colour_gradient(low="red",high="blue") + 
  labs(x="Rich Factor",y="KEGG pathway",title="",
       colour=expression(Pvalue),size="Number") +
  theme(plot.title = element_text(hjust = 0.5))
p<-p+facet_grid(Group~., scale = "free")
ggsave(p,file='./NS_vs_Nor_KEGG_bubble.pdf')
ggsave(p,file='./NS_vs_Nor_KEGG_bubble.png')

write.table(df,'./NS_vs_Nor_KEGG_bubble.xls', sep ='\t',row.names = F, quote =FALSE)
