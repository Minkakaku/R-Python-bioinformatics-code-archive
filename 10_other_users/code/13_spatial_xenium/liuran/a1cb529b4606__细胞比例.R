library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

dir.create("RCTD/Ratio")

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")


all_data_sp_rctd_anno <- readRDS("RCTD/xenium_rctd.rds")

out_gene_umi=as.matrix(all_data_sp_rctd_anno@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD/Ratio/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)



# 柱状图 ------
## 所有样本 ------
data <- read.table("RCTD/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,6,8,10)]
colnames(data) <- c("Cells","Group","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data <- data %>% filter(Annotation != "Unknown")


mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD/Ratio/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="RCTD/Ratio/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("RCTD/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="RCTD/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD/Ratio/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file="RCTD/Ratio/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")


# 限制样本顺序和分组顺序
order_sample <- read.table("RCTD/Ratio/order_sample.txt",header=F,sep="\t",quote="",check.names=F)
order_group <- read.table("RCTD/Ratio/order_group.txt",header=F,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = order_sample[,1])
data4plot2$Group <- factor(data4plot2$Group, levels = order_group[,1])


plot_list <- c()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  # scale_fill_manual(values = plan)+
  ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  # scale_fill_manual(values = plan)+
  ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2.5, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD/Ratio/ratio_barplot.pdf",width=12, height=6)
ggsave(p,file="RCTD/Ratio/ratio_barplot.png",width=12, height=6,dpi=300)



## Air ------
data <- read.table("RCTD/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,6,8,10)]
colnames(data) <- c("Cells","Group","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data <- data %>% filter(Annotation != "Unknown")
data <- data %>% filter(Group %in% c("AirMale", "AirFemale"))
table(data$Group)

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD/Ratio/Cell_type_stat_Air.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="RCTD/Ratio/Cell_type_stat_percent_Air.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("RCTD/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD/Ratio/Cell_type_stat_group_Air.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file="RCTD/Ratio/Cell_type_stat_percent_group_Air.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
order_sample <- read.table("RCTD/Ratio/order_sample.txt",header=F,sep="\t",quote="",check.names=F)
order_group <- read.table("RCTD/Ratio/order_group.txt",header=F,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = order_sample[,1])
data4plot2$Group <- factor(data4plot2$Group, levels = order_group[,1])

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2.5, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD/Ratio/ratio_barplot_Air.pdf",width=12, height=6)
ggsave(p,file="RCTD/Ratio/ratio_barplot_Air.png",width=12, height=6,dpi=300)



## Hyperoxia ------
data <- read.table("RCTD/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,6,8,10)]
colnames(data) <- c("Cells","Group","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data <- data %>% filter(Annotation != "Unknown")
data <- data %>% filter(Group %in% c("HyperoxiaMale", "HyperoxiaFemale"))
table(data$Group)

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD/Ratio/Cell_type_stat_Hyperoxia.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="RCTD/Ratio/Cell_type_stat_percent_Hyperoxia.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("RCTD/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD/Ratio/Cell_type_stat_group_Hyperoxia.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file="RCTD/Ratio/Cell_type_stat_percent_group_Hyperoxia.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
order_sample <- read.table("RCTD/Ratio/order_sample.txt",header=F,sep="\t",quote="",check.names=F)
order_group <- read.table("RCTD/Ratio/order_group.txt",header=F,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = order_sample[,1])
data4plot2$Group <- factor(data4plot2$Group, levels = order_group[,1])

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2.5, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD/Ratio/ratio_barplot_Hyperoxia.pdf",width=12, height=6)
ggsave(p,file="RCTD/Ratio/ratio_barplot_Hyperoxia.png",width=12, height=6,dpi=300)



# # 双向图 ------
# info<-read.table(file="RCTD/Ratio/Cell_type_stat_percent_group.txt",header=T,sep="\t",comment.char = "",check.names=F)
# info <- melt(info)
# colnames(info) <- c("celltype","group","percent")
# 
# sample_map<-read.table(file="RCTD/Ratio/sample_map_group.txt",header=T,sep="\t",comment.char = "")
# orders <- read.table("RCTD/Ratio/order.txt",header=F,sep="\t",quote="",check.names=F)
# info$celltype <- factor(info$celltype,levels=unique(orders[,1]))
# info$group <- factor(info$group,levels=unique(sample_map[,2]))
# 
# info1 <- info[which(info$group %in% c("AirMale","AirFemale","HyperoxiaMale","HyperoxiaFemale")),]
# 
# data1 <- lapply(unique(orders[,1]),function(x){
#   n <- c()
#   tmp <- info1[which(info1$celltype==x),]
#   value <- tmp[which(tmp$group=="24M"),]$percent/tmp[which(tmp$group=="4M"),]$percent
#   n <- append(n,x)
#   n <- append(n,value)
#   n
# }) %>% do.call("rbind", .) 
# 
# data1 <- as.data.frame(data1)
# colnames(data1) <- c("celltype","value")
# data1$value <- as.numeric(data1$value)
# data1$log2 <- log2(data1$value)
# write.table(data1,file="RCTD/Ratio/log2_barplot.txt",quote=F,sep="\t",row.names = F)
# 
# 
# color <- rep("#ae4531", nrow(data1))
# color[which(data1$log2 < 0)] <- "#2f73bb"
# data1$color<-color
# data1$log2 <- as.numeric(data1$log2)
# up <- data1[which(data1$log2 > 0),]
# down<- data1[which(data1$log2 < 0),] #主要为了后面方便添加X轴的标签
# 
# 
# p <- ggplot(data=data1,aes(x=reorder(celltype,log2),log2,fill=color))+
#   geom_bar(stat = "identity")+
#   theme_classic()+
#   ylim(-1.5,1.5)+ # 限制X轴范围
#   coord_flip()+
#   scale_fill_manual(values = c("#2f73bb","#ae4531"))+
#   #geom_hline(yintercept=0,linetype=1)+
#   geom_segment(aes(y = 0, yend = 0,x = 0, xend = nrow(data1)+1))+ #geom_segment比genom_line更为灵活，因为可以定义起始位置
#   geom_segment(aes(y = 1, yend = 1,x = 0, xend = nrow(data1)+1),linetype=2,color="grey")+
#   geom_segment(aes(y = -1, yend = -1,x = 0, xend = nrow(data1)+1),linetype=2,color="grey")+
#   theme(
#     legend.position = "none", #消除图例
#     plot.title = element_text(hjust = 0.5), #标题居中
#     axis.line.y = element_blank(), #删除Y轴
#     axis.title.y = element_blank(), #删除Y轴标题
#     axis.ticks.y = element_blank(), #删除Y轴刻度
#     axis.text.y = element_blank(), #删除Y轴标签
#     axis.text.x=element_text(size=15,color="black"),
#     legend.text=element_text(size=15,color="black"),
#     axis.title.x=element_text(size=15,face = "bold",color="black")
#   )+
#   labs(y="log2(Model/Control)")+
#   geom_text(data=up,aes(x = celltype, y = -0.1, label = celltype), hjust=1,size = 5)+ #添加score是正的label
#   geom_text(data=down,aes(x = celltype, y = 0.1, label = celltype), hjust=0,size = 5) #添加score是负的label
# ggsave(p,file="RCTD/Ratio/log2_barplot.pdf",width=8,height=5)
# ggsave(p,file="RCTD/Ratio/log2_barplot.png",width=8,height=5)






# 箱线图（带显著性标记）------
library(ggpubr)

## 普通版 所有样本 ------
data<-read.table(file = "RCTD/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"
data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by="SampleID")

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="RCTD/Ratio/diff_boxplot.pdf",width=15,height=5)
ggsave(p,file="RCTD/Ratio/diff_boxplot.png",width=15,height=5)




## 普通版 Air ------
data<-read.table(file = "RCTD/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"
data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by="SampleID")

# filter
plot_data <- plot_data %>% filter(Group %in% c("AirMale", "AirFemale"))

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="RCTD/Ratio/diff_boxplot_Air.pdf",width=15,height=5)
ggsave(p,file="RCTD/Ratio/diff_boxplot_Air.png",width=15,height=5)




## 普通版 Hyperoxia ------
data<-read.table(file = "RCTD/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"
data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by="SampleID")

# filter
plot_data <- plot_data %>% filter(Group %in% c("HyperoxiaMale", "HyperoxiaFemale"))

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="RCTD/Ratio/diff_boxplot_Hyperoxia.pdf",width=15,height=5)
ggsave(p,file="RCTD/Ratio/diff_boxplot_Hyperoxia.png",width=15,height=5)


