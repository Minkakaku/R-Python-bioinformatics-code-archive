# conda activate scMetabolism
library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds")
levels(adata)

#### 分为>0 和 =0
marker_file1 <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = "\t",head=FALSE)
marker_file2 <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
marker1 <- marker_file1$V1
marker2 <- marker_file2$HALLMARK_GLYCOLYSIS
#转换成list
gene<-list(marker1,marker2)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'select_geneset'
colnames(adata@meta.data)[15]<-'GLYCOLYSIS'

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_g_select_geneset_barplot/"

# mydata<-FetchData(adata,vars=c('THBS1','select_geneset','HALLMARK_GLYCOLYSIS',"cell_type","group"))
# mydata<-mydata[mydata$THBS1!=0,]
# xlim <- range(mydata$THBS1)
p = VlnPlot(object = adata, features = 'THBS1',pt.size = 0)+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
ggsave(paste0(outdir,'mono_cell_type-THBS1_VlnPlot.pdf'),width = 7,height = 7)
ggsave(paste0(outdir,'mono_cell_type-THBS1_VlnPlot.png'),width = 7,height = 7)
p = VlnPlot(object = adata, features = 'THBS1',pt.size = 0,group.by='thbs1_group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
ggsave(paste0(outdir,'mono-THBS1_VlnPlot.pdf'),width = 7,height = 7)
ggsave(paste0(outdir,'mono-THBS1_VlnPlot.png'),width = 7,height = 7)

thbs1.expr <- FetchData(adata, vars = "THBS1")[[1]]
adata$thbs1_group <- ifelse(thbs1.expr > 0,"Positive", "Negative") %>%factor(levels = c("Negative","Positive"))
# adata$thbs1_group <- cut(thbs1.expr, breaks = c(-Inf, 0, 1, Inf), labels = c("THBS1_zero", "THBS1_low", "THBS1_high"), right = TRUE) %>%  factor() 

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_g_select_geneset_barplot/"
i='GLYCOLYSIS'
p = VlnPlot(object = adata, features = i,pt.size = 0, group.by = 'thbs1_group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
# ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
# ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot.png'),width = 7,height = 7)
p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=F)+ylim(-0.1,0.2)
# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot_sig_ana.pdf'),width = 10,height = 10)
ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot_sig_ana.png'),width = 10,height = 10)

