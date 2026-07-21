library(data.table)
library(readxl)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ggsci)
library(ggh4x)
library(ropls)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(reshape2)

############################## Parameters ######################################
# workflow = "Metabolomics" # Must be one of "Metabolomics", "Lipidomics" 
# organism = "hsa" # Must be organism IDs from KEGG 
# sample_info = "input/FWMS20252134/samples.csv"
# comparison = "input/FWMS20252134/comparison.txt"
# results_dir = "Results/"
# p.method = "none" # 是否进行pvalue校正 "BH","BY","holm"
# p.value = 0.05 # 差异阈值
# logfc = 0 # 差异阈值
# VIP = 1 # 差异阈值
# diff.method = "T test" # One of "T test" or "Paired T test" 

p.method = snakemake@params[["p_method"]]
p.value = snakemake@params[["p_value"]]
logfc = snakemake@params[["logfc"]]
VIP = snakemake@params[["VIP"]] 
diff.method = snakemake@params[["diff_method"]]

workflow = snakemake@params[["workflow"]]
organism = snakemake@params[["organism"]]
data.final.table.file = snakemake@input[["data_final_table"]]
sample_info = snakemake@input[["sample_info"]]
comparison = snakemake@input[["comparison"]]
results_dir = snakemake@params[["results_dir"]]

############################ Functions defined #################################
strip <- strip_themed(background_y = elem_list_rect(fill = brewer.pal(9,"Set3")),text_x = elem_list_text(face=2,size=12),background_x = elem_list_rect(fill = pal_npg()(10)))
strip1 <- strip_themed(background_x = elem_list_rect(fill = pal_npg()(10)[3:10]),text_x = elem_list_text(face=2,size=12))

significance_B=function(log2FC,intensity){
  df=length(log2FC)/1000
  if(df<5) df = 5
  intensity=log2(intensity)
  a=(intensity-min(intensity)+1)*0.5
  t=quantile(log2FC,probs = c(0.5,0.8413))
  sapply(1:length(log2FC),function(x){
    z=abs((log2FC[x]-t[[1]])/(t[[2]]-t[[1]]))*a[x]
    2*pt(-z,df = df)
  },simplify = T)
}

my.diff=function(data,x,orthoI=NA,paired=F,diff.method,p.adjust = "BH", v=1, f=0.01, fc=0,type="Metabolomics"){
  data=data %>% filter(Group %in% c(x$Case,x$Control))
  sample=data %>% distinct(Sample,Group)
  data.matrix=as.matrix(dcast(data.table(data),`Alignment ID`~Sample,value.var = "value") %>% column_to_rownames("Alignment ID"))
  
  if (diff.method=="significance_B"){
    data.diff=data %>% group_by(`Alignment ID`) %>% arrange(Group,Sample) %>% mutate(log2fc=mean(value[Group==x$Case])-mean(value[Group==x$Control])) %>% dplyr::select(-value) %>% distinct() %>% ungroup() %>%
      mutate(p.value=significance_B(log2fc,as.numeric(Intensity))) %>% ungroup()
    sacurine.plsda = NULL
    data.vip=data.frame(V1=unique(data.diff$`Alignment ID`),VIP = 1.1)
    colnames(data.vip)[1]="Alignment ID"
  } else {
    #------------------- OPLS-DA -------------------#
    ci=7;if(nrow(sample)<7) ci=nrow(sample)
    sacurine.plsda <- opls(t(data.matrix[,as.character(sample$Sample)]),sample$Group,orthoI = orthoI, predI =1, crossvalI=ci) 
    if(nrow(sacurine.plsda@modelDF)==0){sacurine.plsda <- opls(t(data.matrix[,as.character(sample$Sample)]),sample$Group,orthoI = 1, predI =1, crossvalI=ci) }
    data.vip=data.frame(VIP=getVipVn(sacurine.plsda)) %>% rownames_to_column("Alignment ID")
 
    #------------------- T-test -------------------#
    if(paired==T){
      data$Ind = unlist(lapply(str_split(data$Sample, pattern = "_"),function(x) x[[1]]))
      data.diff=data %>% dcast(`Alignment ID`+`Metabolite name`+`Category1`+`Category2`+`HMDB ID`+Ind~Group, value.var="value") %>% group_by(`Alignment ID`) %>%
        mutate(p.value=t.test(x=cur_data()[[x$Case]],y=cur_data()[[x$Control]],paired=paired)$p.value,
              log2fc=mean(cur_data()[[x$Case]])-mean(cur_data()[[x$Control]])) %>% ungroup()
    }else{
      data.diff=data %>% group_by(`Alignment ID`) %>% arrange(Group,Sample) %>% mutate(p.value=t.test(value~Group)$p.value) %>%
        mutate(log2fc=mean(value[Group==x$Case])-mean(value[Group==x$Control])) %>% ungroup()
    }

  }
  if(type=="Metabolomics"){
    data.diff=data.diff %>% distinct(`Alignment ID`,`Metabolite name`, `HMDB_SuperClass`,`HMDB_Class`,`HMDB_SubClass`, `HMDB ID`,Intensity, log2fc, p.value) %>% 
      mutate(fdr=p.adjust(p.value, method= p.adjust)) %>% 
      merge(data.vip, by="Alignment ID") %>% 
      mutate(regulation=case_when(fdr<f&VIP>v&log2fc>fc~'up',
                                  fdr<f&VIP>v&log2fc< -fc~'down',TRUE~'unchanged')) 
    
    data.diff.long=subset(data,select=c("Sample","Group","Metabolite name","value")) %>% 
      merge(data.diff %>% filter(regulation!='unchanged') %>% 
      filter(nchar(`Metabolite name`)< 50) %>% group_by(regulation) %>%  
      slice_min(p.value,n=15),by="Metabolite name") %>% 
      filter(Group %in% c(x$Case,x$Control)) %>% 
      group_by(`Metabolite name`) %>% 
      mutate(value=value-mean(value))
    
  }
  if (type=="Lipidomics"){
    data.diff=data.diff %>% distinct(`Alignment ID`,`Metabolite name`,Abbreviation,Categories,`LION ID`,`LION name`,Name,Intensity,log2fc,p.value) %>% 
      mutate(fdr=p.adjust(p.value,method= p.adjust)) %>%
      merge(data.vip,by="Alignment ID") %>% 
      mutate(regulation=case_when(fdr<f&VIP>v&log2fc>fc~'up',fdr<f&VIP>v&log2fc< -fc~'down',TRUE~'unchanged'))
    
    data.diff.long=subset(data,select=c("Sample","Group","Metabolite name","value")) %>% 
      merge(data.diff %>% filter(regulation!='unchanged') %>% 
      filter(nchar(`Metabolite name`)< 50) %>% group_by(regulation) %>% 
      slice_max(abs(log2fc),n= 50) %>% 
      slice_min(p.value,n=25),by="Metabolite name") %>% 
      filter(Group %in% c(x$Case,x$Control)) %>% 
      group_by(`Metabolite name`) %>% 
      mutate(value=value-mean(value))
  }

  return(list(model=sacurine.plsda,matrix = data.matrix, data=data.diff, data.diff.long=data.diff.long ))
}
plada.plot=function(sacurine.plsda,data.matrix,s.color="red"){
  print(sacurine.plsda@modelDF)
  p1=ggplot(as.data.frame(sacurine.plsda@modelDF[,c(3,5)]) %>% rownames_to_column("Component") %>% filter(Component!='sum') %>% melt(id.vars=1) %>% mutate(Component=factor(Component,levels=unique(Component))) %>%mutate(variable=factor(variable,levels=c("R2Y","Q2"))),
            aes(Component,value*100,fill=variable))+geom_col(position = position_dodge2(width = 1))+theme_classic(base_size = 14)+scale_fill_manual(values=c(R2Y=pal_npg()(2)[1],Q2=pal_npg()(2)[2]))+labs(title = "Model Overview",y="% of total variance")+theme(plot.title = element_text(hjust = 0.5))
  p2=ggplot(data.frame(sacurine.plsda@suppLs$permMN[,c(2,3,7)]) %>% melt(id.vars="sim"),aes(sim,value,color=variable,shape=variable))+geom_point(size=3)+geom_smooth(method='lm',se=F)+theme_classic(base_size = 14)+
    labs(title = "Permutation testing",x="Similarity(Y,Yperm)",subtitle = paste0("pR2Y = ", sacurine.plsda@summaryDF$pR2Y," , pQ2 = ",sacurine.plsda@summaryDF$pQ2))+scale_color_manual(values=c(R2Y.cum.=pal_npg()(2)[1],Q2.cum.=pal_npg()(2)[2]))+theme(plot.title = element_text(hjust = 0.5))
  p3=ggplot(data.frame(P1.cov=cov(t(data.matrix),sacurine.plsda@scoreMN[,1],),VIP=sacurine.plsda@vipVn,P1.cor=cor(t(data.matrix),sacurine.plsda@scoreMN)[,1]),
            aes(P1.cov,P1.cor,color=VIP))+geom_point(size=2)+scale_color_material(s.color)+theme_classic(base_size = 14)+labs(title = 'Loading S-plot (OPLS-DA)',x="P1 Covariance",y="P1 Correlation")+theme(plot.title = element_text(hjust = 0.5))
  p4=ggplot(data.frame(P1=sacurine.plsda@scoreMN[,1],O1=sacurine.plsda@orthoScoreMN[,1],Group=sacurine.plsda@suppLs$y),
            aes(P1,O1,color=Group,shape=Group))+geom_point(size=3)+scale_color_manual(values=setNames(pal_npg()(2)[1:2],nm = unique(sacurine.plsda@suppLs$y)))+theme_classic(base_size = 14)+stat_ellipse()+
    labs(title = "Score (OPLS-DA) plot",x=sprintf("Component 1 (%.2f%%)",sacurine.plsda@modelDF["p1","R2X"]*100),y=sprintf("Orthogonal Component 1 (%.2f%%)",sacurine.plsda@modelDF["o1","R2X"]*100))+theme(plot.title = element_text(hjust = 0.5))
  ggarrange(p1,p2,p3,p4)
  return(ggarrange(p1,p2,p3,p4))
}
volcano.plot=function(data.diff,fdr=0.01,fc=0,xmax=as.integer(max(abs(data.diff$log2fc))),ymax=as.integer(max(abs(-log10(data.diff$fdr)))),p.adjust="none"){
  p=ggplot(data.diff, aes(log2fc, -log10(fdr)))+
    geom_point(aes(color=regulation, size = VIP))+
    scale_size(range=c(0.1, 4))+ 
    scale_color_manual(name = "Regulation",values=c(up=pal_npg()(2)[1],down= pal_npg()(2)[2],unchanged='grey'),labels=(data.frame(table(data.diff$regulation)) %>% mutate(label=paste(Var1,Freq,sep=": ")))$label)+
    theme_classic(base_size = 14)+xlim(-xmax, xmax) + ylim(0, ymax) 
  if (p.adjust =="none"){p1=p+labs(y="-log10(p.value)")}else{p1=p+labs(y="-log10(adj.pvalue)")}
  if (fc==0){p1=p1+geom_hline(yintercept = -log10(fdr), linetype=2)}else{p1=p1+geom_hline(yintercept = -log10(fdr), linetype=2)+geom_vline(xintercept=c(-fc, fc),linetype=2)}
  return(p1)
}

######################## Differential analysis  ################################
group.data = fread(comparison)
if(workflow=="Metabolomics"){
  data.final.table = fread(data.final.table.file)
}
if(workflow=="Lipidomics"){
  data.final.table = read_excel(data.final.table.file)
}

sample = fread(sample_info)
if(all(sample$Sample == "-") == T){
  sample = sample %>% mutate(Sample = ID)
} else {
  sample = sample 
}

data.long = data.final.table %>% melt(measure.vars = sample$Sample, variable.name = "Sample") %>% merge(sample,by="Sample")

if(diff.method=="Paired-T test"){paired=T} else{paired=F}
diff.data=lapply(split(group.data, group.data$Label), function(x){my.diff(data.long, x, orthoI = NA,diff.method = diff.method, paired = paired, p.adjust = p.method,
                                                                                               v = VIP, f = p.value, fc = logfc, type = workflow)})
data.diff.combine=data.frame()
for (x in names(diff.data)){
  temp=diff.data[[x]]$data[,c("Alignment ID","log2fc","p.value","fdr","VIP","regulation")] %>% column_to_rownames(var ="Alignment ID")
  colnames(temp) = paste(colnames(temp),x, sep="_")
  if(nrow(data.diff.combine)==0){
    data.diff.combine=temp
  }else{
    data.diff.combine=cbind(data.diff.combine, temp)
  }
}
data.diff.combine=data.diff.combine %>% rownames_to_column(var = "Alignment ID")

values = list()
for (x in group.data$Label){
  i1=paste(x,"table",sep=".")
  i2=paste(x,"plot1",sep=".")
  i3=paste(x,"plot2",sep=".")
  i4=paste(x,"plot3",sep=".")
  i5=paste(x,"table.list",sep=".")
  if(nrow(diff.data[[x]]$data)==0){
    values[[i2]]=NULL
    values[[i3]]=NULL
    values[[i4]]=NULL
    values[[i5]]=NULL
  }else{
    if(diff.method=="significance_B"){
      values[[i2]]= ggplot()+theme(panel.background = element_blank())
      values[[i5]]=NULL
    }else{
      values[[i2]]=plada.plot(diff.data[[x]]$model,diff.data[[x]]$matrix)
      values[[i5]]=list(
        oplsda.table1=as.data.frame(diff.data[[x]]$model@modelDF[,c(3,5)]) %>% rownames_to_column("Component") %>% filter(Component!='sum') %>% melt(id.vars=1) %>% mutate(Component=factor(Component,levels=unique(Component))) %>%mutate(variable=factor(variable,levels=c("R2Y","Q2"))),
        oplsda.table2=data.frame(diff.data[[x]]$model@suppLs$permMN[,c(2,3,7)]) %>% melt(id.vars="sim"),
        oplsda.table3=data.frame(P1.cov=cov(t(diff.data[[x]]$matrix),diff.data[[x]]$model@scoreMN[,1],),VIP=diff.data[[x]]$model@vipVn,P1.cor=cor(t(diff.data[[x]]$matrix),diff.data[[x]]$model@scoreMN)[,1]),
        oplsda.table4=data.frame(P1=diff.data[[x]]$model@scoreMN[,1],O1=diff.data[[x]]$model@orthoScoreMN[,1],Group=diff.data[[x]]$model@suppLs$y) %>% rownames_to_column(var = "Sample")
      )
    }
    values[[i3]]=volcano.plot(data.diff=diff.data[[x]]$data,fdr=p.value,fc=logfc,p.adjust=p.method)
    
    data.diff.long=diff.data[[x]]$data.diff.long
   
    if(nrow(data.diff.long)>0){
      c.angle <- 60; s.value <- 10
      if (workflow =="Metabolomics"){
        p6=ggplot(data.diff.long)+geom_tile(aes(Sample,reorder(`Metabolite name`,log2fc),fill=value),linewidth=0.1,color='white')+
          facet_grid2(HMDB_SuperClass~Group,scales = 'free',space = 'free',switch = "y",strip = strip)+
          scale_fill_gradientn(colors = c(pal_material("blue",reverse = T)(s.value),'white',pal_material("red")(s.value)))+theme_test(base_size=14)+
          theme(axis.text.x = element_text(angle = c.angle,hjust=1),legend.position ="top" ,legend.justification = 0,axis.title.y=element_blank(),
                strip.text.y.left=element_text(angle = 0))
        p7=ggplot(data.diff.long,aes(log2fc,reorder(`Metabolite name`,log2fc)))+
          geom_segment(aes(x=0,xend=log2fc,y=reorder(`Metabolite name`,log2fc),yend=reorder(`Metabolite name`,log2fc)),linewidth=0.4,linetype=3)+
          geom_point(aes(size=VIP,color=-log10(fdr)))+scale_color_material("pink",limits=c(1,5),oob=scales::squish)+theme_test(base_size = 14)+
          facet_grid2(HMDB_SuperClass~'Log2FC',scales = 'free',space = 'free',strip=strip1)+scale_size(rang=c(1,4))+geom_vline(xintercept = 0,color='grey')+
          theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y=element_blank(),
                legend.position = "top",legend.justification = 1,strip.text.y=element_blank())+guides(size=guide_legend(ncol = 2))
      }
      if (workflow =="Lipidomics"){
        p6=ggplot(data.diff.long)+geom_tile(aes(Sample,reorder(`Metabolite name`,log2fc),fill=value),linewidth=0.1,color='white')+
          facet_grid2(Categories~Group,scales = 'free',space = 'free',switch = "y",strip = strip)+
          scale_fill_gradientn(colors = c(pal_material("blue",reverse = T)(s.value),'white',pal_material("red")(s.value)),limits=c(-2,2),oob=scales::squish)+theme_test(base_size=14)+
          theme(axis.text.x = element_text(angle = c.angle,hjust=1),legend.position = "top",legend.justification = 0,axis.title.y=element_blank(),
                strip.text.y.left=element_text(angle = 0))
        p7=ggplot(data.diff.long,aes(log2fc,reorder(`Metabolite name`,log2fc)))+
          geom_segment(aes(x=0,xend=log2fc,y=reorder(`Metabolite name`,log2fc),yend=reorder(`Metabolite name`,log2fc)),linewidth=0.4,linetype=3)+
          geom_point(aes(size=VIP,color=-log10(fdr)))+scale_color_material("pink",limits=c(1,5),oob=scales::squish)+theme_test(base_size=14)+
          facet_grid2(Categories~'Log2FC',scales = 'free',space = 'free',strip=strip1)+scale_size(rang=c(1,4))+geom_vline(xintercept = 0,color='grey')+
          theme(axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y=element_blank(),legend.position = "top",legend.justification = 1,
                strip.text.y=element_blank())+guides(size=guide_legend(ncol = 2))
      }
      values[[i4]]=ggarrange(p6,p7,widths = c(4,1),align = 'h')}
  }
}

########################## Export Results #####################################
if(workflow=="Metabolomics"){temp.dir=paste(results_dir,"3-Differentially abundant metabolites",sep="/")}else{temp.dir=paste(results_dir,"3-Differentially abundant lipids",sep="/")}
if(dir.exists(temp.dir)) unlink(temp.dir,recursive = T,force = T)
dir.create(temp.dir)

diff.table=data.frame()
for(x in names(diff.data)){
  group = unlist(stringr::str_split(x,"-vs-"))
  data.diff.long.temp=data.long %>% filter(Group %in% group) %>% distinct(`Alignment ID`, Sample, value) %>% merge(diff.data[[x]]$data,by="Alignment ID")
  if (workflow =="Metabolomics"){
    values[[paste(x,"diff.table",sep=".")]]=reshape2::dcast(as.data.frame(data.diff.long.temp), `Alignment ID`+`Metabolite name` + HMDB_SuperClass + HMDB_Class + HMDB_SubClass + `HMDB ID` + log2fc + p.value + fdr + VIP + regulation ~ Sample) %>% arrange(regulation)
    fwrite(values[[paste(x,"diff.table",sep=".")]],file=paste0(temp.dir,"/Table 3a. Differentially abundant metabolites list of ",x,".xls"),sep = "\t")
  }
  if (workflow =="Lipidomics"){
    values[[paste(x,"diff.table",sep=".")]]=reshape2::dcast(as.data.frame(data.diff.long.temp), `Alignment ID`+`Metabolite name` + Abbreviation + Categories+`LION ID`+`LION name` + Name + log2fc + p.value + fdr + VIP + regulation ~ Sample) %>% arrange(regulation)
    fwrite(values[[paste(x,"diff.table",sep=".")]],file=paste0(temp.dir,"/Table 3a. Differentially abundant lipids list of ",x,".xls"),sep = "\t")
  }
  
  t=data.frame(Label=x,
               `up-regulated`=sum(values[[paste(x,"diff.table",sep=".")]]$regulation=='up'),
               `down-regulated`=sum(values[[paste(x,"diff.table",sep=".")]]$regulation=='down')
  )
  diff.table=rbind(diff.table,t)
  # saveRDS(diff.data[[x]]$model, file = paste0(temp.dir, "/OPLSDA ",x, ".RData"))
  s=sample %>% filter(Group %in% group)
  if(max(nchar(s$Sample))>20){width1=14}else{width1=10}
  if(max(nchar(data.diff.long.temp$`Metabolite name`))>30){width2=16}else{width2=14}
  if(nrow(s)>30 & nrow(s)<50){width2=width2+nrow(s)*0.15}
  if(nrow(s)>50){width2=width2+nrow(s)*0.08}
  
  if(max(nchar(s$Group))>20){width=12}else{width=8}
  # ggsave(values[[paste(x,"plot1",sep=".")]],filename = paste0(temp.dir,"/Figure 3a. OPLS-DA plot of ",x,".pdf"),width = width1,height=8)
  # ggsave(values[[paste(x,"plot1",sep=".")]],filename = paste0(temp.dir,"/Figure 3a. OPLS-DA plot of ",x,".png"),width = width1,height=8)
  ggsave(values[[paste(x,"plot2",sep=".")]],filename = paste0(temp.dir,"/Figure 3b. Volcano plot of ",x,".pdf"),width = 6,height=6)
  ggsave(values[[paste(x,"plot2",sep=".")]],filename = paste0(temp.dir,"/Figure 3b. Volcano plot of ",x,".png"),width = 6,height=6)
  ggsave(values[[paste(x,"plot3",sep=".")]],filename = paste0(temp.dir,"/Figure 3c. Heatmap of ",x,".pdf"),width = width2,height=10)
  ggsave(values[[paste(x,"plot3",sep=".")]],filename = paste0(temp.dir,"/Figure 3c. Heatmap of ",x,".png"),width = width2,height=10)
  diff.table = diff.table
  if(is.null(values[[paste(x,"table.list",sep=".")]]) == F){
    ggsave(values[[paste(x,"plot1",sep=".")]],filename = paste0(temp.dir,"/Figure 3a. OPLS-DA plot of ",x,".pdf"),width = width1,height=8)
    ggsave(values[[paste(x,"plot1",sep=".")]],filename = paste0(temp.dir,"/Figure 3a. OPLS-DA plot of ",x,".png"),width = width1,height=8)
    
    writexl::write_xlsx(values[[paste(x,"table.list",sep=".")]],path = paste0(temp.dir,"/Table 3c. OPLS-DA results of ",x,".xlsx"))
    saveRDS(diff.data[[x]]$model, file = paste0(temp.dir, "/OPLSDA ",x, ".RData"))
  }
}
fwrite(as.data.frame(merge(data.final.table,data.diff.combine,by="Alignment ID")), file=paste0(temp.dir, "/Table 3b. Combined DEMs list.xls"), sep ="\t")
write.table(diff.table, file=paste0(temp.dir, "/diff.stats.txt"), sep = "\t", row.names = F, quote = F)
