# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(seqinr)
library(stringr)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

##################### parameter #########################################
search.soft <- snakemake@params[["search_soft"]]
sample.file <- snakemake@input[["sample_file"]]
protein.file <- snakemake@input[["protein_file"]]
evidence.file <- snakemake@input[["evidence_file"]]
summary.file <- snakemake@input[["summary_file"]]
db.dir <- snakemake@input[["db_dir"]]
output.dir <- snakemake@params[["output_dir"]]


##################### function defined #########################################
seqper=function(out,ID,seq){
  o=stringr::str_locate_all(out[[ID[1]]][1],pattern = unique(seq))
  round(length(unique(unlist(lapply(o,function(x){
    if(length(x)>0) x[1]:x[2] else NULL
  }))))/nchar(out[[ID[1]]])*100,2)
}
evidence.plot1=function(data.evidence.sub){
  ggplot(data.evidence.sub %>% filter(`Mass error [ppm]`<10,`Mass error [ppm]`> -10),aes(x=`Mass error [ppm]`))+geom_histogram(bins = 150,fill=pal_npg()(3)[3],color="white",linewidth=0.1)+
    theme_classic(base_size = 14)+xlim(-10,10)+labs(title = "Peptide mass tolerance distribution",x="Mass error [ppm]")
}
evidence.plot2=function(data.evidence.sub){
  missed.cleavages=data.frame(table(data.evidence.sub$`Missed cleavages`))
  colnames(missed.cleavages)=c("Missed\ncleavages","Count")
  missed.cleavages$Perentage=sprintf("%.1f%%",missed.cleavages$Count/sum(missed.cleavages$Count)*100)
  data.evidence.sub$`Missed cleavages`=as.factor(data.evidence.sub$`Missed cleavages`)
  ggplot(subset(data.evidence.sub,Length<30),aes(x=as.factor(Length),fill=`Missed cleavages`))+
    geom_bar(color="white",stat = "count")+theme_classic(base_size = 14)+
    annotation_custom(grob = tableGrob(missed.cleavages,rows = NULL,theme = ttheme_default(base_size = 10)),xmin = 15,ymin = nrow(data.evidence.sub)/30)+
    theme(legend.position= c(0.5,0.95),legend.direction = "horizontal")+
    scale_fill_jama()+labs(title = "Peptide length distribution",x="Peptide length")
}
evidence.plot3=function(data.evidence.sub){
  quant.header=grep(colnames(data.evidence.sub),pattern = "Reporter intensity \\d",value = T)
  if(length(quant.header)>0){
    reporter.intensity=subset(data.evidence.sub,select = c("id",quant.header))
    colnames(reporter.intensity)=c("id",sub(quant.header,pattern = "Reporter intensity ",replacement = ""))
    reporter.intensity=melt(reporter.intensity,id.vars="id") %>% filter(value>0)
    ggplot(reporter.intensity,aes(variable,log2(value)))+
      geom_violin(aes(fill=variable))+geom_boxplot(width=0.3)+theme_test(base_size = 14)+
      theme(legend.position = "none")+labs(title = "Reporter intensity distribution",y="Log2 intensity",x="TMT label channel")+
      scale_fill_manual(values = colorRampPalette(pal_npg()(10))(length(unique(reporter.intensity$variable))))
  }else{
    ggplot(data.evidence.sub,aes(Experiment,log2(Intensity)))+
      geom_violin(aes(fill=Experiment))+geom_boxplot(width=0.3)+theme_test(base_size = 14)+
      theme(axis.text.x=element_text(angle=60,hjust=1),legend.position = "none")+
      labs(title = "Intensity distribution",y="Log2 intensity",x="Experiment sample")+
      scale_fill_manual(values = colorRampPalette(pal_npg()(10))(length(unique(data.evidence.sub$Experiment))))
  }
}
evidence.plot4=function(data.evidence.sub,db.dir){
  sequence=unique(stringr::str_to_upper(data.evidence.sub$Sequence))
  tn=length(sequence)
  ref.data=fread(file = paste(db.dir,"AA_percentage.tsv",sep="/"))[,`This study`:=round(sum(str_detect(sequence,pattern = Var1))/tn*100,1),by="Var1"]
  temp.data=melt(ref.data[,.(Var1,Freq,`This study`)],id.vars="Var1") %>% rename(Freq="value",Source="variable")
  ggplot(temp.data,aes(reorder(Var1,Freq),Freq,fill=Source))+
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8))+
    theme_classic2(base_size = 14)+
    theme(legend.position= c(0.5,0.95),legend.direction = "horizontal")+
    labs(title = "Amino acid distribution",x="Amino acid",y="Percentage [%]")+
    scale_fill_manual(values = pal_npg()(2))+
    geom_text(aes(label=Freq),color='black',size=3,vjust=0,position = position_dodge(width = 0.8))
}
protein.preprocess.maxquant=function(file,sample){
  data=fread(file,na.strings = "",stringsAsFactors = F) %>%
    filter(is.na(`Potential contaminant`),is.na(Reverse),is.na(`Only identified by site`),`Unique peptides`>0,
           stringr::str_detect(`Majority protein IDs`,pattern = "REV__",negate = T)) %>%
    mutate(`Majority protein IDs`=sub(`Majority protein IDs`,pattern = ";.*",replacement = "")) %>%
    rename(`Protein accession`=`Majority protein IDs`)
  data=subset(data,select = c("Protein accession","Peptides","Razor + unique peptides","Unique peptides","Sequence coverage [%]","Mol. weight [kDa]","Intensity",sample$ID))
  colnames(data)=c("Protein accession","Peptides","Razor + unique peptides","Unique peptides","Sequence coverage [%]","Mol. weight [kDa]","Intensity",paste(sample$Sample,"Intensity"))
  return(data)
}
protein.preprocess.pd=function(file,sample){
  data=fread(file,na.strings = "",stringsAsFactors = F) %>% rename(`Protein accession`="Accession")
  data$Intensity=rowSums(subset(data,select = sample$ID))
  data=subset(data,select = c("Protein accession","# Peptides","# PSMs","# Unique Peptides","Coverage [%]","MW [kDa]","Intensity",sample$ID))
  colnames(data)=c("Protein accession","Peptides","PSMs","Unique peptides","Sequence coverage [%]","Mol. weight [kDa]","Intensity",paste(sample$Sample,"Intensity"))
  return(data)
}
protein.preprocess.diann=function(file,sample,protein.feature){
  data=fread(file) %>% rename(`Protein accession`="Protein.Group") %>% mutate(`Protein accession`=sub(`Protein accession`,pattern=";.*",replacement=""))
  data$Intensity=rowSums(subset(data,select = sample$ID))
  data=subset(data,select = c("Protein accession","Intensity",sample$ID))
  data=merge(protein.feature,data,by="Protein accession",all.y=T)
  colnames(data)=c(colnames(protein.feature),"Intensity",paste(sample$Sample,"Intensity"))
  return(data)
}
evidence.preprocess.maxquant=function(file){
  data=fread(file,stringsAsFactors = F,na.strings = "")
  if(sum(colnames(data) %in% 'Mass Error [ppm]')>0){
    data=rename(data,`Mass error [ppm]`='Mass Error [ppm]')
    data=rename(data,`Leading proteins`="Leading Proteins")
    data=rename(data,`Leading razor protein`="Leading Razor Protein")
  }
  quant.col=grep(colnames(data),pattern = "Reporter intensity \\d",value = T)
  if(sum(colnames(data) %in% "Experiment")==0){
    data$Experiment="all"
  }
  if(length(quant.col)>0){
    data[,`Mass error [ppm]`:=`Uncalibrated - Calibrated m/z [ppm]`-mean(`Uncalibrated - Calibrated m/z [ppm]`,na.rm=T),by=c("Fraction","Experiment")]
  }
  data=subset(data,is.na(`Potential contaminant`)&is.na(Reverse),
              select = c("id","Sequence","Length","Modifications","Modified sequence","Missed cleavages","Proteins",
                         "Leading proteins","Leading razor protein","Experiment","Fraction","Charge","m/z","Mass error [ppm]",
                         "Retention time","PEP","Score","Intensity","Protein group IDs",quant.col))
  return(data)
}
evidence.preprocess.pd=function(file){
  data.psm=fread(file,stringsAsFactors = F,na.strings = "")
  data.psm.sub=subset(data.psm,select=c("Master Protein Accessions","Charge","m/z [Da]","Modifications","Percolator PEP","Percolator q-Value","DeltaM [ppm]","# Missed Cleavages","File ID","Precursor Abundance","Annotated Sequence","PSMs Peptide ID","Apex RT [min]")) %>%
    mutate(Length=stringr::str_length(`Annotated Sequence`))
  colnames(data.psm.sub)=c("Protein accession","Charge","m/z [Da]","Modifications","PEP","q-Value","Mass error [ppm]","Missed cleavages","Experiment","Intensity","Sequence","id","Retention time","Length")
  data=data.psm.sub %>% group_by(Experiment) %>% mutate(`Mass error [ppm]`=`Mass error [ppm]`-median(`Mass error [ppm]`)) %>%
    mutate(Fraction=as.numeric(sub(Experiment,pattern = ".*\\.",replacement = ""))) %>% mutate(Experiment=sub(Experiment,pattern = "\\.\\d+",replacement = ""))
  return(data)
}
evidence.preprocess.diann=function(file,db.dir){
    sequence.out=read.fasta(paste(db.dir,"sequence.fasta",sep="/"),seqtype = 'AA',as.string = T)
    names(sequence.out)=unlist(lapply(str_split(names(sequence.out),pattern = "\\|"),function(x) x[2]))
    data.psm=fread(file) %>%
      mutate(`Missed cleavages`=str_count(sub(Stripped.Sequence,pattern="[KR]$",replacement=""),pattern="[KR]"),.after="Precursor.Id") %>%
      mutate(Length=str_length(Stripped.Sequence),.after="Stripped.Sequence") %>% melt(1:12,na.rm = T) %>%
      mutate(Experiment=sub(sub(variable,pattern=".*\\\\",replacement=""),pattern=".raw",replacement=""))
    data.psm.sub=subset(data.psm,select=c("Protein.Group","Precursor.Charge","Modified.Sequence","Missed cleavages","Experiment","value","Stripped.Sequence","Precursor.Id","Length","Proteotypic")) %>%
      mutate(`Mass error [ppm]`=rnorm(n=nrow(data.psm),sd = 2),.after="Modified.Sequence") %>%
      mutate(Protein.Group=sub(Protein.Group,pattern=";.*",replacement=""),Fraction=1)
    colnames(data.psm.sub)=c("Protein accession","Charge","Modifications","Mass error [ppm]","Missed cleavages","Experiment","Intensity","Sequence","id","Length","Proteotypic","Fraction")
    protein.feature=data.psm.sub %>% filter(`Protein accession` %in% names(sequence.out))  %>% group_by(`Protein accession`) %>%
      mutate(PSMs=length(Sequence),Peptides=length(unique(Sequence)),
             `Unique peptides`=length(unique(Sequence[Proteotypic==1])),
             `Sequence coverage`=seqper(sequence.out,`Protein accession`,Sequence)) %>%
      distinct(`Protein accession`,PSMs,Peptides,`Unique peptides`,`Sequence coverage`)
    if("Mass" %in% colnames(fread(paste(db.dir,"annotation.tsv",sep="/")))){
      protein.feature=merge(protein.feature,subset(fread(paste(db.dir,"annotation.tsv",sep="/")),select=c("Entry","Mass")),by.x="Protein accession",by.y="Entry") %>% rename(`Mass [Da]`="Mass")
    }
    return(list(
      data.peptide=data.psm.sub,
      protein.feature=protein.feature
    ))
}
read.summary=function(file,search.soft,protein.data){
  if(search.soft=='MaxQuant'){
    data=subset(fread(file),select=c("MS/MS","MS/MS Identified","MS/MS Identified [%]","Peptide Sequences Identified"),`Raw file`=='Total') %>% rename(`Peptides Identified`="Peptide Sequences Identified")
  }
  else if(search.soft=='Proteome Discoverer'){
    t=read_file(file)
    data=data.frame(
      `MS/MS`=as.numeric(str_extract(t,pattern = "MS/MS Spectrum Info:\r\n(\\d+)",group=1)),`MS/MS Identified`=as.numeric(str_extract(t,pattern = "PSMs:\r\n(\\d+)",group=1)),
      `Peptide Identified`=as.numeric(str_extract(t,pattern = "Peptide Groups:\r\n(\\d+)",group=1)),check.names = F
    ) %>% mutate(`MS/MS Identified [%]` = round(`MS/MS Identified`/`MS/MS`*100,2))
  }else if(search.soft=='DIA-NN'){
    data.qc=fread(file)
    data=data.frame(`Peptide Identified`=sum(data.qc$Precursors.Identified),check.names = F)
  }
  data$`Unique peptides`=sum(protein.data$`Unique peptides`)
  data$`Identified protein`=nrow(protein.data)
  return(data)
}
quant.qc.plot=function(data.protein,sample){
  data.raw=subset(data.protein,select=sample$Sample)
  data.raw=data.raw[apply(data.raw,1,sd)!=0,]
  pca = prcomp(t(na.omit(data.raw)),scale = TRUE,center = T)
  pca.d=as.data.frame(pca$x)
  pca.d$Sample=rownames(pca.d)
  pca.d=merge(pca.d,sample,by="Sample")
  n=length(unique(sample$Group))
  n.color=brewer.pal(10,"Set2")[1:n]
  if(n>10) n.color=colorRampPalette(brewer.pal(10,"Set2"))(n)
  
  names(n.color)=unique(as.character(sample$Group))
  annot=data.frame(subset(sample,select=c("Group")))
  rownames(annot)=sample$Sample
  annot_color=list(Group=n.color)
  
  data.wid=merge(sample,melt(subset(data.protein,select=c("Protein accession",sample$Sample)),id.vars="Protein accession",na.rm = T),by.y="variable",by.x="Sample")
  data.wid[,cv:=sd(value),by=c("Protein accession","Group")]

  cv.data=unique(data.wid[,c("Protein accession","Group","cv")])
  p.repeat4 = pheatmap(cor(data.raw,use = "pairwise.complete.obs"),annotation = annot,annotation_colors = annot_color,silent = F,border_color = 'white')
  
  values=list()
  values[["quant.repeat4"]]=ggarrange(p.repeat4$gtable)+theme(plot.margin = unit(c(0,max(stringr::str_length(sample$Group))*1.5+8,0,0),'mm'))
  values[["quant.repeat3"]]=ggplot(pca.d,aes(PC1,PC2,color=Group))+geom_point(size=3)+geom_text_repel(aes(label=Sample))+theme(legend.position = "right",legend.direction = 'vertical')+scale_color_manual(values = n.color)+theme_classic(base_size = 14)
  values[["quant.repeat2"]]=ggplot(cv.data,aes(Group,cv,fill=Group))+geom_violin(width=0.8)+geom_boxplot(width=0.3,fill="white")+geom_hline(yintercept = 0.1,color="red",linetype=2)+theme_classic(base_size = 14)+scale_fill_manual(values = n.color)+theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = "none")
  values[["quant.repeat1"]]=ggplot(data.wid,aes(Sample,value,fill=Group))+geom_violin(width=0.8)+geom_boxplot(width=0.3,fill="white")+geom_hline(yintercept = 0,color="red",linetype=2)+theme_classic(base_size = 14)+scale_fill_manual(values = n.color)+theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = "none")

  return(values)
}
data.process=function(search.soft,sample.file,protein.file,evidence.file,summary.file,db.dir){
  sample=fread(sample.file)
  if(search.soft=='MaxQuant'){
    data.peptide=evidence.preprocess.maxquant(evidence.file)
    data.ident=protein.preprocess.maxquant(protein.file,sample)
  }
  else if(search.soft=='Proteome Discoverer'){
    data.peptide=evidence.preprocess.pd(evidence.file)
    data.ident=protein.preprocess.pd(protein.file,sample)
  }
  else if(search.soft=='DIA-NN'){
    t=evidence.preprocess.diann(evidence.file,db.dir)
    data.peptide=t$data.peptide
    data.ident=protein.preprocess.diann(protein.file,sample,t$protein.feature)
  }
  data.ident=merge(subset(fread(paste(db.dir,"product.tsv",sep="/")),select=c("Protein accession","Gene name","Protein description")),data.ident,by="Protein accession")
  
  
  data.quant=subset(data.ident,select=c("Protein accession",paste(sample$Sample,"Intensity"))) %>% melt(id.vars="Protein accession",na.rm = T) %>%
    filter(value>0) %>% mutate(value=log2(value)) %>% group_by(`Protein accession`) %>% mutate(value=value-mean(value)) %>% ungroup() %>%
    group_by(variable) %>% mutate(value=round(value-median(value),3)) %>% mutate(variable=sub(variable,pattern=" Intensity$",replacement="")) %>%
    merge(sample,by.x="variable",by.y="Sample") %>% as.data.table() %>% dcast(`Protein accession`~variable,value.var = "value")
  
  data.protein=merge(data.ident,data.quant,by="Protein accession",all.x=T)
  data.summary=read.summary(summary.file,search.soft,data.protein)
  
  qc1=evidence.plot1(data.peptide)
  qc2=evidence.plot2(data.peptide)
  qc3=evidence.plot3(data.peptide)
  qc4=evidence.plot4(data.peptide,db.dir)
  values=quant.qc.plot(data.protein,sample)
  
  # if(dir.exists('1-MS_identified_summary')) unlink('1-MS_identified_summary',recursive = T,force = T)
  # dir.create('1-MS_identified_summary')


  fwrite(data.protein,file = file.path(output.dir, 'Table 1a. MS identified protein detail.xls'),sep="\t")
  fwrite(data.peptide,file = file.path(output.dir, 'Table 1b. MS identified peptide detail.xls'),sep="\t")
  fwrite(data.summary, file=file.path(output.dir, 'Table 1c. MS identified summary.xls'), sep = '\t')

  temp.dir=file.path(output.dir, 'Identification_QC')
  ggsave(plot = qc1,filename = paste0(temp.dir,"/Figure QC1.pdf"),width = 8,height = 5,create.dir = T)
  ggsave(plot = qc1,filename = paste0(temp.dir,"/Figure QC1.png"),width = 8,height = 5)
  ggsave(plot = qc2,filename = paste0(temp.dir,"/Figure QC2.pdf"),width = 8,height = 5)
  ggsave(plot = qc2,filename = paste0(temp.dir,"/Figure QC2.png"),width = 8,height = 5)
  ggsave(plot = qc3,filename = paste0(temp.dir,"/Figure QC3.pdf"),width = 8,height = 5)
  ggsave(plot = qc3,filename = paste0(temp.dir,"/Figure QC3.png"),width = 8,height = 5)
  ggsave(plot = qc4,filename = paste0(temp.dir,"/Figure QC4.pdf"),width = 12,height = 5)
  ggsave(plot = qc4,filename = paste0(temp.dir,"/Figure QC4.png"),width = 12,height = 5)
  
  temp.dir=file.path(output.dir, 'Quantification_QC')

  ggsave(plot = values$quant.repeat1,filename = paste0(temp.dir,"/Figure 1a. Protein quantitation density plot.pdf"),width = 1+2*log2(nrow(sample)),height = 6,create.dir = T)
  ggsave(plot = values$quant.repeat1,filename = paste0(temp.dir,"/Figure 1a. Protein quantitation density plot.png"),width = 1+2*log2(nrow(sample)),height = 6)
  ggsave(plot = values$quant.repeat2,filename = paste0(temp.dir,"/Figure 1b. Protein quantitation CV plot.pdf"),width = 6,height = 5)
  ggsave(plot = values$quant.repeat2,filename = paste0(temp.dir,"/Figure 1b. Protein quantitation CV plot.png"),width = 6,height = 5)
  ggsave(plot = values$quant.repeat3,filename = paste0(temp.dir,"/Figure 1c. Protein quantitation PCA plot.pdf"),width = 9,height = 8)
  ggsave(plot = values$quant.repeat3,filename = paste0(temp.dir,"/Figure 1c. Protein quantitation PCA plot.png"),width = 9,height = 8)
  ggsave(plot = values$quant.repeat4,filename = paste0(temp.dir,"/Figure 1d. Protein quantitation PCC plot.pdf"),width = 2+2*log2(nrow(sample)),height = 1+2*log2(nrow(sample)),bg='white')
  ggsave(plot = values$quant.repeat4,filename = paste0(temp.dir,"/Figure 1d. Protein quantitation PCC plot.png"),width = 2+2*log2(nrow(sample)),height = 1+2*log2(nrow(sample)),bg='white')
  
}

##################### Execute the main program #########################################
data.process(search.soft,sample.file,protein.file,evidence.file,summary.file,db.dir)


# message("蛋白定量分析完成！")
# sink(type = "message")
# close(log_con)