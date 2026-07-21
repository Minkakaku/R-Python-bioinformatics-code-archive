library(data.table)
library(dplyr)
library(seqinr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(ggrepel)

##################### parameter #########################################
search.soft <- snakemake@params[["search_soft"]]
sample.file <- snakemake@input[["sample_file"]]
ptm.file <- snakemake@input[["ptm_file"]]
db.dir <- snakemake@input[["db_dir"]]
motified.list.file <- snakemake@params[["motified_list_file"]]
modified.type <- snakemake@params[["modified_type"]]
normalized_file <- snakemake@params[["normalized_file"]]
output.dir <- snakemake@params[["output_dir"]]

##################### function defined #########################################
ptm.file.preprocess.maxquant = function(ptm_file,sample,product,normalized_file){
  data=fread(file=ptm_file,na.strings = "",stringsAsFactors = F,integer64 = "numeric") %>%
    rename(`Protein accession`= Protein) %>% filter(is.na(`Potential contaminant`),is.na(Reverse)) %>%
    rename(`Localization probability` = `Localization prob`) %>% filter(`Localization probability` > 0.75) %>%
    filter(grep(`Protein accession`,pattern = "REV__",invert = T))
  l=c("Protein accession","Position","Amino acid","Localization probability","PEP","Score","Intensity",
      grep(colnames(data),pattern = "Probabilities",value = T),"Charge","Mass error [ppm]",sample$ID)
  m=c("Protein accession","Position","Amino acid","Localization probability","PEP","Score","Intensity",
      "Modified sequence","Charge","Mass error [ppm]",paste(sample$Sample,"Intensity"))
  data=subset(data,select=l);colnames(data)=m
  data.quant=melt(subset(data,select=c("Protein accession","Position",paste(sample$Sample,"Intensity"))),c("Protein accession","Position"),na.rm = T) %>%
    filter(value>0) %>% mutate(value=log2(value)) %>%
    group_by(variable) %>% mutate(value=value-median(value)+15) %>% ungroup() %>%
    group_by(`Protein accession`,`Position`) %>% mutate(value=value-mean(value)) %>% ungroup() %>%
    mutate(variable=sub(variable,pattern=" Intensity$",replacement="")) %>% as.data.table() %>%
    dcast(`Protein accession`+`Position`~variable,value.var = "value")
  
  if(file.exists(normalized_file)){
    ptm.data.long = melt(subset(data.quant,select=c("Protein accession","Position",sample$Sample)),id.vars = c("Protein accession","Position"),na.rm = T) %>% mutate(variable = as.character(variable))
    protein.data.long = subset(fread(normalized_file,na.strings = ""),select=c("Protein accession",sample$Sample)) %>% melt(id.vars = "Protein accession",na.rm = T) %>% mutate(variable = as.character(variable))
    ptm.data.nor = merge(ptm.data.long,protein.data.long,by = c("Protein accession","variable"),suffixes = c(".ptm",".pro"),all.x=T) %>% mutate(value = round(value.ptm-value.pro,3))
    data.quant = dcast(ptm.data.nor,`Protein accession`+Position ~ variable,value.var = "value")
  }
  
  data = merge(merge(product,data,by="Protein accession"),data.quant,by=c("Protein accession","Position"))
  return(data)
}
ptm.file.preprocess.dia = function(ptm_file,sample,product,sequence.out,id,normalized_file){
  id=paste0("[A-Z]\\(",id,"\\)")
  data=fread(file=ptm_file,na.strings = "",stringsAsFactors = F,integer64 = "numeric") %>% dplyr::rename(`Protein accession` = "Protein.Group") %>%
    mutate(`Protein accession`= sub(`Protein accession`,pattern=";.*",replacement="")) %>% filter(`Protein accession` %in% names(sequence.out))
  data.ptm=subset(data,grepl(Modified.Sequence,pattern=id))
  data.sub=data.ptm %>% distinct(`Protein accession`,Stripped.Sequence)
  pep.position=rbindlist(lapply(split(data.sub,data.sub$`Protein accession`),function(x){
    sitelocation(sequence.out[[x$`Protein accession`[1]]][1],x$Stripped.Sequence)
  }),idcol = "Protein accession")
  data.ptm.mod=data.ptm %>% distinct(Modified.Sequence) %>% mutate(seq=gsub(gsub(unique(Modified.Sequence),pattern=id,replacement="#"),pattern="\\(UniMod:\\d+\\)",replacement=""))
  site.position=merge(data.ptm.mod,sitelocation0(unique(data.ptm.mod$seq),'#'),by="seq") %>% dplyr::select(Modified.Sequence,start) %>% dplyr::rename(Peptide.position="start")
  data.temp=merge(data.ptm,pep.position,by=c("Protein accession","Stripped.Sequence")) %>% merge(site.position,by="Modified.Sequence") %>%
    mutate(Position=start+Peptide.position-1,`Amino acid`=str_sub(Stripped.Sequence,Peptide.position,Peptide.position))
  data.temp$Intensity=rowSums(subset(data.temp,select=grep(pattern = ".raw",colnames(data.temp),value = T)),na.rm = T)
  data.temp=data.temp %>% arrange(-Intensity) %>% distinct(`Protein accession`,Position,.keep_all = T)
  data=subset(data.temp,select=c("Protein accession","Amino acid","Position","Modified.Sequence","Intensity",sample$ID))
  colnames(data)=c("Protein accession","Amino acid","Position","Modified sequence","Intensity",paste(sample$Sample,"Intensity"))
  data.quant=melt(subset(data,select=c("Protein accession","Position",paste(sample$Sample,"Intensity"))),c("Protein accession","Position"),na.rm = T) %>%
    filter(value>0) %>% mutate(value=log2(value)) %>%
    group_by(variable) %>% mutate(value=value-median(value)+15) %>% ungroup() %>%
    group_by(`Protein accession`,`Position`) %>% mutate(value=value-mean(value)) %>% ungroup() %>%
    mutate(variable=sub(variable,pattern=" Intensity$",replacement="")) %>% as.data.table() %>%
    dcast(`Protein accession`+`Position`~variable,value.var = "value")

  if(file.exists(normalized_file)){
    ptm.data.long = melt(subset(data.quant,select=c("Protein accession","Position",sample$Sample)),id.vars = c("Protein accession","Position"),na.rm = T) %>% mutate(variable = as.character(variable))
    protein.data.long = subset(fread(normalized_file,na.strings = ""),select=c("Protein accession",sample$Sample)) %>% melt(id.vars = "Protein accession",na.rm = T) %>% mutate(variable = as.character(variable))
    ptm.data.nor = merge(ptm.data.long,protein.data.long,by = c("Protein accession","variable"),suffixes = c(".ptm",".pro"),all.x=T) %>% mutate(value = round(value.ptm-value.pro,3))
    data.quant = dcast(ptm.data.nor,`Protein accession`+Position ~ variable,value.var = "value")
  }
  
  data = merge(merge(product,data,by="Protein accession"),data.quant,by=c("Protein accession","Position"))
  return(data)
}
sitelocation=function(sequence,peptide){
  o=str_locate_all(sequence,peptide)
  o=lapply(o,function(x)data.frame(x)[1,])
  names(o)=peptide
  rbindlist(o,idcol = "Stripped.Sequence")
}
sitelocation0=function(sequence,peptide){
  o=str_locate_all(sequence,peptide)
  o=lapply(o,data.frame)
  names(o)=sequence
  rbindlist(o,idcol = "seq")
}
##################### Execute the main program #########################################
motified.list=fread(motified.list.file)
product=subset(fread(paste(db.dir,"product.tsv",sep="/")),select=c("Protein accession","Gene name","Protein description"))
sample=fread(sample.file)

## data quantification
if(search.soft=='MaxQuant') data=ptm.file.preprocess.maxquant(ptm.file,sample,product,normalized_file)
if(search.soft=='DIA-NN'){
  sequence.out <- read.fasta(paste(db.dir,"sequence.fasta",sep="/"),seqtype = 'AA',as.string = T)
  names(sequence.out)=unlist(lapply(str_split(names(sequence.out),pattern = "\\|"),function(x) x[2]))
  data=ptm.file.preprocess.dia(ptm.file,sample,product,sequence.out,motified.list[Name==modified.type]$Accession,normalized_file)
}

summary.table <- data.frame(
  `Modified peptides` = length(unique(data$`Modified sequence`)),
  `Identified protein` = length(unique(data$`Protein accession`)),
  `Identified sites` = nrow(data),check.names = F
)

fwrite(summary.table,file=file.path(output.dir, 'Table 1b. MS identified modification summary.xls'),sep="\t")

## qc plot
data.raw1 = t(na.omit(subset(data,select=sample$Sample)));data.raw1 <- data.raw1[, which(apply(data.raw1, 2, var) != 0)]
pca = prcomp(data.raw1,scale = TRUE,center = T);pca.d=as.data.frame(pca$x);pca.d$Sample=rownames(pca.d);pca.d=merge(pca.d,sample,by="Sample")
n=length(unique(sample$Group));if(n<=8){n.color=brewer.pal(8,"Set2")[1:n]}else{n.color=colorRampPalette(brewer.pal(8,"Set2"))(n)}

names(n.color)=unique(as.character(sample$Group));annot=data.frame(subset(sample,select=c("Group")))
rownames(annot)=sample$Sample;annot_color=list(Group=n.color)

data.wid=subset(data,select=c("Protein accession","Amino acid","Position",sample$Sample))
data.wid=melt(data.wid,id.vars=c("Protein accession","Amino acid","Position"),na.rm = T)
data.wid=merge(sample,data.wid,by.y="variable",by.x="Sample")
data.wid[,cv:=sd(value),by=c("Protein accession","Amino acid","Position","Group")] 
cv.data=unique(data.wid[,c("Protein accession","Amino acid","Position","Group","cv")])            
p.repeat4 = pheatmap(cor(subset(data,select=sample$Sample),use = "pairwise.complete.obs"),annotation = annot,annotation_colors = annot_color,silent = F,border_color = 'white')

p4=ggarrange(p.repeat4$gtable)+theme(plot.margin = unit(c(0,max(stringr::str_length(sample$Group))*1.5+8,0,0),'mm'))
p3=ggplot(pca.d,aes(PC1,PC2,color=Group))+geom_point(size=3)+geom_text_repel(aes(label=Sample))+theme(legend.position = "right",legend.direction = 'vertical')+scale_color_manual(values = n.color)+theme_classic(base_size = 14)
p2=ggplot(cv.data[cv<1],aes(Group,cv,fill=Group))+geom_violin(width=0.8)+geom_boxplot(width=0.3,fill="white")+geom_hline(yintercept = 0,color="red",linetype=2)+theme_classic(base_size = 14)+scale_fill_manual(values = n.color)+theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = "none")
p1=ggplot(data.wid[abs(value)<5,],aes(Sample,value,fill=Group))+geom_violin(width=0.8)+geom_boxplot(width=0.3,fill="white")+geom_hline(yintercept = 0,color="red",linetype=2)+theme_classic(base_size = 14)+scale_fill_manual(values = n.color)+theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = "none")

## output

fwrite(data,file=file.path(output.dir, 'Table 1a. MS identified modification site detail.xls'),sep="\t")

temp.dir=file.path(output.dir, 'Quantification_QC')
ggsave(plot = p1,filename = paste0(temp.dir,"/Figure 1a. Protein quantitation density plot.pdf"),width = 1+2*log2(nrow(sample)),height = 6,create.dir = T)
ggsave(plot = p1,filename = paste0(temp.dir,"/Figure 1a. Protein quantitation density plot.png"),width = 1+2*log2(nrow(sample)),height = 6)
ggsave(plot = p2,filename = paste0(temp.dir,"/Figure 1b. Protein quantitation CV plot.pdf"),width = 6,height = 5)
ggsave(plot = p2,filename = paste0(temp.dir,"/Figure 1b. Protein quantitation CV plot.png"),width = 6,height = 5)
ggsave(plot = p3,filename = paste0(temp.dir,"/Figure 1c. Protein quantitation PCA plot.pdf"),width = 9,height = 8)
ggsave(plot = p3,filename = paste0(temp.dir,"/Figure 1c. Protein quantitation PCA plot.png"),width = 9,height = 8)
ggsave(plot = p4,filename = paste0(temp.dir,"/Figure 1d. Protein quantitation PCC plot.pdf"),width = 2+2*log2(nrow(sample)),height = 1+2*log2(nrow(sample)),bg='white')
ggsave(plot = p4,filename = paste0(temp.dir,"/Figure 1d. Protein quantitation PCC plot.png"),width = 2+2*log2(nrow(sample)),height = 1+2*log2(nrow(sample)),bg='white')
