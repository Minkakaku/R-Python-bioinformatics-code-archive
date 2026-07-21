library(data.table)
library(tibble)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggsci)
library(stringr)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(impute)
library(ellipse)
library(xlsx)
library(reshape2)

############################## Parameters ######################################
workflow = snakemake@params[["workflow"]]
organism = snakemake@params[["organism"]]
data.file = snakemake@input[["data_file"]]
istd.file = snakemake@input[["istd_file"]]
pos.png = snakemake@input[["pos_png"]]
neg.png = snakemake@input[["neg_png"]]
sample_info = snakemake@input[["sample_info"]]
db_dir = snakemake@params[["db_dir"]]
pca.scale = snakemake@params[["pca_scale"]]
add.ellipse = TRUE
results_dir = snakemake@params[["results_dir"]]

############################ Functions defined #################################
convert_integer64_to_numeric <- function(dt){
  for (i in 1:ncol(dt)) {
    if (inherits(dt[[i]], "integer64")) {
      dt[[i]] <- as.numeric(dt[[i]])
    }}
  return(dt)
}
get_ellipse <- function(df, level = 0.90) { # Generate ellipse data for each group
  cov_matrix <- cov(df[, 1:2])
  center <- colMeans(df[, 1:2])
  ellipse_points <- as.data.frame(ellipse(cov_matrix, centre = center, level = level))
  return(ellipse_points)
}
diff.classy.plot1=function(data, font = 14,color.pal="Set3",  angle = 75, margin =15, margin2 = 10, type ="Metabolomics", label.font=4){
  if(type=="Metabolomics"){
    n = brewer.pal.info[color.pal,"maxcolors"]
    if (length(unique(data$HMDB_SuperClass))>n){
      colors = colorRampPalette(ggsci.colors[[color.pal]][1:n])(length(unique(data$HMDB_SuperClass)))
    }else{
      colors = ggsci.colors[[color.pal]][1:length(unique(data$HMDB_SuperClass))]
    }
    names(colors)=unique(data$HMDB_SuperClass)
    data=data %>% arrange(-Freq)
    data$HMDB_SuperClass = factor(data$HMDB_SuperClass, levels=unique(data$HMDB_SuperClass))
    p1 = ggplot(data,aes(reorder(HMDB_Class,-Freq),Freq,fill=HMDB_SuperClass))+geom_col()+facet_grid(.~HMDB_SuperClass,scales = 'free',space = 'free')+scale_fill_manual(values = colors)+
      theme_classic(base_size = font)+scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      theme(strip.text = element_blank(),axis.text.x=element_text(angle = angle,hjust=1),axis.text = element_text(color='black'), plot.margin=unit(c(margin2,5,5,margin),'mm'),
            legend.position = "top",legend.title.position = "top",legend.background=element_blank(), legend.box.background=element_blank())+
      geom_text(aes(label=Freq),size=label.font, vjust=0)+labs(x = "HMDB Category",y="Number of compounds")+
      guides(fill=guide_legend(ncol=3))
  }
  if(type=="Lipidomics"){
    colors = ggsci.colors[[color.pal]][1:length(unique(data$Categories))]
    names(colors)=unique(data$Categories)
    p1 = ggplot(data, aes(reorder(Name,Freq),Freq,fill=Categories))+geom_col()+
      facet_grid(.~Categories,scales = 'free',space = 'free')+scale_fill_manual(values = colors)+theme_classic(base_size = font)+scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      theme(strip.text = element_blank(),axis.text.x=element_text(angle = angle,hjust=1),axis.text = element_text(color='black'), plot.margin=unit(c(margin2,5,5,margin),'mm'))+
      geom_text(aes(label=Freq),size=label.font, vjust=0)+theme(legend.background=element_blank(), legend.box.background=element_blank())+labs(x = "Lipid category",y="Number of compounds")+
      theme(legend.position = "top",legend.title.position = "top")+
      guides(fill=guide_legend(ncol = 2))
  }
  return(p1)
}
diff.classy.plot2=function(data, font=14, color.pal="Set2", angle = 60, margin = 15, margin2 = 10,label.font=4){
  colors = ggsci.colors[[color.pal]][1:length(unique(data$Pathway.level.1))]
  names(colors)=unique(data$Pathway.level.1)
  p2 = ggplot(data,
              aes(reorder(Pathway.level.2,Freq),Freq,fill=Pathway.level.1))+geom_col()+
    facet_grid(.~Pathway.level.1,scales = 'free',space = 'free')+scale_fill_manual(values = colors)+theme_classic(base_size = 14)+scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme(strip.text = element_blank(),axis.text.x=element_text(angle = angle,hjust=1),axis.text = element_text(color='black'), plot.margin=unit(c(margin2,5,5,margin),'mm'))+
    geom_text(aes(label=Freq),vjust=0,size=label.font)+theme(legend.position = "top",legend.title.position ="top",legend.background=element_blank(), legend.box.background=element_blank())+labs(x = "KEGG pathway",y="Number of compounds", fill ="Pathway Level 1")+
    guides(fill=guide_legend(ncol = 2))
  return(p2)
}
ggsci.colors <- list(
  "Set1"=brewer.pal(9,"Set1"),
  "Set2"=brewer.pal(8,"Set2"),
  "Set3"=brewer.pal(12,"Set3"),
  "Pastel1"=brewer.pal(9,"Pastel1"),
  "Pastel2"=brewer.pal(8,"Pastel2"),
  "Paired"=brewer.pal(12,"Paired"),
  "Dark2"=brewer.pal(8,"Dark2"),
  "Accent"=brewer.pal(8,"Accent"),
  "npg"=pal_npg()(10),
  "aaas"=pal_aaas()(10),
  "nejm"=pal_nejm()(8),
  "lancet"=pal_lancet()(9),
  "jama"=pal_jama()(7),
  "jco"=pal_jco()(10),
  "ucscgb"=pal_ucscgb()(26),
  "d3_category10"=pal_d3("category10")(10),
  "d3_category20"=pal_d3("category20")(20),
  "d3_category20c"=pal_d3("category20c")(20),
  "igv"=pal_igv()(51),
  "uchicago"=pal_uchicago()(9),
  "rickandmorty"=pal_rickandmorty()(12),
  "simpsons"=pal_simpsons()(16)
)
########################### Data Quantification ################################
if (workflow =="Metabolomics"){
  ## 读取ISTD内标文件
  a <- str_split(istd.file, pattern="\\.")[[1]][2]
  if (a == "xlsx"){
    istd.data=read_xlsx(istd.file) %>% dplyr::rename(`Metabolite name`="Molecule Name") %>% as.data.frame()
  }else{
    istd.data=fread(file=istd.file) %>% dplyr::rename(`Metabolite name`="Molecule Name") %>% as.data.frame()
  }
  ## 读取Sample文件（包含样本名修改以及样本删减等情况）
  sample=fread(file=sample_info)
  if(all(sample$Sample == "-") == T){
    sample = sample %>% mutate(Sample = ID)
  } else {
    sample = sample 
  }
  
  istd.data[,sample$ID]=convert_integer64_to_numeric(istd.data[,sample$ID])
  ## 读取data文件
  a <- str_split(data.file, pattern="\\.")[[1]][2]
  if (a == "xlsx"){
    data=readxl::read_excel(data.file) %>% as.data.frame()
  }else{
    data=fread(file=data.file) %>% as.data.frame()
  }
  data[,sample$ID]=convert_integer64_to_numeric(data[,sample$ID])
  data=data %>% dplyr::rename(`Metabolite name` = "name",`Adduct type`="adduct",`Average Mz`="mz",
                              `Average Rt(s)`="rt",Score="score",Source ="source",Level="Levels",
                              Inchikey="inchikey",`HMDB ID`="hmdb_accession",HMDB_SuperClass="hmdb_super_class",HMDB_Class="hmdb_class",HMDB_SubClass="hmdb_sub_class") %>% 
    mutate(HMDB_SuperClass=case_when(HMDB_SuperClass==""~"Unclassified",TRUE~HMDB_SuperClass),
           HMDB_Class=case_when(HMDB_Class==""~"Unclassified",TRUE~HMDB_Class),
           HMDB_SubClass=case_when(HMDB_SubClass==""~"Unclassified",TRUE~HMDB_SubClass)) %>% distinct()
  data$`Alignment ID`=paste("Meta",1:nrow(data),sep="_")
  hmdb=fread(paste0(db_dir,"/hmdb.txt"))
  kegg.data=fread(paste0(db_dir,"/KEGG_compound_pathway_annotation.tsv")) %>% 
    filter(`Pathway ID` %in% sub(fread(paste0("https://rest.kegg.jp//list/pathway/",organism),header = F)$V1,pattern=organism,replacement="map"))

  tmp = data %>% dplyr::select(`Alignment ID`,sample$ID)
  tmp = merge(melt(tmp,measure.vars=sample$ID, variable.name="ID"),sample,by="ID") %>%
    group_by(`Alignment ID`) %>% mutate(Intensity = sum(value,na.rm = T)) %>% ungroup() %>%
    dcast(`Alignment ID`+Intensity~Sample)
  meta.data= merge(data %>% dplyr::select(-sample$ID),tmp, by="Alignment ID") 
  data.o=istd.data %>% melt(measure.vars = sample$ID, variable.name="ID") %>% merge(sample,by="ID") 

  ## 数据标准化
  data.nor=melt(meta.data, measure.vars = sample$Sample,id.vars="Alignment ID") %>% filter(value>0) %>% mutate(value=log2(value)) %>%
    group_by(`Alignment ID`) %>% filter(length(`Alignment ID`)>nrow(sample)*2/3) %>% ungroup() %>% merge(sample, by.x = "variable", by.y= "Sample")
}

if(workflow == "Lipidomics"){
  ## 读取data文件
  a <- str_split(data.file, pattern="\\.")[[1]][2]
  if (a == "xlsx"){
    data=readxl::read_excel(data.file) %>% as.data.frame()
  }else{
    data=fread(file=data.file) %>% as.data.frame()
  }
  inlist=fread(paste0(db_dir,"/internal.standard.for.normalization.txt"))
  inlist=inlist$`Internal Standard`
  lion.anno=fread(paste0(db_dir, "/LION_data.txt"))
  lipid.category=fread(paste0(db_dir,"/Lipid_category.txt"))
  
  data=data %>% dplyr::rename(`Metabolite name` = "Molecule Name", Abbreviation = "Molecule List Name")
  sample=fread(file=sample_info)
  if(all(sample$Sample == "-") == T){
    sample = sample %>% mutate(Sample = ID)
  } else {
    sample = sample 
  }
  
  data=data %>% mutate(`Alignment ID`=paste("Lipid",1:length(`Metabolite name`),sep="_"))
  meta.data=merge(lipid.category, merge(lion.anno, data, by ="Metabolite name", all.y=T), by="Abbreviation") %>% filter(!`Metabolite name` %in% inlist)
  
  data = merge(melt(data,measure.vars=sample$ID, variable.name="ID"),sample,by="ID") %>%
    dcast(`Metabolite name`+Abbreviation+`Precursor Mz`+`Product Mz`+`Explicit Retention Time`~Sample)
  meta.data = merge(melt(meta.data,measure.vars=sample$ID,variable.name="ID"),sample, by="ID") %>%
    group_by(`Alignment ID`) %>% mutate(Intensity = sum(value,na.rm = T)) %>% ungroup() %>%
    dcast(`Metabolite name`+`Alignment ID`+`中文名`+Abbreviation+Categories+Name+`LION name`+`LION ID`+ChEBI+LIPID.MAPS.ID+SwissLipids.ID+HMDB.ID+KEGG.ID+PubChem.CID+`Precursor Mz`+`Product Mz`+`Explicit Retention Time`+Intensity~Sample)
  data.o=merge(sample,data %>% filter(`Metabolite name` %in% inlist) %>% reshape2::melt(measure.vars=sample$Sample),by.y="variable",by.x="Sample")
  
  ## 数据标准化
  data.o=merge(sample, data %>% filter(`Metabolite name` %in% inlist) %>% melt(measure.vars = sample$Sample), by.y="variable",by.x="Sample") %>%
    merge(lipid.category,by="Abbreviation") %>% group_by(Abbreviation) %>% mutate(Nor.value=log2(value)-mean(log2(value))) %>% ungroup() %>% distinct()
  data.l=merge(sample,melt(meta.data, measure.vars=sample$Sample) %>% filter(value>0) %>% mutate(value=log2(value)), by.x="Sample", by.y="variable") %>% group_by(`Metabolite name`) %>% filter(length(value)>nrow(sample)*2/3) %>% ungroup %>% as.data.table()
  data.nor=data.l %>% merge(subset(data.o, select=c("Abbreviation", "Sample","Nor.value")), by=c("Abbreviation","Sample"), all.x=T) %>% mutate(Nor.value=case_when(is.na(Nor.value)~0, TRUE~Nor.value)) %>% mutate(value=value-Nor.value) %>% dplyr::rename(variable ="Sample")
}

## KNN缺失值填补
int=15
set.seed(123)
data.impute=data.nor %>% group_by(variable) %>% mutate(value=value-median(value)+int) %>% ungroup() %>% reshape2::dcast(`Alignment ID`~variable,value.var = "value") %>% column_to_rownames("Alignment ID")
data.impute=merge(melt(meta.data, measure.vars = sample$Sample) %>% dplyr::select(-value),impute.knn(as.matrix(data.impute))$data %>% melt(),by.x=c("Alignment ID","variable"),by.y=c("Var1","Var2")) %>%
  mutate(value = round(value,3))
data.long=merge(sample,data.impute,by.x="Sample",by.y="variable") %>% as.data.table()
data.wid=dcast(data.long,`Alignment ID`~Sample,value.var = "value")
data.matrix=as.matrix(subset(data.wid,select=c(sample$Sample)))
rownames(data.matrix)=data.wid$`Alignment ID`
if (workflow=="Metabolomics"){
  data.final.table=reshape2::dcast(data.long,`Metabolite name`+HMDB_SuperClass+HMDB_Class+HMDB_SubClass+`HMDB ID`+Inchikey+`Alignment ID`+`Average Rt(s)`+`Average Mz`+`Adduct type`+Score+Source+Level+Intensity~Sample) %>% distinct()
}
if (workflow=="Lipidomics"){
  data.final.table=reshape2::dcast(data.long,`Alignment ID`+`Metabolite name`+ Name + `中文名` + Abbreviation + Categories + `LION ID` + `LION name` + ChEBI + LIPID.MAPS.ID+SwissLipids.ID+HMDB.ID+KEGG.ID+PubChem.CID+`Precursor Mz`+`Product Mz`+`Explicit Retention Time`+Intensity~Sample) %>% distinct()
}

## 绘制QC图
data.pca=prcomp(t(data.matrix),scale. = pca.scale,center = TRUE)
n=length(unique(sample$Group))
if(n<3){
  n.color=brewer.pal(n,"Set1")[1:n]
}else if(n<=brewer.pal.info["Set1","maxcolors"]){
  n.color=brewer.pal(n,"Set1")
}else{
  n.color=colorRampPalette(brewer.pal(brewer.pal.info["Set1","maxcolors"],"Set1"))(n)
}
names(n.color)=unique(as.character(sample$Group))
qc.plot=ggplot(data.o,aes(Sample,value,fill=Group))+geom_col()+facet_grid(`Metabolite name`~Group,scales = 'free',space = 'free_x')+
  scale_fill_manual(values=n.color)+theme_test(base_size = 14)+
  theme(axis.text.x=element_text(angle = 60,hjust=1),legend.position = 'none',strip.text.y=element_text(angle = 0))

annot=data.frame(Group = sample$Group,row.names = sample$Sample)
annot_color=list(Group=n.color)
names(annot_color$Group)=unique(sample$Group)
quant.qc1=ggplot(data.long,aes(Sample,value,fill=Group))+geom_violin()+
    geom_boxplot(fill='white',width=0.3,outlier.shape = NA)+theme_classic(base_size = as.numeric(12))+
    theme(axis.text.x = element_text(angle = 45,hjust=1))+scale_fill_manual(values=n.color)

if(min(table(sample$Group)) < 2) {
  quant.qc2 = ggplot()+theme(panel.background = element_blank())
} else {
   quant.qc2=ggplot(data.long %>% group_by(`Alignment ID`,Group) %>% mutate(CV=sd(value)/abs(mean(value))) %>% distinct(`Alignment ID`,Group,CV))+
      stat_ecdf(aes(CV,color=Group),linewidth=1)+theme_classic(base_size = as.numeric(12))+scale_x_continuous(expand = expansion(mult = 0))+scale_color_manual(values=n.color)

}
if(add.ellipse==TRUE){
  if(min(table(sample$Group))>3){
    quant.qc3=ggplot(merge(data.pca$x,sample,by.y="Sample",by.x="row.names"),aes(PC1,PC2,color=Group))+
        geom_point(size=3)+theme_classic(base_size = as.numeric(12))+
        geom_text_repel(aes(label=Row.names))+scale_color_manual(values=n.color)+
        stat_ellipse(aes(fill=Group),type = "norm",geom = "polygon",alpha=0.2, color = NA, level = 0.95)+
        labs(x=sprintf('PC1 (%.2f%%)',summary(data.pca)[6]$importance[2,1]*100),
            y=sprintf('PC2 (%.2f%%)',summary(data.pca)[6]$importance[2,2]*100))
  }else if(min(table(sample$Group))!= 1){
    df=merge(data.pca$x,sample,by.x = "row.names",by.y = "Sample") %>% dplyr::rename(Sample = Row.names) %>% as.data.frame() %>% dplyr::select(PC1,PC2, everything())
    # Generate ellipse data for each group
    ellipses <- do.call(rbind, lapply(split(df, df$Group), function(group_df) {
      e <- get_ellipse(group_df)
      e$Group <- unique(group_df$Group)
      return(e)
    }))
    quant.qc3=ggplot(df, aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_classic(base_size = as.numeric(12))+
      geom_text_repel(aes(label=Sample),hjust = 0.5,vjust = -1,size = 3,show.legend = F)+
      geom_path(data = ellipses, aes(group = Group), linewidth = 0.5) +  # 手动添加椭圆
      scale_color_manual(values=n.color) +
      labs(x=sprintf('PC1 (%.2f%%)',summary(data.pca)[6]$importance[2,1]*100),
           y=sprintf('PC2 (%.2f%%)',summary(data.pca)[6]$importance[2,2]*100))
  }else {
    quant.qc3=ggplot(merge(data.pca$x,sample,by.y="Sample",by.x="row.names"),aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_classic(base_size = as.numeric(12))+
      geom_text_repel(aes(label=Row.names))+scale_color_manual(values=n.color)+ 
      labs(x=sprintf('PC1 (%.2f%%)',summary(data.pca)[6]$importance[2,1]*100),
           y=sprintf('PC2 (%.2f%%)',summary(data.pca)[6]$importance[2,2]*100),colour = "Group",shape = "Group")
  }
}else{
  quant.qc3=ggplot(merge(data.pca$x,sample,by.y="Sample",by.x="row.names"),aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_classic(base_size = as.numeric(12))+
    geom_text_repel(aes(label=Row.names))+scale_color_manual(values=n.color)+ 
    labs(x=sprintf('PC1 (%.2f%%)',summary(data.pca)[6]$importance[2,1]*100),
         y=sprintf('PC2 (%.2f%%)',summary(data.pca)[6]$importance[2,2]*100),colour = "Group",shape = "Group")
}
p=pheatmap(cor(data.matrix,method = 'pearson',use = 'pairwise.complete.obs'),border_color = 'white',display_numbers = FALSE ,number_color = 'black',
           color=colorRampPalette(c("steelblue","white","red"))(100),annotation = annot,annotation_colors = annot_color,silent = T,show_rownames = F,
           treeheight_row = 30,treeheight_col = 30)
quant.qc4=ggarrange(p$gtable)+theme(plot.margin = unit(c(0,max(stringr::str_length(annot$Group))*1.5+8,0,0),'mm'))

## Metabolite Classification plots
if(workflow=="Metabolomics"){
  data.classify=meta.data %>% group_by(HMDB_SuperClass,HMDB_Class) %>% summarize(Freq=length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`,collapse = " ; ")) 
  classify.plot1=diff.classy.plot1(data.classify %>% ungroup() %>% slice_max(Freq,n=40, with_ties = FALSE), font=14, color.pal = "Set3", type = workflow, angle=75, margin=18,margin2=12,label.font=4)
  data.classify2=kegg.data %>% merge(meta.data[,c("Metabolite name","HMDB ID")], by = "HMDB ID") %>% ungroup() %>% distinct(`Metabolite name`, `HMDB ID`,`Pathway level 1`,`Pathway level 2`) %>%
    group_by(`Pathway level 1`,`Pathway level 2`) %>% summarize(Freq=length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`, collapse = " ; "), HMDBs = paste(`HMDB ID`, collapse = ";"))
  colnames(data.classify2)=gsub(" ",".", colnames(data.classify2))
  classify.plot2=diff.classy.plot2(data.classify2 %>% slice_max(Freq,n=40), font=14, color.pal = "Set2", angle=60, margin=18, margin2=15)
}
if(workflow=="Lipidomics"){
  data.classify=meta.data %>% group_by(Name, Categories) %>% summarize(Freq=length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`, collapse = " ; "))
  classify.plot1=diff.classy.plot1(data.classify %>% ungroup() %>% slice_max(Freq,n=40, with_ties = FALSE), font=14, color.pal = "Set3", type = workflow, angle=60, margin=18,margin2=12,label.font=4)
}
########################## Export Results #####################################
if(dir.exists(results_dir)) unlink(results_dir,recursive = T,force = T)
dir.create(results_dir)

temp.dir=paste(results_dir, "1-MS_identified_summary",sep="/")
if(dir.exists(temp.dir)) unlink(temp.dir,recursive = T,force = T)
dir.create(temp.dir)

if(workflow =="Metabolomics"){
  fwrite(meta.data,file=paste0(temp.dir,"/Table 1.1 MS identified metabolites detail_raw.xls"),sep = "\t")
  fwrite(data.final.table, file=paste0(temp.dir,"/Table 1.2 MS identified metabolites detail_normalized.xls"),sep = "\t")
}
if(workflow =="Lipidomics"){
  xlsx::write.xlsx(meta.data, file=paste0(temp.dir,"/Table 1.1 MS identified lipids detail_raw.xlsx"),row.names = F, showNA = F)
  xlsx::write.xlsx(data.final.table, file=paste0(temp.dir,"/Table 1.2 MS identified lipids detail_normalized.xlsx"),row.names = F, showNA = F)
}

temp.dir=paste(results_dir,"1-MS_identified_summary/Identification_QC",sep="/")
if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
if(!is.null(qc.plot)){
  ggsave(plot=qc.plot,filename=paste0(temp.dir, "/QC plot.pdf"), width=10, height=12)
  ggsave(plot=qc.plot,filename=paste0(temp.dir, "/QC plot.png"), width=10, height=12)
}
if(workflow=="Metabolomics"){
  file.copy(from = pos.png, to = paste(temp.dir,"/Pos TIC plot.png", sep = "/"),overwrite = T)
  file.copy(from = neg.png, to = paste(temp.dir,"/Neg TIC plot.png", sep = "/"),overwrite = T)
}

temp.dir=paste(results_dir,"1-MS_identified_summary/Quantification_QC",sep="/")
if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
if(workflow=="Metabolomics"){t="Metabolites"}else{t="Lipids"}

ggsave(plot = quant.qc1,filename = paste0(temp.dir,"/Figure 1a. ",t," quantitation density plot.pdf"),width = 8,height = 6)
ggsave(plot = quant.qc1,filename = paste0(temp.dir,"/Figure 1a. ",t," quantitation density plot.png"),width = 8,height = 6)

if(min(table(sample$Group)) != 1){
  ggsave(plot = quant.qc2,filename = paste0(temp.dir,"/Figure 1b. ",t," quantitation CV plot.pdf"),width = 8,height = 6)
  ggsave(plot = quant.qc2,filename = paste0(temp.dir,"/Figure 1b. ",t," quantitation CV plot.png"),width = 8,height = 6)
}

ggsave(plot = quant.qc3,filename = paste0(temp.dir,"/Figure 1c. ",t," quantitation PCA plot.pdf"),width = 8,height = 6)
ggsave(plot = quant.qc3,filename = paste0(temp.dir,"/Figure 1c. ",t," quantitation PCA plot.png"),width = 8,height = 6)
ggsave(plot = quant.qc4,filename = paste0(temp.dir,"/Figure 1d. ",t," quantitation PCC plot.pdf"),width = 8,height = 6,bg = "white")
ggsave(plot = quant.qc4,filename = paste0(temp.dir,"/Figure 1d. ",t," quantitation PCC plot.png"),width = 8,height = 6,bg = "white")

if(workflow=="Metabolomics"){temp.dir=paste(results_dir,"2-Metabolite Functional Classification",sep="/")}
if(workflow=="Lipidomics"){temp.dir=paste(results_dir,"2-Lipid Functional Classification",sep="/")}
if(!dir.exists(temp.dir)){dir.create(path = temp.dir)}
ggsave(plot=classify.plot1,filename=paste0(temp.dir,"/Figure 2a. ",t," Classification by Category.pdf"), width=16, height=8)
ggsave(plot=classify.plot1,filename=paste0(temp.dir,"/Figure 2a. ",t," Classification by Category.png"), width=16, height=8)
if(workflow=="Metabolomics"){
  ggsave(plot=classify.plot2,filename=paste0(temp.dir,"/Figure 2b. Metabolites Classification by Pathway.pdf"), width=12, height=8)
  ggsave(plot=classify.plot2,filename=paste0(temp.dir,"/Figure 2b. Metabolites Classification by Pathway.png"), width=12, height=8)
  fwrite(as.data.frame(data.classify), file=paste0(temp.dir,"/Table 2a. Metabolite Classification by Category.xls"),sep = "\t")
  fwrite(as.data.frame(data.classify2), file=paste0(temp.dir,"/Table 2b. Metabolite Classification by Pathway.xls"),sep = "\t")
}
if(workflow=="Lipidomics"){
  fwrite(as.data.frame(subset(meta.data, select=c("Metabolite name","Name","LION name","LION ID","Categories","ChEBI"))),file=paste0(temp.dir,"/Table 2. Lipid Classification.xls"),sep = "\t")
}

