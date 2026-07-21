# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)

##################### parameter #########################################
db.dir="C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/hum"
out.dir="C:/Users/Administrator/Desktop/WORK/20260128公共数据蛋白组"
diff.data.file="C:/Users/Administrator/Desktop/WORK/20260128公共数据蛋白组/DEG/Table 2. DEPs list of FAM134B_EBSS-vs-EBSS.xls"
##################### function defined #########################################
enrich.test=function(x){
  a=x[1]
  b=x[2]-a
  c=x[3]-a
  d=x[4]-a-b-c
  p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value
}
enrich.demo=function(go.data,diff.list){
  
  go.data.diff=merge(go.data,subset(diff.list,select=c("Protein accession","Gene name")),by="Protein accession")
  temp.s=go.data.diff %>% group_by(Terms) %>% mutate(`Gene names`=paste(`Gene name`,collapse = ";"),`Proteins`=paste(`Protein accession`,collapse = ';')) %>%
    ungroup() %>% distinct(Terms,`Gene names`,Proteins)
  
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","backgroupd")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background`=length(unique(go.data$`Protein accession`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  go.data.state=subset(go.data.state,mapping>1&`Fold enrich`>1)
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value),3)
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich","p.value","fdr"))
  go.data.state=merge(go.data.state,temp.s,by="Terms")
  return(go.data.state)
}
enrich.combine=function(go.data,diff.list){
  a=split.data.frame(diff.list,diff.list$Regulation)
  enrich.out=lapply(a,function(x) enrich.demo(go.data,x))
  enrich.out.data=rbindlist(enrich.out,idcol = "Regulation")
}

##################### Execute the main program #########################################
compare.name=sub(sub(diff.data.file,pattern=".*Table 2. DEPs list of ",replacement=""),pattern=".xls$",replacement="")
data.diff=fread(diff.data.file) %>% filter(Regulation!='unchanged') %>% distinct(`Protein accession`,Regulation,.keep_all = T)
go.data=fread(file=paste(db.dir,"GO_annotation.tsv",sep="/"))
nog.data=fread(file=paste(db.dir,"NOG_annotation.tsv",sep="/"))
pfam.data=fread(file=paste(db.dir,"Pfam_annotation.tsv",sep="/"))
kegg.data=fread(file=paste(db.dir,"KEGG_annotation.tsv",sep="/"))
annot.data=fread(file=paste(db.dir,"annotation.tsv",sep="/"))

## annotation data combine
pfam.data$Category="Pfam";kegg.data$Category="KEGG"
gmt.data=rbindlist(list(
  subset(go.data,select=c("Protein accession","Category","Terms")),
  subset(pfam.data,select=c("Protein accession","Category","Terms")),
  subset(kegg.data,select=c("Protein accession","Category","Terms"))
)) %>% distinct_all()

## enrichment analysis
enrich.out=lapply(split(gmt.data,gmt.data$Category),function(x) enrich.combine(x,data.diff))

## data output
fwrite(rbindlist(enrich.out,idcol = "Category"),file = file.path(out.dir, paste0('Table 4. Functional enrichment of ',compare.name,'.xls')),sep="\t")

## plot output
lapply(names(enrich.out),function(x){
  p=ggplot(enrich.out[[x]][p.value<0.05&nchar(as.character(Terms))<60,.SD[frankv(p.value, ties.method = "dense") < 20],by="Regulation"],aes(y=reorder(Terms,`Fold enrich`),x=log2(`Fold enrich`),color=p.value))+
    geom_point(aes(size=mapping))+theme_test(base_size = 14)+ylab("Terms")+facet_grid(.~Regulation)+scale_color_material(reverse = T,palette = 'red')+
    scale_size(range = c(3,6))+scale_x_continuous(expand = expansion(mult = c(0.2,0.2)))+scale_y_discrete(expand = expansion(mult = c(0.02,0.02)))
  ggsave(plot = p,filename = file.path(out.dir, paste0("Figure 4. ",x," enrichment of ",compare.name,".pdf")),width = 10,height = 8)
  ggsave(plot = p,filename = file.path(out.dir, paste0("Figure 4. ",x," enrichment of ",compare.name,".png")),width = 10,height = 8)
})


# message("差异蛋白功能富集分析完成！")
# sink(type = "message")
# close(log_con)
