library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(tidyr)
library(pheatmap)
library(ggpubr)

##################### parameter #########################################
enrich.out.dir=snakemake@params[["input_dir"]]
output.dir=snakemake@params[["out_dir"]]

##################### function defined #########################################
enrich.cluster=function(enrich.data,topn=40,th.pvalue=0.05){
  data=subset(enrich.data,select=c("Compare","Terms","Regulation","p.value"))
  if(nrow(data)<3){return(NULL)}
  data$p.value=-log10(data$p.value+1e-276)
  data[,maxv:=max(p.value),by="Terms"]
  data.wid=data.table::dcast(data.table(data),Terms+maxv~Regulation+Compare,value.var = "p.value",fill = 0,sep=":") %>%
    arrange(-maxv) %>% filter(nchar(as.vector(Terms))<50) %>%
    filter(maxv > -log10(th.pvalue)) %>% head(topn)
  if(ncol(data.wid)<4){return(NULL)}
  if(nrow(data.wid)<3){return(NULL)}
  annot=data.frame(ID=colnames(data.wid)[c(-1,-2)])
  rownames(annot)=annot$ID
  annot=separate(annot,col = "ID",sep = ":",into = c("Regulation","Compare"))
  data.wid=as.data.frame(data.wid)
  rownames(data.wid)=data.wid$Terms
  scalev='row'
  if(ncol(data.wid) == 4){scalev="none"}
  p=pheatmap(data.wid[c(-1,-2)],scale = scalev,annotation = annot,silent = T,cluster_cols = F,border_color = 'white',main = enrich.data$Category[1],show_colnames = F)
  return(p)
}

##################### Execute the main program #########################################
files=dir(enrich.out.dir,pattern = "Table 4. Functional enrichment of ",full.names = T)
file.list=setNames(files,sub(sub(files,pattern=".*Table 4. Functional enrichment of ",replacement=""),pattern=".xls$",replacement=""))
data=rbindlist(lapply(file.list,fread),idcol = "Compare")

plot.list=lapply(split(data,data$Category),function(x) enrich.cluster(x,topn = 40,th.pvalue = 0.05))

for(i in names(plot.list)){
  ggsave(plot.list[[i]]$gtable,filename=file.path(output.dir, paste0("Figure 5. Enrichment cluster of ",i,".pdf")),width = 9,height = 8,create.dir = T)
  ggsave(plot.list[[i]]$gtable,filename=file.path(output.dir, paste0("Figure 5. Enrichment cluster of ",i,".png")),width = 9,height = 8)
}
