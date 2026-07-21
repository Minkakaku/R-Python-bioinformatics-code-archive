# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(data.table)
library(dplyr)
library(tidyr)
library(RCurl)
library(visNetwork)
library(ggsci)

##################### parameter #########################################
db.dir="C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/mum"
out.dir="C:/Users/Administrator/Desktop/WORK/20260128公共数据蛋白组"
diff.data.file="C:/Users/Administrator/Desktop/WORK/20260128公共数据蛋白组/DEG/Table 2. DEPs list of EBSS-vs-Dox.xls"

##################### function defined #########################################
get.ppi=function(annot.data,diff.list){
  diff.list=diff.list %>% mutate(score= -log10(adj.P.Val)*abs(log2FC)) %>% slice_max(score,n=500)
  data=merge(subset(diff.list,select = c("Protein accession","Regulation")),
             subset(annot.data,!is.na(`Cross-reference (STRING)`),select=c("Entry","Cross-reference (STRING)")),
             by.x="Protein accession",by.y="Entry")
  
  data=separate_rows(data,`Cross-reference (STRING)`,sep = ";")
  data=subset(data,select=c("Cross-reference (STRING)","Regulation")) %>% distinct(`Cross-reference (STRING)`,.keep_all = T)
  colnames(data)=c("Label","Regulation")
  data$ID=c(1:nrow(data))
  string_api_url = "https://cn.string-db.org/api/"
  output_format = "tsv"
  method = "network"
  uri=paste(string_api_url,output_format,method,sep="/")
  species=names(sort(table(sub(data$Label,pattern = "\\..*",replacement = "")),decreasing = T)[1])
  params = list(
    "identifiers"=paste(collapse ="%0d",data$Label),
    "species"=species,
    "required_score"=700
  )
  response=postForm(uri = uri,.params = params)
  out.data=fread(response)
  out.data=merge(out.data,data,by.x="stringId_A",by.y="Label")
  out.data=merge(out.data,data,by.x="stringId_B",by.y="Label")
  edge.data=subset(out.data,select=c("ID.x","ID.y"))
  colnames(edge.data)=c("from","to")
  node.data1=subset(out.data,select=c("ID.x","preferredName_A","Regulation.x"))
  node.data2=subset(out.data,select=c("ID.y","preferredName_B","Regulation.y"))
  colnames(node.data1)=c("id","label","group")
  colnames(node.data2)=c("id","label","group")
  node.data=data.table(rbind(node.data1,node.data2))
  node.data[,value:=length(label),by=id]
  node.data=unique(node.data)
  return(list(edge=edge.data,node=node.data))
}
ppi.top=function(ppi.data,topn){
  node=ppi.data$node
  edge=ppi.data$edge
  node=node[order(-value)]
  topid=head(node,topn)$id
  edge2=subset(edge, from %in% topid | to %in% topid)
  allid=unique(c(edge2$from,edge2$to))
  node=subset(node,id %in% allid)
  edge=subset(edge, from %in% allid & to %in% allid)
  return(list(node=node,edge=edge))
}

##################### Execute the main program #########################################

## load data
compare.name=sub(sub(diff.data.file,pattern=".*Table 2. DEPs list of ",replacement=""),pattern=".xls$",replacement="")
data.diff=fread(diff.data.file) %>% filter(Regulation!='unchanged')
annot.data=fread(paste(db.dir,"annotation.tsv",sep="/"))

ppi.data=ppi.top(get.ppi(annot.data,data.diff),50)

ppi.plot=visNetwork(ppi.data$node,ppi.data$edge) %>% visInteraction(hover = TRUE) %>%
  visIgraphLayout(layout = 'layout_with_kk') %>%
  visExport(background = "white",float = "left", label = "Save network",type = 'png') %>%
  visNodes(font=list(size=30,color="black")) %>%
  visGroups(groupname = 'up-regulated',color="#D95F02") %>%
  visGroups(groupname = 'down-regulated',color="#1B9E77") %>%
  visOptions(highlightNearest = TRUE,width = "100%",height = "800px")

visSave(ppi.plot, file = file.path(out.dir, paste0("Figure 6. PPI of ",compare.name,".html")),selfcontained = T,background = 'white')


# message("差异蛋白互作网络分析完成！")
# sink(type = "message")
# close(log_con)
