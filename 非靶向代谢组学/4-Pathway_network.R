library(data.table)
library(dplyr)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(igraph)
library(ggraph)
library(qgraph)

############################## Parameters ######################################
# workflow = "Metabolomics" # Must be one of "Metabolomics", "Lipidomics" 
# db_dir = "DB/HMDB_20250114"
# results_dir = "Results"

workflow = snakemake@params[["workflow"]]
db_dir = snakemake@params[["db_dir"]]
results_dir = snakemake@params[["results_dir"]]

############################ Functions defined #################################
network.plot = function(data.kegg, gsea.enrich, font=4){
  edge=data.kegg %>% filter(regulation!='unchanged',Terms %in% gsea.enrich$Terms[gsea.enrich$p.value<0.05]) %>% dplyr::select(Terms,`Metabolite name`) %>% distinct_all()
  node.p=gsea.enrich %>% filter(Terms %in% edge$Terms) %>% dplyr::select(Terms,Count,NES) %>% mutate(shape="pathway")
  colnames(node.p)=c("id","size","fill","shape")
  node.m=data.kegg %>% filter(`Metabolite name` %in% edge$`Metabolite name`) %>% dplyr::select(`Metabolite name`,VIP,log2fc) %>% mutate(shape="metabolite") %>% distinct_all() 
  colnames(node.m)=c("id","size","fill","shape")
  node=rbind(node.m,node.p) %>% mutate(color=case_when(fill>0~"up",TRUE~'down')) %>% distinct(id,.keep_all = T)
  
  ig=graph_from_data_frame(edge,directed = F,vertices = node)
  e=get.edgelist(ig,names = FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(ig),area=8*(vcount(ig)^2),repulse.rad=(vcount(ig)^3.1))
  ggraph(ig,layout = l)+geom_edge_link(color='grey')+
    geom_node_point(aes(size=size,fill=fill,color=color,shape=shape))+
    geom_node_text(aes(label=stringr::str_wrap(name, width=40)),size=as.numeric(font),repel = T)+
    scale_shape_manual(values = c(metabolite=21,pathway=23))+
    scale_size(range=c(2,8))+theme_void()+
    scale_fill_gradientn(colours = c(pal_npg()(3)[3],"white",pal_npg()(3)[1]),limits=c(-1,1),oob=scales::squish)+
    scale_color_manual(values = c(up=pal_npg()(3)[1],down=pal_npg()(3)[3]))  
}
l.network.plot=function(network,data,lion.data, font=4){
  node=data.frame(ChEBI=unique(c(network$source,network$target))) %>% merge(subset(lion.data,select=c("LION name","ChEBI")) %>% separate_rows(ChEBI,sep="\\|") %>% distinct(),by="ChEBI") %>%
    merge(data,by="LION name") %>% distinct(ChEBI,.keep_all = T) %>% dplyr::select(ChEBI,`Metabolite name`,log2fc,VIP,regulation) %>% distinct()
  network=network %>% filter(source %in% node$ChEBI & target %in% node$ChEBI)
  
  ig=graph_from_data_frame(network,directed = F,vertices = node)
  e=get.edgelist(ig,names = FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(ig),area=8*(vcount(ig)^2),repulse.rad=(vcount(ig)^3.1))
  ggraph(ig,layout = l)+geom_edge_link(aes(color=reaction_type))+
    geom_node_point(shape=21,aes(fill=regulation,size=VIP))+
    geom_node_text(aes(label=`Metabolite name`),size=as.numeric(font))+
    scale_fill_manual(values = c(up=pal_npg()(3)[1],unchanged='grey',down=pal_npg()(3)[3]))+
    theme_void() + theme(plot.margin=unit(c(5,5,5,5),'mm'))
}

############################ Network plot #####################################
values = list()
if (workflow =="Metabolomics"){
  files = list.files(path = paste(results_dir,"4-Functional Annotation and Enrichment",sep="/"),pattern = "Table 4a")
  data.kegg = lapply(files, function(x) fread(paste(results_dir,"4-Functional Annotation and Enrichment",x,sep="/")))
  names(data.kegg) = gsub("Table 4a.KEGG Annotation of ","",gsub(".xls","", files))
  
  files = list.files(path = paste(results_dir,"4-Functional Annotation and Enrichment",sep="/"),pattern = "Table 4c")
  gsea.enrich = lapply(files, function(x) fread(paste(results_dir,"4-Functional Annotation and Enrichment",x,sep="/")))
  names(gsea.enrich) = gsub("Table 4c.MSEA Enrichment of ","", gsub(".xls","",files))
  com = names(data.kegg)
  for (x in names(data.kegg)){
    if(!is.null(gsea.enrich[[x]])){
      values[[paste(x,"network.plot",sep=".")]]=network.plot(data.kegg[[x]],gsea.enrich[[x]], font=4)
    }
  }
}
if(workflow =="Lipidomics"){
  temp.dir=paste(results_dir,"3-Differentially abundant lipids",sep="/")
  files = list.files(path = temp.dir, pattern = "Table 3a.")
  diff.list = lapply(files, function(x) fread(paste(temp.dir,x,sep="/")))
  names(diff.list) = gsub("Table 3a. Differentially abundant lipids list of ","",gsub(".xls","", files))
  
  lion.anno=fread(paste0(db_dir, "/LION_data.txt"))
  network=fread(paste0(db_dir,"/Network.txt"),sep="\t")
  for (x in names(diff.list)){
    values[[paste(x,"network.plot",sep=".")]]=l.network.plot(network, diff.list[[x]], lion.anno %>% separate_rows(ChEBI,sep =" \\| "),font=4)
  }
  com = names(diff.list)
}

############################ Export plots #####################################
temp.dir=paste(results_dir,"5-Pathway Network Analysis",sep="/")
if(!dir.exists(temp.dir)) unlink(temp.dir,recursive = T,force = T)
dir.create(path = temp.dir)

for (x in com){
  if(!is.null(values[[paste(x,"network.plot",sep=".")]])){
    if(workflow=="Metabolomics"){
      ggsave(values[[paste(x,"network.plot",sep=".")]], filename = paste0(temp.dir,"/Figure 5. Pathway network plot of ", x,".pdf"),width=12, height=12)
      ggsave(values[[paste(x,"network.plot",sep=".")]], filename = paste0(temp.dir,"/Figure 5. Pathway network plot of ", x,".png"),width=12, height=12, bg="white")
    }
    if(workflow=="Lipidomics"){
      ggsave(values[[paste(x,"network.plot",sep=".")]], filename = paste0(temp.dir,"/Figure 5. Lipids reactome network plot of ", x,".pdf"),width=14, height=14)
      ggsave(values[[paste(x,"network.plot",sep=".")]], filename = paste0(temp.dir,"/Figure 5. Lipids reactome network plot of ", x,".png"),width=14, height=14, bg="white")
    }
  }
}

