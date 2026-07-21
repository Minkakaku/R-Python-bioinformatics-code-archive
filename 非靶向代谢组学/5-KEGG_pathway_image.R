library(data.table)
library(dplyr)
library(XML)
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyverse)

############################## Parameters ######################################
workflow = snakemake@params[["workflow"]]
kegg_dir = snakemake@params[["kegg_dir"]]
results_dir = snakemake@params[["results_dir"]]

############################ Functions defined #################################
pathway_color_plot=function(kegg.dir,kegg.data,out.dir){
  tag=kegg.data$`Pathway ID`[1]
  oufile=paste0(out.dir,"/image/",tag,".png")
  xmlfile=paste0(kegg.dir,tag,".xml")
  pngfile=paste0(kegg.dir,tag,".png")
  htmlfile=paste0(kegg.dir,tag,".html")
  out=rbindlist(lapply(getNodeSet(xmlParse(xmlfile),path = '//entry'),function(x){
    id=xmlAttrs(x)['name']
    if(xmlAttrs(x)['type']=='compound'){
      rbindlist(lapply(getNodeSet(x,path = 'graphics'),function(y){
          as.data.frame(t(xmlAttrs(y)))
      }),fill = T) %>%  mutate(cpd_id=sub("cpd:","",id))
    }
  }),fill = T) %>% separate_rows(cpd_id,sep=" ") %>% distinct(cpd_id,x,y,width,height,type) %>% merge(kegg.data,by="cpd_id") %>%
    mutate(x=as.numeric(x),y= -as.numeric(y),width=as.numeric(width),height=as.numeric(height))
  img=png::readPNG(pngfile)
  w=ncol(img);h=nrow(img)
  p=ggplot(data=out %>% distinct(x,y,.keep_all = T),aes(x=x+0.5,y = y))+
    annotation_custom(grid::rasterGrob(img,interpolate = T))+
    theme_void(base_family = 'DejaVu Serif')+geom_point(aes(fill=log2fc),shape = 21,stroke = 0.3)+
    # geom_text(aes(label=`Gene name`),size=1.5)+
    scale_fill_gradientn(colours = c("blue","white","red"),limits=c(-1,1),oob=scales::squish)+
    scale_x_continuous(expand = expansion(mult = 0),limits = c(0,w))+
    scale_y_continuous(expand = expansion(mult = 0),limits = c(-h,0))+
    theme(legend.position=c(0.8,0.95),legend.direction = 'horizontal',
          legend.title = element_blank(),legend.text = element_text(size = 4),
          legend.key.height = unit(3,'mm'),legend.key.width = unit(4,'mm'))
  ggsave(p,filename=oufile,width = w,height = h,units = 'px',dpi = 300)
  file.copy(htmlfile,out.dir)
}

############################ KEGG pathway images ###############################
if(workflow =="Metabolomics"){
  files = list.files(path = paste(results_dir,"3-Differentially abundant metabolites",sep="/"), pattern = "Table 3a")
  diff.list = lapply(files, function(x) fread(paste(results_dir,"3-Differentially abundant metabolites",x,sep="/")))
  names(diff.list) = gsub("Table 3a. Differentially abundant metabolites list of ","",gsub(".xls","", files))
  
  files = list.files(path = paste(results_dir,"4-Functional Annotation and Enrichment",sep="/"),pattern = "Table 4a")
  kegg.data = lapply(files, function(x) fread(paste(results_dir,"4-Functional Annotation and Enrichment",x,sep="/")))
  names(kegg.data) = gsub("Table 4a.KEGG Annotation of ","",gsub(".xls","", files))
  
  files = list.files(path = paste(results_dir,"4-Functional Annotation and Enrichment",sep="/"),pattern = "Table 4b")
  kegg.enrich = lapply(files, function(x) fread(paste(results_dir,"4-Functional Annotation and Enrichment",x,sep="/")))
  names(kegg.enrich) = gsub("Table 4b.KEGG Enrichment of ","", gsub(".xls","",files))
  
  files = list.files(path = paste(results_dir,"4-Functional Annotation and Enrichment",sep="/"),pattern = "Table 4c")
  gsea.enrich = lapply(files, function(x) fread(paste(results_dir,"4-Functional Annotation and Enrichment",x,sep="/")))
  names(gsea.enrich) = gsub("Table 4c.MSEA Enrichment of ","", gsub(".xls","",files))
  
  temp.dir=paste(results_dir,"6-KEGG pathway image",sep ="/")
  if(dir.exists(temp.dir)) unlink(temp.dir,recursive = T,force = T)
  dir.create(path=temp.dir) 
  total =length(diff.list)
  for (x in names(diff.list)){
    temp.dir2=paste(temp.dir,x,sep = "/")
    if(!dir.exists(temp.dir2)){dir.create(path=temp.dir2)}
    temp.dirx = paste(temp.dir2,"image",sep = "/")
    if(!dir.exists(temp.dirx)){dir.create(path=temp.dirx)}
    
    kegg.res=kegg.enrich[[x]] %>% filter(p.value < 0.05)
    gsea.res=gsea.enrich[[x]] %>% filter(p.value < 0.05)
    terms = unique(c(as.character(kegg.res$Terms),gsea.res$Terms))
    if(length(terms)>0){
      data.kegg = kegg.data[[x]] %>% filter(Terms %in% terms)
      data.kegg = subset(data.kegg, select=c("Pathway ID","cpd_id","Metabolite name","log2fc")) %>% distinct()
      lapply(split(data.kegg,data.kegg$`Pathway ID`),function(x){
        pathway_color_plot(kegg_dir,x,temp.dir2)
      })
    }
  }
}

