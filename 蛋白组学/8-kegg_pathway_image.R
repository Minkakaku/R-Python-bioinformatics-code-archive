# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(XML)

##################### parameter #########################################
db.dir=snakemake@params[["db_dir"]]
diff.data.file=snakemake@input[["dep_results"]]
enrich.file=snakemake@input[["enrich_file"]]
out.dir=snakemake@params[["out_dir"]]

##################### function defined #########################################
pathway_color_plot=function(kegg.dir,kegg.data,out.dir){
  tag=kegg.data$`Pathway ID`[1]
  oufile=paste0(out.dir,"/image/",tag,".png")
  xmlfile=paste0(kegg.dir,tag,".xml")
  pngfile=paste0(kegg.dir,tag,".png")
  htmlfile=paste0(kegg.dir,tag,".html")
  out=rbindlist(lapply(getNodeSet(xmlParse(xmlfile),path = '//entry'),function(x){
    id=xmlAttrs(x)['name']
    if(xmlAttrs(x)['type']=='gene'|xmlAttrs(x)['type']=='ortholog'){
      rbindlist(lapply(getNodeSet(x,path = 'graphics'),function(y){
        if(xmlAttrs(y)['type']=="line"){
          as.data.frame(t(xmlAttrs(y))) %>% separate(col = "coords",into = c("x","y","width","height"),sep = ',') %>%
            mutate(x=as.numeric(x),width=as.numeric(width),y=as.numeric(y)) %>%
            mutate(width=width-x,height=16) %>% mutate(x=x+width/2)
        }else{
          as.data.frame(t(xmlAttrs(y)))
        }
      }),fill = T) %>% mutate(`KEGG ID`=id)
    }
  }),fill = T) %>% separate_rows(`KEGG ID`,sep=" ") %>% distinct(`KEGG ID`,x,y,width,height,type) %>% merge(kegg.data,by="KEGG ID") %>%
    mutate(x=as.numeric(x),y= -as.numeric(y),width=as.numeric(width),height=as.numeric(height))
  img=png::readPNG(pngfile)
  w=ncol(img);h=nrow(img)
  p=ggplot(data=out %>% distinct(x,y,.keep_all = T),aes(x=x+0.5,y = y))+
    annotation_custom(grid::rasterGrob(img,interpolate = T))+
    theme_void(base_family = 'DejaVu Serif')+geom_tile(aes(width=width,height=height,fill=log2FC))+
    geom_text(aes(label=`Gene name`),size=1.5)+
    scale_fill_gradientn(colours = c("blue","white","red"),limits=c(-1,1),oob=scales::squish)+
    scale_x_continuous(expand = expansion(mult = 0),limits = c(0,w))+
    scale_y_continuous(expand = expansion(mult = 0),limits = c(-h,0))+
    theme(legend.position=c(0.8,0.95),legend.direction = 'horizontal',
          legend.title = element_blank(),legend.text = element_text(size = 4),
          legend.key.height = unit(3,'mm'),legend.key.width = unit(4,'mm'))
  ggsave(p,filename=oufile,width = w,height = h,units = 'px',dpi = 300)
  file.copy(htmlfile,out.dir)
}


##################### Execute the main program #########################################

## load data
compare.name=sub(sub(diff.data.file,pattern=".*Table 2. DEPs list of ",replacement=""),pattern=".xls$",replacement="")
data.kegg.enrich=fread(enrich.file) %>% filter(p.value<0.05)
data.kegg.diff=merge(fread(file=paste(db.dir,"KEGG_annotation.tsv",sep="/")),fread(diff.data.file) %>%
                       filter(Regulation!='unchanged'),by="Protein accession") %>% filter(Terms %in% data.kegg.enrich$Terms) %>%
  select(`Pathway ID`,`KEGG ID`,log2FC,`Gene name`) %>% distinct(`Pathway ID`,`KEGG ID`,`Gene name`,.keep_all = T)

temp.dirx=paste(out.dir,compare.name,sep="/")
if(dir.exists(temp.dirx)){unlink(temp.dirx,force = T,recursive = T)}
dir.create(path = temp.dirx)
dir.create(path = paste0(temp.dirx,"/image"))

## plot kegg pathway
kegg.dir=paste(db.dir,"kegg.dir/",sep="/")
lapply(split(data.kegg.diff,data.kegg.diff$`Pathway ID`),function(x){pathway_color_plot(kegg.dir,x,temp.dirx)})

# message("KEGG pathway image分析完成！")
# sink(type = "message")
# close(log_con)