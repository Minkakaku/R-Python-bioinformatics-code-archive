library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)

##################### parameter #########################################
#db.dir="C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/human"
#out.dir="C:/Users/Administrator/Desktop/WORK/武靖/富集圈图+和弦图/demo"
#diff.data.file="Table 2.DEPs list of HD-vs-CON.xls"

# 提取对比组名称
#compare.name = gsub(".*Table 2.DEPs list of ","",diff.data.file)
#compare.name = gsub(".xls$","",compare.name)

# 创建文件夹
#dbs <- c("GO","KEGG","Reactome","NOG","InterPro","Wikipathways")
#invisible(lapply(file.path(out.dir, dbs), dir.create, showWarnings = FALSE))


##################### enrich test #####################
enrich.test=function(x){
  a=x[1]; b=x[2]-a; c=x[3]-a; d=x[4]-a-b-c
  fisher.test(matrix(c(a,b,c,d),2), alternative="greater")$p.value
}


##################### enrich.single #####################
enrich.single=function(db.data, diff.data, DBtype="Other"){
  
  db.diff = merge(db.data, diff.data[,c("Protein accession","Gene name","log2FC")],
                  by="Protein accession")
  
  term.info = db.diff %>%
    group_by(Terms) %>%
    summarise(
      `Gene names` = paste(unique(`Gene name`), collapse=";"),
      Proteins = paste(unique(`Protein accession`), collapse=";"),
      up_count = sum(log2FC>0),
      down_count = sum(log2FC<0),
      hit_n = n_distinct(`Protein accession`)
    )
  
  mapping = data.frame(table(db.diff$Terms))
  background = data.frame(table(db.data$Terms))
  
  state = merge(mapping, background, by="Var1")
  colnames(state)=c("Terms","mapping","background")
  
  state$`All mapping`=length(unique(db.diff$`Protein accession`))
  state$`All background`=length(unique(db.data$`Protein accession`))
  
  state$`Fold enrich` = round(
    (state$mapping/state$background) /
      (state$`All mapping`/state$`All background`),2)
  
  state = subset(state, mapping>1 & `Fold enrich`>1)
  
  state$p.value = signif(apply(state[,c("mapping","background","All mapping","All background")],1,enrich.test),3)
  state$fdr = signif(p.adjust(state$p.value,"BH"),3)
  
  final = merge(state, term.info, by="Terms")
  
  if(DBtype=="GO"){
    ontology = db.data %>% distinct(Terms,Category) %>%
      mutate(ONTOLOGY = recode(Category,
                               "Biological process"="BP",
                               "Molecular function"="MF",
                               "Cellular component"="CC"))
    final = merge(final, ontology[,c("Terms","ONTOLOGY")], by="Terms", all.x=TRUE)
  }
  
  return(final)
}


##################### plot functions #####################
plot_GO_facet=function(go_df, top_n=10, wrap_width=40){
  go_df$logFDR = -log10(go_df$fdr)
  topGO = go_df %>% group_by(ONTOLOGY) %>% arrange(fdr) %>% slice(1:top_n)
  topGO$Terms = str_wrap(topGO$Terms, wrap_width)
  topGO$ONTOLOGY=factor(topGO$ONTOLOGY,levels=c("BP","MF","CC"))
  
  ggplot(topGO,
         aes(x=`Fold enrich`, y=reorder(Terms,`Fold enrich`)))+
    geom_point(aes(size=mapping,color=logFDR))+
    scale_color_gradient(low="#C7E4E2",high="#F08686")+
    scale_size(range=c(3,10))+
    facet_wrap(~ONTOLOGY, ncol=1, scales="free_y", strip.position="right")+
    theme_bw(14)
}

plot_DB_bubble=function(df, title_name="Enrichment", top_n=20, wrap_width=50){
  df$logFDR=-log10(df$fdr)
  df_top=df %>% arrange(fdr) %>% slice(1:top_n)
  df_top$Terms=str_wrap(df_top$Terms, wrap_width)
  
  ggplot(df_top,
         aes(x=`Fold enrich`, y=reorder(Terms,`Fold enrich`)))+
    geom_point(aes(size=mapping,color=logFDR))+
    scale_color_gradient(low="#C7E4E2",high="#F08686")+
    scale_size(range=c(3,10))+
    labs(title=title_name)+
    theme_bw(14)
}

plot_DB_bar=function(df,title_name="Enrichment",top_n=20,wrap_width=50){
  df$logFDR=-log10(df$fdr)
  df_top=df %>% arrange(fdr) %>% slice(1:top_n)
  df_top$Terms=str_wrap(df_top$Terms,wrap_width)
  
  ggplot(df_top,
         aes(x=mapping, y=reorder(Terms,mapping), fill=logFDR))+
    geom_col(width=0.7)+
    scale_fill_gradient(low="#C7E4E2",high="#F08686")+
    labs(title=title_name, x="Gene Count (mapping)")+
    theme_bw(14)
}


##################### load diff and annotation #####################
data.diff = fread(diff.data.file) %>%
  filter(Regulation!="unchanged") %>%
  distinct(`Protein accession`, .keep_all=TRUE)


##################### load all background DBs #####################
go.data       = fread(file.path(db.dir,"GO_annotation.tsv"))
kegg.data     = fread(file.path(db.dir,"KEGG_annotation.tsv"))
reactome.data = fread(file.path(db.dir,"Reactome_annotation.tsv"))
nog.data      = fread(file.path(db.dir,"NOG_annotation.tsv"))

ipr.data      = fread(file.path(db.dir,"InterPro_annotation.tsv"))
setnames(ipr.data,"Protein_accession","Protein accession")
ipr.data$Category="InterPro"

wiki.data     = fread(file.path(db.dir,"Wikipathways_annotation.tsv"))
setnames(wiki.data,"Protein_accession","Protein accession")
wiki.data$Category="Wikipathways"


##################### run enrich #####################
res.go       = enrich.single(go.data, data.diff,"GO")
res.kegg     = enrich.single(kegg.data, data.diff)
res.reactome = enrich.single(reactome.data, data.diff)
res.nog      = enrich.single(nog.data, data.diff)
res.ipr      = enrich.single(ipr.data, data.diff)
res.wiki     = enrich.single(wiki.data, data.diff)


##################### save results + plots #####################
save_enrich=function(res, name){
  fwrite(res, file.path(out.dir,name,paste0(name,"_enrichment_",compare.name,".xls")), sep="\t")
  p1=plot_DB_bubble(res, paste0(name," Enrichment"))
  ggsave(file.path(out.dir,name,paste0(name,"_bubble_",compare.name,".pdf")), p1, width=8,height=6)
  p2=plot_DB_bar(res, paste0(name," Enrichment"))
  ggsave(file.path(out.dir,name,paste0(name,"_bar_",compare.name,".pdf")), p2, width=8,height=6)
}

save_enrich(res.kegg,"KEGG")
save_enrich(res.reactome,"Reactome")
save_enrich(res.nog,"NOG")
save_enrich(res.ipr,"InterPro")
save_enrich(res.wiki,"Wikipathways")

# GO 特殊输出
fwrite(res.go, file.path(out.dir,"GO",paste0("GO_enrichment_",compare.name,".xls")), sep="\t")
ggsave(file.path(out.dir,"GO",paste0("GO_bubble_",compare.name,".pdf")), plot_GO_facet(res.go), width=10,height=12)
ggsave(file.path(out.dir,"GO",paste0("GO_bar_",compare.name,".pdf")), plot_GO_bar(res.go), width=10,height=12)
