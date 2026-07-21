# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)

sub.go.list=c("GO:0009507","GO:0005737","GO:0005856","GO:0005783","GO:0005576","GO:0005794","GO:0005739","GO:0016604","GO:0005777","GO:0005886","GO:0005774")

classif.go=function(go.data,diff.list){
    go.data=merge(go.data,diff.list,by="Protein accession")
    go.data=subset(go.data,`GO levels`==2,select=c("Category","Terms","Regulation"))
    data=subset(as.data.frame(table(go.data)),Freq>0&nchar(as.character(Terms))<60)
    data=data[order(data$Category,data$Regulation,data$Freq),]
    data$Terms=factor(data$Terms,levels = unique(data$Terms))
    data
}

classif.nog=function(nog.data,diff.list){
    nog.data=merge(nog.data,diff.list,by="Protein accession")
    nog.data=subset(nog.data,select=c("Terms","Regulation"))
    data=subset(as.data.frame(table(nog.data)),Freq>0)
    data=data[order(data$Regulation,data$Freq),]
    data$Terms=factor(data$Terms,levels = unique(data$Terms))
    data
}

classif.subcell=function(go.data,diff.list){
    go.data=merge(go.data,diff.list,by="Protein accession")
    go.data=subset(go.data,ID %in% sub.go.list,select=c("Terms","Regulation"))
    data=subset(as.data.frame(table(go.data)),Freq>0)
    data=data[order(data$Regulation,data$Freq),]
    data$Terms=factor(data$Terms,levels = unique(data$Terms))
    return(data)
}

classif.kegg=function(kegg.data,diff.list){
    kegg.data=merge(kegg.data,diff.list,by="Protein accession")
    kegg.data=subset(kegg.data,select=c("Pathway level 1","Pathway level 2","Regulation"))
    data=subset(data.table(table(kegg.data)),N>0) %>% dplyr::rename(Freq=N)
    data=data[order(data$`Pathway level 1`,data$Regulation,data$Freq),]
    data$`Pathway level 2`=factor(data$`Pathway level 2`,levels = unique(data$`Pathway level 2`))
    return(data)
}

message("读取蛋白注释文件...")
go.data=fread(paste(snakemake@params[["db_dir"]],'GO_annotation.tsv',sep="/"))
kegg.data=fread(paste(snakemake@params[["db_dir"]],'KEGG_annotation.tsv',sep="/"))
nog.data=fread(paste(snakemake@params[["db_dir"]],'NOG_annotation.tsv',sep="/"))

message("读取差异表达蛋白数据...")
diff.data <- fread(snakemake@input[["dep_results"]])
diff.list=subset(diff.data, Regulation!='unchanged') %>% distinct(`Protein accession`,Regulation,.keep_all = T)

goclass.data=classif.go(go.data,diff.list)
nogclass.data=classif.nog(nog.data,diff.list)
subcellclass.data=classif.subcell(go.data,diff.list)
keggclass.data=classif.kegg(kegg.data,diff.list)

fwrite(goclass.data,file=snakemake@output[["go_classif"]],sep = "\t")
fwrite(nogclass.data,file=snakemake@output[["kog_classif"]],sep = "\t")
fwrite(subcellclass.data,file=snakemake@output[["subcell_classif"]],sep = "\t")
fwrite(keggclass.data,file=snakemake@output[["kegg_classif"]],sep = "\t")

message("绘图...")
goclass.plot=ggplot(goclass.data,aes(x=Terms,y=Freq,fill=Category))+
                  geom_bar(stat="identity")+theme_test(base_size = 14)+
                  facet_grid(.~Regulation,scales = 'free_y',space = 'free')+coord_flip()+theme(legend.position = "top",legend.justification = c(1,1))+
                  geom_text(aes(label=Freq),hjust=0)+scale_y_continuous(expand = expansion(mult = c(0,0.2)))
ggsave(goclass.plot,filename=snakemake@output[["go_classif_plot"]],width = 8,height = 8)
ggsave(goclass.plot,filename=sub(snakemake@output[["go_classif_plot"]],pattern=".png",replacement=".pdf"),width = 8,height = 8)

nogclass.plot=ggplot(nogclass.data,aes(x=Terms,y=Freq,fill=Terms))+
                  geom_bar(stat="identity")+theme_test(base_size = 14)+facet_grid(.~Regulation,scales = 'free_y',space = 'free')+
                  coord_flip()+theme(legend.position = "none")+geom_text(aes(label=Freq),hjust=0)+
                  scale_y_continuous(expand = expansion(mult = c(0,0.2)))
ggsave(nogclass.plot,filename=snakemake@output[["kog_classif_plot"]],width = 8,height = 8)
ggsave(nogclass.plot,filename=sub(snakemake@output[["kog_classif_plot"]],pattern=".png",replacement=".pdf"),width = 8,height = 8)

subcellclass.plot=ggplot(subcellclass.data,aes(x=Terms,y=Freq,fill=Terms))+
                  geom_bar(stat="identity")+theme_test(base_size = 14)+
                  facet_grid(.~Regulation,scales = 'free_y',space = 'free')+coord_flip()+theme(legend.position = "none")+
                  geom_text(aes(label=Freq),hjust=0)+scale_y_continuous(expand = expansion(mult = c(0,0.2)))
ggsave(subcellclass.plot,filename=snakemake@output[["subcell_classif_plot"]],width = 8,height = 8)
ggsave(subcellclass.plot,filename=sub(snakemake@output[["subcell_classif_plot"]],pattern=".png",replacement=".pdf"),width = 8,height = 8)

keggclass.plot=ggplot(keggclass.data,aes(x=`Pathway level 2`,y=Freq,fill=`Pathway level 1`))+
                  geom_bar(stat="identity")+theme_test(base_size = 14)+
                  facet_grid(.~Regulation,scales = 'free_y',space = 'free')+
                  coord_flip()+theme(legend.position = "top",legend.justification = c(1,1))+
                  guides(fill=guide_legend(nrow = 3))+geom_text(aes(label=Freq),hjust=0)+
                  scale_y_continuous(expand = expansion(mult = c(0,0.2)))
ggsave(keggclass.plot,filename=snakemake@output[["kegg_classif_plot"]],width = 8,height = 8)
ggsave(keggclass.plot,filename=sub(snakemake@output[["kegg_classif_plot"]],pattern=".png",replacement=".pdf"),width = 8,height = 8)


# message("功能分类分析完成！")
# sink(type = "message")
# close(log_con)