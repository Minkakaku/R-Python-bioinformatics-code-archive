library(readxl)
library(dplyr)
library(ggsci)
library(ggplot2)
library(data.table)
library(RColorBrewer)


enrich.test=function(x){
  a=x[1]
  b=x[2]-a
  c=x[3]-a
  d=x[4]-a-b-c
  p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value
}
enrich.demo=function(go.data,diff.list){
  
  go.data.diff=merge(go.data,subset(diff.list,select=c("Protein accession","Gene name")),by="Protein accession")
  temp.s=go.data.diff %>% group_by(Terms,Category) %>% mutate(`Gene names`=paste(`Gene name`,collapse = ";"),`Proteins`=paste(`Protein accession`,collapse = ';')) %>%
    ungroup() %>% distinct(Terms,`Gene names`,Proteins,Category)
  
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","backgroupd")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background`=length(unique(go.data$`Protein accession`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  go.data.state=subset(go.data.state,mapping>1&`Fold enrich`>1)
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value,method = "BH"),3)
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich","p.value","fdr"))
  go.data.state=merge(go.data.state,temp.s,by="Terms")
  return(go.data.state)
}



bar.plot = function(enrich.data){
  # p1 = ggplot(enrich.data, aes(log2(`Fold enrich`), reorder(Terms, -p.value), fill = p.value))+
  p1 = ggplot(enrich.data, aes(log2(`Fold enrich`), reorder(Terms, mapping), fill = p.value))+
  geom_col(width = 0.8, color = "black")+
    theme_bw(base_size = 14)+ # 建议改成 theme_bw，分面边框更好看
    labs(x = "log2(Fold enrich)", y = "")+
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60))+
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu"))+
    theme(axis.text = element_text(size = 12.5),
          strip.text = element_text(size = 14, face = "bold"))+ # 设置分面标签的字体
    # scale_x_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    # 下面这行就是分面的灵魂代码！
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") 
  
  return(p1)
}

enrich = function(data, diff.list, topn = 20) {
  enrich.data = enrich.demo(data, diff.list) %>%
  mutate(Category=case_when(
    grepl("biological", Category, ignore.case = TRUE) ~ "BP",
    grepl("cellular", Category, ignore.case = TRUE) ~ "CC",
    grepl("molecular", Category, ignore.case = TRUE) ~ "MF",
    TRUE ~ Category # 如果有不匹配的，保留原样
  )) %>%
  arrange(p.value)
  df = enrich.data %>% filter(p.value < 0.05) %>% group_by(Category) %>% slice_min(p.value, n = topn)
  
  # p1 = dot.plot(df)
  p2 = bar.plot(df)
  
  return(list(
    enrich.data = enrich.data,
    # dotplot = p1,
    barplot = p2
  ))
}

data_raw_down <- read_excel(file.choose())
# 提取 Region 列的所有不重复的字符
unique(data_raw_down$Region)
data_Common <- data_raw_down %>% filter(Region == "Common")
data_H_Unique <- data_raw_down %>% filter(Region == "H-vs-C Unique")
data_L_Unique <- data_raw_down %>% filter(Region == "L-vs-C Unique")
go_tada <- fread(file.choose())

# venn_down <- enrich(data = go_tada,
#                     diff.list = data_raw_down,
#                     topn = 10)

res_common <- enrich(data = go_tada, diff.list = data_Common, topn = 10)
res_common$barplot # 在右下角显示 Common 的图
ggsave(filename = "GO_Barplot_Common.pdf", plot = res_common$barplot, width = 10, height = 8)
ggsave(filename = "GO_Barplot_Common.png", plot = res_common$barplot, width = 10, height = 8)
# venn_down$barplot
# enrich_data <- venn_down$enrich.data
