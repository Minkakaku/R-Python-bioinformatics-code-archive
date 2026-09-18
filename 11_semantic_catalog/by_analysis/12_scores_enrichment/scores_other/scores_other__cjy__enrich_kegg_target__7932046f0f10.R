# ==========================================
# ==========================================
library(readxl)
library(dplyr)
library(ggsci)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stringr) 

# ==========================================
# ==========================================
target_keywords <- c(
  "ferroptosis",
  "apoptosis",
  "hif-1",
  "pi3k-akt",
  "foxo",
  "jak-stat",
  "tgf-beta",
  "ecm-receptor"
)

# ==========================================
enrich.test=function(x){
  a=x[1]
  b=x[2]-a
  c=x[3]-a
  d=x[4]-a-b-c
  p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value
}

enrich.demo=function(kegg.data, diff.list){
  kegg.data.diff = merge(kegg.data, subset(diff.list, select=c("Protein accession","Gene name")), by="Protein accession")
  
  temp.s = kegg.data.diff %>% 
    group_by(Terms) %>% 
    mutate(`Gene names`=paste(`Gene name`,collapse = ";"), `Proteins`=paste(`Protein accession`,collapse = ';')) %>%
    ungroup() %>% 
    distinct(Terms, `Gene names`, Proteins)
  
  mapping = data.frame(table(kegg.data.diff$Terms))
  background = data.frame(table(kegg.data$Terms))
  kegg.data.state = merge(mapping, background, by="Var1")
  colnames(kegg.data.state) = c("Terms", "mapping", "backgroupd")
  kegg.data.state$`All mapping` = length(unique(kegg.data.diff$`Protein accession`))
  kegg.data.state$`All background` = length(unique(kegg.data$`Protein accession`))
  
  # 计算比例与富集倍数
  kegg.data.state$`Fold enrich` = round(kegg.data.state$mapping/kegg.data.state$backgroupd/kegg.data.state$`All mapping`*kegg.data.state$`All background`, 2)
  kegg.data.state$ProteinsRatio = kegg.data.state$mapping / kegg.data.state$`All mapping`
  
  # P值与FDR计算
  kegg.data.state = subset(kegg.data.state, mapping>1 & `Fold enrich`>1)
  kegg.data.state$p.value = signif(apply(kegg.data.state[2:5], 1, enrich.test), 3)
  kegg.data.state$fdr = signif(p.adjust(kegg.data.state$p.value, method = "BH"), 3)
  
  # 拼装最终输出数据
  kegg.data.state = subset(kegg.data.state, select=c("Terms", "mapping", "Fold enrich", "ProteinsRatio", "p.value", "fdr"))
  kegg.data.state = merge(kegg.data.state, temp.s, by="Terms")
  
  return(kegg.data.state)
}

# ==========================================
# ==========================================
# 柱状图
bar.plot.integrated = function(enrich.data){
  p = ggplot(enrich.data, aes(x = ProteinsRatio, y = reorder(Terms, ProteinsRatio), fill = p.value))+
    geom_col(width = 0.8, color = "black")+
    theme_bw(base_size = 14)+ 
    labs(x = "Proteins Ratio", y = "", title = "")+ 
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60))+
    scale_fill_gradient(low = "red", high = "blue")+ # 最显著为红，不显著为蓝
    theme(axis.text = element_text(size = 12.5))+ 
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) 
  return(p)
}

# 气泡图
dot.plot.integrated = function(enrich.data){
  p = ggplot(enrich.data, aes(x = ProteinsRatio, y = reorder(Terms, ProteinsRatio))) +
    geom_point(aes(size = mapping, color = p.value)) +
    theme_bw(base_size = 14) + 
    labs(x = "Proteins Ratio", y = "", title = "", color = "p.value", size = "Protein number") + 
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60)) +
    # scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) + 
    scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+
    scale_size_continuous(range = c(3, 8)) +
    theme(axis.text = element_text(size = 12.5)) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) 
  return(p)
}
# ==========================================
# ==========================================
enrich = function(data, diff.list, targets) {
  enrich.data = enrich.demo(data, diff.list) 
  
  df_filtered = enrich.data %>% 
    filter(p.value < 0.05) %>% 
    arrange(p.value)
  
  search_pattern <- paste(tolower(targets), collapse = "|")
  matched_idx <- which(grepl(search_pattern, tolower(df_filtered$Terms)))
  
  if(length(matched_idx) == 0) {
    df_selected <- head(df_filtered, 10)
  } else {
    max_target_rank <- max(matched_idx)
    # 最少画 10 个，最多画到目标位置，不超过 20 个
    cutoff <- min(max(10, max_target_rank), 20)
    df_selected <- head(df_filtered, cutoff)
  }
  
  return(list(
    enrich.data = df_selected,
    # bar_integrated = bar.plot.integrated(df_selected),
    dot_integrated = dot.plot.integrated(df_selected) 
  ))
}

# ==========================================
# ==========================================
# 1. target proteins
data_raw <- fread(file.choose())

# 2. kegg annotation
kegg_tada <- fread(file.choose())

# 3. 执行富集与智能挑选
res_kegg <- enrich(data = kegg_tada, diff.list = data_raw, targets = target_keywords)

# ======== 查看并保存柱状图 ========
# print(res_kegg$bar_integrated)
# ggsave(filename = "KEGG_Barplot.pdf", plot = res_kegg$bar_integrated, width = 9, height = 8)
# ggsave(filename = "KEGG_Barplot.png", plot = res_kegg$bar_integrated, width = 9, height = 8)
# 
# ======== 查看并保存气泡图 ========
print(res_kegg$dot_integrated)
ggsave(filename = "KEGG_Dotplot.pdf", plot = res_kegg$dot_integrated, width = 9, height = 8)
ggsave(filename = "KEGG_Dotplot.png", plot = res_kegg$dot_integrated, width = 9, height = 8)

# ======== 导出画在图上的明细表 ========
plot_data_export <- res_kegg$enrich.data %>%
  select(Terms, p.value, fdr, ProteinsRatio, `Fold enrich`, mapping, `Gene names`, Proteins) %>%
  arrange(p.value)

write.csv(plot_data_export, file = "KEGG_Plotted_Pathways_Detail.csv", row.names = FALSE)