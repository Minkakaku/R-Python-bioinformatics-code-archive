library(readxl)
library(dplyr)
library(ggsci)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stringr) 
# ==========================================
# ==========================================
target_terms <- c(
  "mitochondria-associated ER membrane",
  "mitochondrial intermembrane space",
  "lysosomal lumen",
  "endosome membrane",
  "endocytic vesicle lumen",
  "membrane raft",
  "organelle envelope"
)

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
  
  # 计算 Fold enrich 和 ProteinsRatio
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  go.data.state$ProteinsRatio = go.data.state$mapping / go.data.state$`All mapping`
  
  go.data.state=subset(go.data.state,mapping>1&`Fold enrich`>1)
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value,method = "BH"),3)
  
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich", "ProteinsRatio", "p.value","fdr"))
  go.data.state=merge(go.data.state,temp.s,by="Terms")
  return(go.data.state)
}

# ==========================================
# 绘图函数 1：分面图 (Facet)
# ==========================================
bar.plot.facet = function(enrich.data){
  p = ggplot(enrich.data, aes(x = ProteinsRatio, y = reorder(Terms, ProteinsRatio), fill = p.value))+
    geom_col(width = 0.8, color = "black")+
    theme_bw(base_size = 14)+ 
    labs(x = "Proteins Ratio", y = "", title = "")+ 
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60))+
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu"))+
    theme(axis.text = element_text(size = 12.5),
          strip.text = element_text(size = 14, face = "bold"))+ 
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") 
  return(p)
}

# ==========================================
# 绘图函数 2：整合图 (Integrated, 不分面)
# ==========================================
bar.plot.integrated = function(enrich.data){
  # 整合图统一按照 ProteinsRatio 从大到小排序
  p = ggplot(enrich.data, aes(x = ProteinsRatio, y = reorder(Terms, ProteinsRatio), fill = p.value))+
    geom_col(width = 0.8, color = "black")+
    theme_bw(base_size = 14)+ 
    labs(x = "Proteins Ratio", y = "", title = "")+ 
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60))+
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu"))+
    theme(axis.text = element_text(size = 12.5))+ 
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) 
  return(p)
}

# ==========================================
# 绘图函数 3：气泡图分面 (Dot Plot - Faceted)
# ==========================================
dot.plot.facet = function(enrich.data){
  p = ggplot(enrich.data, aes(x = ProteinsRatio, y = reorder(Terms, ProteinsRatio))) +
    geom_point(aes(size = mapping, color = p.value)) + 
    theme_bw(base_size = 14) + 
    labs(x = "Proteins Ratio", y = "", title = "", color = "p.value", size = "Protein number") + 
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60)) +
    # scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) + 
    scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+
    scale_size_continuous(range = c(3, 8)) + 
    theme(axis.text = element_text(size = 12.5),
          strip.text = element_text(size = 14, face = "bold")) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") 
  return(p)
}

# ==========================================
# 绘图函数 4：气泡图整合 (Dot Plot - Integrated)
# ==========================================
dot.plot.integrated = function(enrich.data){
  p = ggplot(enrich.data, aes(x = ProteinsRatio, y = reorder(Terms, ProteinsRatio))) +
    geom_point(aes(size = mapping, color = p.value)) +
    theme_bw(base_size = 14) + 
    labs(x = "Proteins Ratio", y = "", title = "", color = "p.value", size = "Protein number") + 
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60)) +
    # scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
    scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+
    theme(axis.text = element_text(size = 12.5)) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) 
  return(p)
}

# ==========================================
# ==========================================
enrich = function(data, diff.list, targets) {
  enrich.data = enrich.demo(data, diff.list) %>%
    mutate(Category=case_when(
      grepl("biological", Category, ignore.case = TRUE) ~ "BP",
      grepl("cellular", Category, ignore.case = TRUE) ~ "CC",
      grepl("molecular", Category, ignore.case = TRUE) ~ "MF",
      TRUE ~ Category 
    ))
  
  # 挑选逻辑
  df_selected = enrich.data %>% 
    filter(p.value < 0.05) %>% 
    arrange(p.value) %>%          # 1. 首先按 p 值进行全排序
    group_by(Category) %>%        # 2. 分别在 BP, CC, MF 中进行挑选
    group_modify(~ {
      # 找出目标通路在当前 Category 中的排位索引
      matched_idx <- which(tolower(.x$Terms) %in% tolower(targets))
      
      if(length(matched_idx) == 0) {
        # 如果这个 Category 里根本没有目标通路，直接卡死前 10
        head(.x, 10)
      } else {
        # 如果有目标通路，找到最靠后的那一个的排位
        max_target_rank <- max(matched_idx)
        # 计算最终截断位置：至少10，最高到 target 的排位，但绝不能超过 20
        cutoff <- min(max(10, max_target_rank), 20)
        head(.x, cutoff)
      }
    }) %>% 
    ungroup()
  
  # 分别生成两张图
  # p_facet <- bar.plot.facet(df_selected)
  # p_integrated <- bar.plot.integrated(df_selected)
  
  return(list(
    enrich.data = df_selected,
    # plot_facet = p_facet,
    # plot_integrated = p_integrated,
    dot_facet = dot.plot.facet(df_selected),          
    dot_integrated = dot.plot.integrated(df_selected) 
  ))
}

# ==========================================
# 运行
# ==========================================
data_raw <- fread(file.choose())
go_tada <- fread(file.choose())

# 传入 target_terms 列表
res_raw <- enrich(data = go_tada, diff.list = data_raw, targets = target_terms)

# # 1. 打印并在右下角查看分面图
# print(res_raw$plot_facet)
# ggsave(filename = "GO_Barplot_Facet.pdf", plot = res_raw$plot_facet, width = 10, height = 8)
# ggsave(filename = "GO_Barplot_Facet.png", plot = res_raw$plot_facet, width = 10, height = 8)
# 
# # 2. 打印并在右下角查看整合不分面图
# print(res_raw$plot_integrated)
# # 注意：整合图因为堆叠在一起，高度可能需要适当调高，以免文字拥挤（这里设了 height = 9）
# ggsave(filename = "GO_Barplot_Integrated.pdf", plot = res_raw$plot_integrated, width = 10, height = 9)
# ggsave(filename = "GO_Barplot_Integrated.png", plot = res_raw$plot_integrated, width = 10, height = 9)
# 
#气泡图
print(res_raw$dot_facet)
print(res_raw$dot_integrated)
ggsave(filename = "GO_Dotplot_Facet.pdf", plot = res_raw$dot_facet, width = 10, height = 8)
ggsave(filename = "GO_Dotplot_Facet.png", plot = res_raw$dot_facet, width = 10, height = 8)
ggsave(filename = "GO_Dotplot_Integrated.pdf", plot = res_raw$dot_integrated, width = 10, height = 9)
ggsave(filename = "GO_Dotplot_Integrated.png", plot = res_raw$dot_integrated, width = 10, height = 9)



# ==========================================
# 导出通路明细表
# ==========================================
plot_data <- res_raw$enrich.data
plot_data_export <- plot_data %>%
  select(Category, Terms, p.value, fdr, ProteinsRatio, `Fold enrich`, mapping, `Gene names`, Proteins) %>%
  arrange(Category, p.value) 
write.csv(plot_data_export, file = "GO_Plotted_Pathways_Detail.csv", row.names = FALSE)
# library(writexl)
# write_xlsx(plot_data_export, path = "GO_Plotted_Pathways_Detail.xlsx")