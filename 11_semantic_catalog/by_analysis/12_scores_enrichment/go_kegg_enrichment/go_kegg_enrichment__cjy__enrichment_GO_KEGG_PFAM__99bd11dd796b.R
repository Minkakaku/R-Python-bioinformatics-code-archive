# ==========================================
# 1. 环境准备与加载你的核心函数
# ==========================================
library(tidyverse)
library(openxlsx)
library(ggsci)
library(ggplot2)
library(data.table)
library(RColorBrewer)

# ------------------------------------------
# 核心计算与绘图函数保持不变
# ------------------------------------------
enrich.test=function(x){
  a=x[1]
  b=x[2]-a
  c=x[3]-a
  d=x[4]-a-b-c
  p.value=fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")$p.value
}

enrich.demo=function(go.data,diff.list){
  go.data.diff=merge(go.data,subset(diff.list,select=c("Protein accession","Gene name")),by="Protein accession")
  temp.s=go.data.diff %>% group_by(Terms) %>% mutate(`Gene names`=paste(`Gene name`,collapse = ";"),`Proteins`=paste(`Protein accession`,collapse = ';')) %>%
    ungroup() %>% distinct(Terms,`Gene names`,Proteins)
  
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","backgroupd")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background`=length(unique(go.data$`Protein accession`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  
  # 注意：给 Fold enrich 加上反引号
  go.data.state=subset(go.data.state, mapping>1 & `Fold enrich`>1) 
  
  go.data.state$p.value=signif(apply(go.data.state[2:5],1,enrich.test),3)
  go.data.state$fdr=signif(p.adjust(go.data.state$p.value,method = "BH"),3)
  go.data.state=subset(go.data.state,select=c("Terms","mapping","Fold enrich","p.value","fdr"))
  go.data.state=merge(go.data.state,temp.s,by="Terms")
  return(go.data.state)
}

dot.plot = function(enrich.data){
  p1 = ggplot(enrich.data,aes(log2(`Fold enrich`),reorder(Terms,-p.value),fill = p.value))+
    geom_point(aes(size = mapping),shape = 21)+
    scale_size(range = c(3,8))+
    theme_bw(base_size = 14)+
    labs(x = "log2(Fold enrich)",y = "")+
    scale_fill_gradientn(colors = brewer.pal(11,"RdBu"))+
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60))+
    theme(axis.text = element_text(size = 12.5))+
    scale_x_continuous(expand = expansion(mult = c(0.2,0.2)))
  return(p1)
}

bar.plot = function(enrich.data){
  p1 = ggplot(enrich.data,aes(log2(`Fold enrich`),reorder(Terms,-p.value),fill = p.value))+
    geom_col(width = 0.8,color = "black")+
    theme_classic(base_size = 14)+
    labs(x = "log2(Fold enrich)",y = "")+
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 60))+
    scale_fill_gradientn(colors = brewer.pal(11,"RdBu"))+
    theme(axis.text = element_text(size = 12.5))+
    scale_x_continuous(expand = c(0,0))
  return(p1)
}

enrich = function(data, diff.list, topn = 20){
  enrich.data = enrich.demo(data,diff.list) %>% arrange(p.value)
  df = enrich.data %>% filter(p.value < 0.05) %>% slice_min(p.value,n = topn)
  
  if(nrow(df) == 0) return(NULL) 
  
  p1 = dot.plot(df)
  p2 = bar.plot(df)
  
  return(list(enrich.data = enrich.data, dotplot = p1, barplot = p2))
}


# ==========================================
# 2. 交互式读取并标准化背景注释数据库
# ==========================================
load_db_safe <- function(file_path) {
  df <- fread(file_path, check.names = FALSE)
  colnames(df) <- gsub("\\.", " ", colnames(df))
  idx <- grep("accession", colnames(df), ignore.case = TRUE)
  if(length(idx) > 0) colnames(df)[idx[1]] <- "Protein accession"
  return(df %>% select(`Protein accession`, Terms))
}

# 弹窗选择三个数据库
message("请在弹出的窗口中选择【Pfam_annotation.tsv】文件...")
db_pfam <- load_db_safe(file.choose())

message("请在弹出的窗口中选择【KEGG_annotation.tsv】文件...")
db_kegg <- load_db_safe(file.choose())

message("请在弹出的窗口中选择【GO_annotation.tsv】文件...")
db_go   <- load_db_safe(file.choose())

db_list <- list("Pfam" = db_pfam, "KEGG" = db_kegg, "GO" = db_go)


# ==========================================
# ==========================================
# 3. 交互式读取大表并提取分类名单 (修复 Bug 版)
# ==========================================
message("请在弹出的窗口中选择【包含差异分组的 Excel 大表】...")
venn_detail <- read.xlsx(file.choose(), sheet = 1)

# 标准化大表的列名，把点换成空格
colnames(venn_detail) <- gsub("\\.", " ", colnames(venn_detail))

# 【核心修复】：精准包容 Excel 里的命名格式
venn_detail <- venn_detail %>%
  mutate(Category = case_when(
    Region == "Common" | (str_detect(Region, "H-vs-C") & str_detect(Region, "L-vs-C")) ~ "Common",
    str_detect(Region, "H-vs-C") ~ "H_Unique",
    str_detect(Region, "L-vs-C") ~ "L_Unique",
    TRUE ~ "Other"
  ))

# 提取 diff_sets，确保列名完全符合 enrich.demo 的要求
diff_sets <- list(
  "Common"   = venn_detail %>% filter(Category == "Common")   %>% select(`Protein accession`, `Gene name`),
  "H_Unique" = venn_detail %>% filter(Category == "H_Unique") %>% select(`Protein accession`, `Gene name`),
  "L_Unique" = venn_detail %>% filter(Category == "L_Unique") %>% select(`Protein accession`, `Gene name`)
)

# 打印提取结果，让你确认没问题再往下走
message("--- 提取结果核对 ---")
message("Common 数量: ", nrow(diff_sets[["Common"]]))
message("H_Unique 数量: ", nrow(diff_sets[["H_Unique"]]))
message("L_Unique 数量: ", nrow(diff_sets[["L_Unique"]]))
message("--------------------")
# ==========================================
# 4. 自动化循环分析并保存结果
# ==========================================
out_dir <- file.path(getwd(), "Enrichment_Results")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("开始执行批量富集分析...")

for(db_name in names(db_list)) {
  for(set_name in names(diff_sets)) {
    
    message(paste0("正在分析: [", set_name, "] 分组在 [", db_name, "] 库中的富集..."))
    
    current_db <- db_list[[db_name]]
    current_set <- diff_sets[[set_name]]
    
    if(nrow(current_set) == 0) {
      message("   --> 跳过: 名单为空")
      next
    }
    
    # 运行你的核心分析函数
    res <- tryCatch({
      enrich(data = current_db, diff.list = current_set, topn = 20)
    }, error = function(e) {
      message("   --> 运行出错: ", e$message)
      return(NULL)
    })
    
    if(!is.null(res)) {
      prefix <- file.path(out_dir, paste0(set_name, "_", db_name))
      write.xlsx(res$enrich.data, paste0(prefix, "_Enrichment_Table.xlsx"))
      ggsave(paste0(prefix, "_DotPlot.pdf"), res$dotplot, width = 10, height = 7)
      ggsave(paste0(prefix, "_BarPlot.pdf"), res$barplot, width = 10, height = 7)
    } else {
      message(paste0("   --> 跳过: ", set_name, " 在 ", db_name, " 中无显著富集 (P<0.05) 的数据"))
    }
  }
}

message("全部富集分析完成！结果已存入文件夹：", out_dir)