# ==========================================
# 1. 环境准备与加载你的核心函数 (保持完全不变)
# ==========================================
library(tidyverse)
library(openxlsx)
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
  temp.s=go.data.diff %>% group_by(Terms) %>% mutate(`Gene names`=paste(`Gene name`,collapse = ";"),`Proteins`=paste(`Protein accession`,collapse = ';')) %>%
    ungroup() %>% distinct(Terms,`Gene names`,Proteins)
  
  mapping=data.frame(table(go.data.diff$Terms))
  background=data.frame(table(go.data$Terms))
  go.data.state=merge(mapping,background,by="Var1")
  colnames(go.data.state)=c("Terms","mapping","backgroupd")
  go.data.state$`All mapping`=length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background`=length(unique(go.data$`Protein accession`))
  go.data.state$`Fold enrich`=round(go.data.state$mapping/go.data.state$backgroupd/go.data.state$`All mapping`*go.data.state$`All background`,2)
  
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
# 2. 读取数据库
# ==========================================
load_db_safe <- function(file_path) {
  df <- fread(file_path, check.names = FALSE)
  colnames(df) <- gsub("\\.", " ", colnames(df))
  idx <- grep("accession", colnames(df), ignore.case = TRUE)
  if(length(idx) > 0) colnames(df)[idx[1]] <- "Protein accession"
  return(df %>% select(`Protein accession`, Terms))
}

db_pfam <- load_db_safe(file.choose())
db_kegg <- load_db_safe(file.choose())
go_raw <- fread(file.choose(), check.names = FALSE)
colnames(go_raw) <- gsub("\\.", " ", colnames(go_raw))

# 标准化 accession 列名
idx <- grep("accession", colnames(go_raw), ignore.case = TRUE)
if(length(idx) > 0) colnames(go_raw)[idx[1]] <- "Protein accession"

# 找表示“分类”的那一列 (Category)
cat_idx <- grep("Category", colnames(go_raw), ignore.case = TRUE)
if(length(cat_idx) == 0) {
  stop("🚨 错误: 找不到 GO 分类列！你的 GO 表格中必须有一列标明了它是 BP、CC 还是 MF。")
}
cat_colname <- colnames(go_raw)[cat_idx[1]]
message(paste0("检测到 GO 分类列: [", cat_colname, "]，开始切割 BP, CC, MF..."))

# 将 GO 数据库拆分成 3 份
db_go_bp <- go_raw %>% filter(str_detect(toupper(.data[[cat_colname]]), "BIOLOGICAL|\\bBP\\b")) %>% select(`Protein accession`, Terms)
db_go_cc <- go_raw %>% filter(str_detect(toupper(.data[[cat_colname]]), "CELLULAR|\\bCC\\b")) %>% select(`Protein accession`, Terms)
db_go_mf <- go_raw %>% filter(str_detect(toupper(.data[[cat_colname]]), "MOLECULAR|\\bMF\\b")) %>% select(`Protein accession`, Terms)

# 组装数据库
db_list <- list(
  "Pfam"  = db_pfam, 
  "KEGG"  = db_kegg, 
  "GO_BP" = db_go_bp, 
  "GO_CC" = db_go_cc, 
  "GO_MF" = db_go_mf
)


# ==========================================
# 3. 交互式读取大表并提取分类名单
# ==========================================
message("请在弹出的窗口中选择【包含差异分组的 Excel 大表】...")
venn_detail <- read.xlsx(file.choose(), sheet = 1)

colnames(venn_detail) <- gsub("\\.", " ", colnames(venn_detail))

venn_detail <- venn_detail %>%
  mutate(Category = case_when(
    Region == "Common" | (str_detect(Region, "H-vs-C") & str_detect(Region, "L-vs-C")) ~ "Common",
    str_detect(Region, "H-vs-C") ~ "H_Unique",
    str_detect(Region, "L-vs-C") ~ "L_Unique",
    TRUE ~ "Other"
  ))

diff_sets <- list(
  "Common"   = venn_detail %>% filter(Category == "Common")   %>% select(`Protein accession`, `Gene name`),
  "H_Unique" = venn_detail %>% filter(Category == "H_Unique") %>% select(`Protein accession`, `Gene name`),
  "L_Unique" = venn_detail %>% filter(Category == "L_Unique") %>% select(`Protein accession`, `Gene name`)
)


# ==========================================
# 4. 自动化循环分析并保存结果 
# ==========================================
out_dir <- file.path(getwd(), "Enrichment_Results_0513")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# 循环会自动跑 3x5=15 次
for(db_name in names(db_list)) {
  for(set_name in names(diff_sets)) {
    
    message(paste0("正在分析: [", set_name, "] 分组在 [", db_name, "] 库中的富集..."))
    
    current_db <- db_list[[db_name]]
    current_set <- diff_sets[[set_name]]
    
    if(nrow(current_set) == 0) {
      message("   --> 跳过: 名单为空")
      next
    }
    
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
