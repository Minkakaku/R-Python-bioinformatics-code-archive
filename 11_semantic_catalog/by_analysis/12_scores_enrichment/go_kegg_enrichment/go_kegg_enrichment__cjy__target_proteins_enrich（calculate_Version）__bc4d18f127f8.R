# ============================================================
# 目的蛋白列表富集分析：一对一版本
# 输入：
#   1. 项目鉴定蛋白整体大表
#   2. Pfam 注释表
#   3. KEGG 注释表
#   4. GO 注释表
#   5. 目的蛋白列表
# 输出：
#   KEGG、GO、PFAM的富集表、dotplot、barplot
#   可选输出 GO_BP/GO_CC/GO_MF 分面图
# ============================================================

library(tidyverse)
library(openxlsx)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(readxl)

# 是否输出 GO 分面图
draw_go_facet <- TRUE

# ----------------------------
# 1. 核心富集函数
# ----------------------------
enrich.test <- function(x){
  a <- x[1]
  b <- x[2] - a
  c <- x[3] - a
  d <- x[4] - a - b - c
  fisher.test(matrix(c(a,b,c,d), nrow = 2), alternative = "greater")$p.value
}

enrich.demo <- function(go.data, diff.list){
  
  go.data.diff <- merge(
    go.data,
    subset(diff.list, select = c("Protein accession", "Gene name")),
    by = "Protein accession"
  )
  
  temp.s <- go.data.diff %>%
    group_by(Terms) %>%
    mutate(
      `Gene names` = paste(`Gene name`, collapse = ";"),
      Proteins = paste(`Protein accession`, collapse = ";")
    ) %>%
    ungroup() %>%
    distinct(Terms, `Gene names`, Proteins)
  
  mapping <- data.frame(table(go.data.diff$Terms))
  background <- data.frame(table(go.data$Terms))
  
  go.data.state <- merge(mapping, background, by = "Var1")
  colnames(go.data.state) <- c("Terms", "mapping", "background")
  
  go.data.state$`All mapping` <- length(unique(go.data.diff$`Protein accession`))
  go.data.state$`All background` <- length(unique(go.data$`Protein accession`))
  
  go.data.state$`Fold enrich` <- round(
    go.data.state$mapping / go.data.state$background /
      go.data.state$`All mapping` * go.data.state$`All background`,
    2
  )
  
  go.data.state <- subset(go.data.state, mapping > 1 & `Fold enrich` > 1)
  
  go.data.state$p.value <- signif(apply(go.data.state[2:5], 1, enrich.test), 3)
  go.data.state$fdr <- signif(p.adjust(go.data.state$p.value, method = "BH"), 3)
  
  go.data.state <- subset(
    go.data.state,
    select = c("Terms", "mapping", "Fold enrich", "p.value", "fdr")
  )
  
  go.data.state <- merge(go.data.state, temp.s, by = "Terms")
  
  return(go.data.state)
}

dot.plot <- function(enrich.data){
  ggplot(enrich.data, aes(log2(`Fold enrich`), reorder(Terms, -p.value), fill = p.value)) +
    geom_point(aes(size = mapping), shape = 21) +
    scale_size(range = c(3,8)) +
    theme_bw(base_size = 14) +
    labs(x = "log2(Fold enrich)", y = "") +
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu")) +
    # scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+  ##马卡龙色
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 80)) +
    theme(axis.text = element_text(size = 12.5)) +
    scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))
}

bar.plot <- function(enrich.data){
  ggplot(enrich.data, aes(log2(`Fold enrich`), reorder(Terms, -p.value), fill = p.value)) +
    geom_col(width = 0.8, color = "black") +
    # theme_classic(base_size = 14) +  #没右上外框
    theme_bw(base_size = 14) +
    labs(x = "log2(Fold enrich)", y = "") +
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 80)) +
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu")) +
    # scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+
    theme(axis.text = element_text(size = 12.5)) +
    scale_x_continuous(expand = c(0, 0))
}

# ----------------------------
# 2. 新增：GO 分面图函数
# ----------------------------
go_facet_bar.plot <- function(go_enrich_df){
  ggplot(go_enrich_df, aes(log2(`Fold enrich`), reorder(Terms, -p.value), fill = p.value)) +
    geom_col(width = 0.8, color = "black") +
    facet_grid(
      GO_class ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    theme_bw(base_size = 14) +
    labs(x = "log2(Fold enrich)", y = "") +
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu")) +
    # scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 80)) +
    theme(
      axis.text = element_text(size = 10),
      strip.text.y = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = "black")
    ) +
    scale_x_continuous(expand = c(0, 0))
}

go_facet_dot.plot <- function(go_enrich_df){
  ggplot(go_enrich_df, aes(log2(`Fold enrich`), reorder(Terms, -p.value), fill = p.value)) +
    geom_point(aes(size = mapping), shape = 21) +
    facet_grid(
      GO_class ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_size(range = c(3, 8)) +
    theme_bw(base_size = 14) +
    labs(x = "log2(Fold enrich)", y = "") +
    scale_fill_gradientn(colors = brewer.pal(11, "RdBu")) +
    # scale_color_gradientn(colors = c("red", "orange", "yellow", "green3", "skyblue", "darkblue"))+
    scale_y_discrete(label = function(x) stringr::str_wrap(x, width = 80)) +
    theme(
      axis.text = element_text(size = 10),
      strip.text.y = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = "black")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))
}


enrich <- function(data, diff.list, topn = 10){
  enrich.data <- enrich.demo(data, diff.list) %>% arrange(p.value)
  # df <- enrich.data %>% filter(p.value < 0.05) %>% slice_min(p.value, n = topn)
  df <- enrich.data %>%
    slice_min(p.value, n = topn,with_ties = FALSE)  ## 不筛选，直接选指定topn terms,with_ties = FALSE严格topn,pvalue一样不显示
  if(nrow(df) == 0) return(NULL)
  
  list(
    enrich.data = enrich.data,
    plot.data = df,
    dotplot = dot.plot(df),
    barplot = bar.plot(df)
  )
}

# ----------------------------
# 3. 文件读取
# ----------------------------
standard_colnames <- function(df){
  colnames(df) <- gsub("\\.", " ", colnames(df))
  idx <- grep("accession", colnames(df), ignore.case = TRUE)
  if(length(idx) > 0) colnames(df)[idx[1]] <- "Protein accession"
  return(df)
}

load_annotation <- function(file_path){
  df <- fread(file_path, check.names = FALSE)
  df <- standard_colnames(df)
  
  if(!all(c("Protein accession", "Terms") %in% colnames(df))){
    stop("注释表必须包含 Protein accession 和 Terms 两列")
  }
  
  df %>% select(`Protein accession`, Terms) %>% distinct()
}

# ----------------------------
# 4. 读取输入文件
# ----------------------------
## 蛋白鉴定背景表
all_protein <- fread(file.choose()) %>% standard_colnames()
##target proteins
target_protein <- fread(file.choose()) %>% standard_colnames()
##Pfam 注释表
db_pfam <- load_annotation(file.choose())
##KEGG 注释表
db_kegg <- load_annotation(file.choose())
##GO 注释表，需包含 Category 列
go_raw <- fread(file.choose(), check.names = FALSE) %>% standard_colnames()

if(!all(c("Protein accession", "Terms", "Category") %in% colnames(go_raw))){
  stop("GO 注释表必须包含 Protein accession、Terms、Category 三列")
}

# ----------------------------
# 5. 整理背景和目的蛋白列表
# ----------------------------
if(!all(c("Protein accession", "Gene name") %in% colnames(all_protein))){
  stop("项目鉴定蛋白整体大表必须包含 Protein accession 和 Gene name")
}

if(!("Protein accession" %in% colnames(target_protein))){
  stop("目的蛋白列表必须包含 Protein accession")
}

target_protein <- target_protein %>%
  select(`Protein accession`) %>%
  distinct() %>%
  left_join(all_protein %>% select(`Protein accession`, `Gene name`) %>% distinct(),
            by = "Protein accession") %>%
  mutate(`Gene name` = ifelse(is.na(`Gene name`), `Protein accession`, `Gene name`))

background_accession <- unique(all_protein$`Protein accession`)

db_pfam <- db_pfam %>% filter(`Protein accession` %in% background_accession)
db_kegg <- db_kegg %>% filter(`Protein accession` %in% background_accession)

go_raw <- go_raw %>% filter(`Protein accession` %in% background_accession)

db_go <- go_raw %>%
  select(`Protein accession`, Terms) %>%
  distinct()

db_go_bp <- go_raw %>%
  filter(str_detect(toupper(Category), "BIOLOGICAL|\\bBP\\b")) %>%
  select(`Protein accession`, Terms) %>%
  distinct()

db_go_cc <- go_raw %>%
  filter(str_detect(toupper(Category), "CELLULAR|\\bCC\\b")) %>%
  select(`Protein accession`, Terms) %>%
  distinct()

db_go_mf <- go_raw %>%
  filter(str_detect(toupper(Category), "MOLECULAR|\\bMF\\b")) %>%
  select(`Protein accession`, Terms) %>%
  distinct()

##输出控制，不需要输出的注释掉---------------------
db_list <- list(
  # Pfam = db_pfam,
  # KEGG = db_kegg,
  GO = db_go
)

db_facet_list <- list(
  BP = db_go_bp,
  CC = db_go_cc,
  MF = db_go_mf
)


# ----------------------------
# 6. KEGG / Pfam / 不分面 GO 富集绘图
# 根据 db_list 中包含的数据库自动循环分析
# ----------------------------

out_dir <- file.path(getwd(), "L_down")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for(db_name in names(db_list)){
  
  message("正在分析：", db_name)
  
  res <- tryCatch({
    enrich(db_list[[db_name]], target_protein, topn = 30)
  }, error = function(e){
    message("运行出错：", e$message)
    return(NULL)
  })
  
  if(is.null(res)){
    message(db_name, " 无显著富集结果，跳过绘图")
    next
  }
  
  prefix <- file.path(out_dir, db_name)
  
  write.xlsx(
    res$enrich.data,
    paste0(prefix, "_Enrichment_Table.xlsx")
  )
  
  ggsave(
    paste0(prefix, "_DotPlot.pdf"),
    res$dotplot,
    width = 10,
    height = 7
  )
  
  ggsave(
    paste0(prefix, "_BarPlot.pdf"),
    res$barplot,
    width = 10,
    height = 7
  )
}
# ----------------------------
# 7. 只生成 GO 分面图
# ----------------------------

out_dir <- file.path(getwd(), "L_down")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

go_facet_data <- list()

for(db_name in names(db_facet_list)){
  message("正在分析：", db_name)
  
  res <- enrich(db_facet_list[[db_name]], target_protein, topn = 10)
  
  if(is.null(res)){
    message(db_name, " 无显著富集结果，跳过")
    next
  }
  
  tmp <- res$plot.data
  tmp$GO_class <- db_name
  
  go_facet_data[[db_name]] <- tmp
}

go_facet_df <- bind_rows(go_facet_data)

write.xlsx(
  go_facet_df,
  file.path(out_dir, "Target_GO_Facet_Plot_Data.xlsx")
)

p_go_dot <- go_facet_dot.plot(go_facet_df)
p_go_bar <- go_facet_bar.plot(go_facet_df)

ggsave(file.path(out_dir, "GO_Facet_DotPlot.pdf"), p_go_dot, width = 10, height = 8)
ggsave(file.path(out_dir, "GO_Facet_BarPlot.pdf"), p_go_bar, width = 10, height = 8)