# ============================================================
# GO 分面富集图专用版
# 适用：GO_BP / GO_CC / GO_MF 分面图
# 输入：
#   1. 项目鉴定蛋白整体大表
#   2. 目标蛋白列表
#   3. GO 注释表：Protein accession + Terms + Category
# 输出：
#   GO_Facet_DotPlot.pdf
#   GO_Facet_BarPlot.pdf
#   GO_Facet_Plot_Data.xlsx
# ============================================================

library(tidyverse)
library(openxlsx)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(readxl)

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
  df <- enrich.data %>%
    slice_min(p.value, n = topn, with_ties = FALSE)
  
  if(nrow(df) == 0) return(NULL)
  
  list(
    enrich.data = enrich.data,
    plot.data = df
  )
}

# ----------------------------
# 2. 文件读取函数
# ----------------------------
standard_colnames <- function(df){
  colnames(df) <- gsub("\\.", " ", colnames(df))
  idx <- grep("accession", colnames(df), ignore.case = TRUE)
  if(length(idx) > 0) colnames(df)[idx[1]] <- "Protein accession"
  return(df)
}

# ----------------------------
# 3. 读取输入文件
# ----------------------------
## 项目鉴定蛋白整体大表
all_protein <- fread(file.choose()) %>% standard_colnames()

## 目标蛋白列表
target_protein <- fread(file.choose()) %>% standard_colnames()

## GO 注释表
go_raw <- fread(file.choose(), check.names = FALSE) %>% standard_colnames()

if(!all(c("Protein accession", "Terms", "Category") %in% colnames(go_raw))){
  stop("GO 注释表必须包含 Protein accession、Terms、Category 三列")
}

# ----------------------------
# 4. 整理背景和目标蛋白
# ----------------------------
if(!all(c("Protein accession", "Gene name") %in% colnames(all_protein))){
  stop("项目鉴定蛋白整体大表必须包含 Protein accession 和 Gene name")
}

if(!("Protein accession" %in% colnames(target_protein))){
  stop("目标蛋白列表必须包含 Protein accession")
}

target_protein <- target_protein %>%
  select(`Protein accession`) %>%
  distinct() %>%
  left_join(
    all_protein %>% select(`Protein accession`, `Gene name`) %>% distinct(),
    by = "Protein accession"
  ) %>%
  mutate(`Gene name` = ifelse(is.na(`Gene name`), `Protein accession`, `Gene name`))

background_accession <- unique(all_protein$`Protein accession`)

go_raw <- go_raw %>%
  filter(`Protein accession` %in% background_accession)

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

db_facet_list <- list(
  BP = db_go_bp,
  CC = db_go_cc,
  MF = db_go_mf
)

# ----------------------------
# 5. GO 分面富集和输出
# ----------------------------
out_dir <- file.path(getwd(), "GO_Facet_Enrichment")
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
  file.path(out_dir, "GO_Facet_Plot_Data.xlsx")
)

p_go_dot <- go_facet_dot.plot(go_facet_df)
p_go_bar <- go_facet_bar.plot(go_facet_df)

ggsave(
  file.path(out_dir, "GO_Facet_DotPlot.pdf"),
  p_go_dot,
  width = 10,
  height = 8
)

ggsave(
  file.path(out_dir, "GO_Facet_BarPlot.pdf"),
  p_go_bar,
  width = 10,
  height = 8
)