# 对于Cosine.R 脚本的 MP=4 的结果进行后续分析
library(Seurat)
library(tidyverse)
library(data.table)
library(fgsea)
library(ggplot2)
library(stringr)
library(tibble)
# devtools::install_github("immunogenomics/presto") # 加速差异wilcox分析
library(presto)


setwd("D:/Project/AE2507033447-08/shouhou17/")

meta = read.csv(file = "immune_and_Endo_meta.csv")
All <- readRDS("D:/Project/AE2507033447-08/shouhou6/Main_mapping.rds")
filter_tumor_EC <- subset(All,cells = meta$barcode)

# 先根据MP聚类的内容将细胞的ECType分出来
meta1 <- meta |> mutate( ECType = case_when(
                                   sub_cell_type %in% c("Artery EC",'CD2+ Hypoxia EC','Lymphatic EC','NOX5+ Hypoxia EC','PLCXD3+ Vein EC', 'SELE+ Vein EC', 'Tip EC') ~ "EC1",
                                   sub_cell_type %in% c("Neu_CCL4","Neu_PI3","Neu_IFIT2","Neu_PROK2","Neu_SSH2") ~ "EC3",
                                   sub_cell_type %in% c("Mac_CXCL10","Mac_SPP1","Mac_CCL4","Mac_MACRO",'Mac_APOC1','Mac_IL1B',"cDC2","Mono_S100A12","mregDC","pDC",'Cycling') ~ "EC4",
                                   sub_cell_type %in% c("Naive_B","CD4_Tfh","GC_B","Memory_B","Activated_B","IgM Plasma","IgA Plasma","IgG Plasma",
                                                         "CD4_Th17","Mast cells","FDC","CD8_Tc","CD8_Tm","ILC","CD4_Tn",
                                                         "CD4_Tcm","CD4_Tem","Treg","CD4_Tex","CD8_Tex","NK","NKT") ~ "EC2"
) 
)

# identical(rownames(filter_tumor_EC@meta.data),meta1$barcode)
# [1] TRUE
filter_tumor_EC@meta.data <- meta1
rownames(filter_tumor_EC@meta.data) <- meta1$barcode

Idents(filter_tumor_EC) <- "ECType"

# 参考原文的计算差异基因
# 分各个ECType计算差异基因，使用差异基因进行后续的GSEA分析然后绘图
outdir <- "D:/Project/AE2507033447-08/shouhou17/EC_GSEA/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# 读入GSEA的gmt文件
gmt_files <- c(
  "D:/gmt/c5.go.bp.v2026.1.Hs.symbols.gmt",
  "D:/gmt/c5.go.cc.v2026.1.Hs.symbols.gmt",
  "D:/gmt/c5.go.mf.v2026.1.Hs.symbols.gmt",
  "D:/gmt/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
)
read_one_gmt <- function(gmt_file) {
  x <- fgsea::gmtPathways(gmt_file)
  return(x)
}

pathway_list <- lapply(gmt_files, read_one_gmt)
pathways <- do.call(c, pathway_list)

# 去掉重复 pathway 名
pathways <- pathways[!duplicated(names(pathways))]
message("Total pathways: ", length(pathways))
names(pathways)[1:5]

# =========================
# 函数
# =========================
# 对差异基因排序进行后续的GSEA分析
make_rank_vector <- function(deg) {
  deg <- deg %>%
  mutate(
    across(all_of(fc_col), as.numeric),
    p_val_adj = as.numeric(p_val_adj)
  ) %>%
  filter(!is.na(.data[[fc_col]])) %>%
  filter(!is.na(p_val_adj)) %>%
  filter(!is.na(gene))
  
  deg$p_val_adj[deg$p_val_adj == 0] <- 1e-300
  
  # GSEA 排序分数
  # 正值代表该 ECType 上调，负值代表其它 ECType 更高
  # 同一个 gene 只保留绝对值最大的
  deg <- deg %>%
    arrange(desc(all_of(fc_col))) %>%
    distinct(gene, .keep_all = TRUE)
  
  ranks <- deg$avg_log2FC
  names(ranks) <- deg$gene
  ranks <- sort(ranks, decreasing = TRUE)
  return(ranks)
}
# 美化通路名称
clean_term_name <- function(x) {
  x <- as.character(x)
  
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GOCC_", "", x)
  x <- gsub("^GOMF_", "", x)
  x <- gsub("^GO_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("^REACTOME_", "", x)
  
  x <- gsub("_", " ", x)
  x <- stringr::str_to_sentence(tolower(x))
  
  return(x)
}
# 绘图
plot_gsea_bar_vertical <- function(
    gsea_df,
    ec_type,
    outdir,
    top_n = 5,
    padj_cutoff = 0.05,
    only_positive_NES = TRUE
) {
  
  if (nrow(gsea_df) == 0) {
    message("[SKIP] ", ec_type, ": no GSEA result")
    return(NULL)
  }
  
  df <- gsea_df %>%
    as.data.frame() %>%
    filter(!is.na(NES)) %>%
    filter(!is.na(padj))
  
  if (only_positive_NES) {
    df <- df %>% filter(NES > 0)
  }
  
  if (nrow(df) == 0) {
    message("[SKIP] ", ec_type, ": no valid positive NES result")
    return(NULL)
  }
  
  sig <- df %>%
    filter(padj < padj_cutoff)
  
  if (nrow(sig) >= top_n) {
    plot_df <- sig %>%
      arrange(desc(NES)) %>%
      head(top_n)
  } else {
    plot_df <- df %>%
      arrange(padj, desc(NES)) %>%
      head(top_n)
  }
  
  plot_df <- plot_df %>%
    arrange(desc(NES)) %>%
    mutate(
      pathway_clean = clean_term_name(pathway),
      
      # 不再强行截太短，先保留最多 100 字符
      pathway_clean = ifelse(
        nchar(pathway_clean) > 100,
        paste0(substr(pathway_clean, 1, 100), "..."),
        pathway_clean
      ),
      
      # 每 37 个字符换行
      pathway_label = stringr::str_wrap(pathway_clean, width = 37),
      
      pathway_id = paste0("P", seq_len(n())),
      pathway_id = factor(pathway_id, levels = pathway_id),
      log10padj = -log10(padj + 1e-300)
    )
  
  # 文字从接近 x 轴的位置开始
  text_y_start <- max(plot_df$NES, na.rm = TRUE) * 0.03
  library(shadowtext)
  p <- ggplot(plot_df, aes(x = pathway_id, y = NES, fill = log10padj)) +
    geom_col(width = 0.85) +
  shadowtext::geom_shadowtext(
    aes(
      y = text_y_start,
      label = pathway_label
    ),
    angle = 90,
    hjust = 0,
    vjust = 0.5,
    color = "white",
    bg.color = "black",
    bg.r = 0.12,
    fontface = "bold",
    size = 9,
    lineheight = 0.9
  ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08))
    ) +
    scale_fill_gradient(
      low = "#B8E3D5",
      high = "#1F8FA3",
      name = expression(-Log[10](P[adj]))
    ) +
    labs(
      title = ec_type,
      x = NULL,
      y = "NES"
    ) +
    theme_classic(base_size = 28) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 28),
      legend.position = "top",
      legend.title = element_text(size = 28),
      legend.text = element_text(size = 20)
    )
  
  ggsave(
    filename = file.path(outdir, paste0(ec_safe, "_GSEA_NES_barplot_vertical.pdf")),
    plot = p,
    width = max(16, nrow(plot_df) * 0.85),
    height = 10
    )
  
  ggsave(
    filename = file.path(outdir, paste0(ec_safe, "_GSEA_NES_barplot_vertical.png")),
    plot = p,
    width = max(16, nrow(plot_df) * 0.85),
    height = 10,
    dpi = 300
  )
  
  return(p)
}

# 进行差异分析及后续分析
group_col <- "ECType"
ec_types <- unique(filter_tumor_EC@meta.data[[group_col]])
ec_types <- ec_types[!is.na(ec_types)]
ec_types <- as.character(ec_types)

all_deg_list <- list()
all_gsea_list <- list()
# 差异分析参数
logfc_cutoff <- 0
min_pct <- 0.1
fc_col <- 'avg_log2FC'
obj = filter_tumor_EC

for (ec_safe in ec_types) {
  
  message("======================================")
  message("[ECType] ", ec_safe)
  message("======================================")
  
  ec_outdir <- file.path(outdir, ec_safe)
  dir.create(ec_outdir, recursive = TRUE, showWarnings = FALSE)
  
  cells_ec <- rownames(obj@meta.data)[obj@meta.data[[group_col]] == ec_safe]
  cells_other <- rownames(obj@meta.data)[obj@meta.data[[group_col]] != ec_safe]
  
  if (length(cells_ec) < 3 || length(cells_other) < 3) {
    message("[SKIP] ", ec_safe, ": cells too few")
    next
  }
  
  # -----------------------------
  # 1. DEG: 当前 ECType vs 其它 ECType，MAST分析很慢，而且有部分亚群没有结果
  # -----------------------------
  deg <- FindMarkers(
    object = obj,
    ident.1 = ec_safe,
    ident.2 = setdiff(ec_types, ec_safe),
    group.by = group_col,
    test.use = "wilcox",
    logfc.threshold = logfc_cutoff  ,
    min.pct = min_pct
  )
  
  deg <- deg %>%
    rownames_to_column("gene")
  
  deg <- deg %>%
    arrange(p_val_adj, desc(.data[[fc_col]]))
  
  fwrite(
    deg,
    file = file.path(ec_outdir, paste0(ec_safe, "_vs_other_DEG.xls")),
    sep = "\t"
  )
  
  all_deg_list[[ec_safe]] <- deg
  
  # -----------------------------
  # 2. 生成 GSEA rank
  # -----------------------------
  ranks <- make_rank_vector(deg)
  
  if (length(ranks) < 10) {
    message("[SKIP] ", ec, ": too few ranked genes")
    next
  }
  
  fwrite(
    data.frame(gene = names(ranks), score = ranks),
    file = file.path(ec_outdir, paste0(ec_safe, "_rank_for_GSEA.rnk")),
    sep = "\t",
    col.names = FALSE
  )
  
  # -----------------------------
  # 3. fgsea
  # -----------------------------
  gsea_res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 5,
    maxSize = 500,
    # nperm = 1000
  )
  
  gsea_res <- gsea_res %>%
    as.data.frame() %>%
    arrange(padj, desc(abs(NES)))
  
  fwrite(
    gsea_res,
    file = file.path(ec_outdir, paste0(ec_safe, "_fgsea_result.xls")),
    sep = "\t"
  )
  
  all_gsea_list[[ec_safe]] <- gsea_res
  
  # -----------------------------
  # 4. 绘图
  # -----------------------------
  p2 <- plot_gsea_bar_vertical(
    gsea_df = gsea_res,
    ec_type = ec_safe,
    outdir = ec_outdir,
    top_n = 10,
    padj_cutoff = 0.05,
    only_positive_NES = TRUE
  )
}
