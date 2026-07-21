# 从gmt文件中匹配通路，根据每个diff差异文件对每个通路进行GSEA分析
# 每个通路都绘制一张图
# 每个文件的结果保存在一个文件夹内
library(dplyr)
library(stringr)
library(fgsea)
library(ggplot2)
library(data.table)

# =========================
# 1. pathwaylist
# =========================
# pathwaylist <- list(
#   "VEGF" = c(
#     "VEGF signaling pathway",
#     "HIF-1 signaling pathway",
#     "response to hypoxia",
#     "angiogenesis",
#     "blood vessel development",
#     "vasculature development",
#     "endothelial cell proliferation",
#     "endothelial cell migration",
#     "vascular permeability"
#   ),
#   
#   "AMPK" = c(
#     "AMPK signaling pathway",
#     "mTOR signaling pathway",
#     "autophagy",
#     "mitophagy",
#     "oxidative phosphorylation",
#     "electron transport chain",
#     "respiratory electron transport chain",
#     "ATP synthesis coupled electron transport",
#     "generation of precursor metabolites and energy",
#     "TCA cycle",
#     "citrate cycle",
#     "carbon metabolism",
#     "pyruvate metabolism",
#     "glucose metabolic process",
#     "fatty acid oxidation"
#   ),
#   
#   "PPARα" = c(
#     "PPAR signaling pathway",
#     "fatty acid degradation",
#     "fatty acid oxidation",
#     "fatty acid beta-oxidation",
#     "lipid oxidation",
#     "lipid catabolic process",
#     "fatty acid metabolic process",
#     "fatty acid catabolic process",
#     "peroxisome",
#     "acyl-CoA oxidase activity",
#     "acyl-CoA metabolic process",
#     "long-chain fatty acid metabolic process",
#     "bile acid metabolic process"
#   ),
#   
#   "PPARγ" = c(
#     "PPAR signaling pathway",
#     "lipid transport",
#     "lipid localization",
#     "lipid storage",
#     "lipid droplet organization",
#     "triglyceride metabolic process",
#     "triglyceride biosynthetic process",
#     "glycerolipid metabolic process",
#     "phospholipid metabolic process",
#     "cholesterol transport",
#     "cholesterol homeostasis",
#     "lipid and atherosclerosis",
#     "adipogenesis",
#     "fat cell differentiation"
#   )
# )
# 
# # =========================
# # 2. 文件
# # =========================
All_diff_files <- list.files(
  "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/shouhou8/Diff_and_enrich/",
  pattern = "_diffgenes\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
# 
gmtfile <- "/mnt/NexData_file/home/wmx/gmt/GO_Reactome_Hallmark_Kegg_msigdb.v2026.1.Mm.symbols.gmt"
# 
out_base <- "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/shouhou8/GSEA_pathwaylist"
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
# 
# # =========================
# # 3. 工具函数
# # =========================
safe_name <- function(x) {
  x <- gsub("[/ ]", "_", x)
  x <- gsub("α", "alpha", x)
  x <- gsub("γ", "gamma", x)
  x <- gsub("[()]", "", x)
  x <- gsub("[,:;|]", "_", x)
  x <- gsub("__+", "_", x)
  x
}
# 
# normalize_text <- function(x) {
#   x <- tolower(x)
#   x <- gsub("_", " ", x)
#   x <- gsub("-", " ", x)
#   x <- gsub("\\s+", " ", x)
#   x <- trimws(x)
#   x
# }
# 
# escape_regex <- function(x) {
#   gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", x)
# }
# 
# # 边界匹配：term 前后不能是字母/数字
# is_term_match <- function(gmt_name, term) {
#   gmt_name2 <- normalize_text(gmt_name)
#   term2 <- normalize_text(term)
#   term2 <- escape_regex(term2)
#   
#   pattern <- paste0("(^|[^a-z0-9])", term2, "([^a-z0-9]|$)")
#   str_detect(gmt_name2, regex(pattern, ignore_case = TRUE))
# }
# 
# # =========================
# # 4. 读取 GMT
# # =========================
# gmt_list <- fgsea::gmtPathways(gmtfile)
# 
# gmt_names <- names(gmt_list)
# 
# # =========================
# # 5. 从 GMT 中匹配 pathwaylist
# # =========================
# selected_pathways <- list()
# pathway_anno <- data.frame()
# 
# for (main_name in names(pathwaylist)) {
#   
#   terms <- pathwaylist[[main_name]]
#   
#   for (term in terms) {
#     
#     hits <- gmt_names[sapply(gmt_names, function(x) is_term_match(x, term))]
#     
#     if (length(hits) > 0) {
#       for (h in hits) {
#         selected_pathways[[h]] <- gmt_list[[h]]
#         
#         pathway_anno <- rbind(
#           pathway_anno,
#           data.frame(
#             main_name = main_name,
#             keyword = term,
#             pathway = h,
#             stringsAsFactors = FALSE
#           )
#         )
#       }
#     }
#   }
# }
# 
# pathway_anno <- pathway_anno %>%
#   distinct(main_name, keyword, pathway, .keep_all = TRUE)
# 
# selected_pathways <- pathway_anno[c(1,2,4,8,26,27,35,37,47,57,58,59,73,77,81,82,87,89,90,92,96,98,101,102,108,109,110,113,116,122,123,129,157,158,160,164,166,167,173,174,180,181,184,185,187,188,190,195,196,197,200),]
# 
# # 获取挑选的这些基因的ID
# library(msigdbr)
# library(dplyr)
# 
# # mouse GO BP
# msig_go_bp <- msigdbr(
#   species = "Mus musculus",
#   category = "C5",
#   subcategory = "GO:BP"
# )
# 
# msig_go_mf <- msigdbr(
#   species = "Mus musculus",
#   category = "C5",
#   subcategory = "GO:MF"
# )
# 
# msig_go_cc <- msigdbr(
#   species = "Mus musculus",
#   category = "C5",
#   subcategory = "GO:CC"
# )
# 
# go_map <- bind_rows(msig_go_bp, msig_go_mf, msig_go_cc) %>%
#   select(gs_name, gs_exact_source) %>%
#   distinct()
# # 查看列名
# colnames(msig_go)
# 
# # 构建 gs_name -> GO ID 映射
# go_map <- msig_go %>%
#   select(gs_name, gs_exact_source) %>%
#   distinct()
# 
# head(go_map)
# selected_pathways <- selected_pathways %>%
#   left_join(
#     go_map,
#     by = c("pathway" = "gs_name")
#   ) %>%
#   mutate(
#     ID = ifelse(!is.na(gs_exact_source), gs_exact_source, pathway)
#   )
# 
# write.csv(
#   selected_pathways,
#   file.path(out_base, "matched_pathways_from_gmt.csv"),
#   row.names = FALSE
# )

sp <- read.csv(file = "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/shouhou8/matched_pathways_from_gmt.csv")


# =========================
#  自动识别列名
# =========================
find_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA)
  hit[1]
}

# =========================
# 从 diff 文件构建 geneList
# =========================
make_gene_list <- function(diff_file) {
  
  df <- read.csv(diff_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  gene_col <- find_first_col(df, c("gene", "Gene", "genes", "names", "Symbol"))
  fc_col <- find_first_col(df, c("log2FC", "avg_log2FC", "avg_logFC", "logfoldchanges"))
  
  if (is.na(gene_col) || is.na(fc_col)) {
    message("[SKIP] cannot find gene/logFC column: ", diff_file)
    return(NULL)
  }
  
  df$gene_tmp <- as.character(df[[gene_col]])
  df$logFC_tmp <- suppressWarnings(as.numeric(df[[fc_col]]))
  
  df <- df %>%
    filter(
      !is.na(gene_tmp),
      gene_tmp != "",
      !is.na(logFC_tmp),
      is.finite(logFC_tmp)
    )
  
  if (nrow(df) == 0) return(NULL)
  
  df$score <- df$logFC_tmp
  
  df <- df %>%
    filter(!is.na(score), is.finite(score)) %>%
    group_by(gene_tmp) %>%
    slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  gene_list <- df$score
  names(gene_list) <- df$gene_tmp
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

# =========================
# 单个 diff 文件做 GSEA
# =========================
run_gsea_one_file <- function(
    diff_file,
    selected_pathways,
    pathway_anno,
    out_base,
    minSize = 5,
    maxSize = 500
) {
  
  message("========== GSEA: ", basename(diff_file), " ==========")
  
  file_stub <- tools::file_path_sans_ext(basename(diff_file))
  out_dir <- file.path(out_base, safe_name(file_stub))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  gene_list <- make_gene_list(diff_file)
  
  if (is.null(gene_list) || length(gene_list) < 10) {
    message("[SKIP] invalid gene list: ", diff_file)
    return(NULL)
  }
  
  # 只保留在差异基因排序向量中至少有 minSize 个基因的通路
  pathway_use <- selected_pathways[
    sapply(selected_pathways, function(gs) {
      sum(gs %in% names(gene_list)) >= minSize
    })
  ]
  
  if (length(pathway_use) == 0) {
    message("[SKIP] no pathway has enough genes: ", diff_file)
    return(NULL)
  }
  
  fgsea_res <- fgsea::fgsea(
    pathways = pathway_use,
    stats = gene_list,
    minSize = minSize,
    maxSize = maxSize
  )
  
  fgsea_res <- as.data.frame(fgsea_res)
  
  if (nrow(fgsea_res) == 0) {
    message("[SKIP] empty fgsea result: ", diff_file)
    return(NULL)
  }
  
  fgsea_res <- fgsea_res %>%
    left_join(pathway_anno, by = "pathway") %>%
    arrange(padj)
  
  # leadingEdge 是 list，保存 csv 前转成字符串
  if ("leadingEdge" %in% colnames(fgsea_res)) {
    fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = "/")
  }
  
  write.csv(
    fgsea_res,
    file.path(out_dir, paste0(file_stub, "_GSEA_result.csv")),
    row.names = FALSE
  )
  
  # 每个 main_name 也保存一个结果
  for (main_name in unique(fgsea_res$main_name)) {
    tmp <- fgsea_res %>% filter(main_name == !!main_name)
    write.csv(
      tmp,
      file.path(out_dir, paste0(file_stub, "_", safe_name(main_name), "_GSEA_result.csv")),
      row.names = FALSE
    )
  }
  
  # =========================
  # 每个通路画一张图
  # =========================
  for (i in seq_len(nrow(fgsea_res))) {
    
    pathway_name <- fgsea_res$pathway[i]
    main_name <- fgsea_res$main_name[i]
    keyword <- fgsea_res$keyword[i]
    NES <- round(fgsea_res$NES[i], 3)
    pval <- signif(fgsea_res$pval[i], 3)
    padj <- signif(fgsea_res$padj[i], 3)
    
    p <- fgsea::plotEnrichment(
      pathway_use[[pathway_name]],
      gene_list
    ) +
      labs(
        title = pathway_name,
        subtitle = paste0(
          file_stub,
          " | NES=", NES,
          " | p=", pval,
          " | padj=", padj
        ),
        x = "Rank",
        y = "Enrichment score"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5)
      )
    
    fig_dir <- file.path(out_dir, safe_name(main_name))
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    
    out_pdf <- file.path(
      fig_dir,
      paste0(
        safe_name(file_stub), "__",
        safe_name(main_name), "__",
        safe_name(pathway_name),
        ".pdf"
      )
    )
    
    out_png <- file.path(
      fig_dir,
      paste0(
        safe_name(file_stub), "__",
        safe_name(main_name), "__",
        safe_name(pathway_name),
        ".png"
      )
    )
    
    ggsave(out_pdf, p, width = 7, height = 5)
    ggsave(out_png, p, width = 7, height = 5, dpi = 300)
    
    message("[SAVE] ", out_pdf)
  }
  
  return(fgsea_res)
}

message("[INFO] matched pathways: ", length(unique(sp$pathway)))
print(table(sp$main_name))



# =========================
# 批量运行
# =========================

gmt_list <- fgsea::gmtPathways(gmtfile)
gene_list <- gmt_list[names(gmt_list) %in% unique(sp$pathway)]

for (ff in All_diff_files) {
  
  res <- run_gsea_one_file(
    diff_file = ff,
    selected_pathways = gene_list,
    pathway_anno = sp,
    out_base = out_base,
    minSize = 5,
    maxSize = 500
  )
  
}

