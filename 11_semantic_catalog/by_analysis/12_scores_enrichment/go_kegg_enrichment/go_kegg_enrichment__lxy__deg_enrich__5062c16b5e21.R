#!/usr/bin/env Rscript
#=====================================================================================
# 单细胞各细胞亚型批量差异基因 + GO/KEGG富集分析脚本
# 命名参数模式，无需记忆参数顺序，原分析阈值完全保持不变
# 依赖包：Seurat, tidyverse, clusterProfiler, org.Mm.eg.db/org.Hs.eg.db, enrichplot, EnhancedVolcano, ggsci, ggrepel, argparse
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2.rds --species human -g Group -c annotation2 --ident1 Mild --ident2 Normal -o Mild_vs_Normal_DEG &
#=====================================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(EnhancedVolcano)
  library(ggsci)
  library(ggrepel)
  library(argparse)
})

# -------------------------- argparse 参数定义（替换原optparse） --------------------------
parser <- ArgumentParser(
  description = "批量每个细胞亚型做 ident1 vs ident2 差异基因，批量GO/KEGG富集分析"
)

parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "【必填】Seurat对象rds文件路径", metavar = "FILE")
parser$add_argument("--species", type = "character", required = TRUE,
                    help = "【必填】物种，仅支持 human / mouse", metavar = "STR")
parser$add_argument("-g", "--group_col", type = "character", required = TRUE,
                    help = "【必填】meta中分组列名", metavar = "STR")
parser$add_argument("-c", "--celltype_col", type = "character", required = TRUE,
                    help = "【必填】meta中细胞亚型列名", metavar = "STR")
parser$add_argument("--ident1", type = "character", required = TRUE,
                    help = "【必填】实验组（差异：ident1 VS ident2）", metavar = "STR")
parser$add_argument("--ident2", type = "character", required = TRUE,
                    help = "【必填】对照组（差异：ident1 VS ident2）", metavar = "STR")
parser$add_argument("-o", "--outdir", type = "character", default = "./deg_result",
                    help = "【可选】结果输出目录，默认 ./deg_result", metavar = "DIR")

args <- parser$parse_args()

rds_path     <- args$input
species_in   <- tolower(args$species)
group_col    <- args$group_col
celltype_col <- args$celltype_col
ident1       <- args$ident1
ident2       <- args$ident2
out_dir      <- args$outdir

# 物种合法性校验
if(!species_in %in% c("human","mouse")){
  cat("\n错误：--species 只能填写 human 或 mouse\n")
  quit(status=1)
}

# 匹配注释库与KEGG物种编号
if(species_in == "mouse"){
  OrgDb_pkg <- "org.Mm.eg.db"
  kegg_org  <- "mmu"
}else{
  OrgDb_pkg <- "org.Hs.eg.db"
  kegg_org  <- "hsa"
}
suppressPackageStartupMessages(library(OrgDb_pkg, character.only = TRUE))
OrgDb <- get(OrgDb_pkg)

# -------------------------- 读取Seurat对象并校验元数据 --------------------------
if(!file.exists(rds_path)){
  cat(sprintf("\n错误：输入文件不存在 -> %s\n", rds_path))
  quit(status=1)
}
cat(sprintf("正在读取Seurat对象：%s\n", rds_path))
obj <- readRDS(rds_path)

# 全局一次性合并Assay5所有分层，只做一次，后续子集不再出现维度错乱
if (inherits(obj[["RNA"]], "Assay5")) {
  message("检测到Assay5分层对象，全局一次性合并layers，规避子集维度bug")
  obj <- JoinLayers(obj)
  # 及时回收大内存
  gc(reset = TRUE)
}

meta_df <- obj@meta.data
if(!group_col %in% colnames(meta_df)){
  cat(sprintf("\n错误：元数据不存在分组列：%s\n", group_col))
  quit(status=1)
}
if(!celltype_col %in% colnames(meta_df)){
  cat(sprintf("\n错误：元数据不存在细胞亚型列：%s\n", celltype_col))
  quit(status=1)
}

group_all <- unique(meta_df[[group_col]])
if(!ident1 %in% group_all || !ident2 %in% group_all){
  cat(sprintf("\n错误：分组列不存在 ident1=%s 或 ident2=%s\n现有分组：%s\n",
              ident1, ident2, paste(group_all, collapse = ",")))
  quit(status=1)
}

# 创建输出文件夹
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("结果输出目录：%s\n", normalizePath(out_dir)))

# 获取全部细胞亚型
cell_types <- unique(as.character(meta_df[[celltype_col]]))
names(cell_types) <- cell_types

#=================================================================================
# 1. 批量各亚型 FindMarkers（原脚本参数完全未改动）
#=================================================================================
message("\n===== 开始批量计算各细胞亚型差异基因 =====")
deg_list <- imap_dfr(cell_types, function(ctype, name){
  message(sprintf("处理细胞亚型：%s", name))
  sub_obj <- subset(obj, !!sym(celltype_col) == name)
  Idents(sub_obj) <- group_col
  sub_idents <- unique(Idents(sub_obj))

  # 当前亚型缺少一组则跳过
  if(!all(c(ident1, ident2) %in% sub_idents)){
    message(sprintf("亚型 %s 分组不全，跳过", name))
    return(NULL)
  }

  de_res <- tryCatch({
    FindMarkers(
      sub_obj,
      ident.1 = ident1,
      ident.2 = ident2,
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      only.pos = FALSE
    )
  }, error = function(e) {
    message(sprintf("亚型 %s FindMarkers运行失败，跳过，错误信息：%s", name, e$message))
    return(NULL)
  })

  if (is.null(de_res) || nrow(de_res) == 0) return(NULL)
  de_res$gene <- rownames(de_res)
  de_res$celltype <- name
  return(de_res)
})

# 保存总差异基因表
deg_csv <- file.path(out_dir, "ALL_celltype_DEGs_ident1_vs_ident2.csv")
write.csv(deg_list, deg_csv, row.names = FALSE)

# 拆分列表用于富集
if(is.null(deg_list) || nrow(deg_list) == 0){
  message("未找到任何差异基因，富集分析终止")
  quit(status = 0)
}
deg_split <- deg_list %>% group_by(celltype) %>% group_split()
names(deg_split) <- unique(deg_list$celltype)

#=================================================================================
# 2. 批量 GO + KEGG 富集：分上调、下调分别富集，原富集参数完全原样保留
#=================================================================================
message("\n===== 开始批量GO/KEGG富集分析（上调、下调分开富集） =====")

walk(deg_split, function(de_tab){
  ctype_name <- unique(de_tab$celltype)
  # 清洗细胞名，作为合法文件夹名
  safe_ct <- gsub("[^a-zA-Z0-9_\\-]", "_", ctype_name)
  ct_folder <- file.path(out_dir, safe_ct)
  dir.create(ct_folder, recursive = TRUE, showWarnings = FALSE)

  # 1. 当前亚型全部差异基因总表存入自己文件夹
  deg_csv_path <- file.path(ct_folder, paste0(safe_ct, "_DEG.csv"))
  write.csv(de_tab, deg_csv_path, row.names = FALSE)
  message(sprintf("亚型 %s 总差异表已写入：%s", ctype_name, deg_csv_path))

  # -------------------------- 拆分上调、下调基因 --------------------------
  # 上调基因完整表格 + 基因向量
  up_de_tab <- de_tab %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0.25 & pct.1 > 0.1)
  up_genes <- pull(up_de_tab, gene)
  
  # 下调基因完整表格 + 基因向量
  down_de_tab <- de_tab %>%
    filter(p_val_adj < 0.05 & avg_log2FC < -0.25 & pct.2 > 0.1)
  down_genes <- pull(down_de_tab, gene)
  
  # 单独保存上调、下调DEG子表
  write.csv(up_de_tab, file.path(ct_folder, "significant_up_DEG.csv"), row.names = FALSE)
  write.csv(down_de_tab, file.path(ct_folder, "significant_down_DEG.csv"), row.names = FALSE)
  message(sprintf("  上调DEG表、下调DEG表已分别保存至 %s", ct_folder))

  # 定义批量富集函数，复用逻辑，减少重复代码
  run_enrichment <- function(gene_vec, enrich_subfolder, title_tag){
    if(length(gene_vec) < 10){
      message(sprintf("  %s：基因数量不足10个，跳过富集", title_tag))
      return()
    }
    enrich_dir <- file.path(ct_folder, enrich_subfolder)
    dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

    # GO富集
    go_obj <- enrichGO(
      gene = gene_vec,
      OrgDb = OrgDb,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )

    # KEGG富集
  kegg_obj <- tryCatch({
    eg_map <- bitr(
      gene_vec,
      OrgDb = OrgDb,
      fromType = "SYMBOL",
      toType = "ENTREZID"
    )
    entrez_ids <- discard(eg_map$ENTREZID, is.na)
    
    if(length(entrez_ids) > 0) {
      enrichKEGG(
        gene = entrez_ids,
        organism = kegg_org,
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        use_internal_data = FALSE
      )
    } else {
      NULL
    }
  }, error = function(e) {
    message(sprintf("  KEGG富集失败: %s", e$message))
    message("  跳过KEGG，继续GO分析")
    return(NULL)
  })

    # 输出GO结果
    go_df <- as.data.frame(go_obj)
    if(nrow(go_df) > 0){
      write.csv(go_df, file.path(enrich_dir, "GO_table.csv"), row.names = FALSE)
      ggsave(file.path(enrich_dir, "GO_dotplot.pdf"),
             dotplot(go_obj, showCategory = 10, label_format = 50),
             width = 8, height = 6)
      ggsave(file.path(enrich_dir, "GO_barplot.pdf"),
             barplot(go_obj, showCategory = 10, label_format = 50),
             width = 8, height = 6)
    }

  if(!is.null(kegg_obj) && nrow(as.data.frame(kegg_obj)) > 0){
    # 输出KEGG结果
    kegg_obj@result$Description <- str_remove(kegg_obj@result$Description, " - Mus musculus \\(house mouse\\)$")
    kegg_obj@result$Description <- str_remove(kegg_obj@result$Description, " - Homo sapiens \\(human\\)$")
    kegg_df <- as.data.frame(kegg_obj)
    if(nrow(kegg_df) > 0){
      write.csv(kegg_df, file.path(enrich_dir, "KEGG_table.csv"), row.names = FALSE)
      ggsave(file.path(enrich_dir, "KEGG_dotplot.pdf"),
             dotplot(kegg_obj, showCategory = 10, label_format = 50),
             width = 8, height = 6)
      ggsave(file.path(enrich_dir, "KEGG_barplot.pdf"),
             barplot(kegg_obj, showCategory = 10, label_format = 50),
             width = 8, height = 6)
    }
    message(sprintf("  %s富集完成，路径：%s", title_tag, enrich_dir))
      } else {
    message("  KEGG无显著结果或分析失败")
  }
  }

  # 分别执行上调、下调富集
  message(sprintf("开始处理亚型：%s", ctype_name))
  run_enrichment(up_genes, "up_regulate", "上调基因(ident1高表达)")
  run_enrichment(down_genes, "down_regulate", "下调基因(ident2高表达)")
})


message("\n✅ 全部分析流程执行完毕！")
message(sprintf("最终结果目录：%s", normalizePath(out_dir)))
