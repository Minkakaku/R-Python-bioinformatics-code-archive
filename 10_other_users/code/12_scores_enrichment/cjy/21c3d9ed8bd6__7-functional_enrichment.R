# ============================================================
# Step 7: 功能富集 (DB + Fisher's Exact Test, base R)
# - 差异前体蛋白 GO 富集
# - 差异前体蛋白 KEGG 富集
# - 气泡图展示
# ============================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

cat("========== Step 7: Functional Enrichment (DB-based Fisher's Test) ==========\n")

cfg <- load_config(snakemake@params[["config"]])
db <- load_db_annotation(cfg$annotation$db_dir)

compare_info <- read_compare_info(snakemake@input[["compare_info"]])
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- 获取universe（全部鉴定蛋白）----------------
pep_info_file <- file.path(cfg$subdirs$preprocessing, "peptide_info.tsv")
if (file.exists(pep_info_file)) {
  pep_info_all <- fread(pep_info_file)
  universe <- unique(trimws(unlist(strsplit(na.omit(pep_info_all$Protein.Ids), ";"))))
  universe <- universe[universe != ""]
} else {
  universe <- unique(db$go$`Protein accession`)
}
cat("  Universe size:", length(universe), "\n")

top_n <- cfg$enrichment$top_n %||% 20
pcut <- cfg$enrichment$pvalue_cutoff %||% 0.05

# ---------- 纯 base R Fisher 富集函数 ------------------
fisher_enrich <- function(query_proteins, term2gene_df, term2name_df,
                          bg_proteins, pvalue_cutoff = 0.05, top_n = 20,
                          prot2gene = NULL) {
  # prot2gene: optional named vector gene=protein, for Genes column
  query_proteins <- intersect(query_proteins, bg_proteins)
  N <- length(bg_proteins)
  k <- length(query_proteins)

  terms <- unique(term2gene_df$Term)
  results <- lapply(terms, function(term) {
    in_term <- term2gene_df$Gene[term2gene_df$Term == term]
    in_term_bg <- intersect(in_term, bg_proteins)
    in_term_q <- intersect(in_term, query_proteins)
    M <- length(in_term_bg)
    x <- length(in_term_q)
    if (M == 0 || x == 0) {
      return(NULL)
    }
    mat <- matrix(c(x, k - x, M - x, N - k - M + x), nrow = 2)
    pval <- tryCatch(fisher.test(mat, alternative = "greater")$p.value, error = function(e) 1)
    gene_ratio <- paste0(x, "/", k)
    hit_prots <- paste(sort(in_term_q), collapse = ";")
    hit_genes <- if (!is.null(prot2gene)) {
      gs <- prot2gene[in_term_q]
      gs <- gs[!is.na(gs) & gs != ""]
      paste(sort(unique(gs)), collapse = ";")
    } else {
      hit_prots
    }
    data.frame(
      ID = term,
      Description = ifelse(term %in% term2name_df$Term,
        term2name_df$Name[term2name_df$Term == term][1], term
      ),
      GeneRatio = gene_ratio,
      Count = x,
      Proteins = hit_prots,
      Protein.Group = hit_genes,
      pvalue = pval,
      stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }
  results$p.adjust <- p.adjust(results$pvalue, method = "BH")
  results <- results[results$pvalue < pvalue_cutoff, ]
  if (nrow(results) == 0) {
    return(NULL)
  }
  results <- results[order(results$pvalue), ]
  head(results, top_n)
}

# ---------- 构建 GO term2gene / term2name -----------------
build_go_t2g <- function(go_db) {
  # go_db 有 "ID", "Terms", "Protein accession", "Category", "GO Level"
  t2g <- data.frame(
    Term = go_db$ID,
    Gene = go_db$`Protein accession`,
    stringsAsFactors = FALSE
  )
  t2n <- unique(data.frame(
    Term = go_db$ID,
    Name = go_db$Terms,
    stringsAsFactors = FALSE
  ))
  cat_map <- unique(data.frame(
    Term = go_db$ID,
    Category = go_db$Category,
    stringsAsFactors = FALSE
  ))
  # GO Level 列（如果存在）
  level_col <- if ("GO Level" %in% colnames(go_db)) "GO Level" else if ("Level" %in% colnames(go_db)) "Level" else NULL
  level_map <- if (!is.null(level_col)) {
    unique(data.frame(
      Term = go_db$ID,
      GO_Level = go_db[[level_col]],
      stringsAsFactors = FALSE
    ))
  } else {
    NULL
  }
  list(t2g = t2g, t2n = t2n, cat_map = cat_map, level_map = level_map)
}

build_kegg_t2g <- function(kegg_db) {
  t2g <- data.frame(
    Term = kegg_db$`Pathway ID`,
    Gene = kegg_db$`Protein accession`,
    stringsAsFactors = FALSE
  )
  t2n <- unique(data.frame(
    Term = kegg_db$`Pathway ID`,
    Name = kegg_db$Terms,
    stringsAsFactors = FALSE
  ))
  # Level1 / Level2 映射
  level_map <- unique(data.frame(
    Term = kegg_db$`Pathway ID`,
    Level1 = kegg_db$`Pathway level 1`,
    Level2 = kegg_db$`Pathway level 2`,
    stringsAsFactors = FALSE
  ))
  list(t2g = t2g, t2n = t2n, level_map = level_map)
}

# 预构建（全 DB）
go_maps <- if (!is.null(db$go)) build_go_t2g(db$go) else NULL
kegg_maps <- if (!is.null(db$kegg)) build_kegg_t2g(db$kegg) else NULL

# ---------- 对每个比较组执行富集分析 -----------------------
for (i in seq_len(nrow(compare_info))) {
  label <- compare_info$Label[i]
  cat(sprintf("\n--- [%d/%d] %s ---\n", i, nrow(compare_info), label))

  diff_file <- file.path(cfg$subdirs$differential, paste0("Table02_diff_", label, ".tsv"))
  out_go <- file.path(cfg$subdirs$enrichment, paste0("Fig16_GO_enrichment_", label, ".pdf"))
  out_kegg <- file.path(cfg$subdirs$enrichment, paste0("Fig17_KEGG_enrichment_", label, ".pdf"))
  go_table_file <- file.path(cfg$subdirs$enrichment, paste0("Table03_GO_", label, ".tsv"))
  kegg_table_file <- file.path(cfg$subdirs$enrichment, paste0("Table04_KEGG_", label, ".tsv"))

  # 确保输出目录存在
  dir.create(dirname(out_go), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(diff_file)) {
    cat("  [WARN] Diff table not found:", diff_file, "\n")
    for (f in c(out_go, out_kegg)) save_empty_plot(f, "Diff table not found")
    write.table(data.frame(), go_table_file, sep = "\t", row.names = FALSE)
    write.table(data.frame(), kegg_table_file, sep = "\t", row.names = FALSE)
    next
  }

  diff_df <- fread(diff_file)
  sig_df <- diff_df[diff_df$Status != "Not Significant", ]

  if (nrow(sig_df) == 0) {
    cat("  No significant peptides, skipping.\n")
    for (f in c(out_go, out_kegg)) save_empty_plot(f, paste0("No sig. peptides: ", label))
    write.table(data.frame(), go_table_file, sep = "\t", row.names = FALSE)
    write.table(data.frame(), kegg_table_file, sep = "\t", row.names = FALSE)
    next
  }

  sig_proteins <- unique(trimws(unlist(strsplit(na.omit(sig_df$Protein.Ids), ";"))))
  sig_proteins <- sig_proteins[sig_proteins != ""]
  cat("  Significant proteins:", length(sig_proteins), "\n")

  # 构建蛋白→基因名映射（用于 Genes 列）
  prot2gene_map <- NULL
  gene_col <- if ("Protein.Group" %in% colnames(sig_df)) "Protein.Group" else NULL
  if (!is.null(gene_col)) {
    prot2gene_map <- local({
      m <- character(0)
      pid_col <- "Protein.Ids"
      for (ri in seq_len(nrow(sig_df))) {
        prots <- trimws(strsplit(as.character(sig_df[[pid_col]][ri]), ";")[[1]])
        gn <- as.character(sig_df[[gene_col]][ri])
        prots <- prots[prots != ""]
        if (length(prots) > 0 && !is.na(gn) && gn != "") m[prots] <- gn
      }
      m
    })
  }

  if (length(sig_proteins) < 3) {
    cat("  Too few proteins, skipping.\n")
    for (f in c(out_go, out_kegg)) save_empty_plot(f, "Too few proteins for enrichment")
    write.table(data.frame(), go_table_file, sep = "\t", row.names = FALSE)
    write.table(data.frame(), kegg_table_file, sep = "\t", row.names = FALSE)
    next
  }

  # === GO 富集 ===
  cat("  GO enrichment...\n")
  if (!is.null(go_maps)) {
    go_res <- fisher_enrich(sig_proteins, go_maps$t2g, go_maps$t2n,
      universe,
      pvalue_cutoff = pcut, top_n = top_n * 3,
      prot2gene = prot2gene_map
    )
    if (!is.null(go_res) && nrow(go_res) > 0) {
      go_res$Ontology <- go_maps$cat_map$Category[
        match(go_res$ID, go_maps$cat_map$Term)
      ]
      go_res$Ontology <- gsub(
        "Biological process", "BP",
        gsub(
          "Molecular function", "MF",
          gsub("Cellular component", "CC", go_res$Ontology)
        )
      )
      # 加入 GO Level
      if (!is.null(go_maps$level_map)) {
        go_res$GO_Level <- go_maps$level_map$GO_Level[
          match(go_res$ID, go_maps$level_map$Term)
        ]
      }
      # 列顺序：ID | Description | Ontology | GO_Level | GeneRatio | Count | ...
      first_cols <- intersect(
        c(
          "ID", "Description", "Ontology", "GO_Level",
          "GeneRatio", "Count", "Proteins", "Protein.Group",
          "pvalue", "p.adjust"
        ),
        colnames(go_res)
      )
      go_res <- go_res[, c(first_cols, setdiff(colnames(go_res), first_cols))]
      write.table(go_res, go_table_file, sep = "\t", row.names = FALSE, quote = FALSE)

      go_top <- do.call(rbind, lapply(c("BP", "MF", "CC"), function(ont) {
        sub <- go_res[!is.na(go_res$Ontology) & go_res$Ontology == ont, ]
        if (nrow(sub) > 0) head(sub[order(sub$pvalue), ], min(top_n, nrow(sub)))
      }))

      if (!is.null(go_top) && nrow(go_top) > 0) {
        go_top$GeneRatio_num <- sapply(go_top$GeneRatio, function(x) {
          p <- as.numeric(strsplit(x, "/")[[1]])
          if (length(p) == 2 && p[2] > 0) p[1] / p[2] else 0
        })
        go_top$Ontology <- factor(go_top$Ontology, levels = c("BP", "CC", "MF"))
        p16 <- ggplot(go_top, aes(
          x = GeneRatio_num,
          y = reorder(Description, GeneRatio_num)
        )) +
          geom_point(aes(size = Count, color = pvalue), alpha = 0.85) +
          scale_color_gradientn(
            colors = PVAL_GRADIENT, name = "P-value",
            trans = "log10"
          ) +
          scale_size_continuous(range = c(3, 8), name = "Count") +
          facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
          labs(title = paste0("GO Enrichment: ", label), x = "Gene Ratio", y = NULL) +
          theme_sci() +
          theme(
            strip.text.y = element_text(angle = 0, size = 13, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold")
          )
        save_plot(p16, out_go,
          width = 12,
          height = adaptive_height(nrow(go_top), 0.35, 6, 14)
        )
      } else {
        save_empty_plot(out_go, "No significant GO terms")
      }
    } else {
      write.table(data.frame(), go_table_file, sep = "\t", row.names = FALSE)
      save_empty_plot(out_go, "No significant GO terms")
    }
  } else {
    write.table(data.frame(), go_table_file, sep = "\t", row.names = FALSE)
    save_empty_plot(out_go, "GO database not available")
  }

  # === KEGG 富集 ===
  cat("  KEGG enrichment...\n")
  if (!is.null(kegg_maps)) {
    kegg_res <- fisher_enrich(sig_proteins, kegg_maps$t2g, kegg_maps$t2n,
      universe,
      pvalue_cutoff = pcut, top_n = top_n,
      prot2gene = prot2gene_map
    )
    if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
      # 加入 Pathway Level1 / Level2
      if (!is.null(kegg_maps$level_map)) {
        kegg_res$Pathway_Level1 <- kegg_maps$level_map$Level1[
          match(kegg_res$ID, kegg_maps$level_map$Term)
        ]
        kegg_res$Pathway_Level2 <- kegg_maps$level_map$Level2[
          match(kegg_res$ID, kegg_maps$level_map$Term)
        ]
      }
      # 列顺序：ID | Description | Pathway_Level1 | Pathway_Level2 | GeneRatio | Count | ...
      first_cols <- intersect(
        c(
          "ID", "Description", "Pathway_Level1", "Pathway_Level2",
          "GeneRatio", "Count", "Proteins", "Protein.Group",
          "pvalue", "p.adjust"
        ),
        colnames(kegg_res)
      )
      kegg_res <- kegg_res[, c(first_cols, setdiff(colnames(kegg_res), first_cols))]
      write.table(kegg_res, kegg_table_file, sep = "\t", row.names = FALSE, quote = FALSE)
      kegg_res$GeneRatio_num <- sapply(kegg_res$GeneRatio, function(x) {
        p <- as.numeric(strsplit(x, "/")[[1]])
        if (length(p) == 2 && p[2] > 0) p[1] / p[2] else 0
      })
      p17 <- ggplot(kegg_res, aes(
        x = GeneRatio_num,
        y = reorder(Description, GeneRatio_num)
      )) +
        geom_point(aes(size = Count, color = pvalue), alpha = 0.85) +
        scale_color_gradientn(
          colors = PVAL_GRADIENT, name = "P-value",
          trans = "log10"
        ) +
        scale_size_continuous(range = c(3, 8), name = "Count") +
        labs(
          title = paste0("KEGG Pathway Enrichment: ", label),
          x = "Gene Ratio", y = NULL
        ) +
        theme_sci() +
        theme(axis.text.y = element_text(size = 12, face = "bold"))
      save_plot(p17, out_kegg,
        width = 11,
        height = adaptive_height(nrow(kegg_res), 0.35, 5, 12)
      )
    } else {
      write.table(data.frame(), kegg_table_file, sep = "\t", row.names = FALSE)
      save_empty_plot(out_kegg, "No significant KEGG pathways")
    }
  } else {
    write.table(data.frame(), kegg_table_file, sep = "\t", row.names = FALSE)
    save_empty_plot(out_kegg, "KEGG database not available")
  }
}

cat("\n========== Step 7 Complete ==========\n")
