# =========================
# CD8_Tm_C2 vs Other T cells
# DEG + volcano + GO/KEGG
# =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(stringr)
})


scRNA = readRDS("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/T_cells_R/CD8_T/sub_CD8/CD8T_sub_harmony_umap_tsne_cell_type.rds")
scRNA@meta.data <- scRNA@meta.data[, !grepl("cluster", colnames(scRNA@meta.data), ignore.case = TRUE)]
Idents(scRNA) <- 'cell_type'

scRNA@meta.data <- scRNA@meta.data %>%
  mutate(
    group = as.character(group),
    Treatment = case_when(
      substr(group, start = nchar(group), stop = nchar(group)) == "A" ~ "After",
      substr(group, start = nchar(group), stop = nchar(group)) == "B" ~ "Before",
      TRUE ~ NA_character_
    ),
    Response = case_when(
      substr(group, start = 1, stop = nchar(group) - 2) == "non_pcr" ~ "Non_PCR",
      substr(group, start = 1, stop = nchar(group) - 2) == "pcr" ~ "PCR",
      TRUE ~ NA_character_
    ),
    Cluster =  case_when(
      RNA_snn_res.0.8 == "0" ~ "CD8_Tm_C1",
      RNA_snn_res.0.8 == "3" ~ "CD8_Tm_C2",
      RNA_snn_res.0.8 == "5" ~ "CD8_Tm_C2",
      RNA_snn_res.0.8 == "6" ~ "CD8_Tm_C3",
      RNA_snn_res.0.8 == "7" ~ "CD8_Tm_C3",
      TRUE ~ cell_type
    )
  )
outdir =  "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou11/"
# =========================
# 0. 输出目录
# =========================
deg_outdir <- file.path(outdir, "CD8_Tm_C2_vs_Other_DEG_Enrich")
dir.create(deg_outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1. 构建比较分组
# =========================
scRNA@meta.data <- scRNA@meta.data %>%
  mutate(
    Cluster = as.character(Cluster),
    DEG_group = case_when(
      Cluster == "CD8_Tm_C2" ~ "CD8_Tm_C2",
      !is.na(Cluster) & Cluster != "CD8_Tm_C2" ~ "Other_CD8_Tm",
      TRUE ~ NA_character_
    )
  )

table(scRNA$DEG_group, useNA = "ifany")

cells_use <- rownames(scRNA@meta.data)[
  !is.na(scRNA@meta.data$DEG_group)
]

scRNA_deg <- subset(scRNA, cells = cells_use)

Idents(scRNA_deg) <- "DEG_group"

# =========================
# 2. 差异基因分析
# =========================
deg <- FindMarkers(
  object = scRNA_deg,
  ident.1 = "CD8_Tm_C2",
  ident.2 = "Other_CD8_Tm",
  assay = DefaultAssay(scRNA_deg),
  slot = "data",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  only.pos = FALSE
)

deg <- deg %>%
  rownames_to_column("gene") %>%
  mutate(
    p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
    p_val = ifelse(is.na(p_val), 1, p_val),
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up in CD8_Tm_C2",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down in CD8_Tm_C2",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

write.csv(
  deg,
  file.path(deg_outdir, "CD8_Tm_C2_vs_Other_DEG_all.csv"),
  row.names = FALSE,
  quote = FALSE
)

deg_sig <- deg %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)

write.csv(
  deg_sig,
  file.path(deg_outdir, "CD8_Tm_C2_vs_Other_DEG_sig.csv"),
  row.names = FALSE,
  quote = FALSE
)

# =========================
# 3. 定义 functional marker genes
# =========================

marker_gene_list <- list(
  "Naive" = c(
    "CCR7", "SELL", "TCF7", "LEF1", "IL7R", "LTB"
  ),
  
  "Activation/Effector" = c(
    "GZMB", "GZMH", "GZMK", "PRF1", "NKG7", "GNLY",
    "IFNG", "CCL3", "CCL4", "CCL5", "CX3CR1", "KLRG1"
  ),
  
  "Cytotoxicity" = c(
    "NKG7", "GNLY", "GZMB", "GZMH", "GZMA", "PRF1",
    "KLRD1", "KLRB1", "KLRK1", "FCGR3A"
  ),
  
  "Exhaustion" = c(
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
    "TOX", "ENTPD1", "LAYN", "CXCL13"
  ),
  
  "TCR signaling" = c(
    "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2",
    "LCK", "LAT", "ZAP70", "ITK"
  ),
  
  "Adhesion" = c(
    "ITGAE", "ITGA1", "ITGAL", "ITGB1", "CD44",
    "CXCR6", "CD69", "SELL"
  ),
  
  "Chemokine/Chemokine receptor" = c(
    "CCR7", "CXCR3", "CXCR4", "CXCR6", "CCR5",
    "CCL3", "CCL4", "CCL5", "XCL1", "XCL2"
  ),
  
  "IFNγ-response" = c(
    "IFNG", "STAT1", "IRF1", "GBP1", "GBP2",
    "CXCL9", "CXCL10", "CXCL11"
  ),
  
  "NF-κB signaling" = c(
    "NFKBIA", "NFKBIZ", "TNFAIP3", "JUN", "FOS",
    "RELA", "BIRC3"
  )
)


marker_anno <- data.frame(
  functional_group = rep(names(marker_gene_list), lengths(marker_gene_list)),
  gene = unlist(marker_gene_list, use.names = FALSE),
  stringsAsFactors = FALSE
) %>%
  group_by(gene) %>%
  summarise(
    functional_group = paste(unique(functional_group), collapse = "; "),
    .groups = "drop"
  )


drop_cols <- intersect(
  c(
    "functional_group",
    "functional_group.x",
    "functional_group.y",
    "is_functional_gene",
    "functional_label"
  ),
  colnames(deg)
)

deg <- deg[, !colnames(deg) %in% drop_cols, drop = FALSE]

deg <- dplyr::left_join(deg, marker_anno, by = "gene")

deg <- deg %>%
  mutate(
    is_functional_gene = !is.na(functional_group),
    functional_label = ifelse(
      is.na(functional_group),
      NA,
      sapply(strsplit(functional_group, "; "), `[`, 1)
    )
  )

# =========================
# 4. 上下调分类
# =========================

deg <- deg %>%
  mutate(
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Up-regulated",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down-regulated",
      TRUE ~ "NS"
    ),
    sig = factor(
      sig,
      levels = c("Up-regulated", "Down-regulated", "NS")
    )
  )

# =========================
# 5. 需要标注的 functional genes
# =========================
# 当前版本：所有出现在 deg 里的 functional genes 都标注

# label_genes <- deg %>%
#   filter(
#     is_functional_gene,
#     !is.na(avg_log2FC),
#     !is.na(neg_log10_padj)
#   )

# 如果标注太密，可以改用下面这个，只标显著 functional genes
label_genes <- deg %>%
  filter(
    is_functional_gene,
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    !is.na(avg_log2FC),
    !is.na(neg_log10_padj)
  )

write.csv(
  label_genes,
  file.path(deg_outdir, "functional_genes_labeled_in_volcano.csv"),
  row.names = FALSE
)

# =========================
# 6. 颜色设置
# =========================

sig_cols <- c(
  "Up-regulated" = "#D73027",
  "Down-regulated" = "#4575B4",
  "NS" = "grey80"
)

func_cols <- c(
  "Naive" = "#4DBBD5",
  "Activation/Effector" = "#E64B35",
  "Cytotoxicity" = "#F39B7F",
  "Exhaustion" = "#7E6148",
  "TCR signaling" = "#3C5488",
  "Adhesion" = "#00A087",
  "Chemokine/Chemokine receptor" = "#8491B4",
  "IFNγ-response" = "#91D1C2",
  "NF-κB signaling" = "#DC0000",
  "TRM markers" = "#A73030"
)

# =========================
# 7. 坐标范围
# =========================

x_max <- max(abs(deg$avg_log2FC), na.rm = TRUE)
x_lim <- ceiling(x_max)

if (x_lim < 1) {
  x_lim <- 1
}

y_max <- max(deg$neg_log10_padj, na.rm = TRUE)
y_arrow <- y_max * 1.08
y_text  <- y_max * 1.15

# =========================
# 8. 绘制火山图
# =========================

p_volcano <- ggplot(deg, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(
    aes(color = sig),
    size = 0.9,
    alpha = 0.75
  ) +
  geom_point(
    data = label_genes,
    aes(color = sig),
    size = 1.8,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-0.25, 0.25),
    linetype = "dashed",
    color = "grey80",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey80",
    linewidth = 0.4
  ) +
  geom_label_repel(
    data = label_genes,
    aes(label = gene, fill = functional_label),
    color = "black",
    size = 3.0,
    label.size = 0.15,
    label.padding = unit(0.12, "lines"),
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    segment.linewidth = 0.3,
    min.segment.length = 0,
    show.legend = TRUE
  ) +
  annotate(
    "segment",
    x = -x_lim * 0.95,
    xend = -0.12,
    y = y_arrow,
    yend = y_arrow,
    arrow = arrow(length = unit(0.18, "cm"), ends = "first"),
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = 0.12,
    xend = x_lim * 0.95,
    y = y_arrow,
    yend = y_arrow,
    arrow = arrow(length = unit(0.18, "cm"), ends = "last"),
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = -x_lim * 0.55,
    y = y_text,
    label = "Other CD8+ T cells",
    size = 4.2,
    color = "black"
  ) +
  annotate(
    "text",
    x = x_lim * 0.55,
    y = y_text,
    label = "CD8_Tm_C2",
    size = 4.2,
    color = "black"
  ) +
  scale_color_manual(
    values = sig_cols,
    name = "Regulation"
  ) +
  scale_fill_manual(
    values = func_cols,
    name = "Functional type",
    na.value = "grey90"
  ) +
  scale_x_continuous(
    limits = c(-x_lim, x_lim),
    breaks = pretty(c(-x_lim, x_lim), n = 5)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.22))
  ) +
  labs(
    title = NULL,
    x = expression(Log[2]~"Fold change"),
    y = expression(-Log[10]~"adjusted p value")
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    plot.margin = margin(25, 15, 10, 10)
  )

print(p_volcano)
ggsave(
  file.path(deg_outdir, "CD8_Tm_C2_vs_Other_volcano_functional_genes_label.pdf"),
  p_volcano,
  width = 12,
  height = 8
)

ggsave(
  file.path(deg_outdir, "CD8_Tm_C2_vs_Other_volcano_functional_genes_label.png"),
  p_volcano,
  width = 12,
  height = 8,
  dpi = 300
)
# =========================
# 9. 准备富集基因
# 默认用 CD8_Tm_C2 上调基因
# =========================
gene_up <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  pull(gene) %>%
  unique()

gene_down <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC < -0.25) %>%
  pull(gene) %>%
  unique()


# SYMBOL 转 ENTREZID
gene_up_entrez <- bitr(
  gene_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_down_entrez <- bitr(
  gene_down,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


# =========================
# 10. GO 富集分析：BP
# =========================
ego_up <- NULL

if (nrow(gene_up_entrez) >= 5) {
  
  ego_up <- enrichGO(
    gene = unique(gene_up_entrez$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  ego_up_df <- as.data.frame(ego_up)
  
  write.csv(
    ego_up_df,
    file.path(deg_outdir, "CD8_Tm_C2_up_GO_ALL_enrichment.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  if (nrow(ego_up_df) > 0) {
    
    ego_up_top10 <- ego_up_df %>%
      dplyr::filter(!is.na(ONTOLOGY)) %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        ONTOLOGY = factor(
          ONTOLOGY,
          levels = c("BP", "CC", "MF")
        ),
        Description = factor(
          Description,
          levels = rev(unique(Description))
        ),
        GeneRatio_num = sapply(GeneRatio, function(x) {
          eval(parse(text = x))
        }),
        log10_padj = -log10(p.adjust + 1e-300)
      )
    
    write.csv(
      ego_up_top10,
      file.path(deg_outdir, "CD8_Tm_C2_up_GO_ALL_each_ontology_top10.csv"),
      row.names = FALSE,
      quote = FALSE
    )
    
    p_go <- ggplot(
      ego_up_top10,
      aes(
        x = GeneRatio_num,
        y = Description
      )
    ) +
      geom_point(
        aes(
          size = Count,
          color = log10_padj
        ),
        alpha = 0.9
      ) +
      facet_grid(
        ONTOLOGY ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      scale_color_gradient(
        low = "#4575B4",
        high = "#D73027",
        name = expression(-log[10]~adjusted~P)
      ) +
      scale_size_continuous(
        range = c(2.5, 7),
        name = "Gene count"
      ) +
      labs(
        title = "GO enrichment: CD8_Tm_C2 up genes",
        x = "GeneRatio",
        y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black"),
        strip.text.y = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_up_GO_ALL_each_ontology_top10_dotplot.pdf"),
      p_go,
      width = 12,
      height = 9
    )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_up_GO_ALL_each_ontology_top10_dotplot.png"),
      p_go,
      width = 12,
      height = 9,
      dpi = 300
    )
    
    print(p_go)
  }
}

# =========================
# 11. KEGG 富集分析
# =========================
ekegg_up <- NULL

if (nrow(gene_up_entrez) >= 5) {
  
  ekegg_up <- enrichKEGG(
    gene = unique(gene_up_entrez$ENTREZID),
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  ekegg_up_df <- as.data.frame(ekegg_up)
  
  if (nrow(ekegg_up_df) > 0) {
    ekegg_up <- setReadable(
      ekegg_up,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    ekegg_up_df <- as.data.frame(ekegg_up)
  }
  
  write.csv(
    ekegg_up_df,
    file.path(deg_outdir, "CD8_Tm_C2_up_KEGG_enrichment.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  if (nrow(ekegg_up_df) > 0) {
    
    p_kegg <- dotplot(
      ekegg_up,
      showCategory = 5,
      font.size = 12
    ) +
      ggtitle("KEGG enrichment: CD8_Tm_C2 up genes") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_up_KEGG_dotplot_top5.pdf"),
      p_kegg,
      width = 7.5,
      height = 5
    )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_up_KEGG_dotplot_top5.png"),
      p_kegg,
      width = 7.5,
      height = 5,
      dpi = 300
    )
  }
}

# =========================
# 12. 下调基因也做一份 GO/KEGG
# =========================
ego_down <- NULL
ekegg_down <- NULL

if (nrow(gene_down_entrez) >= 5) {
  
  ego_down <- enrichGO(
    gene = unique(gene_down_entrez$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  ego_down_df <- as.data.frame(ego_down)
  
  write.csv(
    ego_down_df,
    file.path(deg_outdir, "CD8_Tm_C2_down_GO_BP_enrichment.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  if (nrow(ego_down_df) > 0) {
    
    ego_down_top10 <- ego_down_df %>%
      dplyr::filter(!is.na(ONTOLOGY)) %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        ONTOLOGY = factor(
          ONTOLOGY,
          levels = c("BP", "CC", "MF")
        ),
        Description = factor(
          Description,
          levels = rev(unique(Description))
        ),
        GeneRatio_num = sapply(GeneRatio, function(x) {
          eval(parse(text = x))
        }),
        log10_padj = -log10(p.adjust + 1e-300)
      )
    
    write.csv(
      ego_down_top10,
      file.path(deg_outdir, "CD8_Tm_C2_down_GO_ALL_each_ontology_top10.csv"),
      row.names = FALSE,
      quote = FALSE
    )
    
    p_go <- ggplot(
      ego_down_top10,
      aes(
        x = GeneRatio_num,
        y = Description
      )
    ) +
      geom_point(
        aes(
          size = Count,
          color = log10_padj
        ),
        alpha = 0.9
      ) +
      facet_grid(
        ONTOLOGY ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      scale_color_gradient(
        low = "#4575B4",
        high = "#D73027",
        name = expression(-log[10]~adjusted~P)
      ) +
      scale_size_continuous(
        range = c(2.5, 7),
        name = "Gene count"
      ) +
      labs(
        title = "GO enrichment: CD8_Tm_C2 down genes",
        x = "GeneRatio",
        y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black"),
        strip.text.y = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_down_GO_ALL_each_ontology_top10_dotplot.pdf"),
      p_go,
      width = 12,
      height = 9
    )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_down_GO_ALL_each_ontology_top10_dotplot.png"),
      p_go,
      width = 12,
      height = 9,
      dpi = 300
    )
    
    print(p_go)
  
  ekegg_down <- enrichKEGG(
    gene = unique(gene_down_entrez$ENTREZID),
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  ekegg_down_df <- as.data.frame(ekegg_down)
  
  if (nrow(ekegg_down_df) > 0) {
    ekegg_down <- setReadable(
      ekegg_down,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    ekegg_down_df <- as.data.frame(ekegg_down)
  }
  
  write.csv(
    ekegg_down_df,
    file.path(deg_outdir, "CD8_Tm_C2_down_KEGG_enrichment.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  if (nrow(ekegg_down_df) > 0) {
    
    p_kegg_down <- dotplot(
      ekegg_down,
      showCategory = 5,
      font.size = 12
    ) +
      ggtitle("KEGG enrichment: CD8_Tm_C2 down genes") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_down_KEGG_dotplot_top5.pdf"),
      p_kegg_down,
      width = 7.5,
      height = 5
    )
    
    ggsave(
      file.path(deg_outdir, "CD8_Tm_C2_down_KEGG_dotplot_top5.png"),
      p_kegg_down,
      width = 7.5,
      height = 5,
      dpi = 300
    )
  }
  }
}

