suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(data.table)
  library(dplyr)
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-f","--file"), type="character"),
  make_option(c("-o","--out"),  type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

# ========================
# 1. read marker genes
# ========================
deg <- fread(opt$file)
rabbit_genes <- unique(deg$names)

# ========================
# 2. Rabbit → Human ortholog
# ========================
rabbit_mart <- useMart("ensembl", dataset="ocuniculus_gene_ensembl")

rabbit_map <- getBM(
  attributes = c("external_gene_name","ensembl_gene_id"),
  filters = "external_gene_name",
  values = rabbit_genes,
  mart = rabbit_mart
)

rabbit2human <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hsapiens_homolog_associated_gene_name"
  ),
  filters = "ensembl_gene_id",
  values = rabbit_map$ensembl_gene_id,
  mart = rabbit_mart
)

human_symbols <- unique(
  rabbit2human$hsapiens_homolog_associated_gene_name
)
human_symbols <- human_symbols[human_symbols != ""]

# ========================
# 3. Human SYMBOL → ENTREZ
# ========================
gene_df <- bitr(
  human_symbols,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)

human_entrez <- unique(gene_df$ENTREZID)

# ========================
# 4. GO enrichment
# ========================
ego <- enrichGO(
  gene = human_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  readable = TRUE
)

if (!is.null(ego)) {
  write.table(
    as.data.frame(ego),
    file = file.path(opt$out,"GO_enrichment.xls"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
}


if (nrow(ego) > 0) {
  p1 <- dotplot(ego, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", label_format = 80) +
    facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(file.path(opt$out, "GO_dotplot.pdf"), p1, width = 12, height = 9)
  ggsave(file.path(opt$out, "GO_dotplot.png"), p1, width = 12, height = 9)

  p2 <- barplot(ego, x = "Count", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", label_format = 80) +
    facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(file.path(opt$out, "GO_barplot.pdf"), p2, width = 12, height = 9)
  ggsave(file.path(opt$out, "GO_barplot.png"), p2, width = 12, height = 9)
}

# ========================
# 5. KEGG enrichment
# ========================
KEGG <- enrichKEGG(
  gene = human_entrez,
  organism = "hsa",
  use_internal_data = TRUE
)

if (!is.null(KEGG)) {
  write.table(
    as.data.frame(KEGG),
    file = file.path(opt$out,"KEGG_enrichment.xls"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
}

if (nrow(KEGG) > 0) {
  p3 <- dotplot(KEGG, x = "GeneRatio", color = "p.adjust", showCategory = 10, label_format = 100)
  ggsave(file.path(opt$out, "KEGG_dotplot.pdf"), p3, width = 12, height = 5)
  ggsave(file.path(opt$out, "KEGG_dotplot.png"), p3, width = 12, height = 5)

  p4 <- barplot(KEGG, x = "Count", color = "p.adjust", showCategory = 10, label_format = 80)
  ggsave(file.path(opt$out, "KEGG_barplot.pdf"), p4, width = 12, height = 5)
  ggsave(file.path(opt$out, "KEGG_barplot.png"), p4, width = 12, height = 5)
}

