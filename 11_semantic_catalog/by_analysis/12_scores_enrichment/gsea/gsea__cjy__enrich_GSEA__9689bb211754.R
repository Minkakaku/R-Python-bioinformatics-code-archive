if (length(snakemake@log) > 0) {
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")
}

work_dir <- getwd()
to_abs <- function(path) {
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(work_dir, path)
}

options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggridges)
# if (snakemake@params[["species"]] == "homo_sapiens") {
#   orgdb <- "org.Hs.eg.db"
#   spe1 <- "hsa"
#   spe2 <- "human"
#   library(org.Hs.eg.db)
#   library(KEGG.db)
#   pathway_id2description <- read.table(snakemake@params[["kegg_hsa_pathway_id2description"]], sep="\t", header=TRUE)
# } else if (snakemake@params[["species"]] == "mus_musculus") {
#   orgdb <- "org.Mm.eg.db"
#   spe1 <- "mmu"
#   spe2 <- "mouse"
#   library(org.Mm.eg.db)
#   library(KEGG.db)
#   pathway_id2description <- read.table(snakemake@params[["kegg_mmu_pathway_id2description"]], sep="\t", header=TRUE)
# }  
library(pathview)
library(ReactomePA)
library(data.table)
library(stringr)
library(jsonlite)

use_scientific_notation <- function(values) {
  values <- suppressWarnings(as.numeric(as.character(values)))
  finite_values <- values[is.finite(values)]
  length(finite_values) > 0 && max(finite_values) < 0.001
}

with_ontology_levels <- function(plot_df, ontology_levels = NULL) {
  if (!is.null(ontology_levels) && "ONTOLOGY" %in% colnames(plot_df)) {
    ontology_values <- unique(na.omit(as.character(plot_df$ONTOLOGY)))
    if (all(ontology_values %in% ontology_levels)) {
      plot_df$ONTOLOGY <- factor(plot_df$ONTOLOGY, levels = ontology_levels)
    }
  }
  plot_df
}

write_enrichment_bubble_json <- function(plot_df, output_file) {
  bubble_df <- plot_df %>% arrange(desc(Description))
  if (nrow(bubble_df) == 0) {
    return(invisible(NULL))
  }

  count_values <- suppressWarnings(as.numeric(as.character(bubble_df$Count)))
  gene_ratio_values <- suppressWarnings(as.numeric(as.character(bubble_df$GeneRatio)))
  size_legend_values <- as.list(round(quantile(count_values, probs = c(0.25, 0.5, 0.75, 1.0), na.rm = TRUE)))
  bubble_json <- list(
    width = 1020, height = 680, outerBorder = TRUE, margin = list(top = 30, right = 220, bottom = 70, left = 340), color = list("#B71B1B", "#4682B4"),
    dataset = list(list(source = list(pathway = as.list(as.character(bubble_df$Description)), geneRatio = as.list(round(gene_ratio_values, 6)), count = as.list(count_values), pAdjust = as.list(bubble_df$p.adjust)))),
    series = list(list(type = "pathwayEnrichmentBubble", datasetIndex = 0, encode = list(pathway = "pathway", geneRatio = "geneRatio", count = "count", pAdjust = "pAdjust"), border = list(show = TRUE), label = list(show = FALSE))),
    xAxis = list(list(type = "scaleLinear", position = "bottom", show = TRUE, name = "GeneRatio", axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"), nameTextStyle = list(color = "#000", fontFamily = "sans-serif", fontSize = "14", fontWeight = "bold"))),
    yAxis = list(list(type = "scaleBand", position = "left", show = TRUE, axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"))),
    legend = list(), sizeLegend = list(show = TRUE, title = "Count", orient = "vertical", count = 4, value = size_legend_values, left = 820, top = 220, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"),
    colorbar = list(show = TRUE, title = "p.adjust", width = 20, height = 160, site = "right", left = 820, top = 40, fontSize = 12, tickCount = 5, useScientificNotation = use_scientific_notation(bubble_df$p.adjust))
  )
  write(jsonlite::toJSON(bubble_json, auto_unbox = TRUE, pretty = TRUE), output_file)
}

write_faceted_enrichment_bubble_json <- function(plot_df, output_file, ontology_levels = NULL) {
  facet_bubble_df <- with_ontology_levels(plot_df, ontology_levels)
  if (nrow(facet_bubble_df) == 0) {
    return(invisible(NULL))
  }

  facet_bubble_df <- facet_bubble_df %>% group_by(ONTOLOGY) %>% arrange(desc(Description), .by_group = TRUE) %>% ungroup()
  count_values <- suppressWarnings(as.numeric(as.character(facet_bubble_df$Count)))
  gene_ratio_values <- suppressWarnings(as.numeric(as.character(facet_bubble_df$GeneRatio)))
  max_gene_ratio <- max(gene_ratio_values, na.rm = TRUE)
  if (!is.finite(max_gene_ratio)) {
    max_gene_ratio <- 0
  }
  p_adjust_values <- suppressWarnings(as.numeric(as.character(facet_bubble_df$p.adjust)))
  p_adjust_range <- range(p_adjust_values, na.rm = TRUE)
  if (!all(is.finite(p_adjust_range))) {
    p_adjust_range <- c(0, 1)
  }
  fb_size_values <- as.list(round(quantile(count_values, probs = c(0.33, 1.0), na.rm = TRUE)))

  facet_bubble_json <- list(
    width = 920, height = 680, outerBorder = TRUE, margin = list(top = 30, right = 140, bottom = 70, left = 310), color = list("#e76b66", "#5e9bd3"),
    dataset = list(list(source = list(description = as.list(as.character(facet_bubble_df$Description)), geneRatio = as.list(round(gene_ratio_values, 6)), pAdjust = as.list(facet_bubble_df$p.adjust), count = as.list(count_values), ontology = as.list(as.character(facet_bubble_df$ONTOLOGY))))),
    series = list(list(type = "facetedEnrichmentBubble", datasetIndex = 0, encode = list(x = "geneRatio", y = "description", size = "count", color = "pAdjust", tooltip = list("description", "geneRatio", "pAdjust", "count", "ontology")), facet = list(field = "ontology", direction = "column", wrap = 1), sort = list(enabled = TRUE, field = "geneRatio", order = "desc"), border = list(show = TRUE), label = list(show = FALSE), facetTitleStyle = list(fontSize = "12", fontWeight = "bold"))),
    xAxis = list(list(type = "scaleLinear", position = "bottom", show = TRUE, name = "GeneRatio", min = 0, max = round(max_gene_ratio * 1.1, 4), tickCount = 5, axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisGrid = list(show = TRUE), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"), nameTextStyle = list(color = "#000", fontFamily = "sans-serif", fontSize = "14", fontWeight = "bold"))),
    yAxis = list(list(type = "scaleBand", position = "left", show = TRUE, name = "", axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisGrid = list(show = TRUE), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"))),
    legend = list(), sizeLegend = list(show = TRUE, title = "Count", orient = "vertical", value = fb_size_values, right = 20, top = 250, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"),
    colorbar = list(show = TRUE, title = "p.adjust", width = 18, height = 90, site = "right", right = 20, top = 100, fontSize = 12, tickCount = 3, colors = list("#e76b66", "#5e9bd3"), reverse = FALSE, domain = as.list(p_adjust_range), useScientificNotation = use_scientific_notation(facet_bubble_df$p.adjust))
  )
  write(jsonlite::toJSON(facet_bubble_json, auto_unbox = TRUE, pretty = TRUE), output_file)
}


# Input
## DESeq2 results
DESeq2 <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE) %>% 
  dplyr::filter(padj != "NA") # remove genes with extreme log2FoldChange value caused by some strange samples
## Cutoff values for enrichment analysis plot
eggnog_annotation <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, quote = "", check.names=FALSE)

# read featureCounts annotation and create gene2symbol dictionary
annot_tsv_path <- snakemake@input[["featurecounts_annotation"]]
annot_tsv <- read.table(annot_tsv_path, sep="\t", header=TRUE, quote="", comment.char="", check.names=FALSE)
get_symbol <- function(eggnog, gtf) {
  res <- ifelse(!is.na(eggnog) & eggnog != "NA" & eggnog != "-", eggnog,
                ifelse(!is.na(gtf) & gtf != "NA" & gtf != "-", gtf, "-"))
  return(res)
}
annot_tsv$GeneSymbol <- get_symbol(annot_tsv$symbol_eggnog, annot_tsv$symbol_gtf)
gene2symbol_dict <- setNames(annot_tsv$GeneSymbol, annot_tsv$Gene)

map_genes_to_symbols <- function(gene_string, gene2symbol_dict) {
  genes <- unlist(strsplit(as.character(gene_string), "/"))
  symbols <- gene2symbol_dict[genes]
  symbols[is.na(symbols)] <- "-"
  return(paste(symbols, collapse="/"))
}

plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
## Contrast
contrast <- gsub("_diffexp.xls", "", basename(snakemake@input[[1]]))
## Output directory
outdir <- file.path(snakemake@params[["outdir"]], contrast)
go_relation <- as.data.frame(fread(snakemake@params[["go_pathway_id2description"]]))
kegg_relation <- as.data.frame(fread(snakemake@params[["kegg_pathway_id2description"]]))
kegg_pathway_dir <- to_abs(snakemake@params[["kegg_pathway_dir"]])
dir.create(kegg_pathway_dir, recursive = TRUE, showWarnings = FALSE)

# Extract gene Ensembl ID and log2FoldChange
# Convert Ensembl ID to Entrez ID
gsea_df <- DESeq2[,c("gene","log2FoldChange")]
gsea_df <- gsea_df %>% 
  dplyr::arrange(desc(log2FoldChange))
geneList <- gsea_df[, "log2FoldChange"]
names(geneList) <- as.character(gsea_df[, "gene"])

# if (snakemake@params[["species"]] == "homo_sapiens") { 
#   gene_map <- select(org.Hs.eg.db, keys=genename, keytype="ENSEMBL", columns=c("ENTREZID"))
# } else if (snakemake@params[["species"]] == "mus_musculus") {
#   gene_map <- select(org.Mm.eg.db, keys=genename, keytype="ENSEMBL", columns=c("ENTREZID"))
# }
# colnames(gene_map)[1] <- "gene"
# gsea_df <- inner_join(gene_map, DESeq2, by = "gene")
# gsea_df <- gsea_df[, c("ENTREZID", "log2FoldChange")]
# gsea_df <- na.omit(gsea_df)
# gsea_df <- dplyr::distinct(gsea_df, ENTREZID, .keep_all=TRUE)
# gsea_df <- gsea_df %>% 
#   dplyr::arrange(desc(log2FoldChange))

# # Data format for GSEA
# geneList <- gsea_df[, "log2FoldChange"]
# names(geneList) <- as.character(gsea_df[, "ENTREZID"])

# Perform GSEA enrichment analysis and plot results

# 1. GSEA GO
# 1.1 gseGO
gene_ids <- eggnog_annotation$gene_id
eggnog_lines_with_go <- eggnog_annotation$GOs!= ""
eggnog_annoations_go <- strsplit(eggnog_annotation[eggnog_lines_with_go,]$GOs, ",")
gene2go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go],
                                 times = sapply(eggnog_annoations_go, length)),
                      term = unlist(eggnog_annoations_go))
dim(gene2go)
gene2go$term[gene2go$term=="-"] <- NA
gene2go<-na.omit(gene2go)
dim(gene2go)
term2gene<-gene2go[,c(2,1)]

GO_gseresult <- GSEA(geneList,
                     pvalueCutoff = plot_padj_cutoff,
                     pAdjustMethod = "BH",
                     TERM2GENE = term2gene)
GO_gseresult@result <- merge(GO_gseresult@result,go_relation,by = "ID") %>% select(-Description.x) %>% rename(Description="Description.y")
GO_gseresult@result <- GO_gseresult@result[order(GO_gseresult@result$pvalue),]
GO_gseresult@result$GeneSymbol <- sapply(GO_gseresult@result$core_enrichment, map_genes_to_symbols, gene2symbol_dict=gene2symbol_dict)
GO_gseresult@result <- GO_gseresult@result %>% relocate(GeneSymbol, .after = core_enrichment)

# GO_gseresult <- gseGO(
#   geneList = geneList,
#   OrgDb = orgdb,
#   keyType = "ENTREZID",
#   ont = "ALL",
#   pvalueCutoff = plot_padj_cutoff,
#   eps = 0.0
# )
# GO_gseresult_read <- DOSE::setReadable(GO_gseresult, 
#                                        OrgDb = orgdb, 
#                                        keyType = 'ENTREZID') #ENTREZID to gene Symbol
dotp <- dotplot(GO_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.pdf", sep="_")), 
       plot = dotp,
       device = "pdf",
       width = 12, 
       height = 10)
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.png", sep="_")), 
       plot = dotp,
       device = "png",
       width = 12, 
       height = 10)
write_enrichment_bubble_json(dotp$data, file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.json", sep="_")))

dotp2 <- dotplot(GO_gseresult, split = "ONTOLOGY", font.size = 14, showCategory = 10) +
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot2.pdf", sep="_")),
       plot = dotp2, 
       device = "pdf",
       width = 12, 
       height = 10)
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot2.png", sep="_")), 
       plot = dotp2, 
       device = "png",
       width = 12, 
       height = 10)
write_faceted_enrichment_bubble_json(dotp2$data, file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot2.json", sep="_")), c("biological_process", "molecular_function", "cellular_component"))

write.table(GO_gseresult@result, 
            file = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_results.xls", sep="_")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
save(GO_gseresult, 
     file = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_results.Rdata", sep="_")))

# 1.2 GSEA GO plots
gene_ids <- eggnog_annotation$gene_id
eggnog_lines_with_go <- eggnog_annotation$GOs!= ""
eggnog_annoations_go <- strsplit(eggnog_annotation[eggnog_lines_with_go,]$GOs, ",")
gene2go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go],
                                 times = sapply(eggnog_annoations_go, length)),
                      term = unlist(eggnog_annoations_go))
dim(gene2go)
gene2go$term[gene2go$term=="-"] <- NA
gene2go<-na.omit(gene2go)
dim(gene2go)
gene2go = merge(gene2go,go_relation,by.x = "term", by.y = "ID")[c("Description","gene")]
colnames(gene2go) = c("term","gene")
term2gene<-gene2go

GO_gseresult <- GSEA(geneList,
                     pvalueCutoff = plot_padj_cutoff,
                     pAdjustMethod = "BH",
                     TERM2GENE = term2gene)
GO_gseresult@result <- GO_gseresult@result[order(GO_gseresult@result$pvalue),]

# ## Ridgeplot
# ridgep <- ridgeplot(GO_gseresult, showCategory = 10) +
#   scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
# ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_ridgeplot.pdf", sep="_")), 
#        plot = ridgep, 
#        device = "pdf",
#        width = 13, 
#        height = 10)
# ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_ridgeplot.png", sep="_")), 
#        plot = ridgep, 
#        device = "png",
#        width = 13, 
#        height = 10)

## GSEAplot
### merged GSEAplot for top5
gseap <- gseaplot2(GO_gseresult, 1:5, color = "red", pvalue_table = T)
png(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_gseaplot_top5.png", sep="_")), 
    width = 16, 
    height = 15, 
    units = "in", 
    res=100)
gseap
dev.off()
pdf(file = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_gseaplot_top5.pdf", sep="_")), 
    width = 16, 
    height = 15)
gseap
dev.off()

### single GSEAplot for top20
go_list <- GO_gseresult@result$ID
if (length(go_list) < 20) {
  go_num <- c(1:length(go_list))
} else {
  go_num <- c(1:20)
}

gseaplot_dir <- file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_gseaplot", sep="_"))
system(paste("mkdir", gseaplot_dir))

setwd(gseaplot_dir)
for (i in go_num) {
  go <- go_list[i]
  go <- go_relation[go_relation$Description==go,"ID"][1]
  p <- gseaplot2(GO_gseresult, geneSetID = i, color = "red", pvalue_table = T)
  
  png(filename = paste0(paste(contrast, "GSEA_GO_enrich_gseaplot", go, sep="_"), ".png"), 
      width = 12, 
      height = 10, 
      units = "in", 
      res=100)
  print(p)
  dev.off()
  
  pdf(file = paste0(paste(contrast, "GSEA_GO_enrich_gseaplot", go, sep="_"), ".pdf"), 
      width = 12, 
      height = 10)
  print(p)
  dev.off()
}
zip_file <-  paste0("../", paste(contrast, "GSEA_GO_enrich_gseaplot", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_GO_enrich_gseaplot", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)

# 2. GSEA KEGG
# 2.1 gseKEGG
gene_ids <- eggnog_annotation$gene_id
eggnog_lines_with_kegg <- eggnog_annotation$KEGG_Pathway!= ""
eggnog_annoations_kegg <- strsplit(eggnog_annotation[eggnog_lines_with_kegg,]$KEGG_Pathway, ",")
gene2kegg <- data.frame(gene = rep(gene_ids[eggnog_lines_with_kegg],
                                 times = sapply(eggnog_annoations_kegg, length)),
                      term = unlist(eggnog_annoations_kegg))
dim(gene2kegg)
gene2kegg$term[gene2kegg$term=="-"] <- NA
gene2kegg<-na.omit(gene2kegg)
dim(gene2kegg)
gene2kegg[grep("^map", gene2kegg[,"term"]),"term"]  <- NA
gene2kegg<-na.omit(gene2kegg)
dim(gene2kegg)
term2gene<-gene2kegg[,c(2,1)]

KEGG_gseresult <- GSEA(geneList,
                     pvalueCutoff = plot_padj_cutoff,
                     pAdjustMethod = "BH",
                     TERM2GENE = term2gene)
KEGG_gseresult@result <- merge(KEGG_gseresult@result,kegg_relation,by = "ID") %>% select(-Description.x) %>% rename(Description="Description.y")
KEGG_gseresult@result <- KEGG_gseresult@result[order(KEGG_gseresult@result$pvalue),]
KEGG_gseresult@result$GeneSymbol <- sapply(KEGG_gseresult@result$core_enrichment, map_genes_to_symbols, gene2symbol_dict=gene2symbol_dict)
KEGG_gseresult@result <- KEGG_gseresult@result %>% relocate(GeneSymbol, .after = core_enrichment)

dotp <- dotplot(KEGG_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.pdf", sep="_")), 
       plot = dotp,
       device = "pdf",
       width = 11, 
       height = 10)
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.png", sep="_")), 
       plot = dotp,
       device = "png",
       width = 11, 
       height = 10)
write_enrichment_bubble_json(dotp$data, file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.json", sep="_")))
# KEGG_gseresult <- gseKEGG(
#   geneList = geneList,
#   organism  = spe1, #hsa mmu
#   keyType = "ncbi-geneid",
#   pvalueCutoff = plot_padj_cutoff,
#   eps = 0.0,
#   use_internal_data = F
# )

# KEGG_gseresult@result <- merge(KEGG_gseresult@result, pathway_id2description, by.x="ID", by.y="pathway_id") %>% 
#   dplyr::select(-Description) %>% 
#   dplyr::rename(Description = pathway_description) %>%
#   dplyr::arrange(p.adjust)
# KEGG_gseresult@result <- KEGG_gseresult@result[, c(1,11,2,3,4,5,6,7,8,9,10)]
# rownames(KEGG_gseresult@result) <- KEGG_gseresult@result$ID

# KEGG_gseresult_read <- DOSE::setReadable(KEGG_gseresult, 
#                                          OrgDb = orgdb, 
#                                          keyType = 'ENTREZID') #ENTREZID to gene Symbol
write.table(KEGG_gseresult@result, 
            file = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_results.xls", sep="_")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
save(KEGG_gseresult, 
     file = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_results.Rdata", sep="_")))

# 2.2 GSEA KEGG plots
gene_ids <- eggnog_annotation$gene_id
eggnog_lines_with_kegg <- eggnog_annotation$KEGG_Pathway!= ""
eggnog_annoations_kegg <- strsplit(eggnog_annotation[eggnog_lines_with_kegg,]$KEGG_Pathway, ",")
gene2kegg <- data.frame(gene = rep(gene_ids[eggnog_lines_with_kegg],
                                 times = sapply(eggnog_annoations_kegg, length)),
                      term = unlist(eggnog_annoations_kegg))
dim(gene2kegg)
gene2kegg$term[gene2kegg$term=="-"] <- NA
gene2kegg<-na.omit(gene2kegg)
dim(gene2kegg)
gene2kegg[grep("^map", gene2kegg[,"term"]),"term"]  <- NA
gene2kegg<-na.omit(gene2kegg)
dim(gene2kegg)
gene2kegg = merge(gene2kegg,kegg_relation,by.x = "term", by.y = "ID")[c("Description","gene")]
colnames(gene2kegg) = c("term","gene")
term2gene<-gene2kegg

KEGG_gseresult <- GSEA(geneList,
                     pvalueCutoff = plot_padj_cutoff,
                     pAdjustMethod = "BH",
                     TERM2GENE = term2gene)
KEGG_gseresult@result <- KEGG_gseresult@result[order(KEGG_gseresult@result$pvalue),]
# ## Ridgeplot
# ridgep <- ridgeplot(KEGG_gseresult, showCategory = 10) +
#   scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
# ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_ridgeplot.pdf", sep="_")), 
#        plot = ridgep, 
#        device = "pdf",
#        width = 13, 
#        height = 10)
# ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_ridgeplot.png", sep="_")), 
#        plot = ridgep, 
#        device = "png",
#        width = 13, 
#        height = 10)

## GSEAplot
### merged GSEAplot for top5
gseap <- gseaplot2(KEGG_gseresult, 1:5, color = "red", pvalue_table = T)
png(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_gseaplot_top5.png", sep="_")), 
    width = 16, 
    height = 15, 
    units = "in", 
    res=100)
gseap
dev.off()
pdf(file = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_gseaplot_top5.pdf", sep="_")), 
    width = 16, 
    height = 15)
gseap
dev.off()

### single GSEAplot for top20
KEGG_list <- KEGG_gseresult@result$ID
if (length(KEGG_list) < 20) {
  KEGG_num <- c(1:length(KEGG_list))
} else {
  KEGG_num <- c(1:20)
}

gseaplot_dir <- file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_gseaplot", sep="_"))
system(paste("mkdir", gseaplot_dir))

setwd(gseaplot_dir)
for (i in KEGG_num) {
  KEGG <- KEGG_list[i]
  KEGG <- kegg_relation[kegg_relation$Description==KEGG,"ID"][1]
  p <- gseaplot2(KEGG_gseresult, geneSetID = i, color = "red", pvalue_table = T)
  
  png(filename = paste0(paste(contrast, "GSEA_KEGG_enrich_gseaplot", KEGG, sep="_"), ".png"), 
      width = 12, 
      height = 10, 
      units = "in", 
      res=100)
  print(p)
  dev.off()
  
  pdf(file = paste0(paste(contrast, "GSEA_KEGG_enrich_gseaplot", KEGG, sep="_"), ".pdf"), 
      width = 12, 
      height = 10)
  print(p)
  dev.off()
}
zip_file <-  paste0("../", paste(contrast, "GSEA_KEGG_enrich_gseaplot", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_KEGG_enrich_gseaplot", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)

### pathview
run_pathview <- function(genelist, pathway_list, spe, n1=1, n2=1, max_success=20){
  pathway_list <- unique(pathway_list[!is.na(pathway_list) & pathway_list != ""])
  success <- max(n2 - 1, 0)

  for (pathway_id in pathway_list) {
    if (success >= max_success) {
      break
    }

    ok <- tryCatch({
      withCallingHandlers({
        pathview(
          gene.data     = genelist,
          pathway.id    = pathway_id,
          species       = "ko" , # hsa mmu
          gene.idtype   = "kegg",
          kegg.native   = T,
          new.signature = F,
          kegg.dir      = kegg_pathway_dir,
          limit         = list(gene=2, cpd=1),
          bins          = list(gene=20, cpd=10))
      }, warning = function(w) {
        message("Pathview warning for ", pathway_id, ": ", conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      file.exists(paste0(pathway_id, ".pathview.png"))
    }, error = function(e) {
      message("Pathview error for ", pathway_id, ": ", conditionMessage(e))
      FALSE
    })

    if (isTRUE(ok)) {
      success <- success + 1
    }
  }

  success
}

gene_ids <- eggnog_annotation$gene_id
eggnog_lines_with_kegg <- eggnog_annotation$KEGG_Pathway!= ""
eggnog_annoations_kegg <- strsplit(eggnog_annotation[eggnog_lines_with_kegg,]$KEGG_Pathway, ",")
gene2kegg <- data.frame(gene = rep(gene_ids[eggnog_lines_with_kegg],
                                 times = sapply(eggnog_annoations_kegg, length)),
                      term = unlist(eggnog_annoations_kegg))
dim(gene2kegg)
gene2kegg$term[gene2kegg$term=="-"] <- NA
gene2kegg<-na.omit(gene2kegg)
dim(gene2kegg)
gene2kegg[grep("^map", gene2kegg[,"term"]),"term"]  <- NA
gene2kegg<-na.omit(gene2kegg)
dim(gene2kegg)
term2gene<-gene2kegg[,c(2,1)]

KEGG_gseresult <- GSEA(geneList,
                     pvalueCutoff = plot_padj_cutoff,
                     pAdjustMethod = "BH",
                     TERM2GENE = term2gene)
KEGG_gseresult@result <- merge(KEGG_gseresult@result,kegg_relation,by = "ID") %>% select(-Description.x) %>% rename(Description="Description.y")
KEGG_gseresult@result <- KEGG_gseresult@result[order(KEGG_gseresult@result$pvalue),]

gene_ids <- eggnog_annotation$gene_id
eggnog_lines_with_keggid <- eggnog_annotation$KEGG_ko!= ""
eggnog_annoations_keggid <- strsplit(eggnog_annotation[eggnog_lines_with_keggid,]$KEGG_ko, ",")
gene2keggid <- data.frame(gene = rep(gene_ids[eggnog_lines_with_keggid],
                                 times = sapply(eggnog_annoations_keggid, length)),
                      term = unlist(eggnog_annoations_keggid))
dim(gene2keggid)
gene2keggid$term[gene2keggid$term=="-"] <- NA
gene2keggid<-na.omit(gene2keggid)
dim(gene2keggid)

log2fc_df = DESeq2[,c("gene","log2FoldChange")]
log2fc_df_koid = merge(gene2keggid,log2fc_df,by = "gene")
log2fc_df_koid = log2fc_df_koid[order(log2fc_df_koid$gene),c("gene","term","log2FoldChange")]
genelist = log2fc_df_koid[,"log2FoldChange"]
names(genelist) <- gsub("^ko:", "", log2fc_df_koid$term)

KEGG_list <- KEGG_gseresult@result$ID
pathview_dir <- file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"))
dir.create(pathview_dir, recursive = TRUE, showWarnings = FALSE)
setwd(pathview_dir)

run_pathview(genelist, KEGG_list, "ko", n1=1, n2=1)

if (length(list.files(".")) == 0) {
  writeLines("No pathview plots were generated for this GSEA result.", "README.txt")
}
zip_file <-  paste0("../", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)


