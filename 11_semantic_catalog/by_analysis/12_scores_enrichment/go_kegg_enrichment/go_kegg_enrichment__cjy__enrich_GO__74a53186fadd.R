if (length(snakemake@log) > 0) {
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")
}

options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
# if (snakemake@params[["species"]] == "homo_sapiens") {
#   orgdb <- "org.Hs.eg.db"
#   library(org.Hs.eg.db)
# } else if (snakemake@params[["species"]] == "mus_musculus") {
#   orgdb <- "org.Mm.eg.db"
#   library(org.Mm.eg.db)
# }  
library(DOSE)
library(data.table)
library(stringr)
library(dplyr)
library(jsonlite)

mixedToFloat <- function(x) {
    x <- sapply(x, as.character)
    is.integer  <- grepl("^-?\\d+$", x)
    is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
    is.float <- grepl("^-?\\d+\\.\\d+$", x)
    is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
    stopifnot(all(is.integer | is.fraction | is.float | is.mixed))
    
    numbers <- strsplit(x, "[ /]")
    
    ifelse(is.integer,  as.numeric(sapply(numbers, `[`, 1)),
           ifelse(is.float,    as.numeric(sapply(numbers, `[`, 1)),
                  ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                             as.numeric(sapply(numbers, `[`, 2)),
                             as.numeric(sapply(numbers, `[`, 1)) +
                             as.numeric(sapply(numbers, `[`, 2)) /
                             as.numeric(sapply(numbers, `[`, 3)))))
}

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

write_enrichment_bar_json <- function(plot_df, output_file) {
  bar_df <- plot_df %>% arrange(desc(Description))
  if (nrow(bar_df) == 0) {
    return(invisible(NULL))
  }

  bar_json <- list(
    width = 800, height = 700, margin = list(top = 30, right = 100, bottom = 60, left = 280), color = list("#B71B1B", "#4682B4"),
    dataset = list(list(source = list(pathway = as.list(as.character(bar_df$Description)), count = as.list(bar_df$Count), pAdjust = as.list(bar_df$p.adjust)))),
    xAxis = list(list(type = "scaleLinear", position = "bottom", show = TRUE, name = "Count", axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, fontSize = "11"), nameTextStyle = list(fontSize = "14", fontWeight = "bold"))),
    yAxis = list(list(type = "scaleBand", position = "left", show = TRUE, axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, fontSize = "11"))),
    series = list(list(type = "pathwayEnrichmentBar", datasetIndex = 0, encode = list(pathway = "pathway", count = "count", pAdjust = "pAdjust"))),
    legend = list(list(show = TRUE, title = "p.adjust", right = 30, top = 50, width = 20, height = 150, nameTextStyle = list(fontSize = "12"), textStyle = list(fontSize = "10"), useScientificNotation = use_scientific_notation(bar_df$p.adjust)))
  )
  write(jsonlite::toJSON(bar_json, auto_unbox = TRUE, pretty = TRUE), output_file)
}

write_faceted_enrichment_bar_json <- function(plot_df, output_file, ontology_levels = NULL) {
  facet_bar_df <- with_ontology_levels(plot_df, ontology_levels)
  if (nrow(facet_bar_df) == 0) {
    return(invisible(NULL))
  }

  facet_bar_df <- facet_bar_df %>% group_by(ONTOLOGY) %>% arrange(desc(Description), .by_group = TRUE) %>% ungroup()
  facet_bar_json <- list(
    width = 800, height = 500, margin = list(top = 40, right = 80, bottom = 60, left = 240), color = list("#B71B1B", "#4682B4"),
    dataset = list(list(source = list(pathway = as.list(as.character(facet_bar_df$Description)), count = as.list(facet_bar_df$Count), pAdjust = as.list(facet_bar_df$p.adjust), ontology = as.list(as.character(facet_bar_df$ONTOLOGY))))),
    xAxis = list(list(type = "scaleLinear", show = TRUE, axisGrid = list(show = TRUE), name = "Count")),
    yAxis = list(list(type = "scaleBand", show = TRUE, axisGrid = list(show = TRUE))),
    series = list(list(type = "bar", barType = "facetedEnrichment", datasetIndex = 0, encode = list(y = "pathway", x = "count", color = "pAdjust"), facet = list(field = "ontology", direction = "column", wrap = 1), label = list(show = FALSE), border = list(show = TRUE))),
    legend = list(list(show = TRUE, useScientificNotation = use_scientific_notation(facet_bar_df$p.adjust)))
  )
  write(jsonlite::toJSON(facet_bar_json, auto_unbox = TRUE, pretty = TRUE), output_file)
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
## DEG results
DEG <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1, check.names=FALSE)
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

go_relation <- as.data.frame(fread(snakemake@params[["go_pathway_id2description"]]))

## Cutoff values for enrichment analysis plot
plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
plot_q_cutoff <- as.double(snakemake@params[["enrich_plot_q_cutoff"]])
## Contrast
contrast <- gsub("_diffexp.DEG.xls", "", basename(snakemake@input[[1]]))
## Output directory
outdir <- file.path(snakemake@params[["outdir"]], contrast)


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

# eggnog_annotation[eggnog_annotation==""] <- NA
# gterms <- eggnog_annotation %>%
#   select(gene_id, GOs) %>%
#   na.omit()
# gene2go <- data.frame(term = character(),
#                       gene = character())
# for (row in 1:nrow(gterms)) {
#   gene_terms <- str_split(gterms[row,"GOs"], ",", simplify = FALSE)[[1]]
#   gene_id <- gterms[row, "gene_id"][[1]]
#   tmp <- data_frame(gene = rep(gene_id, length(gene_terms)),
#                     term = gene_terms)
#   gene2go <- rbind(gene2go, tmp)
# }
# dim(gene2go)
# gene2go$term[gene2go$term=="-"] <- NA
# gene2go<-na.omit(gene2go)
# dim(gene2go)
# term2gene<-gene2go[,c(2,1)]


# Extract genes from DEG
need_DEG <- DEG[,c("log2FoldChange","pvalue","padj","regulate")]
gene_up <- rownames(need_DEG[with(need_DEG, regulate == "Up"),])
gene_down <- rownames(need_DEG[with(need_DEG, regulate == "Down"),])

# Convert gene ENSEMBL ID to ENTREZ ID
gene_up_entrez <- unique(gene_up)
gene_down_entrez <- unique(gene_down)
gene_all_entrez <- unique(c(gene_up_entrez ,gene_down_entrez ))

# Perform GO enrichment analysis and plot results

sets = list("all", "up", "down")

for (set in sets) {
  
  if (set == "all") {
    gene_list = gene_all_entrez
  } else if (set == "up") {
    gene_list = gene_up_entrez
  } else if (set == "down") {
    gene_list = gene_down_entrez
  }

  go_enrich_results <- enricher(gene=gene_list,
                                pvalueCutoff = plot_padj_cutoff,
                                pAdjustMethod = "BH",
                                TERM2GENE = term2gene)

  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
  if (is.null(go_enrich_results)) {

    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, set, paste(contrast, set, "GO_enrich_results.xls", sep="_"))
    headers <- glue("ONTOLOGY\tID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tGeneSymbol\tCount\n")
    
    fileConn<-file(results_xls)
    writeLines(headers, fileConn)
    close(fileConn)
    
    # generate empty plots
    text = "差异表达基因过少，无法进行富集分析！"
    p <- ggplot() + 
      annotate("text", x = 4, y = 25, size=8, label = text) + 
      theme_void() + 
      theme(panel.background = element_rect(fill = 'white'))
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    # ggsave(p, 
    #        filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.pdf", sep="_")), 
    #        device = "pdf",
    #        width = 10, 
    #        height = 10)
    # ggsave(p, 
    #        filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.png", sep="_")), 
    #        device = "png",
    #        width = 4, 
    #        height = 4)
    
  } else {
    # GO enrichment
    go_enrich_results <- enricher(gene=gene_list,
                                  pvalueCutoff = plot_padj_cutoff,
                                  pAdjustMethod = "BH",
                                  TERM2GENE = term2gene)
    go_enrich_results@result <- merge(go_enrich_results@result,go_relation,by = "ID") %>% select(-Description.x) %>% rename(Description="Description.y")
    go_enrich_results@result <- go_enrich_results@result[order(go_enrich_results@result$pvalue),]
    go_enrich_results@result$GeneSymbol <- sapply(go_enrich_results@result$geneID, map_genes_to_symbols, gene2symbol_dict=gene2symbol_dict)
    go_enrich_results@result <- go_enrich_results@result %>% relocate(GeneSymbol, .after = geneID)
    write.table(go_enrich_results@result, 
                file = file.path(outdir, set, paste(contrast, set, "GO_enrich_results.xls", sep="_")),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(go_enrich_results, 
         file = file.path(outdir, set, paste(contrast, set, "GO_enrich_results.Rdata", sep="_")))
    
    # Plots
    ## barplot
    go_enrich_plot = head(go_enrich_results@result,30)
    go_enrich_plot = go_enrich_plot[order(go_enrich_plot$Count,decreasing = TRUE),]
    go_enrich_plot$Description = factor(go_enrich_plot$Description, levels = rev(go_enrich_plot$Description))
    barp <- ggplot(go_enrich_plot,aes(x=Count,y=Description,fill=p.adjust)) +
                geom_bar(stat="identity", width=0.8) +
                theme_bw() +
                ylab("") +
                scale_fill_gradient(high = "blue", low = "red") +
      # geom_rect(ymin = 0, ymax = 12, xmin = 0, xmax = nrow(go_enrich_df) + 1, fill = NA, color = "black") +
      # geom_text(aes(colour = "black", label = Description), size = 3.5, y = -0.2, vjust = 0.4, hjust = 1) +
      # scale_colour_manual(values = "black", guide = F) + 
                theme(axis.title = element_text(face = "plain"),
                      axis.text = element_text(face = "plain"),
                      # panel.grid.minor = element_blank(),
                      # panel.grid.major = element_blank(),
                      legend.title = element_text(size=as.numeric(16), face = "plain"),
                      legend.position="right",
                      legend.text = element_text(size = as.numeric(15)),
                      axis.text.x=element_text(size=as.numeric(15),face = "plain"),
                      axis.text.y=element_text(size=as.numeric(15),face = "plain"),
                      axis.title.x=element_text(size=as.numeric(16),face = "plain"),
                      axis.title.y=element_text(size=as.numeric(16))) +
                scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) + 
                theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 12, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 12, 
           height = 10)
    write_enrichment_bar_json(barp$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.json", sep="_")))
    ### if ont = "ALL"
    go_enrich_plot = rbind(head(go_enrich_results@result[go_enrich_results@result$ONTOLOGY == "biological_process",],10),
                           head(go_enrich_results@result[go_enrich_results@result$ONTOLOGY == "molecular_function",],10),
                           head(go_enrich_results@result[go_enrich_results@result$ONTOLOGY == "cellular_component",],10))
    go_enrich_plot = go_enrich_plot[order(go_enrich_plot$Count,decreasing = TRUE),]
    go_enrich_plot$Description = factor(go_enrich_plot$Description, levels = rev(go_enrich_plot$Description))
    go_enrich_plot$ONTOLOGY = factor(go_enrich_plot$ONTOLOGY, levels = c("biological_process","molecular_function","cellular_component"))
    barp <- ggplot(go_enrich_plot,aes(x=Count,y=Description,fill=p.adjust)) +
                geom_bar(stat="identity", width=0.8) +
                theme_bw() +
                ylab("") +
                scale_fill_gradient(high = "blue", low = "red") +
      # geom_rect(ymin = 0, ymax = 12, xmin = 0, xmax = nrow(go_enrich_df) + 1, fill = NA, color = "black") +
      # geom_text(aes(colour = "black", label = Description), size = 3.5, y = -0.2, vjust = 0.4, hjust = 1) +
      # scale_colour_manual(values = "black", guide = F) + 
                theme(axis.title = element_text(face = "plain"),
                      axis.text = element_text(face = "plain"),
                      # panel.grid.minor = element_blank(),
                      # panel.grid.major = element_blank(),
                      legend.title = element_text(size=as.numeric(16), face = "plain"),
                      legend.position="right",
                      legend.text = element_text(size = as.numeric(15)),
                      axis.text.x=element_text(size=as.numeric(15),face = "plain"),
                      axis.text.y=element_text(size=as.numeric(15),face = "plain"),
                      axis.title.x=element_text(size=as.numeric(16),face = "plain"),
                      axis.title.y=element_text(size=as.numeric(16))) +
                scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) + 
                theme(plot.margin = unit(c(1,1,1,1),'lines')) +
                facet_grid(ONTOLOGY~., scale = "free")
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.pdf", sep="_")), 
           device = "pdf",
           width = 12, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.png", sep="_")), 
           device = "png",
           width = 12, 
           height = 10)
    write_faceted_enrichment_bar_json(barp$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.json", sep="_")), c("biological_process", "molecular_function", "cellular_component"))
    
    ## dotplot 
    go_enrich_plot = head(go_enrich_results@result,30)
    go_enrich_plot$GeneRatio <- mixedToFloat(go_enrich_plot$GeneRatio)
    go_enrich_plot = go_enrich_plot[order(go_enrich_plot$GeneRatio,decreasing = TRUE),]
    go_enrich_plot$Description = factor(go_enrich_plot$Description, levels = rev(go_enrich_plot$Description))
    dotp <- ggplot(go_enrich_plot) +
                geom_point(aes(x=GeneRatio, y=Description, color=p.adjust, size = Count), shape = 19)+
                theme_bw() +
                ylab("") +
                scale_color_gradient(high = "blue", low = "red") +
      # geom_rect(ymin = 0, ymax = 12, xmin = 0, xmax = nrow(go_enrich_df) + 1, fill = NA, color = "black") +
      # geom_text(aes(colour = "black", label = Description), size = 3.5, y = -0.2, vjust = 0.4, hjust = 1) +
      # scale_colour_manual(values = "black", guide = F) + 
                theme(axis.title = element_text(face = "plain"),
                      axis.text = element_text(face = "plain"),
                      # panel.grid.minor = element_blank(),
                      # panel.grid.major = element_blank(),
                      legend.title = element_text(size=as.numeric(16), face = "plain"),
                      legend.position="right",
                      legend.text = element_text(size = as.numeric(15)),
                      axis.text.x=element_text(size=as.numeric(15),face = "plain"),
                      axis.text.y=element_text(size=as.numeric(15),face = "plain"),
                      axis.title.x=element_text(size=as.numeric(16),face = "plain"),
                      axis.title.y=element_text(size=as.numeric(16))) +
                scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) + 
                theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 12, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 12, 
           height = 10)
    write_enrichment_bubble_json(dotp$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.json", sep="_")))
    ### if ont = "ALL"
    go_enrich_plot = rbind(head(go_enrich_results@result[go_enrich_results@result$ONTOLOGY == "biological_process",],10),
                           head(go_enrich_results@result[go_enrich_results@result$ONTOLOGY == "molecular_function",],10),
                           head(go_enrich_results@result[go_enrich_results@result$ONTOLOGY == "cellular_component",],10))
    go_enrich_plot$GeneRatio <- mixedToFloat(go_enrich_plot$GeneRatio)
    go_enrich_plot = go_enrich_plot[order(go_enrich_plot$GeneRatio,decreasing = TRUE),]
    go_enrich_plot$Description = factor(go_enrich_plot$Description, levels = rev(go_enrich_plot$Description))
    go_enrich_plot$ONTOLOGY = factor(go_enrich_plot$ONTOLOGY, levels = c("biological_process","molecular_function","cellular_component"))
    dotp <- ggplot(go_enrich_plot) +
                geom_point(aes(x=GeneRatio, y=Description, color=p.adjust, size = Count), shape = 19)+
                theme_bw() +
                ylab("") +
                scale_color_gradient(high = "blue", low = "red") +
      # geom_rect(ymin = 0, ymax = 12, xmin = 0, xmax = nrow(go_enrich_df) + 1, fill = NA, color = "black") +
      # geom_text(aes(colour = "black", label = Description), size = 3.5, y = -0.2, vjust = 0.4, hjust = 1) +
      # scale_colour_manual(values = "black", guide = F) + 
                theme(axis.title = element_text(face = "plain"),
                      axis.text = element_text(face = "plain"),
                      # panel.grid.minor = element_blank(),
                      # panel.grid.major = element_blank(),
                      legend.title = element_text(size=as.numeric(16), face = "plain"),
                      legend.position="right",
                      legend.text = element_text(size = as.numeric(15)),
                      axis.text.x=element_text(size=as.numeric(15),face = "plain"),
                      axis.text.y=element_text(size=as.numeric(15),face = "plain"),
                      axis.title.x=element_text(size=as.numeric(16),face = "plain"),
                      axis.title.y=element_text(size=as.numeric(16))) +
                scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) + 
                theme(plot.margin = unit(c(1,1,1,1),'lines')) +
                facet_grid(ONTOLOGY~., scale = "free")
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.pdf", sep="_")), 
           device = "pdf",
           width = 12, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.png", sep="_")), 
           device = "png",
           width = 12, 
           height = 10)
    write_faceted_enrichment_bubble_json(dotp$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.json", sep="_")), c("biological_process", "molecular_function", "cellular_component"))
    
    ## emapplot: Enrichment Map
    # pt <- pairwise_termsim(go_enrich_results)
    
    # ### Use trycatch funtion to automatically decrease the value of the parameter showCategory 
    # ###   in order to fix the issue of too few enriched items
    # pass <- TRUE
    # cat_num <- 30 # default value
    # text = "无法生成emapplot富集图！"
    # emapp <- ggplot() +
    #   annotate("text", x = 4, y = 25, size=8, label = text) +
    #   theme_void() +
    #   theme(panel.background = element_rect(fill = 'white'))

    # while (pass & cat_num>0) {
    #   tryCatch({
    #     emapp <- emapplot(pt, node_label = 'category', 
    #                       showCategory = cat_num)  # 'category', 'group', 'all' and 'none'.
    #     pass <- FALSE
    #   }, error = function(e) {
    #     assign("cat_num", cat_num - 5, env=globalenv())
    #   })
    # }
    
    # ggsave(emapp, 
    #        filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.pdf", sep="_")), 
    #        device = "pdf",
    #        width = 13, 
    #        height = 10)
    # ggsave(emapp, 
    #        filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.png", sep="_")), 
    #        device = "png",
    #        width = 13, 
    #        height = 10)
  }
}
