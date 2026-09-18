if (length(snakemake@log) > 0) {
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")
}

options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(jsonlite)

use_scientific_notation <- function(values) {
  values <- suppressWarnings(as.numeric(as.character(values)))
  finite_values <- values[is.finite(values)]
  length(finite_values) > 0 && max(finite_values) < 0.001
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

# Input
## DEG results
DEG <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1, check.names=FALSE)
## Cutoff values for enrichment analysis plot
plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
plot_q_cutoff <- as.double(snakemake@params[["enrich_plot_q_cutoff"]])
## Contrast
contrast <- gsub("_diffexp.DEG.xls", "", basename(snakemake@input[[1]]))
## Output directory
outdir <- file.path(snakemake@params[["outdir"]], contrast)

# Extract genes from DEG
need_DEG <- DEG[,c("log2FoldChange","pvalue","padj","regulate")]
gene_up <- rownames(need_DEG[with(need_DEG, regulate == "Up"),])
gene_down <- rownames(need_DEG[with(need_DEG, regulate == "Down"),])

# Convert gene ENSEMBL ID to ENTREZ ID
gene_up_entrez <- unique(as.character(na.omit(bitr(gene_up,
                                                   fromType="ENSEMBL",
                                                   toType="ENTREZID",
                                                   OrgDb="org.Hs.eg.db")[,2])))
gene_down_entrez <-  unique(as.character(na.omit(bitr(gene_down,
                                                      fromType="ENSEMBL",
                                                      toType="ENTREZID",
                                                      OrgDb="org.Hs.eg.db")[,2])))
gene_all_entrez <- unique(c(gene_up_entrez, gene_down_entrez ))

# Perform DO enrichment analysis and plot results

sets = list("all", "up", "down")

for (set in sets) {
  
  if (set == "all") {
    gene_list = gene_all_entrez
  } else if (set == "up") {
    gene_list = gene_up_entrez
  } else if (set == "down") {
    gene_list = gene_down_entrez
  }

  do_enrich_results <- enrichDO(gene  = gene_list,
                                  pvalueCutoff = plot_padj_cutoff,
                                  qvalueCutoff = plot_q_cutoff,
                                  readable = TRUE)

  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
  if (is.null(do_enrich_results)) {
    
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, set, paste(contrast, set, "DO_enrich_results.xls", sep="_"))
    headers <- glue("ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\n")
    
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
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_emapplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_emapplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    
  } else {
    # DO enrichment
    do_enrich_results <- enrichDO(gene  = gene_list,
                                  pvalueCutoff = plot_padj_cutoff,
                                  qvalueCutoff = plot_q_cutoff,
                                  readable = TRUE)
    write.table(do_enrich_results@result, 
                file = file.path(outdir, set, paste(contrast, set, "DO_enrich_results.xls", sep="_")),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(do_enrich_results, 
         file = file.path(outdir, set, paste(contrast, set, "DO_enrich_results.Rdata", sep="_")))
    
    # Plots
    ## barplot
    barp <- barplot(do_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp,
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_barplot.pdf", sep="_")),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(barp,
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_barplot.png", sep="_")),
           device = "png",
           width = 13,
           height = 10)
    write_enrichment_bar_json(barp$data, file.path(outdir, set, paste(contrast, set, "DO_enrich_barplot.json", sep="_")))
    
    
    ## dotplot
    dotp <- enrichplot::dotplot(do_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp,
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_dotplot.pdf", sep="_")),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(dotp,
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_dotplot.png", sep="_")),
           device = "png",
           width = 13,
           height = 10)
    write_enrichment_bubble_json(dotp$data, file.path(outdir, set, paste(contrast, set, "DO_enrich_dotplot.json", sep="_")))
    
    ## emapplot: Enrichment Map
    pt <- pairwise_termsim(do_enrich_results)
    
    ### Use trycatch funtion to automatically decrease the value of the parameter showCategory 
    ###   in order to fix the issue of too few enriched items
    pass <- TRUE
    cat_num <- 30 # default value
    text = "无法生成emapplot富集图！"
    emapp <- ggplot() +
      annotate("text", x = 4, y = 25, size=8, label = text) +
      theme_void() +
      theme(panel.background = element_rect(fill = 'white'))

    while (pass & cat_num>0) {
      tryCatch({
        emapp <- emapplot(pt, node_label = 'category', 
                          showCategory = cat_num)  # 'category', 'group', 'all' and 'none'.
        pass <- FALSE
      }, error = function(e) {
        assign("cat_num", cat_num - 5, env=globalenv())
      })
    }
    
    ggsave(emapp,
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_emapplot.pdf", sep="_")),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(emapp,
           filename = file.path(outdir, set, paste(contrast, set, "DO_enrich_emapplot.png", sep="_")),
           device = "png",
           width = 13,
           height = 10)
  }
}
