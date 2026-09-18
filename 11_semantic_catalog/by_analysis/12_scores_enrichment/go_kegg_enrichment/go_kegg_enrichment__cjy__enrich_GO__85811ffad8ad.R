log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(jsonlite)

orgdb <- "org.Rn.eg.db"
library(org.Rn.eg.db)

library(DOSE)

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

write_empty_json <- function(path) {
  write(toJSON(list(width=800, height=600, dataset=list(), series=list()), auto_unbox=TRUE, pretty=TRUE), path)
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

write_enrichment_bar_json <- function(plot_df, path) {
  bar_df <- plot_df %>% arrange(desc(Description))
  if (nrow(bar_df) == 0) {
    write_empty_json(path)
    return()
  }
  bar_json <- list(
    width=800, height=700, margin=list(top=30, right=100, bottom=60, left=280), color=list("#B71B1B", "#4682B4"),
    dataset=list(list(source=list(
      pathway=as.list(as.character(bar_df$Description)),
      count=as.list(bar_df$Count),
      pAdjust=as.list(bar_df$p.adjust)
    ))),
    xAxis=list(list(type="scaleLinear", position="bottom", show=TRUE, name="Count", axisLine=list(show=TRUE, lineStyle=list(width=1)), axisTick=list(show=TRUE, lineStyle=list(width=1)), axisLabel=list(show=TRUE, fontSize="11"), nameTextStyle=list(fontSize="14", fontWeight="bold"))),
    yAxis=list(list(type="scaleBand", position="left", show=TRUE, axisLine=list(show=TRUE, lineStyle=list(width=1)), axisTick=list(show=TRUE, lineStyle=list(width=1)), axisLabel=list(show=TRUE, fontSize="11"))),
    series=list(list(type="pathwayEnrichmentBar", datasetIndex=0, encode=list(pathway="pathway", count="count", pAdjust="pAdjust"))),
    legend=list(list(show=TRUE, title="p.adjust", right=30, top=50, width=20, height=150, nameTextStyle=list(fontSize="12"), textStyle=list(fontSize="10"), useScientificNotation=use_scientific_notation(bar_df$p.adjust)))
  )
  write(toJSON(bar_json, auto_unbox=TRUE, pretty=TRUE), path)
}

write_faceted_enrichment_bar_json <- function(plot_df, path, ontology_levels = NULL) {
  facet_bar_df <- with_ontology_levels(plot_df, ontology_levels)
  if (nrow(facet_bar_df) == 0) {
    write_empty_json(path)
    return()
  }
  facet_bar_df <- facet_bar_df %>% group_by(ONTOLOGY) %>% arrange(desc(Description), .by_group=TRUE) %>% ungroup()
  facet_bar_json <- list(
    width=800, height=500, margin=list(top=40, right=80, bottom=60, left=240), color=list("#B71B1B", "#4682B4"),
    dataset=list(list(source=list(
      pathway=as.list(as.character(facet_bar_df$Description)),
      count=as.list(facet_bar_df$Count),
      pAdjust=as.list(facet_bar_df$p.adjust),
      ontology=as.list(as.character(facet_bar_df$ONTOLOGY))
    ))),
    xAxis=list(list(type="scaleLinear", show=TRUE, axisGrid=list(show=TRUE), name="Count")),
    yAxis=list(list(type="scaleBand", show=TRUE, axisGrid=list(show=TRUE))),
    series=list(list(type="bar", barType="facetedEnrichment", datasetIndex=0, encode=list(y="pathway", x="count", color="pAdjust"), facet=list(field="ontology", direction="column", wrap=1), label=list(show=FALSE), border=list(show=TRUE))),
    legend=list(list(show=TRUE, useScientificNotation=use_scientific_notation(facet_bar_df$p.adjust)))
  )
  write(toJSON(facet_bar_json, auto_unbox=TRUE, pretty=TRUE), path)
}

write_enrichment_bubble_json <- function(plot_df, path) {
  bubble_df <- plot_df %>% arrange(desc(Description))
  if (nrow(bubble_df) == 0) {
    write_empty_json(path)
    return()
  }
  count_values <- suppressWarnings(as.numeric(as.character(bubble_df$Count)))
  gene_ratio_values <- suppressWarnings(as.numeric(as.character(bubble_df$GeneRatio)))
  size_values <- as.list(round(quantile(count_values, probs=c(0.25, 0.5, 0.75, 1.0), na.rm=TRUE)))
  bubble_json <- list(
    width=1020, height=680, outerBorder=TRUE, margin=list(top=30, right=220, bottom=70, left=340), color=list("#B71B1B", "#4682B4"),
    dataset=list(list(source=list(
      pathway=as.list(as.character(bubble_df$Description)),
      geneRatio=as.list(round(gene_ratio_values, 6)),
      count=as.list(count_values),
      pAdjust=as.list(bubble_df$p.adjust)
    ))),
    series=list(list(type="pathwayEnrichmentBubble", datasetIndex=0, encode=list(pathway="pathway", geneRatio="geneRatio", count="count", pAdjust="pAdjust"), border=list(show=TRUE), label=list(show=FALSE))),
    xAxis=list(list(type="scaleLinear", position="bottom", show=TRUE, name="GeneRatio", axisLine=list(show=TRUE, lineStyle=list(width=1)), axisTick=list(show=TRUE, lineStyle=list(width=1)), axisLabel=list(show=TRUE, color="#000", fontFamily="sans-serif", fontSize="11"), nameTextStyle=list(color="#000", fontFamily="sans-serif", fontSize="14", fontWeight="bold"))),
    yAxis=list(list(type="scaleBand", position="left", show=TRUE, axisLine=list(show=TRUE, lineStyle=list(width=1)), axisTick=list(show=TRUE, lineStyle=list(width=1)), axisLabel=list(show=TRUE, color="#000", fontFamily="sans-serif", fontSize="11"))),
    legend=list(), sizeLegend=list(show=TRUE, title="Count", orient="vertical", count=4, value=size_values, left=820, top=220, bordered=TRUE, color="#ffffff", fontSize=12, site="right"),
    colorbar=list(show=TRUE, title="p.adjust", width=20, height=160, site="right", left=820, top=40, fontSize=12, tickCount=5, useScientificNotation=use_scientific_notation(bubble_df$p.adjust))
  )
  write(toJSON(bubble_json, auto_unbox=TRUE, pretty=TRUE), path)
}

write_faceted_enrichment_bubble_json <- function(plot_df, path, ontology_levels = NULL) {
  facet_bubble_df <- with_ontology_levels(plot_df, ontology_levels)
  if (nrow(facet_bubble_df) == 0) {
    write_empty_json(path)
    return()
  }
  facet_bubble_df <- facet_bubble_df %>% group_by(ONTOLOGY) %>% arrange(desc(Description), .by_group=TRUE) %>% ungroup()
  count_values <- suppressWarnings(as.numeric(as.character(facet_bubble_df$Count)))
  gene_ratio_values <- suppressWarnings(as.numeric(as.character(facet_bubble_df$GeneRatio)))
  max_gene_ratio <- max(gene_ratio_values, na.rm=TRUE)
  if (!is.finite(max_gene_ratio)) {
    max_gene_ratio <- 0
  }
  p_adjust_values <- suppressWarnings(as.numeric(as.character(facet_bubble_df$p.adjust)))
  p_adjust_range <- range(p_adjust_values, na.rm=TRUE)
  if (!all(is.finite(p_adjust_range))) {
    p_adjust_range <- c(0, 1)
  }
  size_values <- as.list(round(quantile(count_values, probs=c(0.33, 1.0), na.rm=TRUE)))
  facet_bubble_json <- list(
    width=920, height=680, outerBorder=TRUE, margin=list(top=30, right=140, bottom=70, left=310), color=list("#e76b66", "#5e9bd3"),
    dataset=list(list(source=list(
      description=as.list(as.character(facet_bubble_df$Description)),
      geneRatio=as.list(round(gene_ratio_values, 6)),
      pAdjust=as.list(facet_bubble_df$p.adjust),
      count=as.list(count_values),
      ontology=as.list(as.character(facet_bubble_df$ONTOLOGY))
    ))),
    series=list(list(type="facetedEnrichmentBubble", datasetIndex=0, encode=list(x="geneRatio", y="description", size="count", color="pAdjust", tooltip=list("description", "geneRatio", "pAdjust", "count", "ontology")), facet=list(field="ontology", direction="column", wrap=1), sort=list(enabled=TRUE, field="geneRatio", order="desc"), border=list(show=TRUE), label=list(show=FALSE), facetTitleStyle=list(fontSize="12", fontWeight="bold"))),
    xAxis=list(list(type="scaleLinear", position="bottom", show=TRUE, name="GeneRatio", min=0, max=round(max_gene_ratio * 1.1, 4), tickCount=5, axisLine=list(show=TRUE, lineStyle=list(width=1)), axisTick=list(show=TRUE, lineStyle=list(width=1)), axisGrid=list(show=TRUE), axisLabel=list(show=TRUE, color="#000", fontFamily="sans-serif", fontSize="11"), nameTextStyle=list(color="#000", fontFamily="sans-serif", fontSize="14", fontWeight="bold"))),
    yAxis=list(list(type="scaleBand", position="left", show=TRUE, name="", axisLine=list(show=TRUE, lineStyle=list(width=1)), axisTick=list(show=TRUE, lineStyle=list(width=1)), axisGrid=list(show=TRUE), axisLabel=list(show=TRUE, color="#000", fontFamily="sans-serif", fontSize="11"))),
    legend=list(), sizeLegend=list(show=TRUE, title="Count", orient="vertical", value=size_values, right=20, top=250, bordered=TRUE, color="#ffffff", fontSize=12, site="right"),
    colorbar=list(show=TRUE, title="p.adjust", width=18, height=90, site="right", right=20, top=100, fontSize=12, tickCount=3, colors=list("#e76b66", "#5e9bd3"), reverse=FALSE, domain=as.list(p_adjust_range), useScientificNotation=use_scientific_notation(facet_bubble_df$p.adjust))
  )
  write(toJSON(facet_bubble_json, auto_unbox=TRUE, pretty=TRUE), path)
}

# Extract genes from DEG
need_DEG <- DEG[,c("log2FoldChange","pvalue","padj","regulate")]
gene_up <- rownames(need_DEG[with(need_DEG, regulate == "Up"),])
gene_down <- rownames(need_DEG[with(need_DEG, regulate == "Down"),])

# Convert gene ENSEMBL ID to ENTREZ ID
gene_up_entrez <- unique(as.character(na.omit(bitr(gene_up,
                                                   fromType="ENSEMBL",
                                                   toType="ENTREZID",
                                                   OrgDb=orgdb)[,2]))) #"org.Hs.eg.db" "org.Mm.eg.db"
gene_down_entrez <- unique(as.character(na.omit(bitr(gene_down,
                                                     fromType="ENSEMBL",
                                                     toType="ENTREZID",
                                                     OrgDb=orgdb)[,2]))) #"org.Hs.eg.db" "org.Mm.eg.db"
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
  
  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
  if (length(gene_list) < 2) {
    
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, set, paste(contrast, set, "GO_enrich_results.xls", sep="_"))
    headers <- glue("ONTOLOGY\tID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\n")
    
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
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    write_empty_json(file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.json", sep="_")))
    write_empty_json(file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.json", sep="_")))
    write_empty_json(file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.json", sep="_")))
    write_empty_json(file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.json", sep="_")))
    
  } else {
    # GO enrichment
    go_enrich_results <- enrichGO(gene = gene_list,
                                  OrgDb = orgdb,
                                  ont   = "ALL"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                                  pvalueCutoff  = plot_padj_cutoff,
                                  qvalueCutoff  = plot_q_cutoff,
                                  readable      = TRUE)
    write.table(go_enrich_results@result, 
                file = file.path(outdir, set, paste(contrast, set, "GO_enrich_results.xls", sep="_")),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(go_enrich_results, 
         file = file.path(outdir, set, paste(contrast, set, "GO_enrich_results.Rdata", sep="_")))
    
    # Plots
    ## barplot
    barp1 <- barplot(go_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp1, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(barp1, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)
    ### if ont = "ALL"
    barp2 <- barplot(go_enrich_results, split = "ONTOLOGY", font.size = 14, showCategory = 10) +
      facet_grid(ONTOLOGY~., scale = "free") +    
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp2, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(barp2, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.png", sep="_")), 
           device = "png",
           width = 13, 
            height = 10)
    
    ## dotplot 
    dotp1 <- enrichplot::dotplot(go_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp1, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(dotp1, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)
    ### if ont = "ALL"
    dotp2 <- enrichplot::dotplot(go_enrich_results, font.size = 14, split = 'ONTOLOGY', showCategory = 10) +
      facet_grid(ONTOLOGY~., scale = "free") + 
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp2, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(dotp2, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## emapplot: Enrichment Map
    pt <- pairwise_termsim(go_enrich_results)
    
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
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(emapp, 
           filename = file.path(outdir, set, paste(contrast, set, "GO_enrich_emapplot.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)

    write_enrichment_bar_json(barp1$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot.json", sep="_")))
    write_faceted_enrichment_bar_json(barp2$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_barplot2.json", sep="_")), c("biological_process", "molecular_function", "cellular_component"))
    write_enrichment_bubble_json(dotp1$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot.json", sep="_")))
    write_faceted_enrichment_bubble_json(dotp2$data, file.path(outdir, set, paste(contrast, set, "GO_enrich_dotplot2.json", sep="_")), c("biological_process", "molecular_function", "cellular_component"))
  }
}
