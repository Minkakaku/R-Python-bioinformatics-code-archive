log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)

# Input
## snv results
snv <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE)
## Cutoff values for enrichment analysis plot
plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
plot_q_cutoff <- as.double(snakemake@params[["enrich_plot_q_cutoff"]])
## Output directory
outdir <- snakemake@params[["outdir"]]

# Convert gene SYMBOL ID to ENTREZ ID
gene_entrez <- unique(as.character(na.omit(bitr(snv$Hugo_Symbol,
                                                   fromType="SYMBOL",
                                                   toType="ENTREZID",
                                                   OrgDb="org.Hs.eg.db")[,2])))

# Perform DO enrichment analysis and plot results
gene_list = gene_entrez

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
    results_xls <- file.path(outdir, "DO_enrich_results.xls")
    headers <- glue("ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\n")
    
    fileConn<-file(results_xls)
    writeLines(headers, fileConn)
    close(fileConn)
    
    # generate empty plots
    text = "突变基因过少，无法进行富集分析！"
    p <- ggplot() + 
      annotate("text", x = 4, y = 25, size=8, label = text) + 
      theme_void() + 
      theme(panel.background = element_rect(fill = 'white'))
    ggsave(p, 
           filename = file.path(outdir, "DO_enrich_barplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "DO_enrich_barplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "DO_enrich_dotplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "DO_enrich_dotplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "DO_enrich_emapplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "DO_enrich_emapplot.png"), 
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
                file = file.path(outdir, "DO_enrich_results.xls"),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(do_enrich_results, 
         file = file.path(outdir, "DO_enrich_results.Rdata"))
    
    # Plots
    ## barplot
    barp <- barplot(do_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp,
           filename = file.path(outdir, "DO_enrich_barplot.pdf"),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(barp,
           filename = file.path(outdir, "DO_enrich_barplot.png"),
           device = "png",
           width = 13,
           height = 10)
    
    
    ## dotplot
    dotp <- enrichplot::dotplot(do_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp,
           filename = file.path(outdir, "DO_enrich_dotplot.pdf"),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(dotp,
           filename = file.path(outdir, "DO_enrich_dotplot.png"),
           device = "png",
           width = 13,
           height = 10)
    
    ## emapplot: Enrichment Map
    pt <- pairwise_termsim(do_enrich_results)
    
    ### Use trycatch funtion to automatically decrease the value of the parameter showCategory 
    ###   in order to fix the issue of too few enriched items
    pass <- TRUE
    cat_num <- 30 # default value
    while (pass) {
      tryCatch({
        emapp <- emapplot(pt, node_label = 'category', 
                          showCategory = cat_num)  # 'category', 'group', 'all' and 'none'.
        pass <- FALSE
      }, error = function(e) {
        assign("cat_num", cat_num - 5, env=globalenv())
      })
    }
    
    ggsave(emapp,
           filename = file.path(outdir, "DO_enrich_emapplot.pdf"),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(emapp,
           filename = file.path(outdir, "DO_enrich_emapplot.png"),
           device = "png",
           width = 13,
           height = 10)
}
