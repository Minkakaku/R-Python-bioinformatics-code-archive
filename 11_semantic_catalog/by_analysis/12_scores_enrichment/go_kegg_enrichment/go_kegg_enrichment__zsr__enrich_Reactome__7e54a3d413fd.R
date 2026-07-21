log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
if (snakemake@params[["species"]] == "homo_sapiens") {
  orgdb <- "org.Hs.eg.db"
  spe <- "human"
  library(org.Hs.eg.db)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  orgdb <- "org.Mm.eg.db"
  spe <- "mouse"
  library(org.Mm.eg.db)
}
library(DOSE)
library(ReactomePA)

# Input
## DEG results
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
                                                   OrgDb=orgdb)[,2]))) #"org.Hs.eg.db" "org.Mm.eg.db"

# Perform Reactome enrichment analysis and plot results
gene_list = gene_entrez

sets = list("all", "up", "down")

reactome_enrich_results <- enrichPathway(gene  = gene_list,
                                         organism  = spe, #human mouse
                                         pvalueCutoff = plot_padj_cutoff,
                                         qvalueCutoff = plot_q_cutoff,
                                         readable = TRUE)

  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
if (is.null(reactome_enrich_results)) {
    
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, "Reactome_enrich_results.xls")
    headers <- glue("ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\n")
    
    fileConn<-file(results_xls)
    writeLines(headers, fileConn)
    close(fileConn)
    
    # generate empty plots
    text = "突变基因基因过少，无法进行富集分析！"
    p <- ggplot() + 
      annotate("text", x = 4, y = 25, size=8, label = text) + 
      theme_void() + 
      theme(panel.background = element_rect(fill = 'white'))
    ggsave(p, 
           filename = file.path(outdir, "Reactome_enrich_barplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "Reactome_enrich_barplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "Reactome_enrich_dotplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "Reactome_enrich_dotplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "Reactome_enrich_emapplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "Reactome_enrich_emapplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    
  } else {
    # Reactome enrichment
    reactome_enrich_results <- enrichPathway(gene  = gene_list,
                                             organism  = spe, #human mouse
                                             pvalueCutoff = plot_padj_cutoff,
                                             qvalueCutoff = plot_q_cutoff,
                                             readable = TRUE)
    write.table(reactome_enrich_results@result, 
                file = file.path(outdir, "Reactome_enrich_results.xls"),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(reactome_enrich_results, 
         file = file.path(outdir, "Reactome_enrich_results.Rdata"))
    
    # Plots
    ## barplot
    barp <- barplot(reactome_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp,
           filename = file.path(outdir, "Reactome_enrich_barplot.pdf"),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(barp,
           filename = file.path(outdir, "Reactome_enrich_barplot.png"),
           device = "png",
           width = 13,
           height = 10)
    
    ## dotplot 
    dotp <- enrichplot::dotplot(reactome_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp,
           filename = file.path(outdir, "Reactome_enrich_dotplot.pdf"),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(dotp,
           filename = file.path(outdir, "Reactome_enrich_dotplot.png"),
           device = "png",
           width = 13,
           height = 10)
    
    ## emapplot: emapplot: Enrichment Map
    pt <- pairwise_termsim(reactome_enrich_results)
    
    ### Use trycatch funtion to automatically decrease the value of the parameter showCategory 
    ###   in order to fix the issue of too few enriched items
    pass <- TRUE
    cat_num <- 30 # default value
    text = "无法生成emapplot富集图！"
    emapp <- ggplot() + 
      annotate("text", x = 4, y = 25, size=8, label = text) + 
      theme_void() + 
      theme(panel.background = element_rect(fill = 'white'))

    while (pass & cat_num > 0) {
      tryCatch({
        emapp <- emapplot(pt, node_label = 'category', 
                          showCategory = cat_num)  # 'category', 'group', 'all' and 'none'.
        pass <- FALSE
      }, error = function(e) {
        assign("cat_num", cat_num - 5, env=globalenv())
      })
    }
    
    ggsave(emapp,
           filename = file.path(outdir, "Reactome_enrich_emapplot.pdf"),
           device = "pdf",
           width = 13,
           height = 10)
    ggsave(emapp,
           filename = file.path(outdir, "Reactome_enrich_emapplot.png"),
           device = "png",
           width = 13,
           height = 10)
}
