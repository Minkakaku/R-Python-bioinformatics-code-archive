log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
if (snakemake@params[["species"]] == "homo_sapiens") {
  orgdb <- "org.Hs.eg.db"
  library(org.Hs.eg.db)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  orgdb <- "org.Mm.eg.db"
  library(org.Mm.eg.db)
}  
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
                                                   OrgDb=orgdb)[,2]))) #"org.Hs.eg.db" "org.Mm.eg.db"

# Perform GO enrichment analysis and plot results
gene_list = gene_entrez
  
# if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
if (length(gene_list) < 2) {
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, "GO_enrich_results.xls")
    headers <- glue("ONTOLOGY\tID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\n")
    
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
           filename = file.path(outdir, "GO_enrich_barplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_barplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_barplot2.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_barplot2.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_dotplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_dotplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_dotplot2.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_dotplot2.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_emapplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "GO_enrich_emapplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    
  } else {
    # GO enrichment
    go_enrich_results <- enrichGO(gene = gene_list,
                                  OrgDb = orgdb,
                                  ont   = "ALL"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                                  pvalueCutoff  = plot_padj_cutoff,
                                  qvalueCutoff  = plot_q_cutoff,
                                  readable      = TRUE)
    write.table(go_enrich_results@result, 
                file = file.path(outdir, "GO_enrich_results.xls"),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(go_enrich_results, 
         file = file.path(outdir, "GO_enrich_results.Rdata"))
    
    # Plots
    ## barplot
    barp <- barplot(go_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp, 
           filename = file.path(outdir, "GO_enrich_barplot.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, "GO_enrich_barplot.png"), 
           device = "png",
           width = 13, 
           height = 10)
    ### if ont = "ALL"
    barp <- barplot(go_enrich_results, split = "ONTOLOGY", font.size = 14, showCategory = 10) +
      facet_grid(ONTOLOGY~., scale = "free") +    
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp, 
           filename = file.path(outdir, "GO_enrich_barplot2.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, "GO_enrich_barplot2.png"), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## dotplot 
    dotp <- enrichplot::dotplot(go_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp, 
           filename = file.path(outdir, "GO_enrich_dotplot.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, "GO_enrich_dotplot.png"), 
           device = "png",
           width = 13, 
           height = 10)
    ### if ont = "ALL"
    dotp <- enrichplot::dotplot(go_enrich_results, font.size = 14, split = 'ONTOLOGY', showCategory = 10) +
      facet_grid(ONTOLOGY~., scale = "free") + 
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp, 
           filename = file.path(outdir, "GO_enrich_dotplot2.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, "GO_enrich_dotplot2.png"), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## emapplot: Enrichment Map
    pt <- pairwise_termsim(go_enrich_results)
    
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
           filename = file.path(outdir, "GO_enrich_emapplot.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(emapp, 
           filename = file.path(outdir, "GO_enrich_emapplot.png"), 
           device = "png",
           width = 13, 
           height = 10)
}


