log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

work_dir <- getwd()

options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
if (snakemake@params[["species"]] == "homo_sapiens") {
  orgdb <- "org.Hs.eg.db"
  spe <- "hsa"
  library(org.Hs.eg.db)
  library(KEGG.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_hsa_pathway_id2description"]], sep="\t", header=TRUE)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  orgdb <- "org.Mm.eg.db"
  spe <- "mmu"
  library(org.Mm.eg.db)
  library(KEGG.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_mmu_pathway_id2description"]], sep="\t", header=TRUE)
}  
library(DOSE)
library(pathview)

# Input
## snv results
snv <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE)
## Cutoff values for enrichment analysis plot
plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
plot_q_cutoff <- as.double(snakemake@params[["enrich_plot_q_cutoff"]])
## Output directory
outdir <- snakemake@params[["outdir"]]
tmpdir <- R.utils::getAbsolutePath(snakemake@params[["tmpdir"]])

# Convert gene SYMBOL ID to ENTREZ ID
gene_entrez <- unique(as.character(na.omit(bitr(snv$Hugo_Symbol,
                                                   fromType="SYMBOL",
                                                   toType="ENTREZID",
                                                   OrgDb=orgdb)[,2]))) #"org.Hs.eg.db" "org.Mm.eg.db"

# Perform KEGG enrichment analysis and plot results
gene_list = gene_entrez

kegg_enrich_results <- enrichKEGG(gene  = gene_list,
                                  organism  = spe, #hsa mmu
                                  pvalueCutoff = plot_padj_cutoff,
                                  qvalueCutoff = plot_q_cutoff,
                                  use_internal_data = T)

  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
if (is.null(kegg_enrich_results)) {
    
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir,"KEGG_enrich_results.xls")
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
           filename = file.path(outdir, "KEGG_enrich_barplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "KEGG_enrich_barplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "KEGG_enrich_dotplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "KEGG_enrich_dotplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, "KEGG_enrich_emapplot.pdf"), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, "KEGG_enrich_emapplot.png"), 
           device = "png",
           width = 4, 
           height = 4)
    
    # generate empty pathview zip file
    pathview_zip <- file.path(outdir, "KEGG_enrich_pathway.zip")
    cmd2 <- glue("echo | zip -q > {pathview_zip} && zip -dq {pathview_zip} -")
    system(cmd2, intern=TRUE)
    
  } else {
    # KEGG enrichment
    kegg_enrich_results <- enrichKEGG(gene  = gene_list,
                                      organism  = spe, #hsa mmu
                                      pvalueCutoff = plot_padj_cutoff,
                                      qvalueCutoff = plot_q_cutoff,
                                      use_internal_data = T)
    
    kegg_enrich_results@result <- merge(kegg_enrich_results@result, pathway_id2description, by.x="ID", by.y="pathway_id") %>% 
      dplyr::select(-Description) %>% 
      dplyr::rename(Description = pathway_description) %>%
      dplyr::arrange(p.adjust)
    #kegg_enrich_results@result <- kegg_enrich_results@result[, c(1,9,2,3,4,5,6,7,8)] 
    kegg_enrich_results@result <- kegg_enrich_results@result %>% relocate(geneID, .after = last_col())  %>% mutate(Count =as.numeric(sub("(\\d+)/\\d+","\\1",GeneRatio)) )
    rownames(kegg_enrich_results@result) <- kegg_enrich_results@result$ID
    
    kk_read <- DOSE::setReadable(kegg_enrich_results, 
                                 OrgDb = orgdb, 
                                 keyType = 'ENTREZID') #ENTREZID to gene Symbol
    write.table(kk_read@result, 
                file = file.path(outdir, "KEGG_enrich_results.xls"),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(kegg_enrich_results, 
         file = file.path(outdir, "KEGG_enrich_results.Rdata"))
    
    # Plots
    ## barplot
    barp <- barplot(kegg_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp, 
           filename = file.path(outdir, "KEGG_enrich_barplot.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, "KEGG_enrich_barplot.png"), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## dotplot 
    dotp <- enrichplot::dotplot(kegg_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp, 
           filename = file.path(outdir, "KEGG_enrich_dotplot.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, "KEGG_enrich_dotplot.png"), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## emapplot: Enrichment Map
    pt <- pairwise_termsim(kegg_enrich_results)
    
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
           filename = file.path(outdir, "KEGG_enrich_emapplot.pdf"), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(emapp, 
           filename = file.path(outdir, "KEGG_enrich_emapplot.png"), 
           device = "png",
           width = 13, 
           height = 10)
    
#    ## pathview
#    genelist <- as.numeric(DEG[, "log2FoldChange"])
#    names(genelist) <- row.names(DEG)
#    pathview_dir <- file.path(outdir, "KEGG_enrich_pathway")
#    system(paste("mkdir", pathview_dir))
#    
#    pathway_list <- kk_read@result$ID
#    problematic_pathways <- c("hsa05416", "hsa05206", "mmu05206", "hsa04723", "mmu04723", "hsa04215", "mmu04215", "hsa01212", "mmu01212", 
#                              "hsa01040", "hsa01250", "hsa01230", "hsa01200", "hsa00534", "hsa00533", "hsa00511", "hsa01232", "hsa00515", "hsa00531")
#    sel = which(pathway_list %in% problematic_pathways)
#    if (length(sel)) {
#      pathway_list <- pathway_list[-sel]
#    }
#    
#    setwd(pathview_dir)
#    if (length(pathway_list) < 20) {
#      select_pathway <- pathway_list[1:length(pathway_list)]
#    } else {
#      select_pathway <- pathway_list[1:20]
#    }
#    
#    pv.out <- pathview(
#      gene.data     = genelist,
#      pathway.id    = select_pathway,
#      species       = spe , # hsa mmu
#      gene.idtype   = "ENSEMBL",
#      kegg.native   = T,
#      new.signature = F,
#      limit         = list(gene=2, cpd=1),
#      bins          = list(gene=20, cpd=10)
#    )
#    
#    # pick the first pathview plot as the example plot in summary report
#    # save it to results/tmp/pathview_example.png if it doesn't exist
#    if (set == "all") {
#      if (file.exists(file.path(tmpdir, "pathview_example.png"))) {
#        system(paste("rm -f ", file.path(tmpdir, "pathview_example.png")), intern=TRUE)
#      }
#      dir.create(tmpdir, showWarnings = FALSE)
#      first_pathway_id <- pathway_list[1]
#      system(paste("cp", paste0(first_pathway_id, ".pathview.png"), file.path(tmpdir, "pathview_example.png")), intern=TRUE)
#    }
#    
#    # compress the folder of pathview plots
#    zip_file <-  "../KEGG_enrich_pathway.zip"
#    zip_dir <- "../KEGG_enrich_pathway"
#    system(paste("zip", zip_file, "./*"), intern=TRUE)
#    system(paste("rm -rf", zip_dir), intern=TRUE)
#    setwd(work_dir)
}
