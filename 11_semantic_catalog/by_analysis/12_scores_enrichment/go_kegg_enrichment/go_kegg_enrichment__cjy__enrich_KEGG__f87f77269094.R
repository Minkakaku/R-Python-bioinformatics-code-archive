log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

work_dir <- getwd()

options(pathview.download.method = 'auto')
options(clusterProfiler.download.method = 'auto')
options(timeout = 1200)
options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(jsonlite)

reset_dir <- function(path) {
  unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

kegg_pathway_dir <- normalizePath(snakemake@params[["kegg_pathway_dir"]], mustWork = FALSE)

#if ("package:KEGG.db" %in% search()) {
#  detach("package:KEGG.db",unload = TRUE)
#}
#if ("KEGG.db" %in% installed.packages()) {
#  remove.packages("KEGG.db")
#}

contrast <- gsub("_diffexp.DEG.xls", "", basename(snakemake@input[[1]]))
user_lib <- snakemake@params[["kegg_r_lib"]]
.libPaths(c(user_lib, .libPaths()))
library(KEGG.db, lib.loc=user_lib)

if (snakemake@params[["species"]] == "homo_sapiens") {
  orgdb <- "org.Hs.eg.db"
  spe <- "hsa"
  library(org.Hs.eg.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_hsa_pathway_id2description"]], sep="\t", header=TRUE)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  orgdb <- "org.Mm.eg.db"
  spe <- "mmu"
  library(org.Mm.eg.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_mmu_pathway_id2description"]], sep="\t", header=TRUE)
}  
library(DOSE)
library(pathview)

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
tmpdir <- R.utils::getAbsolutePath(snakemake@params[["tmpdir"]])

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
gene_all_entrez <- unique(c(gene_up_entrez, gene_down_entrez))

# Perform KEGG enrichment analysis and plot results

sets = list("all", "up", "down")

run_pathview <- function(genelist, pathway_list, spe, n1, n2){
  while (n2 <= 20 && n1 <= length(pathway_list)) {
    result <- tryCatch({
      pathview(
        gene.data     = genelist,
        pathway.id    = pathway_list[n1],
        species       = spe, # hsa mmu
        gene.idtype   = "ENSEMBL",
        kegg.native   = T,
        new.signature = F,
        kegg.dir = kegg_pathway_dir,
        limit         = list(gene=2, cpd=1),
        bins          = list(gene=20, cpd=10))
      TRUE
    }, error = function(e){
      print("Handling with error:")
      print(str(e))
      FALSE
    }, warning = function(w){
      print("Handling with warning:")
      print(str(w))
      FALSE
    })
    n1 <- n1 + 1
    if (isTRUE(result)) {
      n2 <- n2 + 1
    }
  }
  return(n2)
}

for (set in sets) {
  
  if (set == "all") {
    gene_list = gene_all_entrez
  } else if (set == "up") {
    gene_list = gene_up_entrez
  } else if (set == "down") {
    gene_list = gene_down_entrez
  }

  kegg_enrich_results <- enrichKEGG(gene  = gene_list,
                                      organism  = spe, #hsa mmu
                                      pvalueCutoff = plot_padj_cutoff,
                                      qvalueCutoff = plot_q_cutoff,
                                      use_internal_data = T)

  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
  if (is.null(kegg_enrich_results) || nrow(kegg_enrich_results@result) == 0) {
    
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, set, paste(contrast, set, "KEGG_enrich_results.xls", sep="_"))
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
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.pdf", sep="_")), 
           device = "pdf",
           width = 10, 
           height = 10)
    ggsave(p, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.png", sep="_")), 
           device = "png",
           width = 4, 
           height = 4)
    
    # generate empty pathview zip file
    pathview_zip <- file.path(outdir, set, paste(contrast, set, "KEGG_enrich_pathway.zip", sep="_"))
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
                file = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_results.xls", sep="_")),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(kegg_enrich_results, 
         file = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_results.Rdata", sep="_")))
    
    # Plots
    ## barplot
    barp <- barplot(kegg_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + 
      theme(plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## dotplot 
    dotp <- enrichplot::dotplot(kegg_enrich_results, font.size = 14, showCategory = 30) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
      theme(legend.key.size = unit(10, "pt"),
            plot.margin = unit(c(1,1,1,1),'lines'))
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)
    
    ## emapplot: Enrichment Map
    pt <- pairwise_termsim(kegg_enrich_results)
    
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
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.pdf", sep="_")), 
           device = "pdf",
           width = 13, 
           height = 10)
    ggsave(emapp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.png", sep="_")), 
           device = "png",
           width = 13, 
           height = 10)

    # --- JSON GENERATION (FROM PLOT DATA) ---
    # 1. Barplot JSON
    bar_df <- barp$data %>% arrange(desc(Description))
    if(nrow(bar_df) > 0) {
      bar_json <- list(
        width = 800, height = 700, margin = list(top = 30, right = 100, bottom = 60, left = 280), color = list("#B71B1B", "#4682B4"),
        dataset = list(list(source = list(pathway=as.list(as.character(bar_df$Description)), count=as.list(bar_df$Count), pAdjust=as.list(bar_df$p.adjust)))),
        xAxis = list(list(type = "scaleLinear", position = "bottom", show = TRUE, name = "Count", axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, fontSize = "11"), nameTextStyle = list(fontSize = "14", fontWeight = "bold"))),
        yAxis = list(list(type = "scaleBand", position = "left", show = TRUE, axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, fontSize = "11"))),
        series = list(list(type = "pathwayEnrichmentBar", datasetIndex = 0, encode = list(pathway = "pathway", count = "count", pAdjust = "pAdjust"))),
        legend = list(list(show = TRUE, title = "p.adjust", right = 30, top = 50, width = 20, height = 150, nameTextStyle = list(fontSize = "12"), textStyle = list(fontSize = "10"), useScientificNotation = (max(bar_df$p.adjust) < 0.001)))
      )
      write(toJSON(bar_json, auto_unbox=T, pretty=T), file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.json", sep="_")))
    }

    # 2. Dotplot JSON
    bubble_df <- dotp$data %>% arrange(desc(Description))
    if(nrow(bubble_df) > 0) {
      size_legend_values <- as.list(round(quantile(bubble_df$Count, probs = c(0.25, 0.5, 0.75, 1.0))))
      bubble_json <- list(
        width = 1020, height = 680, outerBorder = TRUE, margin = list(top = 30, right = 220, bottom = 70, left = 340), color = list("#B71B1B", "#4682B4"),
        dataset = list(list(source = list(pathway=as.list(as.character(bubble_df$Description)), geneRatio=as.list(round(bubble_df$GeneRatio, 6)), count=as.list(bubble_df$Count), pAdjust=as.list(bubble_df$p.adjust)))),
        series = list(list(type = "pathwayEnrichmentBubble", datasetIndex = 0, encode = list(pathway = "pathway", geneRatio = "geneRatio", count = "count", pAdjust = "pAdjust"), border = list(show = TRUE), label = list(show = FALSE))),
        xAxis = list(list(type = "scaleLinear", position = "bottom", show = TRUE, name = "GeneRatio", axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"), nameTextStyle = list(color = "#000", fontFamily = "sans-serif", fontSize = "14", fontWeight = "bold"))),
        yAxis = list(list(type = "scaleBand", position = "left", show = TRUE, axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"))),
        legend = list(), sizeLegend = list(show = TRUE, title = "Count", orient = "vertical", count = 4, value = size_legend_values, left = 820, top = 220, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"),
        colorbar = list(show = TRUE, title = "p.adjust", width = 20, height = 160, site = "right", left = 820, top = 40, fontSize = 12, tickCount = 5, useScientificNotation = (max(bubble_df$p.adjust) < 0.001))
      )
      write(toJSON(bubble_json, auto_unbox=T, pretty=T), file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.json", sep="_")))
    }
    
    ## pathview
    genelist <- as.numeric(DEG[, "log2FoldChange"])
    names(genelist) <- row.names(DEG)
    pathview_dir <- file.path(outdir, set, paste(contrast, set, "KEGG_enrich_pathway", sep="_"))
    reset_dir(pathview_dir)
    
    pathway_list <- kk_read@result$ID
    setwd(pathview_dir)
    run_pathview(genelist, pathway_list, spe, n1=1, n2=1)

    # pick the first pathview plot as the example plot in summary report
    # save it to results/tmp/pathview_example.png if it doesn't exist
    if (set == "all") {
      if (file.exists(file.path(tmpdir, "pathview_example.png"))) {
        system(paste("rm -f ", file.path(tmpdir, "pathview_example.png")), intern=TRUE)
      }
      dir.create(tmpdir, showWarnings = FALSE)
      first_pathway_id <- pathway_list[1]
      system(paste("cp", paste0(first_pathway_id, ".pathview.png"), file.path(tmpdir, "pathview_example.png")), intern=TRUE)
    }
    
    # compress the folder of pathview plots
    zip_file <-  paste0("../", paste(contrast, set, "KEGG_enrich_pathway", sep="_"), ".zip")
    zip_dir <- paste0("../", paste(contrast, set, "KEGG_enrich_pathway", sep="_"))
    system(paste("zip", zip_file, "./*"), intern=TRUE)
    system(paste("rm -rf", zip_dir), intern=TRUE)
    setwd(work_dir)
  }
}
