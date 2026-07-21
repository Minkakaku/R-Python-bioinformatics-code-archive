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
library(stringr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(data.table)
# if (snakemake@params[["species"]] == "homo_sapiens") {
#   orgdb <- "org.Hs.eg.db"
#   spe <- "hsa"
#   library(org.Hs.eg.db)
#   library(KEGG.db)
#   pathway_id2description <- read.table(snakemake@params[["kegg_hsa_pathway_id2description"]], sep="\t", header=TRUE)
# } else if (snakemake@params[["species"]] == "mus_musculus") {
#   orgdb <- "org.Mm.eg.db"
#   spe <- "mmu"
#   library(org.Mm.eg.db)
#   library(KEGG.db)
#   library(KEGG.mmu.db)
#   pathway_id2description <- read.table(snakemake@params[["kegg_mmu_pathway_id2description"]], sep="\t", header=TRUE)
# }  
library(DOSE)
library(pathview)
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

kegg_relation <- as.data.frame(fread(snakemake@params[["kegg_pathway_id2description"]]))

## Cutoff values for enrichment analysis plot
plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
plot_q_cutoff <- as.double(snakemake@params[["enrich_plot_q_cutoff"]])
## Contrast
contrast <- gsub("_diffexp.DEG.xls", "", basename(snakemake@input[[1]]))
## Output directory
outdir <- file.path(snakemake@params[["outdir"]], contrast)
tmpdir <- R.utils::getAbsolutePath(snakemake@params[["tmpdir"]])
kegg_pathway_dir <- to_abs(snakemake@params[["kegg_pathway_dir"]])
dir.create(kegg_pathway_dir, recursive = TRUE, showWarnings = FALSE)



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


# Extract genes from DEG
need_DEG <- DEG[,c("log2FoldChange","pvalue","padj","regulate")]
gene_up <- rownames(need_DEG[with(need_DEG, regulate == "Up"),])
gene_down <- rownames(need_DEG[with(need_DEG, regulate == "Down"),])

# Convert gene ENSEMBL ID to ENTREZ ID
gene_up_entrez <- unique(gene_up)
gene_down_entrez <- unique(gene_down)
gene_all_entrez <- unique(c(gene_up_entrez ,gene_down_entrez))

# Perform KEGG enrichment analysis and plot results

sets = list("all", "up", "down")

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

for (set in sets) {
  
  if (set == "all") {
    gene_list = gene_all_entrez
  } else if (set == "up") {
    gene_list = gene_up_entrez
  } else if (set == "down") {
    gene_list = gene_down_entrez
  }

  kegg_enrich_results <- enricher(gene=gene_list,
                                  pvalueCutoff = plot_padj_cutoff,
                                  pAdjustMethod = "BH",
                                  TERM2GENE = term2gene)

  # if there are too few DEGs as input to enrichment analysis, just generate some empty tables and plots
  if (is.null(kegg_enrich_results)) {
    
    library(glue)
    library(ggplot2)
    library(showtext)
    showtext_auto()
    
    # generate empty results table
    results_xls <- file.path(outdir, set, paste(contrast, set, "KEGG_enrich_results.xls", sep="_"))
    headers <- glue("ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tGeneSymbol\tCount\n")
    
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
    # ggsave(p, 
    #        filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.pdf", sep="_")), 
    #        device = "pdf",
    #        width = 10, 
    #        height = 10)
    # ggsave(p, 
    #        filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.png", sep="_")), 
    #        device = "png",
    #        width = 4, 
    #        height = 4)
    
    # generate empty pathview zip file
    pathview_zip <- file.path(outdir, set, paste(contrast, set, "KEGG_enrich_pathway.zip", sep="_"))
    cmd2 <- glue("echo | zip -q > {pathview_zip} && zip -dq {pathview_zip} -")
    system(cmd2, intern=TRUE)
    
  } else {
    # KEGG enrichment
    kegg_enrich_results <- enricher(gene=gene_list,
                                    pvalueCutoff = plot_padj_cutoff,
                                    pAdjustMethod = "BH",
                                    TERM2GENE = term2gene)
    kegg_enrich_results@result <- merge(kegg_enrich_results@result,kegg_relation,by = "ID") %>% select(-Description.x) %>% rename(Description="Description.y")
    kegg_enrich_results@result <- kegg_enrich_results@result[order(kegg_enrich_results@result$pvalue),]
    kegg_enrich_results@result$GeneSymbol <- sapply(kegg_enrich_results@result$geneID, map_genes_to_symbols, gene2symbol_dict=gene2symbol_dict)
    kegg_enrich_results@result <- kegg_enrich_results@result %>% relocate(GeneSymbol, .after = geneID)
    # kegg_enrich_results@result <- merge(kegg_enrich_results@result, pathway_id2description, by.x="ID", by.y="pathway_id") %>% 
    #   dplyr::select(-Description) %>% 
    #   dplyr::rename(Description = pathway_description) %>%
    #   dplyr::arrange(p.adjust)
    # #kegg_enrich_results@result <- kegg_enrich_results@result[, c(1,9,2,3,4,5,6,7,8)] 
    # kegg_enrich_results@result <- kegg_enrich_results@result %>% relocate(geneID, .after = last_col())  %>% mutate(Count =as.numeric(sub("(\\d+)/\\d+","\\1",GeneRatio)) )
    # rownames(kegg_enrich_results@result) <- kegg_enrich_results@result$ID
    
    # kk_read <- DOSE::setReadable(kegg_enrich_results, 
    #                              OrgDb = orgdb, 
    #                              keyType = 'ENTREZID') #ENTREZID to gene Symbol
    write.table(kegg_enrich_results@result, 
                file = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_results.xls", sep="_")),
                row.names = FALSE,
                col.names = TRUE,
                sep = "\t") 
    save(kegg_enrich_results, 
         file = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_results.Rdata", sep="_")))
    
    # Plots
    ## barplot
    kegg_enrich_plot = head(kegg_enrich_results@result,30)
    kegg_enrich_plot = kegg_enrich_plot[order(kegg_enrich_plot$Count,decreasing = TRUE),]
    kegg_enrich_plot$Description = factor(kegg_enrich_plot$Description, levels = rev(kegg_enrich_plot$Description))
    barp <- ggplot(kegg_enrich_plot,aes(x=Count,y=Description,fill=p.adjust)) +
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
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.pdf", sep="_")), 
           device = "pdf",
           width = 11, 
           height = 10)
    ggsave(barp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.png", sep="_")), 
           device = "png",
           width = 11, 
           height = 10)
    write_enrichment_bar_json(barp$data, file.path(outdir, set, paste(contrast, set, "KEGG_enrich_barplot.json", sep="_")))
    
    ## dotplot
    kegg_enrich_plot = head(kegg_enrich_results@result,30)
    kegg_enrich_plot$GeneRatio <- mixedToFloat(kegg_enrich_plot$GeneRatio)
    kegg_enrich_plot = kegg_enrich_plot[order(kegg_enrich_plot$GeneRatio,decreasing = TRUE),]
    kegg_enrich_plot$Description = factor(kegg_enrich_plot$Description, levels = rev(kegg_enrich_plot$Description))
    dotp <- ggplot(kegg_enrich_plot) +
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
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.pdf", sep="_")), 
           device = "pdf",
           width = 11, 
           height = 10)
    ggsave(dotp, 
           filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.png", sep="_")), 
           device = "png",
           width = 11, 
           height = 10)
    write_enrichment_bubble_json(dotp$data, file.path(outdir, set, paste(contrast, set, "KEGG_enrich_dotplot.json", sep="_")))
    
    # ## emapplot: Enrichment Map
    # pt <- pairwise_termsim(kegg_enrich_results)
    
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
    #        filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.pdf", sep="_")), 
    #        device = "pdf",
    #        width = 13, 
    #        height = 10)
    # ggsave(emapp, 
    #        filename = file.path(outdir, set, paste(contrast, set, "KEGG_enrich_emapplot.png", sep="_")), 
    #        device = "png",
    #        width = 13, 
    #        height = 10)
    
    ## pathview

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

    log2fc_df = DEG[rownames(DEG) %in% gene_list,c("log2FoldChange","pvalue")]
    log2fc_df$gene = rownames(log2fc_df)
    log2fc_df_koid = merge(gene2keggid,log2fc_df,by = "gene")
    log2fc_df_koid = log2fc_df_koid[order(log2fc_df_koid$gene),c("gene","term","log2FoldChange")]
    genelist = log2fc_df_koid[,"log2FoldChange"]
    names(genelist) <- gsub("^ko:", "", log2fc_df_koid$term)

    pathview_dir <- file.path(outdir, set, paste(contrast, set, "KEGG_enrich_pathway", sep="_"))
    dir.create(pathview_dir, recursive = TRUE, showWarnings = FALSE)
    
    pathway_list <- kegg_enrich_results@result$ID
    setwd(pathview_dir)
    run_pathview(genelist, pathway_list, "ko", n1=1, n2=1)

    # pick the first pathview plot as the example plot in summary report
    # save it to results/tmp/pathview_example.png if it doesn't exist
    if (set == "all") {
      if (file.exists(file.path(tmpdir, "pathview_example.png"))) {
        system(paste("rm -f ", file.path(tmpdir, "pathview_example.png")), intern=TRUE)
      }
      dir.create(tmpdir, showWarnings = FALSE)
      first_pathway_id <- pathway_list[1]
      first_pathway_png <- paste0(first_pathway_id, ".pathview.png")
      if (file.exists(first_pathway_png)) {
        file.copy(first_pathway_png, file.path(tmpdir, "pathview_example.png"), overwrite=TRUE)
      }
    }
    
    # compress the folder of pathview plots
    if (length(list.files(".")) == 0) {
      writeLines("No pathview plots were generated for this gene set.", "README.txt")
    }
    zip_file <-  paste0("../", paste(contrast, set, "KEGG_enrich_pathway", sep="_"), ".zip")
    zip_dir <- paste0("../", paste(contrast, set, "KEGG_enrich_pathway", sep="_"))
    system(paste("zip", zip_file, "./*"), intern=TRUE)
    system(paste("rm -rf", zip_dir), intern=TRUE)
    setwd(work_dir)
  }
}
