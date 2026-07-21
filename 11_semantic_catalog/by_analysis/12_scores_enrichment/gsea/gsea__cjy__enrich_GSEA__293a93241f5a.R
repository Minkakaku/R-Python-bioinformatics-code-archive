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
library(dplyr)
library(ggplot2)
library(ggridges)
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

contrast <- gsub("_diffexp.xls", "", basename(snakemake@input[[1]]))
user_lib <- snakemake@params[["kegg_r_lib"]]
.libPaths(c(user_lib, .libPaths()))
library(KEGG.db, lib.loc=user_lib)

if (snakemake@params[["species"]] == "homo_sapiens") {
  orgdb <- "org.Hs.eg.db"
  spe1 <- "hsa"
  spe2 <- "human"
  library(org.Hs.eg.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_hsa_pathway_id2description"]], sep="\t", header=TRUE)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  orgdb <- "org.Mm.eg.db"
  spe1 <- "mmu"
  spe2 <- "mouse"
  library(org.Mm.eg.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_mmu_pathway_id2description"]], sep="\t", header=TRUE)
}  
library(pathview)
library(ReactomePA)

# Input
## DESeq2 results
DESeq2 <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE) %>% 
  dplyr::filter(padj != "NA") # remove genes with extreme log2FoldChange value caused by some strange samples
## Cutoff values for enrichment analysis plot
plot_padj_cutoff <- as.double(snakemake@params[["enrich_plot_padj_cutoff"]])
## Contrast
contrast <- gsub("_diffexp.xls", "", basename(snakemake@input[[1]]))
## Output directory
outdir <- file.path(snakemake@params[["outdir"]], contrast)

# Extract gene Ensembl ID and log2FoldChange
# Convert Ensembl ID to Entrez ID
genename <- as.character(DESeq2[, "gene"])
if (snakemake@params[["species"]] == "homo_sapiens") { 
  gene_map <- select(org.Hs.eg.db, keys=genename, keytype="ENSEMBL", columns=c("ENTREZID"))
} else if (snakemake@params[["species"]] == "mus_musculus") {
  gene_map <- select(org.Mm.eg.db, keys=genename, keytype="ENSEMBL", columns=c("ENTREZID"))
}
colnames(gene_map)[1] <- "gene"
gsea_df <- inner_join(gene_map, DESeq2, by = "gene")
gsea_df <- gsea_df[, c("ENTREZID", "log2FoldChange")]
gsea_df <- na.omit(gsea_df)
gsea_df <- dplyr::distinct(gsea_df, ENTREZID, .keep_all=TRUE)
gsea_df <- gsea_df %>% 
  dplyr::arrange(desc(log2FoldChange))

# Data format for GSEA
geneList <- gsea_df[, "log2FoldChange"]
names(geneList) <- as.character(gsea_df[, "ENTREZID"])

# Perform GSEA enrichment analysis and plot results

# 1. GSEA GO
# 1.1 gseGO
GO_gseresult <- gseGO(
  geneList = geneList,
  OrgDb = orgdb,
  keyType = "ENTREZID",
  ont = "ALL",
  pvalueCutoff = plot_padj_cutoff,
  eps = 0.0
)
GO_gseresult_read <- DOSE::setReadable(GO_gseresult, 
                                       OrgDb = orgdb, 
                                       keyType = 'ENTREZID') #ENTREZID to gene Symbol
write.table(GO_gseresult_read@result, 
            file = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_results.xls", sep="_")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
save(GO_gseresult, 
     file = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_results.Rdata", sep="_")))

# 1.2 GSEA GO plots
## Dotplot
dotp1 <- dotplot(GO_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.pdf", sep="_")), 
       plot = dotp1,
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.png", sep="_")), 
       plot = dotp1,
       device = "png",
       width = 13, 
       height = 10)

dotp2 <- dotplot(GO_gseresult, split = "ONTOLOGY", font.size = 14, showCategory = 10) +
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot2.pdf", sep="_")),
       plot = dotp2, 
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot2.png", sep="_")), 
       plot = dotp2, 
       device = "png",
       width = 13, 
       height = 10)

## Ridgeplot
ridgep <- ridgeplot(GO_gseresult, showCategory = 10) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_ridgeplot.pdf", sep="_")), 
       plot = ridgep, 
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_ridgeplot.png", sep="_")), 
       plot = ridgep, 
       device = "png",
       width = 13, 
       height = 10)

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
reset_dir(gseaplot_dir)

setwd(gseaplot_dir)
for (i in go_num) {
  go <- go_list[i]
  p <- gseaplot2(GO_gseresult, geneSetID = i, color = "red", pvalue_table = T)
  
  png(filename = paste0(paste(contrast, "GSEA_GO_enrich_gseaplot", go, sep="_"), ".png"), 
      width = 13, 
      height = 10, 
      units = "in", 
      res=100)
  print(p)
  dev.off()
  
  pdf(file = paste0(paste(contrast, "GSEA_GO_enrich_gseaplot", go, sep="_"), ".pdf"), 
      width = 13, 
      height = 10)
  print(p)
  dev.off()
}
zip_file <-  paste0("../", paste(contrast, "GSEA_GO_enrich_gseaplot", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_GO_enrich_gseaplot", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)


# GSEA GO JSON GENERATION
# 1. GSEA GO Dotplot JSON
bubble_df <- dotp1$data %>% arrange(desc(Description))
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
  write(toJSON(bubble_json, auto_unbox=T, pretty=T), file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.json", sep="_")))
}

# 2. GSEA GO Faceted Dotplot JSON
facet_bubble_df <- dotp2$data
if(nrow(facet_bubble_df) > 0) {
  facet_bubble_df$ONTOLOGY <- factor(facet_bubble_df$ONTOLOGY, levels = c("BP", "CC", "MF"))
  facet_bubble_df <- facet_bubble_df %>% group_by(ONTOLOGY) %>% arrange(desc(Description), .by_group=TRUE) %>% ungroup()
  fb_size_values <- as.list(round(quantile(facet_bubble_df$Count, probs = c(0.33, 1.0))))
  facet_bubble_json <- list(
    width = 920, height = 680, outerBorder = TRUE, margin = list(top = 30, right = 140, bottom = 70, left = 310), color = list("#e76b66", "#5e9bd3"),
    dataset = list(list(source = list(description = as.list(as.character(facet_bubble_df$Description)), geneRatio = as.list(round(facet_bubble_df$GeneRatio, 6)), pAdjust = as.list(facet_bubble_df$p.adjust), count = as.list(facet_bubble_df$Count), ontology = as.list(as.character(facet_bubble_df$ONTOLOGY))))),
    series = list(list(type = "facetedEnrichmentBubble", datasetIndex = 0, encode = list(x = "geneRatio", y = "description", size = "count", color = "pAdjust", tooltip = list("description", "geneRatio", "pAdjust", "count", "ontology")), facet = list(field = "ontology", direction = "column", wrap = 1), sort = list(enabled = TRUE, field = "geneRatio", order = "desc"), border = list(show = TRUE), label = list(show = FALSE), facetTitleStyle = list(fontSize = "12", fontWeight = "bold"))),
    xAxis = list(list(type = "scaleLinear", position = "bottom", show = TRUE, name = "GeneRatio", min = 0, max = round(max(facet_bubble_df$GeneRatio) * 1.1, 4), tickCount = 5, axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisGrid = list(show = TRUE), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"), nameTextStyle = list(color = "#000", fontFamily = "sans-serif", fontSize = "14", fontWeight = "bold"))),
    yAxis = list(list(type = "scaleBand", position = "left", show = TRUE, name = "", axisLine = list(show = TRUE, lineStyle = list(width = 1)), axisTick = list(show = TRUE, lineStyle = list(width = 1)), axisGrid = list(show = TRUE), axisLabel = list(show = TRUE, color = "#000", fontFamily = "sans-serif", fontSize = "11"))),
    legend = list(), sizeLegend = list(show = TRUE, title = "Count", orient = "vertical", value = fb_size_values, right = 20, top = 250, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"),
    colorbar = list(show = TRUE, title = "p.adjust", width = 18, height = 90, site = "right", right = 20, top = 100, fontSize = 12, tickCount = 3, colors = list("#e76b66", "#5e9bd3"), reverse = FALSE, domain = as.list(range(facet_bubble_df$p.adjust)), useScientificNotation = (max(facet_bubble_df$p.adjust) < 0.001))
  )
  write(toJSON(facet_bubble_json, auto_unbox=T, pretty=T), file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot2.json", sep="_")))
}

# 2. GSEA KEGG
# 2.1 gseKEGG

print(spe1)

KEGG_gseresult <- gseKEGG(
  geneList = geneList,
  organism  = spe1, #hsa mmu
  keyType = "ncbi-geneid",
  pvalueCutoff = plot_padj_cutoff,
  eps = 0.0,
  use_internal_data = T
)

KEGG_gseresult@result <- merge(KEGG_gseresult@result, pathway_id2description, by.x="ID", by.y="pathway_id") %>% 
  dplyr::select(-Description) %>% 
  dplyr::rename(Description = pathway_description) %>%
  dplyr::arrange(p.adjust)
KEGG_gseresult@result <- KEGG_gseresult@result[, c(1,11,2,3,4,5,6,7,8,9,10)]
rownames(KEGG_gseresult@result) <- KEGG_gseresult@result$ID

KEGG_gseresult_read <- DOSE::setReadable(KEGG_gseresult, 
                                         OrgDb = orgdb, 
                                         keyType = 'ENTREZID') #ENTREZID to gene Symbol
write.table(KEGG_gseresult_read@result, 
            file = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_results.xls", sep="_")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
save(KEGG_gseresult, 
     file = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_results.Rdata", sep="_")))

# 2.2 GSEA KEGG plots
## Dotplot
dotp_kegg <- dotplot(KEGG_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.pdf", sep="_")), 
       plot = dotp_kegg,
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.png", sep="_")), 
       plot = dotp_kegg,
       device = "png",
       width = 13, 
       height = 10)

## Ridgeplot
ridgep <- ridgeplot(KEGG_gseresult, showCategory = 10) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_ridgeplot.pdf", sep="_")), 
       plot = ridgep, 
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_ridgeplot.png", sep="_")), 
       plot = ridgep, 
       device = "png",
       width = 13, 
       height = 10)

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
reset_dir(gseaplot_dir)

setwd(gseaplot_dir)
for (i in KEGG_num) {
  KEGG <- KEGG_list[i]
  p <- gseaplot2(KEGG_gseresult, geneSetID = i, color = "red", pvalue_table = T)
  
  png(filename = paste0(paste(contrast, "GSEA_KEGG_enrich_gseaplot", KEGG, sep="_"), ".png"), 
      width = 13, 
      height = 10, 
      units = "in", 
      res=100)
  print(p)
  dev.off()
  
  pdf(file = paste0(paste(contrast, "GSEA_KEGG_enrich_gseaplot", KEGG, sep="_"), ".pdf"), 
      width = 13, 
      height = 10)
  print(p)
  dev.off()
}
zip_file <-  paste0("../", paste(contrast, "GSEA_KEGG_enrich_gseaplot", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_KEGG_enrich_gseaplot", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)


# GSEA KEGG Dotplot JSON
bubble_df <- dotp_kegg$data %>% arrange(desc(Description))
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
  write(toJSON(bubble_json, auto_unbox=T, pretty=T), file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.json", sep="_")))
}

### pathview
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

pathview_dir <- file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"))
reset_dir(pathview_dir)
setwd(pathview_dir)

genelist <- as.numeric(DESeq2[, "log2FoldChange"])
names(genelist) <- DESeq2[, "gene"]
run_pathview(genelist, KEGG_list, spe1, n1=1, n2=1)

zip_file <-  paste0("../", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)

# 3. GSEA Reactome
# 3.1 gseReactome
Reactome_gseresult <- gsePathway(
  geneList = geneList,
  organism  = spe2, #human mouse
  pvalueCutoff = plot_padj_cutoff,
  eps = 0.0
)
Reactome_gseresult_read <- DOSE::setReadable(Reactome_gseresult, 
                                             OrgDb = orgdb, 
                                             keyType = 'ENTREZID') #ENTREZID to gene Symbol
write.table(Reactome_gseresult_read@result, 
            file = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_results.xls", sep="_")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
save(Reactome_gseresult, 
     file = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_results.Rdata", sep="_")))

# 3.2 GSEA Reactome plots
## Dotplot
dotp_reactome <- dotplot(Reactome_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_dotplot.pdf", sep="_")), 
       plot = dotp_reactome,
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_dotplot.png", sep="_")), 
       plot = dotp_reactome,
       device = "png",
       width = 13, 
       height = 10)

## Ridgeplot
ridgep <- ridgeplot(Reactome_gseresult, showCategory = 10) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100))
ggsave(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_ridgeplot.pdf", sep="_")), 
       plot = ridgep, 
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_ridgeplot.png", sep="_")), 
       plot = ridgep, 
       device = "png",
       width = 13, 
       height = 10)

## GSEAplot
### merged GSEAplot for top5
gseap <- gseaplot2(Reactome_gseresult, 1:5, color = "red", pvalue_table = T)
png(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_gseaplot_top5.png", sep="_")), 
    width = 16, 
    height = 15, 
    units = "in", 
    res=100)
gseap
dev.off()
pdf(file = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_gseaplot_top5.pdf", sep="_")), 
    width = 16, 
    height = 15)
gseap
dev.off()

### single GSEAplot for top20
Reactome_list <- Reactome_gseresult@result$ID
if (length(Reactome_list) < 20) {
  Reactome_num <- c(1:length(Reactome_list))
} else {
  Reactome_num <- c(1:20)
}

gseaplot_dir <- file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_gseaplot", sep="_"))
reset_dir(gseaplot_dir)

setwd(gseaplot_dir)
for (i in Reactome_num) {
  Reactome <- Reactome_list[i]
  p <- gseaplot2(Reactome_gseresult, geneSetID = i, color = "red", pvalue_table = T)
  
  png(filename = paste0(paste(contrast, "GSEA_Reactome_enrich_gseaplot", Reactome, sep="_"), ".png"), 
      width = 13, 
      height = 10, 
      units = "in", 
      res=100)
  print(p)
  dev.off()
  
  pdf(file = paste0(paste(contrast, "GSEA_Reactome_enrich_gseaplot", Reactome, sep="_"), ".pdf"), 
      width = 13, 
      height = 10)
  print(p)
  dev.off()
}
zip_file <-  paste0("../", paste(contrast, "GSEA_Reactome_enrich_gseaplot", sep="_"), ".zip")
zip_dir <- paste0("../", paste(contrast, "GSEA_Reactome_enrich_gseaplot", sep="_"))
system(paste("zip", zip_file, "./*"), intern=TRUE)
system(paste("rm -rf", zip_dir), intern=TRUE)
setwd(work_dir)


# GSEA Reactome Dotplot JSON
bubble_df <- dotp_reactome$data %>% arrange(desc(Description))
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
  write(toJSON(bubble_json, auto_unbox=T, pretty=T), file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_dotplot.json", sep="_")))
}
