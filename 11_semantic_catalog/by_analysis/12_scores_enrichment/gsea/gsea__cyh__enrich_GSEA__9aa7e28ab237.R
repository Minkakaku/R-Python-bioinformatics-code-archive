log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

work_dir <- getwd()

options(stringsAsFactors = F) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(ggridges)
if (snakemake@params[["species"]] == "homo_sapiens") {
  orgdb <- "org.Hs.eg.db"
  spe1 <- "hsa"
  spe2 <- "human"
  library(org.Hs.eg.db)
  library(KEGG.db)
  pathway_id2description <- read.table(snakemake@params[["kegg_hsa_pathway_id2description"]], sep="\t", header=TRUE)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  orgdb <- "org.Mm.eg.db"
  spe1 <- "mmu"
  spe2 <- "mouse"
  library(org.Mm.eg.db)
  library(KEGG.db)
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
dotp <- dotplot(GO_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.pdf", sep="_")), 
       plot = dotp,
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "GO", paste(contrast, "GSEA_GO_enrich_dotplot.png", sep="_")), 
       plot = dotp,
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
system(paste("mkdir", gseaplot_dir))

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

# 2. GSEA KEGG
# 2.1 gseKEGG
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
dotp <- dotplot(KEGG_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.pdf", sep="_")), 
       plot = dotp,
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_dotplot.png", sep="_")), 
       plot = dotp,
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
system(paste("mkdir", gseaplot_dir))

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

### pathview
pathview_dir <- file.path(outdir, "KEGG", paste(contrast, "GSEA_KEGG_enrich_pathway", sep="_"))
system(paste("mkdir", pathview_dir))
setwd(pathview_dir)

genelist <- as.numeric(DESeq2[, "log2FoldChange"])
names(genelist) <- DESeq2[, "gene"]

problematic_pathways <- c("hsa05206", "mmu05206", "hsa04723", "mmu04723", "hsa04215", "mmu04215")
sel = which(KEGG_list %in% problematic_pathways)
if (length(sel)) {
  KEGG_list <- KEGG_list[-sel]
}

if (length(KEGG_list) < 20) {
  select_pathway <- KEGG_list[1:length(KEGG_list)]
} else {
  select_pathway <-KEGG_list[1:20]
}

pv.out <- pathview(
  gene.data     = genelist,
  pathway.id    = select_pathway,
  species       = spe1 , # hsa mmu
  gene.idtype   = "ENSEMBL",
  kegg.native   = T,
  new.signature = F,
)

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
dotp <- dotplot(Reactome_gseresult, font.size = 14, showCategory = 30) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1, 1, 1, 1), 'lines'))
ggsave(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_dotplot.pdf", sep="_")), 
       plot = dotp,
       device = "pdf",
       width = 13, 
       height = 10)
ggsave(filename = file.path(outdir, "Reactome", paste(contrast, "GSEA_Reactome_enrich_dotplot.png", sep="_")), 
       plot = dotp,
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
system(paste("mkdir", gseaplot_dir))

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
