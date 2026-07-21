#!/usr/bin/env Rscript
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(KEGG.db)
library(DOSE)
library(ReactomePA)

args = commandArgs(trailingOnly=T)

input.mc = args[1]      # ./results/05.diff_methy/FWGN20252361.Diff_region.sig.addGene
outdir = args[2]        # ./results/05.diff_methy/
prefix = args[3]        # FWGN20252361
species = args[4]       # hsa

if (species == "homo_sapiens"){
    spe = "hsa"
    spe2 = "human"
}

orgdb = "org.Hs.eg.db"
plot_padj_cutoff = 1.0
plot_q_cutoff = 1.0

df = read_tsv(input.mc, col_names = F)
gset = as.vector(unique(df$X8))
gene_list = unique(as.character(na.omit(bitr(
  gset,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = orgdb)[, 2])))

###################### GO ###################################
go_enrich_results <- enrichGO(gene = gene_list,
                              OrgDb = orgdb,
                              ont   = "ALL"  ,     #One of "BP", "MF"  "CC"  "ALL"
                              pvalueCutoff = plot_padj_cutoff, 
                              qvalueCutoff  = plot_q_cutoff,
                              readable      = TRUE)
write_tsv(go_enrich_results@result, 
          file = file.path(outdir, paste(prefix, "GO_enrich_results.xls", sep = "_")))
save(go_enrich_results, 
     file = file.path(outdir, paste(prefix, "GO_enrich_results.Rdata", sep = "_")))

# barplot
barp <- barplot(go_enrich_results, font.size = 14, showCategory = 30, color="pvalue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(barp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_barplot.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(barp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_barplot.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

# barplot2
barp <- barplot(go_enrich_results, split = "ONTOLOGY", font.size = 14, showCategory = 10, color="pvalue") +
  facet_grid(ONTOLOGY~., scale = "free") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(barp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_barplot2.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(barp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_barplot2.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

# dotplot
dotp <- enrichplot::dotplot(go_enrich_results, font.size = 14, showCategory = 30, color="pvalue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_dotplot.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_dotplot.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

# dotplot2
dotp <- enrichplot::dotplot(go_enrich_results, font.size = 14, split = 'ONTOLOGY', showCategory = 10, color="pvalue") +
  facet_grid(ONTOLOGY~., scale = "free") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_dotplot2.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "GO_enrich_dotplot2.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

## emapplot: Enrichment Map
#pt <- pairwise_termsim(go_enrich_results)
#pass <- TRUE
#cat_num <- 30
#while (pass) {
#  tryCatch({
#    emapp <- emapplot(pt, node_label = 'category', color="pvalue",
#                      showCategory = cat_num, color="pvalue")  # 'category', 'group', 'all' and 'none'.
#    pass <- FALSE
#  }, error = function(e) {
#    assign("cat_num", cat_num - 5, env=globalenv())
#  })
#}
#
#ggsave(emapp, 
#       filename = file.path(outdir, paste(prefix, "GO_enrich_emapplot.pdf", sep = "_")), 
#       device = "pdf", 
#       width = 13, 
#       height = 10)

##################### KEGG #########################
kegg_hsa_pathway_id2description = "/mnt/NexData_file/home/zsr/1.develop/1.WGBS/1.FWGN20252361/2.hsa_wgbs_203/resources/KEGG_database/KEGG_hsa_pathway_id2description.txt"
kegg_db = "/mnt/NexData_file/home/zsr/1.develop/1.WGBS/1.FWGN20252361/2.hsa_wgbs_203/resources/KEGG_database/KEGG_hsa_database_20230310.tar.gz"

pathway_id2description = read.table(kegg_hsa_pathway_id2description, sep = "\t", header = T)
kegg_enrich_results = enrichKEGG(
  gene = gene_list, 
  organism = spe, 
  pvalueCutoff = plot_padj_cutoff, 
  qvalueCutoff = plot_q_cutoff, 
  use_internal_data = T
)

kegg_enrich_results@result = merge(kegg_enrich_results@result, pathway_id2description, by.x = "ID", by.y = "pathway_id") %>% 
  dplyr::select(-Description) %>%
  dplyr::rename(Description = pathway_description) %>% 
  dplyr::arrange(p.adjust)

kegg_enrich_results@result = kegg_enrich_results@result %>% 
  dplyr::relocate(geneID, .after = last_col()) %>% 
  dplyr::mutate(Count =as.numeric(sub("(\\d+)/\\d+","\\1",GeneRatio)))

rownames(kegg_enrich_results@result) = kegg_enrich_results@result$ID

kk_read = DOSE::setReadable(kegg_enrich_results, OrgDb = orgdb, keyType = "ENTREZID")

write_tsv(kk_read@result, 
            file = file.path(outdir, paste(prefix, "KEGG_enrich_results.xls", sep = "_")))
save(kegg_enrich_results, 
     file = file.path(outdir, paste(prefix, "KEGG_enrich_results.Rdata", sep = "_")))

## barplot
barp <- barplot(kegg_enrich_results, font.size = 14, showCategory = 30, color="pvalue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(barp, 
       filename = file.path(outdir, paste(prefix, "KEGG_enrich_barplot.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(barp, 
       filename = file.path(outdir, paste(prefix, "KEGG_enrich_barplot.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

## dotplot
dotp <- enrichplot::dotplot(kegg_enrich_results, font.size = 14, showCategory = 30, color="pvalue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "KEGG_enrich_dotplot.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "KEGG_enrich_dotplot.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

### emapplot: Enrichment Map
#pt = pairwise_termsim(kegg_enrich_results)
#pass = TRUE
#cat_num = 30
#while (pass) {
#  tryCatch({
#    emapp <- emapplot(pt, node_label = 'category', color="pvalue", 
#                      showCategory = cat_num)  # 'category', 'group', 'all' and 'none'.
#    pass = FALSE
#  }, error = function(e){
#    assign("cat_num", cat_num - 5, env=globalenv())
#  })
#}
#
#ggsave(emapp, 
#       filename = file.path(outdir, paste(prefix, "KEGG_enrich_emapplot.pdf", sep = "_")), 
#       device = "pdf", 
#       width = 13, 
#       height = 10)

######################## Reactome #########################
reactome_enrich_results <- enrichPathway(gene  = gene_list,
                                         organism  = spe2, #human mouse
                                         pvalueCutoff = plot_padj_cutoff,
                                         qvalueCutoff = plot_q_cutoff,
                                         readable = TRUE)
write_tsv(reactome_enrich_results@result, 
          file = file.path(outdir, paste(prefix, "Reactome_enrich_results.xls", sep = "_")))
save(reactome_enrich_results, 
     file = file.path(outdir, paste(prefix, "Reactome_enrich_results.Rdata", sep = "_")))

## barplot
barp <- barplot(reactome_enrich_results, font.size = 14, showCategory = 30, color="pvalue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(barp,
       filename = file.path(outdir, paste(prefix, "Reactome_enrich_barplot.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(barp,
       filename = file.path(outdir, paste(prefix, "Reactome_enrich_barplot.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

## dotplot
dotp <- enrichplot::dotplot(reactome_enrich_results, font.size = 14, showCategory = 30, color="pvalue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  theme(legend.key.size = unit(10, "pt"),
        plot.margin = unit(c(1,1,1,1),'lines'))
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "Reactome_enrich_dotplot.png", sep = "_")), 
       device = "png", 
       width = 13, 
       height = 10)
ggsave(dotp, 
       filename = file.path(outdir, paste(prefix, "Reactome_enrich_dotplot.pdf", sep = "_")), 
       device = "pdf", 
       width = 13, 
       height = 10)

### emapplot
#pt = pairwise_termsim(reactome_enrich_results)
#pass <- TRUE
#cat_num <- 30
#text = "无法生成emapplot富集图！"
#emapp <- ggplot() +
#  annotate("text", x = 4, y = 25, size=8, label = text) +
#  theme_void() +
#  theme(panel.background = element_rect(fill = 'white'))
#
#while (pass & cat_num > 0) {
#  tryCatch({
#    emapp <- emapplot(pt, node_label = 'category', color="pvalue", 
#                      showCategory = cat_num)  # 'category', 'group', 'all' and 'none'.
#    pass <- FALSE
#  }, error = function(e) {
#    assign("cat_num", cat_num - 5, env=globalenv())
#  })
#}
#
#ggsave(emapp, 
#       filename = file.path(outdir, paste(prefix, "Reactome_enrich_emapplot.pdf", sep = "_")), 
#       device = "pdf", 
#       width = 13, 
#       height = 10)
