library(dplyr)
library(data.table)
library(stringr)   
library(ggrepel)   
library(gridExtra) 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
diff.data <- fread(file.choose())
gsea.enrich <- fread(file.choose())
kegg.data = fread("KEGG_annotation.tsv") %>%
  distinct(`Protein accession`, Terms)

gmt.list <- split(kegg.data$`Protein accession`, kegg.data$Terms)
geneList <- sort(setNames(diff.data$log2FC, diff.data$`Protein accession`) %>% na.omit(), decreasing = TRUE)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
gsea.enrich$Global_Rank <- 1:nrow(gsea.enrich) 

c_keywords <- c("Ferroptosis", "Apoptosis", "HIF-1 signaling", "PI3K-Akt", "FoxO")
gsea.enrich$Priority <- 3 

gsea.enrich$Priority[str_detect(gsea.enrich$pathway, regex(paste(c_keywords, collapse = "|"), ignore_case = TRUE))] <- 1
gsea.enrich$Priority[str_detect(gsea.enrich$pathway, regex("Oxidative phosphorylation", ignore_case = TRUE)) & gsea.enrich$pval < 0.05] <- 1

gsea_sorted <- gsea.enrich[order(gsea.enrich$Priority, gsea.enrich$Global_Rank), ]

kegg.gsea.sig <- head(gsea_sorted, 10)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
calc_ES <- function(geneList, pathway, p = 1) {
  hits <- names(geneList) %in% pathway
  Nh <- sum(hits)
  Nm <- length(geneList) - Nh
  Phit <- ifelse(hits, abs(geneList)^p / sum(abs(geneList[hits])^p), 0)
  Pmiss <- ifelse(!hits, 1 / Nm, 0)
  RES <- cumsum(Phit - Pmiss)
  ES <- ifelse(max(RES) > abs(min(RES)), max(RES), min(RES))
  list(ES = ES, DF = data.frame(rank = seq_along(geneList), score = geneList, hit = as.integer(hits), runningES = RES))
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
gsea.plot <- function(pathway_genes, geneList, gsea.res) {
  gsea.res <- gsea.res[1, ]
  pathway_name <- gsea.res$pathway
  es_result <- calc_ES(geneList, pathway_genes[[pathway_name]])
  max_es <- es_result$ES
  df <- es_result$DF
  yn <- round(max_es / 10, 2)
  
  # Annotation table position
  if (max_es > 0) {
    xmin <- 0.6 * max(df$rank)
    xmax <- 0.9 * max(df$rank)
  } else {
    xmin <- 0.1 * max(df$rank)
    xmax <- 0.4 * max(df$rank)
  }
  
  ggplot(df, aes(x = rank, y = runningES)) +
    geom_tile(aes(y = 0, fill = score), height = 2 * yn) +
    scale_fill_gradientn(colors = c("steelblue", "white", "red"), limits = c(-1, 1), oob = scales::squish) +
    geom_line(size = 1.1, color = "#1f78b4") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(data = df[df$hit == 1, ], aes(x = rank, xend = rank, y = -yn, yend = yn), color = "black", alpha = 0.6) +
    geom_point(data = df[which.max(abs(df$runningES)), ], aes(x = rank, y = max_es), size = 4, color = "firebrick") +
    geom_text_repel(
      data = df[which.max(abs(df$runningES)), ],
      aes(x = rank, y = max_es, label = paste0("ES: ", round(max_es, 3))),
      nudge_y = yn, size = 5, color = "firebrick", segment.color = "firebrick", segment.size = 0.5, show.legend = FALSE
    ) +
    annotation_custom(
      grob = tableGrob(
        gsea.res[, c("pval", "NES", "size")] %>% mutate(pval = signif(pval, 2), NES = round(NES, 3)),
        rows = NULL
      ),
      xmin = xmin, xmax = xmax, ymin = 0.5 * max_es, ymax = max_es
    ) +
    labs(x = "Rank in Ordered Gene List", y = "Running Enrichment Score", title = pathway_name) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(), axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), plot.title = element_text(size = 18, face = "bold"),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if (!dir.exists("result/GSEA_plot")) {
  dir.create("result/GSEA_plot", recursive = TRUE)
}

if (nrow(kegg.gsea.sig) > 0) {
  lapply(split(kegg.gsea.sig, kegg.gsea.sig$pathway), function(x) {
    tryCatch({
      p = gsea.plot(pathway_genes = gmt.list, geneList = geneList, gsea.res = x)
      ggsave(plot = p, paste0("GSEA_plot_", x$pathway[1], ".png"), width = 8, height = 5)
      ggsave(plot = p, paste0("GSEA_plot_", x$pathway[1], ".pdf"), width = 8, height = 5)
    }, error = function(e) {
      message("跳过通路 ", x$pathway[1], "，原因为：", e$message)
    })
  })
} else {
  message("未筛选到符合条件的通路。")
}