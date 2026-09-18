library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(fgsea)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
get_gsdata <- function(x, geneSetID) {
  if (length(geneSetID) == 1) {
    gsdata <- enrichplot:::gsInfo(x, geneSetID)
    return(gsdata)
  } 
  lapply(geneSetID, enrichplot:::gsInfo, object = x) |>
    yulab.utils::rbindlist()
}

gseaplot <- function(gseadata, pathwayID, customtitle = FALSE, showNES = TRUE, showPvalue = TRUE,
                     showFDR = TRUE, showCase = TRUE, showCon = TRUE, ES_linetype = "solid", ES_linecolor = "#41AB5D",
                     ES_linewidth = 1.5, fz_linetype = "solid", fz_linecolor = "grey", fz_linewidth = 1.5, p1ytitlefontface = "bold",
                     p1ytextsize = 12, p1ytextfontface = "plain", p1title = "AAAAAAAAAAA", p1titlesize = 20, p1ytitlesize = 16,
                     p1titlefontface = "bold", p1titlecolor = "black", p1titlehjust = 0.5, p1titlevjust = 0.5, NESx = 0.85, NESy = 0.9,
                     p1bordersize = 1.5, NESsize = 5, NEScolor = "black", NESfontface = "italic", Pvaluex = 0.85, Pvaluey = 0.8, Pvaluesize = 5,
                     Pvaluecolor = "black", Pvaluefontface = "italic", FDRx = 0.85, FDRy = 0.7, FDRsize = 5, FDRcolor = "black",
                     FDRfontface = "italic", NESdigit = 2, Pvaluedigit = 2, FDRdigit = 2, p2_linetype = 1, p2_linecolor = "black",
                     p2_linewidth = 0.8, rec_height_ratio = 0.3, rec_height_alpha = 0.8, rec_col_high = "Reds", rec_col_low = "Blues",
                     p3ytitlesize = 16, p3ytitlefontface = "plain", p3xtitlesize = 16, p3xtitlefontface = "bold",
                     p3xtextsize = 12, p3xtextfontface = "plain", p3bordersize = 1, casename = "Case",
                     conname = "Con", Casex = 0.1, Casey = 0.6, Casesize = 5, Casecolor = "red", Casefontface = "italic",
                     Conx = 0.9, Cony = -0.6, Consize = 5, Concolor = "blue", Confontface = "italic",
                     p1ratio = 1, p2ratio = 0.2, p3ratio=0.5){
  
  maintitle = gseadata@result[pathwayID, "Description"]
  nes = gseadata@result[pathwayID, "NES"]
  enrichmentScore = gseadata@result[pathwayID, "enrichmentScore"]
  pvalue = gseadata@result[pathwayID, "pvalue"]
  fdr = gseadata@result[pathwayID, "p.adjust"]
  
  if (!customtitle) {
    Maintitle = maintitle
  } else {
    Maintitle = p1title
  }
  if (showNES) {
    NES = paste0("NES = ", round(nes, NESdigit))
    if (enrichmentScore >= 0) {
      NESx = length(gseadata@geneList) * NESx
      NESy = enrichmentScore * NESy
    } else {
      NESx = length(gseadata@geneList) * (1-NESx)
      NESy = enrichmentScore * (NESy-0.2)
    }
  } else {
    NES = enrichmentScore
    NESsize = 0
    NESx = length(gseadata@geneList) * 2
    NESy = enrichmentScore * (NESy-0.2)
    NEScolor = "white"
  }
  
  if (showPvalue) {
    Pvalue = paste0("P value = ", format(pvalue, scientific = TRUE, digits = Pvaluedigit))
    if (enrichmentScore >= 0) {
      Pvaluex = length(gseadata@geneList) * Pvaluex
      Pvaluey = enrichmentScore * Pvaluey
    } else {
      Pvaluex = length(gseadata@geneList) * (1-Pvaluex)
      Pvaluey = enrichmentScore * Pvaluey
    }
  } else {
    Pvalue = enrichmentScore
    Pvaluesize = 0
    Pvaluex = length(gseadata@geneList) * 2
    Pvaluey = enrichmentScore * Pvaluey
    Pvaluecolor = "white"
  }
  
  if (showFDR) {
    FDR = paste0("FDR = ", format(fdr, scientific = TRUE, digits = FDRdigit))
    if (enrichmentScore >= 0) {
      FDRx = length(gseadata@geneList) * FDRx
      FDRy = enrichmentScore * FDRy
    } else {
      FDRx = length(gseadata@geneList) * (1-FDRx)
      FDRy = enrichmentScore * (FDRy+0.2)
    }
  } else {
    FDR = enrichmentScore
    FDRsize = 0
    FDRx = length(gseadata@geneList) * 2
    FDRy = enrichmentScore * (FDRy+0.2)
    FDRcolor = "white"
  }
  
  gsdata <- get_gsdata(gseadata, pathwayID)
  
  # p1 部分
  p1 <- ggplot(gsdata, aes(x = x)) + xlab(NULL) +
    geom_hline(yintercept = 0, color = fz_linecolor, linewidth = fz_linewidth, linetype = fz_linetype) +
    geom_line(aes(y = runningScore), color= ES_linecolor, linewidth = ES_linewidth, linetype = ES_linetype) +
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.25,"cm"),
          panel.border = element_rect(linewidth = p1bordersize, fill = NA),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = p1ytextsize, face = p1ytextfontface),
          axis.title.y = element_text(size = p1ytitlesize, face = p1ytitlefontface),
          plot.title=element_text(face=p1titlefontface, color=p1titlecolor, size=p1titlesize, hjust=p1titlehjust, vjust=p1titlevjust)
    ) +
    ylab("Enrichment score") + labs(title = Maintitle) +
    scale_x_continuous(expand=c(0,0), limits = c(0,length(gseadata@geneList))) + 
    scale_y_continuous(expand=c(0,0)) +
    annotate(geom = "text", x = NESx, y = NESy, size = NESsize, label = NES, color = NEScolor, fontface = NESfontface) +
    annotate(geom = "text", x = Pvaluex, y = Pvaluey, size = Pvaluesize, label = Pvalue, color = Pvaluecolor, fontface = Pvaluefontface) +
    annotate(geom = "text", x = FDRx, y = FDRy, size = FDRsize, label = FDR, color = FDRcolor, fontface = FDRfontface)
  
  v <- seq(1, sum(gsdata$position), length.out=9)
  inv <- findInterval(rev(cumsum(gsdata$position)), v)
  if (min(inv) == 0) inv <- inv + 1
  col <- c(rev(brewer.pal(5, rec_col_low)), brewer.pal(5, rec_col_high))
  
  # p2 部分
  p2 <- ggplot(gsdata, aes(x = x)) +
    geom_linerange(aes(ymin=ymin, ymax=ymax), linewidth = p2_linewidth, linetype = as.numeric(p2_linetype), color = p2_linecolor) +
    xlab(NULL) + ylab(NULL) + theme_classic(11) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_blank(),
          axis.line.y = element_line(color = "white"),
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_y_continuous(expand=c(0,0)) 
  
  ymin <- min(p2$data$ymin)
  yy <- (max(p2$data$ymax) - min(p2$data$ymin)) * rec_height_ratio + min(p2$data$ymin)
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, xmax = xmax, col = col[unique(inv)])
  
  p2 <- p2 + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=I(col)),
                       data=d, alpha=rec_height_alpha, inherit.aes=FALSE) +
    scale_x_continuous(expand=c(0,0),limits = c(0,rev(xmax)[1])) 
  
  if (showCase) {
    Casex = length(gseadata@geneList) * Casex
    Casey = max(abs(gseadata@geneList)) * Casey
  } else {
    Casesize = 0
    Casex = length(gseadata@geneList) * 2
    Casey = max(abs(gseadata@geneList)) * 0.6
    Casecolor = "white"
  }
  
  if (showCon) {
    Conx = length(gseadata@geneList) * Conx
    Cony = max(abs(gseadata@geneList)) * Cony
  } else {
    Consize = 0
    Conx = length(gseadata@geneList) * 2
    Cony = max(abs(gseadata@geneList)) * -0.6
    Concolor = "white"
  }
  
  df2 <- p1$data 
  df2$y <- p1$data$geneList[df2$x]
  df2$col <- rep(col[unique(inv)], rev(table(inv)))
  
  # p3 部分
  p3 <- ggplot() + 
    geom_segment(data=df2, aes(x=x, xend=x, y=y, yend=0), color=df2$col) +
    theme_bw() + coord_cartesian(expand = 0) +
    theme(axis.ticks.y = element_blank(),
          axis.line = element_line(colour = "black", linewidth = p3bordersize),
          axis.ticks.length = unit(0.25,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = p3xtextsize, face = p3xtextfontface),
          axis.title.x = element_text(size = p3xtitlesize, face = p3xtitlefontface),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = p3ytitlesize, face = p3ytitlefontface),
          plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm")
    ) +
    scale_y_continuous(expand=c(0,0),limits = c(max(abs(gseadata@geneList))*-1.2,max(abs(gseadata@geneList))*1.2)) +
    scale_x_continuous(expand=c(0,0),limits = c(0,length(gseadata@geneList))) +
    xlab("Ranked in Ordered Dataset") + ylab("Ranked List Metric") +
    annotate(geom = "text", x = Casex, y = Casey, size = Casesize, label = casename, color = Casecolor, fontface = Casefontface) +
    annotate(geom = "text", x = Conx, y = Cony, size = Consize, label = conname, color = Concolor, fontface = Confontface)
  
  # 强行锁定 patchwork::plot_layout 规避因命名空间冲突报错
  p.final = p1 + p2 + p3 + patchwork::plot_layout(nrow = 3, heights = c(p1ratio, p2ratio, p3ratio))
  
  return(p.final)
}

# ==============================================================================
# 步骤一本地 KEGG 字典
# ==============================================================================
kegg_anno_raw <- fread(file.choose())

my_kegg_gmt_raw <- data.frame(
  gs_name = kegg_anno_raw$Terms,
  gene    = kegg_anno_raw$`Protein accession`
)

dictionary_translated <- bitr(my_kegg_gmt_raw$gene, 
                              fromType = "UNIPROT", 
                              toType = "ENTREZID", 
                              OrgDb = org.Hs.eg.db)

my_kegg_gmt_final <- merge(my_kegg_gmt_raw, dictionary_translated, by.x = "gene", by.y = "UNIPROT")
my_kegg_gmt <- data.frame(
  gs_name = my_kegg_gmt_final$gs_name,
  gene    = my_kegg_gmt_final$ENTREZID 
)

my_kegg_gmt$gene <- as.character(my_kegg_gmt$gene)

# ==============================================================================
# 步骤二：（差异表达矩阵总表）提取log2
# ==============================================================================
deg_total <- fread(file.choose()) 

gene_ids <- bitr(deg_total$`Protein accession`, 
                 fromType = "UNIPROT", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

deg_with_id <- merge(deg_total, gene_ids, by.x = "Protein accession", by.y = "UNIPROT")

gene_list <- deg_with_id$log2FC          
names(gene_list) <- deg_with_id$ENTREZID
names(gene_list) <- as.character(names(gene_list))
gene_list <- sort(gene_list, decreasing = TRUE)

# ==============================================================================
# 步骤三：（交集基因验证）
# ==============================================================================
common_genes <- intersect(names(gene_list), my_kegg_gmt$gene)
length(common_genes)

# ==============================================================================
# 步骤四：运行本地 GSEA 运算
# ==============================================================================
gsea_res <- GSEA(geneList = gene_list, 
                 TERM2GENE = my_kegg_gmt, 
                 pvalueCutoff = 1) 

save(gsea_res, file = "my_gsea_result.Rdata")

# ==============================================================================
# 步骤五：target proteins---------
# ==============================================================================
target_proteins <- fread(file.choose()) 

targets_translated <- bitr(target_proteins$`Protein accession`, 
                           fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
my_targets_numeric <- as.character(unique(targets_translated$ENTREZID))

# target proteins gsea enrich data-----
gsea_base_df <- as.data.frame(gsea_res@result)
gsea_base_df$Matched_Count <- 0
gsea_base_df$Matched_Genes <- ""

for (i in 1:nrow(gsea_base_df)) {
  p_id <- gsea_base_df$ID[i]
  pathway_all_genes <- as.character(my_kegg_gmt$gene[my_kegg_gmt$gs_name == p_id])
  intersect_genes <- intersect(my_targets_numeric, pathway_all_genes)
  
  if(length(intersect_genes) > 0){
    symbols_df <- bitr(intersect_genes, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)
    gsea_base_df$Matched_Genes[i] <- paste(symbols_df$SYMBOL, collapse = "|")
  }
  gsea_base_df$Matched_Count[i] <- length(intersect_genes)
}

# ==============================================================================
# ==============================================================================
# target proteins pathway
target_pathways <- gsea_base_df[gsea_base_df$Matched_Count > 0, ]

c_keywords <- c("Ferroptosis", "Apoptosis", "HIF-1", "PI3K-Akt", "FoxO", "Oxidative phosphorylation")
target_pathways$Priority <- ifelse(
  str_detect(target_pathways$ID, regex(paste(c_keywords, collapse = "|"), ignore_case = TRUE)), 
  1, 
  2
)

# 3. 先按 Priority 从小到大排（1排在2前面），再按靶点数量从大到小排（-Matched_Count）
target_pathways <- target_pathways[order(target_pathways$Priority, -target_pathways$Matched_Count), ]

final_plots_df <- head(target_pathways, 10)



# ==============================================================================
# GESA-----------------
# ==============================================================================
if(nrow(final_plots_df) > 0) {
  for (i in 1:nrow(final_plots_df)) {
    p_id <- final_plots_df$ID[i]
    pathway_genes_from_anno <- as.character(my_kegg_gmt$gene[my_kegg_gmt$gs_name == p_id])
    
    p <- plotEnrichment(pathway = pathway_genes_from_anno, stats = gene_list) +
      labs(title = p_id, 
           subtitle = paste0("NES = ", round(final_plots_df$NES[i], 2), 
                             "   pvalue = ", format(final_plots_df$pvalue[i], digits = 3),
                             "\n(Driven by ", final_plots_df$Matched_Count[i], " target proteins: ", final_plots_df$Matched_Genes[i], ")")) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 9))
    
    safe_filename <- gsub("[[:punct:]|[:space:]]", "_", p_id)
    ggsave(filename = paste0("GSEA_", safe_filename, ".pdf"), plot = p, width = 6, height = 4.5)
  }
}