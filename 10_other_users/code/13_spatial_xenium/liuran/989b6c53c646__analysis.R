
# 260120 - analysis--------------
all_data <- readRDS(file.path("data_for_analysis", "all_data_h5ad_ROI_0107.rds"))


### setting dir ------------
results_outdir <- "20260120_analysis"
dir.create(results_outdir)

plan <- c(
  "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
  "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead",
  "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c",
  "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1",
  "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E",
  "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d",
  "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9",
  "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6",
  "#005995", "#fa625f", "#f17d80"
)

### dotplot --------
pdf(file.path(results_outdir, str_glue("Dimplot all_cell group.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "group", label = TRUE, reduction = "umap", cols = plan) + NoLegend()
dev.off()


pdf(file.path(results_outdir, str_glue("Dimplot all_cell sample.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sample", label = TRUE, reduction = "umap", cols = plan) + NoLegend()
dev.off()


pdf(file.path(results_outdir, str_glue("Dimplot all_cell ROI.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "ROI", label = TRUE, reduction = "umap", cols = plan) + NoLegend()
dev.off()






all_data_filt <- all_data[, all_data@meta.data$ROI != "Other"]

DefaultAssay(all_data_filt) <- "RNA"

all_data_filt <-
  all_data_filt %>%
  NormalizeData(scale.factor = floor(mean(all_data_filt$nCount_RNA))) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(
    vars.to.regress = c(
      "nCount_RNA"
    )
  )

all_data_filt <- all_data_filt %>% RunPCA()

pdf(file.path(results_outdir, "pca plot.pdf"), height = 6, width = 8)
ElbowPlot(all_data_filt, ndims = 30)
dev.off()


dims <- 13

library(future)
plan(multisession, workers = 24)
options(future.globals.maxSize = 50 * 1024^3) # 10Gb

all_data_filt <- RunUMAP(all_data_filt, dims = 1:dims)

pdf(file.path(results_outdir, str_glue("Dimplot all_cell sample - filt.pdf")), height = 5, width = 6)
DimPlot(all_data_filt, group.by = "sample", label = TRUE, cols = plan) + NoLegend()
dev.off()

pdf(file.path(results_outdir, str_glue("Dimplot all_cell sample - filt - split by sample.pdf")), height = 10, width = 20)
DimPlot(all_data_filt, group.by = "sample", cols = plan, split.by = "sample", ncol = 6)
dev.off()


pdf(file.path(results_outdir, str_glue("Dimplot all_cell group - filt.pdf")), height = 5, width = 6)
DimPlot(all_data_filt, group.by = "group", label = TRUE, reduction = "umap", cols = plan) + NoLegend()
dev.off()

pdf(file.path(results_outdir, str_glue("Dimplot all_cell group - filt - split by group.pdf")), height = 10, width = 15)
DimPlot(all_data_filt, group.by = "group", reduction = "umap", cols = plan, split.by = "group", ncol = 3)
dev.off()


pdf(file.path(results_outdir, str_glue("Dimplot all_cell ROI - filt.pdf")), height = 5, width = 6)
DimPlot(all_data_filt, group.by = "ROI", label = TRUE, reduction = "umap", cols = plan) + NoLegend()
dev.off()

pdf(file.path(results_outdir, str_glue("Dimplot all_cell ROI - filt - split by ROI.pdf")), height = 5, width = 15)
DimPlot(all_data_filt, group.by = "ROI", reduction = "umap", cols = plan, split.by = "ROI")
dev.off()





file.copy("ROI_analysis - 2025-12-31/spatial plot gene/spatial plot - App.pdf", results_outdir)
file.copy("ROI_analysis - 2025-12-31/spatial plot gene/spatial plot - Mapt.pdf", results_outdir)



sp_genes <- c("Mapt")
sp_genes <- c("App")
sp_cells <- c("DG", "CA3", "CA1")
annotation_col <- "ROI"

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project) { # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")

    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>%
      right_join(dat_expr, by = "cell_id")


    all_cells <- dat[[annotation_col]] %>% unique()
    other_cells <- all_cells %>% setdiff(sp_cells)
    color_for_cells <- c(
      structure(names = sp_cells, plan[3:5]),
      structure(names = other_cells, rep("gray95", times = length(other_cells))),
      Others = "gray95"
    )

    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene) {
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])

        dat[[annotation_col]][dat[[annotation_col]] %in% other_cells] <- "Others"
        dat[[gene]][dat[[annotation_col]] == "Others"] <- NA

        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(
            fill = .data[[gene]],
            # color = .data[[annotation_col]],
            group = cell_id
          )) +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          # scale_fill_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })

# export
pdf(file.path(results_outdir, str_glue("spatial plot - App - filt.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()

pdf(file.path(results_outdir, str_glue("spatial plot - Mapt - filt.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()




pdf(file.path(results_outdir, str_glue("Vlnplot by ROI - Tau and App - filt.pdf")), height = 5, width = 6)
VlnPlot(all_data_filt, features = "Mapt", group.by = "ROI", pt.size = 0, cols = plan)
VlnPlot(all_data_filt, features = "App", group.by = "ROI", pt.size = 0, cols = plan)
dev.off()




mapt_cells <- WhichCells(all_data_filt, expression = Mapt > 0)
app_cells <- WhichCells(all_data_filt, expression = App > 0)

all_data_mapt <- subset(all_data_filt, cells = mapt_cells)
all_data_app <- subset(all_data_filt, cells = app_cells)


pdf(file.path(results_outdir, str_glue("Vlnplot by ROI - Tau and App - expression_filt.pdf")), height = 5, width = 6)
VlnPlot(all_data_mapt, features = "Mapt", group.by = "ROI", pt.size = 0, cols = plan)
VlnPlot(all_data_app, features = "App", group.by = "ROI", pt.size = 0, cols = plan)
dev.off()






# 260421 - analysis ----------
## diff analysis ----------
# M vs K    B vs M    D vs M    H vs M    B vs D    B vs H    D vs H
## setting dir ------------
results_outdir <- "20260421_analysis/diff_analysis_CA1"
dir.create(results_outdir, recursive = T)



library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library("clusterProfiler")
# library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
library(gtools)



# 读取数据
all_data <- readRDS(file.path("data_for_analysis", "all_data_h5ad_ROI_0107.rds"))


groups <- unique(all_data@meta.data$group)

all_combns_ls <- permutations(
  n = length(groups),
  r = 2,
  v = seq_along(groups),
  repeats.allowed = FALSE
) %>% apply(1, c, simplify = FALSE)



### DG -----------
all_combns_ls %>%
  walk(function(x) {
    grp <- groups[x]
    
    diff_markers <- FindMarkers(
      JoinLayers(all_data_DG),
      ident.1 = grp[1],
      ident.2 = grp[2],
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      group.by = "group"
    )
    if (nrow(diff_markers) == 0) {
      cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
      return(NULL)
    }
    
    diffmarkers_up <- diff_markers %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC))
    diffmarkers_down <- diff_markers %>%
      filter(avg_log2FC < 0) %>%
      arrange(desc(avg_log2FC))
    
    diff_outdir <- file.path(results_outdir, str_glue("DG diff analysis/DG {grp[1]} vs {grp[2]}"))
    dir.create(diff_outdir, recursive = TRUE)
    write.csv(diffmarkers, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker.csv")))
    write.csv(diffmarkers_up, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker up.csv")))
    write.csv(diffmarkers_down, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker down.csv")))
    
    
    # GO up
    enrichGO_up_ls <- enrichGO(
      rownames(diffmarkers_up),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichGO_up_ls) > 0) {
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched GO up terms found. Skipping plot generation.\n"))
    }
    
    
    # GO down
    enrichGO_down_ls <- enrichGO(
      rownames(diffmarkers_down),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichGO_down_ls) > 0) {
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched GO down terms found. Skipping plot generation.\n"))
    }
    
    
    # KEGG up
    de_up_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_up),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_up_ls <- enrichKEGG(
      de_up_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichKEGG_up_ls) > 0) {
      enrichKEGG_up_dot <- enrichplot::dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched KEGG up terms found. Skipping plot generation.\n"))
    }
    
    
    # KEGG down
    de_down_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_down),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_down_ls <- enrichKEGG(
      de_down_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichKEGG_down_ls) > 0) {
      enrichKEGG_down_dot <- enrichplot::dotplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_down_bar <- barplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched KEGG down terms found. Skipping plot generation.\n"))
    }
  })












### CA3 -----------
all_combns_ls %>%
  walk(function(x) {
    grp <- groups[x]
    
    diff_markers <- FindMarkers(
      JoinLayers(all_data_CA3),
      ident.1 = grp[1],
      ident.2 = grp[2],
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      group.by = "group"
    )
    if (nrow(diff_markers) == 0) {
      cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
      return(NULL)
    }
    
    diffmarkers_up <- diff_markers %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC))
    diffmarkers_down <- diff_markers %>%
      filter(avg_log2FC < 0) %>%
      arrange(desc(avg_log2FC))
    
    diff_outdir <- file.path(results_outdir, str_glue("CA3 diff analysis/CA3 {grp[1]} vs {grp[2]}"))
    dir.create(diff_outdir, recursive = TRUE)
    write.csv(diffmarkers, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker.csv")))
    write.csv(diffmarkers_up, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker up.csv")))
    write.csv(diffmarkers_down, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker down.csv")))
    
    
    # GO up
    enrichGO_up_ls <- enrichGO(
      rownames(diffmarkers_up),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichGO_up_ls) > 0) {
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched GO up terms found. Skipping plot generation.\n"))
    }
    
    
    # GO down
    enrichGO_down_ls <- enrichGO(
      rownames(diffmarkers_down),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichGO_down_ls) > 0) {
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched GO down terms found. Skipping plot generation.\n"))
    }
    
    
    # KEGG up
    de_up_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_up),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_up_ls <- enrichKEGG(
      de_up_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichKEGG_up_ls) > 0) {
      enrichKEGG_up_dot <- enrichplot::dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched KEGG up terms found. Skipping plot generation.\n"))
    }
    
    
    # KEGG down
    de_down_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_down),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_down_ls <- enrichKEGG(
      de_down_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichKEGG_down_ls) > 0) {
      enrichKEGG_down_dot <- enrichplot::dotplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_down_bar <- barplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched KEGG down terms found. Skipping plot generation.\n"))
    }
  })







### CA1 -----------
all_combns_ls %>%
  walk(function(x) {
    grp <- groups[x]
    
    diff_markers <- FindMarkers(
      JoinLayers(all_data_CA1),
      ident.1 = grp[1],
      ident.2 = grp[2],
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      group.by = "group"
    )
    if (nrow(diff_markers) == 0) {
      cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
      return(NULL)
    }
    
    diffmarkers_up <- diff_markers %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC))
    diffmarkers_down <- diff_markers %>%
      filter(avg_log2FC < 0) %>%
      arrange(desc(avg_log2FC))
    
    diff_outdir <- file.path(results_outdir, str_glue("CA1 diff analysis/CA1 {grp[1]} vs {grp[2]}"))
    dir.create(diff_outdir, recursive = TRUE)
    write.csv(diffmarkers, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker.csv")))
    write.csv(diffmarkers_up, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker up.csv")))
    write.csv(diffmarkers_down, file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} group diffmarker down.csv")))
    
    
    # GO up
    enrichGO_up_ls <- enrichGO(
      rownames(diffmarkers_up),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichGO_up_ls) > 0) {
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched GO up terms found. Skipping plot generation.\n"))
    }
    
    
    # GO down
    enrichGO_down_ls <- enrichGO(
      rownames(diffmarkers_down),
      OrgDb = org.Mm.eg.db,
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL"
    )
    write.csv(enrichGO_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichGO_down_ls) > 0) {
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down dotplot.pdf")),
        height = 6,
        width = 8
      )
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichGO_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} Go down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched GO down terms found. Skipping plot generation.\n"))
    }
    
    
    # KEGG up
    de_up_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_up),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_up_ls <- enrichKEGG(
      de_up_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_up_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichKEGG_up_ls) > 0) {
      enrichKEGG_up_dot <- enrichplot::dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_up_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_up_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG up barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched KEGG up terms found. Skipping plot generation.\n"))
    }
    
    
    # KEGG down
    de_down_genes_id_for_enrichmenet <-
      clusterProfiler::bitr(
        rownames(diffmarkers_down),
        OrgDb = "org.Mm.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrichKEGG_down_ls <- enrichKEGG(
      de_down_genes_id_for_enrichmenet,
      organism = "mmu",
      keyType = "ncbi-geneid",
      # keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
      # use_internal_data = TRUE
    )
    enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(enrichKEGG_down_ls@result, file = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down enriment table.csv")), row.names = FALSE)
    
    
    if (nrow(enrichKEGG_down_ls) > 0) {
      enrichKEGG_down_dot <- enrichplot::dotplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      enrichKEGG_down_bar <- barplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50
      )
      ggsave(
        enrichKEGG_down_dot,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down dotplot.pdf")),
        height = 6,
        width = 8
      )
      ggsave(
        enrichKEGG_down_bar,
        filename = file.path(diff_outdir, str_glue("{grp[1]} vs {grp[2]} KEGG down barplot.pdf")),
        height = 6,
        width = 8
      )
    } else {
      cat(paste0(grp[1], " vs ", grp[2], " has No enriched KEGG down terms found. Skipping plot generation.\n"))
    }
  })






readme <- "DG diff analysis: 
B vs D has No enriched KEGG down terms found. Skipping plot generation.
B vs H has No enriched KEGG up terms found. Skipping plot generation.
B vs K has No enriched KEGG up terms found. Skipping plot generation.
B vs K has No enriched KEGG down terms found. Skipping plot generation.
D vs B has No enriched KEGG up terms found. Skipping plot generation.
D vs H has No enriched KEGG up terms found. Skipping plot generation.
D vs H has No enriched KEGG down terms found. Skipping plot generation.
D vs M has No enriched KEGG up terms found. Skipping plot generation.
D vs M has No enriched KEGG down terms found. Skipping plot generation.
H vs B has No enriched KEGG down terms found. Skipping plot generation.
H vs D has No enriched KEGG up terms found. Skipping plot generation.
H vs D has No enriched KEGG down terms found. Skipping plot generation.
H vs K has No enriched KEGG up terms found. Skipping plot generation.
H vs M has No enriched KEGG up terms found. Skipping plot generation.
H vs M has No enriched KEGG down terms found. Skipping plot generation.
K vs B has No enriched KEGG up terms found. Skipping plot generation.
K vs B has No enriched KEGG down terms found. Skipping plot generation.
K vs H has No enriched KEGG down terms found. Skipping plot generation.
K vs M has No enriched KEGG up terms found. Skipping plot generation.
K vs M has No enriched KEGG down terms found. Skipping plot generation.
M vs D has No enriched KEGG up terms found. Skipping plot generation.
M vs D has No enriched KEGG down terms found. Skipping plot generation.
M vs H has No enriched KEGG up terms found. Skipping plot generation.
M vs H has No enriched KEGG down terms found. Skipping plot generation.
M vs K has No enriched KEGG up terms found. Skipping plot generation.
M vs K has No enriched KEGG down terms found. Skipping plot generation.

CA1 diff analysis: 
B vs K has No enriched KEGG down terms found. Skipping plot generation.
H vs K has No enriched KEGG down terms found. Skipping plot generation.
H vs M has No enriched KEGG down terms found. Skipping plot generation.
K vs B has No enriched KEGG up terms found. Skipping plot generation.
K vs H has No enriched KEGG up terms found. Skipping plot generation.
M vs H has No enriched KEGG up terms found. Skipping plot generation.
"
writeLines(readme, file.path(results_outdir, "readme.txt"))








# 260512 - metadata 导出 ---------
meta_data_outdir <- "20260512_annotation"
dir.create(meta_data_outdir)

annotation_metadata <- all_data[[]] %>% 
  rownames_to_column("cell_id") %>% 
  dplyr::select(cell_id, annotation, fov)
project_groups <- split(annotation_metadata, annotation_metadata$fov)


walk(names(project_groups), function(proj) {
  proj_data <- project_groups[[proj]] %>% 
    dplyr::select(cell_id, group = annotation) %>% 
    mutate(cell_id = str_remove(cell_id, "^\\d+_"))
  
  filename <- paste0("data/", proj, "_cell_annotation.csv")
  
  write.csv(proj_data, file = filename, row.names = FALSE)
})











# 260530 gene plot ----
library(Seurat)
library(tidyverse)
library(ggplot2)

results_outdir <- "20260530_results"
dir.create(results_outdir)

all_data <- readRDS("data_for_analysis/all_data_h5ad_ROI_0107.rds")
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
cells_boundary_ls_filt <- cells_boundary_ls[!names(cells_boundary_ls) %in% c("B1", "B4")]
names(cells_boundary_ls_filt)

# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")


sp_genes <- c("Trem2")
sp_p_ls <- 
  # cells_boundary_ls %>% 
  cells_boundary_ls[!names(cells_boundary_ls) %in% c("B1", "B4")] %>% 
  # .[1] %>% 
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    dat <- 
      cells_boundary %>% 
      right_join(dat_expr, by = "cell_id")
    
    p_ls <- 
      sp_genes %>% 
      set_names(., .) %>% 
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        p <- dat %>% 
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(aes(fill = .data[[gene]], group = cell_id)) + 
          scale_y_reverse() + 
          coord_equal() + 
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") + 
          theme_void() + 
          labs(
            title = project
          )
        p
      })
    p_ls
  })

# export
pdf(file.path(results_outdir, str_glue("spatial plot - {sp_genes}.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()





### 气泡图 ------------
pdf(file.path(results_outdir, str_glue("Dot plot by group - {sp_genes}.pdf")), height = 6, width = 5)
DotPlot(all_data, features = sp_genes, dot.scale = 10, group.by = "group")
dev.off()


pdf(file.path(results_outdir, str_glue("Dot plot by annotation - {sp_genes}.pdf")), height = 6, width = 5)
DotPlot(all_data, features = sp_genes, dot.scale = 10, group.by = "annotation")
dev.off()




plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

### 小提琴图 ------------
pdf(file.path(results_outdir, str_glue("Vlnplot by group - {sp_genes}.pdf")), height = 5, width = 10)
VlnPlot(all_data, features = sp_genes, group.by = "group", pt.size = 0, cols = plan)
dev.off()


pdf(file.path(results_outdir, str_glue("Vlnplot by annotation - {sp_genes}.pdf")), height = 5, width = 12)
VlnPlot(all_data, features = sp_genes, group.by = "annotation", pt.size = 0, cols = plan)
dev.off()






data <- all_data[[]] %>% rownames_to_column("Cell")
dat_filt <- data[, c("Cell", "sample", "group", "annotation")]





# 260713 gene plot ----
library(Seurat)
library(tidyverse)
library(ggplot2)

results_outdir <- "20260713_results"
dir.create(results_outdir)

all_data <- readRDS("data_for_analysis/all_data_h5ad_ROI_0107.rds")
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
cells_boundary_ls_filt <- cells_boundary_ls[!names(cells_boundary_ls) %in% c("B1", "B4")]
names(cells_boundary_ls_filt)

# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")


sp_genes <- c("Nos2", "Cd86", "Tnf", "Map2", "Gap43", "Arg1", "Mrc1", "Cd163", "Lpl", "Cd68", "Opn3", "Opn4")
sp_p_ls <- 
  # cells_boundary_ls %>% 
  cells_boundary_ls_filt %>%
  # cells_boundary_ls[!names(cells_boundary_ls) %in% c("B1", "B4")] %>% 
  # .[1] %>% 
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    dat <- 
      cells_boundary %>% 
      right_join(dat_expr, by = "cell_id")
    
    p_ls <- 
      sp_genes %>% 
      set_names(., .) %>% 
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        p <- dat %>% 
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(aes(fill = .data[[gene]], group = cell_id)) + 
          scale_y_reverse() + 
          coord_equal() + 
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") + 
          theme_void() + 
          labs(
            title = project
          )
        p
      })
    p_ls
  })

# export
pdf(file.path(results_outdir, str_glue("spatial plot - sp genes.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()


for (gene in sp_genes) {
  pdf(file.path(results_outdir, str_glue("spatial plot - {gene}.pdf")), height = 18, width = 18)
  
  cells_boundary_ls_filt %>%
    imap(function(cells_boundary, project){
      dat_expr <- FetchData(all_data, vars = gene) %>% rownames_to_column("cell_id")
      dat <- cells_boundary %>% right_join(dat_expr, by = "cell_id")
      
      # q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
      q_val_rng <- quantile(dat[[gene]], c(0, 1))
      dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
      
      p <-
        dat %>%
        ggplot(aes(x = vertex_x, y = vertex_y)) +
        geom_polygon(aes(fill = .data[[gene]], group = cell_id)) +
        scale_y_reverse() +
        coord_equal() +
        scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
        theme_void() +
        labs(title = project)
      
      print(p)
    })
  
  dev.off()
}







### 气泡图 ------------
pdf(file.path(results_outdir, str_glue("Dot plot by group.pdf")), height = 6, width = 8)
DotPlot(all_data, features = sp_genes, dot.scale = 10, group.by = "group", cols = c("lightgrey", "red")) +
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))
dev.off()


pdf(file.path(results_outdir, str_glue("Dot plot by annotation.pdf")), height = 6, width = 8)
DotPlot(all_data, features = sp_genes, dot.scale = 10, group.by = "annotation", cols = c("lightgrey", "red")) +
  labs(x = "", y = "") + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))
dev.off()




plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

### 小提琴图 ------------
pdf(file.path(results_outdir, str_glue("Vlnplot by group.pdf")), height = 8, width = 6)
VlnPlot(all_data, 
        features = sp_genes, 
        group.by = "group", 
        pt.size = 0, 
        cols = plan,
        stack = T,
        flip = T) +
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 0, vjust = 0, hjust = 0.5),
        strip.text.y = element_text(hjust = 0))
dev.off()


pdf(file.path(results_outdir, str_glue("Vlnplot by annotation.pdf")), height = 8, width = 8)
VlnPlot(all_data, 
        features = sp_genes, 
        group.by = "annotation", 
        pt.size = 0, 
        cols = plan,
        stack = T,
        flip = T) +
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1),
        strip.text.y = element_text(hjust = 0))
dev.off()



