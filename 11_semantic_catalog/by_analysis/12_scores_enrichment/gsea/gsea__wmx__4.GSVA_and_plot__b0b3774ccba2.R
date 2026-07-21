# 读取下列的rds，进行整理，Cluster列为paste0('Cluster',leiden列)
# 对其进行基因list的打分
# 对GSVA的结果进行展示 UMAP和dotplot


library(Seurat)
library(dplyr)
library(ggplot2)
library(GSVA,lib.loc = '/mnt/NexData_file/share/scRNA-tools/envs/pip_GSVA/lib/R/library/')
library(GSEABase)
library(patchwork)
library(viridis)

geneset <- list(
  N1 = c(
    "Ifnb1","Il1b","Tnf","Cxcl9","Cxcl10","Cxcl11",
    "Ccl3","Ccl4","Ccl5","Nos2","Icam1","Stat1",
    "Isg15","Irf1","Irf7","Socs1"
  ),
  N2 = c(
    "Arg1","Mmp9","Vegfa","Tgfb1","Il10","Cxcl1",
    "Cxcl2","Cxcl5","Ccl2","Ccl17","S100a8","S100a9",
    "Osm","Chil3","Lgals3","Ptgs2","Hp"
  )
)

rds_files <- list.files(
  path = "/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou6/GSVA/",
  pattern = "count.rds$",
  full.names = TRUE
)

outdir <- "/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou6/GSVA/GSVA_score_plot"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (rds in rds_files) {
  
  message("========== Processing: ", rds)
  
  obj <- readRDS(rds)
  DefaultAssay(obj) <- "RNA"
  
  sample_name <- gsub("_count$", "", tools::file_path_sans_ext(basename(rds)))
  obj_outdir <- file.path(outdir, sample_name)
  dir.create(obj_outdir, showWarnings = FALSE, recursive = TRUE)
  
  obj$Cluster <- paste0("Cluster", obj$leiden)
  obj$Cluster <- factor(obj$Cluster, levels = paste0("Cluster", sort(unique(obj$leiden))))
  
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  expr <- as.matrix(expr)
  
  geneset_use <- lapply(geneset, function(x) intersect(x, rownames(expr)))
  geneset_use <- geneset_use[lengths(geneset_use) > 0]
  
  print(lapply(geneset_use, length))
  
  param <- GSVA::gsvaParam(
    exprData = expr,
    geneSets = geneset_use,
    kcdf = "Gaussian"
  )
  
  gsva_res <- GSVA::gsva(param)
  gsva_res <- t(gsva_res)
  colnames(gsva_res) <- paste0(colnames(gsva_res), "_GSVA")
  
  obj@meta.data <- cbind(obj@meta.data, gsva_res[rownames(obj@meta.data), , drop = FALSE])
  
  score_cols <- colnames(gsva_res)
  
  write.table(
    obj@meta.data,
    file = file.path(obj_outdir, paste0(sample_name, "_metadata_with_GSVA_score.xls")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  # write.csv(
  #   gsva_res,
  #   file = file.path(obj_outdir, paste0(sample_name, "_GSVA_score_matrix.csv"))
  # )
  
  if ("umap" %in% names(obj@reductions)) {
    
    p_umap <- FeaturePlot(
      obj,
      features = score_cols,
      reduction = "umap",
      cols = c("#F2F2F2", "#2166AC", "#B2182B"),
      ncol = 2,
      pt.size = 0.35
    ) &
      theme_classic(base_size = 14) &
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    
    ggsave(file.path(obj_outdir, paste0(sample_name, "_GSVA_umap.pdf")),
           p_umap, width = 10, height = 5)
    ggsave(file.path(obj_outdir, paste0(sample_name, "_GSVA_umap.png")),
           p_umap, width = 10, height = 5, dpi = 300)
  }
  
  dot_df <- obj@meta.data %>%
    select(Cluster, all_of(score_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(score_cols),
      names_to = "Signature",
      values_to = "Score"
    ) %>%
    group_by(Cluster, Signature) %>%
    summarise(
      avg_score = mean(Score, na.rm = TRUE),
      .groups = "drop"
    )
  
    p_dot <- ggplot(dot_df, aes(x = Signature, y = Cluster)) +
      geom_point(aes(color = avg_score), size = 6) +
      scale_color_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0
      ) +
      theme_bw(base_size = 14) +
      labs(
        x = NULL,
        y = NULL,
        color = "Mean GSVA\nscore"
    )+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle(paste0(sample_name, " GSVA score"))
  
  ggsave(file.path(obj_outdir, paste0(sample_name, "_GSVA_dotplot.pdf")),
         p_dot, width = 6, height = 5)
  ggsave(file.path(obj_outdir, paste0(sample_name, "_GSVA_dotplot.png")),
         p_dot, width = 6, height = 5, dpi = 300)
  
  saveRDS(obj, file.path(obj_outdir, paste0(sample_name, "_GSVA_score.rds")))
  
  message("[DONE] ", sample_name)
}
