# 260710 spatial plot ----------
library(Seurat)
library(tidyverse)
library(patchwork)
# library(sf)
library(ggplot2)
library(harmony)

all_data <- readRDS("data_for_analysis/all_data_h5ad.rds")
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")


### spaial plot of some gene ------------
sp_gene_output <- "260710-spatial gene plot"
dir.create(sp_gene_output)

all_sp_genes <- c(
  "ACAN", "ACKR1", "ACTA2", "APOE", "C3", "C7", "CADM1", "CALML5", "CCL2", "CCL5",
  "CD3D", "CD69", "CD74", "CDH19", "CHI3L1", "CLDND1", "CNN1", "COL1A1", "COL6A5", "CPA3",
  "CRABP1", "CSMD1", "CSTA", "CXCL1", "CYP26B1", "DCD", "DGKZ", "DNAJB1", "DUSP2", "DUSP4",
  "FBN1", "FFAR2", "FN1", "FTH1", "FTL", "GNLY", "GPR183", "GTF2IRD1", "HIGD1B", "HSPA1A",
  "HSPA1B", "ID3", "IGHD", "IL32", "IL4I1", "IL7R", "ISG15", "ITGA10", "ITGB1", "JCHAIN",
  "KLF2", "KLK8", "LAYN", "LTB", "LY6G6C", "LYVE1", "LYZ", "MGST1", "MMP1", "MT1X",
  "MUCL1", "NFKBIZ", "NKG7", "NR4A2", "PCLAF", "PCOLCE2", "PHLDA1", "PLA2G2A", "PMEL", "PROK2",
  "PTGDS", "PTPRC", "RGS13", "RORA", "S100A11", "S100A12", "S100A2", "S100A8", "SAT1", "SLC4A10",
  "SLC7A1", "SOX4", "SPARC", "SPINK5", "SRGN", "TFF3", "THBS4", "TIMP1", "TIMP3", "TNN",
  "TOB1", "TRAC", "TRIM33", "TXNIP", "TYRP1", "VCAN", "VWF", "WIF1", "XCL1", "XIST"
)


for (sp_gene in all_sp_genes) {
  sp_p_ls <-
    cells_boundary_ls_raw %>%
    imap(function(cells_boundary, project) {
      dat_expr <- FetchData(all_data, vars = sp_gene) %>% rownames_to_column("cell_id")

      dat <-
        cells_boundary %>%
        right_join(dat_expr, by = "cell_id")

      p_ls <-
        sp_gene %>%
        set_names(., .) %>%
        map(function(gene) {
          # q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
          q_val_rng <- quantile(dat[[gene]], c(0, 1))
          
          dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])

          p <- dat %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = .data[[gene]], group = cell_id)) +
            scale_y_reverse() +
            coord_equal() +
            scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
            theme_void() +
            labs(title = project)
          p
        })
      p_ls
    })

  pdf(file.path(sp_gene_output, str_glue("spatial plot - {sp_gene}.pdf")), height = 18, width = 18)
  print(sp_p_ls)
  dev.off()
}


### 气泡图 ------------
# theme_90_angle <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))

pdf(file.path(sp_gene_output, str_glue("Dot plot by sample - sp_gene.pdf")), height = 15, width = 25)
DotPlot(all_data, features = all_sp_genes, dot.scale = 10, cols = c("lightgrey", "red"), group.by = "sample") +
  labs(x = "", y = "") +
  theme_bw() + # 改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)
  ) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) # 边框粗细
dev.off()

pdf(file.path(sp_gene_output, str_glue("Dot plot by group - sp_gene.pdf")), height = 5, width = 25)
DotPlot(all_data, features = all_sp_genes, dot.scale = 10, cols = c("lightgrey", "red"), group.by = "group") +
  labs(x = "", y = "") +
  theme_bw() + # 改变ggplot2的主题
  # coord_flip() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)
  ) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) # 边框粗细
dev.off()



pdf(file.path(sp_gene_output, str_glue("Dot plot by annotation - sp_gene.pdf")), height = 6, width = 25)
DotPlot(all_data, features = all_sp_genes, dot.scale = 10, cols = c("lightgrey", "red"), group.by = "annotation") +
  labs(x = "", y = "") +
  theme_bw() + # 改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)
  ) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) # 边框粗细
dev.off()




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
plan <- c(
  "#EA7070", "#FDC4B6", "#E59572", "#2694AB", "#96CEB4", "#600473",
  "#7F95D1", "#DEE2D1", "#EF5285", "#A5DFF9", "#E8A0B8", "#FFEEAD",
  "#737495", "#BCCF3D", "#DE4307", "#F29C2B", "#F6D04D", "#8BC24C",
  "#FED95C", "#D9534F", "#FFAD60", "#57D1C9", "#ED5485", "#FFE869",
  "#FA6E57", "#F69E53", "#248888", "#999999", "#E7475E", "#F0D879",
  "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#F1AC9D", "#F06966",
  "#6ABE83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", "#60C5BA",
  "#FEEE7D", "#FFC300", "#02C9C9", "#FBF579", "#4695D6", "#005995",
  "#FA625F", "#F17D80", "#FFB3A1", "#C6E2B6", "#A3D9D7", "#F7D794",
  "#DA7297", "#89C2D9", "#F4A261", "#B0C4DE", "#E76F51", "#2A9D8F",
  "#E9C46A", "#264653", "#F4A492", "#A8DADC", "#457B9D", "#F08080"
)

### 小提琴图 ------------
pdf(file.path(sp_gene_output, str_glue("Vlnplot by group - sample.pdf")), height = 20, width = 10)
VlnPlot(all_data,
  pt.size  = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
  stack    = T, # 是否堆叠
  flip     = T, # 翻转
  # adjust = 5, # 调整粗细
  fill.by  = "ident",
  cols     = plan,
  group.by = "sample",
  features = all_sp_genes
) +
  labs(x = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, color = "black", angle = 45, vjust = 1),
    strip.text.y = element_text(hjust = 0)
  )
dev.off()

pdf(file.path(sp_gene_output, str_glue("Vlnplot by group - group.pdf")), height = 20, width = 5)
VlnPlot(all_data,
        pt.size  = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
        stack    = T, # 是否堆叠
        flip     = T, # 翻转
        # adjust = 5, # 调整粗细
        fill.by  = "ident",
        cols     = plan,
        group.by = "group",
        features = all_sp_genes
) +
  labs(x = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, color = "black", angle = 45, vjust = 1),
    strip.text.y = element_text(hjust = 0)
  )
dev.off()


pdf(file.path(sp_gene_output, str_glue("Vlnplot by annotation - sp_gene.pdf")), height = 20, width = 8)
VlnPlot(all_data,
        pt.size  = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
        stack    = T, # 是否堆叠
        flip     = T, # 翻转
        # adjust = 5, # 调整粗细
        fill.by  = "ident",
        cols     = plan,
        group.by = "annotation",
        features = all_sp_genes
) +
  labs(x = "") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1),
    strip.text.y = element_text(hjust = 0)
  )
dev.off()





