pacman::p_load(plotthis,qs,tidyverse,FNN,dbscan,openxlsx,ggpubr)

source("/mnt/NexData_file/home/lyh/Script/Xenium/CNs_ana/CNs_V2/all_CNs.R")
source("/mnt/NexData_file/home/lyh/Script/R/R_fuc/color.R")

## ----------------------------
## 0. User settings
## ----------------------------
meta_path <- "/mnt/NexData_file/home/lyh/Xenium_aftersale/Sunhandong/20260402_Roe/meta.qs"
boundary_path <- "/mnt/NexData_file/home/hyj/Rstudio/XE2510316144-02_human_breast/R_OUTPUT/data_for_analysis/cells_boundary_ls.rds"
output_dir <- "/mnt/NexData_file/home/lyh/Xenium_aftersale/Sunhandong/20260408_CNs/"
output_csv <- file.path(output_dir, "add_CNs.csv")
region_col <- "sample"
cluster_col <- "first_type"
x_col <- "x_centroid"
y_col <- "y_centroid"
group_col <- "group"
neighbor_mode <- "knn"   # choose "knn" or "radius"
k_neighbors <- 20
radius_value <- 40
n_clusters <- 8
random_seed <- 0
kmeans_nstart <- 20
kmeans_iter_max <- 100
write_row_names <- TRUE

## ----------------------------
## 3. Read input metadata and setwd
## ----------------------------
#setwd
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)
getwd()
if (!file.exists(meta_path)) {
  stop("meta_path does not exist: ", meta_path)
}
ext <- tolower(tools::file_ext(meta_path))
if (ext == "qs") {
  meta_df<- qs::qread(meta_path)
} else if (ext == "csv") {
 meta_df <- read.csv(meta_path, row.names = 1, check.names = FALSE)
} else if (ext == "rds") {
  meta_df <- readRDS(meta_path)
} else if (ext %in% c("xlsx", "xls")) {
  meta_df <- openxlsx::read.xlsx(meta_path, sheet = 1)
} else {
  stop("不支持的文件格式: ", ext)
}
boundary_df <- readRDS(boundary_path)

## ----------------------------
## 4. Run neighborhood analysis
## ----------------------------
cn_res <- run_cn_analysis(
  df = meta_df,
  region_col = region_col,
  cluster_col = cluster_col,
  x_col = x_col,
  y_col = y_col,
  neighbor_mode = neighbor_mode,
  k_neighbors = k_neighbors,
  radius_value = radius_value,
  n_clusters = n_clusters,
  random_seed = random_seed,
  kmeans_nstart = kmeans_nstart,
  kmeans_iter_max = kmeans_iter_max
)

result_df <- cn_res$result_df
window_df <- cn_res$window_df
feature_matrix <- cn_res$feature_matrix
kmeans_res <- cn_res$kmeans_res
cn_col <- cn_res$cn_col
cn_col
table(result_df[[cn_col]])
result_df[,cn_col] <- paste0("CNs_",result_df[,cn_col])
head(result_df[, c(region_col, cluster_col, x_col, y_col, cn_col)])

#save result
write.csv(result_df  %>% select(cn_col),output_csv)

#
dir.create("plot")

#plot heatmap
p <- cn_exp_heatmap(result_df,cluster_col=cluster_col,cn_col=cn_col)
plot_df <- p$data
n_row <- if (is.factor(plot_df$CN)) nlevels(plot_df$CN) else dplyr::n_distinct(plot_df$CN)
n_col <- if (is.factor(plot_df$CellType)) nlevels(plot_df$CellType) else dplyr::n_distinct(plot_df$CellType)
plot_h <- max(2.5, 8 * n_row * 0.8 / n_col + 1.2)
plot_w <- 8
options(repr.plot.width=plot_w , repr.plot.height=plot_h ) 
p
ggsave("plot/CNs_heatmap.pdf", p, width = plot_w, height = plot_h)
system(paste("convert -density 300 ","plot/CNs_heatmap.pdf","plot/CNs_heatmap.png"))

#set colors
unique_celltype <- unique(result_df[,cluster_col])
unique_cns <- unique(result_df[,cn_col])
color_celltype <- setNames(color6[seq_len(length(unique_celltype))],nm=unique_celltype)                    
colors_cns <- setNames(color3[seq_len(length(unique_cns))],nm=unique_cns)

Sum <- table(result_df[[group_col]], result_df[[cn_col]])
Sum_pct <- as.data.frame(t(Sum)) %>% group_by(Var2) %>% mutate( Freq_pct =Freq/sum(Freq)*100) %>% dplyr::select(-"Freq") 
Sum_pct_1 <- Sum_pct %>% pivot_wider(names_from = "Var2", values_from = "Freq_pct") %>% column_to_rownames("Var1") %>% t() %>% as.data.frame()
write.xlsx(Sum,"plot/CNs_sum_group.xlsx",rowNames=T)
write.xlsx(Sum_pct_1,"plot/CNs_sum_group_pct.xlsx",rowNames=T)
####plot
p <- plotthis::TrendPlot(Sum_pct, x = "Var2", y = "Freq_pct", group_by = "Var1",scale_y = TRUE)+theme_bw()+
  scale_fill_manual(name="CNs",values = colors_cns)  +theme(legend.position = "right",  
                     panel.grid = element_blank(), 
                     axis.text.y = element_text(color = "black",),  # Bold for y-axis text
                     axis.text.x = element_text(color = "black",angle = 45,hjust = 1))+labs(x="",y="Percentage",title = "")
p
n_row <- length(unique(result_df[[group_col]]))
options(repr.plot.width=3*n_row, repr.plot.height=6) 
ggsave("plot/CNs_Stacked_group.pdf",p,height = 6,width = 3*n_row)
ggsave("plot/CNs_Stacked_group.png",p,height = 6,width =3*n_row)

#count sample
plot.order <- paste0("CNs_",seq_len(length(unique_cns)))
metadata.df <- result_df %>% tibble::rownames_to_column(var = "Cell")
cell_num.df <- metadata.df %>% dplyr::count(!!sym(region_col),!!sym(group_col), name = "Num")
ratio.df <- metadata.df %>% dplyr::count(!!sym(region_col),!!sym(group_col),!!sym(cn_col), name = "Sum") %>%
  tidyr::complete(nesting(!!sym(region_col),!!sym(group_col)),!!sym(cn_col), fill = list(Sum = 0)) %>%
  left_join(cell_num.df, by = c(region_col,group_col)) %>%
  dplyr::mutate(Percentage = Sum / Num)
# Create wide format data
wide.df <- ratio.df %>% dplyr::select(!!sym(region_col),!!sym(group_col),!!sym(cn_col),Percentage) %>% pivot_wider(names_from = !!sym(cn_col), values_from = Percentage)
write.xlsx(wide.df,"plot/CNs_Proportion_sample.xlsx")

p <- ggplot(ratio.df,aes(x =!!sym(cn_col), y = Percentage, fill =!!sym(group_col))) +
    geom_boxplot(outlier.shape = NA, width = 0.25, position=position_dodge(width=0.5), size = 0.1, colour = "#FFFFFF") +
    geom_point(pch = 21, position = position_jitterdodge(), colour = "white", size = 0.5) +
    theme_classic() +
    scale_fill_manual(values = color4) +
    stat_compare_means(aes(group =!!sym(group_col)), method = "t.test", 
                        symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")), 
                        label = "p.signif",size = 3) +
    labs(title = "Comparison of CNs Proportion", x = "CNs", y = "CNs Proportion") + theme_classic() +    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.title = element_text(size=10),legend.key.size = unit(0.4,"cm"))
ggsave("plot/Cell_Proportion_BoxPlot_sign.pdf",height = 5,width =4*n_row)
ggsave("plot/Cell_Proportion_BoxPlot_sign.png",height = 5,width =4*n_row,dpi=300)

#spatial plot
dir.create("plot/spatial_plot")
join_df <- result_df[,c(region_col,cluster_col,cn_col)] %>%  rownames_to_column("cell_id") 
join_df$cell_id <- gsub("^.*_","",join_df$cell_id)
parallel::mclapply(names(boundary_df), function(i) {
  tmp_df <- boundary_df[[i]]
  tmp_df$cell_id <- gsub("^.*_", "", tmp_df$cell_id)
  tmp_df <- left_join(tmp_df,join_df %>% filter(!!sym(region_col) == i),by = c("cell_id" = "cell_id"))

  p <- ggplot(tmp_df, aes(x = vertex_x, y = vertex_y)) +
    geom_polygon(
      aes(fill = CN_k20_n8, color = first_type, group = cell_id),
      linewidth = 0.1
    ) +
    scale_fill_manual(name = "CNs", values = colors_cns) +
    scale_color_manual(name = "Cell type", values = color_celltype) +
    guides(
      color = guide_legend(
        override.aes = list(fill = NA, linewidth = 0.5)
      )
    ) +
    coord_equal() +
    scale_y_reverse() +
    theme_void()

  pdf_name <- file.path("plot/spatial_plot", paste0(i, "_spatial.pdf"))
  png_name <- file.path("plot/spatial_plot", paste0(i, "_spatial.png"))

  ggsave(
    filename = pdf_name,
    plot = p,
    width = 10,
    height = 10,
    units = "in",
    limitsize = FALSE
  )

  cmd <- sprintf('convert -density 300 "%s" "%s"', pdf_name, png_name)
  system(cmd)

  return(i)
}, mc.cores = 8)
