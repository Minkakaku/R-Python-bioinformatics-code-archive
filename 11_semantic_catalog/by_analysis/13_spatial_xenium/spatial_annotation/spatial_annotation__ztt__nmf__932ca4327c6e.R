library(Seurat)
library(tidyverse)


# all_data <- readRDS('all_data.rds') 
# all_data <- readRDS('all_data_fixed2.rds')
all_data <- readRDS('all_data_fixed.rds')

sample_match_table <-
  data.frame(
    project = all_data$project %>% unique
  ) %>%
  mutate(
    Sample = str_c("Sample", as.numeric(factor(project)))
  )

sample_match_vars  <- structure(sample_match_table$project, names = sample_match_table$Sample)
sample_match_vars2 <- structure(names(sample_match_vars), names = sample_match_vars)



cell_freq <-
  all_data[[]] %>%
  
  # mutate(
  #   # 使用final_annotation的注释结果  
  #   final_annotation_v1 = final_annotation,
  #   final_annotation    = final_annotation2
  # )  %>%
  
  mutate(major_annotation2 = major_annotation) %>%
  mutate(major_annotation2 = fct_recode(
    major_annotation2,
    'Myeloid cells' = 'Granulocytes', 
    'T/NK cells' = 'B cells',
    Stroma = 'Fibroblasts',
    Stroma = 'Endothelial cells'
  )) %>%
  
  group_by(major_annotation2, final_annotation, project, Group) %>%
  summarise(n = n()) %>%
  group_by(major_annotation2, project, Group) %>%
  mutate(n_total = sum(n)) %>%
  mutate(freq = n/n_total)

cell_freq %>%
  
  ggplot(aes(x = project, y = freq, fill = final_annotation)) +
  geom_col()


## NMF -------
library(NMF)
library(ComplexHeatmap)

cell_freq_wide <- 
  cell_freq %>%
  ungroup() %>%
  dplyr::select(-n, -n_total, -major_annotation2) %>%
  pivot_wider(
    names_from  = "final_annotation",
    values_from = "freq"
  ) 

have_na <- apply(cell_freq_wide, 2, function(x) any(is.na(x)))
if(sum(have_na) > 0){
  cell_freq_wide <- cell_freq_wide[, !have_na]
}


# NMF analysis ------

nrank=5
mtx = cell_freq_wide[, 3:ncol(cell_freq_wide)]
nmf.res = nmf(t(mtx), rank=nrank, nrun=20, seed=123456)

# saveRDS(nmf.res, 'nmf.res.rds')

# W------
w = basis(nmf.res)
rownames(w) = colnames(mtx)
colnames(w) = paste0('NMF', seq(1,nrank))

w_ = w %>% t %>% scale
set.seed(123456)
p=ComplexHeatmap::Heatmap(
  w_, 
  name         = " ",
  width        = 11, 
  height       = 4,
  column_km    = nrank, 
  cluster_rows = F, 
  clustering_method_columns = 'single',
  col = circlize::colorRamp2(c(-2, 0, 2), c("#234B88", "white", "#8A2D58"))
) 
heatmap = ComplexHeatmap::draw(p)


# or
pheatmap::pheatmap(
  t(w),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  scale = "row",
  clustering_method = "single",
  border_color = NA,
  cluster_rows = FALSE,
  # cluster_cols = FALSE
)



# export
pdf(
  "nmf heatmap - NMF basis - 2.pdf",
  height = 4, 
  width = 10
)
print(heatmap)
dev.off()



# H -----------
h = coef(nmf.res)
colnames(h) = cell_freq_wide$project
rownames(h) = paste0('NMF', seq(1,nrank))

cell_freq_wide %>%
  dplyr::select(project, Group) %>%
  group_by(Group) %>%
  nest() %>%
  deframe() %>%
  map(unlist) %>%
  imap(function(sample, name){
    pheatmap::pheatmap(
      h[, sample],
      color  = colorRampPalette(c("navy", "white", "firebrick3"))(50)
    )
  })

# read into new gorup
sample_info <- readxl::read_excel("sample_info update 0903.xlsx")


order_samples <- 
  # cell_freq_wide %>% 
  # dplyr::select(project, Group) %>%
  # arrange(Group) %>%
  # column_to_rownames("project") 
  sample_info %>%
  dplyr::select(project, Group) %>%
  arrange(Group) %>%
  column_to_rownames("project")
  
pheatmap::pheatmap(
  h[, rownames(order_samples)],
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  annotation_col = order_samples,
  scale = "row",
  border_color = NA,
  labels_col = sample_match_vars2[colnames(h)],
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

h_ = h %>% t %>% scale %>% t
set.seed(123456)
p2=ComplexHeatmap::Heatmap(
  h_[, rownames(order_samples)], 
  col = circlize::colorRamp2(c(-2, 0, 2), c("#234B88", "white", "#8A2D58")),
  width = 11,
  height = 4,
  column_labels = sample_match_vars2[colnames(h)],
  # column_km = nrank, 
  top_annotation = HeatmapAnnotation(
    Group =  order_samples$Group, 
    col = list(Group =c(
                   P   = "#71C3A3",
                   WT  = "#87CEEB",
                   MUT = "#F8766C"
                 ))
  ),
  cluster_rows = F, 
  cluster_columns = T, 
  column_split  = order_samples$Group,
  # clustering_method_columns = 'single'
) 
heatmap2 = ComplexHeatmap::draw(p2)


# export
pdf(
  "nmf heatmap - NMF coef - 2.pdf",
  height = 4, 
  width = 10
)
print(heatmap2)
dev.off()









