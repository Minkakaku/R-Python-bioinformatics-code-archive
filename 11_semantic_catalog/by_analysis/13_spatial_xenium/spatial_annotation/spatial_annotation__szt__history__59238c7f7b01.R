net_ls %>%
iwalk(function(net, name){
all_nodes <- attr(V(net), "names")
# export pdf
pdf(
file.path(cpdb_plot_out_dir, str_glue("{name} network plot - single node.pdf")),
height = 10,
width  = 10
)
all_nodes %>%
walk(function(x){ # browser()
coords        <- layout_in_circle(net)  # 设置网络布局
# 调整节点位置的线条角度
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
angles <- atan(coords[igraph::V(net),2]/coords[igraph::V(net),1])
loop.angle <- ifelse(coords[igraph::V(net), 1]>0, -angles, pi-angles)
loop.edge.idx <- which(edge.start[,2] == edge.start[,1])
igraph::E(net)$loop.angle <- 0
igraph::E(net)$loop.angle[loop.edge.idx] <- loop.angle[edge.start[loop.edge.idx, 1]]
# 只保留当前细胞发出的线
E(net)$color[! str_starts(attr(E(net), 'vnames' ), fixed(x))] <- NA
E(net)$label[! str_starts(attr(E(net), 'vnames' ), fixed(x))] <- NA
plot(
net,
edge.curved=0.3,
edge.arrow.size=0.3,
# label.cex = 15,
layout = coords,
margin = rep(0.2, 4),
main = name
)
})
dev.off()
})
## heatmap --------
p_cpdb_heatmap_ls2 <-
network_dat_ls %>%
imap(function(network_dat, name){
#browser()
dat <- network_dat %>%
pivot_wider(
names_from  = "TARGET",
values_from = "counts",
values_fill = 0
) %>%
column_to_rownames("SOURCE")
# pheatmap::pheatmap(
#   dat,
#   scale = "none",
#   color = colorRampPalette(c("#124E8B", "#FDD8B8", "#8C0C51"))(50),
#   method = "ward.D2",
#   border_color = NA,
#   treeheight_row = 0,
#   treeheight_col = 0,
#   # main = name,
#   silent = F,
#   na_col = "white"
# )
cluster_row <- hclust(dist(dat), method = "ward.D2")
cluster_col <- hclust(dist(t(dat)), method = "ward.D2")
dat %>%
rownames_to_column("SOURCE") %>%
pivot_longer(
-SOURCE,
names_to = "TARGET",
values_to = "counts"
) %>%
mutate(
SOURCE = factor(SOURCE, levels = rev(cluster_row$labels[cluster_row$order])),
TARGET = factor(TARGET, levels = cluster_col$labels[cluster_col$order]),
) %>%
ggplot(aes(x = TARGET, y = SOURCE, fill = counts)) +
geom_tile() +
scale_fill_gradientn(colours = c("#124E8B", "#FDD8B8", "#8C0C51")) +
theme_minimal() +
coord_fixed() +
ggtitle(name) +
theme(
# axis.title = element_text(),
axis.ticks.length.x = unit(0, 'mm'),
axis.text.x = element_text(
angle = -90,
hjust = 0,
vjust = 0.5,
size = 15
),
axis.text.y = element_text(size = 15),
axis.title  = element_text(size = 20)
)
})
p_cpdb_heatmap_combn_ls2 <-
p_cpdb_heatmap_ls2 %>%
enframe(name = "group", value = "p") %>%
separate(group, into = c("Group", "group"), sep = "_") %>%
dplyr::select(-group) %>%
group_by(Group) %>%
nest() %>%
deframe() %>%
map(1) %>%
map(function(x){ # browser()
rng_ls <- x %>% map('data') %>% map('counts') %>% map(range)
min_val <- rng_ls %>% map_dbl(1) %>% min()
max_val <- rng_ls %>% map_dbl(2) %>% max()
x %>%
map(~ .x + scale_fill_gradientn(colours = c("#124E8B", "#FDD8B8", "#8C0C51"), limits = c(min_val, max_val))) %>%
patchwork::wrap_plots()
})
pdf(
file.path(cpdb_plot_out_dir, str_glue("heatmap plot.pdf")),
height = 12,
width  = 12
)
p_cpdb_heatmap_combn_ls2 %>%
walk(function(p){
cowplot::plot_grid(p) %>% print()
})
dev.off()
# 解析 mean文件和pvalues文件
inter_df_mean_ls <-
all_files_ls %>%
map(function(all_files){
all_files  %>%
str_subset('_analysis_means_') %>%
read_tsv() %>%
dplyr::select(interacting_pair, where(is.numeric)) %>%
pivot_longer(
-interacting_pair,
names_to = 'cell_pair',
values_to = 'mean'
) %>%
dplyr::filter(!is.na(mean)) %>%
dplyr::filter(!cell_pair == 'rank') %>%
mutate(help_col = str_c(interacting_pair, cell_pair)) %>%
dplyr::select(-c(1, 2))
})
inter_df_p_ls <-
all_files_ls %>%
map(function(all_files){ # browser()
all_files  %>%
str_subset('pvalues') %>%
read_tsv() %>%
dplyr::select(interacting_pair, where(is.numeric)) %>%
pivot_longer(
-interacting_pair,
names_to = 'cell_pair',
values_to = 'p'
) %>%
# mutate(p = p.adjust(p)) %>%
dplyr::filter(p < 0.05) %>% # <------ p值: 没有筛选的作用
mutate(help_col = str_c(interacting_pair, cell_pair))
})
inter_df_inter_score_ls <-
all_files_ls %>%
map(function(all_files){ # browser()
all_files  %>%
str_subset('_interaction_scores_') %>%
read_tsv() %>%
dplyr::select(interacting_pair, where(is.numeric)) %>%
pivot_longer(
-interacting_pair,
names_to = 'cell_pair',
values_to = 'score'
) %>%
# mutate(p = p.adjust(p)) %>%
# dplyr::filter(p < 0.001) %>% # <------ p值: 没有筛选的作用
mutate(help_col = str_c(interacting_pair, cell_pair))
})
dot_dat_ls <-
inter_df_p_ls %>%
imap(function(inter_df_p, name){ # browser()
inter_df_mean <- inter_df_mean_ls[[name]]
inter_df_inter_score <- inter_df_inter_score_ls[[name]]
inter_df_p %>%
inner_join(inter_df_mean, by = 'help_col') %>%
inner_join(inter_df_inter_score[c('score', 'help_col')], by = 'help_col') %>%
mutate(p = str_replace(p,"^0$","1e-12")) %>%
group_by(interacting_pair) %>%
# dplyr::filter(p <= 0.05) %>%
arrange(p) %>%
mutate('-log10(pvalue)' = -log10(as.numeric(p))) %>%
mutate(mean = log2(mean))
})
# dotplot
p_dot_ls <-
dot_dat_ls %>%
imap(function(dot_dat, name){
dot_dat  %>%
as.data.frame() %>%
mutate(
cell_pair2 =
cell_pair %>%
str_split("\\|") %>%
map(~ str_remove_all(.x, "\\s.*$")) %>%
map_chr(str_c, collapse = "|")
) %>%
# slice_sample(prop = 0.1) %>%
ggplot(aes(
x = cell_pair2,
y = interacting_pair,
color = `-log10(pvalue)`
)) +
geom_point(shape = 22, aes(fill = mean, size = score), stroke = 2) +
theme_classic() +
labs(
x = NULL,
y = NULL,
title = name,
fill = "Log2 mean(Molecule 1, Molecule 2)",
# size = "-log10(p value)"
size = "Interaction scores"
) +
scale_fill_gradientn(colors = c("#fce5d7", "#f5ae8b",  "#d00f19" ,  "#6e020f")) +
scale_size(range = c(0.1, 5)) +
scale_color_gradient(low = "#386cb044", high = "#386cb0") +
theme(axis.text.x = element_text(
angle = 90,
hjust = 1,
vjust = 0.5
))
})
pdf(
file.path(cpdb_plot_out_dir, str_glue("dotplot.pdf")),
height = 30,
width  = 50
)
p_dot_ls %>% walk(print)
dev.off()
# dotplot - split ------
# ligand
p_dot_ls2 <-
dot_dat_ls %>%
enframe(name = "group") %>%
unnest(value) %>%
mutate(group_ligand = str_remove(cell_pair, "\\|.*$")) %>%
group_by(group_ligand) %>%
nest() %>%
deframe() %>%
imap(function(dat, name){ # browser()
dat %>%
as.data.frame() %>%
# mutate(
#   cell_pair2 =
#     cell_pair %>%
#     str_split("\\|") %>%
#     map(~ str_remove_all(.x, "\\s.*$")) %>%
#     map_chr(str_c, collapse = "|")
# ) %>%
# slice_sample(prop = 0.1) %>%
ggplot(aes(
x = cell_pair,
y = interacting_pair,
color = `-log10(pvalue)`
)) +
geom_point(shape = 22, aes(fill = mean, size = score), stroke = 2) +
ggpubr::theme_pubclean() +
labs(
x = NULL,
y = NULL,
title = str_glue("{name} - ligand"),
fill = "Log2 mean(Molecule 1, Molecule 2)",
# size = "-log10(p value)"
size = "Interaction scores"
) +
facet_wrap(~ group, ncol = 3) +
scale_fill_gradientn(colors = c("#fce5d7", "#f5ae8b",  "#d00f19" ,  "#6e020f")) +
scale_size(range = c(0.1, 5)) +
scale_color_gradient(low = "#386cb044", high = "#386cb0") +
theme(
axis.text.x = element_text(
angle = 90,
hjust = 1,
vjust = 0.5
),
legend.position = "right"
)
})
# receptor
p_dot_ls3 <-
dot_dat_ls %>%
enframe(name = "group") %>%
unnest(value) %>%
mutate(group_receptor = str_remove(cell_pair, "^.*\\|")) %>%
group_by(group_receptor) %>%
nest() %>%
deframe() %>%
imap(function(dat, name){ # browser()
dat %>%
as.data.frame() %>%
# mutate(
#   cell_pair2 =
#     cell_pair %>%
#     str_split("\\|") %>%
#     map(~ str_remove_all(.x, "\\s.*$")) %>%
#     map_chr(str_c, collapse = "|")
# ) %>%
# slice_sample(prop = 0.1) %>%
ggplot(aes(
x = cell_pair,
y = interacting_pair,
color = `-log10(pvalue)`
)) +
geom_point(shape = 22, aes(fill = mean, size = score), stroke = 2) +
ggpubr::theme_pubclean() +
labs(
x = NULL,
y = NULL,
title = str_glue("{name} - receptor"),
fill = "Log2 mean(Molecule 1, Molecule 2)",
# size = "-log10(p value)"
size = "Interaction scores"
) +
facet_wrap(~ group, ncol = 3) +
scale_fill_gradientn(colors = c("#fce5d7", "#f5ae8b",  "#d00f19" ,  "#6e020f")) +
scale_size(range = c(0.1, 5)) +
scale_color_gradient(low = "#386cb044", high = "#386cb0") +
theme(
axis.text.x = element_text(
angle = -90,
hjust = 0,
vjust = 0.5
),
legend.position = "right"
)
})
pdf(
file.path(cpdb_plot_out_dir, str_glue("dotplot2.pdf")),
height = 20,
width  = 15
)
p_dot_ls2 %>% walk(print)
dev.off()
##
pdf(
file.path(cpdb_plot_out_dir, str_glue("dotplot3.pdf")),
height = 20,
width  = 15
)
p_dot_ls3 %>% walk(print)
dev.off()
colnames(All_T@meta.data)
table(All_T$group)
All_T_meta <-
All_T[[]] %>%
rownames_to_column("cell_id")
All_T_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>%
ggplot(aes(x = group, y = freq, fill = subtype)) +
geom_col() +
ggsci::scale_fill_igv() +
theme(
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
p_col_proj <-
All_T_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>%
ggplot(aes(x = group, y = freq, fill = subtype)) +
geom_col() +
ggsci::scale_fill_igv() +
theme(
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
dr_dir
dr_dir <- 'proportion'
ggsave(
p_col_proj,
filename = file.path("T_subtype", "cell_freq_group.pdf"),
height = 6,
width  = 8
)
ggsave(
p_col_proj,
filename = file.path("T_subtype", "cell_freq_group.pdf"),
height = 6,
width  = 4
)
ggsave(
p_col_proj,
filename = file.path("T_subtype", "cell_freq_group.pdf"),
height = 6,
width  = 6
)
# cell props barplot  --------
Fibroblast_meta <-
Fibroblast[[]] %>%
rownames_to_column("cell_id")
Fibroblast_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>%
ggplot(aes(x = group, y = freq, fill = subtype)) +
geom_col() +
ggsci::scale_fill_igv() +
theme(
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
p_col_proj <-
Fibroblast_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>%
ggplot(aes(x = group, y = freq, fill = subtype)) +
geom_col() +
ggsci::scale_fill_igv() +
theme(
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
ggsave(
p_col_proj,
filename = file.path("Fibroblast_subtype", "cell_freq_group.pdf"),
height = 6,
width  = 6
)
Fibroblast_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count))
Fibroblast_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>% write.csv('Fibroblast_subtype/proportion.csv')
Fibroblast_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>% write.csv('Fibroblast_subtype/proportion.csv',row.names = FALSE)
All_T_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>% write.csv('T_subtype/proportion.csv',row.names = FALSE)
Macrophage_meta <-
Macrophage[[]] %>%
rownames_to_column("cell_id")
p_col_proj <-
Macrophage_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>%
ggplot(aes(x = group, y = freq, fill = subtype)) +
geom_col() +
ggsci::scale_fill_igv() +
theme(
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
# export
ggsave(
p_col_proj,
filename = file.path("Macrophage_subtype", "cell_freq_group.pdf"),
height = 6,
width  = 6
)
Macrophage_meta %>%
group_by(group, subtype) %>%
summarise(count = n()) %>%
mutate(freq = count/sum(count)) %>% write.csv('Macrophage_subtype/proportion.csv',row.names = FALSE)
table(all_data$annotation)
??Vlnplot
Duct_cells_1 <- subset(all_data,annotation == 'Duct_cells_1')
VlnPlot(Duct_cells_1,c('Csf1','Vegfa'),group.by = 'group')
Idents(all_data)<- 'annotation'
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1')
# normalize ------
all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1')
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE)
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE,y.max = 2)
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE)
p<- VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE)
plot_data <- VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE,return.data = TRUE)
expression_data <- all_data[which(meta_data$group == "Duct_cells_1"), c("Csf1", "Vegfa")]
expression_data <- all_data[which(all_data@meta.data$group == "Duct_cells_1"), c("Csf1", "Vegfa")]
expression_data <- all_data[which(all_data$annotation == "Duct_cells_1"), c("Csf1", "Vegfa")]
p
??FetchData
subset(all_data,annotation=='Duct_cells_1') %>% FetchData(vars = c("Csf1", "Vegfa",slot = 'data') -> temp
()
)))))\
subset(all_data,annotation=='Duct_cells_1') %>% FetchData(vars = c("Csf1", "Vegfa",slot = 'data')) -> temp
head(temp)
Duct_cells_1 @meta.data %>% colnames()
head(Duct_cells_1)
Duct_cells_1
temp %>%
rownames_to_column('cell_id') %>%
left_join(Duct_cells_1[[]] %>%
select(group) %>%
column_to_rownames('cell_id'),by = 'cell_id')
temp %>%
rownames_to_column('cell_id')
Duct_cells_1[[]] %>%
select(group)
Duct_cells_1[[]] %>%
select(group) %>%
column_to_rownames('cell_id')
Duct_cells_1[[]] %>%
select(group) %>% head()
temp %>%
rownames_to_column('cell_id') %>%
left_join(Duct_cells_1[[]] %>%
select(group) %>%
rownames_to_column('cell_id'),by = 'cell_id')
getwd()
temp %>%
rownames_to_column('cell_id') %>%
left_join(Duct_cells_1[[]] %>%
select(group) %>%
rownames_to_column('cell_id'),by = 'cell_id') %>%
column_to_rownames('cell_id') %>%
write.csv('vlnplot_data.csv')
savehistory("~/history.R")
