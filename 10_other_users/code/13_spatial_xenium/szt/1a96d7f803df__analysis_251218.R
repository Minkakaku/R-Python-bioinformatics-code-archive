
# 1. 细胞比例图 ------------
library(Seurat)
library(tidyverse)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

dir.create("results_20251218")

Fig1 <- file.path("results_20251218", "Fig1")
dir.create(Fig1)

all_data <- readRDS("all_data_DC_annotation5_1208.rds")
Idents(all_data) <- 'annotation5'

non_na_cells <- which(!is.na(all_data@meta.data$annotation5))
all_data_filt <- all_data[, non_na_cells]

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


data <- all_data_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 44)]
colnames(data) <- c("Cells","Sample","Group","Annotation")

data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(Fig1, "Cell_type_stat.txt"),quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix) <- mymtrix[,"Sample"]
mymtrix <- mymtrix[,-1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file = file.path(Fig1, "Cell_type_stat_percent.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file=file.path(Fig1, "Cell_type_stat_group.txt"),quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file=file.path(Fig1, "Cell_type_stat_percent_group.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt",header=T,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))

# orders <- read.table("order.txt",header=F,sep="\t",quote="",check.names=F)
# data4plot1$Celltype <- factor(data4plot1$Celltype,levels=orders[,1])
# data4plot2$Celltype <- factor(data4plot2$Celltype,levels=orders[,1])


plot_list <- c()
plot_list[[1]] <- ggplot(data4plot1, aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) # 两张图各自的宽度
ggsave(p,file = file.path(Fig1, "ratio_barplot.pdf"), width=12, height=6)
ggsave(p,file = file.path(Fig1, "ratio_barplot.png"), width=12, height=6, dpi=300)









# 2. Proximity analysis -----------
source("spatial analysis.R")
### 1. input data ------------
cells_position_ls <- readRDS(file.path("R_OUTPUT", "data_for_analysis", "cells_position_ls.rds"))

All_type_Proximity_out_dir <- "results_20251218/Fig2 Proximity analysis"
dir.create(All_type_Proximity_out_dir)

# all_data_cn_meta
All_type_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    proj_meta_sub <-
      all_data_filt[[]] %>%  
      dplyr::select(annotation5, group, sample) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = annotation5,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

All_type_meta_pos_subset <- All_type_meta_pos %>% rownames_to_column("cell_id")

saveRDS(All_type_meta_pos_subset, file = file.path(All_type_Proximity_out_dir, "Subtype_meta_pos_subset.rds"))

All_type_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))


# summarized by group and celltype
pa_groups <-
  with(
    All_type_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

# pa_groups2 <-
#   All_type_meta_pos_subset %>%
#   left_join(All_type_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
#   with(
#     str_c(Classification1, group, Slide_ID, CN, sep = '..')
#   )
# 
# pa_groups3 <-
#   All_type_meta_pos_subset %>%
#   left_join(All_type_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
#   with(
#     str_c(Classification1, group, Slide_ID, CN, sep = '..')
#   )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
# res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
# res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 


# res_pa_df2  <-
#   res_pa_df2 %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
#   separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
# res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))
# 
# 
# res_pa_df3  <-
#   res_pa_df3 %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
#   separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
# res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

# p_res_pa_ls_cn10 <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN)) %>%
#   imap(function(dat, name){
#     dat %>%
#       # dplyr::filter(Target == "Epithelial cells") %>%
#       ggplot(aes(x = Others, y = cell_props, fill = group)) +
#       geom_boxplot(width = 0.45, outlier.colour = NA) +
#       # geom_jitter(width = 0.25) +
#       facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#       ggpubr::theme_pubclean() +
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#       ggtitle(name)
#   })
# 
# p_res_pa_ls_cn15 <-
#   res_pa_df3 %>%
#   with(split(res_pa_df3, CN)) %>%
#   imap(function(dat, name){
#     dat %>%
#       # dplyr::filter(Target == "Epithelial cells") %>%
#       ggplot(aes(x = Others, y = cell_props, fill = group)) +
#       geom_boxplot(width = 0.45, outlier.colour = NA) +
#       # geom_jitter(width = 0.25) +
#       facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#       ggpubr::theme_pubclean() +
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#       ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis.pdf"),
  height = 25,
  width  = 12
)
p_res_pa%>% print()
dev.off()


# pdf(
#   file.path(All_type_Proximity_out_dir, "Proximity analysis - CN10.pdf"),
#   height = 25,
#   width  = 12
# )
# p_res_pa_ls_cn10 %>% print()
# dev.off()
# 
# 
# pdf(
#   file.path(All_type_Proximity_out_dir, "Proximity analysis - CN15.pdf"),
#   height = 25,
#   width  = 12
# )
# p_res_pa_ls_cn15 %>% print()
# dev.off()




# 10 50 100 150 200
res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, #### <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, #### <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, #### <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, #### <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, #### <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))




proximity_cells <- c(
  "Duct_cells_1", "CD8_naive", "CD8_Tex", "DC"
  # "GME"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){ # browser()
    cell_idx <- All_type_meta_pos_subset$Classification1 == cell
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- All_type_meta_pos_subset[cell_idx, ]
    # 2 res_pa
    res_pa_sub_ls <- res_pa_ls %>% map(~ .x[, cell_idx])
    # distance diff
    diff_idx_ls <- seq_along(res_pa_sub_ls) %>% rev() %>% { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>% apply(1, c, simplify = FALSE)
    # diff_idx_ls
    res_pa_sub_diff_ls <- vector(length = 5, mode = "list")
    for(idx in diff_idx_ls){
      res_pa_sub_diff_ls[[ idx[1] ]] <- res_pa_sub_ls[[ idx[1] ]] - res_pa_sub_ls[[ idx[2] ]]
    }
    res_pa_sub_diff_ls[[1]] <- res_pa_sub_ls[[1]]
    names(res_pa_sub_diff_ls) <- names(res_pa_sub_ls)
    
    res_pa_sub_count <-
      res_pa_sub_diff_ls %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_count <-
      res_pa_sub_count %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count
    )
  })


p_proximity_cell_count_distribution_ls <- 
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){ # browser()
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })



res_pa_count_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    res_pa_count_tidy <-
      1:nrow(res_pa_count_raw) %>%  
      map(function(x){
        count <- res_pa_count_raw$count[x] 
        if(count < 1) return(NULL) # <- count是median count，有可能是0.5
        xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
        xx
      }) %>%
      purrr::discard(is.null) %>%
      bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    
    res_pa_count_tidy
  })


# only selecet two macrohpages
p_proximity_cell_count_distribution_ls2 <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      # mutate(distance = as.numeric(distance)) %>%
      mutate(distance = factor(distance, levels = c("10", "50", "100", "150", "200"))) %>%
      # dplyr::filter(Cell %in% c("C1QC Macrophage", "SELENOP Macrophage")) %>%
      ggplot(aes(x = Cell, y = count, fill = distance)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ Slide_ID, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })

# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
    #   ggplot(aes(x = distance, fill = Group)) + #, color = Group)) +
    #   geom_density(alpha = 0.6) +
    #   # geom_smooth(aes(y = count), method = "loess") +
    #   facet_wrap(~ Cell, scales = "free_y") +
    #   theme_bw() +
    #   scale_fill_manual(values = grp_cols_sp) +
    #   labs(
    #     x = 'Distance(um)',
    #     y = 'Median count'
    #   )
    
    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      scale_color_manual(values = plan) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })



# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -1.pdf"), width = 14, height = 24)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -2.pdf"), width = 14, height = 24)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - smooth_lines.pdf"), width = 18, height = 12)
p_proximity_ls %>% walk(print)
dev.off()


# radar plot 

p_proximity_radar_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%
      
      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = group)) +
      geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
      geom_hline(yintercept = c(10, 50, 100, 150, 200), color = "gray") +
      geom_point(position = position_jitterdodge(
        jitter.width = 0,
        jitter.height = 10,
        dodge.width = 0.75
      )) +
      geom_text(
        data = data.frame(x = 0, y = c(10, 50, 100, 150, 200)),
        aes(x = x, y = y, label = paste0(y, "um")),
        inherit.aes = FALSE,
        vjust = -1.5
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1,1,1,6, 'mm')
      ) +
      scale_size(range = c(3, 12)) +
      scale_y_continuous(limits = c(-100, 220)) +
      scale_color_manual(values = plan)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - radarplot.pdf"), width = 20, height = 15)
p_proximity_radar_ls %>% walk(print)
dev.off()












# 3. 所有亚群细胞比例图 ------------
library(Seurat)
library(tidyverse)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

# dir.create("results_20251218")

Fig3 <- file.path("results_20251218", "Fig3")
dir.create(Fig3)

# all_data <- readRDS("all_data_DC_annotation5_1208.rds")


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


data <- all_data_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 45)]
colnames(data) <- c("Cells","Sample","Group","Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(Fig3, "Cell_type_stat.txt"),quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix) <- mymtrix[,"Sample"]
mymtrix <- mymtrix[,-1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file = file.path(Fig3, "Cell_type_stat_percent.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file=file.path(Fig3, "Cell_type_stat_group.txt"),quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file=file.path(Fig3, "Cell_type_stat_percent_group.txt"),quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt",header=T,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))

# orders <- read.table("order.txt",header=F,sep="\t",quote="",check.names=F)
# data4plot1$Celltype <- factor(data4plot1$Celltype,levels=orders[,1])
# data4plot2$Celltype <- factor(data4plot2$Celltype,levels=orders[,1])


plot_list <- c()
plot_list[[1]] <- ggplot(data4plot1, aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) # 两张图各自的宽度
ggsave(p,file = file.path(Fig3, "ratio_barplot.pdf"), width=12, height=6)
ggsave(p,file = file.path(Fig3, "ratio_barplot.png"), width=12, height=6, dpi=300)







# 4. 巨噬亚群细胞比例图 ------------
library(Seurat)
library(tidyverse)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

# dir.create("results_20251218")

Fig4 <- file.path("results_20251218", "Fig4")
dir.create(Fig4)


Macrophage <- readRDS("Macrophage_subtype/Second/Macrophage_2.rds")

head(Macrophage)
unique(Macrophage$subtype)


data <- Macrophage[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 45)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(Fig4, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(Fig4, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(Fig4, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(Fig4, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

# orders <- read.table("order.txt", header = F, sep = "\t", quote = "", check.names = F)
# data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
# data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])


plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Sample", y = "Percent(%)", title = "") +
  scale_fill_manual(values = plan) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2)) + NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Group", y = "", title = "") +
  scale_fill_manual(values = plan) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(widths = c(3, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(Fig4, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(Fig4, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)












# 5. T cell 比例图 ------------
library(Seurat)
library(tidyverse)
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

# dir.create("results_20251218")

Fig5 <- file.path("results_20251218", "Fig5")
dir.create(Fig5)


sc_Tcell <- readRDS("T_subtype/third/T_subtype.rds")

head(sc_Tcell)
unique(sc_Tcell$subtype)


data <- sc_Tcell[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 46)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(Fig5, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(Fig5, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(Fig5, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(Fig5, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

# orders <- read.table("order.txt", header = F, sep = "\t", quote = "", check.names = F)
# data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
# data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])


plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Sample", y = "Percent(%)", title = "") +
  scale_fill_manual(values = plan) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2)) + NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Group", y = "", title = "") +
  scale_fill_manual(values = plan) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(widths = c(3, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(Fig5, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(Fig5, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)











# 6. Macrophage 拟时序 ------------
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)


monocle_outdir <- "results_20251218/Monocle_Macrophage"
dir.create(monocle_outdir)
# expr_dat_monocle <- as(LayerData(JoinLayers(sc_fib), "counts"), "sparseMatrix")
# Macrophagell<-subset(Macrophage,subtype == 'M2'|subtype == 'M1'|subtype == 'Cycling M2')

count <- GetAssayData(object = JoinLayers(Macrophage)[["RNA"]], layer = "counts")
fdata <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))
pdata <- Macrophage@meta.data
#构建monocle对象
cds <- newCellDataSet(count,
                      featureData = new("AnnotatedDataFrame", data =fdata),
                      phenoData = new("AnnotatedDataFrame", data = pdata),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#标准化sacle及筛选，即使原先的Seurat对象已经归一化过，官方推荐在monocle中重新进行归一化处理
#标准化细胞之间的mRNA的差异
cds <- estimateSizeFactors(cds)
#离散度值可以帮助我们进行后续的差异分析
cds <- estimateDispersions(cds)

#在fData(cds)中添加“num_cell_expressed”列，可用于过滤细胞
cds<- detectGenes(cds, min_expr = 0.1)

head(pData(cds))
head(fData(cds))

expressed_genes<-row.names(subset(fData(cds),num_cells_expressed >= 10))
length(expressed_genes)

#向Epidata(cds)中添加“UMI”列
pdata$UMI <- Matrix::colSums(exprs(cds))



### 筛选高变基因(离散度高的gene) ----------
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit) 
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-1.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

cds=orderCells(cds,root_state =3) 
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-1.pdf"))
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-1.pdf"), height = 10, width = 25)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-1.pdf"), height = 10, width = 10)
print(p)
dev.off()



### 高变基因 -----------
# group1_data<- group1_data[['RNA']] %>% as('Assay5')
# group1_data <- group1_data%>%NormalizeData(verbose=FALSE)%>%FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# group1_data <- FindVariableFeatures(group1_data, selection.method = "vst", nfeatures = 2000)
ordergene <- VariableFeatures(Macrophage)
cds <- setOrderingFilter(cds, ordergene)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-2.pdf"), height = 10, width = 25)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-2.pdf"), height = 10, width = 10)
print(p)
dev.off()



### 差异基因 ----------
deg.cluster <- FindAllMarkers(Macrophage)
express_genes <- subset(deg.cluster,p_val_adj <0.05)$gene
express_genes <- unique(express_genes)
cds <- setOrderingFilter(cds, express_genes)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_subtype-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-3.pdf"), height = 10, width = 25)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-3.pdf"), height = 10, width = 10)
print(p)
dev.off()


##寻找拟时相关基因
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res = BEAM(cds[ordergene,], branch_point = 1, progenitor_method = 'duplicate')
BEAM_res = BEAM_res[, c("gene_short_name","pval","qval")]
BEAM_res = BEAM_res[order(BEAM_res$qval),]
sum(BEAM_res$pval < 0.05)
tmp1 <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, pval < 0.05)),],
                                    branch_point = 2,
                                    num_clusters = 3, #这些基因被分成几个group
                                    cores = 1,
                                    #branch_labels = c("Cell fate 1", "Cell fate 2"),
                                    #hmcols = NULL, #默认值
                                    #hmcols = colorRampPalette(c("navy","white","firebrick3"))(50),
                                    # hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                    hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(60),
                                    branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                    use_gene_short_name = T,
                                    show_rownames =T,
                                    return_heatmap = T)
write.csv(BEAM_res,'./monocle/group2_branch2.csv')
pdf("./monocle/group2_branch2_heatmap.pdf",width =10,height = 10)
tmp1$ph_res
dev.off()

plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE) 
plot_cell_trajectory(cds,color_by="subtype", size=1, show_backbone=TRUE) + scale_color_manual(values = plan)
plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap("~subtype", nrow = 1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = plan)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype") + scale_color_manual(values = plan)

# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

pdf(file.path(monocle_outdir, "cell trajectory.pdf"), width = 8, height = 4)
plot_cell_trajectory(cds,color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds,color_by = "subtype", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - - split.pdf"), width = 15, height = 8)
plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap("~subtype", nrow = 1) + theme(legend.position = "right") + scale_color_manual(values = plan) 
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = plan)
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - annotation - split by group.pdf"), width = 16, height = 4)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("CON", "IRE_4", "IRE_7"))
  plot_cell_trajectory(cds, color_by = "subtype",size=1,show_backbone=TRUE) +
    scale_color_manual(values = plan) +
    facet_wrap(~ group) +
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(monocle_outdir, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype") +
  scale_color_manual(values = plan) +
  facet_wrap(~ subtype)
dev.off()

# cell barplot 
p_col_cell_state <-
  cds %>% 
  pData() %>%
  dplyr::select(subtype, State) %>%
  group_by(subtype, State) %>%
  mutate(subtype = as.character(subtype)) %>%
  mutate(subtype = factor(subtype, levels = c(
    'M0',
    'M1',
    'M2'
  ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = plan) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_outdir, "cell stat barplot.pdf"), height = 4, width = 8)

# de analysis
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 1, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(BEAM_res, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)  
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)


# BEAM analysis
BEAM_res <- BEAM(cds[ordergene,],branch_point =1,progenitor_method = 'duplicate')
# BEAM_res <- BEAM(cds[express_genes_monocle,], branch_point = 2, cores = 2) 

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res,file.path(monocle_out_dir, "BEAM_res1.csv"), row.names = F)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-20)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 3, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了
#选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  arrange(qval) %>%
  slice_head(n = 100) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 1, 
                                 num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave(file.path(monocle_out_dir,"BEAM_heatmap1.pdf"), p$ph_res, width = 10, height = 10)














# 7. T cell 拟时序 ------------
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)

AllT <- readRDS("T_subtype/third/T_subtype.rds")



monocle_outdir <- "results_20251218/Monocle_AllT"
dir.create(monocle_outdir)
# expr_dat_monocle <- as(LayerData(JoinLayers(sc_fib), "counts"), "sparseMatrix")
# Macrophagell<-subset(Macrophage,subtype == 'M2'|subtype == 'M1'|subtype == 'Cycling M2')

count <- GetAssayData(object = JoinLayers(AllT)[["RNA"]], layer = "counts")
fdata <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))
pdata <- AllT@meta.data
#构建monocle对象
cds <- newCellDataSet(count,
                      featureData = new("AnnotatedDataFrame", data =fdata),
                      phenoData = new("AnnotatedDataFrame", data = pdata),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#标准化sacle及筛选，即使原先的Seurat对象已经归一化过，官方推荐在monocle中重新进行归一化处理
#标准化细胞之间的mRNA的差异
cds <- estimateSizeFactors(cds)
#离散度值可以帮助我们进行后续的差异分析
cds <- estimateDispersions(cds)

#在fData(cds)中添加“num_cell_expressed”列，可用于过滤细胞
cds <- detectGenes(cds, min_expr = 0.1)

head(pData(cds))
head(fData(cds))

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)

#向Epidata(cds)中添加“UMI”列
pdata$UMI <- Matrix::colSums(exprs(cds))



### 筛选高变基因(离散度高的gene) ----------
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit) 
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype, ncol = 4)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-1.pdf"), height = 12, width = 25)
print(p)
dev.off()

# cds = orderCells(cds, root_state =3) 
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-1.pdf"))
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-1.pdf"), height = 10, width = 10)
print(p)
dev.off()



### 高变基因 -----------
# group1_data<- group1_data[['RNA']] %>% as('Assay5')
# group1_data <- group1_data%>%NormalizeData(verbose=FALSE)%>%FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# group1_data <- FindVariableFeatures(group1_data, selection.method = "vst", nfeatures = 2000)
ordergene <- VariableFeatures(AllT)
cds <- setOrderingFilter(cds, ordergene)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-2.pdf"), height = 10, width = 10)
print(p)
dev.off()

cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-2.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype, ncol = 4)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-2.pdf"), height = 12, width = 25)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-2.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-2.pdf"), height = 10, width = 10)
print(p)
dev.off()



### 差异基因 ----------
deg.cluster <- FindAllMarkers(AllT)
express_genes <- subset(deg.cluster,p_val_adj < 0.05)$gene
express_genes <- unique(express_genes)
cds <- setOrderingFilter(cds, express_genes)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_subtype-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype, ncol = 4)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-3.pdf"), height = 12, width = 25)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-3.pdf"), height = 10, width = 10)
print(p)
dev.off()


##寻找拟时相关基因
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res = BEAM(cds[ordergene,], branch_point = 1, progenitor_method = 'duplicate')
BEAM_res = BEAM_res[, c("gene_short_name","pval","qval")]
BEAM_res = BEAM_res[order(BEAM_res$qval),]
sum(BEAM_res$pval < 0.05)
tmp1 <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, pval < 0.05)),],
                                    branch_point = 2,
                                    num_clusters = 3, #这些基因被分成几个group
                                    cores = 1,
                                    #branch_labels = c("Cell fate 1", "Cell fate 2"),
                                    #hmcols = NULL, #默认值
                                    #hmcols = colorRampPalette(c("navy","white","firebrick3"))(50),
                                    # hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                    hmcols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(60),
                                    branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                    use_gene_short_name = T,
                                    show_rownames =T,
                                    return_heatmap = T)
write.csv(BEAM_res,'./monocle/group2_branch2.csv')
pdf("./monocle/group2_branch2_heatmap.pdf",width =10,height = 10)
tmp1$ph_res
dev.off()

plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE) 
plot_cell_trajectory(cds,color_by="subtype", size=1, show_backbone=TRUE) + scale_color_manual(values = plan)
plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap("~subtype", nrow = 1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = plan)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype") + scale_color_manual(values = plan)

# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

pdf(file.path(monocle_outdir, "cell trajectory.pdf"), width = 8, height = 4)
plot_cell_trajectory(cds,color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds,color_by = "subtype", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - - split.pdf"), width = 15, height = 8)
plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap("~subtype", nrow = 1) + theme(legend.position = "right") + scale_color_manual(values = plan) 
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = plan)
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - annotation - split by group.pdf"), width = 16, height = 4)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("CON", "IRE_4", "IRE_7"))
  plot_cell_trajectory(cds, color_by = "subtype",size=1,show_backbone=TRUE) +
    scale_color_manual(values = plan) +
    facet_wrap(~ group) +
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(monocle_outdir, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype") +
  scale_color_manual(values = plan) +
  facet_wrap(~ subtype)
dev.off()

# cell barplot 
p_col_cell_state <-
  cds %>% 
  pData() %>%
  dplyr::select(subtype, State) %>%
  group_by(subtype, State) %>%
  mutate(subtype = as.character(subtype)) %>%
  mutate(subtype = factor(subtype, levels = c(
    'M0',
    'M1',
    'M2'
  ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = plan) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_outdir, "cell stat barplot.pdf"), height = 4, width = 8)

# de analysis
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 1, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(BEAM_res, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)  
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)


# BEAM analysis
BEAM_res <- BEAM(cds[ordergene,],branch_point =1,progenitor_method = 'duplicate')
# BEAM_res <- BEAM(cds[express_genes_monocle,], branch_point = 2, cores = 2) 

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res,file.path(monocle_out_dir, "BEAM_res1.csv"), row.names = F)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-20)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 3, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了
#选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  arrange(qval) %>%
  slice_head(n = 100) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 1, 
                                 num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave(file.path(monocle_out_dir,"BEAM_heatmap1.pdf"), p$ph_res, width = 10, height = 10)





































