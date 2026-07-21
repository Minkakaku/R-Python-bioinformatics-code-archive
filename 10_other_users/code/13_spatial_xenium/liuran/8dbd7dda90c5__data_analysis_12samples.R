
library(tidyverse)
library(Seurat)
library(patchwork)


dir.create("12samples_analysis")
setwd("12samples_analysis")



# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)


list.files("/data/xenium/FWSP20253423")

# library --------------
# copy jar file.
bf_cache_dir <- tools::R_user_dir("RBioFormats", which = "cache")
dir.create(bf_cache_dir)
bf_jar_path <- list.files(system.file("java", package = "RBioFormats"), pattern = "jar", full.names = TRUE)
bf_ver <- read.dcf(file.path(system.file(package = "RBioFormats"), "DESCRIPTION"), "BioFormats")
bf_cache_name <- sprintf("bioformats_package_%s.jar", bf_ver)
file.copy(bf_jar_path, file.path(bf_cache_dir, bf_cache_name), overwrite = TRUE)

options(timeout = 10 ^ 4)
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
library(RBioFormats)
source("xenium.R")
options(timeout = 60)


# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()

# options(timeout = 10000)

# setting variable   --------------------
data_dir <- "/data/xenium/FWSP20253423"
anno_file <- "~/data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "SAM000626 xenium-human-gene-expression-panel-gene-list.xlsx"


# reading data -------------
# ## data path and project name  -----
# parse_data_path <-
#   dir(data_dir, full.names = TRUE) %>% 
#   parseProjName_xenium(pattern = "output[-]?") %>%
#   mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
#   identity()
# 
# ## add id for parse_data_path --------
# parse_data_path <- 
#   parse_data_path %>%
#   mutate(.id = seq_along(.proj))



## image w/h ratio ------
# used to fix width and heigth ratio in In Situ Spatial plot.
# 这段代码的主要目的是计算每个图像的长宽比 (width-to-height ratio)，并将其存储在列表 image_wh_r_ls 中
# 现在可以用 coord_equal 调整宽高比，这段代码可以不跑
# image_wh_r_ls <-
#   with(parse_data_path, split(.path, .proj)) %>%
#   map(function(x){
#     imag_wh <-
#       RBioFormats::read.metadata(
#         file.path(x, "morphology.ome.tif")
#       ) %>%
#       .[[1]] %>%
#       .$coreMetadata %>%
#       .[c("sizeX", "sizeY")]
# 
#     image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
#     image_wh_r
#   })


## read into annotation -------
# if(file.exists(anno_file)){
#   anno_data <- readr::read_csv(anno_file)
#   User_columns <- colnames(anno_data)
# }else{
#   anno_data <- NULL
#   User_columns <- "project"
# }
# 
# 
# ## seurat object ------
# all_data_ls <-
#   parse_data_path %>% 
#   dplyr::select(.counts_path, .proj, .id) %>%
#   pmap(function(.counts_path, .proj, .id){
#     # merge anno_data
#     if(!is.null(anno_data)){
#       row_idx <- which(anno_data$project == .proj)
#       col_idx <- anno_data %>% colnames() 
#       
#       # 没有找到当前样本的注释数据
#       #<<<<<< 所以Annotation.csv也有筛选文件的作用 >>>>>#
#       if(length(row_idx) < 1){
#         message(str_glue("Not found sample {.proj} in column: project of Annotation table.") )
#         return(NULL)
#       }
#       # 没有任何数据，报错
#       if(length(col_idx) < 1) stop("Not found any data in Annotation table.")
#       
#       meta_dat <- anno_data %>% .[row_idx, col_idx] %>% as.list()
#     }else{
#       meta_dat <- list(project = .proj)
#     }
#     
#     # 构建Seurat
#     sc <-
#       Read10X(.counts_path) %>% 
#       .[["Gene Expression"]] %>%
#       CreateSeuratObject(
#         project = .proj,
#         min.cells    = 5,
#         min.features = 5
#       ) %>%
#       AddMetaData(meta_dat) %>%
#       RenameCells(add.cell.id = .id)
#     sc
#   }) 
# 
# all_data_ls <- all_data_ls %>% purrr::discard(is.null)
# names(all_data_ls) <- map_chr(all_data_ls, "project.name")
# 
# 
# # parse ROI----------------
# # sample_df <-  
# #   data.frame(
# #     sample = c("kxz1","kxz6","kxz7","kcz2","kcz3","kcz8","gxz10","gxz11","gxz12","gcz13","gcz14","gcz15"),
# #     group  = c("AirMale","AirMale","AirMale","AirFemale","AirFemale","AirFemale","HyperoxiaMale","HyperoxiaMale","HyperoxiaMale","HyperoxiaFemale","HyperoxiaFemale","HyperoxiaFemale")
# #   )
# 
# ### sample_df中的内容：sample为样本名，group为对应FOV名
# sample_df <- 
#   data.frame(
#     sample = c("B1-1","B1-2","B1-3","B2-1","B2-2","B2-3","B3-1","B3-2","B3-3","BD1-1","BD1-2","BD1-3","BD2-1","BD2-2","BD2-3","BD3-1","N6","BD3-2","C1","C2"),
#     group  = c("B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B1-1_B1-2_B1-3_B2-1_B2-2_B2-3_B3-1_B3-2",
#                "B3-3_BD1-1_BD1-2_BD1-3",
#                "B3-3_BD1-1_BD1-2_BD1-3",
#                "B3-3_BD1-1_BD1-2_BD1-3",
#                "B3-3_BD1-1_BD1-2_BD1-3",
#                "BD2-1",
#                "BD2-2_BD2-3_BD3-1_N6",
#                "BD2-2_BD2-3_BD3-1_N6",
#                "BD2-2_BD2-3_BD3-1_N6",
#                "BD2-2_BD2-3_BD3-1_N6",
#                "BD3-2_C1_C2",
#                "BD3-2_C1_C2",
#                "BD3-2_C1_C2")
#   )
# 
# 
# roi_ls <- list.files(
#   'ROI',
#   full.names = TRUE
# ) 
# 
# names(roi_ls) <- 
#   roi_ls %>%
#   basename() %>%
#   str_remove("\\..*")
# 
# 
# roi_coords_ls <-
#   roi_ls %>%
#   map(function(path){
#     # 读取文件内容，跳过前两行
#     dat <- read.csv(path,skip = 2)
#     # 将列名转换为小写
#     colnames(dat) <- dat %>% colnames() %>% str_to_lower()
#     # 添加一个新的列 'cell'，其值与 'cell.id' 列相同
#     dat<-dat %>% mutate(cell=cell.id)%>%
#       # 将 'cell.id' 列转换为行名
#       column_to_rownames('cell.id')
#     dat
#   }) 
# 
# stopifnot(all(names(roi_coords_ls) %in% sample_df$sample))
# roi_coords_ls <- roi_coords_ls[sample_df$sample]
# 
# 
# ### 根据每个样本所在芯片的顺序编号，来对应all_data_ls中的cell id
# .id<-c("1_","1_","1_","1_","1_","1_","1_","1_","2_","2_","2_","2_","3_","4_","4_","4_","4_","5_","5_","5_")
# # .id<- str_c(names(roi_coords_ls), "_")
# # .id
# names(.id)<-names(roi_coords_ls)
# roi_coords_ls <- roi_coords_ls %>%
#   map2(.id, ~ {
#     .x$cell <- str_c(.y, .x$cell)
#     return(.x)
#   })
# 
# 
# xenium_roi_ls <-
#   roi_coords_ls %>%
#   imap(function(roi, name){ # browser()
#     # get Visium Object
#     xenium_name <- sample_df$group[sample_df$sample == name]
#     xenium <- all_data_ls[[xenium_name]]
#     
#     # subset visium
#     cells <- roi$cell
#     xenium_roi <- subset(xenium, cells = cells)
#     xenium_roi
#   })
# 
# 
# # 计算 roi_coords_ls 中每个样本的细胞数
# cell_counts <- sapply(roi_coords_ls, function(x) length(x$cell))
# 
# # 输出每个样本的细胞数
# cell_counts
# 
# 
# # xenium_roi_ls$`F21-6989A`<-all_data_ls$`F21-6989A`
# # xenium_roi_ls$`F22-2468A`<-all_data_ls$`F22-2468A`
# 
# # add sample_id or sample_name -------------------
# xenium_roi_ls <-
#   xenium_roi_ls %>%
#   imap(function(xenium, name){ #browser()
#     xenium$project <- name
#     Idents(xenium) <- 'project'
#     # RenameCells(visium, add.cell.id = name)
#     xenium
#   })

# visium_roi override visium
# all_data_ls <- xenium_roi_ls

all_data_ls_fix <- readRDS("~/data_for_analysis/all_data_ls_fix.rds")
all_data_ls <- all_data_ls_fix

selected_samples <- c("B1-1","B1-2","B1-3","B2-1","B2-2","B2-3","BD1-1","BD1-2","BD1-3","BD2-1","BD2-2","BD2-3")
all_data_ls_filt <- all_data_ls[selected_samples]

all_data_ls <- all_data_ls_filt


cells_position_ls <- readRDS("~/data_for_analysis/cells_position_ls.rds")
cells_boundary_ls <- readRDS("~/data_for_analysis/cells_boundary_ls.rds")


# sample_df_fix <-
#   data.frame(
#     sample = c("B3-1","B3-2","B3-3","BD3-1","BD3-2"),
#     sample_new  = c("IPF3_1","IPF3_2","IPF3_3","Dex3_1","Dex3_2")
#   )
# 
# 
# all_data_ls_fix <- all_data_ls %>%
#   imap(function(data, name) {
#     # 获取正确的组名
#     correct_group <- sample_df_fix$sample_new[sample_df_fix$sample == name]
#     correct_project <- sample_df_fix$sample[sample_df_fix$sample == name]
#     # 检查是否找到正确的组名
#     if (length(correct_group) == 0) {
#       warning(paste("No group found for sample:", name))
#       return(data)
#     }
#     
#     # 添加或更新组名到元数据
#     data@meta.data$sample <- correct_group
#     data@meta.data$project <- correct_project
#     # 返回更新后的数据
#     data
#   })
# 
# 
# sample_df_group_fix <-
#   data.frame(
#     sample = c("B3-1","B3-2","B3-3","BD3-1","BD3-2"),
#     group  = c("G3","G3","G3","Dex_G3","Dex_G3")
#   )
# 
# all_data_ls_fix <- all_data_ls_fix %>%
#   imap(function(data, name) {
#     # 获取正确的组名
#     correct_group <- sample_df_group_fix$group[sample_df_group_fix$sample == name]
#     # correct_project <- sample_df_group_fix$sample[sample_df_group_fix$sample == name]
#     # 检查是否找到正确的组名
#     if (length(correct_group) == 0) {
#       warning(paste("No group found for sample:", name))
#       return(data)
#     }
#     
#     # 添加或更新组名到元数据
#     data@meta.data$group <- correct_group
#     # data@meta.data$project <- correct_project
#     # 返回更新后的数据
#     data
#   })
# 
# 
# 
# 
# # 验证组名是否更新正确
# print(unique(all_data_ls_fix$`B3-1`@meta.data$group))
# print(unique(all_data_ls_fix$`B3-1`@meta.data$project))
# print(unique(all_data_ls_fix$`B3-1`@meta.data$sample))
# 
# all_data_ls <- all_data_ls_fix

# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name  <- names(all_data_ls)




if(proj_count > 1){
  #merge
  all_data_ls[["Merge"]] <- mergeSeurat(all_data_ls)
}


# ## read coords -------------------
# # cell coordinate and cell bundary
# cells_position_ls <- 
#   # with(parse_data_path, split(.path, .proj)) %>%
#   parse_data_path %>%
#   group_by(.proj) %>%
#   nest %>%
#   deframe() %>%
#   map(function(x){
#     xx <-
#       readr::read_csv(
#         file.path(x$.path, "cells.csv.gz"),
#         col_types = list(cell_id = "c")
#       ) %>%
#       mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
#     xx
#   })
# 
# roi_position_ls <-
#   roi_coords_ls %>%
#   imap(function(roi, name){  #browser()
#     # get Visium Object
#     xenium_name <- sample_df$group[sample_df$sample == name]
#     xenium <- cells_position_ls[[xenium_name]]
#     
#     # subset visium
#     cells <- roi$cell
#     position_roi <-  xenium[which(xenium$cell_id %in% cells),]
#     position_roi
#   })
# 
# cells_position_ls <- cells_position_ls[proj_name]
# cells_position_ls <- cells_position_ls %>% purrr::discard(is.null)
# cells_position_ls <- c(roi_position_ls,cells_position_ls)
# 
# 
# cells_boundary_ls <- 
#   # with(parse_data_path, split(.path, .proj)) %>%
#   parse_data_path %>%
#   group_by(.proj) %>%
#   nest %>%
#   deframe() %>%
#   map(function(x){
#     xx <-
#       readr::read_csv(
#         file.path(x$.path, "cell_boundaries.csv.gz"),
#         col_types = list(cell_id = "c")
#       ) %>%
#       mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
#     xx
#   })
# 
# roi_boundary_ls <-
#   roi_coords_ls %>%
#   imap(function(roi, name){  #browser()
#     # get Visium Object
#     xenium_name <- sample_df$group[sample_df$sample == name]
#     xenium <- cells_boundary_ls[[xenium_name]]
#     
#     # subset visium
#     cells <- roi$cell
#     boundary_roi <- xenium[which(xenium$cell_id %in% cells),]
#     boundary_roi
#   })
# 
# cells_boundary_ls <- cells_boundary_ls[proj_name]
# cells_boundary_ls <- cells_boundary_ls %>% purrr::discard(is.null)
# cells_boundary_ls <- c(roi_boundary_ls,cells_boundary_ls)


selected_samples <- c("B1-1","B1-2","B1-3","B2-1","B2-2","B2-3","BD1-1","BD1-2","BD1-3","BD2-1","BD2-2","BD2-3")
cells_position_ls <- cells_position_ls[selected_samples]
cells_boundary_ls <- cells_boundary_ls[selected_samples]


# export data
dir.create("data_for_analysis") 
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))
saveRDS(all_data_ls, file.path("data_for_analysis", "all_data_ls.rds"))
# saveRDS(image_wh_r_ls, file.path("data_for_analysis", "image_wh_r_ls.rds"))
saveRDS(all_data_ls_fix, file.path("data_for_analysis", "all_data_ls_fix.rds"))




# QC ------------------
## QC_feature ----------
qc_metric <- c(
  "nFeature_RNA",
  "nCount_RNA", 
  "log10GenesPerUMI"
)


## add qc metric ----------
all_data_ls <- all_data_ls %>%
  map(function(sc){
    sc[["log10GenesPerUMI"]] <- log10(sc$nFeature_RNA) / log10(sc$nCount_RNA)
    sc
  })

## nRNA nFeature percent.mt ----------
mito_qc_ls <-
  all_data_ls %>% 
  map(function(sc){
    
    # use to determine the filtering criteria
    p1 <- VlnPlot(object = sc,
                  features = qc_metric,
                  group.by  =  "project",
                  ncol = 3,
                  # log = T,
                  pt.size = 0)
    # geom_hline(yintercept = 0.985)
    p1 <-
      p1 & scale_y_continuous(
        breaks = function(x){
          xx <- round(x)
          
          # 5代表分为5个breaks
          # 10代表以10的倍数进行划分breaks
          steps <- if(xx[2] < 100){
            scaled_times <- ceiling(100/xx[2])
            
            floor((scaled_times * xx[2] / 5 ) / 10 ) * 10 / scaled_times
          }else{
            floor((xx[2] / 5 ) / 10 ) * 10
          }
          
          # 最小值都是0
          bk <- seq(0, xx[2], by = steps)
          
          # 第一间隔细分
          bk_minor <- seq(bk[1], bk[2], length.out = 5)
          bk <- c(bk_minor, bk[-(1:2)])
          
          # 最后间隔细分
          bk <- rev(bk)
          bk_minor <- seq(bk[1], bk[2], length.out = 5)
          bk <- c(bk_minor, bk[-(1:2)])
          bk <- rev(bk)
          
          bk
        }
      )
    p1
  })

# if(proj_count > 1){
#   mito_qc_ls[["Merge_vlnplot"]] <- mito_qc_ls[["Merge"]]
#   mito_qc_ls[["Merge"]] <- NULL
#   mito_qc_ls[["Merge_density"]] <-
#     all_data_ls$Merge[[]] %>%
#     dplyr::select(project, all_of(qc_metric)) %>%
#     pivot_longer(
#       -project,
#       names_to = "QC_metric",
#       values_to = "Value"
#     ) %>%
#     mutate(QC_metric = factor(QC_metric, levels = qc_metric)) %>%
#     ggplot(aes(x = Value, fill = project)) +
#     geom_density(alpha = 0.4)+
#     facet_wrap(~QC_metric, scales = "free", ncol = 3)+
#     ggpubr::theme_pubclean() +
#     labs(x = NULL,
#          y = NULL)
# }


pdf(file.path(qc_outdir, "Vlnplot.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp))
    # cat("  \n\n")
  })
dev.off()

pdf(file.path(qc_outdir, "Vlnplot-log.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp) & scale_y_log10())
    # cat("  \n\n")
  })
dev.off()


# export mito_qc.csv
mito_qc_path <- file.path('data', "mito_qc.csv")
if(!file.exists(mito_qc_path)){
  # 上限
  # nFeature_cutoff <- all_data_ls[["Merge"]][["nFeature_RNA"]] %>% unlist %>% {mean(.)+3*sd(.)} %>% round()
  # nCount_cutoff <- all_data_ls[["Merge"]][["nCount_RNA"]] %>% unlist %>% {mean(.)+3*sd(.)} %>% round()
  
  if(proj_count < 2){
    nFeature_cutoff <- all_data_ls[[1]][["nFeature_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
    nCount_cutoff <- all_data_ls[[1]][["nCount_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
  }else{
    nFeature_cutoff <- all_data_ls[["Merge"]][["nFeature_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
    nCount_cutoff <- all_data_ls[["Merge"]][["nCount_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
  }
  
  
  origin_mito_qc <-
    data.frame(
      proj               = proj_name,
      nFeature_RNA       = paste0(">10"), # paste0(">10; <", nFeature_cutoff),
      # nCount_RNA         = paste0(">20"), # paste0(">20; <", nCount_cutoff),
      log10GenesPerUMI   = ">0.6"
    )
  
  readr::write_csv(origin_mito_qc, mito_qc_path)
}

## cell filter -----------

# read into qc criteria
origin_mito_qc <- 
  # readxl::read_excel(mito_qc_path) %>%
  readr::read_csv(mito_qc_path) %>%
  group_by(proj) %>%
  nest %>%
  deframe()

parsed_mito_qc <-
  origin_mito_qc%>%
  map(function(x){
    # browser()
    x %>% map(
      function(col){
        col %>% 
          str_split(pattern = "[ ;]+") %>%
          .[[1]] %>%
          map(parse_qc_criteria)
      }
    )
  })


if(proj_count > 1){
  all_data_ls$Merge <- NULL
}

# befor mito filter
before_filter_qc <- all_data_ls %>%
  map2(names(.), function(sc, proj){
    # befor mito filter
    p1 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "nFeature_RNA",
                         group.by = "project")
    
    
    p2 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "log10GenesPerUMI",
                         group.by = "project")
    
    p <- p1 + p2 +
      plot_annotation(title = str_glue("{proj}: Before.qc.filter")) + 
      plot_layout(nrow=1, guides = "collect") &
      theme(legend.position = "top")
    p
  })

# 去除基因总量异常的cell
all_data_ls_pass <- 
  all_data_ls %>%
  imap(qc_ls = parsed_mito_qc, function(sc, name, qc_ls){
    
    # qc参数
    if(!name %in% names(qc_ls)){
      stop(str_glue("Please check QC table, {name} not found qc criteria..."))
    }
    
    qc <- qc_ls[[name]]
    # browser()
    sc_have_name <- is.null(names(sc))
    qc_have_name <- is.null(names(qc))
    
    if(all(sc_have_name, qc_have_name)){
      stopifnot(all(names(sc) == names(qc)))
    }
    
    combined_idx <-
      qc %>%
      map2(names(.), function(conds, var){#  browser()
        # var: nFeture_RNA percent.mt
        combined_idx_per_var <-
          conds %>% 
          map(~eval_parsed_criteria(.x, others = sc[[var]], verbose = FALSE)) %>%
          purrr::reduce(`&`)
        combined_idx_per_var
      }) %>%
      purrr::reduce(`&`)
    
    if(sum(combined_idx) < 1){
      stop(str_glue("No cell pass the qc criteria, is qc criteria of {name} too much strict?"))
    }
    # return
    sc[, which(combined_idx)]
  })

# after mito filter
after_filter_qc <- all_data_ls_pass %>%
  map2(names(.), function(sc, proj){
    # befor mito filter
    p1 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "nFeature_RNA",
                         group.by = "project") 
    
    p2 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "log10GenesPerUMI",
                         group.by = "project")
    
    p <- p1 + p2 +
      plot_annotation(title = str_glue("{proj}: After.qc.filter")) + 
      plot_layout(nrow = 1, guides = "collect") &
      theme(legend.position = "top")
    p
  })


# print the image
qc_filter <- 
  map2(before_filter_qc, 
       after_filter_qc, 
       function(before, after){
         cowplot::plot_grid(before, after, ncol = 1)
       })





pdf(file.path(qc_outdir, "qc plot - pass.pdf"), height = 12, width = 18)
qc_filter %>%
  walk2(names(.), 
        function(p, name){
          # if(proj_count > 1){
          #   cat(str_glue("####  {name}\n\n"))
          # }
          print(p)
          # cat("  \n\n")
        })
dev.off()





### qc pass ---------
mito_qc_ls_pass <-
  all_data_ls_pass %>% 
  map(function(sc){
    
    # use to determine the filtering criteria
    p1 <- VlnPlot(object = sc,
                  features = qc_metric,
                  group.by  =  "project",
                  ncol = 3,
                  # log = T,
                  pt.size = 0)
    # geom_hline(yintercept = 0.985)
    p1 <-
      p1 & scale_y_continuous(
        breaks = function(x){
          xx <- round(x)
          
          # 5代表分为5个breaks
          # 10代表以10的倍数进行划分breaks
          steps <- if(xx[2] < 100){
            scaled_times <- ceiling(100/xx[2])
            
            floor((scaled_times * xx[2] / 5 ) / 10 ) * 10 / scaled_times
          }else{
            floor((xx[2] / 5 ) / 10 ) * 10
          }
          
          # 最小值都是0
          bk <- seq(0, xx[2], by = steps)
          
          # 第一间隔细分
          bk_minor <- seq(bk[1], bk[2], length.out = 5)
          bk <- c(bk_minor, bk[-(1:2)])
          
          # 最后间隔细分
          bk <- rev(bk)
          bk_minor <- seq(bk[1], bk[2], length.out = 5)
          bk <- c(bk_minor, bk[-(1:2)])
          bk <- rev(bk)
          
          bk
        }
      )
    p1
  })

pdf(file.path(qc_outdir, "Vlnplot-pass.pdf"), height = 5, width = 12)
mito_qc_ls_pass %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp))
    # cat("  \n\n")
  })
dev.off()

pdf(file.path(qc_outdir, "Vlnplot-log-pass.pdf"), height = 5, width = 12)
mito_qc_ls_pass %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp) & scale_y_log10())
    # cat("  \n\n")
  })
dev.off()







## cell filter table ------
cell_count_before <- all_data_ls      %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "Before filter")
cell_count_after  <- all_data_ls_pass %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "After filter")

# print 
cell_count_df <- cell_count_before %>% left_join(cell_count_after, by = "Sample")

# export
cell_count_df %>% 
  rownames_to_column(" ") %>%
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "cell count after filter.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# export data
saveRDS(all_data_ls_pass, file.path("data_for_analysis", "all_data_ls_pass.rds"))


# merge
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])
all_data_sp <- merge(all_data_ls[[1]], all_data_ls[-1])

# normalize ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
# all_data <- ScaleData(all_data, vars.to.regress = 'nCount_Spatial.008um')

# sketch data --------
all_data <-
  SketchData(
    all_data,
    ncells         = 20000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


## dimension reduction - pca ------------
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 2000)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')

all_data <- RunPCA(all_data)

pdf(file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_outdir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


# dimension reduction --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 50Gb


dims <- 13


all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - split.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap", split.by = "project", ncol = 4) + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

# cluster ------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  )

# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.5", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.7", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.7", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----
major_marker_ls <- c(
  "Epithelial"  = "Epcam",
  AT1 = "Ager",
  AT1 = "Pdpn",
  AT1 = "Csf1r",
  AT1 = "Mafb",
  AT1 = "Batf3",
  AT1 = "Krt79",
  AT1 = "Msr1",
  AT2 = "Abca3",
  AT2 = "Wif1",
  AT2 = "Etv5",
  "Basel cell" = "Ngfr",
  "B cell" = "Cd79a",
  "B cell" = "Cd79b",
  "B cell" = "Cd19",
  "B cell" = "Ms4a1",
  "T cells" = "Cd3d",
  "T cells" = "Cd3e",
  "T cells" = "Cd3g",
  "CD4T cell" = "Cd4",
  "CD4T cell" = "Cxcl13",
  "CD4T cell" = "Il21",
  "CD8T cell" = "Cd2",
  "CD8T cell" = "Cd8a",
  "CD8T cell" = "Cd8b1",
  "CD8T cell" = "Il2ra",
  "CD8T cell" = "Ptprc",
  Cilia = "Gsta1",
  Cilia = "Pifo",
  Cilia = "Foxj1",
  # "Club cell" = "PECAM1",
  Cilliated = "Muc5ac",
  "Dendritic cell" = "Ccl17",
  "Dendritic cell" = "Cd83",
  "Dendritic cell" = "Clec9a",
  "Dendritic cell" = "Ido1",
  Endothelial = "Pecam1",
  Endothelial = "Cdh5",
  Endothelial = "Vwf",
  Endothelial = "Plvap",
  Fibroblast = "Col1a1", 
  Fibroblast = "Col1a2",
  Fibroblast = "Pdgfra",
  Fibroblast = "Lum",
  Granulocyte = "Csf3r",
  "Mono/Mac" = "Cd68",
  "Mono/Mac" = "Fcgr1",
  "Mono/Mac" = "Cd163",
  "Mono/Mac" = "Cd14",
  "Mono/Mac" = "Fcgr3",
  "Mono/Mac" = "Vcan",
  "Mono/Mac" = "Mrc1",
  "Mast cell" = "Gata2",
  "Mast cell" = "Kit",
  "Mast cell" = "Ms4a2",
  "Mast cell" = "Tpsab1",
  "Mast cell" = "Tpsb2",
  "NK cell" = "Ncam1",
  "NK cell" = "Klrd1",
  Pericytes = "Pdgfrb",
  Pericytes = "Mcam",
  "Smooth muscle cells" = "Pln",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Tagln"
)

# fibro_markers %in% rownames(all_data)


# dotplot
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.5", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.7", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

# annotation ----
Idents(all_data) <- "sketch_snn_res.0.8"


cell_cluster <- c(
  '0' = "Mono/Mac", 
  '1' = "T cell", 
  '2' = "Fibroblast", 
  '3' = "Endothelial", 
  '4' = "B cell", 
  '5' = "AT2", 
  '6' = "AT1", 
  '7' = "Endothelial", 
  '8' = "Cilia", 
  '9' = "Cilia", 
  '10' = "AT2", 
  '11' = "Other", 
  '12' = "AT1", 
  '13' = "Mono/Mac", 
  '14' = "Dendritic cell", 
  '15' = "Mono/Mac", 
  '16' = "Pericyte", 
  '17' = "Fibroblast", 
  '18' = "SMC", 
  '19' = "Endothelial", 
  '20' = "Endothelial", 
  '21' = "SMC", 
  '22' = "Immune", 
  '23' = "AT2", 
  '24' = "Other", 
  '25' = "Other"
)

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(all_data) <- "sketch"
Idents(all_data) <- "annotation"
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


# find markers -------
Idents(all_data) <- "sketch_snn_res.0.8"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
all_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(mk_output, "marker.csv"), 
    row.names = FALSE
  )


all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  pull("gene")

all_de_genes_for_enrichmenet <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene,
    names = as.character(cluster)
  ))
all_top5_gene_unique <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <- 
  DotPlot(
    all_data, 
    features = all_top5_gene_unique
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(mk_output, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 25
)

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[[1]]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

avg_expr[all_de_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t() %>%
  pheatmap::pheatmap(
    filename = "4. marker gene/all_markers_heatmap.pdf", width = 20, height = 5,
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )


# project data ------
DefaultAssay(all_data) <- "sketch"
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",# "pca",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "RNA"

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

saveRDS(all_data, "data_for_analysis/all_data_projected.rds")
# all_data <- readRDS("data_for_analysis/all_data_projected.rds")

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_full", features = major_marker_ls, dot.scale = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

# 剔除混杂细胞
all_data_filt <- all_data[, !all_data@meta.data$annotation_full %in% c("Other")]

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected - trimmed.pdf")), height = 6, width = 8)
DimPlot(all_data_filt, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "full.harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected - trimmed.pdf")), height = 6, width = 8)
DimPlot(all_data_filt, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()


saveRDS(all_data_filt, 'data_for_analysis/all_data_projected_trimmed.rds')


# cell props barplot  --------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(sample, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation_full)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_col_proj_group <-
  all_data_meta %>%
  group_by(group, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, sample, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_full, y = freq, fill = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(dr_outdir, "cell_freq_project.pdf"),
  height = 6,
  width  = 10
)

ggsave(
  p_col_proj_group,
  filename = file.path(dr_outdir, "cell_freq_project_group.pdf"),
  height = 6,
  width  = 8
)

ggsave(
  p_box_group,
  filename = file.path(dr_outdir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)






# 细胞比例 全片 -------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

all_data_filt <- readRDS("data_for_analysis/all_data_projected_trimmed.rds")

dir.create("2. dimension reduction and cluster/Cell_ratio")

out_gene_umi = as.matrix(all_data_filt@meta.data)
row_name = c("Cell",rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file = "2. dimension reduction and cluster/Cell_ratio/meta_data.txt",sep = "\t",quote = F,row.names = F,col.names = F)


data <- read.table("2. dimension reduction and cluster/Cell_ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1, 7, 6, 19)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Unknown")


mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = "2. dimension reduction and cluster/Cell_ratio/Cell_type_stat.txt", quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i, ] <- round(mymtrix[i, ]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = "2. dimension reduction and cluster/Cell_ratio/Cell_type_stat_percent.txt", quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("2. dimension reduction and cluster/Cell_ratio/sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap) <- c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/Cell_ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = "2. dimension reduction and cluster/Cell_ratio/Cell_type_stat_group.txt", quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i, ] <- round(mymtrix[i, ]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = "2. dimension reduction and cluster/Cell_ratio/Cell_type_stat_percent_group.txt", quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7)+
  labs(x = "Sample", y = "Percent(%)", title = "")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7)+
  labs(x = "Group", y = "", title = "")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 15, colour = "black"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2)+
  plot_layout(widths = c(2, 1)) + # 两张图各自的宽度
  plot_annotation(
    title = "Cell Type Ratio",
    theme = theme(plot.title = element_text(
      size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
    ))
  )
ggsave(p, file = "2. dimension reduction and cluster/Cell_ratio/ratio_barplot.pdf", width = 12, height = 6)
ggsave(p, file = "2. dimension reduction and cluster/Cell_ratio/ratio_barplot.png", width = 12, height = 6, dpi = 300)









# 细胞比例 全片 -------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

all_data_filt <- readRDS("data_for_analysis/all_data_projected_trimmed.rds")

# 创建输出目录
dir.create("2. dimension reduction and cluster/Cell_ratio", recursive = TRUE, showWarnings = FALSE)

# 获取所有ROI
roi_list <- unique(all_data_filt$ROI)
print(paste("发现ROI:", paste(roi_list, collapse = ", ")))

# 为每个ROI创建子目录
for(roi in roi_list) {
  roi_dir <- file.path("2. dimension reduction and cluster/Cell_ratio", roi)
  dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)
  print(paste("创建目录:", roi_dir))
}

# 定义分析函数
analyze_roi_cell_ratio <- function(seurat_obj, roi_name, output_dir) {
  # 筛选当前ROI的细胞
  roi_data <- subset(seurat_obj, subset = ROI == roi_name)
  
  print(paste("分析ROI:", roi_name, "- 细胞数:", ncol(roi_data)))
  
  # 保存metadata
  out_gene_umi = as.matrix(roi_data@meta.data)
  row_name = c("Cell", rownames(out_gene_umi))
  out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
  out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
  write.table(out_gene_umi, 
              file = file.path(output_dir, "meta_data.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  # 读取并处理数据
  data <- read.table(file.path(output_dir, "meta_data.txt"), 
                     header = T, sep = "\t", check.names = F)
  colnames(data)
  data <- data[, c(1, 7, 6, 19)]  # 根据实际列索引调整
  colnames(data) <- c("Cells", "Sample", "Group", "Annotation")
  
  # 过滤未知细胞类型（如果存在）
  # data <- data %>% filter(Annotation != "Unknown")
  
  # 按样本和细胞类型计数
  mydata <- data %>% count(Sample, Annotation)
  mymtrix <- dcast(mydata, Sample ~ Annotation)
  mymtrix[is.na(mymtrix)] <- 0
  mymtrix_out <- t(mymtrix)
  write.table(as.data.frame(mymtrix_out), 
              file = file.path(output_dir, "Cell_type_stat.txt"), 
              quote = F, row.names = T, col.names = F, sep = "\t")
  
  # 计算百分比
  rownames(mymtrix) <- mymtrix[, "Sample"]
  mymtrix <- mymtrix[, -1]
  mymtrix <- as.matrix(mymtrix)
  percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
    mymtrix[i, ] <- round(mymtrix[i, ]/rowSums(mymtrix)[i]*100, 2)
  }))
  rownames(percent) <- rownames(mymtrix)
  data4plot1 <- melt(percent)
  percent <- t(percent)
  DATA <- cbind(rownames(percent), percent)
  DATA <- rbind(c("Sample", colnames(percent)), DATA)
  write.table(DATA, 
              file = file.path(output_dir, "Cell_type_stat_percent.txt"), 
              quote = F, row.names = F, col.names = F, sep = "\t")
  colnames(data4plot1) <- c("Sample", "Celltype", "Abund")
  
  # 按组分析
  final_data <- data
  mydata <- final_data %>% count(Group, Annotation)
  mymtrix <- dcast(mydata, Group ~ Annotation)
  mymtrix[is.na(mymtrix)] <- 0
  mymtrix_out <- t(mymtrix)
  write.table(as.data.frame(mymtrix_out), 
              file = file.path(output_dir, "Cell_type_stat_group.txt"), 
              quote = F, row.names = T, col.names = F, sep = "\t")
  
  rownames(mymtrix) <- mymtrix[, "Group"]
  mymtrix <- mymtrix[, -1]
  mymtrix <- as.matrix(mymtrix)
  percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
    mymtrix[i, ] <- round(mymtrix[i, ]/rowSums(mymtrix)[i]*100, 2)
  }))
  rownames(percent) <- rownames(mymtrix)
  data4plot2 <- melt(percent)
  percent <- t(percent)
  DATA <- cbind(rownames(percent), percent)
  DATA <- rbind(c("Group", colnames(percent)), DATA)
  write.table(DATA, 
              file = file.path(output_dir, "Cell_type_stat_percent_group.txt"), 
              quote = F, row.names = F, col.names = F, sep = "\t")
  colnames(data4plot2) <- c("Group", "Celltype", "Abund")
  
  # 设置颜色方案（如果plan未定义，使用默认颜色）
  if(!exists("plan")) {
    # 创建默认颜色方案
    n_celltypes <- length(unique(data4plot1$Celltype))
    plan <- scales::hue_pal()(n_celltypes)
    names(plan) <- unique(data4plot1$Celltype)
  }
  
  # 绘图
  plot_list <- list()
  plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
    geom_bar(stat = "identity", size = 0.5, width = 0.7) +
    labs(x = "Sample", y = "Percent(%)", title = paste("ROI:", roi_name)) +
    scale_fill_manual(values = plan) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 15, colour = "black"), 
      axis.title = element_text(size = 15, colour = "black", face = "bold"), 
      legend.title = element_text(size = 15, colour = "black", face = "bold"), 
      plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
      strip.background = element_rect(fill = "white", size = 2)
    ) + NoLegend()
  
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
      strip.background = element_rect(fill = "white", size = 2)
    )
  
  # 组合图形
  p <- wrap_plots(plot_list, ncol = 2) +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = paste("Cell Type Ratio -", roi_name),
      theme = theme(plot.title = element_text(
        size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
      ))
    )
  
  # 保存图形
  ggsave(p, file = file.path(output_dir, "ratio_barplot.pdf"), width = 12, height = 6)
  ggsave(p, file = file.path(output_dir, "ratio_barplot.png"), width = 12, height = 6, dpi = 300)
  
  print(paste("ROI", roi_name, "分析完成，结果保存在:", output_dir))
  
  return(list(
    data_sample = data4plot1,
    data_group = data4plot2,
    plot = p
  ))
}

# 对每个ROI执行分析
results <- list()

for(roi in roi_list) {
  roi_dir <- file.path("2. dimension reduction and cluster/Cell_ratio", roi)
  results[[roi]] <- analyze_roi_cell_ratio(all_data_filt, roi, roi_dir)
}

# 创建合并比较图（可选）
create_comparison_plot <- function(results, roi_list) {
  # 合并所有ROI的数据
  all_data <- data.frame()
  for(roi in roi_list) {
    roi_data <- results[[roi]]$data_group
    roi_data$ROI <- roi
    all_data <- rbind(all_data, roi_data)
  }
  
  # 创建比较图
  p_compare <- ggplot(all_data, aes(x = Group, y = Abund, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~ ROI, ncol = 1) +
    labs(x = "Group", y = "Proportion", title = "Cell Type Proportion Comparison Across ROIs") +
    scale_fill_manual(values = plan) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  ggsave(p_compare, 
         file = "2. dimension reduction and cluster/Cell_ratio/ROI_comparison.pdf", 
         width = 10, height = 8)
  ggsave(p_compare, 
         file = "2. dimension reduction and cluster/Cell_ratio/ROI_comparison.png", 
         width = 6, height = 8, dpi = 300)
  
  return(p_compare)
}

# 生成比较图
comparison_plot <- create_comparison_plot(results, roi_list)

# 保存汇总统计
summary_stats <- data.frame()
for(roi in roi_list) {
  roi_data <- subset(all_data_filt, subset = ROI == roi)
  stats <- data.frame(
    ROI = roi,
    Total_Cells = ncol(roi_data),
    Samples = length(unique(roi_data$sample)),
    CellTypes = length(unique(roi_data$annotation_full))
  )
  summary_stats <- rbind(summary_stats, stats)
}

write.csv(summary_stats, 
          "2. dimension reduction and cluster/Cell_ratio/ROI_summary_statistics.csv", 
          row.names = FALSE)











# spatial plot --------------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })




all_data_temp <- all_data


# 添加ROI信息 -------- 
roi_df <-
  list.files('ROI_2', full.names = TRUE, pattern = "csv") %>%
  # setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
  setNames(., str_extract(., "(?<=ROI_2/)[^.]+")) %>%
  map(~ read.csv(.x, col.names = c('Barcode', 'ROI'))) %>%
  enframe("group") %>%
  unnest(value) %>%
  # mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
  # mutate(roi_group = str_c(group, ROI, sep = "_")) %>%
  dplyr::select(Barcode, ROI) %>%
  column_to_rownames("Barcode")

all_data <- 
  all_data %>% 
  AddMetaData(roi_df)

all_data$ROI[is.na(all_data$ROI)] <- "Affected"
all_data_filt <- all_data[, !all_data@meta.data$annotation_full %in% c("Other")]

saveRDS(all_data, 'data_for_analysis/all_data_projected_ROI.rds')



# 单独画ROI Spatial图
all_data_Normal <- subset(all_data_filt, subset = ROI == "Normal")
all_data_Affected <- subset(all_data_filt, subset = ROI == "Affected")
all_data_More_affected <- subset(all_data_filt, subset = ROI == "More_affected")

# spatial plot Normal --------------
all_data_Normal_meta <-
  all_data_Normal[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_Normal_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output, str_glue("spatial plot {grp} - Normal.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })


# spatial plot Normal --------------
all_data_Affected_meta <-
  all_data_Affected[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_Affected_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output, str_glue("spatial plot {grp} - Affected.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })


# spatial plot Normal --------------
all_data_More_affected_meta <-
  all_data_More_affected[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_More_affected_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output, str_glue("spatial plot {grp} - More_affected.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })







# 基因表达图 ------- 
# spatial plot of some gene ------
sp_genes <- c("Adra2a", "Tnfrsf1a")

sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    
    dat <-
      cells_boundary %>%
      right_join(dat_expr, by = "cell_id")
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.0001, 0.9999))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = .data[[gene]], group = cell_id))  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
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
pdf(file.path(sp_output, str_glue("spatial plot - gene.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()







# gene expression plot
expression_data <- FetchData(
  object = all_data_filt,
  vars = c("Adra2a", "Tnfrsf1a", "group")
)

expression_data_long <- expression_data %>%
  pivot_longer(
    cols = c("Adra2a", "Tnfrsf1a"),
    names_to = "gene",
    values_to = "expression"
  )

expression_data_long$group <- factor(
  expression_data_long$group,
  levels = c("G1", "G2", "Dex_G1", "Dex_G2")
)

group_colors <- c(
  "G1" = "#fdc4b6",
  "G2" = "#e59572", 
  "Dex_G1" = "#96ceb4", 
  "Dex_G2" = "#248888"
  )


p <- ggplot(expression_data_long, aes(x = group, y = expression, fill = group)) +
  # 1. 先绘制箱线图
  geom_boxplot(
    alpha = 0.7, 
    outlier.shape = NA,  # 隐藏异常值
    width = 0.6,
    color = "black"      # 添加黑色边框
  ) +
  # 2. 在箱线图上方添加小点
  geom_point(
    position = position_jitter(width = 0.2),
    size = 0.3,          # 点的大小
    alpha = 0.2,         # 透明度
    shape = 16,          # 实心圆点
    color = "black"      # 点的颜色
  ) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    title = "Adra2a and Tnfrsf1a Expression Across Groups",
    x = "Groups",
    y = "Expression Level",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p)

ggsave(
  "Adra2a_Tnfrsf1a_expression_boxplot.png",
  plot = p, 
  width = 10,
  height = 6,
  dpi = 300
)












