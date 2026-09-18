library(tidyverse)
library(Seurat)
library(patchwork)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
# sp_dir <- "spatial plot" # dir.create(sp_dir)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)


list.files("/mnt/NexData_file/home/liuran/projects/3_Xenium/FWSP20250789_Mouse_Heart_Xenium/xenium_output_cellpose")

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
data_dir <- "/mnt/NexData_file/home/liuran/projects/3_Xenium/FWSP20250789_Mouse_Heart_Xenium/xenium_output_cellpose"
anno_file <- "data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "SAM000626 xenium-human-gene-expression-panel-gene-list.xlsx"


# reading data -------------
## data path and project name  -----
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## add id for parse_data_path --------
parse_data_path <- 
  parse_data_path %>%
  mutate(.id = seq_along(.proj))



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
if(file.exists(anno_file)){
  anno_data <- readr::read_csv(anno_file)
  User_columns <- colnames(anno_data)
}else{
  anno_data <- NULL
  User_columns <- "project"
}


## seurat object ------
all_data_ls <-
  parse_data_path %>% 
  dplyr::select(.counts_path, .proj, .id) %>%
  pmap(function(.counts_path, .proj, .id){
    # merge anno_data
    if(!is.null(anno_data)){
      row_idx <- which(anno_data$project == .proj)
      col_idx <- anno_data %>% colnames() 
      
      # 没有找到当前样本的注释数据
      #<<<<<< 所以Annotation.csv也有筛选文件的作用 >>>>>#
      if(length(row_idx) < 1){
        message(str_glue("Not found sample {.proj} in column: project of Annotation table.") )
        return(NULL)
      }
      # 没有任何数据，报错
      if(length(col_idx) < 1) stop("Not found any data in Annotation table.")
      
      meta_dat <- anno_data %>% .[row_idx, col_idx] %>% as.list()
    }else{
      meta_dat <- list(project = .proj)
    }
    
    # 构建Seurat
    sc <-
      Read10X(.counts_path) %>% 
      .[["Gene Expression"]] %>%
      CreateSeuratObject(
        project = .proj,
        min.cells    = 5,
        min.features = 5
      ) %>%
      AddMetaData(meta_dat) %>%
      RenameCells(add.cell.id = .id)
    sc
  }) 

all_data_ls <- all_data_ls %>% purrr::discard(is.null)
names(all_data_ls) <- map_chr(all_data_ls, "project.name")


# parse ROI----------------
# sample_df <-  
#   data.frame(
#     sample = c("kxz1","kxz6","kxz7","kcz2","kcz3","kcz8","gxz10","gxz11","gxz12","gcz13","gcz14","gcz15"),
#     group  = c("AirMale","AirMale","AirMale","AirFemale","AirFemale","AirFemale","HyperoxiaMale","HyperoxiaMale","HyperoxiaMale","HyperoxiaFemale","HyperoxiaFemale","HyperoxiaFemale")
#   )

# ### sample_df中的内容：sample为样本名，group为对应FOV名
# sample_df <-  
#   data.frame(
#     sample = c("kxz1","kxz6","kxz7","kcz2","kcz3","kcz8","gxz10","gxz11","gxz12","gcz13","gcz14","gcz15"),
#     group  = c("1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo","1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo","1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo","1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo","3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo","3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo","3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo","3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo","12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo","12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo","12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo","12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo")
#   )


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
# ### 根据每个样本所在芯片的顺序编号，来对应all_data_ls中的cell id
# .id<-c("1_","1_","1_","3_","3_","1_","2_","3_","3_","2_","2_","2_")
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
# # # 计算 roi_coords_ls 中每个样本的细胞数
# # cell_counts <- sapply(roi_coords_ls, function(x) length(x$cell))
# # 
# # # 输出每个样本的细胞数
# # cell_counts
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
# 
# # visium_roi override visium
# all_data_ls <- xenium_roi_ls



sample_df_fix <-
  data.frame(
    sample = c("D534","D536","D537","D538","D546","D547","s1","S2"),
    # group  = c()
    new_name  = c("IR7_1","IR7_2","IR3_1","IR3_2","IR1_1","IR1_2","S_1","S_2")
  )

all_data_ls_fix <- all_data_ls %>%
  imap(function(data, name) {
    # 大小写不敏感匹配：统一转换为小写
    idx <- which(tolower(sample_df_fix$sample) == tolower(name))
    
    if (length(idx) == 0) {
      warning(paste("No group found for sample:", name))
      return(data)
    } else if (length(idx) > 1) {
      warning(paste("Multiple matches for sample:", name, ". Using first entry."))
      idx <- idx[1]
    }
    
    data@meta.data$project <- sample_df_fix$new_name[idx]
    return(data)
  })


# 创建名称映射向量（大小写敏感匹配）
name_map <- setNames(sample_df_fix$new_name, sample_df_fix$sample)

# 替换列表名称
names(all_data_ls_fix) <- name_map[names(all_data_ls_fix)]

# 验证结果
names(all_data_ls_fix)

saveRDS(all_data_ls_fix, file.path("data_for_analysis", "all_data_ls_fix.rds"))

all_data_ls <- all_data_ls_fix



# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name  <- names(all_data_ls)


if(proj_count > 1){
  #merge
  all_data_ls[["Merge"]] <- mergeSeurat(all_data_ls)
}


## read coords -------------------
# cell coordinate and cell bundary
cells_position_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  map(function(x){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cells.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
    xx
  })

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

name_map <- setNames(sample_df_fix$new_name, tolower(sample_df_fix$sample))
names(cells_position_ls) <- name_map[tolower(names(cells_position_ls))]
names(cells_position_ls)

cells_position_ls <- cells_position_ls[proj_name]
cells_position_ls <- cells_position_ls %>% purrr::discard(is.null)
# cells_position_ls <- c(roi_position_ls,cells_position_ls)


cells_boundary_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  map(function(x){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cell_boundaries.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
    xx
  })

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

name_map <- setNames(sample_df_fix$new_name, tolower(sample_df_fix$sample))
names(cells_boundary_ls) <- name_map[tolower(names(cells_boundary_ls))]
names(cells_boundary_ls)

cells_boundary_ls <- cells_boundary_ls[proj_name]
cells_boundary_ls <- cells_boundary_ls %>% purrr::discard(is.null)
# cells_boundary_ls <- c(roi_boundary_ls,cells_boundary_ls)


# export data
dir.create("data_for_analysis") 
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))
saveRDS(all_data_ls, file.path("data_for_analysis", "all_data_ls.rds"))
# saveRDS(image_wh_r_ls, file.path("data_for_analysis", "image_wh_r_ls.rds"))




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


pdf(file.path("1. QC", "Vlnplot.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp))
    # cat("  \n\n")
  })
dev.off()

pdf(file.path("1. QC", "Vlnplot-log.pdf"), height = 5, width = 12)
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





pdf(file.path("1. QC", "qc plot - pass.pdf"), height = 12, width = 18)
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


## cell filter table ------
cell_count_before <- all_data_ls      %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "Before filter")
cell_count_after  <- all_data_ls_pass %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "After filter")

# print 
cell_count_df <- cell_count_before %>% left_join(cell_count_after, by = "Sample")

# export
cell_count_df %>% 
  rownames_to_column(" ") %>%
  openxlsx::write.xlsx(
    file = file.path("1. QC", "cell count after filter.xlsx"),
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
    ncells         = 10000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


## dimension reduction - pca ------------
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
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
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


dims <- 16


all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
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
major_marker_ls_1 <- c(
  "Adipo" = "Plin1",
  "Adipo" = "Adrb3",
  "B cell" = "Tshr",
  "B cell" = "Pax5",
  "B cell" = "Ms4a1",
  "Cardio" = "Grm1",
  "Cardio" = "Kcnn2",
  "Endocar" = "Plagl1",
  "Endocar" = "Nrg1",
  "Endocar" = "Plvap",
  "Endocar" = "Tmem108",
  "Endocar" = "Smoc1",
  "EC" = "Car4",
  "Epicar" = "Vtcn1",
  "Epicar" = "Slc26a3",
  "Fibro" = "Hmcn2",
  "Fibro" = "C7",
  "LEC" = "Ccl21a",
  "LEC" = "Reln",
  "LEC" = "Mmrn1",
  "LEC" = "Tbx1",
  "Mono" = "Flt3",
  "Mono" = "Xcr1",
  "Schwann" = "Gpr37l1",
  "Schwann" = "Gfra3",
  "SMC" = "Cdh6",
  "SMC" = "Myh11",
  "T cell" = "Themis",
  "T cell" = "Cd3e",
  "T cell" = "Itk",
  "Macro" = "F13a1",
  "Macro" = "Adgre1",
  "Peri" = "Trpc3",
  "Peri" = "Vtn"
)

major_marker_ls <- c(
  "EC" = "Cdh5",
  "EC" = "Pecam1",
  "FB" = "Col1a2",
  "FB" = "Col1a1",
  "FB" = "Pdgfra",
  "SMC" = "Tagln",
  "SMC" = "Myh11",
  "TC" = "Cd3e",
  "TC" = "Cd4",
  "TC" = "Cd8a",
  "BC" = "Cd79a",
  "BC" = "Ms4a1",
  "Cardiomyocytes" = "Tnni3k",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Pln",
  "Cardiomyocytes" = "Myh7",
  "Cardiomyocytes" = "Actn2",
  "MAC" = "Adgre1",
  "MAC" = "Fcgr1",
  "MAC" = "Csf1r",
  "Endocardial" = "Plagl1",
  "Endocardial" = "Nrg1",
  "Endocardial" = "Plvap",
  "Endocardial" = "Smoc1",
  "LECs" = "Ccl21a",
  "LECs" = "Ccl21b",
  "LECs" = "Mmrn1",
  "LECs" = "Reln",
  "LECs" = "Gfra3",
  "LECs" = "Prox1",
  "LECs" = "Lyve1",
  "Schwann cells" = "Sox10",
  "Schwann cells" = "Cnp",
  "Schwann cells" = "Cd59a",
  "Epicardial" = "Wt1",
  "Epicardial" = "Tbx18",
  "Epicardial" = "Aldh1a2",
  "Neutrophil" = "Csf3r",
  "Neutrophil" = "Cxcr2",
  "Neutrophil" = "Ccr1",
  "Pericytes" = "Pdgfrb",
  "Pericytes" = "Cspg4",
  "Pericytes" = "Mcam"
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
# annotation ---- res.0.3细胞注释
Idents(all_data) <- "sketch_snn_res.0.3"
cell_cluster <- c(
  '0'  = "Macro",
  '1'  = "Fibro",
  '2'  = "EC",
  '3'  = "Cardio",
  '4'  = "EC",
  '5'  = "Peri",
  '6'  = "Macro",
  '7'  = "Endocar",
  '8'  = "T cell",
  '9'  = "LEC/Schwann",
  '10' = "Epicar",
  '11' = "EC",
  '12' = "B cell"
)

# annotation ---- res.0.5细胞注释
Idents(all_data) <- "sketch_snn_res.0.5"
cell_cluster <- c(
  '0'  = "Cardiomyocytes",
  '1'  = "FB",
  '2'  = "Cardiomyocytes",
  '3'  = "EC",
  '4'  = "MAC",
  '5'  = "SMC/Pericytes",
  '6'  = "FB",
  '7'  = "MAC",
  '8'  = "MAC",
  '9'  = "Neutrophil",
  '10' = "Cardiomyocytes",
  '11' = "TC",
  '12' = "Endocardial",
  '13' = "LECs/Schwann cells",
  '14' = "Epicardial",
  '15' = "EC",
  '16' = "BC"
  )


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation_res05 <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(all_data) <- "sketch"
Idents(all_data) <- "annotation_res05"
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_res05", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


# find markers -------
Idents(all_data) <- "sketch_snn_res.0.5"
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

# top5 dotplot
p_dot_top5 <-
  DotPlot(
    all_data, 
    features = unique(all_top5_gene)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(mk_output, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 20
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
      annotation_res05_full = "annotation_res05",
      RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
      RNA_snn_res.0.3 = 'sketch_snn_res.0.3',
      RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
      RNA_snn_res.0.5 = 'sketch_snn_res.0.5',
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.7 = 'sketch_snn_res.0.7',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(all_data) <- 'annotation_res05_full'
DefaultAssay(all_data) <- "RNA"

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_res05_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

saveRDS(all_data, "data_for_analysis/all_data_res05_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_full", features = major_marker_ls, dot.scale = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

# # 发现T与plasma标错顺序了
# all_data$annotation_full2 <- all_data$annotation_full
# all_data$annotation_full2[all_data$annotation_full == "T"]      <- "Plasma"
# all_data$annotation_full2[all_data$annotation_full == "Plasma"] <- "T"
# all_data$annotation_full  <- all_data$annotation_full2

saveRDS(all_data, 'data_for_analysis/all_data_res05_annotation.rds')

# cell props barplot  --------

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation_full)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation_full) %>%
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
  p_box_group,
  filename = file.path(dr_outdir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)

# spatial plot --------------

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
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




# SMC/Pericytes subcluster ---------
sub_out_dir <- "2. dimension reduction and cluster/res.0.5/subcluster" # dir.create(fb_out_dir)
dir.create(sub_out_dir)

SMC_Pericyte <- all_data[, all_data$annotation_res05_full == "SMC/Pericytes"]
DefaultAssay(SMC_Pericyte) <- "RNA"
SMC_Pericyte[["sketch"]] <- NULL

SMC_Pericyte <- SMC_Pericyte %>% 
  FindVariableFeatures(nfeatures = 1000) %>% 
  NormalizeData() %>% 
  ScaleData(vars.to.regress = 'nCount_RNA')

## pca -----
SMC_Pericyte <- SMC_Pericyte %>% RunPCA()

pdf(file.path(sub_out_dir, 'pca plot-SMC_Pericyte.pdf'), height = 6, width = 8)
ElbowPlot(SMC_Pericyte, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_SMC_Pericyte <- 
  SMC_Pericyte[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(SMC_Pericyte) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(sub_out_dir, 'SMC_Pericyte dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(SMC_Pericyte[, suset_cells_SMC_Pericyte], dims = 1:5)
dev.off()
pdf(file.path(sub_out_dir, 'SMC_Pericyte dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(SMC_Pericyte[, suset_cells_SMC_Pericyte], dims = 6:10)
dev.off()
pdf(file.path(sub_out_dir, 'SMC_Pericyte dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(SMC_Pericyte[, suset_cells_SMC_Pericyte], dims = 11:15)
dev.off()
pdf(file.path(sub_out_dir, 'SMC_Pericyte dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(SMC_Pericyte[, suset_cells_SMC_Pericyte], dims = 16:20)
dev.off()


## dimension reduction --------
SMC_Pericyte_dims <- 15 #  20 # 20

# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

SMC_Pericyte <- RunUMAP(SMC_Pericyte, dims = 1:SMC_Pericyte_dims)

pdf(file.path(sub_out_dir, str_glue("SMC_Pericyte Dimplot dim{SMC_Pericyte_dims} project.pdf")), height = 6, width = 8)
DimPlot(SMC_Pericyte, group.by = "project", label = TRUE) + NoLegend()
dev.off()

SMC_Pericyte <- SMC_Pericyte %>% 
  harmony::RunHarmony(
    dims.use = 1:SMC_Pericyte_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

SMC_Pericyte <- RunUMAP(
  SMC_Pericyte, 
  dims = 1:SMC_Pericyte_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(sub_out_dir, str_glue("SMC_Pericyte Dimplot dim{SMC_Pericyte_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(SMC_Pericyte, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
SMC_Pericyte <-
  FindNeighbors(
    SMC_Pericyte,
    reduction = reduction.use,
    dims = 1:SMC_Pericyte_dims,
    k.param = 10
  )

SMC_Pericyte <-
  FindClusters(
    SMC_Pericyte,
    resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(sub_out_dir, str_glue("SMC_Pericyte Dimplot dim{SMC_Pericyte_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

## cell marker ----
major_marker_ls <- c(
  "SMC" = "Cdh6",
  "SMC" = "Myh11",
  "SMC" = "Tagln",
  "Peri" = "Trpc3",
  "Peri" = "Vtn",
  "Peri" = "Pdgfrb",
  "Peri" = "Cspg4",
  "Peri" = "Mcam"
)

major_marker_ls2 <- c(
  "EC" = "Cdh5",
  "EC" = "Pecam1",
  "FB" = "Col1a2",
  "FB" = "Col1a1",
  "FB" = "Pdgfra",
  "SMC" = "Tagln",
  "SMC" = "Myh11",
  "TC" = "Cd3e",
  "TC" = "Cd4",
  "TC" = "Cd8a",
  "BC" = "Cd79a",
  "BC" = "Ms4a1",
  "Cardiomyocytes" = "Tnni3k",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Pln",
  "Cardiomyocytes" = "Myh7",
  "Cardiomyocytes" = "Actn2",
  "MAC" = "Adgre1",
  "MAC" = "Fcgr1",
  "MAC" = "Csf1r",
  "Endocardial" = "Plagl1",
  "Endocardial" = "Nrg1",
  "Endocardial" = "Plvap",
  "Endocardial" = "Smoc1",
  "LECs" = "Ccl21a",
  "LECs" = "Ccl21b",
  "LECs" = "Mmrn1",
  "LECs" = "Reln",
  "LECs" = "Gfra3",
  "LECs" = "Prox1",
  "LECs" = "Lyve1",
  "Schwann cells" = "Sox10",
  "Schwann cells" = "Cnp",
  "Schwann cells" = "Cd59a",
  "Epicardial" = "Wt1",
  "Epicardial" = "Tbx18",
  "Epicardial" = "Aldh1a2",
  "Neutrophil" = "Csf3r",
  "Neutrophil" = "Cxcr2",
  "Neutrophil" = "Ccr1",
  "Pericytes" = "Pdgfrb",
  "Pericytes" = "Cspg4",
  "Pericytes" = "Mcam"
)

## dotplot -----
pdf(file.path(sub_out_dir, str_glue("dotplot dims{SMC_Pericyte_dims} {reduction.use} major_marker.pdf")), width = 10, height = 6)
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.05", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.1", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.2", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.3", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.4", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.6", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.8", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


pdf(file.path(sub_out_dir, str_glue("dotplot dims{SMC_Pericyte_dims} {reduction.use} major_marker2.pdf")), width = 20, height = 6)
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.05", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.1", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.2", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.3", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.4", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.6", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(SMC_Pericyte, group.by = "RNA_snn_res.0.8", features = major_marker_ls2) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


## umap plot of some gene ------
Idents(SMC_Pericyte) <- "RNA_snn_res.0.4"
markers <- read.table("2. dimension reduction and cluster/res.0.5/subcluster/marker.txt",header = T,sep="\t")
gene_list <- as.character(markers$gene)
celltype_list <- as.character(markers$celltype)
plot_list <- list()
for (i in (1:length(gene_list))){
  plot_list[[i]] <- FeaturePlot(SMC_Pericyte,
                                raster=FALSE, # 如果细胞数量比较多，可以将这个参数设置为FALSE
                                features = gene_list[i],
                                reduction = "harmonyumap", # 一般有两种选择：tsne和umap
                                pt.size = 0.05, # 点大小
                                label=T # 是否需要展示cluster标签
  )+
    theme_bw()+ #改变ggplot2的主题
    theme( #进一步修改主题
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
      axis.ticks = element_blank(),axis.title= element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
      legend.position = "none", #去掉图例
      plot.subtitle = element_text(hjust = 0.05,vjust=-9,size=15,color="red"),
      plot.title = element_text(vjust=-5,hjust = 0.5,size=17) #改变标题位置和字体大小
    )+
    labs(title =celltype_list[i], subtitle =gene_list[i])+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
}
p2 <- patchwork::wrap_plots(plot_list,ncol=4) # ncol参数可以调整出图时的列数
ggsave(p2, file = "2. dimension reduction and cluster/res.0.5/subcluster/Featureplot_all_RNA_snn_res.0.4.pdf", width = 22, height = 12) 


Idents(SMC_Pericyte) <- "RNA_snn_res.0.3"
markers <- read.table("2. dimension reduction and cluster/res.0.5/subcluster/marker.txt",header = T,sep="\t")
gene_list <- as.character(markers$gene)
celltype_list <- as.character(markers$celltype)
plot_list <- list()
for (i in (1:length(gene_list))){
  plot_list[[i]] <- FeaturePlot(SMC_Pericyte,
                                raster=FALSE, # 如果细胞数量比较多，可以将这个参数设置为FALSE
                                features = gene_list[i],
                                reduction = "harmonyumap", # 一般有两种选择：tsne和umap
                                pt.size = 0.05, # 点大小
                                label=T # 是否需要展示cluster标签
  )+
    theme_bw()+ #改变ggplot2的主题
    theme( #进一步修改主题
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
      axis.ticks = element_blank(),axis.title= element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
      legend.position = "none", #去掉图例
      plot.subtitle = element_text(hjust = 0.05,vjust=-9,size=15,color="red"),
      plot.title = element_text(vjust=-5,hjust = 0.5,size=17) #改变标题位置和字体大小
    )+
    labs(title =celltype_list[i], subtitle =gene_list[i])+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
}
p2 <- patchwork::wrap_plots(plot_list,ncol=4) # ncol参数可以调整出图时的列数
ggsave(p2, file = "2. dimension reduction and cluster/res.0.5/subcluster/Featureplot_all_RNA_snn_res.0.3.pdf", width = 22, height = 12) 


Idents(SMC_Pericyte) <- "RNA_snn_res.0.2"
markers <- read.table("2. dimension reduction and cluster/res.0.5/subcluster/marker.txt",header = T,sep="\t")
gene_list <- as.character(markers$gene)
celltype_list <- as.character(markers$celltype)
plot_list <- list()
for (i in (1:length(gene_list))){
  plot_list[[i]] <- FeaturePlot(SMC_Pericyte,
                                raster=FALSE, # 如果细胞数量比较多，可以将这个参数设置为FALSE
                                features = gene_list[i],
                                reduction = "harmonyumap", # 一般有两种选择：tsne和umap
                                pt.size = 0.05, # 点大小
                                label=T # 是否需要展示cluster标签
  )+
    theme_bw()+ #改变ggplot2的主题
    theme( #进一步修改主题
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
      axis.ticks = element_blank(),axis.title= element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
      legend.position = "none", #去掉图例
      plot.subtitle = element_text(hjust = 0.05,vjust=-9,size=15,color="red"),
      plot.title = element_text(vjust=-5,hjust = 0.5,size=17) #改变标题位置和字体大小
    )+
    labs(title =celltype_list[i], subtitle =gene_list[i])+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
}
p2 <- patchwork::wrap_plots(plot_list,ncol=4) # ncol参数可以调整出图时的列数
ggsave(p2, file = "2. dimension reduction and cluster/res.0.5/subcluster/Featureplot_all_RNA_snn_res.0.2.pdf", width = 22, height = 12) 



## annotation ----
Idents(SMC_Pericyte) <- "RNA_snn_res.0.3"
cell_cluster <- c(
  '0'  = "Peri",
  '1'  = "Peri",
  '2'  = "SMC",
  '3'  = "Peri",
  '4'  = "Peri",
  '5'  = "Peri",
  '6'  = "Peri",
  '7'  = "Peri",
  '8'  = "SMC",
  '9'  = "Peri"
)

SMC_Pericyte <- RenameIdents(SMC_Pericyte, cell_cluster)
SMC_Pericyte$SMC_Peri_annotation_res.0.3 <- Idents(SMC_Pericyte)


pdf(file.path(sub_out_dir, str_glue("Dimplot dim{SMC_Pericyte_dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(SMC_Pericyte, group.by = "SMC_Peri_annotation_res.0.3", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(SMC_Pericyte, group.by = "SMC_Peri_annotation_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(SMC_Pericyte) <- "RNA"
Idents(SMC_Pericyte) <- "SMC_Peri_annotation_res.0.3"
pdf(file.path(sub_out_dir, str_glue("dotplot dims{SMC_Pericyte_dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(SMC_Pericyte, group.by = "SMC_Peri_annotation_res.0.3", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()



## 替换SMC_Pericyte注释 -----
out_gene_umi=as.matrix(SMC_Pericyte@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="2. dimension reduction and cluster/res.0.5/subcluster/Filter_Basicinfo_nGene_nUMI_mito.txt",sep="\t",quote=F,row.names=F,col.names=F)

SMC_Pericyte_subset_metadata <- read.table("2. dimension reduction and cluster/res.0.5/subcluster/Filter_Basicinfo_nGene_nUMI_mito.txt", header=T, sep="\t", check.names=F)
colnames(SMC_Pericyte_subset_metadata)
SMC_Pericyte_subset_metadata <- SMC_Pericyte_subset_metadata[,c(1,43)]
colnames(SMC_Pericyte_subset_metadata) <- c("barcode","annotation_SMC_Pericyte")


all_data@meta.data$Cells <- rownames(all_data@meta.data)
all_data@meta.data$barcode <- all_data@meta.data$Cells
all_data@meta.data <- all_data@meta.data %>%
  left_join(SMC_Pericyte_subset_metadata, by = "barcode") %>%
  mutate(
    annotation_updated = case_when(
      # 在亚群对象中有细化的注释
      !is.na(annotation_SMC_Pericyte) ~ annotation_SMC_Pericyte,
      
      # 保持非T细胞类型的原始注释
      !annotation_res05_full %in% c("SMC/Pericytes") ~ annotation_res05_full,
      
      # 对于未细化的T细胞给予特殊标记
      TRUE ~ paste0(annotation_res05_full, " (unclassified)")
    )
  ) %>%
  select(-c(barcode, annotation_SMC_Pericyte))  # 移除临时列
unique(all_data$annotation_updated)
all_data@meta.data$cell_type_res05 <- all_data@meta.data$annotation_updated


dr_res0.5_output <- "2. dimension reduction and cluster/res.0.5"
Idents(all_data) <- "cell_type_res05"
rownames(all_data@meta.data) <- all_data@meta.data$Cells


all_data@meta.data$group <- factor(all_data@meta.data$group, levels = c("S", "IR1", "IR3", "IR7"))
pdf(file.path(dr_res0.5_output, str_glue("Full Dimplot dim{SMC_Pericyte_dims} cell_type_res05 {reduction.use} split full.harmonyumap.pdf")), height = 6, width = 22)
DimPlot(all_data, group.by = "cell_type_res05", label = TRUE, reduction = 'full.harmonyumap', ncol = 4, split.by = "group", pt.size = 0.1, raster = FALSE) + NoLegend()
dev.off()


DefaultAssay(all_data) <- "RNA"
Idents(all_data) <- "cell_type_res05"
pdf(file.path(dr_res0.5_output, str_glue("Full dotplot dims{SMC_Pericyte_dims} {reduction.use} annotation major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "cell_type_res05", features = major_marker_ls2, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

saveRDS(all_data, "data_for_analysis/all_data_res05_projected_subcluster.rds")



## barplot ------

library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")



# all_data_sp_2.dimension_reduction_anno <- readRDS("2.dimension_reduction/xenium_2.dimension_reduction.rds")
all_data <- readRDS("data_for_analysis/all_data_res05_projected.rds")

dir.create("2. dimension reduction and cluster/res.0.5/Ratio_new")

out_gene_umi=as.matrix(all_data@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="2. dimension reduction and cluster/res.0.5/Ratio_new/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)


data <- read.table("2. dimension reduction and cluster/res.0.5/Ratio_new/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,5,6,39)]
colnames(data) <- c("Cells","Sample","Group","Annotation")

# data <- data %>% filter(Annotation != "Unknown")


mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="2. dimension reduction and cluster/res.0.5/Ratio_new/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="2. dimension reduction and cluster/res.0.5/Ratio_new/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("2. dimension reduction and cluster/res.0.5/Ratio_new/sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio_new/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="2. dimension reduction and cluster/res.0.5/Ratio_new/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

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
write.table(DATA,file="2. dimension reduction and cluster/res.0.5/Ratio_new/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("2. dimension reduction and cluster/res.0.5/Ratio/sample_map_group.txt",header=F,sep="\t",quote="",check.names=F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[,1])
data4plot2$Group <- factor(data4plot2$Group, levels=unique(sample_map[,2]))

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
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
  # ggsci::scale_fill_igv() +
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
  plot_layout(widths = c(2, 1)) + # 两张图各自的宽度
  plot_annotation(
    title = "Cell Type Ratio",
    theme = theme(plot.title = element_text(
      size = 20, face = "bold", hjust = 0.5, margin = margin(b = 1, unit = "pt")
    ))
  )
ggsave(p,file="2. dimension reduction and cluster/res.0.5/Ratio_new/ratio_barplot.pdf",width=12, height=6)
ggsave(p,file="2. dimension reduction and cluster/res.0.5/Ratio_new/ratio_barplot.png",width=12, height=6,dpi=300)



# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = "2. dimension reduction and cluster/res.0.5/Ratio_new/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"
data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="2. dimension reduction and cluster/res.0.5/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "SampleID"))

sample<-read.table(file="2. dimension reduction and cluster/res.0.5/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by=c("Sample" = "SampleID"))

# 限制分组顺序
plot_data$Group <- factor(plot_data$Group, levels=unique(sample_map[,2]))

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="2. dimension reduction and cluster/res.0.5/Ratio_new/diff_boxplot.pdf",width=15,height=5)
ggsave(p,file="2. dimension reduction and cluster/res.0.5/Ratio_new/diff_boxplot.png",width=15,height=5)







## find markers -------
Idents(SMC_Pericyte) <- "RNA_snn_res.0.3"
DefaultAssay(SMC_Pericyte) <- "RNA"

all_markers <-
  FindAllMarkers(
    JoinLayers(SMC_Pericyte),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
all_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(sub_out_dir, "marker.csv"), 
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

# top5 dotplot
p_dot_top5 <-
  DotPlot(
    all_data, 
    features = all_top5_gene, 
    dot.scale = 10
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(sub_out_dir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 18
)


# spatial plot --------------
dir.create("2. dimension reduction and cluster/res.0.5/spatial plot")
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
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
      geom_polygon(aes(fill = cell_type_res05, group = cell_id))  + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = plan) + 
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






# Fibroblast subcluster ---------
fb_out_dir <- "fibroblast_subcluster" # dir.create(fb_out_dir)

sc_fibro <- all_data[, all_data$annotation_full == "Fibroblast"]
DefaultAssay(sc_fibro) <- "RNA"
sc_fibro[['sketch']] <- NULL


sc_fibro <-
  sc_fibro %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_fibro <- sc_fibro %>% RunPCA()

pdf(file.path(fb_out_dir, 'pca plot-fibro.pdf'), height = 6, width = 8)
ElbowPlot(sc_fibro, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_fibro <- 
  sc_fibro[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_fibro) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(fb_out_dir, 'fibro dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 1:5)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 6:10)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 11:15)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 16:20)
dev.off()


## dimension reduction --------
fibro_dims <- 10 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_fibro <- RunUMAP(sc_fibro, dims = 1:fibro_dims)


pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_fibro <- sc_fibro %>% 
  harmony::RunHarmony(
    dims.use = 1:fibro_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_fibro <- RunUMAP(
  sc_fibro, 
  dims = 1:fibro_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'pca' # 'harmony' 'pca'

## cluster ------
sc_fibro <-
  FindNeighbors(
    sc_fibro,
    reduction = reduction.use,
    dims = 1:fibro_dims,
    k.param = 10
  )

sc_fibro <-
  FindClusters(
    sc_fibro,
    resolution = c(0.05, 0.1, 0.2) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_fibro, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_fibro, group.by = "group2")
# dev.off()

## find markers -------
Idents(sc_fibro) <- "RNA_snn_res.0.2"
DefaultAssay(sc_fibro) <- "RNA"

# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_fibro),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
all_markers <-
  FindAllConservedMarkers(
    JoinLayers(sc_fibro),
    grouping.var = "project",
    test.use = "wilcox",
    # logfc.threshold = 0.25,
    # min.pct = 0.1,
    only.pos = TRUE
  ) %>%
  map(rownames_to_column, "gene") %>%
  enframe("cluster") %>%
  unnest(value) %>%
  rowwise() %>%
  mutate(
    avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
    pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
  ) %>%
  mutate(p_val_adj  = minimump_p_val)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_fibro))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(fb_out_dir, "all_marker_table.xlsx"),
    asTable = TRUE
  )

all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # arrange(p_val_adj, .by_group = TRUE) %>%
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

avg_expr[all_de_genes %>% unique, -(7:9)] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )

fibro_genes <- c(
  'VIM', 'THY', 'S100A4', 'DCN',
  'COL1A1', 'ACTA2', 'PDPN', 'FAP'
)

avg_expr[rownames(avg_expr) %in% fibro_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )

# saveRDS(sc_fibro, file = file.path(fb_out_dir, "sc_fibro.rds"))
## remove 6 7 8--------
sc_fibro2 <- sc_fibro[, !sc_fibro$RNA_snn_res.0.2 %in% 6:8]

pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
DimPlot(sc_fibro2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


# heatmap 7-9是小群，剔除
pdf(file.path(fb_out_dir, str_glue("fibro marker heatmapdim dim{fibro_dims} - {reduction.use} - trimmed.pdf")), height = 4, width = 8)
avg_expr[all_de_genes %>% unique %>% setdiff(c('CCDC80', 'TCF21')), -(7:9)] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )
dev.off()

# cell props 
sc_fibro_cell_props_df <-
  sc_fibro2[[]] %>% 
  dplyr::select(project, group, RNA_snn_res.0.2) %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq_project = count/sum(count)) %>%
  group_by(group) %>%
  mutate(freq_group = count/sum(count))

p_sc_fibro_props_project <-
  sc_fibro_cell_props_df %>%
  ggplot(aes(x = project, y = freq_project, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_bmj()

p_sc_fibro_props_group <-
  sc_fibro_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_bmj()

pdf(file.path(fb_out_dir, str_glue("fibro - cell props baplot - trimmed.pdf")), height = 4, width = 8)
p_sc_fibro_props_project %>% print()
p_sc_fibro_props_group   %>% print()
dev.off()


## dotplot ----------
FeaturePlot(
  sc_fibro,
  features = c(
    "ACTA2", "S100A4", "COL1A1", "VIM", 
    "PDGFRA", "PDGFRB", "THY1",
    "DCN", "FAP", "PDPN"
  ),
  reduction = "harmonyumap"
)


DotPlot(
  sc_fibro,
  features = c("ACTA2", "S100A4", "COL1A1", "VIM", "PDGFRA", "PDGFRB", "CD3D", "CD8A"),
  group.by = "RNA_snn_res.0.6"
)


DimPlot(sc_fibro, reduction = "harmonyumap", label = TRUE) + NoLegend() 



# testing analysis below -------------

# spatial plot of some gene ------
sp_genes <- c("Cdh5", "Pecam1", "Col1a2", "Col1a1", "Pdgfra", "Tagln", "Myh11", "Cd3e", "Cd4", "Cd8a", "Cd79a", "Ms4a1", "Tnni3k", "Ttn", "Pln", "Myh7", "Actn2", "Adgre1", "Fcgr1", "Csf1r", "Plagl1", "Nrg1", "Plvap", "Smoc1", "Ccl21a", "Ccl21b", "Mmrn1", "Reln", "Gfra3", "Prox1", "Lyve1", "Sox10", "Cnp", "Cd59a", "Wt1", "Tbx18", "Aldh1a2", "Csf3r", "Cxcr2", "Ccr1", "Pdgfrb", "Cspg4", "Mcam")
sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    all_data <- all_data_ls[[project]]
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
pdf(file.path(str_glue("spatial plot - gene.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()




sp_p_ls$FWSP20250913 %>% 
  iwalk(function(p, gene){
    pdf(file.path(str_glue("spatial plot - {gene}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


# some nFeature_RNA metric -------
all_data_ls %>%
  map(function(x){
    (x$nFeature_RNA > 1000) %>% table
  }) %>%
  bind_rows()


all_data_ls %>%
  map(function(x){
    (x$nFeature_RNA > 500) %>% table
  }) %>%
  bind_rows() %>%
  dplyr::rename(
    "< 500" = 1,
    "> 500" = 2
  )


all_data_ls %>%
  map(function(x){
    summary(x$nFeature_RNA)
  })  %>%
  enframe("sample") %>%
  mutate(value = bind_rows(value)) %>%
  unnest(value)


# Spatial plot by nFeature_RNA --------
sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    all_data_meta <- all_data_ls[[project]][[]] %>% rownames_to_column('cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = nFeature_RNA, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_gradient(low = "gray90", high = "red") +
      theme_void() +
      labs(
        title = project
      )
    p
  })

# export
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("1. QC", str_glue("spatial plot nFeature {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


# two samples analysis ---------
dir.create("two samples analysis")
all_data_ls %>% names

all_data_ls_susbet <- all_data_ls[c("X_202413076_32", "Y_2241905_12")]


all_data_ls_susbet <- 
  all_data_ls_susbet %>%
  map(function(x){
    subset(
      x,
      log10GenesPerUMI < 0.99 & log10GenesPerUMI > 0.9 
    )
  })

## normalize ------
all_data_ls_susbet <-
  all_data_ls_susbet %>%
  map(NormalizeData) %>%
  map(FindVariableFeatures) %>%
  map(ScaleData, vars.to.regress = 'nCount_RNA')



all_data_ls_susbet <-
  all_data_ls_susbet %>%
  map(RunPCA)


all_data_ls_susbet %>%
  iwalk(function(sc, name){# browser()
    
    pdf(file.path("two samples analysis", str_glue('pca plot {name}.pdf')), height = 6, width = 8)
    ElbowPlot(sc, ndims = 30) %>% print()
    dev.off()
    
    
    set.seed(1234)
    suset_cells <-
      sc[[]] %>%
      rownames_to_column("cell_id") %>%
      group_by(project) %>%
      slice_sample(prop = 50000/ncol(sc) ) %>%
      dplyr::pull("cell_id")
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 1-5 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 1:5)
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 6-10 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 6:10)
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 11-15 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 11:15)
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 16-20 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 16:20)
    dev.off()
  })

## dimension reduction -------
two_samples_dims <- c(
  X_202413076_32 = 10,
  Y_2241905_12   = 9
)


all_data_ls_susbet <- 
  all_data_ls_susbet %>%
  imap(function(sc, name){
    dims <- two_samples_dims[name]
    # sc <- RunUMAP(sc, dims = 1:dims)
    
    pdf(file.path("two samples analysis", str_glue("Dimplot dim{dims} project {name}.pdf")), height = 6, width = 8)
    print(DimPlot(sc, group.by = "project", label = TRUE, reduction = "umap") + NoLegend())
    dev.off()
    
    sc
  })


## cluster -------
all_data_ls_susbet <- 
  all_data_ls_susbet %>%
  imap(function(sc, name){
    reduction.use <- 'pca' # 'harmony' 'pca'
    
    dims <- two_samples_dims[name]
    
    # cluster  
    sc <-
      FindNeighbors(
        sc,
        reduction = reduction.use,
        dims = 1:dims,
        k.param = 20
      )
    
    sc <-
      FindClusters(
        sc,
        resolution = c(0.2, 0.4, 0.6, 0.8)
      )
    
    # dimplot
    pdf(file.path("two samples analysis", str_glue("Dimplot dim{dims} cluster - {reduction.use} {name}.pdf")), height = 6, width = 8)
    print(DimPlot(sc, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend())
    print(DimPlot(sc, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend())
    print(DimPlot(sc, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend())
    print(DimPlot(sc, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend())
    dev.off()
    
    sc
  })


## dotplot -------
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  # "Pericytes" = "PDGFRB",
  
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  
  Plasmablast = "MZB1",
  Plasmablast = "XBP1",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Fat cell' = 'ADIPOQ'
)


# dotplot
all_data_ls_susbet %>%
  iwalk(function(sc, name){
    reduction.use <- 'pca' # 'harmony' 'pca'
    dims <- two_samples_dims[name]
    
    pdf(file.path('two samples analysis', str_glue("dotplot dims{dims} major_marker - {reduction.use} {name}.pdf")), width = 25, height = 6)
    print(DotPlot(sc, group.by = "RNA_snn_res.0.2", features = major_marker_ls))
    print(DotPlot(sc, group.by = "RNA_snn_res.0.4", features = major_marker_ls))
    print(DotPlot(sc, group.by = "RNA_snn_res.0.6", features = major_marker_ls))
    print(DotPlot(sc, group.by = "RNA_snn_res.0.8", features = major_marker_ls))
    dev.off()
  })


## find markers ---------

all_markers_ls_susbet <-
  all_data_ls_susbet %>%
  map(function(sc){
    Idents(sc) <- "RNA_snn_res.0.2"
    all_markers <-
      FindAllMarkers(
        sc,
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        only.pos = TRUE
      )
    all_markers
  })



all_de_genes_ls_susbet <-
  all_markers_ls_susbet %>%
  map(function(all_markers){
    all_markers %>% 
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(pct.1 > 0.25) %>%
      arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = 10) %>%
      pull("gene")
  })

all_markers_ls_susbet %>%
  walk(function(all_markers){
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
  })


de_heatmap_ls <-
  all_data_ls_susbet %>%
  imap(function(sc, name){
    Idents(sc) <- "RNA_snn_res.0.2"
    
    de_genes <- all_de_genes_ls_susbet[[name]]
    
    avg_expr <- AverageExpression(sc)[['RNA']]
    
    p <-
      avg_expr[de_genes, ] %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
      pheatmap::pheatmap(
        scale = "column",
        cluster_row = avg_expr[VariableFeatures(sc), ] %>% t %>% dist %>% hclust,
        cluster_col = FALSE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        border_col = NA,
        silent = TRUE
      )
    cowplot::plot_grid(p$gtable)
  })


## annotation -------

annotation_ls <- list(
  "X_202413076_32" = c(
    '0'  = "Fibroblast",
    '1'  = "CD14 monocyte",
    '2'  = "Plasmablast",
    '3'  = "Endothelial",
    '4'  = "CD8T",
    '5'  = "Mast cell",
    '6'  = "Fat cell",
    '7'  = "Others"
  ),
  "Y_2241905_12"   = c(
    '0'  = "CD8T",
    '1'  = "Others",
    '2'  = "Plasmablast",
    '3'  = "Macrophage",
    '4'  = "B",
    '5'  = "CD14-dim monocyte",
    '6'  = "Epithelial",
    '7'  = "Endothelial",
    '8'  = "Fibroblast",
    '9'  = "Proliferating lymphocyte"
  )
)


# anno
all_data_ls_susbet <-
  all_data_ls_susbet %>%
  imap(function(sc, name){
    reduction.use <- "pca"
    cell_cluster <- annotation_ls[[name]]
    dims <- two_samples_dims[name]
    
    Idents(sc) <- "RNA_snn_res.0.2"
    sc <- RenameIdents(sc, cell_cluster)
    sc$annotation <- Idents(sc)
    
    
    
    pdf(file.path("two samples analysis", str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
    print(DimPlot(sc, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend())
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue("dotplot dims{dims} {reduction.use} annotation major_marker.pdf")) , width = 25, height = 6)
    print(DotPlot(sc, group.by = "annotation", features = major_marker_ls))
    dev.off()
    
    sc
  })

## spatial plot -----
sp_p_ls_subset <-
  cells_boundary_ls %>%
  .[names(all_data_ls_susbet)] %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    all_data_meta <- all_data_ls_susbet[[project]][[]] %>% rownames_to_column('cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      # scale_fill_gradient(low = "gray90", high = "red") +
      theme_void() +
      labs(
        title = project
      )
    p
  })

# export
sp_p_ls_subset %>% 
  iwalk(function(p, grp){
    pdf(file.path("two samples analysis", str_glue("spatial plot annotation {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


