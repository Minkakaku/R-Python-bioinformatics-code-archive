library(tidyverse)
library(Seurat)
library(patchwork)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
sp_dir <- "spatial plot" # dir.create(sp_dir)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)


# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
# cd /home/linrui/.cache/R/RBioFormats
# wget --no-check-certificate https://downloads.openmicroscopy.org/bio-formats/7.3.0/artifacts/bioformats_package.jar
# mv bioformats_package.jar bioformats_package_7.3.0.jar
library(RBioFormats)
source("xenium.R")

# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# 
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()


# setting variable   --------------------
data_dir <- "/data/xenium/FWSP20241023"
anno_file <- "data/Annotation.csv"
#ROI_anno_file <- "../data/ROI_Annotation.xlsx"


# reading data -------------
## data path and project name  ---
parse_data_path <-
  dir(data_dir, full.names = TRUE)%>% 
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## add id for parse_data_path ---
parse_data_path <- 
  parse_data_path %>%
  mutate(.id = seq_along(.proj))

## image w/h ratio ---
# used to fix width and heigth ratio in In Situ Spatial plot.
# image_wh_r_ls <-
#   with(parse_data_path, split(.path, .proj)) %>%
#   map(function(x){ #browser()
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


## read into annotation ---
if(file.exists(anno_file)){
  anno_data <- readr::read_csv(anno_file)
  User_columns <- colnames(anno_data)
}else{
  anno_data <- NULL
  User_columns <- "project"
}


# seurat object ---
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

# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name  <- names(all_data_ls)


# parse ROI-------
sample_df <- 
  data.frame(
    sample = c('LCR_C', 'LCR_T', 'YJ', 'JWL_C', 'JWL_T', 'ZSM', 'CHC', 'ZM', 
               'ZZF', 'C164960', 'T164960', 'T165854', 'GJJ_C', 'GJJ_T', 'LCY',
               'ZLZ', 'HQY', 'T414791', 'WXJ_C', 'WXJ_T', 'ysf_416971', 'ZSL',
               'WYH', 'ZYS'),
    group = c('FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023','FWSP20241023')
  )

roi_ls <- list.files(
  'ROI',
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("\\..*")

roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    dat <- read.csv(path,skip = 2)
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    dat<-dat %>% mutate(cell=cell.id)%>%
      column_to_rownames('cell.id')
    dat
  }) 

.id<-c('1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_', '1_')
names(.id)<-names(roi_coords_ls)
roi_coords_ls <- roi_coords_ls %>%
  map2(.id, ~ {
    .x$cell <- str_c(.y, .x$cell)
    return(.x)
  })

xenium_roi_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){ # browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- all_data_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    xenium_roi <- subset(xenium, cells = cells)
    xenium_roi
  })


# add sample_id or sample_name -------
xenium_roi_ls <-
  xenium_roi_ls %>%
  imap(function(xenium, name){ #browser()
    xenium$project <- name
    Idents(xenium) <- 'project'
    # RenameCells(visium, add.cell.id = name)
    xenium
  })

# visium_roi override visium
all_data_ls <- xenium_roi_ls



# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name  <- names(all_data_ls)


if(proj_count > 1){
  #merge
  all_data_ls[["Merge"]] <- mergeSeurat(all_data_ls)
}


# read coords ------------
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

roi_position_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){  #browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- cells_position_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    position_roi <-  xenium[which(xenium$cell_id %in% cells),]
    position_roi
  })

cells_position_ls <- cells_position_ls[proj_name]
cells_position_ls <- cells_position_ls %>% purrr::discard(is.null)
cells_position_ls <- c(roi_position_ls,cells_position_ls)


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

roi_boundary_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){  #browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- cells_boundary_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    boundary_roi <- xenium[which(xenium$cell_id %in% cells),]
    boundary_roi
  })

cells_boundary_ls <- cells_boundary_ls[proj_name]
cells_boundary_ls <- cells_boundary_ls %>% purrr::discard(is.null)
cells_boundary_ls <- c(roi_boundary_ls,cells_boundary_ls)


# export data
dir.create("data_for_analysis") 
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))
saveRDS(all_data_ls, file.path("data_for_analysis", "all_data_ls.rds"))
#saveRDS(image_wh_r_ls, file.path("data_for_analysis", "image_wh_r_ls.rds"))


#sample to group----
sample_df_fix <-
   data.frame(
      sample = c('LCR_C', 'LCR_T', 'YJ', 'JWL_C', 'JWL_T', 'ZSM', 'CHC', 'ZM', 
                 'ZZF', 'C164960', 'T164960', 'T165854', 'GJJ_C', 'GJJ_T', 'LCY',
                 'ZLZ', 'HQY', 'T414791', 'WXJ_C', 'WXJ_T', 'ysf_416971', 'ZSL',
                 'WYH', 'ZYS'),
      group = c('ADPC', 'ADPC', 'ADPC', 'ADPC', 'ADPC', 'CRPC', 'ADPC', 'ADPC',
                'CRPC', 'ADPC', 'ADPC', 'CRPC', 'ADPC', 'ADPC', 'CRPC', 'ADPC',
                'ADPC', 'CRPC', 'ADPC', 'ADPC', 'ADPC', 'BPH', 'BPH', 'BPH')
  )

all_data_ls_fix <- all_data_ls %>%
  imap(function(data, name) {
    # 获取正确的组名
    correct_group <- sample_df_fix$group[sample_df_fix$sample == name]
    correct_project <- sample_df_fix$sample[sample_df_fix$sample == name]
    # 检查是否找到正确的组名
    if (length(correct_group) == 0) {
      warning(paste("No group found for sample:", name))
      return(data)
    }
    
    # 添加或更新组名到元数据
    data@meta.data$group <- correct_group
    data@meta.data$project <- correct_project
    # 返回更新后的数据
    data
  })

# 验证组名是否更新正确
print(unique(all_data_ls_fix$C164960@meta.data$group))
print(unique(all_data_ls_fix$C164960@meta.data$project))

saveRDS(all_data_ls_fix, file.path("data_for_analysis", "all_data_ls_fix.rds"))

all_data_ls <- all_data_ls_fix

all_data_sp_rctd <- merge(all_data_ls[[1]], all_data_ls[-1])
saveRDS(all_data_sp_rctd, "data_for_analysis/all_data_sp_rctd.rds")


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
      #nCount_RNA         = paste0(">20"), # paste0(">20; <", nCount_cutoff),
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
  imap(qc_ls = parsed_mito_qc, function(sc, name, qc_ls){ # browser()
    
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
      map2(names(.), function(conds, var){ #  browser()
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


## cell filter table ------
cell_count_before <- all_data_ls      %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "Before filter")
cell_count_after  <- all_data_ls_pass %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "After filter")

# print 
cell_count_df <- cell_count_before %>% left_join(cell_count_after, by = "Sample")
saveRDS(all_data_ls_pass, file.path("data_for_analysis", "all_data_ls_pass.rds"))


# export
cell_count_df %>% 
  rownames_to_column(" ") %>%
  openxlsx::write.xlsx(
  file = file.path(qc_outdir, "cell count after filter.xlsx"),
  asTable = TRUE,
  withFilter = FALSE
)


# merge
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])


# normalize ------
all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))
all_data <- FindVariableFeatures(all_data)
# all_data <- ScaleData(all_data, vars.to.regress = 'nCount_Spatial.008um')

# sketch data ---
# all_data <-
#   SketchData(
#     all_data,
#     ncells         = 20000L,
#     sketched.assay = "sketch",
#     method         = "LeverageScore",
#     var.name       = "leverage.score",
#     seed = 123L
#   )


# dimension reduction - pca ---------
#dr_dir <- "dimension reduction" # dir.create(dr_dir)

DefaultAssay(all_data) <- "RNA"

# all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))

all_data <- FindVariableFeatures(all_data, nfeatures = 1500)
VariableFeaturePlot(all_data)

all_data <- ScaleData(all_data, vars.to.regress = c('nCount_RNA'))

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
dims <- 15


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE, n.neighbors = 30)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 8, width = 8)
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
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----
# 
# panel_gene_ls <- read.csv("XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv") %>%
#   dplyr::select(gene_name, cell_type) %>%
#   mutate(cell_type = cell_type %>% str_split(";") %>% map(str_trim)) %>%
#   unnest(cell_type) %>%
#   dplyr::filter(cell_type != "")  %>%
#   group_by(cell_type) %>%
#   nest() %>%
#   deframe() %>%
#   map(unlist, use.names = FALSE)
# 
# hepa_panel_gene_ls <-
#   panel_gene_ls[c(
#     "hepatocyte",
#     "fibroblast",
#     "hepatic stellate cell" ,
#     "midzonal region hepatocyte",
#     "centrilobular region hepatocyte",
#     "hepatoblast",
#     "endothelial cell of hepatic sinusoid",
#     "endothelial cell of pericentral hepatic sinusoid"
#   )]
# 


major_marker1 <- c(
  'epithelial cells'  = 'EPCAM',
  'epithelial cells'  = 'KRT5',
  'epithelial cells'  = 'KRT14',
  
  Hepatocyte          = 'ALB',
  Hepatocyte          = 'APOB',
  Hepatocyte          = 'CYP3A5',
  Hepatocyte          = 'CPB2',
  
  'endothelial cells' = 'PECAM1', 
  'endothelial cells' = 'VWF', 
  
  # 'fibroblasts'       = 'COL1A1', 
  # 'fibroblasts'       = 'DCN', 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  

  'T cells'           = 'CD3D', 
  'T cells'           = 'CD3E', 
  'T cells'           = 'CD3G', 
  'T cells'           = 'TRAC', 
    CD4='CD4',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  # "NK" = "GNLY",
  # "NK" = "NKG7",
  # "NK" = "NCAM1",
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  

  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  'mast cells'        = 'KIT',
  'mast cells'        = 'TPSAB1',
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC", 
  'Macrophage'  = 'CSF1R',
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  netrophil = "CSF3R",
  netrophil = "G0S2",
  netrophil = "RETNLG"
)



major_marker2 <- c(
  Epithelial = "EPCAM",
  Epithelial = "KRT8",
  Epithelial = "KRT5",
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "CLDN5",
  Endothelial = "AQP1",
  

  neron = 'TUBB3',
  neron = 'MAP2',
  neron = 'GAP43',
  neron = 'SNAP25',
  neron = 'SYT1',
  
  `Glutamatergic neurons` = "SLC17A7", 
  `GABAergic neurons` = "GAD2", 
  `GABAergic neurons` = "GAD1", 
  MSN = "SYT6", 
  MSN = "OTOF", 
  MSN = "DRD1", 
  MSN = "ADORA2A", 
  `Dopaminergic neurons` = "TH", 
  `Dopaminergic neurons` = "ALDH1A1", 
  `Serotonergic neurons` = "TPH1", 
  `Serotonergic neurons` = "TPH2", 
  `Cholinergic neurons` = "CHAT",
  UBC='EOMES',
  Purkinje='CA8',
  'Cere GC' = 'RBFOX3',
  
  MSC="NT5E",
  MSC="CD44",
  MSC="THY1",
  MSC="ENG",
  # Progenitor = 'SOX2',
  Progenitor = 'NOTCH1',
  Progenitor = 'FAM107A',
  Progenitor = 'NES',
  
  
  'Epen'='FOXJ1',
  'Epen'='DNAH11',
  'Epen'='TMEM212',
  'Epen'='CCDC153',
  
  Oligo = 'MOG',
  Oligo = 'SOX10',
  Oligo = 'OLIG1',
  Oligo = 'CLDN11',
  Oligo = 'MBP1',
  Astrocytes='AQP4',
  Astro = 'AGT',
  Astro = 'RORB',
  Astro = 'SLC7a10',
  Asrto = 'GJA1',
  Astro = 'GFAP',
  Astro = 'ALDH1L1l',
  Astro = 'SLC1A2',
  Astro = 'SLC1A3',
  RG='SOX2', 
  RG= 'FABP7',
  RG='BLBP',
  RG='GLAST',
  RG= 'PAX6',
  RG= 'HOPX',
  Immune='PTPRC',
  
  "Fibroblasts" = "FAP",
  "Fibroblasts" = "TNC",
  "Fibroblasts" = "PDPN",
  "Fibroblasts" = "PDGFRA",
  "Fibroblasts" = "PDGFRB",
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL4A2",
  "Fibroblasts" = "COL5A2",
  "Fibroblasts" = "STC2",
  "Fibroblasts" = "LOXL2",
  # "Fibroblasts" = "PDGFRA",
  # "Fibroblasts" = "LUM",
  # "Fibroblasts" = "DCN",
  # "Fibroblasts" = "VIM",
  # "Fibroblasts" = "COL1A2",
  
  
  # "Pericytes" = "PDGFRB",
  # Pericytes = "ACTA2",
  # Pericytes = 'RGS5',
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "IGHG1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHG3",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  
  Plasmablast = "MZB1",
  Plasmablast = "XBP1",
  
  "T cells" = "CD3D",
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "LCK",
  "T cells" = "TRAC",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  "T cells" = "GZMK",
  
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
  "Macrophage" = "F13A1",
  "Macrophage" = 'SIGLEC9',
  
  "Macrophage" = "MMP19",
  "Macrophage" = "CX3CR1",
  
  "Microglial" = "TMEM119",
  "Microglial" = "P2RY12",
  "Microglial" = "CTSS",
  "Monocytes" = "CD14", 
  "Monocytes" = "CD16",
  "Monocytes" = "FCGR3A",
  "Monocytes" = "S100A9",
  "Monocytes" = "FCN1",
  "Monocytes" = "CASP1",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",

  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  'Neutrophil' = "FCGR3B",
  'Neutrophil' = "IFITM2",
  'Neutrophil' = "CSF3R"
)

major_marker_psa <- c(
  'B cells' = 'CD19',
  'B cells' = 'CD20',
  'B cells' = 'CD79A',
  'B cells' = 'MS4A1',
  'B cells' = 'BLNK',
  'B cells' = 'CD79B',
  'Plasma' = 'JCHAIN',
  'Plasma' = 'IGJ',
  'Plasma' = 'MZB1',
  'Plasma' = 'TNFRSF17',
  'Plasma' = 'IGHG1',
  'Plasma' = 'IGKC',
  'T cells' = 'CD3E',
  'T cells' = 'CD3D',
  'T cells' = 'CD3G',
  'T cells' = 'CD4',
  'T cells' = 'FOXP3',
  'T cells' = 'CD8A',
  'T cells' = 'GZMB',
  'T cells' = 'TRAC',
  'Endothelial' = 'PECAM1',
  'Endothelial' = 'VWF',
  'Endothelial' = 'AQP1',
  'Endothelial' = 'CLEC4M',
  'Endothelial' = 'CLDN5',
  'Endothelial' = 'CDH5',
  'Endothelial' = 'ENG',
  'Endothelial' = 'SELE',
  'Epithelial' = 'EPCAM',
  'Epithelial' = 'KRT19',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'ACTA2',
  'Fibroblasts' = 'CCDC80',
  'Fibroblasts' = 'PDGFRA',
  'Fibroblasts' = 'SFRP4',
  'monocyte' = 'CD14',
  'monocyte' = 'FCGR3A',
  'monocyte' = 'FCN1',
  'monocyte' = 'VCAN',
  'monocyte' = 'LYZ',
  'macrophage' = 'CD74',
  'macrophage' = 'C1QC',
  'macrophage' = 'CD68',
  'macrophage' = 'CSF1R',
  'macrophage' = 'MACRO',
  'DC' = 'CD83',
  'DC' = 'LILRA4',
  'DC' = 'CLEC9A',
  'DC' = 'IDO1',
  'DC' = 'CD1C',
  'DC' = 'FCER1A',
  'NK' = 'GNLY',
  'NK' = 'NKG7',
  'NK' = 'NCAM1',
  'NK' = 'KLRD1',
  'NK' = 'KLRF1',
  'NK' = 'KLRK1',
  'NK' = 'KLRC1',
  'NK' = 'KLRB1',
  'Mast' = 'CPA3',
  'Mast' = 'CTSG',
  'Mast' = 'TPSAB1',
  'Luminal' = 'CD24A',
  'Luminal' = 'KRT8',
  'Luminal' = 'KRT18',
  'Luminal' = 'AR',
  'Immune' = 'PTPRC',
  'basal' = 'TP63',
  'basal' = 'KRT5',
  'basal' = 'KRT14',
  'basal' = 'KRT15',
  'progenitor' = 'CK19',
  'progenitor' = 'CK5',
  'progenitor' = 'ABCG2',
  'progenitor' = 'TCF1'
)



major_marker_ls <- c(
  'Epithelial cells'  = 'EPCAM',
  'Epithelial cells'  = 'KRT5',
  'Epithelial cells'  = 'KRT14',
  
  'Club cells' = 'SCGB1A3',
  #  'Club cells' = 'WFDC2',
  'Club cells' = 'LCN2',
  'Club cells' = 'MMP7',
  'Club cells' = 'KRT4',
  'Club cells' = 'TACSTD2',
  'Club cells' = 'SCGB3A1',
  'Club cells' = 'LTF',
  'Club cells' = 'RARRES1',
  'Club cells' = 'KRT7',
  'Club cells' = 'AGR2',
  
  'Hillock cells' = 'KRT13',
  'Hillock cells' = 'S100A16',
  'Hillock cells' = 'S100A14',
  'Hillock cells' = 'KRT19',
  'Hillock cells' = 'WFDC2',
  
  'Basal cells' = 'TP63',
  #  'Basal cells' = 'KRT14',
  #  'Basal cells' = 'KRT5',
  'Basal cells' = 'TRIM29',
  'Basal cells' = 'SLC14A1',
  'Basal cells' = 'KRT17',
  'Basal cells' = 'KRT15',
  
  'Luminal cells' = 'KRT8',
  'Luminal cells' = 'KRT18',
  'Luminal cells' = 'NKX3-1',
  'Luminal cells' = 'KLK4',
  'Luminal cells' = 'KLK3',
  'Luminal cells' = 'KLK2',
  'Luminal cells' = 'MSMB',
  'Luminal cells' = 'ACPP',
  
  'NE cells' = 'ASCL-1',
  'NE cells' = 'FoxA2',
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  'NE cells' = 'SYP',
  'NE cells' = 'NKX2-1',
  'NE cells' = 'ENO2',
  'NE cells' = 'MYCN',
  'NE cells' = 'POUF3F2',
  'NE cells' = 'LMO3',
  'NE cells' = 'INSM1',
  'NE cells' = 'SOX2',
  'NE cells' = 'EZH2',
  'NE cells' = 'SIAH2',
  
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2',
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "IGF1",
  "Fibroblasts" = "APOD",
  "Fibroblasts" = "FBLN1",
  "Fibroblasts" = "MEG3",
  "Fibroblasts" = "CXCL12",
  
  'B cells' = 'MS4A1',
  'B cells' = 'CD79A',
  'B cells' = 'CD19',
  'B cells' = 'BTG1',
  'B cells' = 'VPREB3',
  'B cells' = 'BANK1',
  'B cells' = 'IGLL5',
  'B cells' = 'CD79B',
  'B cells' = 'IGHM',
  'B cells' = 'CD22',
  
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  "Plasma cells" = "TNFRSF17",
  "Plasma cells" = "JCHAIN",
  "Plasma cells" = "MZB1",
  "Plasma cells" = "XBP1",
  'Plasma cells' = 'SEC11C',
  'Plasma cells' = 'PRDX4',
  'Plasma cells' = 'SPCS2',
  'Plasma cells' = 'SSR3',
  'Plasma cells' = 'SDF2L1',
  'Plasma cells' = 'C19orf10',
  'Plasma cells' = 'MANF',
  'Plasma cells' = 'TMEM258',
  'Plasma cells' = 'DNAJB9',
  'Plasma cells' = 'CD138',
  'Plasma cells' = 'CD38',
  'Plasma cells' = 'BCMA',
  'Plasma cells' = 'IRF4',
  'Plasma cells' = 'BLIMP1',
  
  'T cells' = 'CD3D', 
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'T cells' = 'TRAC', 
  'T cells' = 'NKG7',
  'T cells' = 'KLRC1',
  'T cells' = 'CCR7',
  'T cells' = 'FOXP3',
  'T cells' = 'CLTA4',
  'T cells' = 'CD8B',
  'T cells' = 'XCR6',
  'T cells' = 'IL32',
  'T cells' = 'CXCR6',
  
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC", 
  'Macrophages'  = 'CSF1R',
  'Macrophages' = 'C1QB',
  'Macrophages' = 'GPR34',
  'Macrophages' = 'MS4A4A',
  
  "NK cells" = "GNLY",
  #  "NK cells" = "NKG7",
  "NK cells" = "NCAM1",
  "NK cells" = "GZMB",
  "NK cells" = "KLRD1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  'Monocytes' = 'S100A9',
  'Monocytes' = 'FCN1',
  'Monocytes' = 'S100A8',
  'Monocytes' = 'CSTA',
  'Monocytes' = 'EREG',
  'Monocytes' = 'NEAT1',
  'Monocytes' = 'TKT',
  'Monocytes' = 'THBS1',
  'Monocytes' = 'VCAN',
  'Monocytes' = 'TSPO',
  'Monocytes' = 'LYZ',
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Netrophils' = "CSF3R",
  'Netrophils' = "G0S2",
  'Netrophils' = "RETNLG",
  
  'Mast cells' = 'KIT',
  'Mast cells'= 'TPSAB1',
  'Mast cells' = 'CPA3',
  'Mast cells' = 'HPGDS',
  'Mast cells' = 'RGS13',
  'Mast cells' = 'IL1RL1',
  'Mast cells' = 'GATA2',
  'Mast cells' = 'FOSB',
  'Mast cells' = 'RP11-354E11.2',
  'Mast cells' = 'VWA5A',
  'Mast cells' = 'SLC18A2',
  'Mast cells' = 'HDC',
  'Mast cells' = 'CAPG',
  'Mast cells' = 'FCER1G',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  'Pericytes' = 'ACTA2',
  #  'Pericytes' = 'MYH11',
  'Pericytes' = 'MT1M',
  'Pericytes' = 'FRZB',
  'Pericytes' = 'MT1A',
  'Pericytes' = 'NDUFA4L2',
  'Pericytes' = 'PPP1R14A',
  'Pericytes' = 'MYLK',
  'Pericytes' = 'PHLDA1',
  
  #  'Smooth muscle cells' = 'ACTA2',
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2',
  
  'Endothelials' = 'PECAM1', 
  'Endothelials' = 'ACKR1', 
  'Endothelials' = 'CLDN5',
  'Endothelials' = 'RAMP2',
  'Endothelials' = 'TM4SF1',
  'Endothelials' = 'RNASE1',
  'Endothelials' = 'EGFL7',
  'Endothelials' = 'RAMP3',
  'Endothelials' = 'PLVAP',
  'Endothelials' = 'AQP1',
  'Endothelials' = 'ECSCR',
  'Endothelials' = 'FKBP1A',
  'Endothelials' = 'AC011526.1',
  'Endothelials' = 'EMP1',
  'Endothelials' = 'DARC',
  'Endothelials' = 'VWF',
  'Endothelials' = 'EMCN'
)

# find markers -------
Idents(all_data) <- "RNA_snn_res.0.6"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

top10_markers <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 30) %>% 
  arrange(desc(pct.1)) %>% 
  slice_head(n=10) %>% 
  nest() %>% 
  deframe() %>% 
  map(pull,'gene') %>% 
  imap(function(genes,name){
    re <- genes %>% 
      map(function(gene){
        names(gene) = str_glue("C_{name}")
        gene
      })
    re
  }) %>% 
  unlist() %>% 
  .[!duplicated(.)]

names(top10_markers) <- sub('.*\\.','',names(top10_markers))



# dotplot
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 40, height = 6)
# DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker)
# DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker1)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker2)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_psa)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = top10_markers)+RotatedAxis()
# DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker)
dev.off()

# dotplot

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} major_marker_ls.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_psa)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_psa)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_psa)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_psa)
dev.off()

# dotplot
# pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} hepa_panel_marker.pdf")), width = 20, height = 6)
# hepa_panel_gene_ls %>%
#   iwalk(function(x, y){
#     # (DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
#     (DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
#     # (DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
#     # (DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
#   })  
# dev.off() 










# annotation ----
Idents(all_data) <- "RNA_snn_res.0.6"


cell_cluster <- c(
  '0' = "luminal",#luminal上皮 FOXA1 FOLH1 AMACR
  '1' = "luminal",#肿瘤luminal ACPP NKX3-1
  '2' = "luminal",#肿瘤细胞 EEF1G
  '3' = "luminal",#基底细胞 PROM2
  '4' = "Fibroblast",#肿瘤相关成纤维细胞
  '5' = "luminal",#间充质样肿瘤细胞 EMT特征（GDF15、SPOCK1）
  '6' = "Macrophage",#巨噬细胞,抗炎M2型（MRC1/CD206）
  '7' = "SMC",#平滑肌细胞 MYLK、COL4A2
  '8' = "stroma",#WNT信号活化基质细胞
  '9' = "T cell",#CD3E、CD8A、CXCR4
  '10' = "Endothelial",#CD34、PECAM1
  '11' = "Pericyte",#PDGFRB
  '12' = "stroma",#衰老基质细胞
  '13' = "luminal",#增殖肿瘤细胞
  '14' = "luminal",#Club-like细胞
  '15' = "luminal"#应激肿瘤细胞
)
all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data_rctd, group.by = "annotation", features = major_marker1)+ RotatedAxis()
# dev.off()

saveRDS(all_data,'all_data_annotation.rds')




#热图----
# valid_gene_list <- gene_list[gene_list %in% rownames(xenium)]
# cat("有效基因:", length(valid_gene_list), "/", length(gene_list))
# # 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
# missing_gene <- setdiff(gene_list, valid_gene_list)
# gene <- valid_gene_list
# 
# # 假设 immune_genes 是已分类的命名向量（如 'T_func' = "CD8A", ...）
# ordered_gene <- names(gene)  # 直接提取基因名（忽略分类名）
# # 或按分类组排序（保持组内顺序）：
# ordered_gene <- unname(gene)  # 仅保留基因名，顺序与定义一致
# 
# #xenium <- ScaleData(xenium, vars.to.regress = c('nCount_RNA'))
# H <- DoHeatmap(
#   object = xenium,
#   features = ordered_gene,   # 关键：传入有序基因向量
#   slot = "scale.data",
#   group.by = "annotation",
#   group.colors = all_colors8,
#   disp.min = -2.5,
#   disp.max = 2.5,
#   angle = -45,
#   hjust = 0.8,
#   size = 3,
#   raster = TRUE               # 加速渲染（大数据集必备）[1,6](@ref)
# )+
#   scale_fill_gradientn(
#     colors = c("navy", "white", "firebrick3"),
#     values = scales::rescale(c(-2.5, 0, 2.5))  # 精确匹配disp.min/max[4](@ref)
#   )
# 
# ggsave("result/MarkerHeatmap.pdf",plot = H,width = 12,height = 12,dpi = 300)



H1 <- DoHeatmap(
  object = all_data, 
  features = top10_markers,   # 关键：传入有序基因向量
  slot = "scale.data",    
  group.by = "annotation", 
  group.colors = all_colors8, 
  disp.min = -2.5, 
  disp.max = 2.5,
  angle = -45, 
  hjust = 0.8,
  size = 3,
  raster = TRUE               # 加速渲染（大数据集必备）[1,6](@ref)
)+
  scale_fill_gradientn(
    colors = c("navy", "white", "firebrick3"),
    values = scales::rescale(c(-2.5, 0, 2.5))  # 精确匹配disp.min/max[4](@ref)
  ) 

ggsave("MarkerHeatmap1.png",plot = H1,width = 12,height = 12,dpi = 300)


# # project data ------
# DefaultAssay(all_data) <- "sketch"
# all_data <- 
#   all_data %>% 
#   ProjectData(
#     assay = "RNA",
#     full.reduction = "harmony.full",
#     
#     sketched.assay = "sketch",
#     sketched.reduction = "harmony",# "pca",
#     
#     umap.model = "harmonyumap", 
#     dims = 1:dims, 
#     
#     refdata = list(
#       annotation_full = "annotation",
#       RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
#       RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
#       RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
#       RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
#     )
#   )
# 
# Idents(all_data) <- 'annotation_full'
# DefaultAssay(all_data) <- "RNA"
# 
# pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
# dev.off()
# 
# # saveRDS(all_data, "all_data_projected.rds")
# # all_data <- readRDS("all_data_projected.rds")
# 
# pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "annotation_full", features = major_marker) + RotatedAxis()
# dev.off()
# 
# ll<-subset(all_data, cells = colnames(all_data)[!is.na(all_data$annotation)])
# pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use} split1.pdf")), height = 6, width = 15)
# # DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(ll, group.by = "annotation", label = TRUE, reduction = 'harmonyumap',split.by = 'group') + NoLegend()
# dev.off()
# 

# cell props barplot  --------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(project, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, fill = group)) +
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
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(dr_outdir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)



#addFOV----
#source("/mnt/NexData_file/home/hyj/h5ad2rds(1).R")
library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library(ggplot2)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(arrow)


all_data2 <-  JoinLayers(all_data
                         #,layers = 'count'
                        )
rownames(all_data2@assays$RNA@layers$counts) <- rownames(all_data2)
colnames(all_data2@assays$RNA@layers$counts) <- colnames(all_data2)
xenium <- CreateSeuratObject(
  all_data2@assays$RNA@layers$counts,
  assay = "Xenium")

#saveRDS(xenium,"xenium_pro.rds")
#xenium <- readRDS("xenium_pro.rds")

xenium@meta.data <- all_data@meta.data
###改行名
rownames(xenium@meta.data) <- gsub("1_","",rownames(xenium@meta.data))
xenium@meta.data$cell_id <- rownames(xenium@meta.data)
xenium@meta.data$sample <- xenium@meta.data$project
rownames(xenium@meta.data) <- paste0(xenium@meta.data$sample,"_",xenium@meta.data$cell_id)

umap <- all_data@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)
xenium <- RenameCells(xenium,new.names = colnames(xenium))


has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
data.dir = "/data/xenium/FWSP20241023/output-XETG00305__0050485__FWSP20241023__20250112__042323"
xenium@meta.data$label = xenium@meta.data$sample
centroids = {
  if (has_dt) {
    centroids <- as.data.frame(data.table::fread(file.path(data.dir,
                                                           "cells.csv.gz")))
  } else {
    centroids <- read.csv(file.path(data.dir, "cells.csv.gz"))
  }
  centroids <- data.frame(x = centroids$x_centroid,
                          y = centroids$y_centroid, cell = centroids$cell_id,
                          stringsAsFactors = FALSE)}

segmentations = {
  if (has_dt) {
    segmentations <- as.data.frame(data.table::fread(file.path(data.dir,
                                                               "cell_boundaries.csv.gz")))
  } else {
    segmentations <- read.csv(file.path(data.dir,
                                        "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
  }
  segmentations <- data.frame(x = segmentations$vertex_x,y = segmentations$vertex_y,
                              cell = segmentations$cell_id)}
mols.qv.threshold = 20
# microns = {
#   if (has_dt) {
#     microns <- data.table::fread(file.path(data.dir, "transcripts.csv"))
#     microns <- subset(microns, qv >= mols.qv.threshold)
#   } else {
#     microns <- read.csv(file.path(data.dir, "transcripts.csv"))
#     microns <- subset(microns, qv >= mols.qv.threshold)
#   }
#   microns <- data.frame(x = microns$x_location, y = microns$y_location,
#                         cell = microns$cell_id,
#                         gene = microns$feature_name, stringsAsFactors = FALSE)
# }

microns = {
  microns <- read_parquet(file.path(data.dir, "transcripts.parquet"))
  microns <- subset(microns, qv >= mols.qv.threshold)
  microns <- data.frame(x = microns$x_location, y = microns$y_location,
                        cell = microns$cell_id,
                        gene = microns$feature_name, stringsAsFactors = FALSE)
}

object <- xenium
assay = "Xenium"
xenium_fov <- function(object,label,assay){
  cells_data <- data.frame(cell_id = object@meta.data$cell_id,sample = object@meta.data$sample)
  cells <- subset(cells_data,sample == label)$cell_id
  features <- rownames(object)
  
  segmentations_sub <- subset(segmentations,cell %in% cells)
  centroids_sub <- subset(centroids,cell %in% cells)
  microns_sub <- subset(microns,gene %in% features & cell %in% cells)
  
  segmentations_sub$cell <- paste0(label,"_",segmentations_sub$cell)
  centroids_sub$cell <- paste0(label,"_",centroids_sub$cell)
  microns_sub$cell <- paste0(label,"_",microns_sub$cell)
  
  segmentations.data <- list(centroids = CreateCentroids(centroids_sub),
                             segmentation = CreateSegmentation(segmentations_sub))
  
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation",
                                                            "centroids"), molecules = microns_sub, assay = assay)
  return(coords)
  
}

sample_label <- unique(xenium@meta.data$sample)
for (i in sample_label) {
  fov <- xenium_fov(object = xenium,label = i,assay = "Xenium")
  xenium@images[[i]] <- fov
}


saveRDS(xenium,"xenium_FOV_annotation.rds")
xenium <- readRDS("xenium_FOV_annotation.rds")


Idents(xenium) <- "annotation"
# rownames(xenium@assays$Xenium@layers$counts) <- rownames(xenium)
# colnames(xenium@assays$Xenium@layers$counts) <- colnames(xenium)

xenium <- NormalizeData(xenium)

###测试一个样本----
DefaultBoundary(xenium@images[["JWL_C"]]) <- "segmentation"
p <- ImageFeaturePlot(
  object = xenium,
  fov = "JWL_C",
  nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
  border.color = "white",  # 细胞边界的颜色（深蓝色）
  border.size = 0.1,    # 细胞边界的线宽
  features = 'EPCAM',
  max.cutoff = "q90",
  axes = FALSE,
  cols = c('#003366', "#E7475E"),
  size = 1.5,
  dark.background = T
)
ggsave(p,filename = paste0("sample","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


###Imagedimplot----
Imagedimplot_dir <- 'result/Imagedimplot'
dir.create(Imagedimplot_dir)


# library(RColorBrewer)
# n_groups <- length(unique(xenium$SubtypeT_Fib))  # 获取分组数量
# colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_groups)
#xenium@meta.data$Subtype <- as.factor(xenium@meta.data$Subtype)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928"
)

all_colors8 <- c(
  "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7"
)

sample_label <- unique(xenium@meta.data$sample)


for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.05,
                    cols = all_colors8,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "annotation"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
}


###SegImagefeature----
SegImagefeature_dir <- 'result/SegImagefeature'
dir.create(SegImagefeature_dir)

gene_list <- c('EPCAM')
valid_gene <- gene_list[gene_list %in% rownames(xenium)]
cat("有效基因:", length(valid_gene), "/", length(gene_list))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(gene_list, valid_gene)

gene_cus <- valid_gene

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_cus) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色
        border.size = 0.1,    # 细胞边界的线宽
        features = gene,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = TRUE
      )
      
      ggsave(
        filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}


#20250714-----
result_dir <- "~/xenium/result"
dir.create(result_dir)
###新分组----
# 构建样本名称到分组的映射向量
sample_to_group1 <- c(
  "LCR_C" = "PERI", "LCR_T" = "ADPC", "YJ" = "ADPC", "JWL_C" = "PERI", 
  "JWL_T" = "ADPC", "ZSM" = "CRPC", "CHC" = "ADPC", "ZM" = "ADPC", 
  "ZZF" = "CRPC", "C164960" = "PERI", "T164960" = "ADPC", "T165854" = "CRPC", 
  "GJJ_C" = "PERI", "GJJ_T" = "ADPC", "LCY" = "CRPC", "ZLZ" = "ADPC", 
  "HQY" = "ADPC", "T414791" = "CRPC", "WXJ_C" = "PERI", "WXJ_T" = "ADPC", 
  "ysf_416971" = "ADPC", "ZSL" = "BPH", "WYH" = "BPH", "ZYS" = "BPH"
)

# 假设您的Seurat对象名为xenium
current_samples <- xenium@meta.data$sample

# 通过样本名称匹配分组
xenium@meta.data$group1 <- sample_to_group1[current_samples]

# 验证是否所有样本均匹配成功
if (any(is.na(xenium@meta.data$group1))) {
  warning("部分样本未匹配到分组: ", 
          paste(unique(current_samples[is.na(xenium@meta.data$group1)]), collapse = ", "))
} else {
  print("分组添加成功，无未匹配样本")
}

###改细胞名称----
# 提取当前注释因子水平
current_levels <- levels(xenium@meta.data$annotation)

# 创建新旧名称映射关系（注意大小写匹配）
new_levels <- ifelse(
  current_levels == "Tcell", "T_cell",       # Tcell → T_cell
  ifelse(current_levels == "luminal", "Epithelial", # luminal → Epithelial
         ifelse(current_levels == "Macrophage", "Myeloid", # Macrophage → Myeloid
                current_levels)  # 其他类型保持不变
  ))

# 应用新因子水平
levels(xenium@meta.data$annotation) <- new_levels
unique(xenium@meta.data$annotation)

# 提取 all_data2 中的 harmonyumap 降维对象
harmony_umap_reduc <- all_data2@reductions[["harmonyumap"]]
# 深度复制对象（避免引用传递）
harmony_umap_copy <- harmony_umap_reduc  # 直接赋值在 S4 对象中为深度复制

# 将副本添加到 xenium 的 reductions 列表
xenium@reductions[["harmonyumap"]] <- harmony_umap_copy

# 强制对齐细胞名称
rownames(xenium@reductions$harmonyumap@cell.embeddings) <- colnames(xenium)

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 6)
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(xenium, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} annotation split1.pdf")), height = 7, width = 20)
DimPlot(xenium, group.by = "annotation",split.by = 'group1', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

#####cell props barplot  --------

all_data_meta <-
  xenium[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(sample, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_box_group  <-
  all_data_meta %>%
  group_by(group1, sample, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, fill = group1)) +
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
  filename = file.path(result_dir, "cell_freq_project.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(result_dir, "cell_freq_group1.pdf"),
  height = 6,
  width  = 10
)

#####循环细胞比例箱线图----
library(ggplot2)
library(dplyr)
library(ggpubr)

box_result_dir <- "~/xenium/result/boxplot"
dir.create(box_result_dir)
# 1. 计算每个样本中细胞类型的比例
cell_prop <- all_data_meta %>%
  count(sample, group1, annotation, name = "cell_count") %>%
  group_by(sample, group1) %>%
  mutate(total_cells = sum(cell_count)) %>%
  ungroup() %>%
  mutate(proportion = cell_count / total_cells)
# 2. 获取所有细胞类型列表
cell_types <- unique(all_data_meta$annotation)

# 3. 循环绘制每个细胞类型的箱线图
for (ct in cell_types) {
  # 筛选当前细胞类型数据
  plot_data <- cell_prop %>% 
    filter(annotation == ct) %>%
    arrange(group1)
  
  # 检查数据是否有效
  if (nrow(plot_data) < 3) {
    warning(sprintf("细胞类型 %s 的数据不足，跳过绘图", ct))
    next
  }
  
  # 设置绘图参数
  plot_title <- paste(ct)
  plot_file <- sprintf("cell_proportion_%s.png", gsub(" ", "_", ct))
  
  # 计算统计显著性 (PERI-ADPC-CRPC-BPH之间)
  stat_test <- compare_means(
    proportion ~ group1, 
    data = plot_data,
    method = "wilcox.test",
    ref.group = "PERI"  # 以PERI作为参照组
  )
  
  # 创建箱线图
  p <- ggplot(plot_data, aes(x = group1, y = proportion, fill = group1)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
    scale_fill_manual(values = c("PERI" = "#FFD700", "ADPC" = "#FF6347", 
                                 "CRPC" = "#9370DB", "BPH" = "#3CB371")) +
    labs(
      title = plot_title,
      x = "group",
      y = "prop",
      fill = "group"
    ) +
    theme_pubr() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # 添加统计显著性标记
  if (nrow(stat_test) > 0) {
    p <- p + stat_pvalue_manual(
      stat_test,
      label = "p.signif",
      y.position = max(plot_data$proportion) * seq(1.1, 1.4, length.out = nrow(stat_test))
    )
  }
  
  # 保存图像
  ggsave(
    filename = file.path(box_result_dir, plot_file),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message(sprintf("已生成 %s 的细胞比例箱线图", ct))
}

# saveRDS(xenium,"xenium_FOV_annotation.rds")
# xenium <- readRDS("xenium_FOV_annotation.rds")

###T细胞再分群----
All_T <- subset(xenium, annotation == 'T_cell')

#####PCA --------
DefaultAssay(All_T) <- "Xenium"
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')
All_T <- RunPCA(All_T)

T_out_dir<-'result/T_subtype' 
dir.create(T_out_dir)


#####pca inspector ------
pdf(file.path(T_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

pdf(file.path(T_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 1:5)
dev.off()
pdf(file.path(T_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 6:10)
dev.off()
pdf(file.path(T_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 11:15)
dev.off()
pdf(file.path(T_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(All_T, dims = 16:20)
dev.off()


#####dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


All_T <- RunUMAP(All_T, dims = 1:dims,return.model = TRUE)



pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


All_T <- All_T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

All_T <- RunUMAP(
  All_T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


#####cluster -------
All_T  <-
  FindNeighbors(
    All_T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T <-
  FindClusters(
    All_T,
    resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
  )



# dimplot
pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "Xenium_snn_res.0.05", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "Xenium_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


#####dotplot -------
T_major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='GZMB',
  'CD8'='GZMA',
  'CD8'='GZMK',
  'CD8'='GZMH',
  'CD8'='CCL5',
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCl5',
  'Naive'='CXCR3',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TRM'='CD89',
  'TRM'='CXCR8',
  # 'TEM'='IFNG',
  'TEM'='OCL4',
  'TEM'='OCL5',
  'TEM'='JAK3',
  # 'TCM'='IDO1',
  # 'TCM'='SKIV2L',
  # 'TCM'='TRAF2',
  # 'TCM'='SLC2A1',
  # 'TCM'='CD40LG',
  'TCM'='DUSP2',
  'TCM'='LYAR',
  'TCM'='ZNF683',
  'TCM'='IL32',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'Effector'='GNLY',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  # 'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  #'Exhausted'='PDCD1',
  'CD4'='FOXP3',
  'CD4'='CD4',
  'CD4'='TRDC',
  'CD4'='IL2RA',
  'CD4'='CD25',
  'CD4'='STAT1',
  'CD4'='IL2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  'Tfh'='PDCD1',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  
  
  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)

valid_Tgenes <- T_major_marker_ls[T_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_Tgenes), "/", length(T_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_Tgene <- setdiff(T_major_marker_ls, valid_Tgenes)
missing_Tgene


CD4_major_marker_ls <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  'T cells' = 'CD3D', 
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'T cells' = 'TRAC',
  'CD4'='CD4',
  # 'CD4'='TRDC',
  # 'CD4'='IL2RA',
  # 'CD4'='CD25',
  # 'CD4'='STAT1',
  # 'CD4'='IL2',
  'CD4_Treg'='TNFRSF4',
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='LAYN',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  'CD4_naive'='LEF1',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  'CD4_Th1'='ANXA1',
  'CD4_Th1'='CXCL13',
  'CD4_Th1'='BHLHE40',
  'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  'CD4_Tcm'='TCF7',
  'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2'
)
valid_CD4Tgenes <- CD4_major_marker_ls[CD4_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_CD4Tgenes), "/", length(CD4_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_CD4Tgene <- setdiff(CD4_major_marker_ls, valid_CD4Tgenes)
missing_CD4Tgene


CD8_major_marker_ls <- c(  
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAYN',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  'CD8_Tex'='ENPD1',
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='IFNG',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='NKG7',
  'CD8_Tnaive'='SELL',
  'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CCR7',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Tcm'='IL7R',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)
valid_CD8Tgenes <- CD8_major_marker_ls[CD8_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_CD8Tgenes), "/", length(CD8_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_CD8Tgene <- setdiff(CD8_major_marker_ls, valid_CD8Tgenes)
missing_CD8Tgene



# find markers
Idents(All_T) <- "Xenium_snn_res.0.6"
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

T_top10_markers <- 
  T_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(T_top10_markers) <- sub('.*\\.','',names(T_top10_markers))



#提取deepseek
T_top10_markers1 <-
  T_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)



# dotplot
pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(All_T, group.by = "Xenium_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "Xenium_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "Xenium_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "Xenium_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "Xenium_snn_res.0.6", features = T_top10_markers)+RotatedAxis()
DotPlot(All_T, group.by = "Xenium_snn_res.0.6", features = NK_major_marker_ls)+RotatedAxis()
dev.off()

NK_major_marker_ls <- c('CD8'='FCGR3A',
                        'CD8'='NCAM1'
)

valid_NKgenes <- NK_major_marker_ls[NK_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_NKgenes), "/", length(NK_major_marker_ls))


valid_Tgenes <- T_major_marker_ls[T_major_marker_ls %in% rownames(all_data)]
cat("有效基因:", length(valid_Tgenes), "/", length(T_major_marker_ls))

for (gene in valid_Tgenes) {
  p <- FeaturePlot(All_T, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(T_out_dir, paste0("feature/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

p <- FeaturePlot(All_T, features = "LAMP3", reduction = 'harmonyumap')+RotatedAxis()
ggsave(
  filename = file.path(T_out_dir, paste0("feature/","lamp3", "_feature.png")),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)


#DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = Fibro_major_marker_ls2)+RotatedAxis()

# Idents(All_T) <- "RNA_snn_res.0.8"
# All_T <- FindSubCluster(All_T, cluster = 0, graph.name = "RNA_snn", algorithm = 1,resolution = 0.6)
# Idents(All_T) <- "sub.cluster"
# # find markers
# T_markers_sub <-
#   FindAllMarkers(
#     JoinLayers(All_T),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
# 
# T_top5_markers_sub <- 
#   T_markers_sub %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 30) %>%
#   arrange(desc(pct.1)) %>%
#   slice_head(n = 5) %>%
#   nest() %>%
#   deframe() %>% 
#   map(pull, 'gene') %>%
#   imap(function(genes, name){
#     re <- genes %>%
#       map(function(gene){
#         # browser()
#         names(gene) = str_glue("cluster_{name}")
#         gene
#       })
#     re
#   }) %>%
#   unlist() %>% 
#   .[!duplicated(.)]
# 
# names(T_top5_markers_sub) <- sub('.*\\.','',names(T_top5_markers_sub))
# 
# DimPlot(All_T, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(All_T, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(All_T, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# 
# pdf(file.path(T_out_dir, str_glue("dotplot dims{dims} major_marker_sub.pdf")), width = 35, height = 6)
# DotPlot(All_T, group.by = "sub.cluster", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = CD4_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "sub.cluster", features = T_top5_markers_sub)+theme_90angle
# dev.off()

#####annotation -----------
cell_cluster2 <- c(
  '0' = "CD8",#CD8A
  '1' = "CD4",#
  '2' = "CD8",#CD8A
  '3' = "CD8",#CD8A
  '4' = "CD4",#
  '5' = "CD4",#
  '6' = "CD4",#
  '7' = "Plasma",#
  '8' = "B_cell",#
  '9' = "DC",#
  '10' = "Mast",#
  '11' = "ProT",#增殖性T细胞
  '12' = "DC",
  '13' = "CD4",
  '14' = "unknow"
)


cell_cluster3 <- c(
  '0' = "Tex", # CD8A GZMA KLRG1耗竭性CD8⁺T细胞
  '1' = "Tfh", # CD4 CD40LG KLRB1 滤泡辅助T细胞
  '2' = "TIMP2+CD8", #基质调节T细胞
  '3' = "Trm", # CD8A FOXA1 HSPD1 组织驻留记忆T细胞
  '4' = "Treg", # CD4 FOXP3 CTLA4#调节性T细胞
  '5' = "Th1",#CXCL9 CCL19 IFN-γ激活的Th1细胞
  '6' = "LIPA+CD4",#脂质代谢相关T细胞
  '7' = "Plasma",#
  '8' = "B_cell",#
  '9' = "LAMP3+DC",#
  '10' = "Mast",#
  '11' = "ProT",#PCNA RRM2增殖性T细胞
  '12' = "cDC",
  '13' = "CTL",#GZMB 细胞毒性T细胞
  '14' = "unknow"
)

Idents(All_T) <- "Xenium_snn_res.0.6"

All_T <- RenameIdents(All_T, cell_cluster3)
All_T$subtype3 <- Idents(All_T)

# cellType_marker_ls <- c(
#   'T' = 'CD3E',
#   
#   'CD8' = 'CD8A',
#   'CD8' = 'CD8B',
#   'CD8' = 'GZMH',
#   
#   'CD4' = 'CD4',
#   'CD4' = 'STAT1',
#   
#   'Fiblike_T' = 'COL5A1',
#   'Fiblike_T' = 'CCDC80',
#   
#   'MoT' = 'CD14',
#   'MoT' = 'CD68',
#   
#   "Epi" = 'EPCAM',
#   
#   'Pro_T' = 'MKI67',
#   'Pro_T' = 'TOP2A',
#   
#   'Treg'='FOXP3',
#   'Treg'='CTLA4',
#   
#   'Blike_T'='CD19',
#   'Blike_T'='CD79A'
# )

# pdf(file.path(T_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
# DotPlot(All_T, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
# dev.off()

pdf(file.path(T_out_dir,str_glue("Dimplot dim{dims} annotation3.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype3", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_out_dir, str_glue("Dimplot dim{dims} annotation3 split1.pdf")), height = 7, width = 18)
DimPlot(All_T, group.by = "subtype3",split.by = 'group1', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


#####cell props barplot  --------
All_T_meta <- All_T[[]]

All_T_meta$sample <- factor(All_T_meta$sample, 
                            levels=c(All_T_meta$sample %>% unique())
)
p_col_proj <-
  All_T_meta %>%
  group_by(sample, subtype3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype3)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/T_subtype/cell_freq3_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  All_T_meta %>%
  group_by(group1, subtype3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group1, y = freq, fill = subtype3)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path("result/T_subtype/cell_freq3_group.pdf"),
  height = 6,
  width  = 4
)


saveRDS(All_T, 'result/T_subtype/all_T.rds')

#####volcano-------
T_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(T_markers,'result/T_subtype/T_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(T_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- T_markers[T_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/T_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

T_markers_ls <- list(T_markers)
names(T_markers_ls) <- 'T cells'

# scatter
p_volcanoplot_major_ls <- 
  T_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/T_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


#####enrichment ----------------------------
T_de_genes_for_enrichmenet <- 
  T_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  T_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output<-'result/T_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

T_de_genes_id_for_enrichmenet <-
  T_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  T_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(All_T, 'result/T_subtype/all_T.rds')


###巨噬再分群----
# Myeloids/Macrophages
#all_data$annotation <- str_replace(all_data$annotation, 'Monocytes', 'Macrophages')

Macrophage_out_dir <- "result/Macrophage_subtype" 
dir.create(Macrophage_out_dir)
Macrophage <- subset(xenium, annotation == 'Myeloid')

#####PCA --------
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

#####pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

pdf(file.path(Macrophage_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 1:5)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 6:10)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 11:15)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 16:20)
dev.off()

#####dimension reduction --------
dims <- 16 

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


#####cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.05", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "Xenium_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Macrophage) <- "Xenium_snn_res.0.6"


#####dotplot -------
macrophage_major_marker_ls <- c(
  M1 = 'NOS2',
  M1 = 'IL12',
  M1 = 'FCGR1A',
  M1 = 'FCGR1B',
  M1 = 'FCGR1C',
  M1 = 'CD80',
  M1 = 'CXCR10',
  M1 = 'IL23',
  M1 = 'CXCL9',
  M1 = 'CXCL10',
  M1 = 'CXCL11',
  M1 = 'CD86',
  M1 = 'IL1A',
  M1 = 'IL1B',
  M1 =  'IL6',
  M1 =  'TNF',
  M1 = 'MHCII',
  M1 = 'CCL5',
  M1 = 'IRF5',
  M1 = 'IRF1',
  M1 = 'CD40',
  M1 = 'KYNU',
  M1 = 'CCR7',
  M1 =  'VCAN',
  M1 =  'CD300E',
  M1 =  'INHBA',
  M1 =  'TNIP3',
  M2 = 'ARG1',
  M2 = 'IL10',
  M2 = 'CD32',
  M2 = 'CD163',
  M2 = 'FCER2',
  M2 = 'CD200R1',
  M2 = 'PDCD1LG2',
  M2 = 'CD274',
  M2 = 'MARCO',
  M2 = 'CSF1R',
  M2 = 'MRC1',
  M2 = 'IL1RN',
  M2 = 'Il1R2',
  M2 = 'IL4R',
  M2 = 'CCL4',
  M2 = 'CCL13',
  M2 = 'CCL20',
  M2 = 'CCL17',
  M2 = 'CCL18',
  M2 = 'CCL22',
  M2 = 'CCL24',
  M2 = 'LYVE1',
  M2 = 'VEGFA',
  M2 = 'VEGFB',
  M2 = 'VEGFC',
  M2 = 'VEGFD',
  M2 = 'EGF',
  M2 = 'CTSA',
  M2 = 'CTSB',
  M2 = 'CSTC',
  M2 = 'CTSD',
  M2 = 'TGFB1',
  M2 = 'CD209',
  M2 = 'CD163L1',
  M2 = 'FOLR2',
  M0 = 'CD68',
  M0 = 'C1QA',
  M0 = 'C1QB',
  M0 = 'SELENOP',
  M0 = 'LYZ',
  M0 = 'HLA-DPB1',
  M0 = 'HLA-DPA1',
  M0 = 'AIF1', 
  TAM = 'SPP1', 
  TAM = 'CQC'
)


Myeloids_major_marker_ls <- c(
  # 单核细胞亚群
  'Classical Monocyte' = "CD14",
  'Classical Monocyte' = "S100A8",
  'Non-classical Monocyte' = "FCGR3A",  # CD16
  
  # 巨噬细胞亚群
  'M1-lMacrophage' = "CD86",
  'M1-Macrophage' = "CXCL9",
  'M2-Macrophage' = "CD163",
  'M2-Macrophage' = "MSR1",
  'Pro-angiogenic TAM' = "SPP1",       # 促血管生成型[5](@ref)
  'Anti-inflammatory TAM' = "C1QC",    # 抗炎型[5](@ref)
  'Anti-inflammatory TAM' = "APOE",
  
  # 树突状细胞亚群
  "Dendritic" = "FCER1A",
  "Dendritic" = "JCHAIN",
  'cDC1' = "CLEC9A",                   # 抗原呈递型[2,5](@ref)
  'cDC1' = "XCR1",
  'cDC2' = "CD1C",
  'LAMP3+ DC' = "LAMP3",               # 迁移型[5](@ref)
  
  # 粒细胞亚群
  'Neutrophil' = "CSF3R",
  'Neutrophil' = "CXCR2",
  'Eosinophil' = "SIGLEC8",
  'Basophil' = "CD123",
  
  # 肥大细胞
  'Mast' = "KIT",                 # 干细胞因子受体[2](@ref)
  'Mast' = "TPSAB1",              # 类胰蛋白酶[2](@ref)
  "Mast" = "MS4A2",
  "Mast" = "CPA3",
  
  # 髓系祖细胞
  'CMP' = "CD34",                      # 共同髓系祖细胞[1](@ref)
  'GMP' = "CD33"                       # 粒-单核祖细胞[4](@ref)
)

myeloids_cell_markers <- list(
  "Mast_KIT" = c("KIT", "TPSAB1", "CPA3"),
  "pDC_LILRA4" = c("LILRA4", "GZMB", "IL3RA"),
  "cDC1_CLEC9A" = c("CLEC9A", "FLT3", "IDO1"),
  "cDC2_CD1C" = c("CD1C", "FCER1A", "HLA-DQA1"),
  "cDC3_LAMP3" = c("LAMP3", "CCFR7", "FSCN1"),
  "Mono_CD14" = c("FCN1", "S100A8", "S100A9"),
  "Mono_CD16" = c("FCGR3A", "LST1", "LILRB2"),
  "Macro_INHBA" = c("INHBA","IL1RN", "CCL4"),
  "Macro_NLRP3" = c("NLRP3", "EREG", "IL1B"),
  "Macro_LYVE1" = c("LYVE1", "PLTP", "SEPP1"),
  "Macro_C1QC" = c("C1QC", "C1QA", "APOE")
)


Idents(Macrophage)<-'Xenium_snn_res.0.6'
macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

macrophage_top10_markers <- 
  macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(macrophage_top10_markers) <- sub('.*\\.','',names(macrophage_top10_markers))


#提取deepseek
macrophage_top10_markers1 <- 
  macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) 


# dotplot
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 27, height = 6)
DotPlot(Macrophage, group.by = "Xenium_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "Xenium_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "Xenium_snn_res.0.6", features = Myeloids_major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "Xenium_snn_res.0.6", features = macrophage_top10_markers)+RotatedAxis()
DotPlot(Macrophage, group.by = "Xenium_snn_res.0.6", features = myeloids_cell_markers)+RotatedAxis()
dev.off()

#feature
valid_genes <- major_marker_ls[major_marker_ls %in% rownames(xenium)]
cat("有效基因:", length(valid_genes), "/", length(major_marker_ls))

valid_mye_genes <- Myeloids_major_marker_ls[Myeloids_major_marker_ls %in% rownames(xenium)]
cat("有效基因:", length(valid_mye_genes), "/", length(Myeloids_major_marker_ls))

myeloids_cell_marker <- unique(unlist(myeloids_cell_markers))
valid_myel_genes <- myeloids_cell_marker[myeloids_cell_marker %in% rownames(xenium)]
cat("有效基因:", length(valid_myel_genes), "/", length(myeloids_cell_marker))

for (gene in valid_myel_genes) {
  p <- FeaturePlot(Macrophage, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(Macrophage_out_dir, paste0("feature3/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

#####score--------------------
M1<-c( 'NOS2',
       'IL12',
       'FCGR1A',
       'FCGR1B',
       'FCGR1C',
       'CD80',
       'CXCR10',
       'IL23',
       'CXCL9',
       'CXCL10',
       'CXCL11',
       'CD86',
       'IL1A',
       'IL1B',
       'IL6',
       'TNF',
       'MHCII',
       'CCL5',
       'IRF5',
       'IRF1',
       'CD40',
       'KYNU',
       'CCR7',
       'VCAN',
       'CD300E',
       'INHBA',
       'TNIP3')
M2<-c( 'ARG1',
       'IL10',
       'CD32',
       'CD163',
       'FCER2',
       'CD200R1',
       'PDCD1LG2',
       'CD274',
       'MARCO',
       'CSF1R',
       'MRC1',
       'IL1RN',
       'Il1R2',
       'IL4R',
       'CCL4',
       'CCL13',
       'CCL20',
       'CCL17',
       'CCL18',
       'CCL22',
       'CCL24',
       'LYVE1',
       'VEGFA',
       'VEGFB',
       'VEGFC',
       'VEGFD',
       'EGF',
       'CTSA',
       'CTSB',
       'CSTC',
       'CTSD',
       'TGFB1',
       'CD209',
       'CD163L1',
       'FOLR2')
Neutrophils <- c("FCGR3B","IFITM2","CSF3R")
Plasma <- c("TNFRSF17","JCHAIN","MZB1","XBP1","IGHG1","IGHA1","IGHG3","IGKC")
Tcells <- c("CD3D","CD3E","CD3G","LCK","TRAC","CD4","CD8A","GZMK")
Fibroblasts <- c("SFRP4","PDGFRA","CCDC80","COL1A1","ACTA2",'DCN',"FAP","TNC","PDPN","PDGFRB","COL4A2","COL5A2","STC2","LOXL2","LUM","VIM","COL1A2")
Bcells <- c("CD19","CD79A","MS4A1")
DC <- c("CD83","LILRA4","CLEC9A","IDO1","CD1C","FCER1A")

pathway<-list(M1,M2,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
names(pathway)<-c('M1','M2','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Macrophage <- AddModuleScore(JoinLayers(Macrophage), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(Macrophage,
           features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.6',
           pt.size = 0)

FeaturePlot(Macrophage,
            features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'))

ggsave(
  p,
  filename = 'result/Macrophage_subtype/M1_M2 score.pdf',
  height = 10,
  width  = 15
) 


#### annotation ----------------------------------------------------------------
cell_cluster2 <- c(
  '0' = 'Fibroblast', #PDGFRB, COL5A1, POSTN
  '1' = 'Epithelial', #EPCAM NKX3−1
  '2' = 'Macro', #CD163, MRC1, FOLR2 M2
  '3' = 'Macro', #CD68 M0
  '4' = 'DC',#CD1C,CLEC10A,FCER1A
  '5' = 'DC',#CD83,NR4A1
  '6' = 'unknow',#
  '7' = 'B_cell', # MS4A1 CD79A
  '8' = 'Epithelial',#FOLH1,SPON2 
  '9' = 'Macro',#CXCL9, STAT1 M1
  '10' = 'Macro',#CD68 CCL9 CXCL9 M1
  '11' = 'DC',#CLEC9A, XCR1, IRF8
  '12' = 'MAST' # KIT SIGLEC6
)


Idents(Macrophage) <- "Xenium_snn_res.0.6"

Macrophage <- RenameIdents(Macrophage, cell_cluster2)
Macrophage$subtype2 <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype2", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 20)
DimPlot(Macrophage, group.by = "subtype2",split.by = 'group1', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# Macrophage_cellTye_markers <- c( 
#   "Macrophage" = "CD68",
#   # "Macrophage" = "CD163",
#   "Macrophage" = "C1QA",
#   "Macrophage" = "C1QC", 
#   # 'Macrophage'  = 'CSF1R',
#   
#   "Fibroblasts" = "PDGFRA", 
#   "Fibroblasts" = "CCDC80",
#   "Fibroblasts" = "POSTN",
#   # "Fibroblasts" = "SFRP4",
#   "Fibroblasts" = "COL1A1",
#   "Fibroblasts" = "DCN",
#   "Fibroblasts" = "ACTA2",
#   
#   'vCAFs' = "MCAM", 
#   'vCAFs' = "NOTCH3", 
#   'vCAFs' = "COL18A1",
#   
#   'M2 macrophages' = "CD163",
#   'M2 macrophages' = "CSF1R",
#   'M2 macrophages' = "MRC1",
#   'M2 macrophages' = "CCL24",
#   'M2 macrophages' = "VEGFB",
#   'M2 macrophages' = "CTSA",
#   'M2 macrophages' = "CD209",
#   'M2 macrophages' = "FOLR2",
#   
#   'Epithelials' = "EPCAM",
#   
#   'Epithelial-like macrophages' = 'PGA3',
#   
#   'cDC2s' = 'CD1C',
#   'cDC2s' = 'CLEC10A',
#   
#   'Proliferating' = 'MKI67',
#   'Proliferating' = 'TOP2A',
#   
#   'LAMP3 DCs' = "LAMP3",
#   'LAMP3 DCs' = 'CD83', 
#   
#   'IFN-responsive macrophages' = "IFIT3", 
#   'IFN-responsive macrophages' = "IFIT2",
#   'IFN-responsive macrophages' = "IFIT1",
#   
#   'cDC1s' = 'CLEC9A',
#   'cDC1s' = 'IDO1'
# )

Macrophage_cellTye_markers <- c( 
  "M0-Mac" = "CD68",
  
  # "Macrophage" = "CD163",
  "M2-Mac" = "CD163",
  "M2-Mac" = "MRC1",
  # 'Macrophage'  = 'CSF1R',
  
  "Caf-like_MoMa" = "POSTN ",
  "Caf-like_MoMa" = "COL5A1 ",
  
  "cDC2" = "CD1C", 
  "cDC2" = "FCER1A",
  
  'Epi-like_MoMa' = "EPCAM",
  
  'M1-Mac' = "CXCL9",
  'M1-Mac' = "CXCL10",
  
  'Pro_MoMa' = 'MKI67',
  'Pro_MoMa' = 'TOP2A',
  
  'Mono' = "FCGBP",
  'Mono' = "CDC42",
  
  "MAST" = "KIT",
  "MAST" = "MS4A2",
  
  'Plasma-like_MoMa' = 'MZB1',
  'Plasma-like_MoMa' = 'XBP1',
  
  'LAMP3_DC' = "LAMP3",
  'LAMP3_DC' = 'FSCN1', 
  
  'pDC' = 'CLEC4C',
  'pDC' = 'LILRA4',
  
  'cDC1s' = 'CLEC9A',
  'cDC1s' = 'IDO1'
)

pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(Macrophage, group.by = "subtype", features = Macrophage_cellTye_markers, dot.scale = 10) 
#+ theme_90angle
dev.off()


#####cell props barplot  ----------------------------------------------------------
Macrophage_meta <- Macrophage[[]]

#Macrophage_meta$sample <- factor(Macrophage_meta$sample,levels=c(Macrophage_meta$sample %>% unique()))

p_col_proj <-
  Macrophage_meta %>%
  group_by(sample, subtype2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/Macrophage_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Macrophage_meta %>%
  group_by(group1, subtype2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group1, y = freq, fill = subtype2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Macrophage_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


saveRDS(Macrophage, 'result/Macrophage_subtype/Macrophage.rds')

# ###Only_Macrophage----
# # 生成子集对象（仅保留三种目标亚型）
# target_subtypes <- c('M0-Mac', 'M1-Mac', 'M2-Mac')
# Macrophage_subset <- subset(Macrophage, subtype %in% target_subtypes)
# 
# # 绘制分组UMAP图（三组对比）
# pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation filtered.pdf")), 
#     height = 6, width = 8)  
# DimPlot(
#   Macrophage_subset, 
#   group.by = "subtype",
#   #split.by = 'group3', 
#   reduction = 'harmonyumap',
#   label = TRUE,
#   repel = TRUE
# ) + NoLegend()
# dev.off()
# 
# # 绘制分组UMAP图
# pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation split3_filtered.pdf")), 
#     height = 6, width = 15)  
# DimPlot(
#   Macrophage_subset, 
#   group.by = "subtype",
#   split.by = 'group3', 
#   reduction = 'harmonyumap',
#   label = TRUE,
#   repel = TRUE
# ) + NoLegend()
# dev.off()
# 
# 
# Macrophage_subset_meta <- Macrophage_subset[[]]
# 
# Macrophage_subset_meta$sample <- factor(Macrophage_subset_meta$sample,levels=c(Macrophage_subset_meta$sample %>% unique()))
# 
# p_col_proj <-
#   Macrophage_subset_meta %>%
#   group_by(sample, subtype) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = sample, y = freq, fill = subtype)) + 
#   geom_col() +
#   ggsci::scale_fill_igv() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# 
# 
# # export
# ggsave(
#   p_col_proj,
#   filename = file.path("result/Macrophage_subtype", "Macro_cell_freq_sample.pdf"),
#   height = 6,
#   width  = 8
# )
# 
# p_col_proj <-
#   Macrophage_subset_meta %>%
#   group_by(group3, subtype) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = group3, y = freq, fill = subtype)) + 
#   geom_col() +
#   ggsci::scale_fill_igv() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# 
# ggsave(
#   p_col_proj,
#   filename = file.path(Macrophage_out_dir, "Macro_cell_freq_group3.pdf"),
#   height = 6,
#   width  = 4
# )
# 
# 



##### vocanol -------------------
Macrophage_markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Macrophage_markers,'Macrophage_subtype/Macrophage_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(Macrophage_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Macrophage_markers[Macrophage_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })

volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Macrophage_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


Macrophage_markers_major<-list(Macrophage = Macrophage_markers)
# scatter
p_volcanoplot_major_ls <- 
  Macrophage_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Macrophage_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 25
    )
  })


##### enrrichment ---------------------------------------
Macrophage_de_genes_for_enrichmenet <- 
  Macrophage_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

Macrophage_er_output<-'Macrophage_subtype/enrichment' # dir.create(Macrophage_er_output)
enrich_ls <-
  Macrophage_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Macrophage_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Macrophage_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Macrophage_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Macrophage_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Macrophage_de_genes_id_for_enrichmenet <-
  Macrophage_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Macrophage_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(Macrophage_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Macrophage_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Macrophage_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Macrophage_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Macrophage_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(Macrophage, 'result/Macrophage_subtype/Macrophage.rds')

Macrophage[['subtype3']] <- Macrophage[['subtype2']]



###映射——add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(All_T[['subtype3']] %>% rownames_to_column('cells'), 
                 #Fibroblast[['subtype']] %>% rownames_to_column('cells'),
                 Macrophage[['subtype3']]%>% rownames_to_column('cells'))

meta.data.t <-
  subtype%>%
  right_join(
    xenium[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtype3 = ifelse((as.character(annotation) %in% c('T_cell',"Myeloid")),
                     as.character(subtype3), 
                     as.character(annotation)))%>%
  dplyr::select(cells, subtype3) %>%
  column_to_rownames('cells')


xenium<-AddMetaData(xenium,meta.data.t['subtype3'])
xenium$annotation %>% table()
xenium$subtype3 %>% table()

# xenium$annotation2 <- NULL

# 确认是否加入
table(xenium$subtype3, xenium$annotation, useNA = "ifany")

saveRDS(xenium, 'result/xenium_subtype.rds')

#####cell props barplot  --------
all_data_meta <-
  xenium[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(sample, subtype3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype3)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_box_group  <-
  all_data_meta %>%
  group_by(group1, sample, subtype3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype3, y = freq, fill = group1)) +
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
  filename = file.path(result_dir, "sub3_freq_project.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(result_dir, "sub3_freq_group1.pdf"),
  height = 6,
  width  = 10
)

#####循环细胞比例箱线图----
library(ggplot2)
library(dplyr)
library(ggpubr)

box_result_dir <- "~/xenium/result/boxplot3"
dir.create(box_result_dir)
# 1. 计算每个样本中细胞类型的比例
cell_prop <- all_data_meta %>%
  count(sample, group1, subtype3, name = "cell_count") %>%
  group_by(sample, group1) %>%
  mutate(total_cells = sum(cell_count)) %>%
  ungroup() %>%
  mutate(proportion = cell_count / total_cells)
# 2. 获取所有细胞类型列表
cell_types <- unique(all_data_meta$subtype3)

# 3. 循环绘制每个细胞类型的箱线图
for (ct in cell_types) {
  # 筛选当前细胞类型数据
  plot_data <- cell_prop %>% 
    filter(subtype3 == ct) %>%
    arrange(group1)
  
  # 检查数据是否有效
  if (nrow(plot_data) < 3) {
    warning(sprintf("细胞类型 %s 的数据不足，跳过绘图", ct))
    next
  }
  
  # 设置绘图参数
  plot_title <- paste(ct)
  plot_file <- sprintf("cell_proportion_%s.png", gsub(" ", "_", ct))
  
  # 计算统计显著性 (PERI-ADPC-CRPC-BPH之间)
  stat_test <- compare_means(
    proportion ~ group1, 
    data = plot_data,
    method = "wilcox.test",
    ref.group = "PERI"  # 以PERI作为参照组
  )
  
  # 创建箱线图
  p <- ggplot(plot_data, aes(x = group1, y = proportion, fill = group1)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
    scale_fill_manual(values = c("PERI" = "#FFD700", "ADPC" = "#FF6347", 
                                 "CRPC" = "#9370DB", "BPH" = "#3CB371")) +
    labs(
      title = plot_title,
      x = "group",
      y = "prop",
      fill = "group"
    ) +
    theme_pubr() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  # 添加统计显著性标记
  if (nrow(stat_test) > 0) {
    p <- p + stat_pvalue_manual(
      stat_test,
      label = "p.signif",
      y.position = max(plot_data$proportion) * seq(1.1, 1.4, length.out = nrow(stat_test))
    )
  }
  
  # 保存图像
  ggsave(
    filename = file.path(box_result_dir, plot_file),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message(sprintf("已生成 %s 的细胞比例箱线图", ct))
}



#####Imagedimplot----
Imagedimplot_dir <- 'result/Imagedimplot3'
dir.create(Imagedimplot_dir)


# library(RColorBrewer)
# n_groups <- length(unique(xenium$SubtypeT_Fib))  # 获取分组数量
# colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_groups)
#xenium@meta.data$Subtype <- as.factor(xenium@meta.data$Subtype)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", "yellow", "red", "pink"
)

all_colors8 <- c(
  "#FF6B6B","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7"
)

sample_label <- unique(xenium@meta.data$sample)


for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.05,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "subtype3"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot3/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
}




#enrichment analysis-----------------
debug(FindMarkers)
Idents(all_data) <- "RNA_snn_res.0.4"
B_markers <-
  FindAllMarkers(
    JoinLayers(B),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(B_markers,'B_subtype/B_diff_gene.csv')
# all_markers1 <-
#   FindMarkers(
#     JoinLayers(all_data),
#     ident.1='Group2_CAFs',
#     ident.2='Group1_CAFs',
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
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


er_output<-'./CAFs_degene_enrichment'
library("org.Hs.eg.db")
library(clusterProfiler)
go<-enrichGO(rownames(all_markers), OrgDb =org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
         qvalueCutoff = 0.05,keyType = 'SYMBOL')
write.csv(
  go@result,
  file = file.path(er_output, "GO_enrichment_table.csv"),
  row.names = FALSE
)
dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
p<-enrichplot::dotplot(
  go,
  showCategory = 20,
  label_format = 50
)

ggsave(
  p,
  filename = file.path(er_output, "dotplot.pdf"),
  height = 6,
  width = 8
)

id<-clusterProfiler::bitr(
  rownames(all_markers),
  OrgDb = "org.Hs.eg.db",
  fromType = "SYMBOL",
  toType = "ENTREZID"
) %>%
  .$ENTREZID

library(createKEGGdb)
getOption("clusterProfiler.download.method")
R.utils::setOption("clusterProfiler.download.method","auto")
kegg <- enrichKEGG(
  id,
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = TRUE
)
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(
  kegg@result,
  file = file.path(er_output, "KEGG_enrichment_table.csv"),
  row.names = FALSE
)
p<-enrichplot::dotplot(
  kegg,
  showCategory = 20,
  label_format = 50
)

ggsave(
  p,
  filename = file.path(er_output, "KEGG_dotplot.pdf"),
  height = 6,
  width = 8
)



all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)


enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    dir.create(file.path(er_output, cluster_name))
    
    write.csv(
      x@result,
      file = file.path(er_output, cluster_name, "enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, cluster_name, "dotplot.pdf"),
      height = 6,
      width = 8
    )
  })
# # spatial plot --------------
# 
# all_data_meta <-
#   all_data[[]] %>%
#   rownames_to_column("cell_id")
# cells_boundary_ls<-readRDS(file.path("data_for_analysis", "cells_boundary_ls.rds"))
# cells_position_ls<-readRDS(file.path("data_for_analysis", "cells_position_ls.rds"))
# image_wh_r_ls<-readRDS(file.path("data_for_analysis", "image_wh_r_ls.rds"))
# 
# sp_p_ls1 <-
#   cells_boundary_ls %>%
#   imap(function(cells_boundary, project){ browser()
#     #image_wh_r <- image_wh_r_ls[[project]]
#     
#     dat <-
#       cells_boundary %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = vertex_x, y = vertex_y)) +
#       geom_polygon(aes(fill = annotation, group = cell_id))  +
#       scale_y_reverse() +
#       #coord_fixed(ratio = 1/image_wh_r) +
#       ggsci::scale_fill_igv() +
#       theme_void()
#     p
#   })
# 
# 
# # export image
# sp_dir <- "spatial plot" # dir.create(sp_dir)
# 
# sp_p_ls1 %>% 
#   iwalk(function(p, grp){
#     pdf(
#       file.path(sp_dir, str_glue("spatial plot {grp}.pdf")),
#       height = 12,
#       width  = 12
#     )
#     # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
#     print(p)
#     # cat("\n\n")
#     dev.off()
#   })


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


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


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
  arrange(desc(avg_log2FC)) %>%
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
  # arrange(p_val_adj) %>%
  arrange(desc(avg_log2FC)) %>%
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
marker<-c('ERBB2','ESR1','PGR','PD-L1','CD44','EPCAM','S100A8','SCD','TOP2A','STAT1','RIC8A','IRF9','ISG15','STAT2','JAK1','CDH1','CTNNA1','GNA13','OAS1')
sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation2, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })


sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    expression <- JoinLayers(all_data) %>% FetchData(marker,layer='data') %>% rownames_to_column('cell_id')
    all_data_meta <- all_data_meta %>% inner_join(expression, by = 'cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    marker <-  marker[which(marker %in% colnames(dat))]
      
    p1 <- marker %>% 
      map(function(gene){
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = .data[[gene]], group = cell_id))  +
          scale_y_reverse() +
          coord_fixed(ratio = 1/image_wh_r) +
          # ggsci::scale_fill_igv() +
          scale_fill_gradient(low = "gray90", high = "red") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p1
  })

sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){  browser()
    image_wh_r <- image_wh_r_ls[[project]]
    expression <- JoinLayers(all_data) %>% FetchData(marker,layer='data') %>% rownames_to_column('cell_id')
    all_data_meta <- all_data_meta %>% inner_join(expression, by = 'cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = marker, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
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
    pdf(file.path("spatial plot", str_glue("spatial plot nFeature {grp}.pdf")), height = 18, width = 40)
    print(patchwork::wrap_plots(p,nrow =4, ncol =4))
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
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 10) %>%
      pull("gene")
  })

all_markers_ls_susbet %>%
  walk(function(all_markers){
    all_markers %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
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


