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

# dimension reduction - pca ---------
DefaultAssay(all_data) <- "RNA"

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
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----
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
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker1)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker2)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_psa)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = top10_markers)+RotatedAxis()
dev.off()

# dotplot
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} major_marker_ls.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_psa)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_psa)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_psa)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_psa)
dev.off()


# annotation ----
Idents(all_data) <- "RNA_snn_res.0.6"

cell_cluster <- c(
  '0' = "luminal",#luminal上皮 FOXA1 FOLH1 AMACR
  '1' = "luminal",#肿瘤luminal ACPP NKX3-1
  '2' = "luminal",#肿瘤细胞​ EEF1G
  '3' = "luminal",#基底细胞 PROM2
  '4' = "Fibroblast",#肿瘤相关成纤维细胞
  '5' = "luminal",#​​间充质样肿瘤细胞​ EMT特征（GDF15、SPOCK1）
  '6' = "Macrophage",#巨噬细胞​,抗炎M2型（MRC1/CD206）
  '7' = "SMC",#平滑肌细胞​ MYLK、COL4A2
  '8' = "stroma",#WNT信号活化基质细胞​
  '9' = "Tcell",#CD3E、CD8A、CXCR4
  '10' = "Endothelial",#CD34、PECAM1
  '11' = "Pericyte",#PDGFRB
  '12' = "stroma",#​​衰老基质细胞​
  '13' = "luminal",#增殖肿瘤细胞​
  '14' = "luminal",#​​Club-like细胞​
  '15' = "luminal"#​​应激肿瘤细胞​
)
all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

saveRDS(all_data,'all_data_annotation.rds')

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
#source("h5ad2rds.R")
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


all_data2 <-  JoinLayers(all_data)
rownames(all_data2@assays$RNA@layers$counts) <- rownames(all_data2)
colnames(all_data2@assays$RNA@layers$counts) <- colnames(all_data2)
xenium <- CreateSeuratObject(
  all_data2@assays$RNA@layers$counts,
  assay = "Xenium")

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
                    dark.background = T,
                    border.color = "white", 
                    border.size = 0.1,
                    cols = all_colors8,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "annotation"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
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


#enrichment analysis-----------------
#debug(FindMarkers)
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

#再分群----
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
