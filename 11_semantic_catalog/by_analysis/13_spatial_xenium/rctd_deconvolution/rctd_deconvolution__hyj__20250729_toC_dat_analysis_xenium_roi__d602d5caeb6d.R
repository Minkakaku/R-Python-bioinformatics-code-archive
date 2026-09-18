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




# setting variable   --------------------
data_dir <- "/data/xenium/FWSP20241023"
anno_file <- "data/Annotation.csv"


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
  nmols = 20000,       # 每个图像中显示的最大分子数
  border.color = "white",  # 细胞边界的颜色
  border.size = 0.1,    # 细胞边界
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
        nmols = 20000,       # 每个图像中显示的最大分子数
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
sample_to_group1 <- c(
  "LCR_C" = "PERI", "LCR_T" = "ADPC", "YJ" = "ADPC", "JWL_C" = "PERI", 
  "JWL_T" = "ADPC", "ZSM" = "CRPC", "CHC" = "ADPC", "ZM" = "ADPC", 
  "ZZF" = "CRPC", "C164960" = "PERI", "T164960" = "ADPC", "T165854" = "CRPC", 
  "GJJ_C" = "PERI", "GJJ_T" = "ADPC", "LCY" = "CRPC", "ZLZ" = "ADPC", 
  "HQY" = "ADPC", "T414791" = "CRPC", "WXJ_C" = "PERI", "WXJ_T" = "ADPC", 
  "ysf_416971" = "ADPC", "ZSL" = "BPH", "WYH" = "BPH", "ZYS" = "BPH"
)


current_samples <- xenium@meta.data$sample
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
  current_levels == "Tcell", "T_cell",  
  ifelse(current_levels == "luminal", "Epithelial", 
         ifelse(current_levels == "Macrophage", "Myeloid", 
                current_levels)  # 其他类型保持不变
  ))

# 应用新因子水平
levels(xenium@meta.data$annotation) <- new_levels
unique(xenium@meta.data$annotation)

# 提取 all_data2 中的 harmonyumap 降维对象
harmony_umap_reduc <- all_data2@reductions[["harmonyumap"]]
harmony_umap_copy <- harmony_umap_reduc  
xenium@reductions[["harmonyumap"]] <- harmony_umap_copy
rownames(xenium@reductions$harmonyumap@cell.embeddings) <- colnames(xenium)

pdf(file.path(result_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 6)
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

#####细胞比例箱线图----
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

cell_types <- unique(all_data_meta$annotation)

# 3. 循环绘制每个细胞类型的箱线图
for (ct in cell_types) {
 
  plot_data <- cell_prop %>% 
    filter(annotation == ct) %>%
    arrange(group1)
  
  if (nrow(plot_data) < 3) {
    warning(sprintf("细胞类型 %s 的数据不足，跳过绘图", ct))
    next
  }
  
  plot_title <- paste(ct)
  plot_file <- sprintf("cell_proportion_%s.png", gsub(" ", "_", ct))
  
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
  
  ggsave(
    filename = file.path(box_result_dir, plot_file),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  message(sprintf("已生成 %s ", ct))
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

NK_major_marker_ls <- c('NK'='FCGR3A',
                        'NK'='NCAM1'
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

###Myeloid再分群----
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
  'Pro-angiogenic TAM' = "SPP1",       
  'Anti-inflammatory TAM' = "C1QC",   
  'Anti-inflammatory TAM' = "APOE",
  
  # 树突状细胞亚群
  "Dendritic" = "FCER1A",
  "Dendritic" = "JCHAIN",
  'cDC1' = "CLEC9A",                  
  'cDC1' = "XCR1",
  'cDC2' = "CD1C",
  'LAMP3+ DC' = "LAMP3",              
  
  # 粒细胞亚群
  'Neutrophil' = "CSF3R",
  'Neutrophil' = "CXCR2",
  'Eosinophil' = "SIGLEC8",
  'Basophil' = "CD123",
  
  # 肥大细胞
  'Mast' = "KIT",                 
  'Mast' = "TPSAB1",              
  "Mast" = "MS4A2",
  "Mast" = "CPA3",
  
  # 髓系祖细胞
  'CMP' = "CD34",                      
  'GMP' = "CD33"                       
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


#####cell props barplot  ----------------------------------------------------------
Macrophage_meta <- Macrophage[[]]

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

Macrophage[['subtype3']] <- Macrophage[['subtype2']]

###映射——add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(All_T[['subtype3']] %>% rownames_to_column('cells'), 
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

# 确认
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

#####细胞比例箱线图----
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
  
  message(sprintf("已生成 %s ", ct))
}


#####Imagedimplot----
Imagedimplot_dir <- 'result/Imagedimplot3'
dir.create(Imagedimplot_dir)



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
                    #shuffle.cols = F, 
                    #size = 1.2,
                    group.by = "subtype3"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #透明度
                    
  )
  ggsave(p,filename = paste0("result/Imagedimplot3/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
}
