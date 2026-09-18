library(tidyverse)
library(Seurat)
library(patchwork)

data_dir <- "/data/xenium/FWSP20250437"
data_dir <- "/data/xenium/FWSP20250437"
anno_file <- "data/Annotation.csv"

> anno_data
# A tibble: 3 × 2
  project                                        FOV                                         
  <chr>                                          <chr>                                         
1 12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo   12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo  
2 1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo 1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo
3 3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo  3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo 



# reading data -------------
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>%
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## add id for parse_data_path --------
parse_data_path <-
  parse_data_path %>%
  mutate(.id = seq_along(.proj))


## read into annotation -------
if (file.exists(anno_file)) {
  anno_data <- readr::read_csv(anno_file)
  User_columns <- colnames(anno_data)
} else {
  anno_data <- NULL
  User_columns <- "project"
}


## seurat object ------
all_data_ls <-
  parse_data_path %>%
  dplyr::select(.counts_path, .proj, .id) %>%
  pmap(function(.counts_path, .proj, .id) {
    # merge anno_data
    if (!is.null(anno_data)) {
      row_idx <- which(anno_data$project == .proj)
      col_idx <- anno_data %>% colnames()

      # 没有找到当前样本的注释数据
      #<<<<<< 所以Annotation.csv也有筛选文件的作用 >>>>>#
      if (length(row_idx) < 1) {
        message(str_glue("Not found sample {.proj} in column: project of Annotation table."))
        return(NULL)
      }
      # 没有任何数据，报错
      if (length(col_idx) < 1) stop("Not found any data in Annotation table.")

      meta_dat <- anno_data %>%
        .[row_idx, col_idx] %>%
        as.list()
    } else {
      meta_dat <- list(project = .proj)
    }

    # 构建Seurat
    sc <-
      Read10X(.counts_path) %>%
      .[["Gene Expression"]] %>%
      CreateSeuratObject(
        project      = .proj,
        min.cells    = 5,
        min.features = 5
      ) %>%
      AddMetaData(meta_dat) %>%
      RenameCells(add.cell.id = .id)
    sc
  })

all_data_ls <- all_data_ls %>% purrr::discard(is.null)
names(all_data_ls) <- map_chr(all_data_ls, "project.name")



### sample_df中的内容：sample为样本名，group为对应FOV名
sample_df <-
  data.frame(
    sample = c("kxz1", "kxz6", "kxz7", "kcz2", "kcz3", "kcz8", "gxz10", "gxz11", "gxz12", "gcz13", "gcz14", "gcz15"),
    FOV    = c("1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo", "1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo", "1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo", "1kongnanzuo-6kongnanzuo-7kongnanzuo-2kongnvzuo", "3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo", "3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo", "3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo", "3kongnvzuo-8kongnvzuo-10gaonanzuo-11gaonanzuo", "12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo", "12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo", "12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo", "12gaonanzuo-13gaonvzuo-14gaonvzuo-15gaonvzuo"),
    group  = c("AirMale", "AirMale", "AirMale", "AirFemale", "AirFemale", "AirFemale", "HyperoxiaMale", "HyperoxiaMale", "HyperoxiaMale", "HyperoxiaFemale", "HyperoxiaFemale", "HyperoxiaFemale")
  )


roi_ls <- list.files(
  "ROI",
  full.names = TRUE
)

names(roi_ls) <-
  roi_ls %>%
  basename() %>%
  str_remove("\\..*")


roi_coords_ls <-
  roi_ls %>%
  map(function(path) {
    # 读取文件内容，跳过前两行
    dat <- read.csv(path, skip = 2)
    # 将列名转换为小写
    colnames(dat) <- dat %>%
      colnames() %>%
      str_to_lower()
    # 添加一个新的列 'cell'，其值与 'cell.id' 列相同
    dat <- dat %>% mutate(cell = cell.id) %>%
      # 将 'cell.id' 列转换为行名
      column_to_rownames("cell.id")
    dat
  })

.id <- c("1_", "1_", "1_", "3_", "3_", "1_", "2_", "3_", "3_", "2_", "2_", "2_")
# .id <- str_c(names(roi_coords_ls), "_")
# .id
names(.id) <- names(roi_coords_ls)
roi_coords_ls <- roi_coords_ls %>%
  map2(.id, ~ {
    .x$cell <- str_c(.y, .x$cell)
    return(.x)
  })


xenium_roi_ls <-
  roi_coords_ls %>%
  imap(function(roi, name) { # browser()
    # get Visium Object
    xenium_name <- sample_df$FOV[sample_df$sample == name]
    group_name <- sample_df$group[sample_df$sample == name]
    xenium <- all_data_ls[[xenium_name]]

    # subset visium
    cells <- roi$cell
    xenium_roi <- subset(xenium, cells = cells)
    xenium_roi@project.name <- name # 更新 project.name
    xenium_roi@meta.data$sample <- name
    xenium_roi@meta.data$group <- group_name

    xenium_roi
  })


cell_counts <- sapply(roi_coords_ls, function(x) length(x$cell))


# add sample_id or sample_name -------------------
xenium_roi_ls <-
  xenium_roi_ls %>%
  imap(function(xenium, name) { # browser()
    xenium$project <- name
    Idents(xenium) <- "project"
    # RenameCells(xenium, add.cell.id = name)
    xenium
  })

# visium_roi override visium
all_data_ls <- xenium_roi_ls


# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name <- names(all_data_ls)

### all_data_ls_pass为all_data_ls过滤细胞后
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])


