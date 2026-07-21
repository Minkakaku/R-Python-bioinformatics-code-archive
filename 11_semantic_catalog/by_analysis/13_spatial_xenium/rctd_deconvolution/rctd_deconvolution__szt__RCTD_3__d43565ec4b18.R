
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)

# all_data <- readRDS('result_raw.rds')
# Idents(all_data) <- "cluster"
# 
# cell_cluster <- c(
#   '0'    = "Granulocyte",
#   '1'    = "Granulocyte",
#   '2'    = "Granulocyte",
#   '3'    = "Granulocyte",
#   '4'    = "B cells",
#   '5'    = "Monocyte macrophages",
#   '6'    = "Granulocyte",
#   '7'    = "Granulocyte",
#   '8'    = "B cells",
#   '9'    = "Erythrocytes",
#   '10'    = "Progenitor cells",
#   '11'    = "Chondrocyte",
#   '12'    = "T cells",
#   '13'    = "Osteoblast",
#   '14'    = "Granulocyte",
#   '15'    = "Granulocyte",
#   '16'    = "Megakaryocytes",
#   '17'    = "Erythrocytes",
#   '18'    = "Monocyte macrophages",
#   '19'    = "Endothelial cells",
#   '20'    = "Basophils",
#   '21'    = "Dendritic cells",
#   '22'    = "granulocyte",
#   '23'    = "Reticular cells",
#   '24'    = "B cells",
#   '25'    = "Smooth muscle cells"
# )


# all_data <- RenameIdents(all_data, cell_cluster)
# all_data$annotation <- Idents(all_data)


##spatial--------
# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


# read data -------
sample_df <- 
  # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('4M-7','24M-16'),
    group  = c('4M','24M')
  )

file_list <- 
  list.files(
    "/mnt/NexData_file/nas/outputdata/spaceranger_out/FWSP20250661_new",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("\\..*$", negate = TRUE)

stopifnot(nrow(file_list) == nrow(sample_df))

# TODO: fix bugs
# Samples order in file_list may not identical to sample_df
# names(file_list) <- sample_df$group # basename(file_list)
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]


visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, spatial.dir = image.dir, image.name = "tissue_hires_image.png", 
                assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(spatial.dir, 
                                                                   "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(spatial.dir, 
                                                                         "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
                         type = "centroids", radius = scale.factors[["spot"]], 
                         assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                          molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                          image = image, scale.factors = scale.factors)
        return(visium.fov)
      }
    
    Load10X_Spatial(
      path, 
      bin.size = 8, 
      slice = name, 
      image = Read10X_Image(
        file.path(path, "spatial"),
        spatial.dir = file.path(path, "binned_outputs/square_008um/spatial"),
        assay = 'Spatial.008um'
      )
    )
  })




# add sample_id or sample_name
visim_ls <-
  visim_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 


all_data <- merge(visim_ls[[1]], visim_ls[-1])


roi_df <-
  list.files('ROI', full.names = TRUE, pattern = "csv") %>%
  # setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
  setNames(., str_extract(., "(?<=ROI/)[^.]+")) %>% 
  map(~ read.csv(.x, col.names = c('Barcode', 'ROI'))) %>%
  enframe("group") %>%
  unnest(value) %>%
  select(Barcode,ROI)



roi_group <- matrix(
  c(
    '24M-1', '24M',
    '24M-2', '24M',
    '4M-1', '4M',
    '4M-2', '4M'
  ),
  ncol  = 2,
  byrow = TRUE
) %>% 
  as.data.frame() %>%
  setNames(c("ROI", "group"))

roi_df <- 
  roi_df %>%
  left_join(
    roi_group,
    by = "ROI"
  ) %>%
  mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
  column_to_rownames("Barcode")

all_data <- all_data %>% subset(cells = rownames(roi_df))

all_data <-
  all_data %>%
  AddMetaData(roi_df %>% select(ROI))

saveRDS(all_data, "RCTD3/all_data.rds")

# RCTD -----------------------
library(spacexr)
rctd_out_dir <- 'RCTD3' # dir.create(rctd_out_dir)

#2. single data
all_data_sc <- readRDS('FWSC20250984.25.6.9.rds')
all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$cell_type_new

# 查看 annotation 列的原始值
unique(all_data_sc$annotation)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation <- gsub("/", "_", all_data_sc$annotation)
all_data_sc$annotation <- gsub("fibroblast_mesenchymal stem cell", "synovial fibroblast", all_data_sc$annotation)
all_data_sc$annotation <- gsub("chondrocyte", "Chondrocyte", all_data_sc$annotation)
all_data_sc$annotation <- gsub("synovial intimal fibroblast", "Synovial intimal fibroblast", all_data_sc$annotation)
# 查看替换后的值
unique(all_data_sc$annotation)
saveRDS(all_data_sc, "all_data_sc_250613.rds")

Idents(all_data_sc) <- "annotation"

count_sc_data <- readRDS('adata_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")
counts_sc  <- count_sc_data %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$total_counts

sc_cell_filter_idx <-  cluster_sc != "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
visim_ls4 <-all_data %>% SplitObject(split.by = "ROI")
rctd_query_ls <-
  visim_ls4 %>%
  map(function(visium){# browser()
    counts_hd <- visium %>% LayerData(layer = "counts")
    coords    <- visium %>% GetTissueCoordinates() %>% .[, 1:2]
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(
      query, reference, max_cores = 10,
      UMI_min = 50
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })




saveRDS(rctd_result_ls, 'RCTD3/rctd_result_ls.rds')

#5. addMetaData
visium_rctd_ls <-
  visim_ls4 %>%
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Other" # "Unknown"
    visium
  })


all_data_rctd <- merge(visium_rctd_ls[[1]], visium_rctd_ls[-1])
#6. spatial plot 
# SpatialFeaturePlot(
#   visium_roi_ls[[1]],
#   image.scale = 'hires',
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

image_scale_ls <-
  visium_rctd_ls %>%
  # visium_roi_ls %>%
  map(GetTissueCoordinates) %>%
  map('x') %>%
  map(range) %>%
  map(diff) %>%
  map(log) %>%
  {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}

# RCTD Spatial plot
pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 3 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()


# RCTD Spatial plot use custom_color

# plan <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#FFEEAD",
#           "#D9534F","#FFAD60","#57D1C9","#ED5485","#7F95D1","#FFE869",
#           "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
#           "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
#           "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
#           "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
#           "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
#           "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
#           "#F17D80","#737495")

plan1 <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#ED5485",
          "#D9534F","#FFAD60","#57D1C9","#FFEEAD","#7F95D1","#FFE869",
          "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
          "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
          "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
          "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
          "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
          "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
          "#F17D80","#737495")

plan2 <- c("#5050ff","#ce3d32","#749b58","#cfe685","#466983","#ba6338",
           "#ef5285","#802268","#6bd766","#d595a7","#ffc300","#96ceb4",
           "#c75127","#d58f5c","#7a65a5","#005995")

pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD-1.pdf'), 
  height = 15, 
  width = 15
)

visium_rctd_ls %>%
  imap(function(visium, name){ 
    # if(name == "HC") name_fix <- "NML" else 
    name_fix <- name
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 3 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = plan2, na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
  })
dev.off()



## RCTD remove non-tissue area -------
pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD - only tissue area.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){ # browser()
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    visium$first_type[visium$first_type == "Other"] <- NA
    
    
    image  <- visim_image_he_ls[[name]]
    image  <- RenameCells(image, new.names = paste0(name, "_", Cells(image)))
    
    coords <- visim_coords_ls[[name]]
    coords_fixed <- paste(name, coords, sep = "_")
    
    image_sub  <- subset(image,  cells = coords_fixed)
    visium_sub <- subset(visium, cells = coords_fixed)
    visium_sub@images[[1]] <- image_sub
    
    
    # all_celltype <- visium_sub$first_type %>% unique
    # cell_colors  <- structure(
    #   names = all_celltype,
    #   ggsci::pal_igv()(length(all_celltype))
    # )
    # cell_colors <- c(cell_colors, Other = NA)
    
    # remove others 
    visium_sub <- visium_sub[, !is.na(visium_sub$first_type)]
    
    SpatialDimPlot(
      visium_sub,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = cell_colors) +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()
