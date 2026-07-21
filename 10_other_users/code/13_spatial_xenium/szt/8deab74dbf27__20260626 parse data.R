library(Seurat)
library(tidyverse)

out_dir <- 'data_analysis_cycle2'
dir.create(out_dir)

#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}

##### nestSingleList
nestSingleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  name <- names(x)
  id_ <- seq_along(name)
  stopifnot(all(str_detect(name, sep)))
  name_1 <- str_split(name, pattern = sep) %>% map_chr(~ .x[1])
  name_2 <- str_split(name, pattern = sep) %>% map(~ .x[2:length(.x)]) %>% map_chr(str_c, collapse = sep)
  
  res_ls <- list()
  for(i in id_){    
    n1 <- name_1[i]
    n2 <- name_2[i]
    res <- x[[i]]
    res_ls[[n1]][[n2]] <- res
  }
  res_ls
}

# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


# set default -------------
bin.size     <- '008um' # '016um' 
assay.name   <- str_glue("Spatial.{bin.size}")
image.name   <- "tissue_hires_image.png" # "tissue_lowres_image.png"

stopifnot(bin.size %in% c('008um', '016um' ))

# read data -------
sample_df <- read.csv("data/Annoation.csv")
# data.frame(
#   sample = c('838375-1', '838375-2'),
#   group  = c('Tumor', 'Handover')
# )

file_list <- 
  list.files(
    '/mnt/NexData_file/nas/outputdata/spaceranger_out/FWSP20240200/', # 'data/FWSP20240200', # "/home/ztt/spatial/visium_HD/FWSP20240200",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) 
names(file_list) <- basename(file_list)


visim_ls <-
  file_list %>%
  imap(function(path, name){
    path <- file.path(path)
    spatial_path_total <- file.path(path, 'spatial')
    spatial_path_bin <- file.path(path, "binned_outputs", str_glue("square_{bin.size}"), 'spatial')
    
    
    Read10X_Image <- function (image.dir, data.dir, image.name = "tissue_lowres_image.png", 
        assay = "Spatial", slice = "slice1", filter.matrix = TRUE){
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(data.dir, 
            "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(data.dir, 
            "*tissue_positions*")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
            type = "centroids", radius = scale.factors[["spot"]], 
            assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
            molecules = fov@molecules, assay = fov@assay, key = fov@key, 
            image = image, scale.factors = scale.factors)
        return(visium.fov)
    }
    
    # read hi-res/low-res
    image <- Read10X_Image(
      image.dir  = spatial_path_total,
      data.dir   = spatial_path_bin,
      assay      = assay.name,
      slice      = name,
      image.name = image.name,
    )
    
    # visium
    Load10X_Spatial(
      path,
      bin.size = as.numeric(str_remove(bin.size,"um")), 
      image = image,
      slice = name
    )
  })


# parse ROI ---------------------
#1. read coordinate of  all cells  in all ROIs
roi_ls <- list.files(
  'data/ROI_coordinates_V3/',
  pattern    = "Projection", 
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("-.*$")

roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    dat <- read.csv(path)
    dat <- 
      dat %>%
      dplyr::rename(
        X = X.Coordinate,
        Y = Y.Coordinate
      )
    dat
  }) 

# compute coords
compute_coords <- FALSE

if(compute_coords){
  #2. compute ROI boundary from ROI coordinate
  roi_boundary_ls <-
    roi_coords_ls %>%
    map(function(x){
      idx <- chull(x[c("X", "Y")])
      dat <- x[idx, c("X", "Y")]
      rownames(dat) <- NULL
      dat
    })
  
  
  #3. get all coordinate of all ROIs
  # roi_ls have all coordinates, why select all points' coordinate manually?
  # 3.1 only 008um have coordinate, not 016um
  # 3.2 coordinates manualy computed were not all identical to roi_ls which download from Loupe
  
  roi_coords_recomputed_ls <-
    roi_boundary_ls %>%
    imap(function(bound_coords, name){ # browser()
      visium_name <- sample_df$sample[sample_df$group == name]
      visium <- visim_ls[[visium_name]]
      
      # all points
      coords <- GetTissueCoordinates(visium)
      
      # compute all_points
      # sp_idx: 
      # 0: point exterior to pol 
      # 1: point interior to pol
      # 2: point lies on a line
      # 3: point lies on a vertex
      sp_idx <-
        sp::point.in.polygon(
          point.x = coords$x,
          point.y = coords$y,
          pol.x = bound_coords$X, 
          pol.y = bound_coords$Y
        )
      coords[as.logical(sp_idx), ]
    })
}else{
  roi_coords_recomputed_ls <- 
    roi_coords_ls %>%
    map(function(x){
      colnames(x) <- x %>% colnames() %>% str_to_lower()
      x$cell <- x$barcode # rownames(x)
      x
    })
}


#4. subset Visium Seurat Object
visium_roi_ls <-
  roi_coords_recomputed_ls %>%
  imap(function(roi, name){ # browser()
    # get Visium Object
    visium_name <- sample_df$sample[sample_df$group == name]
    visium <- visim_ls[[visium_name]]
    
    # subset visium
    cells <- roi$cell
    visium_roi <- subset(visium, cells = cells)
    
    # subset image
    image_name <- make.names(str_glue("{visium_name}.{bin.size}"))
    image <- visium_roi@images[[image_name]]
    stopifnot(!is.null(image))
    Key(image) <- suppressWarnings(Key(name))
    
    # change image names
    visium_roi@images[[image_name]] <- NULL
    visium_roi@images[[name]] <- image
    
    visium_roi
  })

# # test spatial plot
# SpatialFeaturePlot(
#   visium_roi_ls[[1]], 
#   image.scale = 'hires', 
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

# normalize
visium_roi_ls <- 
  visium_roi_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 

# add sample_id or sample_name -------------------
visium_roi_ls <-
  visium_roi_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })

# visium_roi override visium ------------
visim_ls <- visium_roi_ls

# export
saveRDS(visim_ls, file.path(out_dir, 'visim_ls.rds') )



