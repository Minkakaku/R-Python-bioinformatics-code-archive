library(Seurat)
library(tidyverse)

# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


# read data -------
sample_df <- read.csv("data/Annoation.csv")
# data.frame(
#   sample = c('838375-1', '838375-2'),
#   group  = c('Tumor', 'Handover')
# )

file_list <- 
  list.files(
    "/home/ztt/spatial/visium_HD/FWSP20240131",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("_[no|old|remove]", negate = TRUE)

names(file_list) <- sample_df$group # basename(file_list)


visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    Load10X_Spatial(path, bin.size = 16, slice = name)
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
  map(~ {DefaultAssay(.x) <- 'Spatial.016um'; .x}) %>%
  map(FindVariableFeatures, nfeatures = 3000)  %>%
  map(NormalizeData)

# LSGI --------
library(RcppML)  # https://github.com/zdebruine/RcppML
library(LSGI)
library(Matrix)

lsgi_out_dir <- "LSGI"
dir.create(lsgi_out_dir)

scan.nmf.mse <- function(obj, ranks = seq(1, 30, 2), tol = 1e-04) {
  # users can customize the scan by changing 'ranks'
  # dat <- obj@assays$RNA@data
  dat <- LayerData(obj, layer = "data")
  errors <- c()
  ranks <- seq(1, 30, 2)
  for (i in ranks) {
    # cat('rank: ', i, '\n')
    mod <- RcppML::nmf(dat, i, tol = 1e-04, verbose = F)
    mse_i <- mse(dat, mod$w, mod$d, mod$h)
    errors <- c(errors, mse_i)
  }
  results <- data.frame(rank = ranks, MSE = errors)
  return(results)
}

sr.nmf <- function(obj, k = 10, tol = 1e-06, assay = "RNA") {
  # dat <- obj@assays$RNA@data
  dat <- LayerData(obj, layer = "data")
  
  nmf_model <- RcppML::nmf(dat, k = k, tol = tol, verbose = F)
  embeddings <- t(nmf_model$h)
  rownames(embeddings) <- colnames(obj)
  colnames(embeddings) <- paste0("nmf_", 1:k)
  loadings <- nmf_model$w
  rownames(loadings) <- rownames(obj)
  obj@reductions$nmf <- CreateDimReducObject(embeddings = embeddings,
                                             loadings = loadings, key = "nmf_", assay = assay)
  return(obj)
}

# p_scan_nmf_ls <- visium_roi_ls %>%
#   imap(function(visium, name){
#     scan.nmf.res <- scan.nmf.mse(obj = visium)
#     ggplot(scan.nmf.res, aes(x = rank, y = MSE)) + geom_point(size = 0.7) +
#       geom_smooth(method = "loess", span = 0.2, color = "black",
#                   linewidth = 1, se = F) + labs(x = "NMF rank", y = "MSE") +
#       theme_classic() + scale_y_continuous(expand = c(0.01, 0)) +
#       theme(aspect.ratio = 1) +
#       ggtitle(name)
#   })

# run nmf
visim_ls_lsgi_ls <-
  visim_ls %>%
  map(function(x){
    subset(x, features = VariableFeatures(x))
  }) %>%
  map(function(visium){
    rank <- c(3, 5, 10)
    names(rank) <- paste0('rank', rank)
    
    rank %>%
      map(function(rank){
        visium <- sr.nmf(obj = visium, k = rank, tol = 1e-05)
        visium
      })
  })

saveRDS(visim_ls_lsgi_ls, file.path(lsgi_out_dir, 'visim_ls_lsgi_ls.rds'))

# run lsgi
lsgi_res_ls <-
  visim_ls_lsgi_ls %>%
  imap(function(visium_ls, name1){
    print(name1)
    
    visium_ls %>%
      imap(function(visium, name2){
        print(name2)
        
        spatial_coords <- GetTissueCoordinates(visium) %>% .[1:2]
        colnames(spatial_coords) <- c("X", "Y")
        
        embeddings <- visium@reductions$nmf@cell.embeddings
        
        divided_number_ls <- DescTools::Divisors(nrow(spatial_coords)) %>% .[[1]]
        divided_number    <- divided_number_ls[which.min(abs(divided_number_ls - 8))]
        
        divided_number <- 2 # 客户指定为2
         
        print(divided_number_ls)
        # return()
        
        lsgi.res <-
          local.traj.preprocessing(
            spatial_coords    = spatial_coords,
            n.grids.scale     = divided_number ,
            embeddings        = embeddings, 
            n.cells.per.meta  = 15
          )
        
        # export
        saveRDS(lsgi.res, file.path(lsgi_out_dir, str_glue('lsgi_res_{name1}_{name2}.rds')))
        
        lsgi.res
      })
  })

