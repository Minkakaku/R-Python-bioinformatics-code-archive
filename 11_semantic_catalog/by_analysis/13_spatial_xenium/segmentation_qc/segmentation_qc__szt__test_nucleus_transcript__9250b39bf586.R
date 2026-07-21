library(Seurat)
library(tidyverse)

setwd("~")

list.files("data/00037G2")

df_nuc_boundary <- arrow::read_parquet("data/00037G2/nucleus_boundaries.parquet")
df_transcripts  <- arrow::read_parquet("data/00037G2/transcripts.parquet")

df_nuc_boundary
# # A tibble: 6 × 4
# cell_id    vertex_x vertex_y label_id
# <chr>         <dbl>    <dbl>    <int>
# 1 aaaaeffo-1    1431.     378.        1
# 2 aaaaeffo-1    1431.     378.        1
# 3 aaaaeffo-1    1429.     378.        1
# 4 aaaaeffo-1    1429.     378.        1
# 5 aaaaeffo-1    1429.     378.        1
# 6 aaaaeffo-1    1428      379.        1

df_transcripts
# # A tibble: 6 × 11
# transcript_id cell_id    overlaps_nucleus feature_name x_location y_location z_location    qv fov_name nucleus_distance codeword_index
# <dbl> <chr>                 <int> <chr>             <dbl>      <dbl>      <dbl> <dbl> <chr>               <dbl>          <int>
# 1       2.82e14 UNASSIGNED                0 ITGAX             219.        738.       15.2  40   L2                   520.            195
# 2       2.82e14 UNASSIGNED                0 TAP2               50.2       726.       15.2  40   L2                   658.            353
# 3       2.82e14 UNASSIGNED                0 ZEB2               35.8       728.       15.4  39.4 L2                   669.            112
# 4       2.82e14 UNASSIGNED                0 AIF1               95.2       750.       15.1  35.8 L2                   608.             21
# 5       2.82e14 UNASSIGNED                0 TRAC              233.        731.       15.1  17.9 L2                   515.            341
# 6       2.82e14 UNASSIGNED                0 ITM2C              56.2       749.       15.3  40   L2                   640.            222


# add cell into ROI
# ROI2Cell <- 
#   ROI_ls %>%
#   imap(function(roi, proj){
#     cell_position <- cells_position_ls[[proj]]
#     if(is.null(cell_position)) stop(str_glue("Not found cell_position for tissue {proj}"))
#     
#     roi_ls <- roi %>% group_by(ID) %>% nest() %>% deframe()
#     roi_ls <-
#       lapply(roi_ls, function(x)
#         sp::point.in.polygon(cell_position$x_centroid, # point
#                              cell_position$y_centroid, # point
#                              x$X, # polygon
#                              x$Y)) %>% # polygon
#       map(as.logical) %>% # only 0 is FALSE, all others is TRUE
#       map(which) %>%
#       map(~ cell_position$cell_id[.x])
#     
#     roi_ls
#   })


Sys.time()
nuc_ls <- df_nuc_boundary %>% group_by(cell_id) %>% nest() %>% deframe() %>% head(10)
nuc_trans_ls <-
  lapply(nuc_ls, function(x)
    sp::point.in.polygon(df_transcripts$x_location, # point
                         df_transcripts$y_location, # point
                         x$X, # polygon
                         x$Y)) %>% # polygon
  map(as.logical) %>% # only 0 is FALSE, all others is TRUE
  map(which) %>%
  map(~ df_transcripts$transcript_id[.x])

Sys.time()





