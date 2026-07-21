library(tidyverse)
library(Seurat)


# Samples order in file_list may not identical to sample_df
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]
sample_match_table <- 
  data.frame(
    from = sample_df$group,
    to = names(file_list)
  )
# from  to
#    HC BP3
#   BP1 BP2
#   BP2 BP1
#   BP3  HC
stopifnot(! all(sample_match_table$from == sample_match_table$to))

# NML -> HC

# visim_ls --------------
visim_ls_images <-
  visim_ls %>%
  map('images') 

visim_ls_images_new <-
  visim_ls_images %>%
  map(function(img){# browser()
    img_name <- names(img)
    img_name_search <- str_remove(img_name, "\\.008um")
    img <- img[[1]]
    
    new_img_name_search <- sample_match_table$to[sample_match_table$from == img_name_search]
    new_img_name <- str_c(new_img_name_search, ".008um")
    
    new_img <- visim_ls_images[[new_img_name_search]][[new_img_name]]
    
    # return
    res <- list(new_img)
    names(res) <- new_img_name
    res
  })
names(visim_ls_images_new) <- sample_match_table$to



# visim_image_he_ls
# visium_rctd_ls

# lowres -> hires
all_data@images <-
  all_data@images %>%
  imap(function(img, img_name){ # browser()
    img_name_search <- str_remove(img_name, "\\..*$")
    new_img <- visim_image_he_ls[[img_name_search]]
    new_img <- new_img %>% RenameCells(., str_c(img_name_search, Cells(.), sep = "_"))
    new_img
  })

saveRDS(all_data, "all_data_projected_v2_hi_res.rds")

# all_data@images$HC.008um@image %>% grid::grid.raster()
# visim_image_he_ls$HC@image %>% grid::grid.raster()


# all_data.full.annotation@images <-
#   all_data.full.annotation@images %>%
#   imap(function(img, img_name){# browser()
#     
#     new_img_name <- sample_match_table$to[sample_match_table$from == img_name]
#     new_img <- all_data.full.annotation@images[[new_img_name]]
#     
#     new_img
#   }) %>% 
#   structure(names = sample_match_table$to)


all_data.full.annotation$project <-
  all_data.full.annotation$project %>%
  fct_recode(!!!structure(
    sample_match_table$from %>% str_replace("NML", "HC"),
    names = sample_match_table$to %>% str_replace("NML", "HC")
  ))

all_data.full.annotation.new.names <-
  all_data.full.annotation %>% 
  Cells() %>%
  str_remove("^.*_s") %>%
  str_c(all_data.full.annotation$project, "_s", ., sep = "")

all_data.full.annotation <-
  all_data.full.annotation  %>%
  RenameCells(new.names = all_data.full.annotation.new.names)




