library(Seurat)
library(tidyverse)
library(infercnv)

options(scipen = 100)


infercnv_in_dir  <- "infercnv/infercnv_files"
infercnv_out_dir <- "infercnv/infercnv_run"

all_infercnv_obj_ls <- readRDS(file.path('infercnv/', 'all_infercnv_obj_ls.rds'))


all_infercnv_obj_ls %>%
  iwalk(function(infer_obj, name){ 
    print(str_glue("{name} running..."))
    print(str_glue("gene count: "))
    print(all_infercnv_obj_ls[[name]]@count.data %>% rowMeans %>% { sum(. > 0.01) })
    
    tryCatch(
      {
        infer_obj = infercnv::run(
          infer_obj, 
          cutoff              = 0.01, #(see infercnv::run documentation)
          out_dir             = file.path(infercnv_out_dir, name), 
          cluster_by_groups   = FALSE, #unsupervised analysis
          HMM                 = FALSE, 
          write_expr_matrix   = TRUE,
          write_phylo         = TRUE,
          num_threads         = 20,
          no_prelim_plot      = TRUE,
          no_plot             = TRUE,
          # inspect_subclusters = FALSE,
          # useRaster           = FALSE,
          output_format     = "pdf",
          denoise             = TRUE
        ) #denoising applies noise reduction for the plot
        
        print(str_glue("{name} done."))
      },
      error = function(e){
        print(str_glue("{name} failed."))
      }
    )
  })

print("Done.")


# undebug(infercnv:::define_signif_tumor_subclusters)

# infercnv_obj <-
#   all_cell_id_used_ls['A246283'] %>%
#   imap(function(cell_id, name){
#     
#     dir.create(file.path(infercnv_in_dir, name))
#     
#     all_data <- all_data_used[, cell_id] %>% JoinLayers()
#     
#     #  export meta
#     all_data[['group_all']] %>% 
#       data.table::fwrite(
#         file = file.path(infercnv_in_dir, name, "meta_data_infer_subset.txt"),
#         row.names = TRUE,
#         col.names = FALSE,
#         sep = "\t"
#       )
#     
#     
#     raw_counts_matrix_infercnv <- LayerData(all_data, "counts")
#     
#     infercnv_obj = infercnv::CreateInfercnvObject(
#       raw_counts_matrix = raw_counts_matrix_infercnv,
#       annotations_file  = file.path(infercnv_in_dir, name, "meta_data_infer_subset.txt"),
#       gene_order_file   = file.path('infercnv', "GRCh38_gene_id_pos.txt"),
#       
#       min_max_counts_per_cell = c(100, Inf), # <--- qc prefilter ------
#       delim             = "\t",
#       ref_group_names   = c("NORMAL")
#     )
#     
#     infercnv_obj
#   })



infercnv::run(
  all_infercnv_obj_ls$A246283,
  cutoff            = 0.01, # (see infercnv::run documentation)
  out_dir           = file.path(infercnv_out_dir, 'A246283'),
  cluster_by_groups = FALSE, # unsupervised analysis
  HMM               = FALSE,
  write_expr_matrix = TRUE,
  write_phylo       = TRUE,
  num_threads       = 20,
  output_format     = "pdf",
  no_prelim_plot    = TRUE,
  # no_plot           = TRUE,
  denoise           = TRUE
)



