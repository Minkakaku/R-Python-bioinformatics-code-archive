library(Seurat)
library(tidyverse)
library(CellTrek)


celltrek_out_dir <- 'CellTrek'
dir.create(celltrek_out_dir)


visim_ls <- readRDS("data_analysis_cycle2/visim_ls.rds")
# brain_st_cortex <- readRDS("brain_st_cortex.rds")

all_data_sc <- readRDS('data/hast_data.rds')

## Rename the cells/spots with syntactically valid names
visim_ls2 <- visim_ls %>%
  map(function(visium){
    visium <- RenameCells(visium, new.names=make.names(Cells(visium)))
    visium
  })
all_data_sc2 <- RenameCells(all_data_sc, new.names=make.names(Cells(all_data_sc)))

## Visualize the ST data
SpatialDimPlot(visim_ls2[[1]])

# We first co-embed ST and scRNA-seq datasets using traint
visim_traint_ls <- 
  visim_ls2 %>%
  map(function(visim){
    CellTrek::traint(
      st_data=visim, 
      sc_data=all_data_sc2,
      sc_assay='RNA',  
      cell_names='RNA_snn_res.0.05' # <-------------
    )
  })


## We can check the co-embedding result to see if there is overlap between these two data modalities
DimPlot(visim_traint_ls[[1]], group.by = "type") 


# After coembedding, we can chart single cells to their spatial locations.
# Here, we use the non-linear interpolation (intp = T, intp_lin=F) 
# approach to augment the ST spots.
visim_celltrek_ls <- 
  visim_ls2 %>%
  map(function(visium){
    celltrek <-
      CellTrek::celltrek(
        st_sc_int=brain_traint,
        int_assay='traint',
        sc_data=all_data_sc2,
        sc_assay = 'RNA', 
        reduction='pca',
        intp=T,
        intp_pnt=5000,
        intp_lin=F,
        nPCs=30,
        ntree=1000, 
        dist_thresh=0.55,
        top_spot=5,
        spot_n=5, # 8um: 2
        repel_r=20,
        repel_iter=20,
        keep_model=T
      )$celltrek
    
    celltrek$cell_type <- factor(celltrek$cell_type, levels=sort(unique(celltrek$cell_type)))
    celltrek
  })






brain_celltrek <- CellTrek::celltrek(
  st_sc_int=brain_traint,
  int_assay='traint',
  sc_data=all_data_sc2,
  sc_assay = 'RNA', 
  reduction='pca',
  intp=T,
  intp_pnt=5000,
  intp_lin=F,
  nPCs=30,
  ntree=1000, 
  dist_thresh=0.55,
  top_spot=5,
  spot_n=5, # 8um: 2
  repel_r=20,
  repel_iter=20,
  keep_model=T
)$celltrek


# After cell charting, we can interactively visualize the CellTrek result using celltrek_vis
brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, levels=sort(unique(brain_celltrek$cell_type)))

CellTrek::celltrek_vis(
  brain_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
  brain_celltrek@images$anterior1@image, brain_celltrek@images$anterior1@scale.factors$lowres
)





