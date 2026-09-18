# RCTD -----------------------
library(spacexr)
rctd_out_dir <- 'RCTD' # dir.create(rctd_out_dir)

# all_data <- readRDS("all_data.full.annotation.rds")

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

#2. single data
# all_data_sc <- readRDS('BP_AD_Pso-v2.rds')
all_data_sc <- readRDS('GSE207848.rds')
all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$clust_annot 

#sc_gene trans
library(biomaRt)
library(BiocFileCache)

mart <- useMart("ensembl","mmusculus_gene_ensembl")

##人类选择hsapiens_gene_ensembl
gene_name<-getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters = "ensembl_gene_id",values = as.data.frame(row.names(all_data_sc)), mart = mart)

gene_name_unique <- gene_name %>%  
  distinct(external_gene_name, .keep_all = TRUE) 

all_data_sub <- subset(all_data_sc,features=gene_name_unique$ensembl_gene_id)

index <- match(gene_name_unique$ensembl_gene_id, rownames(all_data_sub))
gene_name_sorted <- gene_name_unique[order(index), ]


# all_data_sub=RenameGenesSeurat(obj = all_data_sub, 
#                       newnames = gene_name_sorted$external_gene_name)

all_data_sub$annotation <- all_data_sub$clust_annot 
Idents(all_data_sub) <- "annotation"
colSums(counts_sc)%>%as.data.frame()%>%AddMetaData(all_data_sub,.)->all_data_sub
counts_sc  <- all_data_sub %>% JoinLayers() %>% LayerData(layer = "counts") 
rownames(counts_sc) <- gene_name_sorted$external_gene_name
cluster_sc <- as.factor(all_data_sub$annotation)
nUMI_sc    <- all_data_sub$.


sc_cell_filter_idx <-  cluster_sc != "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
rctd_query_ls <-
  visim_ls %>%
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

# saveRDS(rctd_result_ls, 'rctd_result_ls.rds')

#5. addMetaData
visium_rctd_ls <-
  visim_ls %>%
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Other" # "Unknown"
    visium
  })

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

pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){
    # if(name == "HC") name_fix <- "NML" else name_fix <- name
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()
ll<-all_data[,!is.na(all_data$first_type)]
p_sp_plot_anno <-
  Images(ll)[5:6] %>%
  setNames(., c('DZZ_Saline','SYZ_NH4')) %>%
  imap(function(img, img_name){
    ll@images[[img_name]] <- 
      ll@images[[img_name]] %>%
      subset(cells = colnames(ll))
    p <- SpatialDimPlot(
      ll,
      images = img, 
      image.scale = "hires", 
      group.by = "first_type",
      pt.size.factor	= 3,
      crop  = FALSE
    ) +
      #ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values=col_vector)+
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'),
  height = 10, width = 13
)
p_sp_plot_anno %>% walk(print)
dev.off()

saveRDS(visium_rctd_ls,'RCTD.rds')
#映射-----------------
library(tibble)
library(dplyr)
meta.data <-
  rbind(visium_rctd_ls$DZZ_Saline_34_0um[[c("spot_class","first_type","second_type","first_class","second_class","min_score","singlet_score","conv_all","conv_doublet")]],
        visium_rctd_ls$DZZ_Saline_34_100um[[c("spot_class","first_type","second_type","first_class","second_class","min_score","singlet_score","conv_all","conv_doublet")]],
        visium_rctd_ls$SYZ_NH4_44_0um[[c("spot_class","first_type","second_type","first_class","second_class","min_score","singlet_score","conv_all","conv_doublet")]],
        visium_rctd_ls$SYZ_NH4_44_100um[[c("spot_class","first_type","second_type","first_class","second_class","min_score","singlet_score","conv_all","conv_doublet")]])

all_data<-AddMetaData(all_data,meta.data)  

pdf(file.path(rctd_out_dir, str_glue("Dimplot rctd annotation.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "first_type", label = TRUE, reduction = 'harmonyumap',repel = TRUE)+ NoLegend()
dev.off()


##celltype spatial-------------
Idents(all_data) <- "first_type"
cells <- CellsByIdentities(all_data)
Images(all_data)[5:6] %>%
  setNames(.,c("DZZ_Saline","SYZ_NH4")) %>%
  iwalk(function(img, img_name){
    pdf(file.path(rctd_out_dir, str_glue("celltype spatial - {img_name}.pdf")),width =20, height = 20)
    print(SpatialDimPlot(all_data,
                         image.scale = "hires",
                         images = img, 
                         pt.size.factor	= 1,
                         cells.highlight = cells[setdiff(names(cells), "NA")],
                         cols.highlight = c("#FFFF00", "grey50"),
                         facet.highlight = T, combine = T) + NoLegend())
    dev.off()
  })
Images(all_data)[5:6] %>%
  setNames(.,c("DZZ_Saline","SYZ_NH4")) %>%
  iwalk(function(img, img_name){
    tiff(file.path(rctd_out_dir, str_glue("celltype spatial - {img_name}.tiff")),width =600, height = 600)
    print(SpatialDimPlot(all_data,
                         image.scale = "hires",
                         images = img, 
                         pt.size.factor	= 1,
                         cells.highlight = cells[setdiff(names(cells), "NA")],
                         cols.highlight = c("#FFFF00", "grey50"),
                         facet.highlight = T, combine = T) + NoLegend())
    dev.off()
  })



# 遍历每个注释类别并生成 DimPlot  
plots<-map(names(table(all_data$first_type)), function(name) {  
  p <- DimPlot(all_data, group.by = "annotation",  
               cells.highlight = WhichCells(all_data, idents = name),  
               cols.highlight = c("#DC0000B2"), cols = "grey", label = FALSE, reduction = 'harmonyumap') +  
    NoLegend() +  
    labs(title = name)
  return(p)
})  
combined_plot <- patchwork::wrap_plots(plots[1:44])
ggsave(file.path(rctd_out_dir, "celltype_split_umap.pdf"), combined_plot, width = 70, height = 70,limitsize = FALSE)

#prop barplot-------------
coords_enrichment <-
  all_data %>%
  {left_join(
    .[[c('first_type', "project")]] %>% rownames_to_column("cell"),
    GetTissueCoordinates(.),
    by = 'cell'
  )}
props_df <-
  coords_enrichment %>%
  group_by(project, first_type) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq = count/sum(count) * 100) 
p_prop_celltype <-
  props_df %>%
  
  ggplot(aes(x = project, y = freq, fill = first_type)) +
  geom_col() +
  ggpubr::theme_pubclean() +
  ggsci::scale_fill_igv() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Cell Proportion(%)",
    fill = "Cell Type"
  )
# export
ggsave(
  p_prop_celltype,
  filename = file.path(rctd_out_dir,"Celltype Props distribution.pdf"),
  height = 6,
  width  = 7
)


gene_to_check<-c('Slc17a6','Slc32a1','Slc6a5',
                 'Adcyap1','Phox2b','Dbh',
                 'Cck','Th','Tac1',
                 'Pdyn','Pvalb','P2ry1',
                 'Crhr2','Npy2r','Gpr65',
                 'Ascl1','Nmbr','Calca',
                 'Trpv1','Calb1','Calb2',
                 'Tac2','Htra1','Sst',
                 'Vip','Cbln1','Chat','Oprm1')
tiff(
  'RCTD/Spatial gene feature plot_DZZ2.tiff',
  width = 600, height =1000
)
ll@images[5] %>%
  imap(function(images, name){
    SpatialFeaturePlot(
      ll,
      images = name,
      feature = gene_to_check,
      shape = 21,
      alpha = c(0.3, 1),
      pt.size.factor = 50,
      ncol=4
    ) +
      patchwork::plot_annotation(title = name) +
      scale_fill_gradientn(transform='log1p',colors = Seurat:::SpatialColors(100))
  })
dev.off()

pdf(
  file.path(rctd_out_dir,"gene featureplot umap.pdf"), 
  width = 70, 
  height = 40)
FeaturePlot(
  all_data,
  features = gene_to_check,
  reduction = "harmonyumap",
  ncol = 7,
  raster = TRUE
)  
dev.off()


theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(rctd_out_dir, str_glue("dotplot check_gene.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = gene_to_check, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "first_type", features = gene_to_check, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "annotation", features = gene_to_check, dot.scale = 10) + theme_90angle
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