
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)

rctd_outdir_5 <- "RCTD_20250814"
dir.create(rctd_outdir_5)


## spatial--------
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
  list.files('ROI_2', full.names = TRUE, pattern = "csv") %>%
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

saveRDS(all_data, "RCTD_20250814/all_data.rds")


# RCTD -----------------------
library(spacexr)


# 2. single data
all_data_sc <- readRDS('FWSC20250984.25.6.9.rds')
all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$cell_type_new

# 查看 annotation 列的原始值
unique(all_data_sc$annotation)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation <- gsub("/", "_", all_data_sc$annotation)
all_data_sc$annotation <- gsub("fibroblast_mesenchymal stem cell", "Synovial fibroblasts", all_data_sc$annotation)
all_data_sc$annotation <- gsub("chondrocyte", "Chondrocytes", all_data_sc$annotation)
all_data_sc$annotation <- gsub("synovial intimal fibroblast", "Synovial intimal fibroblasts", all_data_sc$annotation)
all_data_sc$annotation <- gsub("Granulocyte", "Granulocytes", all_data_sc$annotation)
all_data_sc$annotation <- gsub("Osteoblast", "Osteoblasts", all_data_sc$annotation)
all_data_sc$annotation <- gsub("Osteoclast", "Osteoclasts", all_data_sc$annotation)


# 查看替换后的值
unique(all_data_sc$annotation)

out_gene_umi=as.matrix(all_data_sc@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20250814/all_data_sc_meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)



# subcluster数据处理
library(tidyverse)
library(reticulate)

use_python("/home/szt/miniforge3/envs/rpy/bin/python")
reticulate::py_config()

sc <- import("scanpy")
adata <- sc$read("RCTD_20250814/result.h5ad")
counts_x <- adata$X %>% t %>% as.matrix()
colnames(counts_x) <- rownames(adata$obs)

scRNA2 <- CreateSeuratObject(
  counts_x,
  meta.data = adata$obs
)

out_gene_umi=as.matrix(scRNA2@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20250814/scRNA2_meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)


# 添加新的注释信息
annotation_new <- read.table(file = "RCTD_20250814/new_annotation.txt", header = T, sep = "\t")
annot_vec <- setNames(annotation_new$annotation_new, annotation_new$Cell)
all_data_sc_barcodes <- rownames(all_data_sc@meta.data)
new_annotations <- annot_vec[all_data_sc_barcodes]
all_data_sc@meta.data$annotation_new <- new_annotations



# 查看 annotation 列的原始值
unique(all_data_sc$annotation_new)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation_new <- gsub("/", "_", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("fibroblast_mesenchymal stem cell", "Synovial fibroblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("synovial intimal fibroblast", "Synovial intimal fibroblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Granulocyte-sub4", "Granulocytes-sub4", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Osteoblast-sub1", "Osteoblasts-sub1", all_data_sc$annotation_new)


# 查看替换后的值
unique(all_data_sc$annotation_new)

saveRDS(all_data_sc, "RCTD_20250814/all_data_sc_250825.rds")



# 单细胞数据处理
Idents(all_data_sc) <- "annotation_new"

count_sc_data <- readRDS('adata_sc_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")
counts_sc  <- count_sc_data %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation_new)
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




saveRDS(rctd_result_ls, 'RCTD_20250814/rctd_result_ls.rds')

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

saveRDS(all_data_rctd, file = "RCTD_20250814/all_data_rctd.rds")
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
  file.path(rctd_outdir_5, 'spatial plot of RCTD.pdf'), 
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
      ggsci::scale_fill_igv(na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
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
  file.path(rctd_outdir_5, 'spatial plot of RCTD-1.pdf'), 
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
      scale_fill_manual(values = plan1, na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
  })
dev.off()



### RCTD remove non-tissue area -------
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










### 细胞比例 -----
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

dir.create("RCTD_20250814/Ratio", recursive = T)

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")


out_gene_umi=as.matrix(all_data_rctd@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20250814/Ratio/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.table("RCTD_20250814/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,6,8)]
colnames(data) <- c("Cells","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data_filt <- data %>% filter(Annotation != "Other")
data <- data_filt

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250814/Ratio/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="RCTD_20250814/Ratio/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

samplemap <- read.table("RCTD_20250814/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
colnames(samplemap)<-c("Sample","Group")
final_data <- data %>% left_join(samplemap, by = "Sample")

# write.table(final_data,file="RCTD_20250814/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250814/Ratio/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file="RCTD_20250814/Ratio/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD_20250814/Ratio/ratio_barplot.pdf",width=15, height=6)
ggsave(p,file="RCTD_20250814/Ratio/ratio_barplot.png",width=15, height=6,dpi=300)









# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = "RCTD_20250814/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"

data_filt <- data %>% filter(Gene != "Other")
data <- data_filt

data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD_20250814/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD_20250814/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by="SampleID")

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="RCTD_20250814/Ratio/diff_boxplot.pdf",width=20,height=5)
ggsave(p,file="RCTD_20250814/Ratio/diff_boxplot.png",width=20,height=5)









# remove Erythrocyte --------------
rctd_outdir_5_remove_Erythrocyte <- "RCTD_20250814/remove_Erythrocyte"
dir.create("RCTD_20250814/remove_Erythrocyte")


all_data <- readRDS("RCTD_20250814/all_data.rds")

library(spacexr)

# 2. single data
all_data_sc <- readRDS("RCTD_20250814/all_data_sc_250825.rds")

# 单细胞数据处理
Idents(all_data_sc) <- "annotation_new"
# all_data_sc <- subset(all_data_sc, subset = !(annotation_new %in% c("Erythrocytes-9", "Erythrocytes-17")))

count_sc_data <- readRDS('adata_sc_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")
counts_sc  <- count_sc_data %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation_new)
nUMI_sc    <- all_data_sc$total_counts

sc_cell_filter_idx <-  cluster_sc != "Unknown"
sc_cell_filter_idx_new <-  !(cluster_sc %in% c("Unknown","Erythrocytes-9","Erythrocytes-17"))

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx_new], 
  droplevels(cluster_sc[sc_cell_filter_idx_new]), 
  nUMI_sc[sc_cell_filter_idx_new]
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


saveRDS(rctd_result_ls, 'RCTD_20250814/remove_Erythrocyte/rctd_result_ls.rds')


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

saveRDS(all_data_rctd, file = "RCTD_20250814/remove_Erythrocyte/all_data_rctd.rds")


#6. spatial plot 
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
  file.path(rctd_outdir_5_remove_Erythrocyte, 'spatial plot of RCTD.pdf'), 
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
      ggsci::scale_fill_igv(na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
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
  file.path(rctd_outdir_5_remove_Erythrocyte, 'spatial plot of RCTD-1.pdf'), 
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
      scale_fill_manual(values = plan1, na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
  })
dev.off()






### 细胞比例 --------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

dir.create("RCTD_20250814/remove_Erythrocyte/Ratio", recursive = T)

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")


out_gene_umi=as.matrix(all_data_rctd@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20250814/remove_Erythrocyte/Ratio/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.table("RCTD_20250814/remove_Erythrocyte/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,6,8)]
colnames(data) <- c("Cells","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data_filt <- data %>% filter(Annotation != "Other")
data <- data_filt

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250814/remove_Erythrocyte/Ratio/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="RCTD_20250814/remove_Erythrocyte/Ratio/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

samplemap <- read.table("RCTD_20250814/remove_Erythrocyte/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
colnames(samplemap)<-c("Sample","Group")
final_data <- data %>% left_join(samplemap, by = "Sample")

# write.table(final_data,file="RCTD_20250814/remove_Erythrocyte/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250814/remove_Erythrocyte/Ratio/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file="RCTD_20250814/remove_Erythrocyte/Ratio/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD_20250814/remove_Erythrocyte/Ratio/ratio_barplot.pdf",width=15, height=6)
ggsave(p,file="RCTD_20250814/remove_Erythrocyte/Ratio/ratio_barplot.png",width=15, height=6,dpi=300)









# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = "RCTD_20250814/remove_Erythrocyte/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"

data_filt <- data %>% filter(Gene != "Other")
data <- data_filt

data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD_20250814/remove_Erythrocyte/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD_20250814/remove_Erythrocyte/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by="SampleID")

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="RCTD_20250814/remove_Erythrocyte/Ratio/diff_boxplot.pdf",width=20,height=5)
ggsave(p,file="RCTD_20250814/remove_Erythrocyte/Ratio/diff_boxplot.png",width=20,height=5)









