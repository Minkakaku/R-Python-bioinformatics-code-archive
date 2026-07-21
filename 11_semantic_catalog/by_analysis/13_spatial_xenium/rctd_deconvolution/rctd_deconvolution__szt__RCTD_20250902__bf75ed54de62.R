
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)

rctd_outdir_6 <- "RCTD_20250902"
dir.create(rctd_outdir_6)


## spatial--------
# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 10Gb


# read data -------
# sample_df <- 
#   # read.csv("data/Annoation.csv")
#   data.frame(
#     sample = c('4M-7','24M-16'),
#     group  = c('4M','24M')
#   )
# 
# file_list <- 
#   list.files(
#     "/mnt/NexData_file/nas/outputdata/spaceranger_out/FWSP20250661_new",
#     pattern = paste(sample_df$sample, collapse = "|"), 
#     full.names = TRUE
#   ) %>%
#   str_subset("\\..*$", negate = TRUE)
# 
# stopifnot(nrow(file_list) == nrow(sample_df))
# 
# # TODO: fix bugs
# # Samples order in file_list may not identical to sample_df
# # names(file_list) <- sample_df$group # basename(file_list)
# names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]
# 
# 
# visim_ls <-
#   file_list %>% 
#   imap(function(path, name){
#     path <- file.path(path, "outs")
#     
#     Read10X_Image <-
#       function (image.dir, spatial.dir = image.dir, image.name = "tissue_hires_image.png", 
#                 assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
#         
#         image <- png::readPNG(source = file.path(image.dir, image.name))
#         scale.factors <- Read10X_ScaleFactors(filename = file.path(spatial.dir, 
#                                                                    "scalefactors_json.json"))
#         coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(spatial.dir, 
#                                                                          "*tissue_positions*parquet")), filter.matrix)
#         fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
#                          type = "centroids", radius = scale.factors[["spot"]], 
#                          assay = assay, key = Key(slice, quiet = TRUE))
#         visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
#                           molecules = fov@molecules, assay = fov@assay, key = fov@key, 
#                           image = image, scale.factors = scale.factors)
#         return(visium.fov)
#       }
#     
#     Load10X_Spatial(
#       path, 
#       bin.size = 8, 
#       slice = name, 
#       image = Read10X_Image(
#         file.path(path, "spatial"),
#         spatial.dir = file.path(path, "binned_outputs/square_008um/spatial"),
#         assay = 'Spatial.008um'
#       )
#     )
#   })
# 
# 
# # add sample_id or sample_name
# visim_ls <-
#   visim_ls %>%
#   imap(function(visium, name){
#     visium$project <- name
#     Idents(visium) <- 'project'
#     RenameCells(visium, add.cell.id = name)
#   })
# 
# 
# visim_ls <- 
#   visim_ls %>%
#   map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
#   map(NormalizeData) 
# 
# 
# all_data <- merge(visim_ls[[1]], visim_ls[-1])
# 
# 
# roi_df <-
#   list.files('ROI_2', full.names = TRUE, pattern = "csv") %>%
#   # setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
#   setNames(., str_extract(., "(?<=ROI/)[^.]+")) %>% 
#   map(~ read.csv(.x, col.names = c('Barcode', 'ROI'))) %>%
#   enframe("group") %>%
#   unnest(value) %>%
#   select(Barcode,ROI)
# 
# 
# 
# roi_group <- matrix(
#   c(
#     '24M-1', '24M',
#     '24M-2', '24M',
#     '4M-1', '4M',
#     '4M-2', '4M'
#   ),
#   ncol  = 2,
#   byrow = TRUE
# ) %>% 
#   as.data.frame() %>%
#   setNames(c("ROI", "group"))
# 
# roi_df <- 
#   roi_df %>%
#   left_join(
#     roi_group,
#     by = "ROI"
#   ) %>%
#   mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
#   column_to_rownames("Barcode")
# 
# all_data <- all_data %>% subset(cells = rownames(roi_df))
# 
# all_data <-
#   all_data %>%
#   AddMetaData(roi_df %>% select(ROI))
# 
# saveRDS(all_data, "RCTD_20250902/all_data.rds")

all_data <- readRDS("data_for_analysis/all_data_projected.rds")

# RCTD -----------------------
library(spacexr)


# 2. single data --------
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
write.table(out_gene_umi,file="RCTD_20250902/all_data_sc_meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)



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
write.table(out_gene_umi,file="RCTD_20250902/scRNA2_meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)


# 添加新的注释信息
annotation_new <- read.table(file = "RCTD_20250902/new_annotation.txt", header = T, sep = "\t")
annot_vec <- setNames(annotation_new$annotation_new, annotation_new$Cell)
all_data_sc_barcodes <- rownames(all_data_sc@meta.data)
new_annotations <- annot_vec[all_data_sc_barcodes]
all_data_sc@meta.data$annotation_new <- new_annotations



# 查看 annotation 列的原始值
unique(all_data_sc$annotation_new)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation_new <- gsub("/", "_", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("fibroblast_mesenchymal stem cell", "Synovial fibroblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub0", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub3", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub5", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub7", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub8", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub9", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Osteoblast-sub1", "Osteoblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Osteoclasts-18", "Osteoclasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("synovial intimal fibroblast", "Synovial intimal fibroblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Granulocyte-sub4", "Granulocytes-sub4", all_data_sc$annotation_new)


# 查看替换后的值
unique(all_data_sc$annotation_new)

saveRDS(all_data_sc, "RCTD_20250902/all_data_sc_250902.rds")



# remove other --------------

# 单细胞数据处理
Idents(all_data_sc) <- "annotation_new"

count_sc_data <- readRDS('adata_sc_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")

count_sc_data <- subset(count_sc_data, cells = rownames(all_data_sc@meta.data))

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
visim_ls4 <-all_data %>% SplitObject(split.by = "sample")
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



rctd_outdir_6_remove_other <- "RCTD_20250902/remove_other"
dir.create("rctd_outdir_6_remove_other")
saveRDS(rctd_result_ls, 'RCTD_20250902/remove_other/rctd_result_ls.rds')

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

saveRDS(all_data_rctd, file = "RCTD_20250902/remove_other/all_data_rctd.rds")
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
  file.path(rctd_outdir_6_remove_other, 'spatial plot of RCTD.pdf'), 
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
  file.path(rctd_outdir_6_remove_other, 'spatial plot of RCTD-1.pdf'), 
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






### 细胞比例 -----
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

dir.create("RCTD_20250902/remove_other/Ratio", recursive = T)


out_gene_umi=as.matrix(all_data_rctd@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20250902/remove_other/Ratio/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.table("RCTD_20250902/remove_other/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,11,27)]
colnames(data) <- c("Cells","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data_filt <- data %>% filter(Annotation != "Other")
data <- data_filt

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250902/remove_other/Ratio/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

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
write.table(DATA,file="RCTD_20250902/remove_other/Ratio/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

samplemap <- read.table("RCTD_20250902/remove_other/Ratio/sample_map_group.txt", header=F, sep="\t", check.names=F)
colnames(samplemap)<-c("Sample","Group")
final_data <- data %>% left_join(samplemap, by = "Sample")

# write.table(final_data,file="RCTD_20250902/remove_other/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250902/remove_other/Ratio/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

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
write.table(DATA,file="RCTD_20250902/remove_other/Ratio/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
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
  plot_layout(widths = c(1.5, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD_20250902/remove_other/Ratio/ratio_barplot.pdf",width=18, height=6)
ggsave(p,file="RCTD_20250902/remove_other/Ratio/ratio_barplot.png",width=18, height=6,dpi=300)



# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = "RCTD_20250902/remove_other/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"

# data_filt <- data %>% filter(Gene != "Other")
# data <- data_filt

data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD_20250902/remove_other/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD_20250902/remove_other/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
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
ggsave(p,file="RCTD_20250902/remove_other/Ratio/diff_boxplot.pdf",width=20,height=5)
ggsave(p,file="RCTD_20250902/remove_other/Ratio/diff_boxplot.png",width=20,height=5)








# 对每种细胞类型进行正态性检验
shapiro_results <- plot_data %>%
  group_by(Gene) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(value)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(value)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "RCTD_20250902/remove_other/Ratio/shapiro_test_results.csv", row.names = FALSE)






# > plot_data$Gene %>% unique()
# [1] B cells-24                         B cells-4                          B cells-8                          Basophils-20                      
# [5] Chondrocytes                       Dendritic cells-21                 Endothelial cells-19               Erythrocytes-17                   
# [9] Erythrocytes-9                     Granulocytes-0                     Granulocytes-1                     Granulocytes-14                   
# [13] Granulocytes-15                    Granulocytes-2                     Granulocytes-3                     Granulocytes-6                    
# [17] Granulocytes-7                     Granulocytes-sub4                  Megakaryocytes-16                  Monocyte macrophages-5            
# [21] Osteoblasts                        Osteoclasts                        Reticular cells-23                 Smooth muscle cells-25            
# [25] Synovial fibroblasts-sub2          Synovial intimal fibroblasts-sub10 T cells-12                        






plot_data_Chondrocytes  <- subset(plot_data, Gene == "Chondrocytes")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Chondrocytes, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 11),
    breaks = seq(0, 11, by = 2)
  ) +
  # 标签和标题
  labs(
    title = "Chondrocytes",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Chondrocytes_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Chondrocytes_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Osteoblasts <- subset(plot_data, Gene == "Osteoblasts")
p_detailed <- ggplot(plot_data_Osteoblasts, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(6, 9),
    breaks = seq(6, 9, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Osteoblasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Osteoblasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Osteoblasts_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Osteoclasts <- subset(plot_data, Gene == "Osteoclasts")
p_detailed <- ggplot(plot_data_Osteoclasts, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1.5),
    breaks = seq(0, 1.5, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Osteoclasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Osteoclasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Osteoclasts_plot.pdf", p_detailed, width = 10, height = 8)










plot_data_Synovial_fibroblasts_sub2 <- subset(plot_data, Gene == "Synovial fibroblasts-sub2")
p_detailed <- ggplot(plot_data_Synovial_fibroblasts_sub2, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.15, 0.4),
    breaks = seq(0.15, 0.4, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Synovial fibroblasts-sub2",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Synovial_fibroblasts_sub2_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Synovial_fibroblasts_sub2_plot.pdf", p_detailed, width = 10, height = 8)











plot_data_Synovial_intimal_fibroblasts_sub10 <- subset(plot_data, Gene == "Synovial intimal fibroblasts-sub10")
p_detailed <- ggplot(plot_data_Synovial_intimal_fibroblasts_sub10, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.4, 0.65),
    breaks = seq(0.4, 0.65, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Synovial intimal fibroblasts-sub10",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Synovial_intimal_fibroblasts_sub10_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Synovial_intimal_fibroblasts_sub10_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_B_cells_24  <- subset(plot_data, Gene == "B cells-24")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_24, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1.8),
    breaks = seq(0, 1.8, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "B cells-24",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/B cells-24_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/B cells-24_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_B_cells_4  <- subset(plot_data, Gene == "B cells-4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_4, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1.8),
    breaks = seq(0, 1.8, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "B cells-4",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/B cells-4_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/B cells-4_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_B_cells_8  <- subset(plot_data, Gene == "B cells-8")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_8, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 4),
    breaks = seq(0, 4, by = 0.8)
  ) +
  # 标签和标题
  labs(
    title = "B cells-8",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/B cells-8_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/B cells-8_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Basophils_20  <- subset(plot_data, Gene == "Basophils-20")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Basophils_20, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Basophils-20",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Basophils_20_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Basophils_20_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Dendritic_cells_21  <- subset(plot_data, Gene == "Dendritic cells-21")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Dendritic_cells_21, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2),
    breaks = seq(0, 2, by = 0.6)
  ) +
  # 标签和标题
  labs(
    title = "Dendritic cells-21",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Dendritic_cells_21_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Dendritic_cells_21_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Endothelial_cells_19  <- subset(plot_data, Gene == "Endothelial cells-19")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Endothelial_cells_19, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(1.5, 4),
    breaks = seq(1.5, 4, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Endothelial cells-19",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Endothelial_cells_19_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Endothelial_cells_19_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Erythrocytes_17  <- subset(plot_data, Gene == "Erythrocytes-17")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Erythrocytes_17, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(9, 12),
    breaks = seq(9, 12, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Erythrocytes-17",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Erythrocytes_17_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Erythrocytes_17_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Erythrocytes_9  <- subset(plot_data, Gene == "Erythrocytes-9")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Erythrocytes_9, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(24, 29),
    breaks = seq(24, 29, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Erythrocytes-9",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Erythrocytes_9_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Erythrocytes_9_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Granulocytes_0  <- subset(plot_data, Gene == "Granulocytes-0")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_0, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(2, 5),
    breaks = seq(2, 5, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-0",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_0_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_0_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_14  <- subset(plot_data, Gene == "Granulocytes-14")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_14, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(7, 10),
    breaks = seq(7, 10, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-14",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_14_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_14_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_15  <- subset(plot_data, Gene == "Granulocytes-15")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_15, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.3, 1.5),
    breaks = seq(0.3, 1.5, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-15",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_15_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_15_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_2  <- subset(plot_data, Gene == "Granulocytes-2")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_2, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-2",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_2_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_2_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_3  <- subset(plot_data, Gene == "Granulocytes-3")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_3, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.7),
    breaks = seq(0, 0.7, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-3",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_3_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_3_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_Granulocytes_6  <- subset(plot_data, Gene == "Granulocytes-6")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_6, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-6",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_6_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_6_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_7  <- subset(plot_data, Gene == "Granulocytes-7")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_7, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(7, 9),
    breaks = seq(7, 9, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-7",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_7_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_7_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Granulocytes_sub4  <- subset(plot_data, Gene == "Granulocytes-sub4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_sub4, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.7),
    breaks = seq(0, 0.7, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-sub4",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_sub4_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Granulocytes_sub4_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Megakaryocytes_16  <- subset(plot_data, Gene == "Megakaryocytes-16")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Megakaryocytes_16, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(5, 12),
    breaks = seq(5, 12, by = 2)
  ) +
  # 标签和标题
  labs(
    title = "Megakaryocytes-16",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Megakaryocytes_16_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Megakaryocytes_16_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_Monocyte_macrophages_5  <- subset(plot_data, Gene == "Monocyte macrophages-5")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Monocyte_macrophages_5, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2),
    breaks = seq(0, 2, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Monocyte macrophages-5",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Monocyte_macrophages_5_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Monocyte_macrophages_5_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Reticular_cells_23  <- subset(plot_data, Gene == "Reticular cells-23")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Reticular_cells_23, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(4, 6),
    breaks = seq(4, 6, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Reticular cells-23",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Reticular_cells_23_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Reticular_cells_23_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Smooth_muscle_cells_25  <- subset(plot_data, Gene == "Smooth muscle cells-25")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Smooth_muscle_cells_25, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Smooth muscle cells-25",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/Smooth_muscle_cells_25_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/Smooth_muscle_cells_25_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_T_cells_12  <- subset(plot_data, Gene == "T cells-12")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_T_cells_12, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "T cells-12",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_other/Ratio/T_cells_12_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_other/Ratio/T_cells_12_plot.pdf", p_detailed, width = 10, height = 8)







saveRDS(plot_data, file = "RCTD_20250902/remove_other/Ratio/plot_data.rds")


















# remove Erythrocyte --------------
rctd_outdir_6_remove_Erythrocyte <- "RCTD_20250902/remove_Erythrocyte"
dir.create("RCTD_20250902/remove_Erythrocyte")


all_data <- readRDS("RCTD_20250902/all_data.rds")

library(spacexr)

# 2. single data
all_data_sc <- readRDS("RCTD_20250902/all_data_sc_250902.rds")

# 单细胞数据处理
Idents(all_data_sc) <- "annotation_new"
all_data_sc <- subset(all_data_sc, 
                      subset = annotation_new != "Erythrocytes-9" &
                               annotation_new != "Erythrocytes-17")

# all_data_sc <- all_data_sc[, all_data_sc$annotation_new != c("Erythrocytes-9", "Erythrocytes-17")]
table(all_data_sc$annotation_new)

count_sc_data <- readRDS('adata_sc_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")

count_sc_data <- subset(count_sc_data, cells = rownames(all_data_sc@meta.data))

counts_sc  <- count_sc_data %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation_new)
nUMI_sc    <- all_data_sc$total_counts

sc_cell_filter_idx <-  cluster_sc != "Unknown"
# sc_cell_filter_idx_new <-  !(cluster_sc %in% c("Unknown","Erythrocytes-9","Erythrocytes-17"))

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
visim_ls4 <-all_data %>% SplitObject(split.by = "sample")
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


saveRDS(rctd_result_ls, 'RCTD_20250902/remove_Erythrocyte/rctd_result_ls.rds')


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

saveRDS(all_data_rctd, file = "RCTD_20250902/remove_Erythrocyte/all_data_rctd.rds")


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
  file.path(rctd_outdir_6_remove_Erythrocyte, 'spatial plot of RCTD.pdf'), 
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
  file.path(rctd_outdir_6_remove_Erythrocyte, 'spatial plot of RCTD-1.pdf'), 
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

dir.create("RCTD_20250902/remove_Erythrocyte/Ratio", recursive = T)

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
write.table(out_gene_umi,file="RCTD_20250902/remove_Erythrocyte/Ratio/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.table("RCTD_20250902/remove_Erythrocyte/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(data)
data <- data[, c(1,5,27)]
colnames(data) <- c("Cells","Sample","Annotation")

# data_filt <- data[data$Annotation != "Other"]
data_filt <- data %>% filter(Annotation != "Other")
data <- data_filt

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250902/remove_Erythrocyte/Ratio/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

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
write.table(DATA,file="RCTD_20250902/remove_Erythrocyte/Ratio/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

samplemap <- read.table("RCTD_20250902/remove_Erythrocyte/Ratio/sample_map_group.txt", header=T, sep="\t", check.names=F)
colnames(samplemap)<-c("Sample","Group")
final_data <- data %>% left_join(samplemap, by = "Sample")

# write.table(final_data,file="RCTD_20250902/remove_Erythrocyte/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20250902/remove_Erythrocyte/Ratio/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

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
write.table(DATA,file="RCTD_20250902/remove_Erythrocyte/Ratio/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan1)+
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
  scale_fill_manual(values = plan1)+
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
ggsave(p,file="RCTD_20250902/remove_Erythrocyte/Ratio/ratio_barplot.pdf",width=15, height=6)
ggsave(p,file="RCTD_20250902/remove_Erythrocyte/Ratio/ratio_barplot.png",width=15, height=6,dpi=300)



# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = "RCTD_20250902/remove_Erythrocyte/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"

data_filt <- data %>% filter(Gene != "Other")
data <- data_filt

data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="RCTD_20250902/remove_Erythrocyte/Ratio/sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="RCTD_20250902/remove_Erythrocyte/Ratio/sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
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
ggsave(p,file="RCTD_20250902/remove_Erythrocyte/Ratio/diff_boxplot.pdf",width=20,height=5)
ggsave(p,file="RCTD_20250902/remove_Erythrocyte/Ratio/diff_boxplot.png",width=20,height=5)











# 对每种细胞类型进行正态性检验
shapiro_results <- plot_data %>%
  group_by(Gene) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(value)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(value)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "RCTD_20250902/remove_Erythrocyte/Ratio/shapiro_test_results.csv", row.names = FALSE)



plot_data_Chondrocytes  <- subset(plot_data, Gene == "Chondrocytes")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Chondrocytes, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 11),
    breaks = seq(0, 11, by = 2)
  ) +
  # 标签和标题
  labs(
    title = "Chondrocytes",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Chondrocytes_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Chondrocytes_plot.pdf", p_detailed, width = 10, height = 8)



plot_data_Osteoblasts <- subset(plot_data, Gene == "Osteoblasts")
p_detailed <- ggplot(plot_data_Osteoblasts, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
   # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(7, 11),
    breaks = seq(7, 11, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Osteoblasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Osteoblasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Osteoblasts_plot.pdf", p_detailed, width = 10, height = 8)





plot_data_Osteoclasts <- subset(plot_data, Gene == "Osteoclasts")
p_detailed <- ggplot(plot_data_Osteoclasts, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.1),
    breaks = seq(0, 2.1, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Osteoclasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Osteoclasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Osteoclasts_plot.pdf", p_detailed, width = 10, height = 8)





plot_data_Synovial_fibroblasts_sub2 <- subset(plot_data, Gene == "Synovial fibroblasts-sub2")
p_detailed <- ggplot(plot_data_Synovial_fibroblasts_sub2, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.2, 0.35),
    breaks = seq(0.2, 0.35, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Synovial fibroblasts-sub2",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Synovial_fibroblasts_sub2_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Synovial_fibroblasts_sub2_plot.pdf", p_detailed, width = 10, height = 8)






plot_data_Synovial_intimal_fibroblasts_sub10 <- subset(plot_data, Gene == "Synovial intimal fibroblasts-sub10")
p_detailed <- ggplot(plot_data_Synovial_intimal_fibroblasts_sub10, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.4, 0.65),
    breaks = seq(0.4, 0.65, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Synovial intimal fibroblasts-sub10",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Synovial_intimal_fibroblasts_sub10_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Synovial_intimal_fibroblasts_sub10_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_B_cells_24  <- subset(plot_data, Gene == "B cells-24")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_24, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "B cells-24",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/B cells-24_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/B cells-24_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_B_cells_4  <- subset(plot_data, Gene == "B cells-4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_4, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(1, 2.8),
    breaks = seq(1, 2.8, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "B cells-4",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/B cells-4_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/B cells-4_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_B_cells_8  <- subset(plot_data, Gene == "B cells-8")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_8, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 8),
    breaks = seq(0, 8, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "B cells-8",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/B cells-8_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/B cells-8_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Basophils_20  <- subset(plot_data, Gene == "Basophils-20")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Basophils_20, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Basophils-20",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Basophils_20_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Basophils_20_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Dendritic_cells_21  <- subset(plot_data, Gene == "Dendritic cells-21")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Dendritic_cells_21, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.6, 2),
    breaks = seq(0.6, 2, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Dendritic cells-21",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Dendritic_cells_21_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Dendritic_cells_21_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Endothelial_cells_19  <- subset(plot_data, Gene == "Endothelial cells-19")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Endothelial_cells_19, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(2.4, 3),
    breaks = seq(2.4, 3, by = 0.1)
  ) +
  # 标签和标题
  labs(
    title = "Endothelial cells-19",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Endothelial_cells_19_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Endothelial_cells_19_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Granulocytes_0  <- subset(plot_data, Gene == "Granulocytes-0")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_0, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(18, 25),
    breaks = seq(18, 25, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-0",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_0_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_0_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_14  <- subset(plot_data, Gene == "Granulocytes-14")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_14, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(8, 11),
    breaks = seq(8, 11, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-14",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_14_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_14_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_15  <- subset(plot_data, Gene == "Granulocytes-15")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_15, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.6, 2),
    breaks = seq(0.6, 2, by = 0.4)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-15",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_15_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_15_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_2  <- subset(plot_data, Gene == "Granulocytes-2")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_2, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-2",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_2_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_2_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_3  <- subset(plot_data, Gene == "Granulocytes-3")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_3, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-3",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_3_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_3_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_Granulocytes_6  <- subset(plot_data, Gene == "Granulocytes-6")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_6, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-6",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_6_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_6_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_7  <- subset(plot_data, Gene == "Granulocytes-7")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_7, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(15, 18),
    breaks = seq(15, 18, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-7",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_7_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_7_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Granulocytes_sub4  <- subset(plot_data, Gene == "Granulocytes-sub4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_sub4, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-sub4",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_sub4_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Granulocytes_sub4_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Megakaryocytes_16  <- subset(plot_data, Gene == "Megakaryocytes-16")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Megakaryocytes_16, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(5, 14),
    breaks = seq(5, 14, by = 2)
  ) +
  # 标签和标题
  labs(
    title = "Megakaryocytes-16",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Megakaryocytes_16_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Megakaryocytes_16_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_Monocyte_macrophages_5  <- subset(plot_data, Gene == "Monocyte macrophages-5")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Monocyte_macrophages_5, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Monocyte macrophages-5",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Monocyte_macrophages_5_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Monocyte_macrophages_5_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Reticular_cells_23  <- subset(plot_data, Gene == "Reticular cells-23")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Reticular_cells_23, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(3, 5),
    breaks = seq(3, 5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Reticular cells-23",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Reticular_cells_23_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Reticular_cells_23_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Smooth_muscle_cells_25  <- subset(plot_data, Gene == "Smooth muscle cells-25")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Smooth_muscle_cells_25, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Smooth muscle cells-25",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Smooth_muscle_cells_25_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/Smooth_muscle_cells_25_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_T_cells_12  <- subset(plot_data, Gene == "T cells-12")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_T_cells_12, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#00FFFF",    # 青色
      "4M" = "#00008B"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "T cells-12",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/T_cells_12_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20250902/remove_Erythrocyte/Ratio/T_cells_12_plot.pdf", p_detailed, width = 10, height = 8)




saveRDS(plot_data, file = "RCTD_20250902/remove_Erythrocyte/Ratio/plot_data.rds")




































# 处理数据：将两次rctd的结果添加到all_data的meta.data中 ----- 
all_data <- readRDS(file = "data_for_analysis/all_data_projected.rds")


### rctd4 ------
meta_data_rctd4 <- read.table("RCTD_4/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(meta_data_rctd4) <- gsub("first_type", "rctd4", colnames(meta_data_rctd4))
meta_data_rctd4$Cell <- gsub("24M_s_", "s_", meta_data_rctd4$Cell)
meta_data_rctd4$Cell <- gsub("4M_s_", "s_", meta_data_rctd4$Cell)

meta_data_rctd4 <- meta_data_rctd4 %>% 
  mutate(new_rownames = paste(ROI, Cell, sep = "_")) %>% 
  column_to_rownames("new_rownames")

head(meta_data_rctd4)

all_data <- all_data %>% 
  AddMetaData(
    meta_data_rctd4 %>% select(rctd4)
  )



### rctd remove other ------
meta_data_rctd_20250902 <- read.table("RCTD_20250902/remove_other/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(meta_data_rctd_20250902) <- gsub("first_type", "rctd_20250902", colnames(meta_data_rctd_20250902))

meta_data_rctd_20250902 <- meta_data_rctd_20250902 %>% 
  column_to_rownames("Cell")

head(meta_data_rctd_20250902)

all_data <- all_data %>% 
  AddMetaData(
    meta_data_rctd_20250902 %>% select(rctd_20250902)
  )
head(all_data@meta.data)



### rctd remove erythrocyte other ------
meta_data_rctd_remove_erythrocyte_20250902 <- read.table("RCTD_20250902/remove_Erythrocyte/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(meta_data_rctd_remove_erythrocyte_20250902) <- gsub("first_type", "rctd_remove_erythrocyte_20250902", colnames(meta_data_rctd_remove_erythrocyte_20250902))

meta_data_rctd_remove_erythrocyte_20250902 <- meta_data_rctd_remove_erythrocyte_20250902 %>% 
  column_to_rownames("Cell")

head(meta_data_rctd_remove_erythrocyte_20250902)

all_data <- all_data %>% 
  AddMetaData(
    meta_data_rctd_remove_erythrocyte_20250902 %>% select(rctd_remove_erythrocyte_20250902)
  )
head(all_data@meta.data)




# 保存数据 ------- 
saveRDS(all_data, file = "data_for_analysis/all_data_projected_rctd_20250902.rds")






# rctd remove other 注释信息导出 ------ 
rctd_remore_other <- 
  all_data[[]] %>% 
  rownames_to_column("Barcode") %>% 
  select(Barcode, rctd_20250902, sample, group)
head(rctd_remore_other)

project_groups <- split(rctd_remore_other, rctd_remore_other$sample)

# 为每个 project 创建 CSV 文件
walk(names(project_groups), function(proj) {
  # 提取当前 project 的数据
  proj_data <- project_groups[[proj]] %>%
    select(Barcode, CellType = rctd_20250902) %>% 
    mutate(Barcode = str_remove(Barcode, paste0(proj, "_")))# 重命名列
  
  # 创建文件名
  filename <- paste0("RCTD_20250902/remove_other/", proj, "_cell_annotations.csv")
  
  # 保存为 CSV
  write.csv(proj_data, file = filename, row.names = FALSE)
  
  # 输出信息
  cat("已保存:", filename, "-", nrow(proj_data), "个细胞\n")
})




# rctd remove other erythrocyte 注释信息导出 -------
rctd_remove_errythrocyte <- 
  all_data[[]] %>% 
  rownames_to_column("Barcode") %>% 
  select(Barcode, rctd_remove_erythrocyte_20250902, sample)

project_groups <- split(rctd_remove_errythrocyte, rctd_remove_errythrocyte$sample)

# 为每个 project 创建 CSV 文件
walk(names(project_groups), function(proj) {
  # 提取当前 project 的数据
  proj_data <- project_groups[[proj]] %>%
    select(Barcode, CellType = rctd_remove_erythrocyte_20250902) %>% 
    mutate(Barcode = str_remove(Barcode, paste0(proj, "_")))# 重命名列
  
  # 创建文件名
  filename <- paste0("RCTD_20250902/remove_Erythrocyte/", proj, "_cell_annotations.csv")
  
  # 保存为 CSV
  write.csv(proj_data, file = filename, row.names = FALSE)
  
  # 输出信息
  cat("已保存:", filename, "-", nrow(proj_data), "个细胞\n")
})



