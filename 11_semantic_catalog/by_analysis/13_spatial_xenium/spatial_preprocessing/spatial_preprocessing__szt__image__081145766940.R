# match table - id
images_ls <-
  file_list %>%
  imap(function(path, name){
    path <- file.path(path, "outs")
    spatial_path <- file.path(path, "binned_outputs", str_glue("square_{bin.size}"), 'spatial')
    
    # read hi-res/low-res
    image <- Read10X_Image(
      spatial_path,
      assay      = assay.name,
      slice      = name,
      image.name = image.name,
    )
  })



#str_remove（barcode，"^.*_s_"）
library(dplyr)
table_match<-all_data@meta.data[which(all_data$orig.ident=='DZZ_Saline'),]%>%
  rownames_to_column('barcode')%>% 
  dplyr::select('barcode')%>% 
  mutate(raw_id= gsub("^.*?_s", "s",.$barcode))%>%
  column_to_rownames('raw_id')

images_ls$`DZZ_Saline`<-subset(images_ls$`DZZ_Saline`,cells=rownames(table_match))
table_match<-table_match[Cells(images_ls$`DZZ_Saline`),]
images_ls$`DZZ_Saline`<-RenameCells(images_ls$`DZZ_Saline`,new.names=table_match)

table_match<-all_data@meta.data[which(all_data$orig.ident=='SYZ_NH4'),]%>%
  rownames_to_column('barcode')%>% 
  dplyr::select('barcode')%>% 
  mutate(raw_id= gsub("^.*?_s", "s",.$barcode))%>%
  column_to_rownames('raw_id')

images_ls$`SYZ_NH4`<-subset(images_ls$`SYZ_NH4`,cells=rownames(table_match))
table_match<-table_match[Cells(images_ls$`SYZ_NH4`),]
images_ls$`SYZ_NH4`<-RenameCells(images_ls$`SYZ_NH4`,new.names=table_match)











'Slc17a6','Slc32a1','Slc6a5','Adcyap1','Phox2b','Dbh','Cck','Th','Tac1','Pdyn','Pvalb','P2ry1','Crhr2','Npy2r','Gpr65','Ascl1','Nmbr','Calca','Trpv1','Calb1','Calb2','Tac2','Htra1','Sst','Vip','Cbln1','Chat','Oprm1'
