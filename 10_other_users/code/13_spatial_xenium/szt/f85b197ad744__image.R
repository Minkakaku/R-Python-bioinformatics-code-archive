# match table - id
table_match<-all_data@meta.data[which(all_data$orig.ident=='WT'),]%>%
  rownames_to_column('barcode')%>% 
  select('barcode')%>% 
  mutate(raw_id= str_sub(.$barcode,5,-1))%>%
  column_to_rownames('raw_id')

images_ls$`A1-WT-101`<-subset(images_ls$`A1-WT-101`,cells=rownames(table_match))
table_match<-table_match[Cells(images_ls$`A1-WT-101`),]
images_ls$`A1-WT-101`<-RenameCells(images_ls$`A1-WT-101`,new.names=table_match)

table_match<-all_data@meta.data[which(all_data$orig.ident=='KA'),]%>%
  rownames_to_column('barcode')%>% 
  select('barcode')%>% 
  mutate(raw_id= str_sub(.$barcode,5,-1))%>%
  column_to_rownames('raw_id')

images_ls$`KA1-9`<-subset(images_ls$`KA1-9`,cells=rownames(table_match))
table_match<-table_match[Cells(images_ls$`KA1-9`),]
images_ls$`KA1-9`<-RenameCells(images_ls$`KA1-9`,new.names=table_match)




                        

  