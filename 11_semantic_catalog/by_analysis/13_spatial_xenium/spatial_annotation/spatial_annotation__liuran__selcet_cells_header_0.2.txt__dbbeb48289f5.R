
library(Seurat)
library(SeuratObject)
all_data_sub <- readRDS("all_data_sub.rds")
#all_data_sub<-UpdateSeuratObject(all_data_sub)
table(all_data_sub@meta.data$annotation_sub)
allCells = names(Idents(all_data_sub))
allType = levels(Idents(all_data_sub))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(all_data_sub)== x ]
  num = round(length(cgCells)*0.2)
  cg=sample(cgCells,num)
  cg
}))

all_data_sub <- subset(all_data_sub, cells = choose_Cells)
table(all_data_sub@meta.data$annotation_sub)

out_gene_umi = as.matrix(all_data_sub@meta.data)
row_name = c("Cell", rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
write.table(out_gene_umi, file = "extract_cell/selcet_cells_header_0.2.txt", sep = "\t", quote = F, row.names = F, col.names = F)
