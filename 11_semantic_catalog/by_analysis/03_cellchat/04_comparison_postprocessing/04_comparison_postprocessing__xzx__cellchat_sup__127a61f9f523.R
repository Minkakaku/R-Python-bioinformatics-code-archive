# options(encoding = "UTF-8")
# conda activate "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/"
# debug:shenhao

library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(argparser)
library(pdftools)
library(stringr)
library(jpeg)
library(rmarkdown)
# install.packages('NMF')
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("jokergoo/ComplexHeatmap")
# devtools::install_github("sqjin/CellChat")
# devtools::install_github("jinworks/CellChat") 这个是新版的，要装这个
# BiocManager::install("ggalluvial")
# library(dior)
# ibrary(SingleR)
# library('getopt')
# library(tidyverse)
# options(stringsAsFactors = FALSE)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="input rds")
argv <- add_argument(argv,"--prefix", help="group name prefix",default='F')
argv <- add_argument(argv,"--outpath", help="the output dir")
argv <- add_argument(argv,"--splitby", help="the data split.by",default='group')
argv <- add_argument(argv,"--subcelltype", help="the cell_type used",default='all')
argv <- add_argument(argv,"--comparison", help="Group or sample for comparison",default='F')
argv <- add_argument(argv,"--species", help="species used,human or mouse",default='human')
argv <- parse_args(argv)


#rds <- argv$rds
#prefix <- argv$prefix
outpath <- argv$outpath
if(!dir.exists(outpath)) {dir.create(outpath, recursive = T)}
#splitby <- argv$splitby
#subcelltype <- unlist(strsplit(argv$subcelltype, ","))
# comparison <- "A:B,C:D,E:F"
comparison <- unlist(strsplit(argv$comparison, ","))
#species <- argv$species

outdir = argv$outpath
load(paste0(outpath,"/tmp_object.RData"))
if(argv$comparison=="F"){
    pathway.union_intersection <- object.list[[1]]@netP$pathways
    }else{
        pathway.union_intersection <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
        }
height_h <- ceiling(length(pathway.union_intersection) / 9) * 10

#提前生成标题
if(argv$comparison=="F"){
  cellchat <- object.list[1]
  title_Number <- paste0("Number of interactions - ",names(object.list)[1])
  title_Weight <- paste0("Weight of interactions - ",names(object.list)[1])
  title_name <- names(object.list)[1]
}else{
  cellchat <- mergeCellChat(object.list,add.names = names(object.list))
  title_Number <- c(paste0("Number of interactions - ",names(object.list)[1]),paste0("Number of interactions - ",names(object.list)[2]))
  title_Weight <- c(paste0("Weight of interactions - ",names(object.list)[1]),paste0("Weight of interactions - ",names(object.list)[2]))
  title_name <- c(names(object.list)[1],names(object.list)[2])
}

for(i in c("outgoing","incoming")){
    if(argv$comparison=="F"){
        ht <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = i, signaling = pathway.union_intersection,title = title_name[1],width = 8, height = height_h)
        }else{
            ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = i, signaling = pathway.union_intersection,title = title_name[1],width = 8, height = height_h)
			ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = i,signaling = pathway.union_intersection,title = title_name[2],width = 8, height = height_h)
			ht <- ht1 + ht2
			}
    pdf(paste0(outdir,i,".pdf"),width = 9,height = 0.5*height_h+1)
	print(ht)
	dev.off()
	}

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
    }
cat("Processed", length(pdf_files), "PDF files.\n")
cp_parts <- unlist(strsplit(comparison[cp], ":"))
if(argv$comparison=="F"){
    params<-list(Group1=names(cellchat_sp))
    file.copy("/mnt/NexData_file/share/scRNA-tools/report/pip_report_cellchat_single.rmd",outdir)
    render(input = paste0(outdir,"pip_report_cellchat_single.rmd"),output_dir=outdir,params=params)
    }else{
        params<-list(Group1=cp_parts[1],Group2=cp_parts[2])
        file.copy("/mnt/NexData_file/share/scRNA-tools/report/pip_report_cellchat.rmd",outdir)
        render(input = paste0(outdir,"pip_report_cellchat.rmd"),output_dir=outdir,params=params)
        }
