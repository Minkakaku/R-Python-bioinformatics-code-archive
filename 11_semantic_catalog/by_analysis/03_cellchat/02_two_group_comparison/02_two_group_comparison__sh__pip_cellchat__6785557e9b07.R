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
argv <- add_argument(argv,"--subgroup", help="the group used",default='all')
argv <- add_argument(argv,"--comparison", help="Group or sample for comparison",default='F')
argv <- add_argument(argv,"--species", help="species used,human or mouse",default='human')
argv <- parse_args(argv)

rds <- argv$rds
prefix <- argv$prefix
outpath <- argv$outpath
if(!dir.exists(outpath)) {dir.create(outpath, recursive = T)}
splitby <- argv$splitby
subcelltype <- unlist(strsplit(argv$subcelltype, ","))
subgroup <- unlist(strsplit(argv$subgroup, ","))
# comparison <- "A:B,C:D,E:F"
comparison <- unlist(strsplit(argv$comparison, ","))
species <- argv$species

#' Run CellChat analysis
#'
#' This function performs CellChat analysis on the given cellchat object.
#'
#' @param cellchat A CellChat object.
#' @param name_ct A character string representing the name of the cell type.
#' @param species A character string indicating the species ("human" or "mouse").
#' @param outdir A character string specifying the output directory.
#' @return A CellChat object with updated analysis results.
runCellchat = function(cellchat, name_ct, species, outdir){
  # showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
  cellchat <- setIdent(cellchat, ident.use = "cell_type") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
  CellChatDB <- if (species == "human") CellChatDB.human else CellChatDB.mouse
  CellChatDB.use <- subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
  cellchat@DB <- CellChatDB.use
  cellchat@idents <- droplevels(cellchat@idents)
  unique(cellchat@idents)
  cellchat@data.signaling 
  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
  cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
  # project gene expression data onto PPI network (optional)
  cellchat <- if (species == "human") smoothData(cellchat, adj = PPI.human) else smoothData(cellchat, adj = PPI.mouse)
  cellchat <- computeCommunProb(cellchat,type="triMean",trim=0,raw.use = TRUE)###这应该是核心的一行代码
  cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
  cellchat <- computeCommunProbPathway(cellchat)
  
  df.net <- subsetCommunication(cellchat)
  relation <- paste0(df.net[,1],"->",df.net[,2])
  df.net <- cbind(relation,df.net)
  colnames(df.net)[1] <- "relation"
  write.table(df.net,paste0(outdir,"/",name_ct,"-net_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  df.netp <- subsetCommunication(cellchat,slot.name="netP")
  write.table(df.netp,paste0(outdir,"/",name_ct,"-netPathways_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}

# 分组生成cellchat对象
process_cellchat_comparison <- function(comparison_sub, cellchat_sp, cp, species,outdir) {
    if(comparison_sub=="F"){cp_parts <- names(cellchat_sp)}else{cp_parts <- unlist(strsplit(comparison_sub, ":"))}
		common_cluster <- lapply(cp_parts, function(x) cellchat_sp[[x]]@meta.data$cell_type)
		common_clusters <- Reduce(intersect, common_cluster)
    meta_data <- lapply(cp_parts, function(x) {
			meta_data_x <- cellchat_sp[[x]]@meta.data
			meta_data_x[meta_data_x$cell_type %in% common_clusters, ]
		})
    cellchats <- lapply(seq_along(cp_parts), function(i) {
        createCellChat(object = subset(cellchat_sp[[cp_parts[i]]],idents=common_clusters), meta = meta_data[[i]], group.by = "cell_type") %>%
            runCellchat(cp_parts[i], species,outdir)
    })
    setNames(cellchats, cp_parts)
}

vertex_cal = function (object, signaling, slot.name = "netP"){
		centr0 <- slot(object, slot.name)$centr[signaling]
		centr0 <- centr0[[1]]
    mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
    mat <- mat[c("outdeg","indeg"),]
		mat <- as.data.frame(t(mat))
		top3_indices <- head(order(mat$outdeg, decreasing = TRUE), 3)
}

Seuset <- readRDS(rds)
# Seuset <- readRDS("/home/sh/sc_pipline/cellchat/test/FW2024-225.rds")
# levels(Seuset@meta.data$group)[2]<-'ADM'
Idents(Seuset)<-Seuset@meta.data$cell_type
if (length(subcelltype) > 1) {
  Seuset<-subset(Seuset,idents=subcelltype) #"CD8TN","CD8TCM","CD8TEM","CD8TEX"
  Seuset@meta.data$cell_type<-Idents(Seuset)
}
# if (length(subgroup) > 1) {
  Seuset<-subset(Seuset,group %in% subgroup)
  print(levels(Seuset$group))
  # Seuset@meta.data$cell_type<-Idents(Seuset)
# }
common_cell_types <- levels(Seuset) #用来判断分细胞类型弦图的输出
width_heat <- 8 #热图的比例
if(argv$comparison != "F"){
	common_cell_types <- c()
	for(i in 1:length(comparison)){
		comparison_sub <- unlist(strsplit(comparison[i], ":"))
		groups <- split(as.character(Seuset$cell_type[which(Seuset$group %in% comparison_sub)]), as.character(Seuset$group[which(Seuset$group %in% comparison_sub)]))
		common_cell_type <- Reduce(intersect, groups)
		common_cell_types<-c(common_cell_types,common_cell_type)
	}
	common_cell_types<-unique(common_cell_types)
	Seuset<-subset(Seuset,idents=common_cell_types)
	width_heat <- 12
}
if(length(common_cell_types)<=6){labeledge<-"T"}else{labeledge<-"F"} #弦图的数字展不展示
# 计算唯一细胞类型的数量
num_celltypes <- length(common_cell_types)
if (num_celltypes <= 4) n <- c(1,num_celltypes) else n <- c(ceiling(num_celltypes / 3),3) #用来计算分细胞类型弦图的排版
Seuset@meta.data$cell_type<-Idents(Seuset)
if(argv$comparison != "F"){
	cellchat_sp <- SplitObject(Seuset,split.by = splitby)
	}else{
		cellchat_sp <- list(all=Seuset)
		}
measure_list <- c("count","weight") #循环省事
heatmap_color <- c("white","red2") #热图的色块

for (cp in 1:length(comparison)) {
  if(argv$comparison=="F"){
    outdir <- paste0(outpath,"/")
    }else{
      outdir<-paste0(outpath,"/",str_replace(comparison[cp],":","_"),"/")
      dir.create(outdir, recursive = T)
      }
  object.list <- process_cellchat_comparison(comparison[cp], cellchat_sp, cp = cp, species = species, outdir = outdir)
  # levels1 <- levels(droplevels(object.list[[1]]@meta$cell_type))
  # levels2 <- levels(droplevels(object.list[[2]]@meta$cell_type))
  # if(!identical(levels1,levels2)){
  #   levels_inter <- intersect(levels1,levels2)
  #   object.list[[1]] <- subsetCellChat(object.list[[1]],idents.use =levels_inter) 这个subset函数会改数据格式
  #   object.list[[2]] <- subsetCellChat(object.list[[2]],idents.use =levels_inter)
  #   }
	save(object.list,file=paste0(outdir,"tmp_object.RData"))
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
	title_list <- list(title_Number,title_Weight)
	## 细胞通讯网络可视化
	#分组细胞类型count数量热图
	for (i in seq_along(measure_list)) {
		gg1 = netVisual_heatmap(object.list[[1]], measure = measure_list[i],color.heatmap = heatmap_color,title.name = title_list[[i]][1])
		if(argv$comparison=="F"){
			gg <- gg1
			}else{
				gg2 = netVisual_heatmap(object.list[[2]], measure = measure_list[i],color.heatmap = heatmap_color,title.name = title_list[[i]][2])
				gg <- gg1+gg2
				}
		if(i==1){label <- "Number"}else{label <- "Weight"}
		pdf(paste0(outdir,"group_",label,"_netVisual_heatmap.pdf"),width = width_heat,height = 7)
		print(gg)
		dev.off()
		}

	# if(length(common_cell_types)<=12){
		#分组细胞类型相互作用数量弦图，以最大count做scale
		count.max <- getMaxWeight(object.list, attribute = c("idents","count"))
		weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
		if(argv$comparison=="F"){
			pdf(paste0(outdir,"Number_of_interactions_netplot.pdf"),width = 8,height = 7)
			netVisual_circle(object.list[[1]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number)
			dev.off()
			#分组细胞类型相互作用权重弦图，以最大weight做scale
			pdf(paste0(outdir,"Weight_of_interactions_netplot.pdf"),width = 8,height = 7)
			netVisual_circle(object.list[[1]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight)
			dev.off()
			}else{
				pdf(paste0(outdir,"Number_of_interactions_netplot.pdf"),width = 12,height = 7)
				par(mfrow = c(1,2), xpd=TRUE)
				for (i in 1:length(object.list)) {
					netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
					}
				dev.off()
				#分组细胞类型相互作用权重弦图，以最大weight做scale
				pdf(paste0(outdir,"Weight_of_interactions_netplot.pdf"),width = 12,height = 7)
				par(mfrow = c(1,2), xpd=TRUE)
				for (i in 1:length(object.list)) {
					netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
					}
				dev.off()
				}
		# }

	# if(length(common_cell_types)<=12){
		if(argv$comparison=="F"){
			mat_count = cellchat[[1]]@net$count
			mat_weight = cellchat[[1]]@net$weight
			pdf(paste0(outdir,"",names(cellchat_sp),"_Number_netVisual_circle.pdf"),width = 5*n[2],height = 5*n[1])
			par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
			for (i in 1:nrow(mat_count)) {
				 mat_count2 = matrix(0,nrow = nrow(mat_count),ncol = ncol(mat_count),dimnames = dimnames(mat_count))
				 mat_count2[i,] = mat_count[i,]
				 # mat_count2<-mat_count2[1:5,1:5]
				 netVisual_circle(mat_count2,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_count),edge.label.cex = 1,title.name = rownames(mat_count)[i])
      }
      dev.off()
			pdf(paste0(outdir,"",names(cellchat_sp),"_Weight_netVisual_circle.pdf"),width = 5*n[2],height = 5*n[1])
			par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
			for (i in 1:nrow(mat_count)) {
				mat_weight2 = matrix(0,nrow = nrow(mat_weight),ncol = ncol(mat_weight),dimnames = dimnames(mat_weight))
				mat_weight2[i,] = mat_weight[i,]
				netVisual_circle(mat_weight2,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_weight),edge.label.cex = 1,title.name = rownames(mat_weight)[i])
      }
      dev.off()
		}else{
			#单细胞类型贝壳图，感觉比较有用
			for(j in 1:2){
				mat_count = cellchat@net[[j]]$count
				mat_weight = cellchat@net[[j]]$weight
				pdf(paste0(outdir,"",names(cellchat@net)[j],"_Number_netVisual_circle.pdf"),width = 5*n[2],height = 5*n[1])
				par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
				for (i in 1:nrow(mat_count)) {
					mat_count2 = matrix(0,nrow = nrow(mat_count),ncol = ncol(mat_count),dimnames = dimnames(mat_count))
					mat_count2[i,] = mat_count[i,]
					netVisual_circle(mat_count2,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_count),edge.label.cex = 1,title.name = rownames(mat_count)[i])
				}
				dev.off()
				pdf(paste0(outdir,"",names(cellchat@net)[j],"_Weight_netVisual_circle.pdf"),width = 5*n[2],height = 5*n[1])
				par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
				for (i in 1:nrow(mat_count)) {
					mat_weight2 = matrix(0,nrow = nrow(mat_weight),ncol = ncol(mat_weight),dimnames = dimnames(mat_weight))
					mat_weight2[i,] = mat_weight[i,]
					netVisual_circle(mat_weight2,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_weight),edge.label.cex = 1,title.name = rownames(mat_weight)[i])
				}
				dev.off()
				}
			}
		# }

	num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	weight.MinMax <- c(min(num.link), max(num.link))
	if(argv$comparison=="F"){
		pdf(paste0(outdir,"Interactions_compare_scatterplot.pdf"),width = 8,height = 7)
		print(netAnalysis_signalingRole_scatter(object.list[[1]], title = title_name[1], weight.MinMax = weight.MinMax))
		dev.off()
	}else{
		gg <- list()
		for (i in 1:length(object.list)) {
  		gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = title_name[i], weight.MinMax = weight.MinMax)
			}
		maxx <- ceiling(max(max(gg[[1]]$data$x),max(gg[[2]]$data$x)) / 5)*5
		maxy <- ceiling(max(max(gg[[1]]$data$y),max(gg[[2]]$data$y)) / 5)*5
		for (i in 1:length(object.list)) {
  		gg[[i]] <- gg[[i]]+ylim(0, maxy)+xlim(0, maxx)
			}
		pdf(paste0(outdir,"Interactions_compare_scatterplot.pdf"),width = 12,height = 7)
		print(patchwork::wrap_plots(plots = gg))
		dev.off()
		}
	## 运行上述代码出现报错 Error in plot.new() : figure margins too large
	# par("mar")
	# par(mar=c(1,1,1,1))
	# 重新运行上面的代码

	if(argv$comparison!="F"){
		#比较交互总数和交互强度分组柱状图，这个高宽可以固定
		cellchat
		gg1 <- compareInteractions(cellchat, show.legend = F, group = c(2, 1),measure = "count")
		gg2 <- compareInteractions(cellchat, show.legend = F, group = c(2, 1), measure = "weight")
		gg <- gg1 + gg2
		pdf(paste0(outdir,"Interactions_compare_barplot.pdf"),width = 8,height = 7)
		print(gg)
		dev.off()

		#相互作用的差异数量或差异作用强度细胞类型热图，与弦图对应，高宽需要计算一下，只跟细胞类型数量有关，14：12/7
		gg1 <- netVisual_heatmap(cellchat, comparison = c(2, 1))
		gg2 <- netVisual_heatmap(cellchat, comparison = c(2, 1), measure = "weight")
		gg <- gg1 + gg2
		pdf(paste0(outdir,"Interactions_diff_heatmap.pdf"),width = 12,height = 7)
		print(gg)
		dev.off()

		#相互作用的差异数量或差异作用强度细胞类型弦图，这个高宽可以固定
		# if(length(common_cell_types)<=12){
	  	pdf(paste0(outdir,"Interactions_diff_netplot.pdf"),width = 12,height = 7)
	  	par(mfrow = c(1,2), xpd=TRUE)
			for(i in c("count","weight")){
				netVisual_diffInteraction(cellchat, comparison = c(2, 1), weight.scale = T,measure = i)
				}
			dev.off()
			# }
		}

	#比较每个细胞群相关的传出传入信号,应该跟信息流flow放在一起+后面添加一下自动计算横纵比例
	# combining all the identified signaling pathways from different datasets 
	if(argv$comparison=="F"){
		pathway.union_intersection <- object.list[[1]]@netP$pathways
		}else{
			pathway.union_intersection <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
			}
	height_h <- ceiling(length(pathway.union_intersection) / 10) * 2.5
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

	select_pathway <- c()
	if(argv$comparison!="F"){
		for (i in seq_along(measure_list)) {
			if(i==1){label <- "Number"}else{label <- "Weight"}
			# Flow 信息流
			gg = rankNet(cellchat,mode = "comparison",measure = measure_list[i],sources.use = NULL,targets.use = NULL, comparison =c(1, 2) ,stacked = T,do.stat = TRUE)
			height_f <- ceiling(nrow(gg$data) / 20)
			pdf(paste0(outdir,"Relative_",label,"_information_flow.pdf"),height= height_f)
			print(gg)
			dev.off()
			select_pathway <- c(select_pathway,gg$data %>% filter(pvalues <= 0.05) %>% rownames())
			gg = rankNet(cellchat,mode = "comparison",measure = measure_list[i],sources.use = NULL,targets.use = NULL, comparison = c(1, 2),stacked = F,do.stat = TRUE)
			pdf(paste0(outdir,label,"_Information_flow.pdf"),height= height_f)
			print(gg)
			dev.off()
			}
		select_pathway <- select_pathway %>% table() %>% as.data.frame() %>% filter(Freq > 1) %>% .$.
		write.table(select_pathway,paste0(outdir,"select_pathway.txt"),sep="\t",quote=F,row.names=F,col.names=T)
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
	}