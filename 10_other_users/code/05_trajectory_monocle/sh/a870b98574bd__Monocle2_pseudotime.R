# Load necessary libraries
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(argparser)
library(rmarkdown)

# Define command line arguments
argv <- arg_parser('')
argv <- add_argument(argv, "--rds", type = "character", help = "Path to input rds")
argv <- add_argument(argv, "--gene", default = 'F', type = "character", help = "Path to gene file (tab-delimited text)")
argv <- add_argument(argv, "--lambda", default = 8, type = "numeric", help = "Lambda value (numeric)")
argv <- add_argument(argv, "--outdir", type = "character", help = "Output path")
argv <- add_argument(argv, "--RootState", default = 'F', type = "numeric", help = "Root state value (numeric)")
argv <- add_argument(argv, "--Reverse", default = 'F', type = "logical", help = "Reverse flag (logical)")
argv <- add_argument(argv, "--celllist", default = NULL, type = "character", help = "Path to cell list file")
argv <- add_argument(argv, "--clusters", default = "F", type = "character", help = "Cluster information")
argv <- add_argument(argv, "--RandomCells", default = "F", help = "Random cells flag (logical)")
argv <- add_argument(argv, "--CellsPercent", default = 0.1, type = "numeric", help = "Percentage of cells (numeric)")
argv <- add_argument(argv, "--MinCellsPerCluster", default = 1000, type = "numeric", help = "Minimum number of cells per cluster (numeric)")
args <- parse_args(argv)
options(stringsAsFactors = FALSE, warn = -1)

# Assign arguments to variables
rdsFile = args$rds
geneFile = args$gene; print(paste0("gene use:", head(geneFile)))
outdir = args$outdir; print(paste0("out path:", outdir))
lambda = args$lambda; print(paste0("lambda use:", lambda))
RootState = args$RootState
Reverse = args$Reverse
celllist = args$celllist
cluster = args$clusters
RandomCells = args$RandomCells; print(paste0("Random Cell:", RandomCells))
cellsPercent = args$CellsPercent; print(paste0("cells Percent:", cellsPercent))
minCellsPerCluster = args$MinCellsPerCluster; print(paste0("min Cells PerCluster:", minCellsPerCluster))

# Create output directory and set working directory
dir.create(outdir)
# dir.create(paste0(outdir, "Pseudotime"))
# setwd(paste0(outdir, "Pseudotime"))

# rm(list = ls())
# 函数名: Psepre
# 参数: seuset - Seurat对象 celllists -抽取的细胞列表 cluster -选择细胞类型 RandomCells -随机抽取
# 返回值: seurat - Seurat对象
# 功能: 对传入的seuset进行处理并返回
Psepre = function(seuset, outdir, celllist, cluster, RandomCells, cellsPercent, minCellsPerCluster) {
	# if(!is.null(celllist)){
	#   newcluster = read.table(celllist,sep="\t",header=F)
	#   realcluster = newcluster[is.finite(match(newcluster[,1],colnames(seuset))),]
	#   seuset = subset(seuset,cells=as.character(realcluster[,1]))
	#   seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]),
	#                     as.character(realcluster[,2]))
	# }else if(!is.null(cluster)){
	#   clusters = unlist(strsplit(cluster, "[,]"))
	#   celllist = names(seuset@active.ident[is.finite(match(seuset@active.ident,clusters))])
	#   seuset = subset(seuset,cells=celllist)
	# }
	if (cluster!="F") { 
		cluster = unlist(strsplit(cluster, ","))
		seuset = subset(seuset, idents = c(cluster))
		seuset@meta.data$cell_type <- droplevels(seuset@meta.data$cell_type)
		Idents(seuset) <- droplevels(Idents(seuset))
		print(paste0("cluster use:", levels(seuset)))
		}
	if (RandomCells != "F") {
		newcelllist = randomSample(seuset@active.ident, cellsPercent = 0.1, minCellsPerCluster = 1000)
		seuset = subset(seuset, cells = as.character(newcelllist))
		write.table(data.frame(Cell = colnames(seuset), Cluster = seuset@active.ident),
								file = paste0(outdir,"/","RandomCells.txt"), sep = "\t", row.names = F, quote = F)
	}
	seuset <- AddMetaData(object = seuset, metadata = seuset@active.ident, col.name = "Cluster")
	return(seuset)
}

# 函数名: setNum
# 参数: number - 传入细胞数 cellsPercent -抽取的比例 minCellsPerCluster -最小细胞数
# 返回值: finalNum - 需要抽取的细胞数量
# 功能: 对传入的细胞数以及参数进行处理并返回最终抽取细胞数量
setNum = function(number, cellsPercent, minCellsPerCluster) {
	selectNum = round(number * cellsPercent)
	if (selectNum >= minCellsPerCluster) {
		finalNum = selectNum
	} else {
		finalNum = min(number, minCellsPerCluster)
	}
	return(finalNum)
}

# Function to randomly sample cells
randomSample = function(cluster, size = NULL, cellsPercent, minCellsPerCluster) {
	samples = c()
	id = unique(cluster)
	for (i in id) {
		clusterInfo = cluster[cluster == i]
		if (is.null(size)) {
			finalNum = setNum(length(clusterInfo), cellsPercent, minCellsPerCluster)
		} else {
			finalNum = size
		}
		insample = sample(names(clusterInfo), finalNum, replace = FALSE)
		samples = c(samples, insample)
	}
	return(samples)
}

# 函数名: ImportSeurat4CDS
# 参数: otherCDS - Seurat对象
# 返回值: monocle_cds - 转换后的Monocle对象
# 功能: 将Seurat对象转换为Monocle对象，以便在Monocle中进行进一步的分析和可视化。
ImportSeurat4CDS = function(otherCDS) {
	data = GetAssayData(otherCDS, slot = "counts", assay = "RNA")
	# data = GetAssayData(otherCDS, slot = "counts", assay = "RNA3") #V5版本替换
	if (class(data) == "data.frame") {
		data = as(as.matrix(data), "sparseMatrix")
	}
	pd = tryCatch(
		{
			pd = new("AnnotatedDataFrame", data = otherCDS@meta.data)
			pd
		}, error = function(e) {
			pData = data.frame(cell_id = colnames(data), row.names = colnames(data))
			pd = new("AnnotatedDataFrame", data = pData)
			message("This Seurat object doesn't provide any meta data")
			pd
		}
	)
	if (length(setdiff(colnames(data), rownames(pd))) > 0) {
		data = data[, rownames(pd)]
	}
	fd = data.frame(gene_short_name = row.names(data), row.names = row.names(data))
	fd = new("AnnotatedDataFrame", data = fd)
	# 这里使用VGAM包的函数定义传入数据的分布类型，其中negbinomial.size()以及negbinomial()适用于UMI数据前者速度更快
	# 通常情况下将会运行第一种，FPKM值适用tobit
	# if (all(data == floor(data))) {
	expressionFamily <- negbinomial.size()
	# }else if(any(data < 0)){
	#   expressionFamily <- uninormal()
	# }else{
	#   expressionFamily <- tobit()
	# }
	# 创建monocle使用object，expressionFamily表达相应的VGAM函数
	monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
																lowerDetectionLimit = 0, expressionFamily = expressionFamily)
	if ("Monocle" %in% names(otherCDS@misc)) {
		otherCDS@misc$Monocle@auxClusteringData$seurat = NULL
		otherCDS@misc$Monocle@auxClusteringData$scran = NULL
		monocle_cds = otherCDS@misc$Monocle
		mist_list = otherCDS
	} else {
		mist_list = otherCDS
	}
	monocle_cds@auxClusteringData$seurat <- mist_list
	return(monocle_cds)
}

# Function to set ggplot theme
set_theme <- function() {
	theme(
		legend.position = "right",
		axis.title.x = element_text(size = 25),
		axis.title.y = element_text(size = 25),
		axis.text.x = element_text(size = 15),
		axis.text.y = element_text(size = 15)
	)
}

# Function to save plots
save_plots <- function(gg, base_name, outdir, ncol, nrow=1) {
	ggsave(paste0(outdir,"/",base_name, ".png"), plot = gg, width = 4 * ncol, height = 7 * nrow, limitsize = FALSE)
	ggsave(paste0(outdir,"/",base_name, ".pdf"), plot = gg, width = 4 * ncol, height = 7 * nrow, limitsize = FALSE)
}

# Function to calculate facet rows
calculate_facet_rows <- function(sample_vector) {
	unique_count <- length(unique(sample_vector))
	rows <- floor(sqrt(unique_count))
	return(rows)
}
calculate_facet_cols <- function(sample_vector, rows = NULL) {
  unique_count <- length(unique(sample_vector))
  if (is.null(rows)) {
    rows <- calculate_facet_rows(sample_vector)
  }
  cols <- ceiling(unique_count / rows)
  return(cols)
}

# 函数名: PsePlot
# 参数: CDS - monocle_cds对象
# 返回值: 画基础图
# 功能: 画各种图
PsePlot <- function(seuCDS, outdir) {
	color_groups <- list(State = "Pseudotime.State", Cluster = "Pseudotime.Cluster", Pseudotime = "Pseudotime.Pseudotime")
	for (color_var in names(color_groups)) {
		gg <- plot_cell_trajectory(seuCDS, color_by = color_var) + set_theme()
		save_plots(gg, color_groups[[color_var]], outdir, 2)
	}
	gg <- plot_cell_trajectory(seuCDS, color_by = "Cluster") +
		facet_wrap(~ Cluster, nrow = calculate_facet_rows(seuCDS$Cluster)) +
		set_theme()
	save_plots(gg, "Pseudotime.ClusterEach", outdir, calculate_facet_cols(seuCDS$Cluster))
	gg <- plot_cell_trajectory(seuCDS, color_by = "State") +
		facet_wrap(~ State, nrow = calculate_facet_rows(seuCDS$State)) +
		set_theme()
	save_plots(gg, "Pseudotime.StateEach", outdir, calculate_facet_cols(seuCDS$State))
	gg <- plot_complex_cell_trajectory(seuCDS, 
		color_by = 'Cluster', 
		show_branch_points = TRUE, 
		cell_size = 0.5, 
		cell_link_size = 0.3) +
		set_theme()
	save_plots(gg, "Pseudotime.tree", outdir,2)
}

# 函数名: PsePlotSample
# 参数: CDS - monocle_cds对象
# 返回值: 画样本图
# 功能: 画样本图
PsePlotSample <- function(seuCDS, outdir) {
	if ("sample" %in% colnames(seuCDS@phenoData@data) | "orig.ident" %in% colnames(seuCDS@phenoData@data)) {
		if ("orig.ident" %in% colnames(seuCDS@phenoData@data)) {
			seuCDS$Sample <- as.factor(seuCDS$orig.ident)
		} else {
			seuCDS$Sample <- as.factor(seuCDS$sample)
		}
		if (!is.null(seuCDS$Sample)) {
			gg <- plot_cell_trajectory(seuCDS, color_by = "Sample") + set_theme()
			save_plots(gg, "Pseudotime.Sample", outdir, 2)
		}
	}
	if ("group" %in% colnames(seuCDS@phenoData@data)) {
		gg <- plot_cell_trajectory(seuCDS, color_by = "group") + set_theme()
		save_plots(gg, "Pseudotime.group", outdir, 2)
	}
	plot_and_save_facet_wrap <- function(color_by, facet_var, base_name) {
		if (facet_var %in% colnames(seuCDS@phenoData@data)) {
			gg <- plot_cell_trajectory(seuCDS, color_by = color_by) +
				facet_wrap(as.formula(paste("~", facet_var)), nrow = calculate_facet_rows(seuCDS[[facet_var]])) +
				set_theme()
			save_plots(gg, paste0("Pseudotime.", base_name), outdir, calculate_facet_cols(seuCDS[[facet_var]]))
			
			gg_complex <- plot_complex_cell_trajectory(seuCDS, color_by = color_by, show_branch_points = TRUE, cell_size = 0.5, cell_link_size = 0.3) +
				facet_wrap(as.formula(paste("~", facet_var)), nrow = calculate_facet_rows(seuCDS[[facet_var]])) +
				scale_size(range = c(0.2, 0.2)) +
				theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
				set_theme()
			save_plots(gg_complex, paste0(base_name, ".tree"), outdir, calculate_facet_cols(seuCDS[[facet_var]]))
		}
	}
	plot_and_save_facet_wrap("Cluster", "Sample", "SampleEach")
	plot_and_save_facet_wrap("Cluster", "group", "groupEach")
}

# Function to generate density plots
generate_density_plots <- function(data, fill_var, facet_var, prefix,outdir) {
		nrow = floor(sqrt(length(unique(data[[facet_var]]))))
		ncol = ceiling(length(unique(data[[facet_var]])) / nrow)
		gg = ggplot(data, aes(x = Pseudotime, fill = !!rlang::sym(fill_var))) +
			geom_density(alpha = 0.5) +
			theme_cowplot() +
			facet_wrap(as.formula(paste("~", facet_var)), nrow = nrow)
		ggsave(paste0(outdir,"/",prefix, "Each.png"), gg, width = 4 * ncol, height = 2.5 * nrow, limitsize = FALSE, bg = "white")
		ggsave(paste0(outdir,"/",prefix, "Each.pdf"), gg, width = 4 * ncol, height = 2.5 * nrow, limitsize = FALSE, bg = "white")
		gg = ggplot(data, aes(x = Pseudotime, fill = !!rlang::sym(fill_var))) +
			geom_density(alpha = 0.5) +
			theme_cowplot()
		ggsave(paste0(outdir,"/",prefix, ".png"), gg, width = 8, height = 5, limitsize = FALSE, bg = "white")
		ggsave(paste0(outdir,"/",prefix, ".pdf"), gg, width = 8, height = 5, limitsize = FALSE, bg = "white")
}

# Function to plot pseudotime density
PseDenPlot = function(forplot,outdir) {
	if ("Sample" %in% colnames(forplot)) {
		if (length(unique(forplot$Sample)) > 1) {
			generate_density_plots(forplot, "Sample", "Sample", "PseudotimeDensity.Sample",outdir)
		} else {
			generate_density_plots(forplot, "Cluster", "Cluster", "PseudotimeDensity.Cluster",outdir)
		}
	}
	if ("group" %in% colnames(forplot)) {
		generate_density_plots(forplot, "group", "group", "PseudotimeDensity.group",outdir)
	}
	gg = ggplot(forplot, aes(x = Pseudotime)) +
		geom_density(alpha = 0.5, fill = "blue") +
		theme_cowplot()
	ggsave(paste0(outdir,"/PseudotimeDensity.png"), gg, width = 8, height = 5, limitsize = FALSE, bg = "white")
	ggsave(paste0(outdir,"/PseudotimeDensity.pdf"), gg, width = 8, height = 5, limitsize = FALSE, bg = "white")
}

# Main function to run the analysis
mainfun = function(seuset, outdir, geneFile, lambda, RootState, Reverse, celllist, cluster, RandomCells,
									cellsPercent, minCellsPerCluster) {
	
	seuset = readRDS(seuset)
	# seuset[["RNA"]] <- as(object = seuset[["RNA"]], Class = "Assay") #V5版本替换
	# DefaultAssay(seuset) <- 'RNA'
	seuset$cell_type = Idents(seuset)
	# seuset <- subset(seuset,subset=sample %in% c('HLH01-LAHS-B','HLH02-LAHS-T','HLH03-LAHS-T','HLH04-LAHS'))
	seuset$group<-factor(seuset$group)
	seuset$sample<-factor(seuset$sample)
	# seuset@active.ident = seuset$cell_type
	if ("cell_type" %in% colnames(seuset@meta.data)) {
		seuset@active.ident = seuset$cell_type
	} else {
		stop("Error: 'cell_type' column not found in seuset@meta.data")
	}
	seuset = Psepre(seuset, outdir, celllist, cluster, RandomCells, cellsPercent, minCellsPerCluster)
	print("Psepre done")
	if (geneFile!='F') {
		genelist = unique(read.table(geneFile, sep = "\t", header = T)[, 1])
	} else {
		# seuset = FindVariableFeatures(seuset, selection.method = 'vst', nfeatures = 2000)
		# genelist = VariableFeatures(seuset)
		# write.table(genelist,file = paste0(outdir,"/","VariableFeatures.txt"), sep = "\t", row.names = F, quote = F)
	}
	# 获取lambda参数
	lambda = lambda*ncol(seuset)
	# 创建monocle对象,标记用于聚类的基因，标化数据，并进行降维,排序，计算拟时间
	# estimateSizeFactors 用于裱花细胞间的差异 lamda越大branch数越小
	# reduceDimension max_components降维空间维数
	# max_components is the number of dimensions you want returned.
	# num_dim is the number of dimension for the PCA preprocessing.
	seuCDS = ImportSeurat4CDS(seuset) %>% estimateSizeFactors() %>% estimateDispersions()
	disp_table <- dispersionTable(seuCDS)
	genelist <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
	write.table(genelist,file = paste0(outdir,"/","dispersionTable.txt"), sep = "\t", row.names = F, quote = F)
	# deg_table <- FindAllMarkers(seuset)
	# genelist <- subset(deg_table,p_val_adj < 0.05)$gene
	# write.table(genelist,file = paste0(outdir,"/","FindAllMarkers.txt"), sep = "\t", row.names = F, quote = F)
	seuCDS = setOrderingFilter(seuCDS,genelist) %>% 
		reduceDimension(reduction_method = "DDRTree", lambda = lambda, num_dim = 6
		, residualModelFormulaStr = "~sample"    #根据样本去批次，分支乱有批次的时候使用
		# , auto_param_selection = F  , ncenter = 5,  #细胞数量较低的时候需要有这行代码 Error in kmeans(t(Z), K, centers = centers) : initial centers are not distinct,
		) %>%
		orderCells()
	print("orderCells done")
	ggdim = plot_pc_variance_explained(seuCDS, 
                           return_all = F,
                           max_components = 30)
	ggsave(paste0(outdir,"/Pseudotime_dim.png"), ggdim, width = 8, height = 5, limitsize = FALSE, bg = "white")
	ggsave(paste0(outdir,"/Pseudotime_dim.pdf"), ggdim, width = 8, height = 5, limitsize = FALSE, bg = "white")
	if (RootState!='F' | Reverse!='F') {
		seuCDS = orderCells(seuCDS, root_state = RootState, reverse = Reverse)
	}
	table = table(seuCDS$Cluster, seuCDS$State) %>% as.data.frame() %>% 
		dcast(Var1 ~ Var2)
	colnames(table)[1] = "Cluster"
	table %>% write.table(paste0(outdir,"/Pseudotime.Statistics.txt"), sep = "\t", quote = F)
	seuCDS@phenoData@data %>% data.frame(CellName = row.names(seuCDS@phenoData@data)) %>%
		write.table(paste0(outdir,"/Pseudotime.Summary.txt"), sep = "\t", quote = F, row.names = F)
	PsePlot(seuCDS, outdir)
	if (length(unique(seuCDS$sample)) > 1 | length(unique(seuCDS$orig.ident)) > 1) {
		PsePlotSample(seuCDS, outdir)
	}
	if ("sample" %in% colnames(seuCDS@phenoData@data)) {
		forplot = data.frame(Sample = seuCDS$sample, Cluster = seuCDS$Cluster, Pseudotime = seuCDS$Pseudotime)
		PseDenPlot(forplot, outdir)
	}
	if ("orig.ident" %in% colnames(seuCDS@phenoData@data)) {
		forplot = data.frame(Sample = seuCDS$orig.ident, Cluster = seuCDS$Cluster, Pseudotime = seuCDS$Pseudotime)
		PseDenPlot(forplot, outdir)
	}
	if ("group" %in% colnames(seuCDS@phenoData@data)) {
		forplot = data.frame(group = seuCDS$group, Cluster = seuCDS$Cluster, Pseudotime = seuCDS$Pseudotime)
		PseDenPlot(forplot, outdir)
	}
	
	# 寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图。
	Time_diff <- differentialGeneTest(seuCDS[genelist,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
	Time_diff <-Time_diff[, c(5, 3, 4)] # 把gene放前面，也可以不改
	Time_genes <- Time_diff %>% filter(pval < 0.05) %>% arrange(qval)
	write.csv(Time_genes, paste0(outdir,"/Pseudotime_diff_all.csv"), row.names = F)
	Time_genes <- Time_genes %>% head(10) %>% pull(gene_short_name) %>% as.character() # 提取Top10gene
# 直接查看一些基因随细胞状态等的表达变化,选择前10个top基因并将其对象取出
	seuCDS_subset <- seuCDS[Time_genes,]
# 可视化：以state/cluster/pseudotime进行
	gg=plot_genes_in_pseudotime(seuCDS_subset, color_by = "group", ncol = 5)
	ggsave(paste0(outdir,"/pseudotime_genes_group.png"),limitsize = FALSE,width= 20,height= 10)
	ggsave(paste0(outdir,"/pseudotime_genes_group.pdf"),limitsize = FALSE,width= 20,height= 10)

	saveRDS(seuCDS, file = paste0(outdir,"/Pseudotime.monocle.rds"))
}

# Run the main function
mainfun(rdsFile, outdir = outdir, geneFile = geneFile, lambda = lambda, RootState = RootState,
				Reverse = Reverse, celllist = celllist, cluster = cluster,
				RandomCells = RandomCells, cellsPercent = cellsPercent,
				minCellsPerCluster = minCellsPerCluster)
