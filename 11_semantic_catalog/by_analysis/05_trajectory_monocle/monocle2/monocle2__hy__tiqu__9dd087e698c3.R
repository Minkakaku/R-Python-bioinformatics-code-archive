# 加载必要的包
library(Seurat)
library(Matrix)

# 函数：读取RDS文件并提取metadata和counts
read_rds_info <- function(rds_path, layer = "counts") {
  # 参数说明：
  # rds_path: RDS文件路径
  # layer: 要提取的counts层，默认为"counts"，也可以是"data"等
  
  # 检查文件是否存在
  if (!file.exists(rds_path)) {
    stop(paste("文件不存在:", rds_path))
  }
  
  # 读取RDS文件
  cat("正在读取RDS文件...\n")
  seu_obj <- readRDS(rds_path)
  cat("文件读取完成。\n")
  
  # 检查对象类型
  if (inherits(seu_obj, "Seurat")) {
    cat("检测到Seurat对象。\n")
    
    # 提取metadata
    metadata <- seu_obj@meta.data
    cat(paste0("metadata维度: ", nrow(metadata), "行 x ", ncol(metadata), "列\n"))
    
    # 提取counts
    cat(paste0("尝试提取", layer, "矩阵...\n"))
    
    # 检查是否有RNA assay
    if ("RNA" %in% names(seu_obj@assays)) {
      # 尝试获取counts
      if (layer == "counts" && !is.null(seu_obj@assays$RNA@counts)) {
        counts_matrix <- seu_obj@assays$RNA@counts
        cat(paste0("成功提取", layer, "矩阵: ", nrow(counts_matrix), 
                  "个基因 x ", ncol(counts_matrix), "个细胞\n"))
      } 
      # 尝试获取data
      else if (layer == "data" && !is.null(seu_obj@assays$RNA@data)) {
        counts_matrix <- seu_obj@assays$RNA@data
        cat(paste0("成功提取", layer, "矩阵: ", nrow(counts_matrix), 
                  "个基因 x ", ncol(counts_matrix), "个细胞\n"))
      } 
      # 尝试获取其他层
      else if (exists(paste0(layer, "s"), seu_obj@assays$RNA)) {
        counts_matrix <- seu_obj@assays$RNA[[paste0(layer, "s")]]
        cat(paste0("成功提取", layer, "矩阵: ", nrow(counts_matrix), 
                  "个基因 x ", ncol(counts_matrix), "个细胞\n"))
      } else {
        cat("警告: 未找到指定的", layer, "层，尝试其他方法...\n")
        
        # 列出可用的层
        cat("可用的层:\n")
        print(names(seu_obj@assays$RNA))
        
        # 尝试使用GetAssayData函数
        counts_matrix <- GetAssayData(seu_obj, layer = layer)
        if (length(dim(counts_matrix)) > 0) {
          cat(paste0("通过GetAssayData提取", layer, "矩阵: ", 
                    nrow(counts_matrix), "个基因 x ", ncol(counts_matrix), "个细胞\n"))
        } else {
          stop(paste("无法提取", layer, "矩阵"))
        }
      }
    } else {
      cat("警告: 未找到RNA assay，尝试其他assay...\n")
      
      # 列出可用的assay
      cat("可用的assay:\n")
      print(names(seu_obj@assays))
      
      # 使用第一个可用的assay
      assay_name <- names(seu_obj@assays)[1]
      cat(paste0("使用assay: ", assay_name, "\n"))
      
      counts_matrix <- GetAssayData(seu_obj, assay = assay_name, layer = layer)
      cat(paste0("提取", assay_name, " assay的", layer, "矩阵: ", 
                nrow(counts_matrix), "个基因 x ", ncol(counts_matrix), "个细胞\n"))
    }
    
  } else if (inherits(seu_obj, "SingleCellExperiment")) {
    cat("检测到SingleCellExperiment对象。\n")
    
    # 提取metadata
    metadata <- colData(seu_obj)
    cat(paste0("metadata维度: ", nrow(metadata), "行 x ", ncol(metadata), "列\n"))
    
    # 提取counts
    cat(paste0("尝试提取", layer, "矩阵...\n"))
    
    if (layer == "counts") {
      counts_matrix <- counts(seu_obj)
    } else {
      counts_matrix <- assay(seu_obj, layer)
    }
    
    if (!is.null(counts_matrix)) {
      cat(paste0("成功提取", layer, "矩阵: ", nrow(counts_matrix), 
                "个基因 x ", ncol(counts_matrix), "个细胞\n"))
    } else {
      stop(paste("无法提取", layer, "矩阵"))
    }
    
  } else {
    cat("未知的对象类型，尝试直接提取...\n")
    
    # 尝试提取metadata
    if ("meta.data" %in% slotNames(seu_obj)) {
      metadata <- seu_obj@meta.data
      cat(paste0("metadata维度: ", nrow(metadata), "行 x ", ncol(metadata), "列\n"))
    } else {
      metadata <- NULL
      cat("警告: 无法提取metadata\n")
    }
    
    # 尝试提取counts
    if ("counts" %in% names(seu_obj)) {
      counts_matrix <- seu_obj$counts
      cat(paste0("成功提取counts矩阵: ", nrow(counts_matrix), 
                "个基因 x ", ncol(counts_matrix), "个细胞\n"))
    } else {
      cat("尝试查找counts矩阵...\n")
      if (exists("counts", where = seu_obj)) {
        counts_matrix <- seu_obj$counts
        cat(paste0("成功提取counts矩阵: ", nrow(counts_matrix), 
                  "个基因 x ", ncol(counts_matrix), "个细胞\n"))
      } else {
        counts_matrix <- NULL
        cat("警告: 无法提取counts矩阵\n")
      }
    }
  }
  
  # 返回结果
  result <- list(
    metadata = metadata,
    counts = counts_matrix,
    object_type = class(seu_obj)[1]
  )
  
  return(result)
}

# 函数：保存提取的数据
save_extracted_data <- function(extracted_data, output_dir = ".") {
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存metadata
  if (!is.null(extracted_data$metadata)) {
    metadata_path <- file.path(output_dir, "metadata.csv")
    write.csv(extracted_data$metadata, metadata_path, row.names = TRUE)
    cat(paste0("metadata已保存到: ", metadata_path, "\n"))
  }
  
  # 保存counts矩阵
  if (!is.null(extracted_data$counts)) {
    counts_path <- file.path(output_dir, "counts_matrix.rds")
    saveRDS(extracted_data$counts, counts_path)
    cat(paste0("counts矩阵已保存到: ", counts_path, "\n"))
    
    # 如果矩阵较小，也可以保存为CSV
    if (nrow(extracted_data$counts) * ncol(extracted_data$counts) < 1e7) {
      # 转换为普通矩阵
      if (inherits(extracted_data$counts, "sparseMatrix")) {
        counts_dense <- as.matrix(extracted_data$counts)
      } else {
        counts_dense <- extracted_data$counts
      }
      
      counts_csv_path <- file.path(output_dir, "counts_matrix.csv")
      write.csv(counts_dense, counts_csv_path, row.names = TRUE)
      cat(paste0("counts矩阵(密集格式)已保存到: ", counts_csv_path, "\n"))
    } else {
      cat("矩阵太大，跳过保存CSV格式。\n")
    }
  }
  
  # 保存信息摘要
  summary_info <- list(
    object_type = extracted_data$object_type,
    metadata_dim = if(!is.null(extracted_data$metadata)) dim(extracted_data$metadata) else NULL,
    counts_dim = if(!is.null(extracted_data$counts)) dim(extracted_data$counts) else NULL
  )
  
  summary_path <- file.path(output_dir, "summary_info.rds")
  saveRDS(summary_info, summary_path)
  cat(paste0("信息摘要已保存到: ", summary_path, "\n"))
  
  # 打印摘要
  cat("\n===== 数据摘要 =====\n")
  cat(paste0("对象类型: ", summary_info$object_type, "\n"))
  if(!is.null(summary_info$metadata_dim)) {
    cat(paste0("metadata维度: ", summary_info$metadata_dim[1], 
              "行 x ", summary_info$metadata_dim[2], "列\n"))
  }
  if(!is.null(summary_info$counts_dim)) {
    cat(paste0("counts矩阵维度: ", summary_info$counts_dim[1], 
              "个基因 x ", summary_info$counts_dim[2], "个细胞\n"))
  }
  cat("===== 结束 =====\n")
}

# 主程序
main <- function() {
  # 设置参数
  rds_path <- "xenium_rctd_2.rds"  # 替换为您的RDS文件路径
  output_dir <- "extracted_data"  # 输出目录
  layer_to_extract <- "counts"  # 要提取的层，可以是"counts"或"data"等
  
  # 提取数据
  extracted_data <- read_rds_info(rds_path, layer = layer_to_extract)
  
  # 保存数据
  save_extracted_data(extracted_data, output_dir = output_dir)
  
  # 返回提取的数据
  return(extracted_data)
}

# 执行主程序
if (interactive()) {
  # 如果是交互式环境，可以注释掉下面这行，手动运行
  # result <- main()
  
  # 或者使用这个快捷函数
  quick_read <- function(file_path) {
    cat("快捷读取函数示例:\n")
    cat("1. 设置文件路径: rds_path <- \"your_file.rds\"\n")
    cat("2. 调用: data_info <- read_rds_info(rds_path)\n")
    cat("3. 查看metadata: head(data_info$metadata)\n")
    cat("4. 查看counts维度: dim(data_info$counts)\n")
  }
  
  cat("脚本已加载。请先设置正确的RDS文件路径，然后调用read_rds_info()函数。\n")
}
