# 整合T和Epith的rds，以便进行后续的细胞通讯分析
library(Seurat)
library(dplyr)
setwd("/mnt//NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/")

# 1. 加载数据
t_data <- readRDS("T_score/T_score.rds")
epith_data <- readRDS("Epith_score/Epith_score.rds")

# 2. 确保元数据有您需要的列
# 假设您原来的数据里已经有了 group 和 cell_type，如果列名不同，需要手动统一
# 例如: t_data$cell_type <- t_data$orig.ident (如果是按样本存的)

# 3. 合并对象
# 使用 merge，会自动处理表达矩阵的交集/并集
combined_obj <- merge(t_data, y = epith_data, 
                      add.cell.ids = c("Tcell", "Epith"), 
                      project = "4T1_CellChat")

# 4. 检查 metadata 是否完好
# 确保所有细胞都有这两个关键列
table(combined_obj$group)
table(combined_obj$cell_type)

dir.create('cellchat/')
saveRDS(combined_obj,'cellchat/T_epith_combined.rds')
