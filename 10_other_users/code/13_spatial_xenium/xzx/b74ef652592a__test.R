library(jsonlite)
library(languageserver)
library(httpgd)

setwd("/mnt/NexData_file/home/xzx/spatial_project/XE2512297838-01/")

dir.create("/mnt/NexData_file/home/xzx/spatial_project/XE2512297838-01/test")
setwd("/mnt/NexData_file/home/xzx/spatial_project/XE2512297838-01/test")

library(ggplot2)
library(tidyverse)
library(Seurat)
library(spacexr)
library(httpgd)
hgd()


install.packages("tidyverse", repos ="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

install.packages("devtools")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
hgd()

# 创建测试数据
df <- data.frame(
  x = 1:10,
  y = (1:10) + rnorm(10)
)

# 绘制简单散点图
ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "red") +
  ggtitle("Test Plot: ggplot2") +
  theme_minimal()
hgd_browse()
