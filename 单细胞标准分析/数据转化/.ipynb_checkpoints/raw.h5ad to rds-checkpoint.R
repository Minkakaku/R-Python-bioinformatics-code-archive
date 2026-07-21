
# 清理所有 reticulate 的 Python 缓存变量
Sys.unsetenv("RETICULATE_PYTHON")
Sys.unsetenv("RETICULATE_CONDA")
Sys.unsetenv("RETICULATE_PYTHON_FALLBACK")
unlink("C:/Users/Administrator/AppData/Local/R/cache/R/reticulate", recursive = TRUE, force = TRUE)
#############此处务必保存环境并重启R#######################

Sys.setenv(RETICULATE_CONDA = "D:/Users/Administrator/anaconda3/Scripts/conda.exe")

library(reticulate)

use_condaenv(
  condaenv = "scanpy_env",
  conda = Sys.getenv("RETICULATE_CONDA"),
  required = TRUE
)
#D:\Users\Administrator\anaconda3\envs\scanpy_env\python.exe

sceasy::convertFormat("adata.h5ad",from="anndata",to="seurat",
                      outFile = "result.rds")
