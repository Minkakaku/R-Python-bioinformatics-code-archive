
suppressPackageStartupMessages({
  library(utils)
})
rda <- "C:/Users/Administrator/Desktop/生信标准工具/Database/omnipath/CellChatDB.mouse.rda"
out <- "C:/Users/Administrator/Desktop/生信标准工具/Database/omnipath/CellChatDB_mouse_interaction.tsv"

# 1) load 进当前环境
e <- new.env(parent = emptyenv())
load(rda, envir = e)
objs <- ls(e)
cat("Objects loaded:", paste(objs, collapse=", "), "\n")

# 2) 找到包含 interaction 的对象（通常是 CellChatDB.mouse）
target <- NULL
for (nm in objs) {
  obj <- get(nm, envir = e)
  if (is.list(obj) && ("interaction" %in% names(obj))) {
    target <- nm
    break
  }
}
if (is.null(target)) {
  stop("No object with $interaction found in the RDA.")
}
cat("Using object:", target, "\n")

db <- get(target, envir = e)
inter <- db$interaction

# 3) 导出 TSV
write.table(inter, file = out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
cat("Saved:", out, "\n")
cat("interaction dim:", paste(dim(inter), collapse=" x "), "\n")
cat("interaction cols (first 40):", paste(head(colnames(inter), 40), collapse=", "), "\n")
