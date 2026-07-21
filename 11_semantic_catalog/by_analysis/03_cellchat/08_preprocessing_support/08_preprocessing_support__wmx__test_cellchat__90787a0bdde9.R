
library(CellChat,lib.loc = "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/lib/R/library/")

load("/mnt/NexData_file/home/yrc/AE2408071595-01/4.15/cellchat/HLH_NC/tmp_object.RData")

pathways <- c("BAFF", "TWEAK", "TRAIL", "APRIL", "LIGHT", "CD30")

lapply(object.list, function(x) {
  grep(paste(pathways, collapse = "|"), x@netP$pathways, value = TRUE)
})


# cellchat_merge <- mergeCellChat(
#   object.list,
#   add.names = names(object.list)
# )
# netVisual_bubble(
#   cellchat_merge,
#   comparison = c(1, 2),
#   signaling = c("BAFF", "TWEAK","TRAIL", "APRIL","LIGHT", "CD30"),
#   angle.x = 45,
#   remove.isolate = FALSE
# )
# 

load("/mnt/NexData_file/home/yrc/AE2408071595-01/4.15/cellchat/HLH_iMCD-TAFRO/tmp_object.RData")

pathways <- c("BAFF", "TWEAK", "TRAIL", "APRIL", "LIGHT", "CD30")

lapply(object.list, function(x) {
  grep(paste(pathways, collapse = "|"), x@netP$pathways, value = TRUE)
})
# $HLH
# character(0)
# 
# $`iMCD-TAFRO`
# [1] "TRAIL"

