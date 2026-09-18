library(CellChat,lib.loc = "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/lib/R/library/")
load("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/cellchat/LARC_N/tmp_object.RData")


cellchat_merge <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

pdf(
  file.path(outdir, paste0(pathway_use, "_group_comparison_bubble.pdf")),
  width = 10,
  height = 6
)

print(
  netVisual_bubble(
    cellchat_merge,
    # signaling = pathway_use,
    sources.use = celltypes_use,
    targets.use = celltypes_use,
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()