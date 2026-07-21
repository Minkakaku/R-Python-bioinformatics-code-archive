library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)
install.packages("gridExtra")


load("A-vs-B_GSEA_GO_enrich_results.Rdata")
load("A-vs-B_GSEA_KEGG_enrich_results.Rdata")
load("A-vs-B_GSEA_Reactome_enrich_results.Rdata")

# pathways <- c(
#   "Neutrophil degranulation",
#   "Complement and coagulation cascades",
#   "Cytokine-cytokine receptor interaction",
#   "IL-17 signaling pathway"
# )

pathways <- c(
  "intermediate filament cytoskeleton",
  "Complement and coagulation cascades",
  "Cytokine-cytokine receptor interaction",
  "Cell Cycle, Mitotic",
  "Neutrophil degranulation"
)



plot_and_save_gsea <- function(gsea_obj, pathway_name, outdir, db){
  
  idx <- which(gsea_obj@result$Description == pathway_name)
  
  if(length(idx) == 0){
    return(NULL)
  }
  
  p <- gseaplot2(
    gsea_obj,
    geneSetID = idx,
    title = pathway_name,
    color = "red",
    pvalue_table = TRUE
  )
  
  # 文件名安全处理（防止空格/特殊符号）
  safe_name <- gsub("[/:*?\"<>| ]", "_", pathway_name)
  
  pdf(file = file.path(outdir, paste0(db, "_", safe_name, ".pdf")),
      width = 10, height = 7)
  
  print(p)
  dev.off()
  
  png(file.path(outdir, paste0(db, "_", safe_name, ".png")),
      width = 2000, height = 1400, res = 200)
  
  print(p)
  dev.off()
  
  return(p)
}
gsea_list <- list(
  GO = GO_gseresult,
  KEGG = KEGG_gseresult,
  Reactome = Reactome_gseresult
)


outdir <- "./GSEA_plots"
dir.create(outdir, showWarnings = FALSE)

# plot_list <- list()
# for(db in names(gsea_list)){
#   
#   obj <- gsea_list[[db]]
#   
#   for(pw in pathways){
#     
#     p <- plot_gsea(obj, pw)
#     
#     if(!is.null(p)){
#       plot_list[[paste(db, pw, sep = "_")]] <- p
#     }
#   }
# }
# 

plot_store <- list()

for(db in names(gsea_list)){
  
  obj <- gsea_list[[db]]
  
  for(pw in pathways){
    
    p <- plot_and_save_gsea(
      gsea_obj = obj,
      pathway_name = pw,
      outdir = outdir,
      db = db
    )
    
    if(!is.null(p)){
      plot_store[[paste(db, pw, sep = "_")]] <- p
    }
  }
}

for(p in plot_list){
  print(p)
}