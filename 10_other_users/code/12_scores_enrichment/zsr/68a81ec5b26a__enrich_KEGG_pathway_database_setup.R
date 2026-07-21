log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type=c("message"))

if (snakemake@params[["species"]] == "homo_sapiens") {
  install.packages(snakemake@params[["kegg_hsa_pathway_database"]], repos=NULL, type = "source", quiet=TRUE)
  library(KEGG.db)
} else if (snakemake@params[["species"]] == "mus_musculus") {
  install.packages(snakemake@params[["kegg_mmu_pathway_database"]], repos=NULL, type = "source", quiet=TRUE)
  library(KEGG.db)
}  

system(paste("touch", snakemake@output[[1]]), intern=TRUE)