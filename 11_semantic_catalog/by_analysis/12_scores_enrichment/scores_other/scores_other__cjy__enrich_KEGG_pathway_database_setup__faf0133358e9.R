log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type=c("message"))

user_lib <- snakemake@params[["kegg_r_lib"]]
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

if (snakemake@params[["species"]] == "homo_sapiens") {
  kegg_database <- snakemake@params[["kegg_hsa_pathway_database"]]
} else if (snakemake@params[["species"]] == "mus_musculus") {
  kegg_database <- snakemake@params[["kegg_mmu_pathway_database"]]
}  

if (!file.exists(file.path(user_lib, "KEGG.db", "DESCRIPTION"))) {
  install.packages(kegg_database, lib=user_lib, repos=NULL, type = "source", quiet=TRUE)
}
library(KEGG.db, lib.loc=user_lib)

system(paste("touch", snakemake@output[[1]]), intern=TRUE)
