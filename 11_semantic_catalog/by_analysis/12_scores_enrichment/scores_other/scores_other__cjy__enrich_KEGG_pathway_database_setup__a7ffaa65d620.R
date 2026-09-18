log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type=c("message"))

user_lib <- snakemake@params[["kegg_r_lib"]]
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

if (!file.exists(file.path(user_lib, "KEGG.db", "DESCRIPTION"))) {
  install.packages(snakemake@params[["kegg_rat_pathway_database"]], lib=user_lib, repos=NULL, type = "source", quiet=TRUE)
}
library(KEGG.db, lib.loc=user_lib)

if (snakemake@params[["species"]] == "rattus_norvegicus") {
  # install.packages(snakemake@params[["kegg_hsa_pathway_database"]], repos=NULL, type = "source", quiet=TRUE)
  library("org.Rn.eg.db")
}


system(paste("touch", snakemake@output[[1]]), intern=TRUE)
