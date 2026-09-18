log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type=c("message"))

install.packages(snakemake@params[["kegg_rat_pathway_database"]], repos=NULL, type = "source", quiet=TRUE)
library(KEGG.db)

system(paste("touch", snakemake@output[[1]]), intern=TRUE)