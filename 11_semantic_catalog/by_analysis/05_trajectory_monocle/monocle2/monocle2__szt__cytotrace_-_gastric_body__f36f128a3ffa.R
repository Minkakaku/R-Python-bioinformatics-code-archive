library(reticulate)
use_python("/home/szt/miniforge3/envs/rpy/bin/python")

library(CytoTRACE)
source("cytotrace-fix bug.R")

out_dir <- "cytotrace-body/"
dir.create(out_dir)


sc_epi2 <- readRDS("monocle - v2 - gastric body - remove distal/sc_epi2.rds")
exp_dat <- sc_epi2 %>% JoinLayers() %>% LayerData(layer = "counts") %>% as.matrix()

cytotrace_results <- CytoTRACE(exp_dat, ncores = 1)

# saveRDS(cytotrace_results, file = file.path(out_dir, "cytotrace_results.rds"))

phenot <- sc_epi2$celltype
phenot <- as.character(phenot)
names(phenot) <- rownames(sc_epi2@meta.data)

# emb <- sc_epi2@reductions[["harmonyumap"]]@cell.embeddings
emb <- data.frame(
  row.names = colnames(sc_epi2),
  x = rep(1, length.out = ncol(sc_epi2)),
  y = rep(1, length.out = ncol(sc_epi2))
)

undebug(plotCytoTRACE)

plotCytoTRACE(cytotrace_results, phenotype = phenot, emb = emb, outputDir = out_dir)
plotCytoGenes(cytotrace_results, numOfGenes = 30, outputDir = out_dir)


# cds <- readRDS("monocle - v2 - gastric antrum - remove distal/cds2.rds")
p_state <- 
  data.frame(
    score = cytotrace_results$CytoTRACE,
    celltype = phenot
  )  %>%
  rownames_to_column("cell_id") %>%
  right_join(
    rownames_to_column(pData(cds)['State'], "cell_id"),
    by = 'cell_id'
  ) %>%
  ggplot(aes(x = State, y = score, fill = State)) +
  geom_boxplot() +
  theme_classic() +
  ggsci::scale_fill_igv()

ggsave(
  p_state,
  filename = file.path(out_dir, "cytotrace score of monocle state.pdf"),
  height = 6,
  width = 8
)





######################## cytotrace2 ######################

# load the CytoTRACE 2 package
library(CytoTRACE2) 


# running CytoTRACE 2 main function - cytotrace2 - with default parameters
cytotrace2_result <- cytotrace2(exp_dat)
# Error in names(cytotrace) <- unlist(sample_names) : 
#   'names' attribute [6376] must be the same length as the vector [7]


phenot <- sc_epi2$annotation
phenot <- as.character(phenot)
names(phenot) <- rownames(sc_epi2@meta.data)


# generate prediction and phenotype association plots with plotData function
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = phenot,
                  expression_data = exp_dat
)


