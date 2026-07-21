# library(tradeSeq)

out_dir_tradeseq <- "test_tradeseq"
dir.create(out_dir_tradeseq)


data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

#0 parallel multisession  --------
# bpparms <- BiocParallel::SerialParam()
bpparms <- BiocParallel::MulticoreParam(worker = 6)

BiocParallel::register(bpparms, default = FALSE) # change backend list
BiocParallel::bpparam() # registed backend

sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE, parallel=TRUE, BPPARAM = bpparms)


#1 Fit negative binomial model  --------
tradeSeq::plotGeneCount(curve = crv, clusters = celltype, 
              title = "Colored by cell type")


# set.seed(5)
# icMat <- evaluateK(
#   counts = counts,
#   sds = crv,
#   k = 3:10,
#   nGenes = 200,
#   verbose = T
# )
### Downstream of any trajectory inference method using pseudotime and cell weights
set.seed(7)
pseudotime <- slingPseudotime(crv, na=FALSE)
cellWeights <- slingCurveWeights(crv)
icMat2 <- tradeSeq::evaluateK(
  counts = counts, 
  pseudotime = pseudotime,
  cellWeights = cellWeights,
  k=3:7, nGenes = 100,
  verbose = FALSE, plot = TRUE
)

tradeSeq::plot_evalutateK_results(icMat2)


#2 fit GAM   --------

# set.seed(6)
# sce <- fitGAM(counts = counts, sds = crv, nknots = 6, verbose = FALSE)


set.seed(7)
# pseudotime <- slingPseudotime(crv, na = FALSE)
# cellWeights <- slingCurveWeights(crv)
sce <- tradeSeq::fitGAM(
  counts = counts, pseudotime = pseudotime,
  cellWeights = cellWeights,
  nknots = 6, verbose = FALSE
)

# if add covariates 
# batch <- factor(rep(c("A", "B"), each = ncol(counts)/2))
# U <- model.matrix(~batch)
# sceBatch <- fitGAM(counts = counts,
#                    U = U,
#                    sds = crv, 
#                    nknots = 6, 
#                    verbose = FALSE)



# # more control in GAM
# control <- gam.control()
# control$maxit <- 1000 #set maximum number of iterations to 1K
# # pass to control argument of fitGAM as below:
# # 
# # gamList <- fitGAM(counts = counts,
# #                   pseudotime = slingPseudotime(crv, na = FALSE),
# #                   cellWeights = slingCurveWeights(crv),
# #                   control = control)


# 3 Within-lineage comparisons  --------
## association
# is testing the null hypothesis that all smoother coefficients are equal to each other.

assoRes <- tradeSeq::associationTest(sce)
head(assoRes)

## start vs end 
startRes <- tradeSeq::startVsEndTest(sce, lineages = FALSE)

oStart <- order(startRes$waldStat, decreasing = TRUE)
startRes %>% head
# waldStat df       pvalue logFClineage1 logFClineage2
# Acin1   4.31946  2 1.153563e-01     0.3641027    0.03131154
# Actb  199.17232  2 0.000000e+00     0.6915002    0.31875690
# Ak2    64.68311  2 8.992806e-15     0.5912243    1.07951355
# Alad  110.38859  2 0.000000e+00    -0.9043989    1.96150328
# Alas1 639.89546  2 0.000000e+00     3.6477817   -0.21519090
# Aldoa  75.92282  2 0.000000e+00     1.1115989    0.12368087
sigGeneStart <- names(sce)[oStart[3]]
tradeSeq::plotSmoothers(sce, counts, gene = sigGeneStart)


tradeSeq::plotGeneCount(crv, counts, gene = sigGeneStart)

# compare the expression for each gene at pseudotime values of 0.8 and 0.1
# customRes <- tradeSeq::startVsEndTest(sce, pseudotimeValues = c(0.1, 0.8))


# 4. Between-lineage comparisons ----------------
# 4.1  Discovering differentiated cell type markers
# one can assess all pairwise comparisons using the pairwise=TRUE argument.
endRes <- tradeSeq::diffEndTest(sce)

o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
tradeSeq::plotSmoothers(sce, counts, sigGene)
tradeSeq::plotGeneCount(crv, counts, gene = sigGene)

# 4.2 different expression patterns

patternRes <- tradeSeq::patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])


tradeSeq::plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][4])
tradeSeq::plotGeneCount(crv, counts, gene = rownames(patternRes)[oPat][4])


# 4.3 if you insterested in genes that have different expression patterns
# but no different expression at the differentiated cell type level. 

# library(ggplot2)

patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes <- patternRes[, c("Gene", "pattern")]

endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes <- endRes[, c("Gene", "end")]

compare <- merge(patternRes, endRes, by = "Gene", all = FALSE)
compare$transientScore <- 
  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2

ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "patternTest Wald Statistic (log scale)",
       y = "diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()

topTransient <- compare[which.max(compare$transientScore), "Gene"]
tradeSeq::plotSmoothers(sce, counts, gene = topTransient)

tradeSeq::plotGeneCount(crv, counts, gene = topTransient)

head(
  compare[order(compare$transientScore, decreasing = TRUE), "Gene"],
  n = 5
)



# 4.4 Early drivers of differentiation
cluster_label <- apply(slingClusterLabels(crv), 1, which.max) 
cluster_label <- colnames(slingClusterLabels(crv))[cluster_label]

tradeSeq::plotGeneCount(curve = crv, counts = counts,
              clusters = apply(slingClusterLabels(crv), 1, which.max),
              models = sce)

earlyDERes <- tradeSeq::earlyDETest(sce, knots = c(1, 2))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])

tradeSeq::plotSmoothers(sce, counts, gene = rownames(earlyDERes)[oEarly][2])
tradeSeq::plotGeneCount(crv, counts, gene = head(rownames(earlyDERes)[oEarly])[2])

# 5 set log2fc --------------
# testing against fold change threshold of 2
start2 <- startVsEndTest(sce, l2fc = log2(2))
# testing against fold change threshold of 1.5
pat2 <- patternTest(sce, l2fc = log2(1.5))


# 6
# 6.1 Extracting fitted values
yhat <- predictCells(models = sce, gene = "Irf8")

ysmooth <- predictSmooth(models = sce, gene = "Irf8", nPoints = 40)



# 6.2 
library(clusterExperiment)
nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                     genes = rownames(counts)[1:100])


clusterLabels <- primaryCluster(clusPat$rsec)


cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes

for (xx in cUniq[1:4]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 2),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("orange", "darkseagreen3"),
                       breaks = c("0", "1"))  
  print(p)
}
