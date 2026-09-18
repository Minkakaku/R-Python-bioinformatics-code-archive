library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
source('color_pallete.R')

# set default ------- 
enrichment_out_dir <- "enrichment_result" 
# dir.create(enrichment_out_dir) 

Luminal4 <- readRDS('result_8um/Luminal_reannotation3/Luminal4_final.rds')
Luminal4$clusters %>% table()

Luminal4$annotation_monocle <- 
  Luminal4$clusters %>% 
  str_replace_all(c('0' = 'LE_0', 
                    '1' = 'LE_1', 
                    '2' = 'LE_2', 
                    '3' = 'LE_3', 
                    '4' = 'LE_4', 
                    '5' = 'LE_5'))

Luminal4$annotation_monocle %>% table()


NE <- readRDS('result_8um/NE_subtype/filtered/NE_filterd.rds')
NE$Spatial.008um_snn_res.0.25 %>% table()
NE$clusters <- as.character(NE$Spatial.008um_snn_res.0.25)
NE$clusters <- factor(NE$clusters, levels = c(0:4))
NE$clusters %>% table()

NE$annotation_monocle <-
  NE$clusters %>%
  str_replace_all(c('0' = 'NE_0',
                    '1' = 'NE_1',
                    '2' = 'NE_2',
                    '3' = 'NE_3',
                    '4' = 'NE_4'))

NE$annotation_monocle %>% table()


Tumor <- merge(Luminal4, NE)
Tumor$annotation_monocle %>% table()

Tumor <- JoinLayers(Tumor)
Tumor <- NormalizeData(Tumor) # <----------------------
Tumor <- FindVariableFeatures(Tumor, nfeatures = 3000)




# DEGs -----------------------
Tumor_DEGs <- 
  FindMarkers(
    Tumor, 
    group.by = 'project', 
    ident.1 = 'NEPC', 
    ident.2 = 'ADPC'
  )

Tumor_DEGs %>% 
  rownames_to_column('gene') %>% 
  openxlsx::write.xlsx(file.path(enrichment_out_dir, 'Tumor_DEGs.xlsx'))

NE_DEGs <- 
  Tumor_DEGs %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.5) %>% 
  dplyr::filter(pct.1 - pct.2 > 0.1)

NE_DEGs %>% 
  rownames_to_column('gene') %>% 
  openxlsx::write.xlsx(file.path(enrichment_out_dir, 'NE_DEGs_to_keep.xlsx'))


LE_DEGs <- 
  Tumor_DEGs %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC < -0.5) %>% 
  dplyr::filter(pct.2 - pct.1 > 0.1)

LE_DEGs %>% 
  rownames_to_column('gene') %>% 
  openxlsx::write.xlsx(file.path(enrichment_out_dir, 'LE_DEGs_to_keep.xlsx'))




# enrichmrnt -------------------------
## NE ---------------------
NE_DEGs_for_enrichment <- rownames(NE_DEGs)
enrich_go_NE <- 
  clusterProfiler::enrichGO(
    NE_DEGs_for_enrichment,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
    ont='ALL',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05,
    keyType = 'SYMBOL'
  )

NE_DEGs_id_for_enrichment <- 
  clusterProfiler::bitr(
    NE_DEGs_for_enrichment,
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrich_kegg_NE <- 
  clusterProfiler::enrichKEGG(
    NE_DEGs_id_for_enrichment,
    organism = "hsa", 
    keyType = "ncbi-geneid", 
    pvalueCutoff = 0.05, 
    pAdjustMethod = "BH", 
    minGSSize = 10,
    maxGSSize = 500, 
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )

enrich_go_NE@result %>% write.csv('enrichment_result/enrich_go_NE.csv')
enrich_kegg_NE@result %>% write.csv('enrichment_result/enrich_kegg_NE.csv')


## LE ---------------------
LE_DEGs_for_enrichment <- rownames(LE_DEGs)
enrich_go_LE <- 
  clusterProfiler::enrichGO(
    LE_DEGs_for_enrichment,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
    ont='ALL',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05,
    keyType = 'SYMBOL'
  )

LE_DEGs_id_for_enrichment <- 
  clusterProfiler::bitr(
    LE_DEGs_for_enrichment,
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrich_kegg_LE <- 
  clusterProfiler::enrichKEGG(
    LE_DEGs_id_for_enrichment,
    organism = "hsa", 
    keyType = "ncbi-geneid", 
    pvalueCutoff = 0.05, 
    pAdjustMethod = "BH", 
    minGSSize = 10,
    maxGSSize = 500, 
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )

enrich_go_LE@result %>% write.csv('enrichment_result/enrich_go_LE.csv')
enrich_kegg_LE@result %>% write.csv('enrichment_result/enrich_kegg_LE.csv')


# dotplot ----------------------
df_NE_go <- fortify(enrich_go_NE, showCategory = 10, split = NULL)
df_NE_go <- dplyr::mutate(df_NE_go, x = eval(parse(text = 'GeneRatio')))

df_LE_go <- fortify(enrich_go_LE, showCategory = 10, split = NULL)
df_LE_go$GeneRatio <- -df_LE_go$GeneRatio
df_LE_go <- dplyr::mutate(df_LE_go, x = eval(parse(text = 'GeneRatio')))

df_NE_kegg <- fortify(enrich_kegg_NE, showCategory = 10, split = NULL)
df_NE_kegg <- dplyr::mutate(df_NE_kegg, x = eval(parse(text = 'GeneRatio')))
df_NE_kegg$ONTOLOGY <- c(rep('KEGG', dim(df_NE_kegg)[1]))

df_LE_kegg <- fortify(enrich_kegg_LE, showCategory = 10, split = NULL)
df_LE_kegg$GeneRatio <- -df_LE_kegg$GeneRatio
df_LE_kegg <- dplyr::mutate(df_LE_kegg, x = eval(parse(text = 'GeneRatio')))
df_LE_kegg$ONTOLOGY <- c(rep('KEGG', dim(df_LE_kegg)[1]))


df <- rbind(
  df_NE_go[c('GeneRatio', 'Description', 'Count', 'pvalue', 'ONTOLOGY')], 
  df_LE_go[c('GeneRatio', 'Description', 'Count', 'pvalue', 'ONTOLOGY')], 
  df_NE_kegg[c('GeneRatio', 'Description', 'Count', 'pvalue', 'ONTOLOGY')], 
  df_LE_kegg[c('GeneRatio', 'Description', 'Count', 'pvalue', 'ONTOLOGY')]
  )
idx <- order(df[['GeneRatio']], decreasing = T)
df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))

p <- 
  ggplot(df, aes(x = GeneRatio, y = Description)) + 
  # 添加从y轴到点的水平线段
  geom_segment(
    aes(x = 0, xend = GeneRatio, y = Description, yend = Description),
    color = ifelse(df$GeneRatio > 0, '#e8a0b8', 'skyblue'),  # 设置线段颜色为浅灰色
    linewidth = 0.5     # 设置线段粗细
  ) + 
  geom_segment(
    aes(x = 0, xend = 0, y = 0, yend = length(unique(df$Description))+1), 
    color = 'darkgrey',
    linewidth = 0.5  # 加粗Y轴
  ) + 
  geom_point(aes(size = Count, color = pvalue)) +
  scale_size_area() +
  scale_color_continuous(low = 'red', high = 'orange') + 
  xlab('') + 
  theme(
    axis.text.y = element_text(color = ifelse(df$ONTOLOGY == 'BP', 'darkgreen', ifelse(df$ONTOLOGY == 'CC', 'purple', ifelse(df$ONTOLOGY == 'MF', 'black', 'red')))), 
    panel.grid.major = element_blank(), #主网格线
    panel.grid.minor = element_blank(), #次网格线
    panel.border = element_rect(color = "darkgrey", size = 1, fill = NA), #边框
    legend.position = c(1.2, 0.65), 
    plot.margin = margin(1,5,1.5,0.5,'cm'), 
    panel.background = element_rect(fill = 'white')
  )
p

p2 <- 
  p + 
  coord_cartesian(clip = 'off', xlim = c(-0.08, 0.14), ylim = c(0.6, length(unique(df$Description))+0.5)) + #在非图形区域绘图,且要定好y轴范围
  annotate('text',x=0.14,y=-2,label='NE',size=5,color='#E377C2FF') + 
  annotate('text',x=-0.06,y=-2,label='LE',size=5,color='#29A2C6') + 
  annotate('text',x=0.04,y=-2,label='<-----Richfactor----->',size=5,color='black') + 
  annotate('text',x=0.2,y=4.7,label='Pathways',size=4,family='sans',color='black') + 
  annotate('text',x=0.2,y=3.5,label='GOBP',size=4,color='darkgreen') + 
  annotate('text',x=0.2,y=2.5,label='GOCC',size=4,color='purple') + 
  annotate('text',x=0.2,y=1.5,label='GOMF',size=4,color='black') + 
  annotate('text',x=0.2,y=0.5,label='KEGG',size=4,color='red')
p2

pdf(file.path(enrichment_out_dir, 'enrichment dotplot diff - NE vs LE.pdf'), width = 7.5, height = 6)
p2
dev.off()




# heatmap --------------------------------需要客户再归纳总结一下功能分几类
pathyway_info <- 
  openxlsx::read.xlsx('enrichment_result/pathway.xlsx') %>% 
  select(pathway, Type) %>% 
  drop_na()
all_description <- 
  c(str_c('KEGG_', enrich_kegg_NE@result$Description), 
    str_c('KEGG_', enrich_kegg_LE@result$Description), 
    str_c('GO_', enrich_go_NE@result$ONTOLOGY, '_', enrich_go_NE@result$Description), 
    str_c('GO_', enrich_go_LE@result$ONTOLOGY, '_', enrich_go_LE@result$Description)) %>% 
  unique() %>% 
  toupper() %>% 
  str_replace_all(' ', '_')

pathyway_in_paper <- all_description[all_description %in% (pathyway_info$pathway)]
pathyway_info[pathyway_info$pathway %in% pathyway_in_paper,]




# circle图 ---------------------
library(circlize)
gene_description <- read_tsv("enrichment_result/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

tf_CistromeDB <- read_tsv("enrichment_result/TF_CistromeDB.txt", col_names = F)
table(tf_gene %in% tf_CistromeDB$X1)

# tf_TRRUST <- openxlsx::read.xlsx('enrichment_result/Human genes included in TRRUST (25052017.xlsx')
# table(tf_gene %in% tf_TRRUST$Gene.symbol)
# 
# tf_HOCOMOCO <- read_tsv("enrichment_result/HUMAN_mono_motifs_HOCOMOCO.tsv")
# table(tf_gene %in% str_replace(tf_HOCOMOCO$`Transcription factor`, 'HUMAN:', ''))
# 
# table(rownames(NE_DEGs) %in% tf_gene)
# table(rownames(NE_DEGs) %in% tf_CistromeDB$X1)
# table(rownames(NE_DEGs) %in% tf_TRRUST$Gene.symbol)
# table(rownames(LE_DEGs) %in% tf_gene)
# table(rownames(LE_DEGs) %in% tf_CistromeDB$X1)
# table(rownames(LE_DEGs) %in% tf_TRRUST$Gene.symbol)
# rownames(NE_DEGs)[rownames(NE_DEGs) %in% tf_gene]
# rownames(NE_DEGs)[rownames(NE_DEGs) %in% tf_CistromeDB$X1]
# rownames(NE_DEGs)[rownames(NE_DEGs) %in% tf_TRRUST$Gene.symbol]
# rownames(NE_DEGs)[rownames(NE_DEGs) %in% str_replace(tf_HOCOMOCO$`Transcription factor`, 'HUMAN:', '')]
# rownames(LE_DEGs)[rownames(LE_DEGs) %in% tf_gene]
# rownames(LE_DEGs)[rownames(LE_DEGs) %in% tf_CistromeDB$X1]
# rownames(LE_DEGs)[rownames(LE_DEGs) %in% tf_TRRUST$Gene.symbol]
# rownames(LE_DEGs)[rownames(LE_DEGs) %in% str_replace(tf_HOCOMOCO$`Transcription factor`, 'HUMAN:', '')]


# Cistrome的tf较全，结合biomart的tf进行tf的筛选，TRRUST的TF NCBI复核不准
# tf的target gene(TFBS)通过TRRUST数据库查找，有对应ref文献链接
tf_list <- 
  c(tf_gene, tf_CistromeDB$X1) %>% 
  unique()

NE_TF <- rownames(NE_DEGs)[rownames(NE_DEGs) %in% tf_list]
LE_TF <- rownames(LE_DEGs)[rownames(LE_DEGs) %in% tf_list]

# NE_TF_ls <- 
#   list(
#     'ASCL1' = c('CHRNA3', 'CHRNA4', 'CHRNA5', 'DKK1', 'IGF2', 'PCSK6', 'TPH1'), 
#     'HES6' = c('ASCL1', 'HTR1A'), 
#     'PPARGC1A' = c('ALDOB', 'APOA5', 'APOC3', 'CYP2C9', 'CYP7A1', 'MFN2', 'MLYCD', 'OTC', 'PCK2'), 
#     'PROX1' = c('AKT2', 'CCNE1', 'CDKN1B', 'CDKN1C', 'CYP7A1', 'IFNG', 'NR1I2', 'NR5A2', 'RGS4', 'SIX3', 'TWIST1')
#   )
# 
# NE_TF_pair <- 
#   NE_TF_ls %>% 
#   enframe(name = 'TF', value = 'TFBS') %>% 
#   unnest(cols = c(TFBS))
# 
# LE_TF_ls <- 
#   list(
#     'NKX3-1' = c('ACTG2', 'AR', 'CCND1', 'ESR1', 'FOS', 'IGF1R', 'KLK3', 'SIX6', 'TMPRSS2', 'VEGFA'), 
#     'ZNF217' = c('CDH1'), 
#     'ELK4' = c('FOS', 'MCL1', 'NR2C2'), 
#     'GRHL2' = c('ZEB1'), 
#     'FOXA1' = c('ABCA1', 'AGR2', 'APOB', 'BCL2', 'BRCA1', 'CCNG2', 'CFTR', 'ESR1', 'GCG', 'HSPA1A', 'KRT7', 'MYC', 'NF1', 'NKX2-1', 'RET', 'RPRM', 'SERPINC1', 'SFTPB', 'TFF1', 'UGT2B15', 'UGT2B17')
#   )
# 
# LE_TF_pair <- 
#   LE_TF_ls %>% 
#   enframe(name = 'TF', value = 'TFBS') %>% 
#   unnest(cols = c(TFBS))


# Spearman’s correlation coefficients筛选TF对 <------------------------
# NE_TF_allgene <- 
#   c(NE_TF_pair$TF, NE_TF_pair$TFBS) %>% 
#   unique()
# LE_TF_allgene <- 
#   c(LE_TF_pair$TF, LE_TF_pair$TFBS) %>% 
#   unique()
# TF_allgene <- 
#   c(NE_TF_pair$TF, LE_TF_pair$TF, NE_TF_pair$TFBS, LE_TF_pair$TFBS) %>% 
#   unique()
# all_TF <- 
#   c(NE_TF_pair$TF, LE_TF_pair$TF) %>% 
#   unique()
# all_TFBS <- 
#   c(NE_TF_pair$TFBS, LE_TF_pair$TFBS) %>% 
#   unique()
all_TF <- c(NE_TF, LE_TF) %>% unique()
all_TFBS <- c(rownames(NE_DEGs), rownames(LE_DEGs)) %>% unique()
all_TFBS <- setdiff(all_TFBS, all_TF)
NE_expr <- LayerData(NE, layer = 'data')
LE_expr <- LayerData(Luminal4, layer = 'data')
# all_expr <- cbind(NE_expr, LE_expr)

NE_TF_expr <- NE_expr[(rownames(NE_expr) %in% all_TF), ]
NE_TFBS_expr <- NE_expr[(rownames(NE_expr) %in% all_TFBS), ]
NE_spearman_coef <- cor(as.matrix(t(NE_TF_expr)), as.matrix(t(NE_TFBS_expr)), method = 'spearman')

LE_TF_expr <- LE_expr[(rownames(LE_expr) %in% all_TF), ]
LE_TFBS_expr <- LE_expr[(rownames(LE_expr) %in% all_TFBS), ]
LE_spearman_coef <- cor(as.matrix(t(LE_TF_expr)), as.matrix(t(LE_TFBS_expr)), method = 'spearman')

# #####如果有NA,且LE和NE的NA位置一致#####
# # NE_TFBS_kept <- all_TFBS[all_TFBS %in% colnames(NE_spearman_coef)]
# # NE_spearman_coef <- NE_spearman_coef[all_TF, NE_TFBS_kept]
# NE_spearman_coef <- NE_spearman_coef[, colSums(is.na(NE_spearman_coef)) != nrow(NE_spearman_coef)]
# 
# # LE_TFBS_kept <- all_TFBS[all_TFBS %in% colnames(LE_spearman_coef)]
# # LE_spearman_coef <- LE_spearman_coef[all_TF, LE_TFBS_kept]
# LE_spearman_coef <- LE_spearman_coef[, colSums(is.na(LE_spearman_coef)) != nrow(LE_spearman_coef)]
# #####如果有NA,且LE和NE的NA位置一致#####

NE_spearman_coef_p <- matrix(
  NA, # 初始填充NA
  nrow = nrow(NE_TF_expr), 
  ncol = nrow(NE_TFBS_expr) 
)
for (i in seq(1, nrow(NE_TF_expr))) {
  for (j in seq(1, nrow(NE_TFBS_expr))) { # browser()
    # print(str_glue('i: {i}    j:{j}'))
    tmp <- cor.test(as.matrix(NE_TF_expr)[i,], as.matrix(NE_TFBS_expr)[j,], method = 'spearm')
    NE_spearman_coef_p[i, j] <- tmp$p.value
  }
}
rownames(NE_spearman_coef_p) <- rownames(NE_TF_expr)
colnames(NE_spearman_coef_p) <- rownames(NE_TFBS_expr)
NE_spearman_coef_p
#NE_spearman_coef_p <- NE_spearman_coef_p[, colSums(is.na(NE_spearman_coef_p)) != nrow(NE_spearman_coef_p)]


LE_spearman_coef_p <- matrix(
  NA, # 初始填充NA
  nrow = nrow(LE_TF_expr), 
  ncol = nrow(LE_TFBS_expr) 
)
for (i in seq(1,nrow(LE_TF_expr))) {
  for (j in seq(1, nrow(LE_TFBS_expr))) { # browser()
    # print(str_glue('i: {i}    j:{j}'))
    tmp <- cor.test(as.matrix(LE_TF_expr)[i,], as.matrix(LE_TFBS_expr)[j,], method = 'spearm')
    LE_spearman_coef_p[i, j] <- tmp$p.value
  }
}
rownames(LE_spearman_coef_p) <- rownames(LE_TF_expr)
colnames(LE_spearman_coef_p) <- rownames(LE_TFBS_expr)
LE_spearman_coef_p
#LE_spearman_coef_p <- LE_spearman_coef_p[, colSums(is.na(LE_spearman_coef_p)) != nrow(LE_spearman_coef_p)]


# spearman_coef_mx <- matrix(
#   NA, # 初始填充NA
#   nrow = dim(NE_spearman_coef_p)[1],
#   ncol = dim(NE_spearman_coef_p)[2]
# )
# for (i in seq(1,dim(NE_spearman_coef_p)[1])) {
#   for (j in seq(1, dim(NE_spearman_coef_p)[2])) { # browser()
#     print(str_glue('i: {i}    j:{j}'))
#     if(NE_spearman_coef_p[i,j] < 0.05 & LE_spearman_coef_p[i,j] < 0.05){
#       if(NE_spearman_coef[i,j] > 0 & LE_spearman_coef[i,j] > 0){
#         spearman_coef_mx[i,j] <- 'Both positive'
#       }else if(NE_spearman_coef[i,j] > 0 & LE_spearman_coef[i,j] <= 0){
#         spearman_coef_mx[i,j] <- 'NE positive-LE negative'
#       }else if(NE_spearman_coef[i,j] <= 0 & LE_spearman_coef[i,j] > 0){
#         spearman_coef_mx[i,j] <- 'LE positive-NE negative'
#       }else if(NE_spearman_coef[i,j] < 0 & LE_spearman_coef[i,j] < 0){
#         spearman_coef_mx[i,j] <- 'Both negative'
#       }
#     }else if(NE_spearman_coef_p[i,j] < 0.05 & LE_spearman_coef_p[i,j] >= 0.05){
#       if(NE_spearman_coef[i,j] > 0){
#         spearman_coef_mx[i,j] <- 'NE positive'
#       }else if(NE_spearman_coef[i,j] < 0){
#         spearman_coef_mx[i,j] <- 'NE negative'
#       }
#     }else if(NE_spearman_coef_p[i,j] >= 0.05 & LE_spearman_coef_p[i,j] < 0.05){
#       if(LE_spearman_coef[i,j] > 0){
#         spearman_coef_mx[i,j] <- 'LE positive'
#       }else if(LE_spearman_coef[i,j] < 0){
#         spearman_coef_mx[i,j] <- 'LE negative'
#       }
#     }else if(NE_spearman_coef_p[i,j] >= 0.05 & LE_spearman_coef_p[i,j] >= 0.05){
#       spearman_coef_mx[i,j] <- 'Non-siganificant'
#     }
#   }
# }
# 
# rownames(spearman_coef_mx) <- rownames(NE_spearman_coef_p)
# colnames(spearman_coef_mx) <- colnames(NE_spearman_coef_p)
# spearman_coef_mx


dat_NE <- (NE_spearman_coef_p < 0.05) * sign(NE_spearman_coef)
dat_NE[dat_NE == 0] <- "NE_Non"
dat_NE[dat_NE == "1"] <- "NE_Pos"
dat_NE[dat_NE == "-1"] <- "NE_Neg"

dat_LE <- (LE_spearman_coef_p < 0.05) * sign(LE_spearman_coef)
dat_LE[dat_LE == 0] <- "LE_Non"
dat_LE[dat_LE == "1"] <- "LE_Pos"
dat_LE[dat_LE == "-1"] <- "LE_Neg"

stopifnot(nrow(dat_NE) == nrow(dat_LE))
stopifnot(ncol(dat_NE) == ncol(dat_LE))
stopifnot(all(rownames(dat_NE) == rownames(dat_LE)))
stopifnot(all(colnames(dat_NE) == colnames(dat_LE)))

spearman_coef_mx2 <- matrix(
  str_c(dat_NE, dat_LE, sep = "-" ),
  nrow = nrow(dat_NE),
  ncol = ncol(dat_NE),
  dimnames = dimnames(dat_NE)
)

df2 <- 
  spearman_coef_mx2 %>% 
  as.data.frame() %>% 
  rownames_to_column('TF') %>% 
  pivot_longer(
    cols = -TF, 
    names_to = 'TFBS', 
    values_to = 'Classification'
  )
df2$label <- 
  str_replace_all(df2$Classification, c(
    'NE_Neg-LE_Neg' = NA, 
    'NE_Neg-LE_Non' = 'N', 
    'NE_Neg-LE_Pos' = NA, 
    'NE_Non-LE_Neg' = 'L', 
    'NE_Non-LE_Non' = NA, 
    'NE_Non-LE_Pos' = 'L', 
    'NE_Pos-LE_Neg' = NA, 
    'NE_Pos-LE_Non' = 'N', 
    'NE_Pos-LE_Pos' = NA
    ))
# df2$TF <- factor(df2$TF, levels = c('ASCL1', 'HES6', 'PPARGC1A', 'PROX1', 'NKX3-1', 'ZNF217', 'ELK4', 'GRHL2', 'FOXA1'))
# df2$TFBS <- factor(df2$TFBS, levels = c(
#   'CHRNA3', 'CHRNA4', 'CHRNA5', 'DKK1', 'IGF2', 'PCSK6', 'TPH1', 
#   'ASCL1', 'HTR1A', 
#   'ALDOB', 'APOA5', 'APOC3', 'CYP2C9', 'CYP7A1', 'MFN2', 'MLYCD', 'OTC', 'PCK2', 
#   'AKT2', 'CCNE1', 'CDKN1B', 'CDKN1C', 'IFNG', 'NR1I2', 'NR5A2', 'RGS4', 'SIX3', 'TWIST1', 
#   'ACTG2', 'AR', 'CCND1', 'ESR1', 'FOS', 'IGF1R', 'KLK3', 'SIX6', 'TMPRSS2', 'VEGFA', 
#   'CDH1', 
#   'MCL1', 'NR2C2', 
#   'ZEB1', 
#   'ABCA1', 'AGR2', 'APOB', 'BCL2', 'BRCA1', 'CCNG2', 'CFTR', 'GCG', 'HSPA1A', 'KRT7', 'MYC', 'NF1', 'NKX2-1', 'RET', 'RPRM', 'SERPINC1', 'SFTPB', 'TFF1', 'UGT2B15', 'UGT2B17'
# ))

df2$TF <- factor(df2$TF, levels = all_TF)
df2$TFBS <- factor(df2$TFBS, levels = all_TFBS)

TF_headtmap <- 
  ggplot(df2, aes(x = TFBS, y = TF, fill = Classification)) + 
  geom_tile() + 
  scale_fill_manual(values = c(
    'NE_Non-LE_Non' = 'grey', 
    'NE_Pos-LE_Non' = 'pink', 
    'NE_Non-LE_Pos' = 'pink', 
    'NE_Neg-LE_Non' = 'skyblue', 
    'NE_Non-LE_Neg' = 'skyblue', 
    'NE_Pos-LE_Pos' = 'red', 
    'NE_Neg-LE_Neg' = 'blue', 
    'NE_Pos-LE_Neg' = 'purple', 
    'NE_Neg-LE_Pos' = 'purple'
  )
  )+ 
  geom_text(aes(label = label)) +
  theme(
    axis.text.y = element_text(color = ifelse(all_TF %in% NE_TF, 'red', 'skyblue')), 
    axis.text.x = element_text(color = ifelse(all_TFBS %in% setdiff(rownames(NE_DEGs), all_TF), 'red', 'skyblue'), 
                               angle = 90)
  ) 

pdf(file.path(enrichment_out_dir, 'TF_headtmap.pdf'), width = 35, height = 3)
TF_headtmap
dev.off()


# ComplexHeatmap::Heatmap(spearman_coef)
# 
# df <- 
#   spearman_coef %>% 
#   as.data.frame() %>% 
#   rownames_to_column('TF') %>% 
#   pivot_longer(
#     cols = -TF, 
#     names_to = 'TFBS', 
#     values_to = 'spearman_coef'
#   )
# row_cluster <- hclust(dist(spearman_coef))
# col_cluster <- hclust(dist(t(spearman_coef)))
# df_clustered <- df
# df_clustered$TF <- factor(df_clustered$TF, levels = rownames(spearman_coef)[row_cluster$order])
# df_clustered$TFBS <- factor(df_clustered$TFBS, levels = colnames(spearman_coef)[col_cluster$order])
# 
# ggplot(df_clustered, aes(x = TFBS, y = TF, fill = spearman_coef)) + 
#   geom_tile() + 
#   scale_fill_gradient2(
#     low = 'blue', 
#     mid = 'white', 
#     high = 'red', 
#     midpoint = 0, 
#     limits = c(-0.4, 0.4), 
#     name = 'spearman_coef'
#   ) + 
#   theme(
#     axis.text.y = element_text(color = ifelse(rownames(spearman_coef) %in% c('ASCL1', 'HES6', 'PPARGC1A', 'PROX1'), 'red', 'skyblue')), 
#     axis.text.x = element_text(color = ifelse(colnames(spearman_coef) %in% c('CHRNA3', 'CHRNA4', 'CHRNA5', 'DKK1', 'IGF2', 'PCSK6', 'TPH1', 
#                                                                              'ASCL1', 'HTR1A', 
#                                                                              'ALDOB', 'APOA5', 'APOC3', 'CYP2C9', 'CYP7A1', 'MFN2', 'MLYCD', 'OTC', 'PCK2', 
#                                                                              'AKT2', 'CCNE1', 'CDKN1B', 'CDKN1C', 'CYP7A1', 'IFNG', 'NR1I2', 'NR5A2', 'RGS4', 'SIX3', 'TWIST1'
#                                                                              ), 'red', 'skyblue'), 
#                                angle = 90)
#   )  




TF_pair_filtered <- 
  df2 %>% 
  filter(Classification %in% c('NE_Pos-LE_Pos', 'NE_Neg-LE_Neg'))
TF_order <- TF_pair_filtered$TF %>% table() %>% sort(decreasing = T)
# TF_pair_filtered$TF <- factor(TF_pair_filtered$TF, levels = names(TF_order))
TF_pair_filtered_sorted <- TF_pair_filtered %>% arrange(factor(TF, levels = names(TF_order)))

mycol1 <- RColorBrewer::brewer.pal(length(unique(TF_pair_filtered_sorted$TF)), 'Set3')
names(mycol1) <- unique(TF_pair_filtered_sorted$TF)
mycol2 <- c(rep('black', length(unique(TF_pair_filtered_sorted$TFBS))))
names(mycol2) <- unique(TF_pair_filtered_sorted$TFBS)

pdf(file.path(enrichment_out_dir, 'TF_chordDiagram.pdf'), width = 12, height = 12)
circos.par(start.degree = 90, # 设置circos图的起始角度为90度
           gap.degree = 1) # 设置相邻扇区之间的间隙为5度
chordDiagram(
  x = TF_pair_filtered_sorted, # 输入数据
  grid.col = c(mycol1, mycol2), # TF和TFBS的颜色
  # grid.border = 'black', 
  transparency = 0, # 透明度
  directional = 1, # 线的方向 1指从绘图数据第一列指向绘图数据第二列
  # direction.type = c("diffHeight+arrows"),
  direction.type = c("diffHeight"),
  link.arr.length = 0.1, 
  link.arr.type = 'curved', 
  big.gap = 8, #只有当TFBS中没有TF中的基因的情况下才生效
  # small.gap = 5,
  diffHeight = -0.04,
  annotationTrack = c("grid", "axis"), # 注释轨道，'grid'指扇形网格
  annotationTrackHeight = c(0.05, 0.1), # 注释轨道高度
  link.sort = 'default', # 调整线的次序，使得线具有最小的交集
  link.decreasing = TRUE,
  preAllocateTracks = list(track.height = 0.2))#调整circle大小

# 添加文本注释
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, cex=0.6,adj = c(-0.5, 0.5))
}, bg.border = NA)
circos.clear()
dev.off()




# radar plot --------------
library(plotly)
# min_radius <-  10
# 
# radius <- seq(10, 110, by = 50)
# theta_circle <- seq(0, 2 * pi, length.out = 360)
# 
# df_circle <-
#   radius %>%
#   setNames(., .) %>%
#   map(function(x)list(x = x * cos(theta_circle), y = x * sin(theta_circle))) %>%
#   map(as.data.frame) %>%
#   enframe('radius') %>%
#   unnest(value)

pathway <- openxlsx::read.xlsx('enrichment_result/pathway.xlsx')
'HALLMARK_DNA_REPAIR' %in% pathway$pathway
'GOBP_SPLICEOSOMAL_SNRNP_ASSEMBLY' %in% pathway$pathway
'*_SPLICEOSOMAL_SNRNP_ASSEMBLY' %in% pathway$pathway

enrich_go_NE_result <- read_csv(file.path(enrichment_out_dir, 'enrich_go_NE.csv'))
enrich_go_NE_pathway <- 
  enrich_go_NE_result %>% 
  filter(pvalue < 0.05 & qvalue < 0.25) %>% 
  mutate(
    path = str_c('GO', ONTOLOGY, '_', toupper(Description))
    ) %>% 
  .$path  %>% 
  str_replace_all(c(' ' = '_'))
enrich_go_NE_pathway %in% pathway$pathway %>% table()

enrich_go_LE_result <- read_csv(file.path(enrichment_out_dir, 'enrich_go_LE.csv'))
enrich_go_LE_pathway <- 
  enrich_go_LE_result %>% 
  filter(pvalue < 0.05 & qvalue < 0.25) %>% 
  mutate(
    path = str_c('GO', ONTOLOGY, '_', toupper(Description))
  ) %>% 
  .$path  %>% 
  str_replace_all(c(' ' = '_'))
enrich_go_LE_pathway %in% pathway$pathway %>% table()

enrich_kegg_NE_result <- read_csv(file.path(enrichment_out_dir, 'enrich_kegg_NE.csv'))
enrich_kegg_NE_pathway <- 
  enrich_kegg_NE_result %>% 
  filter(pvalue < 0.05 & qvalue < 0.25) %>% 
  mutate(
    path = str_c('KEGG', '_', toupper(Description))
  ) %>% 
  .$path  %>% 
  str_replace_all(c(' ' = '_'))
enrich_kegg_NE_pathway %in% pathway$pathway %>% table()

enrich_kegg_LE_result <- read_csv(file.path(enrichment_out_dir, 'enrich_kegg_LE.csv'))
enrich_kegg_LE_pathway <- 
  enrich_kegg_LE_result %>% 
  filter(pvalue < 0.05 & qvalue < 0.25) %>% 
  mutate(
    path = str_c('KEGG', '_', toupper(Description))
  ) %>% 
  .$path  %>% 
  str_replace_all(c(' ' = '_'))
enrich_kegg_LE_pathway %in% pathway$pathway %>% table()

pathway_to_keep_NE <- c(
  enrich_go_NE_pathway[enrich_go_NE_pathway %in% pathway$pathway], 
  enrich_kegg_NE_pathway[enrich_kegg_NE_pathway %in% pathway$pathway]
  )
pathway_to_keep_LE <- c(
  enrich_go_LE_pathway[enrich_go_LE_pathway %in% pathway$pathway], 
  enrich_kegg_LE_pathway[enrich_kegg_LE_pathway %in% pathway$pathway]
)

df_NE <- 
  pathway[pathway$pathway %in% pathway_to_keep_NE, ] %>% 
  dplyr::select(pathway, Type) %>% 
  mutate(group = 'NE')
df_NE
df_LE <- 
  pathway[pathway$pathway %in% pathway_to_keep_LE, ] %>% 
  dplyr::select(pathway, Type) %>% 
  mutate(group = 'LE')
df_LE

df <- 
  rbind(df_NE, df_LE) %>% 
  group_by(Type) %>% 
  mutate(count_Type = n()) %>% 
  group_by(group, Type) %>% 
  mutate(count = n()) %>% 
  mutate(freq = 100*count/count_Type) %>% 
  dplyr::select(-pathway) %>% 
  unique() 

# df %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select(Type, group, freq) %>% 
#   tidyr::complete(Type, group, fill = list(freq = 0)) %>% 
#   arrange(Type) %>%
#   ggplot(aes(x = Type, y = freq, color = group)) + 
#   geom_vline(xintercept = 0:4 , color = "gray") +
#   geom_hline(yintercept = c(0, 50, 100), color = "gray") + 
#   geom_point() +
#   # geom_path(
#   #   aes(group = group)
#   # ) +
#   # geom_path(
#   #   data = ~ mutate(rbind(.[1:2, ], .[1:2, ]), Type = c(0, 0, 1, 1) )%>% print,
#   #   aes(group = group)
#   # ) +
#   geom_text(
#     data = data.frame(x = 0, y = c(0, 50, 100)),
#     aes(x = x, y = y, label = paste0(y, "%")),
#     inherit.aes = FALSE,
#     vjust = -1.5
#   ) +
#   theme_minimal() + 
#   theme(
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     plot.margin = margin(1,1,1,6, 'mm')
#   ) +
#   scale_size(range = c(3, 12)) +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-10, 100)) +
#   coord_polar(theta = 'x', clip = 'off', start = -1/2 * pi)
# 
# 
# df %>% 
#   dplyr::ungroup() %>% 
#   dplyr::select(Type, group, freq) %>% 
#   tidyr::complete(Type, group, fill = list(freq = 0)) %>% 
#   mutate(freq = freq + min_radius) %>% 
#   mutate(
#     theta = scales::rescale(as.numeric(factor(Type)), to = c(0, 2*pi)),
#     xx = freq * cos(theta),
#     yy = freq * sin(theta)
#   ) %>%
#   
#   ggplot(aes(x = xx, y = yy, color = group, fill = group)) +
# 
#   geom_polygon(
#     data = df_circle,
#     aes(x = x, y = y, group = radius),
#     color = "blue",
#     fill = NA,
#     inherit.aes = FALSE
#   ) + 
#   theme_void() + 
#   geom_polygon(alpha = 0.5) +
#   coord_equal()

df2 <- 
  df %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Type, group, freq) %>% 
  tidyr::complete(Type, group, fill = list(freq = 0))

p <- 
  plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  ) %>% 
  add_trace(
    r = df2[df2$group == 'LE', ]$freq,
    theta = df2[df2$group == 'LE', ]$Type,
    name = 'LE'
  )%>% 
  add_trace(
    r = df2[df2$group == 'NE', ]$freq,
    theta = df2[df2$group == 'NE', ]$Type,
    name = 'NE'
  ) %>%
  layout(
    polar = list(
      radialaxis = list(
        range = c(-10, 100),  # 设置半径范围为 -10 到 100
        visible = TRUE
      )
    )
  )


library(reticulate)
use_python("/home/szt/miniforge3/envs/rpy/bin/python")
save_image(p, file.path(enrichment_out_dir, "radar_plot.pdf"), engine = "orca", width = 800, height = 600)




# gene-pair scatterplot --------------------
# Fisher Z变换函数
fisher_z_transform <- function(rho) {
  0.5 * log((1 + rho) / (1 - rho))
}

# 完整的工作流程：从数据到Fisher Z变换的p值
fisher_z_test <- function(x, y, alternative = "two.sided") {
  # 移除缺失值
  complete_cases <- complete.cases(x, y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  n <- length(x_clean)
  
  # 计算Spearman相关系数
  spearman_rho <- cor(x_clean, y_clean, method = "spearm")
  
  # Fisher Z变换
  fisher_z <- fisher_z_transform(spearman_rho)
  
  # 计算标准误
  se <- 1 / sqrt(n - 3)
  
  # 计算Z-score（零假设：总体相关系数为0）
  z_score <- fisher_z / se
  
  # 返回结果
  return(z_score)
}

# 比较两个独立Z-score的差异
compare_two_zscores <- function(x1, y1, x2, y2, alternative = "two.sided") {
  complete_cases1 <- complete.cases(x1, y1)
  x_clean1 <- x1[complete_cases1]
  n1 <- length(x_clean1)
  
  complete_cases2 <- complete.cases(x2, y2)
  x_clean2 <- x2[complete_cases2]
  n2 <- length(x_clean2)
  
  z1 <- fisher_z_test(x1, y1)
  z2 <- fisher_z_test(x2, y2)
  
  # 计算两个Z-score差异的标准误
  se_diff <- sqrt(1/(n1-3) + 1/(n2-3))
  
  # 计算检验统计量
  z_stat <- (z1 - z2) / se_diff
  
  # 计算p值
  if (alternative == "two.sided") {
    p_value <- 2 * (1 - pnorm(abs(z_stat)))
  } else if (alternative == "greater") { # 右侧检验
    p_value <- 1 - pnorm(z_stat)
  } else if (alternative == "less") { # 左侧检验
    p_value <- pnorm(z_stat)
  } else {
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  }
  
  return(list(
    z1 = z1,
    z2 = z2,
    z_difference = z1 - z2,
    z_statistic = z_stat,
    p_value = p_value,
    standard_error = se_diff,
    confidence_interval = c(
      (z1 - z2) - 1.96 * se_diff,
      (z1 - z2) + 1.96 * se_diff
    )
  ))
}

# gene_pairs <- c(
#   'AR', 'KLK3', 'TMPRSS2', 'NKX3.1', 
#   'ASCL1', 'CHGA', 
#   'NEUROD1', 'SYP', 
#   'INSM1', 'ENO2', 
#   'E2F1', 'MKI67', 
#   'MYCN', 'AURKA', 
#   'CDK4', 'CCND1', 
#   'BRD4', 'MYC', 'MYCN', 
#   'KDM1A', 'LSD1', 'SOX2', 
#   'HDAC1', 'REST', 
#   'PTEN', 'CCND1', 
#   'AKT1', 'mTOR', 
#   'RICTOR', 'AKT1', 
#   'CTNNB1', 'AXIN2', 
#   'LEF1', 'SOX2', 
#   'DKK1', 'CTNNB1'
# )
gene_pair_ls <- c(
  'AR' = 'KLK3', 
  'AR' = 'TMPRSS2', 
  'AR' = 'NKX3-1', 
  'ASCL1' = 'CHGA', 
  'NEUROD1' = 'SYP', 
  'INSM1' = 'ENO2', 
  'E2F1' = 'MKI67', 
  'MYCN' = 'AURKA', 
  'CDK4' = 'CCND1', 
  'BRD4' = 'MYC',
  'KDM1A' = 'SOX2', 
  'HDAC1' = 'REST', 
  'PTEN' = 'CCND1', 
  'AKT1' = 'MTOR', 
  'RICTOR' = 'AKT1', 
  'CTNNB1' = 'AXIN2', 
  'LEF1' = 'SOX2', 
  'DKK1' = 'CTNNB1'
)

all_gene_pairs_expr <- Tumor %>% 
  FetchData(vars = unique(c(gene_pair_ls, names(gene_pair_ls))), layer = "data") %>% 
  rownames_to_column('cell_id')

names(gene_pair_ls) <- str_c(names(gene_pair_ls), gene_pair_ls, sep = '_')
meta_expr <- 
  Tumor@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(all_gene_pairs_expr, by = 'cell_id')

meta_expr %>% head


p_ls <- 
  gene_pair_ls %>% 
  imap(function(gene1, gene_pair){ # browser()
    gene2 <- str_split(gene_pair, '_')[[1]][1]
    
    r_ADPC <- 
      cor.test(meta_expr[meta_expr$project == 'ADPC', gene2], meta_expr[meta_expr$project == 'ADPC', gene1], method = 'spearm')$estimate %>% 
      round(3)
    r_ADPC_p <- cor.test(meta_expr[meta_expr$project == 'ADPC', gene2], meta_expr[meta_expr$project == 'ADPC', gene1], method = 'spearm')$p.value
    r_ADPC_p_sig <- 
      case_when(
        r_ADPC_p < 0.001 ~ "***",
        r_ADPC_p < 0.01 ~ "**",
        r_ADPC_p < 0.05 ~ "*", 
        TRUE ~ " ns"
      )
    
    r_NEPC <- 
      cor.test(meta_expr[meta_expr$project == 'NEPC', gene2], meta_expr[meta_expr$project == 'NEPC', gene1], method = 'spearm')$estimate %>% 
      round(3)
    r_NEPC_p <- cor.test(meta_expr[meta_expr$project == 'NEPC', gene2], meta_expr[meta_expr$project == 'NEPC', gene1], method = 'spearm')$p.value
    r_NEPC_p_sig <- 
      case_when(
        r_NEPC_p < 0.001 ~ "***",
        r_NEPC_p < 0.01 ~ "**",
        r_NEPC_p < 0.05 ~ "*",
        TRUE ~ " ns"
      )
    
    result <- compare_two_zscores(
      meta_expr[meta_expr$project == 'ADPC', gene2], 
      meta_expr[meta_expr$project == 'ADPC', gene1], 
      meta_expr[meta_expr$project == 'NEPC', gene2], 
      meta_expr[meta_expr$project == 'NEPC', gene1], 
      alternative = "two.sided")
    
    fisher_z_p <- result$p_value
    
    
    # new version ----------
    p <- 
      meta_expr %>% 
      dplyr::select(any_of(c(gene2, gene1)), project) %>% 
      # dplyr::filter(.data[[gene2]] != 0, .data[[gene1]] != 0) %>%
      ggplot(aes(x = .data[[gene1]], y = .data[[gene2]], color = project)) + 
      geom_point() + 
      # geom_smooth(method = "lm", se = TRUE, aes(fill = project), alpha = 0.1) + 
      geom_smooth(method = "lm", se = FALSE, linetype = 'dashed') + 
      scale_color_manual(values = c(
        'NEPC' = '#E377C2FF', 
        'ADPC' = '#29A2C6'
      )) + 
      theme_bw() + 
      theme(legend.position = "none")
    #############################
    
    p2 <- 
      meta_expr %>% 
      dplyr::select(any_of(gene1), project) %>% 
      
      ggplot(aes(y = project, x = .data[[gene1]], color = project)) + 
      geom_boxplot()  + 
      scale_color_manual(values = c(
        'NEPC' = '#E377C2FF', 
        'ADPC' = '#29A2C6'
      )) + 
      coord_cartesian(clip = 'off', ylim = c(0.6, 3)) + 
      annotate('text',x=0.2,y=2.9,label=str_glue('{r_NEPC}{r_NEPC_p_sig}'),size=3,color='#E377C2FF') + 
      annotate('text',x=1.8,y=2.9,label=str_glue('{r_ADPC}{r_ADPC_p_sig}'),size=3,color='#29A2C6') + 
      annotate('text',x=3.8,y=2.9,label=str_glue('fisher_z_p: {fisher_z_p}'),size=3,color='black') + 
      theme_bw() + 
      theme(legend.position = "none")
    
    p3 <- 
      meta_expr %>% 
      dplyr::select(any_of(gene2), project) %>% 
      
      ggplot(aes(x = project, y = .data[[gene2]], color = project)) + 
      geom_boxplot()  + 
      scale_color_manual(values = c(
        'NEPC' = '#E377C2FF', 
        'ADPC' = '#29A2C6'
      )) + 
      theme_bw() + 
      theme(legend.position = "none")
    
    p4 <- 
      p %>% 
      aplot::insert_top(p2 + theme_void() + theme(legend.position = "none"), height = 0.2) %>% 
      aplot::insert_right(p3 + theme_void() + theme(legend.position = "none"), width = 0.2) 
    
    p4
  })

pdf(file.path(enrichment_out_dir, 'gene-pair scatterplot_patch.pdf'), width = 20, height = 30)
p_ls %>%
  imap(function(p, gene_pair){ # browser()
    p %>% aplot::gglistGrob() %>% ggplotify::as.ggplot() + ggtitle(gene_pair)
  }) %>%
  patchwork::wrap_plots(ncol = 3)
dev.off()

pdf(file.path(enrichment_out_dir, 'gene-pair scatterplot.pdf'), width = 6, height = 6)
p_ls %>%
  imap(function(p, gene_pair){ # browser()
    p
  })
dev.off()

r_ADPC <- 
  cor.test(meta_expr[meta_expr$project == 'ADPC', 'AR'], meta_expr[meta_expr$project == 'ADPC', 'KLK3'], method = 'spearm')$estimate %>% 
  round(3)
r_ADPC_p <- cor.test(meta_expr[meta_expr$project == 'ADPC', 'AR'], meta_expr[meta_expr$project == 'ADPC', 'KLK3'], method = 'spearm')$p.value
r_ADPC_p_sig <- 
  case_when(
    r_ADPC_p < 0.001 ~ "***",
    r_ADPC_p < 0.01 ~ "**",
    r_ADPC_p < 0.05 ~ "*", 
    TRUE ~ " ns"
  )
  
r_NEPC <- 
  cor.test(meta_expr[meta_expr$project == 'NEPC', 'AR'], meta_expr[meta_expr$project == 'NEPC', 'KLK3'], method = 'spearm')$estimate %>% 
  round(3)
r_NEPC_p <- cor.test(meta_expr[meta_expr$project == 'NEPC', 'AR'], meta_expr[meta_expr$project == 'NEPC', 'KLK3'], method = 'spearm')$p.value
r_NEPC_p_sig <- 
  case_when(
    r_NEPC_p < 0.001 ~ "***",
    r_NEPC_p < 0.01 ~ "**",
    r_NEPC_p < 0.05 ~ "*",
    TRUE ~ " ns"
  )

result <- compare_two_zscores(
  meta_expr[meta_expr$project == 'ADPC', 'AR'], 
  meta_expr[meta_expr$project == 'ADPC', 'KLK3'], 
  meta_expr[meta_expr$project == 'NEPC', 'AR'], 
  meta_expr[meta_expr$project == 'NEPC', 'KLK3'], 
  alternative = "two.sided")

fisher_z_p <- result$p_value



p <- 
  meta_expr %>% 
  dplyr::select(AR, KLK3, project) %>% 
  
  ggplot(aes(x = KLK3, y = AR, color = project)) + 
  geom_point() + 
  # geom_smooth(method = "lm", se = TRUE, aes(fill = project), alpha = 0.1) + 
  geom_smooth(method = "lm", se = FALSE, linetype = 'dashed') + 
  scale_color_manual(values = c(
    'NEPC' = '#E377C2FF', 
    'ADPC' = '#29A2C6'
  )) + 
  theme_bw() + 
  theme(legend.position = "none")

p2 <- 
  meta_expr %>% 
  dplyr::select(KLK3, project) %>% 
  
  ggplot(aes(y = project, x = KLK3, color = project)) + 
  geom_boxplot()  + 
  scale_color_manual(values = c(
    'NEPC' = '#E377C2FF', 
    'ADPC' = '#29A2C6'
  )) + 
  coord_cartesian(clip = 'off', ylim = c(0.6, 3)) + 
  annotate('text',x=0.2,y=2.9,label=str_glue('{r_NEPC}{r_NEPC_p_sig}'),size=3,color='#E377C2FF') + 
  annotate('text',x=1.8,y=2.9,label=str_glue('{r_ADPC}{r_ADPC_p_sig}'),size=3,color='#29A2C6') + 
  annotate('text',x=3.8,y=2.9,label=str_glue('fisher_z_p: {fisher_z_p}'),size=3,color='black') + 
  theme_bw() + 
  theme(legend.position = "none")

p3 <- 
  meta_expr %>% 
  dplyr::select(AR, project) %>% 
  
  ggplot(aes(x = project, y = AR, color = project)) + 
  geom_boxplot()  + 
  scale_color_manual(values = c(
    'NEPC' = '#E377C2FF', 
    'ADPC' = '#29A2C6'
  )) + 
  theme_bw() + 1
  theme(legend.position = "none")

p4 <- 
  p %>% 
  aplot::insert_top(p2 + theme_void() + theme(legend.position = "none"), height = 0.2) %>% 
  aplot::insert_right(p3 + theme_void() + theme(legend.position = "none"), width = 0.2) 

p4


