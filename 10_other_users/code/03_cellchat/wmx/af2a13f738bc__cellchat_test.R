# 常见信号通路受体在肝细胞上表达都很低
# 对cellchat只能调一下参数试一下
# 1. computeCommunProb(cellchat,type="truncatedMean",trim=0,raw.use = TRUE)### default triMean  producing fewer but stronger interactions; 
# 2. computeCommunProb(cellchat,type="truncatedMean",trim=0,raw.use = FALSE)### default triMean  producing fewer but stronger interactions; 

library(CellChat)
library(dplyr)
library(Seurat)
library(ggplot2)

load("/mnt/NexData_file/home/wmx/Project/SC2512227719-01/shouhou4/results/06_细胞通讯/CKO_C/tmp_object.RData")

CKO <- object.list$CKO 
tt <- read.delim2("/mnt/NexData_file/home/wmx/Project/SC2512227719-01/shouhou4/results/06_细胞通讯/CKO_C/CKO-net_lr.txt")
tt1 <- tt %>% dplyr::filter(source == "Hepatocytes" | target == "Hepatocytes")

sce <- readRDS("/mnt/NexData_file/home/wmx/Project/SC2512227719-01/shouhou4/for_cellchat.rds")
table(sce$cell_type_sub,sce$cell_type)
subsce <- subset(sce,subset = cell_type  == "Hepatocytes")


# 定义小鼠肝细胞的核心受体列表,这里发现就是没有什么接受信号的表达
hep_receptors_mouse <- unique(c(
  # growth / regeneration
  "Met",
  "Egfr", "Erbb2", "Erbb3",
  "Fgfr1", "Fgfr2", "Fgfr3", "Fgfr4", "Klb",
  "Igf1r", "Igf2r",
  
  # inflammation
  "Il6r", "Il6st",
  "Il1r1", "Il1rap",
  "Tnfrsf1a", "Tnfrsf1b",
  "Ifngr1", "Ifngr2",
  "Ifnar1", "Ifnar2",
  "Osmr", "Lifr",
  
  # fibrosis / development
  "Tgfbr1", "Tgfbr2", "Tgfbr3",
  "Bmpr1a", "Bmpr1b", "Bmpr2",
  "Acvr1", "Acvr1b", "Acvr2a", "Acvr2b",
  "Fzd1", "Fzd2", "Fzd3", "Fzd4", "Fzd5", "Fzd6", "Fzd7",
  "Lrp5", "Lrp6",
  "Notch1", "Notch2", "Notch3",
  
  # ECM / adhesion
  "Itga1", "Itga2", "Itga3", "Itga5", "Itga6",
  "Itgav", "Itgb1",
  "Cd44",
  "Cd36", "Cd47",
  "Cd74",
  "Cxcr2", "Cxcr4",
  
  # vascular / others
  "Flt1", "Kdr", "Nrp1", "Nrp2",
  "Tek", "Tie1",
  "Plxna1", "Plxna2", "Plxna3", "Plxna4",
  "Epha2", "Ephb2", "Ephb4",
  "Ptprz1", "Sdc1", "Sdc2", "Sdc3", "Sdc4", "Ncl"
))
# 检查这些基因在你的 Seurat 对象中是否存在
available_genes <- intersect(hep_receptors_mouse, rownames(sce))

# 绘制点图（按你注释的细胞类型分组）
DotPlot(sce, features = available_genes, group.by = "cell_type") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#先看你之前整理的小鼠 hepatocyte receptor 在 data.signaling 
intersect(hep_receptors_mouse, rownames(CKO@data.signaling))
setdiff(hep_receptors_mouse, rownames(CKO@data.signaling))


# 检测一下调整参数有没有影响
CKO_test <- CKO

CKO_test <- computeCommunProb(
  CKO_test,
  type = "truncatedMean",
  trim = 0.1,
  raw.use = TRUE
)

CKO_test <- filterCommunication(CKO_test, min.cells = 10)
CKO_test <- computeCommunProbPathway(CKO_test)
CKO_test <- aggregateNet(CKO_test)

df.testp <- subsetCommunication(CKO_test, slot.name = "netP")

table(df.testp$target)
df.testp[df.testp$target == "Hepatocytes", ]

# > df.testp[df.testp$target == "Hepatocytes", ]
# source      target pathway_name         prob pval
# 1551    CD8_Tem_C1 Hepatocytes         PARs 1.332341e-04    0
# 2036    CD8_Tem_C3 Hepatocytes         PARs 2.252590e-05    0
# 2519    CD8_Tem_C4 Hepatocytes         PARs 2.630172e-05    0
# 3698     CD8_Tn_C3 Hepatocytes         PARs 1.168667e-06    0
# 4240            DC Hepatocytes          FN1 4.779054e-06    0
# 4661   Hepatocytes Hepatocytes     GALECTIN 8.159288e-03    0
# 4662   Hepatocytes Hepatocytes         PARs 4.653084e-05    0
# 5056 Kupffer cells Hepatocytes          APP 5.666865e-03    0
# 5678    Macrophage Hepatocytes          APP 8.145444e-03    0
# 5679    Macrophage Hepatocytes          FN1 1.038989e-04    0
# 5680    Macrophage Hepatocytes         THBS 6.229800e-04    0
# 6295      Monocyte Hepatocytes          APP 1.125715e-02    0
# 6296      Monocyte Hepatocytes     COLLAGEN 2.614810e-06    0
# 6297      Monocyte Hepatocytes          FN1 3.202122e-03    0
# 6298      Monocyte Hepatocytes         THBS 3.061399e-03    0
# 6867           NKT Hepatocytes         PARs 4.749906e-04    0
# 7376            NK Hepatocytes         PARs 6.235624e-06    0
# 8833    Naive B_C4 Hepatocytes         PARs 2.930005e-07    0
# 9221        Plasma Hepatocytes         PARs 2.418203e-06    0

CKO_test1 <- CKO
CKO_test1 <- smoothData(CKO_test1, adj = PPI.mouse)
CKO_test1 <- computeCommunProb(
  CKO_test1,
  type = "triMean",
  trim = 0.1,
  raw.use = FALSE
)
CKO_test1 <- filterCommunication(CKO_test1, min.cells = 10)
CKO_test1 <- computeCommunProbPathway(CKO_test1)
CKO_test1 <- aggregateNet(CKO_test1)

df.testp1 <- subsetCommunication(CKO_test1, slot.name = "netP")

table(df.testp1$target)
df.testp1[df.testp1$target == "Hepatocytes", ]
# > df.testp1[df.testp1$target == "Hepatocytes", ]
# source      target pathway_name         prob pval
# 242        CD4_Tem Hepatocytes          LCK 1.397246e-04 0.02
# 243        CD4_Tem Hepatocytes        MHC-I 5.327636e-05 0.00
# 689         CD4_Tn Hepatocytes        MHC-I 2.015707e-05 0.00
# 1090       CD8_Tcm Hepatocytes        MHC-I 2.650991e-05 0.00
# 1545    CD8_Tem_C1 Hepatocytes        MHC-I 7.277651e-05 0.00
# 2015    CD8_Tem_C3 Hepatocytes        MHC-I 3.871730e-05 0.01
# 2912     CD8_Tn_C1 Hepatocytes        MHC-I 3.081059e-05 0.00
# 3302     CD8_Tn_C2 Hepatocytes        MHC-I 3.511179e-05 0.00
# 3675     CD8_Tn_C3 Hepatocytes        MHC-I 1.742698e-05 0.00
# 4179            DC Hepatocytes       CEACAM 2.674334e-05 0.00
# 4180            DC Hepatocytes        MHC-I 2.888519e-05 0.00
# 4579   Hepatocytes Hepatocytes        MHC-I 9.022789e-05 0.00
# 4955 Kupffer cells Hepatocytes        MHC-I 3.081274e-05 0.00
# 5594    Macrophage Hepatocytes       CEACAM 3.013083e-04 0.00
# 5595    Macrophage Hepatocytes        MHC-I 1.742216e-05 0.00
# 6217      Monocyte Hepatocytes       CEACAM 1.839665e-04 0.00
# 6218      Monocyte Hepatocytes        MHC-I 2.110597e-05 0.00
# 6766           NKT Hepatocytes        MHC-I 6.429824e-05 0.00
# 7240            NK Hepatocytes        MHC-I 5.222464e-05 0.00
# 8330    Naive B_C3 Hepatocytes        MHC-I 2.909657e-05 0.00
# 9040        Plasma Hepatocytes        MHC-I 1.728986e-05 0.04

CKO_test2 <- CKO
CKO_test2 <- smoothData(CKO_test2, adj = PPI.mouse)
CKO_test2 <- computeCommunProb(
  CKO_test1,
  type = "truncatedMean",
  trim = 0.1,
  raw.use = FALSE
)
CKO_test2 <- filterCommunication(CKO_test2, min.cells = 10)
CKO_test2 <- computeCommunProbPathway(CKO_test2)
CKO_test2 <- aggregateNet(CKO_test2)

df.testp2 <- subsetCommunication(CKO_test2, slot.name = "netP")

table(df.testp2$target)
df.testp2[df.testp2$target == "Hepatocytes", ]

# > df.testp2[df.testp2$target == "Hepatocytes", ]
# source      target pathway_name         prob  pval
# 298         CD4_Tem Hepatocytes         CD40 1.432431e-07 0.000
# 299         CD4_Tem Hepatocytes         CD96 2.364054e-08 0.000
# 300         CD4_Tem Hepatocytes        FASLG 1.017649e-07 0.000
# 926          CD4_Tn Hepatocytes         CD40 2.056972e-08 0.000
# 927          CD4_Tn Hepatocytes         CD96 1.044723e-08 0.030
# 2110     CD8_Tem_C1 Hepatocytes        FASLG 1.014968e-07 0.030
# 2111     CD8_Tem_C1 Hepatocytes         PARs 8.882670e-05 0.000
# 2772     CD8_Tem_C3 Hepatocytes         PARs 1.501738e-05 0.000
# 3434     CD8_Tem_C4 Hepatocytes         PARs 1.753463e-05 0.000
# 3435     CD8_Tem_C4 Hepatocytes         PDGF 1.132737e-08 0.000
# 4529      CD8_Tn_C2 Hepatocytes         CD40 1.280158e-08 0.040
# 5054      CD8_Tn_C3 Hepatocytes         PARs 7.791114e-07 0.000
# 5778             DC Hepatocytes         CD80 4.406151e-11 0.010
# 5779             DC Hepatocytes          CDH 3.172151e-09 0.000
# 5780             DC Hepatocytes          EGF 5.400722e-10 0.000
# 5781             DC Hepatocytes          FN1 4.779054e-06 0.000
# 5782             DC Hepatocytes       NECTIN 7.380599e-09 0.000
# 5783             DC Hepatocytes         PDGF 5.877653e-09 0.000
# 6382    Hepatocytes Hepatocytes          FGF 7.438949e-09 0.000
# 6383    Hepatocytes Hepatocytes     GALECTIN 8.159288e-03 0.000
# 6384    Hepatocytes Hepatocytes         PARs 3.102080e-05 0.000
# 6895  Kupffer cells Hepatocytes          APP 5.666865e-03 0.000
# 6896  Kupffer cells Hepatocytes         CD80 3.461626e-11 0.030
# 6897  Kupffer cells Hepatocytes          CSF 6.445814e-07 0.025
# 6898  Kupffer cells Hepatocytes          EGF 2.507579e-07 0.000
# 6899  Kupffer cells Hepatocytes         PDGF 7.431826e-08 0.000
# 7779     Macrophage Hepatocytes          APP 8.145444e-03 0.000
# 7780     Macrophage Hepatocytes       CEACAM 8.452582e-04 0.000
# 7781     Macrophage Hepatocytes          CSF 2.781201e-07 0.010
# 7782     Macrophage Hepatocytes          EGF 1.454926e-09 0.000
# 7783     Macrophage Hepatocytes          FN1 1.038989e-04 0.000
# 7784     Macrophage Hepatocytes       NECTIN 6.284182e-09 0.000
# 7785     Macrophage Hepatocytes         PDGF 2.679754e-08 0.000
# 7786     Macrophage Hepatocytes         THBS 6.229800e-04 0.000
# 8661       Monocyte Hepatocytes          APP 1.125715e-02 0.000
# 8662       Monocyte Hepatocytes         CD80 5.308561e-11 0.010
# 8663       Monocyte Hepatocytes       CEACAM 5.620151e-04 0.000
# 8664       Monocyte Hepatocytes     COLLAGEN 2.614810e-06 0.000
# 8665       Monocyte Hepatocytes          CSF 2.412182e-07 0.000
# 8666       Monocyte Hepatocytes          FN1 3.202122e-03 0.000
# 8667       Monocyte Hepatocytes       NECTIN 3.214858e-09 0.000
# 8668       Monocyte Hepatocytes         THBS 3.061399e-03 0.000
# 9446            NKT Hepatocytes         PARs 3.167106e-04 0.000
# 9447            NKT Hepatocytes         PDGF 4.431037e-08 0.000
# 10151            NK Hepatocytes         PARs 4.157091e-06 0.000
# 10152            NK Hepatocytes         PDGF 2.060300e-08 0.000
# 11766    Naive B_C3 Hepatocytes       CEACAM 4.685116e-05 0.010
# 12295    Naive B_C4 Hepatocytes         PARs 1.953337e-07 0.000
# 12838        Plasma Hepatocytes         PARs 1.612137e-06 0.000
# 13413          Th17 Hepatocytes          EGF 3.791351e-08 0.000
# 13414          Th17 Hepatocytes           LT 2.557540e-08 0.000
# 14159          Treg Hepatocytes         CD40 7.142864e-08 0.000
# 14160          Treg Hepatocytes         CD96 3.417658e-08 0.010
# 14161          Treg Hepatocytes          EGF 7.523577e-08 0.000
# 14162          Treg Hepatocytes        FASLG 1.109393e-08 0.030
