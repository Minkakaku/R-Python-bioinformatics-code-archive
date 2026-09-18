library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
# library(pheatmap)
# library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyr)
library(stringr)
# library(GSVA)
library(patchwork)
# library(limma)
# library(clusterProfiler)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpubr)

# genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/script/new_geneset.csv",sep = "\t",head=TRUE)
genesets <- list(
  ISC_identity = c("Lgr5", "Olfm4", "Ascl2", "Smoc2", "Axin2", "Sox9"),
  TFRC_iron_metabolism = c("Tfrc", "Slc11a2", "Slc40a1", "Fth1", "Ftl1", "Hmox1", "Ncoa4", "Steap3", "Ireb2", "Aco1"),
  ferroptosis = c("Gpx4", "Slc7a11", "Slc3a2", "Acsl4", "Lpcat3", "Alox15", "Alox5", "Sat1", "Chac1", "Nfe2l2", "Hmox1"),
  proliferation = c("Mki67", "Top2a", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7", "Ccna2", "Ccnb1"),
  DNA_damage = c("Trp53", "Cdkn1a", "Gadd45a", "Bax", "Bbc3", "Mdm2", "Ddit4"),
  oxidative_stress = c("Nqo1", "Gclm", "Gclc", "Hmox1", "Txnrd1", "Srxn1", "Prdx1", "Sod1", "Sod2"),
  apoptosis = c("Bax", "Bbc3", "Pmaip1", "Casp3", "Casp7"),
  senescence = c("Cdkn1a", "Cdkn2a", "Cdkn2b", "Serpine1", "Igfbp7", "Glb1"),
  regeneration = c("Clu", "Anxa1", "Tacstd2", "Ly6a", "Krt19", "Yap1", "Tead1", "Areg")
)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/02.Mission/mapping_raw/Main_mapping.rds")
# adata<-subset(adata,cell_type %in% c("Dendritic cell","Classical monocyte","Non-classical monocyte"))
# gene<-c(genesets$GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION)
# gene <- gene[gene != ""]
# gene<-list(gene)

adata<-AddModuleScore(
    object=adata,
    features=genesets,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[15:23]<-names(genesets)
# colnames(adata@meta.data)[15]<-'GOBP_CELL_ACTIVATION'
write.table(adata@meta.data,"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/geneset_score/all_cells_features.xls",sep="\t")
# Idents(adata)<-adata$annotation
# p <- VlnPlot(adata,features='GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',pt.size = 0)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all.VlnPlot.pdf",width = 10,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all.VlnPlot.png",width = 10,height = 10)
outdir<-"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/geneset_score_cell_type/"
adata$cell_type_group<-paste0(adata$cell_type,adata$group)
adata<-subset(adata,cell_type_group %in% 'B cellA2',invert=TRUE)
for(i in names(genesets)){
	p <- VlnPlot(adata,features=i,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(i,' ~ split')),group.by='ident',prefix=i,map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))
get_sig_ana(p$data,as.formula(paste0(i,' ~ split')),group.by='ident',prefix=paste0(i,'_cell_type'),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.pdf"),width = 21,height = 10)
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.png"),width = 21,height = 10)
}

outdir<-"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/geneset_score_sub_cell_type/"
Idents(adata)<-adata$sub_cell_type
for(i in names(genesets)){
	p <- VlnPlot(adata,features=i,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave(paste0(outdir,i,"sub_cell_type_group.VlnPlot.pdf"),width = 25,height = 10)
ggsave(paste0(outdir,i,"sub_cell_type_group.VlnPlot.png"),width = 25,height = 10)

p <- get_sig_ana(p,as.formula(paste0(i,' ~ split')),group.by='ident',prefix=i,map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))
get_sig_ana(p$data,as.formula(paste0(i,' ~ split')),group.by='ident',prefix=paste0(i,'_sub_cell_type'),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,i,"sub_cell_type_group_sig_ana.VlnPlot.pdf"),width = 36,height = 10)
ggsave(paste0(outdir,i,"sub_cell_type_group_sig_ana.VlnPlot.png"),width = 36,height = 10)
}


#### 线粒体健康评分
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/02.Mission/mapping_raw/Main_mapping.rds")
adata<-subset(adata,sub_cell_type %in% 'Paneth')
mito_genes <- grep("^mt-", rownames(adata), value = TRUE)
genesets <- list(
  Mitochondrial_health_score = mito_genes
)

adata<-AddModuleScore(
    object=adata,
    features=genesets,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[15]<-names(genesets)

write.table(adata@meta.data,"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/Mitochondrial_health_score/all_cells_features.xls",sep="\t")

outdir<-"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/Mitochondrial_health_score/"
for(i in names(genesets)){
	p <- VlnPlot(adata,features=i,group.by='group',pt.size = 0)+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.png"),width = 7,height = 10)

p <- get_sig_ana(p,as.formula(paste0(i,' ~ ident')),group.by='ident',prefix=i,map_sig_level=T,only_sig = F,reset_y=TRUE,method = "wilcox.test")+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))
get_sig_ana(p$data,as.formula(paste0(i,' ~ ident')),group.by='ident',prefix=paste0(i,'_cell_type'),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.png"),width = 7,height = 10)
}


#### B细胞 NF-kappaB signaling pathway
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/02.Mission/mapping_raw/Main_mapping.rds")
adata<-subset(adata,cell_type %in% 'B cell')
genesets <- list(
  'NF-kappa_B_signaling_pathway' = c('Lck','Zap70','Lat','Plcg1','Prkcq','Syk','Lyn','Blnk','Btk','Plcg2','Prkcb','Card10','Card11','Card14','Bcl10','Malt1','Il1b','Il1r1','Myd88','Irak1','Irak4','Traf6','Tnf','Tnfrsf1a','Ripk1','Tradd','Traf2','Traf5','Birc3','Birc2','Eda','Edar','Edaradd','Cyld','Eda2r','Rigi','Trim25','Trim80','Lbp','Cd14','Tlr4','Ly96','Tirap','Ticam2','Ticam1','Cd40lg','Cd40','Traf3','Tnfsf11','Tnfrsf11a','Lta','Ltb','Tnfsf14','Ltbr','Map3k14','Map3k7','Tab1','Tab2','Tab3','Tnfsf13b','Tnfrsf13c','Ikbkg','Chuk','Ikbkb','Parp1','Pias4','Ube2i','Atm','Pidd1','Erc1','Nfkbia','Nfkb1','Rela','Ppp1r13l','Cflar','Xiap','Bcl2l1','Bcl2','Gadd45a','Gadd45b','Gadd45g','Traf1','Bcl2a1a','Bcl2a1b','Bcl2a1c','Bcl2a1d','Nfkb2','Tnfaip3','Ptgs2','Ccl4','Cxcl1','Cxcl2','Cxcl3','Vcam1','Plau','Csnk2a1','Csnk2a2','Csnk2b','Relb','Ccl19','Ccl21b','Ccl21e','Ccl21f','Ccl21d','Ccl21a','Cxcl12','Icam1')
)

adata<-AddModuleScore(
    object=adata,
    features=genesets,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[15]<-names(genesets)

write.table(adata@meta.data,"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/NF-kappa B/all_cells_features.xls",sep="\t")

outdir<-"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/NF-kappa B/"
for(i in names(genesets)){
	p <- VlnPlot(adata,features=i,group.by='group',pt.size = 0)+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.png"),width = 7,height = 10)

p <- get_sig_ana(p,as.formula(paste0(i,' ~ ident')),group.by='ident',prefix=i,map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+ylim(-0.1,0.3)
get_sig_ana(p$data,as.formula(paste0(i,' ~ ident')),group.by='ident',prefix=paste0(i,'_cell_type'),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.png"),width = 7,height = 10)
}



#### 细胞增殖、凋亡
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/02.Mission/mapping_raw/Main_mapping.rds")
genesets <- list(
  Cell_cycle = c('Ccnd1','Ccnd2','Ccnd3','Cdk4','Cdk6','Rb1','Rbl1','Rbl2','Abl1','Hdac2','Hdac1','E2f1','E2f2','E2f3','E2f4','E2f5','Tfdp1','Tfdp2','Gsk3b','Tgfb1','Tgfb2','Tgfb3','Smad2','Smad3','Smad4','Myc','Zbtb17','Cdkn2a','Cdkn2b','Cdkn2c','Cdkn2d','Cdkn1b','Cdkn1c','Cdkn1a','Ccne1','Ccne2','Cdk2','Skp1','Cul1','Rbx1','Skp2','Ccna1','Ccna2','Cdc6','Cdt1','Orc1','Orc2','Orc3','Orc4','Orc5','Orc6','Mcm2','Mcm3','Mcm4','Mcm5','Mcm6','Mcm7','Cdc45','Ticrr','Mtbp','Cdc7','Dbf4','Cdk1','Ccnb1','Ccnb1-ps','Ccnb2','Ccnb3','Cdc25b','Cdc25c','Gm2423','Ywhaq','Ywhaz','Ywhab','Ywhae','Ywhag','Ywhah','Plk1','Wee1','Wee2','Pkmyt1','Ccnh','Cdk7','Anapc1','Anapc2','Cdc27','Anapc4','Anapc5','Cdc16','Anapc7','Cdc23','Anapc10','Anapc11','Cdc26','Anapc13','Anapc15','Anapc16','Cdc20','Pttg1','Espl1','Smc1b','Smc1a','Smc3','Rad21','Stag1','Stag2','Atrx','Ddx11','Nipbl','Mau2','Hdac8','Esco2','Esco1','Pds5b','Pds5a','Wapl','Cdca5','Aurkb','Sgo1','Ppp2ca','Ppp2cb','Ppp2r1a','Ppp2r1b','Ppp2r5d','Ppp2r5b','Ppp2r5a','Ppp2r5c','Ppp2r5e','Ndc80','Knl1','Ttk','Mad1l1','Mad2l1','Mad2l2','Bub1b','Bub3','Bub1','Trip13','Mad2l1bp','Fzr1','Fbxo5','Cdc14b','Cdc14a','Crebbp','Ep300','Prkdc','Atr','Atm','Mdm2','Trp53','Gadd45a','Gadd45b','Gadd45g','Pcna','Sfn','Chek1','Chek2','Cdc25a'),
  DNA_replication = c('Ssbp1','Rnaseh1','Rpa1','Pcna','Dna2','Fen1','Lig1','Pola1','Pola2','Prim1','Prim2','Pold1','Pold2','Pold3','Pold4','Pole','Pole2','Pole3','Pole4','Mcm2','Mcm3','Mcm4','Mcm5','Mcm6','Mcm7','Rpa2','Gm6195','Rpa3','Rfc1','Rfc4','Rfc2','Rfc3','Rfc5','Rnaseh2a','Rnaseh2b','Rnaseh2c'),
  Apoptosis = c('Tnfsf10','Tnfrsf10b','Fasl','Fas','Fadd','Tnf','Tnfrsf1a','Tradd','Cflar','Nol3','Casp8','Casp6','Casp3','Casp7','Bid','Bax','Bak1','Diablo','Septin4','Htra2','Cycs','Cyct','Apaf1','Casp9','Prf1','Gzmb','Tuba1a','Tuba1b','Tuba3a','Tuba4a','Tuba1c','Tuba3b','Tubal3','Tuba8','Mcl1','Actb','Actg1','Sptan1','Lmna','Lmnb1','Lmnb2','Parp1','Parp2','Parp3','Parp4','Dffa','Dffb','Endog','Aifm1','Ern1','Traf2','Itpr1','Itpr2','Itpr3','Capn1','Capn2','Casp12','Eif2ak3','Eif2s1','Atf4','Ddit3','Ctsb','Ctsc','Ctsd','Ctsf','Ctsh','Ctsk','Ctsl','Ctso','Ctss','Ctsw','Ctsz','Birc3','Birc2','Xiap','Birc5','Bcl2l11','Bcl2l1','Bcl2','Daxx','Ripk1','Dab2ip','Map3k5','Mapk10','Mapk8','Mapk9','Bad','Jun','Fos','Trp53','Hrk','Map3k14','Chuk','Ikbkb','Ikbkg','Nfkbia','Nfkb1','Rela','Ptpn13','Gadd45a','Gadd45b','Gadd45g','Traf1','Bcl2a1a','Bcl2a1b','Bcl2a1c','Bcl2a1d','Atm','Pidd1','Bbc3','Pmaip1','Casp2','Ngf','Ntrk1','Il3','Il3ra','Csf2rb','Csf2rb2','Pik3ca','Pik3cd','Pik3cb','Pik3r1','Pik3r2','Pik3r3','Pdpk1','Akt1','Akt2','Akt3','Hras','Kras','Nras','Raf1','Map2k1','Map2k2','Mapk1','Mapk3')
)

adata<-AddModuleScore(
    object=adata,
    features=genesets,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[15:17]<-names(genesets)

outdir<-"/mnt/NexData_file/home/sh/Projects/SC2604012626-01/03.Mission/proliferation_apoptosis/"
write.table(adata@meta.data,paste0(outdir,"all_cells_features.xls"),sep="\t")

for(i in names(genesets)){
	p <- VlnPlot(adata,features=i,group.by='group',pt.size = 0)+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,i,"cell_type_group.VlnPlot.png"),width = 7,height = 10)

ylim1<-min(p$data[i])-0.1
ylim2<-max(p$data[i])+0.1
p <- get_sig_ana(p,as.formula(paste0(i,' ~ ident')),group.by='ident',prefix=i,map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(A1='#1f77b4',A2='#ff7f0e',A3='#279e68'))+ylim(ylim1,ylim2)
get_sig_ana(p$data,as.formula(paste0(i,' ~ ident')),group.by='ident',prefix=paste0(i,'_cell_type'),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,i,"cell_type_group_sig_ana.VlnPlot.png"),width = 7,height = 10)
}
