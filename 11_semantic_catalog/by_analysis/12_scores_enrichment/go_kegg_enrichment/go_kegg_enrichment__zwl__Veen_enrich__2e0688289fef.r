setwd('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require2/Microglia_common')

E_IR_Ctrl = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require1/celltype_group_DEG/Microglia/IR_E_vs_Ctrl_E/all_cluster_marker.xls',sep='\t')
E_IRD_IR = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require1/celltype_group_DEG/Microglia/IR_E_vs_IRD_E/all_cluster_marker.xls',sep='\t')
E_NGFE_Ctrl = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require1/celltype_group_DEG/Microglia/Ctrl_E_vs_NGF_E/all_cluster_marker.xls',sep='\t')

library(ggvenn)
A = list('IR vs Ctrl' = E_IR_Ctrl$Symbol, 'IRD vs IR' = E_IRD_IR$Symbol, 'NGFE vs Ctrl' = E_NGFE_Ctrl$Symbol)
p = ggvenn(A,fill_color= c("#E41A1C", "#1E90FF", "#FF8C00"),fill_alpha = .7,#stroke_linetype = "longdash",
           set_name_size = 6,text_size=4)+theme(plot.background = element_rect('white'))
ggsave("1.Early_Time_commen_gene_Veen.png", width = 6, height = 6, dpi = 300,p)
ggsave("1.Early_Time_commen_gene_Veen.pdf", width = 6, height = 6, dpi = 300,p)

L_IR_Ctrl = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require1/celltype_group_DEG/Microglia/IR_L_vs_Ctrl_L/all_cluster_marker.xls',sep='\t')
L_IRD_IR = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require1/celltype_group_DEG/Microglia/IR_L_vs_IRD_L/all_cluster_marker.xls',sep='\t')
L_NGFE_Ctrl = read.csv('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require1/celltype_group_DEG/Microglia/NGF_L_vs_Ctrl_L/all_cluster_marker.xls',sep='\t')

A = list('IR vs Ctrl' = L_IR_Ctrl$Symbol, 'IRD vs IR' = L_IRD_IR$Symbol, 'NGFE vs Ctrl' = L_NGFE_Ctrl$Symbol)
p = ggvenn(A,fill_color= c("#E41A1C", "#1E90FF", "#FF8C00"),fill_alpha = .7,#stroke_linetype = "longdash",
           set_name_size = 6,text_size=4)+theme(plot.background = element_rect('white'))
ggsave("2.Late_Time_commen_gene_Veen.png", width = 6, height = 6, dpi = 300,p)
ggsave("2.Late_Time_commen_gene_Veen.pdf", width = 6, height = 6, dpi = 300,p)

Early = read.csv('./1.Early_Therapy_gene.csv')
Late = read.csv('./2.Late_Therapy_gene.csv')
library(clusterProfiler)
Earlyg = bitr(Early$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Lateg = bitr(Late$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

for (i in c('Early','Late')) {
  gene = get(paste0(i,'g'))
  
  GO_BP = enrichGO(gene          = gene$ENTREZID,
                   OrgDb         = 'org.Mm.eg.db',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

  
  GO_MF = enrichGO(gene          = gene$ENTREZID,
                   OrgDb         = 'org.Mm.eg.db',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

  
  GO_CC = enrichGO(gene          = gene$ENTREZID,
                   OrgDb         = 'org.Mm.eg.db',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

  
  KEGG <- enrichKEGG(gene         = gene$ENTREZID,
                   organism     = 'mus',
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)
  
  
  for(j in c('GO_BP','GO_CC','GO_MF','KEGG')){
    enrich = get(j)
    if(dim(enrich@result)[1]){
      enrich <- setReadable(enrich, OrgDb = 'org.Mm.eg.db', keyType="ENTREZID")
      write.csv(file = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require2/Microglia_common/3_',i,'_',j,'_enichment.csv'),enrich@result)
      A = barplot(enrich, showCategory=10) 
      ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require2/Microglia_common/3_',i,'_',j,'_enichment.png'),A,w=6,h=4)
      ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require2/Microglia_common/3_',i,'_',j,'_enichment.pdf'),A,w=6,h=4)
    }
  }
}

