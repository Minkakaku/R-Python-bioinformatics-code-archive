library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Mm.eg.db")
# library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 

# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)

library(RBioFormats)
source("xenium.R")


# 读取数据
load("data_for_analysis/all_data_annotation_rctd.rData")

# enrichment---all_data_Air------------------
DefaultAssay(all_data) <- "RNA"
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')
saveRDS(all_data, "data_for_analysis/all_data_annotaion_rctd_scale.rds")

all_data_Air <- all_data[, all_data@meta.data$group %in% c("AirMale", "AirFemale")]
table(all_data_Air$group)

All_type <- as.character(unique(all_data_Air$annotation_full))
print(All_type)

All_type %>%
  setNames(.,rep('all_data_Air', 9)) %>%
  imap(function(celltype,name){ #browser()
    tryCatch({
      object <- get(name)
      Idents(object)<-'annotation_full'
      # print(celltype)
      # print(name)
      # return()
      # 验证亚型隔离
      cat("Processing:", celltype, "| Cells:", sum(Idents(object) == celltype), "\n")
      
      # 使用子集对象提高效率
      subset_obj <- subset(JoinLayers(object), idents = celltype)
      
      diffmarkers <-
        FindMarkers(
          subset_obj,
          ident.1 = 'AirMale',
          ident.2 = 'AirFemale',
          test.use = "wilcox",
          logfc.threshold = 0.25,
          min.pct = 0.1,
          group.by = 'group'
        )
      # ==== 检查空结果 ====
      if (nrow(diffmarkers) == 0) {
        cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
        write.csv(data.frame(), str_glue('DGE_and_enrichment_250730/{name}/{celltype}/NO_DEGs.csv'))
        return(NULL)
      }
      
      dir.create(file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/')), recursive = TRUE)
      write.csv(diffmarkers,file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} group diffmarker.csv')))
      
      enrichGO_ls <- enrichGO(
        rownames(diffmarkers),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_ls@result,file = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} Go enriment table.csv')),row.names = FALSE)
      
      enrichGO_dot <- enrichplot::dotplot(
        enrichGO_ls,
        showCategory = 10,
        label_format = 50)
      enrichGO_bar <- barplot(
        enrichGO_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_dot,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} Go dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichGO_bar,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} Go barplot.pdf')),
        height = 6,
        width = 8)
      
      de_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_ls <- enrichKEGG(
        de_genes_id_for_enrichmenet,
        # organism = "mmu", 
        keyType = "ncbi-geneid",
        keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        use_internal_data = FALSE
        # use_internal_data = TRUE
        )
      enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_ls@result,file = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} KEGG enriment table.csv')),row.names = FALSE)
      
      enrichKEGG_dot <- enrichplot::dotplot(
        enrichKEGG_ls,
        showCategory = 10,
        label_format = 50)
      enrichKEGG_bar <- barplot(
        enrichKEGG_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichKEGG_dot,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} KEGG dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichKEGG_bar,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} KEGG barplot.pdf')),
        height = 6,
        width = 8)
      
      top <- rbind(diffmarkers %>% top_n(n = 25, wt = avg_log2FC),diffmarkers %>% top_n(n = -25, wt = avg_log2FC))
      p<-DoHeatmap(object, label = F,
                   features = rownames(top),
                   group.by = "group",
                   disp.max= 1,
                   disp.min= -1,
                   assay = 'RNA')+ #设置组别颜色
        scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
      pdf(file.path(str_glue("DGE_and_enrichment_250730/{name}/{celltype}/diffgene_heatmap.pdf")), height = 8, width =10)
      print(p)
      dev.off()
      p_volcano <- EnhancedVolcano(diffmarkers,
                                   lab = rownames(diffmarkers),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   title = str_glue('{celltype} volcano plot'),
                                   pCutoff = 10e-2,
                                   FCcutoff = 0,
                                   pointSize = 3.0,
                                   labSize = 6.0)
      ggsave(p_volcano,
             filename = file.path(str_glue("DGE_and_enrichment_250730/{name}/{celltype}/diffgene valcano.pdf")),
             height = 6,
             width = 8)
    }, error = function(e) {
      # ===== 错误处理代码 ===== #
      message("\nERROR in ", celltype, ": ", e$message)
      # 创建专门存放错误对象的目录
      error_dir <- "DGE_and_enrichment_250730/all_data_Air_error_objects/"
      dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
      # 保存引发错误的数据对象
      saveRDS(
        list(celltype = celltype, 
             error_msg = e$message, 
             subset_obj = if(exists("subset_obj")) subset_obj else NULL),
        file.path(error_dir, paste0(celltype, "_error.rds"))
      )
      return(NULL)
    })
  })





# enrichment---all_data_Hyperoxia------------------

all_data_Hyperoxia <- all_data[, all_data@meta.data$group %in% c("HyperoxiaMale", "HyperoxiaFemale")]
table(all_data_Hyperoxia$group)


All_type <- as.character(unique(all_data_Hyperoxia$annotation_full))
print(All_type)

All_type %>%
  setNames(.,rep('all_data_Hyperoxia', 9)) %>%
  imap(function(celltype,name){ #browser()
    tryCatch({
      object <- get(name)
      Idents(object)<-'annotation_full'
      # print(celltype)
      # print(name)
      # return()
      # 验证亚型隔离
      cat("Processing:", celltype, "| Cells:", sum(Idents(object) == celltype), "\n")
      
      # 使用子集对象提高效率
      subset_obj <- subset(JoinLayers(object), idents = celltype)
      
      diffmarkers <-
        FindMarkers(
          subset_obj,
          ident.1 = 'HyperoxiaMale',
          ident.2 = 'HyperoxiaFemale',
          test.use = "wilcox",
          logfc.threshold = 0.25,
          min.pct = 0.1,
          group.by = 'group'
        )
      # ==== 检查空结果 ====
      if (nrow(diffmarkers) == 0) {
        cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
        write.csv(data.frame(), str_glue('DGE_and_enrichment_250730/{name}/{celltype}/NO_DEGs.csv'))
        return(NULL)
      }
      
      dir.create(file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/')), recursive = TRUE)
      write.csv(diffmarkers,file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} group diffmarker.csv')))
      
      enrichGO_ls <- enrichGO(
        rownames(diffmarkers),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_ls@result,file = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} Go enriment table.csv')),row.names = FALSE)
      
      enrichGO_dot <- enrichplot::dotplot(
        enrichGO_ls,
        showCategory = 10,
        label_format = 50)
      enrichGO_bar <- barplot(
        enrichGO_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_dot,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} Go dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichGO_bar,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} Go barplot.pdf')),
        height = 6,
        width = 8)
      
      de_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_ls <- enrichKEGG(
        de_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        use_internal_data = FALSE
        # use_internal_data = TRUE
        )
      enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_ls@result,file = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} KEGG enriment table.csv')),row.names = FALSE)
      
      enrichKEGG_dot <- enrichplot::dotplot(
        enrichKEGG_ls,
        showCategory = 10,
        label_format = 50)
      enrichKEGG_bar <- barplot(
        enrichKEGG_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichKEGG_dot,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} KEGG dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichKEGG_bar,
        filename = file.path(str_glue('DGE_and_enrichment_250730/{name}/{celltype}/{celltype} KEGG barplot.pdf')),
        height = 6,
        width = 8)
      
      top <- rbind(diffmarkers %>% top_n(n = 25, wt = avg_log2FC),diffmarkers %>% top_n(n = -25, wt = avg_log2FC))
      p<-DoHeatmap(object, label = F,
                   features = rownames(top),
                   group.by = "group",
                   disp.max= 1,
                   disp.min= -1,
                   assay = 'RNA')+ #设置组别颜色
        scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
      pdf(file.path(str_glue("DGE_and_enrichment_250730/{name}/{celltype}/diffgene_heatmap.pdf")), height = 8, width =10)
      print(p)
      dev.off()
      p_volcano <- EnhancedVolcano(diffmarkers,
                                   lab = rownames(diffmarkers),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   title = str_glue('{celltype} volcano plot'),
                                   pCutoff = 10e-2,
                                   FCcutoff = 0,
                                   pointSize = 3.0,
                                   labSize = 6.0)
      ggsave(p_volcano,
             filename = file.path(str_glue("DGE_and_enrichment_250730/{name}/{celltype}/diffgene valcano.pdf")),
             height = 6,
             width = 8)
    }, error = function(e) {
      # ===== 错误处理代码 ===== #
      message("\nERROR in ", celltype, ": ", e$message)
      # 创建专门存放错误对象的目录
      error_dir <- "DGE_and_enrichment_250730/all_data_Hyperoxia_error_objects/"
      dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
      # 保存引发错误的数据对象
      saveRDS(
        list(celltype = celltype, 
             error_msg = e$message, 
             subset_obj = if(exists("subset_obj")) subset_obj else NULL),
        file.path(error_dir, paste0(celltype, "_error.rds"))
      )
      return(NULL)
    })
  })




# 单独绘制AT2
Idents(all_data_Hyperoxia)<-'annotation_full'

AT2 <- subset(JoinLayers(all_data_Hyperoxia), idents = "AT2")

diff_markers <- FindMarkers(AT2, 
                            ident.1 = 'HyperoxiaMale', 
                            ident.2 = 'HyperoxiaFemale', 
                            test.use = "wilcox",
                            logfc.threshold = 0.25,
                            min.pct = 0.1,
                            group.by = 'group'
                            )

nrow(diff_markers)

write.csv(diff_markers, file = "all_data_hyperoxia AT2 diffmarker.csv")


AT2_enrichGO <- enrichGO(rownames(diff_markers),
                     OrgDb = org.Mm.eg.db, 
                     ont = "ALL", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05, 
                     qvalueCutoff = 0.05,
                     keyType = 'SYMBOL'
                     )
write.csv(AT2_enrichGO@result, file = "all_data_hyperoxia AT2 GO enrichment table.csv")


AT2_enrichGO_dot <- enrichplot::dotplot(AT2_enrichGO, 
                                        showCategory = 10,
                                        label_format = 50
                                        )
ggsave(AT2_enrichGO_dot, filename = "all_data_hyperoxia AT2 Go dotplot.pdf", width = 8, height = 6)


AT2_enrichGO_bar <- barplot(AT2_enrichGO, 
                            showCategory = 10,
                            label_format = 50
                            )
ggsave(AT2_enrichGO_bar, filename = "all_data_hyperoxia AT2 Go barplot.pdf", width = 8, height = 6)


de_genes_id_for_enrichmenet <- clusterProfiler::bitr(rownames(diff_markers), 
                                                     OrgDb = "org.Mm.eg.db", 
                                                     fromType = "SYMBOL", 
                                                     toType = "ENTREZID"
                                                     ) %>% .$ENTREZID
AT2_enrichKEGG <- enrichKEGG(de_genes_id_for_enrichmenet, 
                             organism = "mmu", 
                             keyType = "ncbi-geneid",
                             # keyType = "kegg", 
                             pvalueCutoff = 0.05, 
                             pAdjustMethod = "BH", 
                             minGSSize = 10,
                             maxGSSize = 500, 
                             qvalueCutoff = 0.2,
                             use_internal_data = FALSE
                             # use_internal_data = TRUE
                             )
AT2_enrichKEGG <- setReadable(AT2_enrichKEGG,  OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(AT2_enrichKEGG@result, file = "all_data_hyperoxia AT2 KEGG enrichment table.csv")



AT2_enrichKEGG_dot <- enrichplot::dotplot(AT2_enrichKEGG,
                                          showCategory = 10, 
                                          label_format = 50
                                          )
ggsave(AT2_enrichKEGG_dot, filename = "all_data_hyperoxia AT2 KEGG dotplot.pdf", width = 8, height = 6)


AT2_enrichKEGG_bar <- barplot(AT2_enrichKEGG, 
                              showCategory = 10,
                              label_format = 50
                              )
ggsave(AT2_enrichKEGG_bar, filename = "all_data_hyperoxia AT2 KEGG barplot.pdf", width = 8, height = 6)


top <- rbind(diff_markers %>% top_n(n = 25, wt = avg_log2FC), diff_markers %>% top_n(n = -25, wt = avg_log2FC))
p <- DoHeatmap(all_data_Hyperoxia, 
               label = F, 
               features = rownames(top),
               group.by = "group",
               disp.max= 1,
               disp.min= -1,
               assay = 'RNA'
               ) + 
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))

pdf("all_data_hyperoxia AT2 diffgene_heatmap.pdf", height = 8, width =10)
print(p)
dev.off()


p_volcano <- EnhancedVolcano(diff_markers, 
                             lab = rownames(diff_markers), 
                             x = 'avg_log2FC', 
                             y = 'p_val_adj', 
                             title = "AT2 volcano plot", 
                             pCutoff = 10e-2, 
                             FCcutoff = 0, 
                             pointSize = 3.0, 
                             labSize = 6.0
                             )
ggsave(p_volcano, filename = "AT2 diffgene valcano.pdf", height = 6, width = 8)









