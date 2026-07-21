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

unique(all_data$sample)


all_data_Air <- all_data[, all_data@meta.data$group %in% c("AirMale", "AirFemale")]
all_data_Hyperoxia <- all_data[, all_data@meta.data$group %in% c("HyperoxiaMale", "HyperoxiaFemale")]
table(all_data_Air$group)
table(all_data_Hyperoxia$group)


all_data@meta.data$group_new <- all_data@meta.data$group

all_data@meta.data$group_new <- gsub("AirFemale", "Air", all_data@meta.data$group_new)
all_data@meta.data$group_new <- gsub("AirMale", "Air", all_data@meta.data$group_new)
all_data@meta.data$group_new <- gsub("HyperoxiaFemale", "Hyperoxia", all_data@meta.data$group_new)
all_data@meta.data$group_new <- gsub("HyperoxiaMale", "Hyperoxia", all_data@meta.data$group_new)

table(all_data@meta.data$group_new)




All_type <- as.character(unique(all_data$annotation_full))
print(All_type)

All_type %>%
  setNames(.,rep('all_data', 9)) %>%
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
          ident.1 = 'Hyperoxia',
          ident.2 = 'Air',
          test.use = "wilcox",
          logfc.threshold = 0.25,
          min.pct = 0.1,
          group.by = 'group_new'
        )
      # ==== 新增：检查空结果 ====
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
                   group.by = "group_new",
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
      error_dir <- "DGE_and_enrichment_250730/all_data_error_objects/"
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







