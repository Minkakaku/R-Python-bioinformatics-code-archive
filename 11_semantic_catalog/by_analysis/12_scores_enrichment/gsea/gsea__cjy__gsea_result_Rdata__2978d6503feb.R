BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("yulab.utils")
library(clusterProfiler)
library(org.Hs.eg.db)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# 1. 读入你的未过滤总表（请确保文件名对应）
deg_total <- data.table::fread(file.choose())

# 2. 用第一列“Protein accession”去完美转换 ID（fromType 必须是 "UNIPROT"）
gene_ids <- bitr(deg_total$`Protein accession`, 
                 fromType = "UNIPROT", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# 3. 把转换后的数字 ID 融合成大表
deg_with_id <- merge(deg_total, gene_ids, by.x = "Protein accession", by.y = "UNIPROT")

# 4. 提取 log2FC 数据（这里以你的 H1 列为例，如果想跑 H2，把这里改成 H2 即可）
gene_list <- deg_with_id$H1

# 5. 用人类 KEGG 专用的 ENTREZID 编号给它命名
names(gene_list) <- deg_with_id$ENTREZID

# 6. GSEA 灵魂排序：从大到小
gene_list <- sort(gene_list, decreasing = TRUE)

# # 7. 运行 GSEA (自动联网下载人类 KEGG 数据库)
# gsea_res <- gseKEGG(geneList = gene_list, 
#                     organism = "hsa", 
#                     pvalueCutoff = 1) # 设为 1 可以保留全部通路结果，方便后续挑选
# 


kegg_anno_raw <- fread(file.choose())
# GSEA 的 TERM2GENE 要求第一列是通路名字，第二列是基因/蛋白ID
my_kegg_gmt <- data.frame(
  gs_name = kegg_anno_raw$Terms,               # 通路名字（如 Glycolysis / Gluconeogenesis）
  gene    = kegg_anno_raw$`Protein accession`  # 质谱蛋白号（如 Q9NUB1）
)
# 7.3 运行通用的 GSEA 函数（不需要联网，速度极快！）
gsea_res <- GSEA(geneList = gene_list, 
                 TERM2GENE = my_kegg_gmt, 
                 pvalueCutoff = 1)

# 9. 当场将其保存为你们苦苦寻找的 .Rdata 原始工程文件！
save(gsea_res, file = "my_gsea_result.Rdata")

