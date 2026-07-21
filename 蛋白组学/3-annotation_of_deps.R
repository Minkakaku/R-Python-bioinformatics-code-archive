# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

# 差异表达蛋白注释脚本
library(data.table)
library(dplyr)

# 读取注释文件
message("读取蛋白注释文件...")
data.go=fread("C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/hum/GO_annotation.tsv")
data.kegg=fread("C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/hum/KEGG_annotation.tsv")
data.pfam=fread("C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/hum/Pfam_annotation.tsv" )
data.nog=fread("C:/Users/Administrator/Desktop/生信标准工具/Database/Proteomics/annotation/hum/NOG_annotation.tsv")

# 读取差异表达蛋白数据
message("读取差异表达蛋白数据...")
diff.data <- fread("Table 2. DEPs list of Dox-vs-EBSS.xls")
diff.list=subset(diff.data, Regulation!='unchanged')
bp.diff=merge(subset(data.go,Category=='Biological process'),diff.list,by="Protein accession")
mf.diff=merge(subset(data.go,Category=='Molecular function'),diff.list,by="Protein accession")
cc.diff=merge(subset(data.go,Category=='Cellular component'),diff.list,by="Protein accession")
kegg.diff=merge(data.kegg,diff.list,by="Protein accession")
pfam.diff=merge(data.pfam,diff.list,by="Protein accession")
kog.diff=merge(data.nog,diff.list,by="Protein accession")

# 写出差异蛋白注释结果文件
message("保存注释结果...")
fwrite(bp.diff,file="bp_annotation.csv",sep="\t")
fwrite(cc.diff,file="cc_annotation.csv",sep="\t")
fwrite(mf.diff,file="mf_annotation.csv",sep="\t")
fwrite(kegg.diff,file="kegg_annotation.csv",sep="\t")
fwrite(pfam.diff,file="pfam_annotation.csv",sep="\t")
fwrite(kog.diff,file="kog_annotation.csv",sep="\t")


# message("差异表达蛋白注释完成！")
# sink(type = "message")
# close(log_con)
