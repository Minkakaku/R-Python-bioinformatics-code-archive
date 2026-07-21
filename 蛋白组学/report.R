# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(dplyr)
library(data.table)

##################### parameter #########################################
source_dir<- snakemake@params[["sourc_dir"]]

sample.file <- snakemake@params[["sample_file"]]
adj <- snakemake@params[['diff_adj']]
fdr.value <- snakemake@params[['diff_fdr']]
fc.value <- snakemake@params[['diff_fc']]

project_id = snakemake@params[['project_id']]
project_name = snakemake@params[['project_name']]
project_institute = snakemake@params[['project_institute']]
customer_name = snakemake@params[['customer_name']]

project.info=data.frame(
  ID=c("项目名称: ","项目编号: ","客户单位: ","客户姓名: "),
  value=c(project_name,project_id,project_institute,customer_name)
)

db.dir <- snakemake@params[["db_dir"]]
quant.type <- snakemake@params[["quant_type"]]
dirname <- snakemake@params[["out_dir"]]

######################### create report ##########################

## load data
db.infor=fread(paste(db.dir,"summary.tsv",sep="/"))
db.name=db.infor$`Database name`
db.time=sub(db.infor$`Create time`,pattern = " .*",replacement = "")
db.count=db.infor$Proteins
sample=fread(sample.file)
files=dir(paste(dirname,"2-Differentially_expressed_protein",sep="/"),pattern = "Table 2. DEPs list of ",full.names = T)

message(paste(dirname,"2-Differentially_expressed_protein",sep="/"))

file.list=setNames(files,sub(sub(files,pattern=".*Table 2. DEPs list of ",replacement=""),pattern=".xls$",replacement=""))
diff.data=diff.table=rbindlist(lapply(file.list,function(x){
  subset(fread(x),Regulation!='unchanged',select=c("Protein accession","Gene name","Protein description","log2FC","adj.P.Val","Regulation"))
}),idcol = "Compare")




diff.table=subset(diff.data,select=c("Regulation","Compare")) %>% table() %>% as.data.table() %>% dcast(Compare~Regulation,value.var = "N",fill = 0)
summary.table=fread(paste(dirname,"1-MS_identified_summary/Table 1c. MS identified summary.xls",sep="/"))

##  copy file
if(quant.type=='TMT'){
  file.copy(from = paste(source_dir,'report.Rmd',sep="/"),to = paste(dirname,"report.Rmd",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'TMT process.png',sep="/"),to = paste(dirname,"TMT process.png",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'TMT.pdf',sep="/"),to = paste(dirname,"Report description.pdf",sep="/"),overwrite = T)
}else if(quant.type=='DIA'){
  file.copy(from = paste(source_dir,'report_v2.Rmd',sep="/"),to = paste(dirname,"report.Rmd",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'Label-free process.png',sep="/"),to = paste(dirname,"TMT process.png",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'Astral.pdf',sep="/"),to = paste(dirname,"Report description.pdf",sep="/"),overwrite = T)
  
  file.copy(from = paste(source_dir,'Astral.png',sep="/"),to = paste(dirname,"Astral.png",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'Astral DIA.png',sep="/"),to = paste(dirname,"Astral DIA.png",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'DIA-NN workflow.png',sep="/"),to = paste(dirname,"DIA-NN workflow.png",sep="/"),overwrite = T)
}else if(quant.type=='PD'){
  file.copy(from = paste(source_dir,'report_v2.bak.Rmd',sep="/"),to = paste(dirname,"report.Rmd",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'Label-free process.png',sep="/"),to = paste(dirname,"TMT process.png",sep="/"),overwrite = T)
  file.copy(from = paste(source_dir,'TMT.pdf',sep="/"),to = paste(dirname,"Report description.pdf",sep="/"),overwrite = T)
}

file.copy(from = paste(source_dir,'logo.png',sep="/"),to = paste(dirname,"logo.png",sep="/"),overwrite = T)
file.copy(from = paste(source_dir,'Data process.png',sep="/"),to = paste(dirname,"Data process.png",sep="/"),overwrite = T)
report=paste(dirname,"report.html",sep="/")

##  render html report
rmarkdown::render(input = paste(dirname,"report.Rmd",sep="/"),output_format = 'html_document',output_file = report)

##  remove file
unlink(paste(dirname,"report.Rmd",sep="/"))
unlink(dir(dirname,pattern = "*.png",full.names = T))

# message("报告生成！")
# sink(type = "message")
# close(log_con)