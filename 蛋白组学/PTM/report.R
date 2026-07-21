library(dplyr)
library(data.table)

##################### parameter #########################################

modified.type <- snakemake@params[["modified_type"]]
if(modified.type=='Phosphorylation (STY)') motif.width=6 else motif.width=10
sample.file <- snakemake@params[["sample_file"]]
adj <- snakemake@params[['diff_adj']]
fdr.value <- snakemake@params[['diff_fdr']]
fc.value <- snakemake@params[['diff_fc']]
db.dir <- snakemake@params[["db_dir"]]
source_dir<- snakemake@params[["sourc_dir"]]
dirname <- snakemake@params[["out_dir"]]

project_id = snakemake@params[['project_id']]
project_name = snakemake@params[['project_name']]
project_institute = snakemake@params[['project_institute']]
customer_name = snakemake@params[['customer_name']]

project.info=data.frame(
  ID=c("项目名称: ","项目编号: ","客户单位: ","客户姓名: ","报告日期: "),
  value=c(project_name,project_id,project_institute,customer_name,format(Sys.Date(),'%Y-%m-%d'))
)

######################### create report ##########################

## load data
db.infor=fread(paste(db.dir,"summary.tsv",sep="/"))
db.name=db.infor$`Database name`
db.time=sub(db.infor$`Create time`,pattern = " .*",replacement = "")
db.count=db.infor$Proteins
sample=fread(sample.file)
files=dir(paste(dirname,"2-Differentially_expressed_protein",sep="/"),pattern = "Table 2. DEPs list of ",full.names = T)
file.list=setNames(files,sub(sub(files,pattern=".*Table 2. DEPs list of ",replacement=""),pattern=".xls$",replacement=""))
diff.data=diff.table=rbindlist(lapply(file.list,function(x){
  subset(fread(x),Regulation!='unchanged',select=c("Protein accession","Amino acid","Position","Gene name","Protein description","log2FC","adj.P.Val","Regulation"))
}),idcol = "Compare")
diff.table=rbind(
  subset(diff.data,select=c("Regulation","Compare")) %>% table() %>% as.data.table() %>% dcast(Compare~Regulation,value.var = "N",fill = 0) %>% mutate(Type="Site",.after="Compare"),
  subset(diff.data %>% distinct(`Protein accession`,Regulation,Compare),select=c("Regulation","Compare")) %>% table() %>% as.data.table() %>% dcast(Compare~Regulation,value.var = "N",fill = 0) %>% mutate(Type="Protein",.after="Compare")
) %>% arrange(Compare)

summary.table=fread(paste(dirname,"1-MS_identified_summary/Table 1b. MS identified modification summary.xls",sep="/"))

##  copy file
file.copy(from = file.path(source_dir,'report_dia_ptm.Rmd'),to = paste(dirname,"report.Rmd",sep="/"),overwrite = T)
file.copy(from = file.path(source_dir,'Astral.png'),to = paste(dirname,"Astral.png",sep="/"),overwrite = T)
file.copy(from = file.path(source_dir,'DIA-NN workflow.png'),to = paste(dirname,"DIA-NN workflow.png",sep="/"),overwrite = T)
file.copy(from = file.path(source_dir,'Astral DIA PTM.png'),to = paste(dirname,"Astral DIA.png",sep="/"),overwrite = T)
file.copy(from = paste(source_dir,'logo.png',sep="/"),to = paste(dirname,"logo.png",sep="/"),overwrite = T)
file.copy(from = paste(source_dir,'PTM Data process.png',sep="/"),to = paste(dirname,"Data process.png",sep="/"),overwrite = T)
file.copy(from = paste(source_dir,'Astral_PTM.pdf',sep="/"),to = paste(dirname,"报告说明文档.pdf",sep="/"),overwrite = T)
report=paste(dirname,"report.html",sep="/")

##  render html report
rmarkdown::render(input = paste(dirname,"report.Rmd",sep="/"),output_format = 'html_document',output_file = report)

##  remove file
unlink(paste(dirname,"report.Rmd",sep="/"))
unlink(dir(dirname,pattern = "*.png",full.names = T))
