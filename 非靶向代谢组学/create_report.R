library(data.table)
library(readxl)
library(dplyr)
library(shiny)
library(shinyFiles)
library(slickR)
library(rmarkdown)

############################## Parameters ######################################
# workflow = "Metabolomics" # Must be one of "Metabolomics", "Lipidomics" 
# results_dir = "C:/Users/css/Desktop/Rstudio/Rshiny/Meta_scripts/Results"
# TICs = TRUE #是否附有TIC图
# sample_info = "input/FWMS20252134/samples.csv"
# p.method = "none" # 是否进行pvalue校正 "BH","BY","holm"
# p.value = 0.05 # 差异阈值
# logfc = 0 # 差异阈值
# VIP = 1 # 差异阈值
# project.info=data.frame(
#   ID=c("Study name: ","Institute: ","Project ID: ","Report date: "),
#   value=c("Quantitative proteomics study of TMT labeled","CWMDA","PJ0001",format(Sys.Date(),'%Y-%m-%d'))
# )

project_id = snakemake@params[['project_id']]
project_name = snakemake@params[['project_name']]
project_institute = snakemake@params[['project_institute']]
customer_name = snakemake@params[['customer_name']]

project.info=data.frame(
  ID=c("项目名称: ","项目编号: ","客户单位: ","客户姓名: "),
  value=c(project_name,project_id,project_institute,customer_name)
)

p.method = snakemake@params[["p_method"]]
p.value = snakemake@params[["p_value"]]
logfc = snakemake@params[["logfc"]]
VIP = snakemake@params[["VIP"]] 
diff.method = snakemake@params[["diff_method"]]

workflow = snakemake@params[["workflow"]]
data.final.table.file = snakemake@input[["data_final_table"]]
sample_info = snakemake@input[["sample_info"]]
results_dir = snakemake@params[["results_dir"]]


############################ Create Report ####################################

## Load data 
sample = fread(sample_info)
if(all(sample$Sample == "-") == T){
  sample = sample %>% mutate(Sample = ID)
} else {
  sample = sample 
}

if(p.method=="none"){my.params=data.frame(p=paste0("p.value < ", p.value), fc=paste0("abs(log2fc) > ", logfc), v=paste0("VIP > ",VIP))}else{
  my.params=data.frame(p=paste0("adj.p.value < ", p.value), fc=paste0("abs(log2fc) > ", logfc), v=paste0("VIP > ",VIP))
}

if(workflow=="Metabolomics"){
  meta.data = fread(paste(results_dir,"1-MS_identified_summary/Table 1.1 MS identified metabolites detail_raw.xls",sep="/"))
  files = list.files(path = paste(results_dir,"3-Differentially abundant metabolites",sep = "/"), pattern = "Table 3a")
  diff.list = lapply(files, function(x) fread(paste(results_dir,"3-Differentially abundant metabolites",x,sep="/")))
  names(diff.list) = gsub("Table 3a. Differentially abundant metabolites list of ","",gsub(".xls","", files))
  diff.table=fread(paste(results_dir,"3-Differentially abundant metabolites/diff.stats.txt",sep = "/"))
}
if(workflow=="Lipidomics"){
  meta.data = read_excel(paste(results_dir,"1-MS_identified_summary/Table 1.1 MS identified lipids detail_raw.xlsx",sep="/"))
  files = list.files(path = paste(results_dir,"3-Differentially abundant lipids",sep = "/"), pattern = "Table 3a")
  diff.list = lapply(files, function(x) fread(paste(results_dir,"3-Differentially abundant lipids",x,sep="/")))
  names(diff.list) = gsub("Table 3a. Differentially abundant lipids list of ","",gsub(".xls","", files))
  diff.table=fread(paste(results_dir,"3-Differentially abundant lipids/diff.stats.txt",sep = "/"))
}

data.diff.table=rbindlist(lapply(diff.list,function(x) subset(x,select =c("Alignment ID","Metabolite name","log2fc","p.value","fdr","VIP","regulation"))),idcol = "Compare") %>% filter(regulation!='unchanged') %>% mutate(log2fc=round(log2fc,3),fdr=signif(fdr,3),p.value=signif(p.value,3),VIP=round(VIP,3))
file.copy(from = 'resources/logo.png',to = paste(results_dir,"logo.png",sep="/"),overwrite = T)

report=paste(results_dir,"report.html",sep="/")
if(workflow=="Metabolomics"){
  file.copy(from = 'resources/meta_dataprocess.png',to = paste(results_dir,"meta_dataprocess.png",sep="/"),overwrite = T)
  file.copy(from = "resources/meta_workflow.png", to = paste(results_dir, "meta_workflow.png", sep = "/"), overwrite = T)
  if(min(table(sample$Group))!= 1){
    file.copy(from = "resources/meta_report.Rmd", to = paste(results_dir, "meta_report.Rmd", sep = "/"), overwrite = T)
  }else{file.copy(from = "resources/meta_report_v2.Rmd", to = paste(results_dir, "meta_report.Rmd", sep = "/"), overwrite = T)}
  file.copy(from = 'resources/代谢组学报告解读文档.pdf', to = paste(results_dir,"代谢组学报告解读文档.pdf",sep ="/"), overwrite =T)
  rmarkdown::render(input = paste(results_dir,"meta_report.Rmd",sep="/"),output_format = 'html_document', output_file = report)
  unlink(paste(results_dir,"meta_report.Rmd",sep="/"))
  unlink(paste(results_dir,"meta_dataprocess.png",sep="/"))
  unlink(paste(results_dir,"meta_workflow.png",sep="/"))
}
if(workflow=="Lipidomics"){
  file.copy(from = 'resources/lipid_dataprocess.png',to = paste(results_dir,"lipid_dataprocess.png",sep="/"),overwrite = T)
  file.copy(from = "resources/lipid_workflow.png", to = paste(results_dir, "lipid_workflow.png", sep = "/"), overwrite = T)
  file.copy(from = "resources/lipid_report.Rmd", to = paste(results_dir, "lipid_report.Rmd", sep = "/"), overwrite = T)
  file.copy(from = 'resources/method.txt',to = paste(results_dir,"method.txt",sep="/"),overwrite = T)
  rmarkdown::render(input = paste(results_dir,"lipid_report.Rmd",sep="/"),output_format = 'html_document', output_file = report)
  unlink(paste(results_dir,"lipid_report.Rmd",sep="/"))
  unlink(paste(results_dir,"method.txt", sep="/"))
  unlink(paste(results_dir,"lipid_dataprocess.png",sep="/"))
  unlink(paste(results_dir,"lipid_workflow.png",sep="/"))
}
unlink(paste(results_dir,"logo.png",sep="/"))


# temp_report=paste(gsub(results_dir,pattern = ".*/",replacement =""),"report.html",sep="_")
# file.copy(from = report,to = paste("www",temp_report,sep="/"),overwrite = T)
