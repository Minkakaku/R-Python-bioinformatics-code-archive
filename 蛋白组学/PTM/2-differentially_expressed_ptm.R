library(limma)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)

# 使用limma包进行差异表达分析
ptm.limma.filter.s=function(data.quant,s,adj="BH"){
  if(min(table(s$Type))==1){
    s0=subset(s,Type==names(which.min(table(s$Type))))
    data.diff=as.data.frame(subset(data.quant,select = c("Protein accession","Amino acid","Position","Gene name","Protein description","Intensity",as.character(s$ID)))) %>%
      melt(id.vars=c("Protein accession","Amino acid","Position","Gene name","Protein description","Intensity",as.character(s0$ID))) %>% dplyr::rename(control.value=as.character(s0$ID)) %>%
      mutate(value=value-control.value) %>% dcast(`Protein accession`+`Amino acid`+Position+`Gene name`+`Protein description`+Intensity~variable,value.var = "value")
    s=subset(s,Type==names(which.max(table(s$Type))))
    MA=subset(data.diff,select=as.character(s$ID))
    MA=subset(MA,rowSums(is.na(MA))<nrow(s)/2)
    fit <- lmFit(MA)
    fit2 <- eBayes(fit)
    out=topTable(fit2,n=Inf,adjust.method = adj) %>% dplyr:rename(log2FC = "logFC")
    data.diff=merge(data.diff,subset(out,select=c("log2FC","P.Value","adj.P.Val")),by="row.names")
    data.diff$log2FC=round(data.diff$log2FC,3)
    data.diff$P.Value=signif(data.diff$P.Value,3)
    data.diff$adj.P.Val=signif(data.diff$adj.P.Val,3)
    data.diff=data.diff[-1]
    return(list(data=data.table(data.diff),sample=s))
  }else{
    data.diff=as.data.frame(subset(data.quant,select = c("Protein accession","Amino acid","Position","Gene name","Protein description","Intensity",s$ID)))   
    MA=subset(data.diff,rowSums(is.na(data.diff))<nrow(s)/2,select=as.character(s$ID))
    design <- model.matrix(~ 0+s$Type)
    compare=levels(s$Type)
    colnames(design) <- compare
    fit <- lmFit(MA, design)
    contrast.matrix <- makeContrasts(contrasts=paste(compare[1],compare[2],sep="-"), levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    out=topTable(fit2,n=Inf,adjust.method = adj) %>% dplyr::rename(log2FC = "logFC")
    data.diff=merge(data.diff,subset(out,select=c("log2FC","P.Value","adj.P.Val")),by="row.names")
    data.diff$log2FC=round(data.diff$log2FC,3)
    data.diff$P.Value=signif(data.diff$P.Value,3)
    data.diff$adj.P.Val=signif(data.diff$adj.P.Val,3)
    data.diff=data.diff[-1]
    return(list(data=data.table(data.diff),sample=s))
  }
}

# 使用significance_B方法进行差异表达分析
significance_filter=function(data.quant,s,adj="BH"){
  data.diff=as.data.frame(subset(data.quant,select = c("Protein accession","Gene name","Position","Amino acid","Protein description","Intensity",s$ID)))
  compare=levels(s$Type)
  data.diff.long=melt(data.diff,id.vars=c("Protein accession","Gene name","Position","Amino acid","Protein description","Intensity"),na.rm = T)%>%
    filter(Intensity >0) %>% merge(s,by.x="variable",by.y="ID") %>% group_by(`Protein accession`,Position)%>%filter(Type %in% compare & unique(length(Type)) > 1) %>%
    mutate(log2FC = mean(value[Type==compare[1]])-mean(value[Type==compare[2]])) %>% ungroup() %>% distinct(`Protein accession`,Position,log2FC,Intensity) %>% 
    mutate(P.Value=significance_B(log2FC,Intensity)) %>% distinct(`Protein accession`,Position,log2FC,P.Value)
  data.diff.long$adj.P.Val=p.adjust(data.diff.long$P.Value,method = adj)
  data.diff=merge(data.diff,data.diff.long,by=c("Protein accession","Position"))
  data.diff$log2FC=round(data.diff$log2FC,3)
  return(list(data=data.table(data.diff),sample=s))
}

# 根据阈值筛选差异表达蛋白
diff.filter=function(data,fdr=0.01,fc=1.5){
  data.diff=data  
  data.diff$Regulation="unchanged"
  data.diff$Regulation[data.diff$adj.P.Val<fdr&data.diff$log2FC>log2(fc)]='up-regulated'
  data.diff$Regulation[data.diff$adj.P.Val<fdr&data.diff$log2FC< -log2(fc)]='down-regulated'
  data.diff$`Gene name`=replace(data.diff$`Gene name`,data.diff$`Gene name`=='---',data.diff$`Protein accession`[data.diff$`Gene name`=='---'])
  return(data.diff)
}

# 计算significance_B值
significance_B=function(log2FC,intensity){
  # df=length(log2FC)/1000
  df = 5
  intensity=log2(intensity)
  a=(intensity-min(intensity)+1)*0.5
  t=quantile(log2FC,probs = c(0.5,0.8413))
  sapply(1:length(log2FC),function(x){
    z=abs((log2FC[x]-t[[1]])/(t[[2]]-t[[1]]))*a[x]
    2*pt(-z,df = df)
  },simplify = T)
}

# 生成火山图
volcano.plot=function(data.diff,fdr=0.01,fc=1.5){
  data.diff.s=merge(data.frame(Var1=c("down-regulated","unchanged","up-regulated"),color=c("steelblue","grey","red")),
                    as.data.frame(table(data.diff$Regulation)),by="Var1",all.x=T) %>% mutate(label=paste(Var1,Freq,sep=": "))
  xmax=quantile(abs(data.diff$log2FC),0.75)+5*IQR(abs(data.diff$log2FC))
  ymax=quantile(-log10(data.diff$adj.P.Val),0.75)+3*IQR(-log10(data.diff$adj.P.Val))
  p2=ggplot(data.diff,aes(x=log2FC,y=-log10(adj.P.Val),color=Regulation))+
    geom_point(size=1)+theme_classic(base_size = 14)+xlim(-xmax,xmax)+ylim(0,ymax)+
    scale_color_manual(values = setNames(data.diff.s$color,data.diff.s$Var1),labels=data.diff.s$label)+
    geom_hline(yintercept = -log10(fdr),linetype=2)+
    geom_vline(xintercept = c(-log2(fc),log2(fc)),linetype=2)
  return(p2)
}

# 生成差异表达蛋白热图
heatmap.plot <- function(data.diff, sample.list, sc="none", colv=TRUE, cn=TRUE, rn=FALSE, svalue=2, base_family='sans') {
  data.diff.filter <- subset(data.diff, Regulation != "unchanged")
  if(nrow(data.diff.filter) < 2) {
    return(NULL)
  }
  
  annot_sample <- data.frame(Type=sample.list$Type)
  rownames(annot_sample) <- sample.list$ID
  mycolor <- c("red", "blue")
  names(mycolor) <- levels(sample.list$Type)
  
  t.plot <- pheatmap(subset(data.diff.filter, select=rownames(annot_sample)), 
                     border_color = 'white', 
                     show_colnames = cn,
                     breaks = seq(-svalue, svalue, length.out=100), 
                     show_rownames = rn, 
                     cluster_cols = colv, 
                     silent = T, 
                     fontfamily=base_family,
                     annotation_col = annot_sample, 
                     cluster_rows = T, 
                     labels_row = data.diff.filter$`Gene name`, 
                     scale = sc,
                     color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100), 
                     annotation_colors = list(Type=mycolor))
  
  return(t.plot)
}

# 生成差异表达蛋白分析报告
generate_diff_report <- function(data.quant, sample.list, compare.group, output_dir, 
                                 diff_method="limma", paired=FALSE, adj="BH",
                                 fdr=0.01, fc=1.5, base_family='sans') {
  # 创建输出目录
  if(!dir.exists(output_dir)) {
    dir.create(path = output_dir, recursive = TRUE)
  }

  # 结果汇总表
  diff_summary <- data.frame()
  
  # 按比较组进行差异分析
  for(group_label in unique(compare.group$Label)) {
    # 提取当前比较组的样本
    s1 <- compare.group$Case[compare.group$Label == group_label][1]
    s2 <- compare.group$Control[compare.group$Label == group_label][1]
    current_samples <- sample.list[sample.list$Group %in% c(s1, s2), c("Sample", "Group")]
    colnames(current_samples) <- c("ID", "Type")
    current_samples$Type <- factor(current_samples$Type, levels = c(s1, s2))
    current_samples <- current_samples[order(current_samples$Type), ]

    # 差异分析
    if(diff_method == 'limma') {
      diff_result <- ptm.limma.filter.s(data.quant, current_samples, adj)
    } else {
      diff_result <- significance_filter(data.quant, current_samples, adj)
    }
    
    # 筛选差异表达蛋白
    diff_data <- diff.filter(diff_result$data, fdr=fdr, fc=fc)
    
    # 生成火山图
    
    volcano_plot <- volcano.plot(diff_data, fdr=fdr, fc=fc)
    
    # 生成热图
    heatmap_plot <- heatmap.plot(diff_data, diff_result$sample)
    
    # 保存结果
    comparison_name <- sub(group_label, pattern = "/", replacement = "-vs-")
    
    # 保存差异表达蛋白列表
    fwrite(diff_data, file=file.path(output_dir, paste0("Table 2. DEPs list of ", comparison_name, ".xls")), sep = "\t")
    
    # 保存火山图
    ggsave(volcano_plot, filename=file.path(output_dir, paste0("Figure 2a. DEPs volcano of ", comparison_name, ".pdf")), width=8, height=8)
    ggsave(volcano_plot, filename=file.path(output_dir, paste0("Figure 2a. DEPs volcano of ", comparison_name, ".png")), width=8, height=8)
    
    # 保存热图
    if(!is.null(heatmap_plot)) {
      ggsave(heatmap_plot$gtable, filename=file.path(output_dir, paste0("Figure 2b. DEPs heatmap of ", comparison_name, ".pdf")), width=6, height=8, bg='white')
      ggsave(heatmap_plot$gtable, filename=file.path(output_dir, paste0("Figure 2b. DEPs heatmap of ", comparison_name, ".png")), width=6, height=8, bg='white')
    }
    
    # 更新结果汇总
    difff_data_sub=subset(diff_data,Regulation!='unchanged')
    t1=table(difff_data_sub %>% select(Regulation)) %>% as.data.frame() %>% mutate(Type="Sites",Compare=group_label)
    t2=table(distinct(difff_data_sub,`Protein accession`,Regulation) %>% select(Regulation)) %>% as.data.frame() %>% mutate(Type="Protein",Compare=group_label)   
    diff_summary=rbind(diff_summary,t1,t2)
  }

  fwrite(diff_summary, file=file.path(output_dir, 'Statistics of DEPs.xls'), sep = '\t')

  # 返回结果汇总
  return(diff_summary)
}

data_quant_file <- snakemake@input[["data_quant"]]
sample_info_file <- snakemake@input[["sample_file"]]
compare_group_file <- snakemake@input[["compare_group"]]
diff_method <- snakemake@params[["diff_method"]]
diff_paired <- snakemake@params[["diff_paired"]]
diff_fc <- snakemake@params[['diff_fc']]
diff_adj <- snakemake@params[['diff_adj']]
diff_fdr <- snakemake@params[['diff_fdr']]
output_dir <- snakemake@params[["output_dir"]]

data.quant = fread(data_quant_file)
sample.list = fread(sample_info_file)
compare.group = fread(compare_group_file)

# 假设已经有了定量数据和样本分组信息
diff_summary <- generate_diff_report(
  data.quant = data.quant,
  sample.list = sample.list,
  compare.group = compare.group,
  output_dir = output_dir,
  diff_method = diff_method,
  paired = diff_paired,
  adj = diff_adj,
  fdr = diff_fdr,
  fc = diff_fc
)


