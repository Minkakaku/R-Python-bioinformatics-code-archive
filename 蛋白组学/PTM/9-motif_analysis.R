library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(seqinr)
library(stringr)
library(rmotifx)
library(ggpubr)
library(ggseqlogo)

##################### parameter #########################################
data_quant_file <- snakemake@input[["data_quant"]]
db.dir=snakemake@params[["db_dir"]]
modified.type=snakemake@params[["modified_type"]]
out.dir=snakemake@params[["out_dir"]]


##################### function defined #########################################
motif.seq.file <- function(file,sites,width){
  t=distinct(file,`Protein accession`,Sequence)
  seqs=setNames(t$Sequence,t$`Protein accession`)
  result_list <- lapply(names(seqs), function(nm) {
    extract_flanking_peptides(nm, seqs[[nm]], sites, width)
  })
  bg.seq.final <- do.call(rbind, result_list)
  ptm.seq <- merge(subset(file,select = c("Protein accession","Position","Amino acid")),
                   bg.seq.final, by = c("Protein accession","Position","Amino acid"))
  return(list(momo.bg = bg.seq.final, momo.fg = ptm.seq))
}
extract_flanking_peptides <- function(seqname, seq, target_aa, window) {
  seq <- as.character(seq)
  target_positions <- unique(unlist(str_locate_all(seq,pattern = target_aa)))
  peptides <- character(0)
  aa_list <- character(0)
  
  for (pos in target_positions[target_positions > 0]) {
    start <- pos - window
    end <- pos + window
    
    # 提取并补齐
    frag <- substring(seq,
                      first = max(start, 1),
                      last  = min(end, nchar(seq)))
    
    prefix_pad <- if (start < 1) paste0(rep("-", 1 - start), collapse = "") else ""
    suffix_pad <- if (end > nchar(seq)) paste0(rep("-", end - nchar(seq)), collapse = "") else ""
    
    padded_frag <- paste0(prefix_pad, frag, suffix_pad)
    
    peptides <- c(peptides, padded_frag)
    aa_list  <- c(aa_list, substr(seq, pos, pos))
  }
  
  data.frame(
    `Protein accession` = seqname,
    Position = target_positions[target_positions > 0],
    `Amino acid` = aa_list,
    frag = peptides,
    stringsAsFactors = FALSE,check.names = F
  )
}
site.model <- function(motif.data,amino,site,len){
  bg <- unique(subset(motif.data$momo.bg,`Amino acid` == site)$frag)
  bg <- as.data.frame(str_split_fixed(bg,pattern = "",n = nchar(bg[1])))
  fg <- unique(subset(motif.data$momo.fg,`Amino acid` == site)$frag)
  fg <- as.data.frame(str_split_fixed(fg,pattern = "",n = nchar(fg[1])))
  colnames(fg)=seq(-1*len,len,by =1)
  colnames(bg)=seq(-1*len,len,by =1)
  enrich.data <- lapply(colnames(fg),function(x){enrich.motif.demo(fg[x],bg[x],amino,site)})
  names(enrich.data) <- seq(-1*len,len,by =1)
  enrich.data.final <- rbindlist(enrich.data,idcol = T)
  colnames(enrich.data.final) <- c("Position","Amino acid","logp")
  enrich.data.final[enrich.data.final$logp == Inf,"logp"] <- max(enrich.data.final$logp[is.finite(enrich.data.final$logp)])
  enrich.data.final[enrich.data.final$logp == -Inf,"logp"] <- min(enrich.data.final$logp[is.finite(enrich.data.final$logp)])
  enrich.data.final$Position = factor(enrich.data.final$Position,levels = seq(-1*len,len))
  return(enrich.data.final)
}
enrich.motif.demo <- function(fg.data,bg.data,amino,site){
  if(all(unique(fg.data[[1]]) %in% site)){
    data.state = data.frame(`Amino acid` = amino,`-logp` = rep(0,length(amino)),check.names = F)
  }else{
    mapping = as.data.frame(table(fg.data[[1]]))
    all_mapping = as.data.frame(table(bg.data[[1]]))
    data = merge(mapping,all_mapping,by = "Var1") %>% mutate(Var1 = as.character(Var1))
    colnames(data) = c("Amino acid","mapping","all mapping")
    data = data %>% mutate(background = nrow(fg.data)-mapping,`all background` = nrow(bg.data)-`all mapping`) %>%
      filter(`Amino acid` %in% amino)
    data$p.value = signif(apply(data[2:5],1,chis.test),3)
    data$`Fold Enrich` = round(data$mapping/data$background/data$`all mapping`*data$`all background`,2)
    data = data %>% mutate(`-logp` = case_when(`Fold Enrich` > 1 ~ -log10(p.value),TRUE ~ log10(p.value)))
    data.state <- subset(data,select = c("Amino acid","-logp"))
    
    if(nrow(data.state) < length(amino)){
      ami <- setdiff(amino,data.state$`Amino acid`)
      ma <- data.frame(`Amino acid` = ami,`-logp` = rep(0,length(ami)),check.names = F)
      data.state <- rbind(data.state,ma)
      data.state <- data.state %>% arrange(`Amino acid`)
    }else{
      data.state <- data.state %>% arrange(`Amino acid`)
    }        
  }
  colnames(data.state) <- c("Amino acid",colnames(fg.data))
  return(data.state)
}
chis.test <- function(x){
  ma <- matrix(as.numeric(x),nrow = 2)
  p.value <- chisq.test(ma)$p.value
  return(p.value)
}
site.heatmap <- function(data,color1,color2,color3,n){
  data = as.data.frame(data)
  data$`Amino acid` = factor(data$`Amino acid`,levels = rev(unique(data$`Amino acid`)))
  p <- ggplot(data,aes(Position,`Amino acid`,fill = logp))+
    geom_tile(color = "grey")+
    scale_fill_gradient2(limits = c((-1)*n,n),midpoint = 0,low = color1,mid = color2,high = color3,oob = scales::squish)+
    guides(fill = guide_colorbar(title = NULL))+
    theme_minimal(base_size = 14)
  return(p)
}
##################### Execute the main program #########################################

## data load
sequence.out <- read.fasta(paste(db.dir,"sequence.fasta",sep="/"),seqtype = "AA",as.string = T)
names(sequence.out)=unlist(lapply(str_split(names(sequence.out),pattern = "\\|"),function(x) x[2]))
sequence <- sapply(sequence.out,"[",1)
sequence <- as.data.frame(sequence) %>% tibble::rownames_to_column("Protein accession")
colnames(sequence) <- c("Protein accession","Sequence")
ptm.seq.temp <- merge(subset(fread(data_quant_file),select = c("Protein accession","Amino acid","Position")), sequence, by = "Protein accession")
sites=unique(ptm.seq.temp$`Amino acid`)
if(modified.type=='Phosphorylation (STY)') width=6 else width=10

## motif peptide 
motif.file = motif.seq.file(ptm.seq.temp,sites,width)

## output

# if(!dir.exists('9-Motif_analysis')) dir.create('9-Motif_analysis')
if(!dir.exists(file.path(out.dir, 'motif_image'))) dir.create(file.path(out.dir, 'motif_image'))

## motif model heatmap
amino = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
lapply(sites,function(x){
  out.data=site.model(motif.file,amino,x,width)
  out.plot=site.heatmap(out.data,color1 = "#39B59D",color2 = "#FDFDFD",color3 = "#E24C4F",n = 20)
  fwrite(out.data,paste0(out.dir, "/Table 7. ",x,".motif.model.xls"),sep="\t")
  ggsave(plot = out.plot,filename = paste0(out.dir,"/Figure 7. Heatmap of ",x,".motif.model.pdf"),width = 6,height = 5)
  ggsave(plot = out.plot,filename = paste0(out.dir,"/Figure 7. Heatmap of ",x,".motif.model.png"),width = 6,height = 5,bg = 'white')
})


result=rbindlist(lapply(sites,function(x){
  fg=motif.file$momo.fg$frag[motif.file$momo.fg$`Amino acid`==x]
  bg=motif.file$momo.bg$frag[motif.file$momo.bg$`Amino acid`==x]
  n=log2(ceiling(length(fg)/1000)+0.3)*50
  tmp = motifx(fg,bg,central.res = x, min.seqs = n, pval.cutoff = 1e-6)
  if(length(tmp) != 0){
    mot = tmp %>%
      mutate(score=round(score,2),fold.increase=round(fold.increase,2)) %>%
      mutate(Logo=paste0('<img src="9-Motif_analysis/motif_image/',motif,'.png" height="40"/>'),.before = "motif")
    lapply(mot$motif,function(y){
      p=ggseqlogo(grep(fg,pattern=y,value = T),seq_type="aa",method="probability")+theme_void()+theme(legend.position = 'none')
      ggsave(plot = p,filename = paste0(out.dir,"/motif_image/",y,".png"),width = 5,height = 1)
    })
    return(mot)
  }
}))

fwrite(result,file.path(out.dir,"Table 7. Motif-x.xls"),sep="\t",quote = F)
