log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(chimeraviz)

# Input
## STAR Fusion results 
star_fusion_results <- snakemake@input[[1]]
## Fusion circos plot output
circos_png <- snakemake@output[[1]]
circos_pdf <- snakemake@output[[2]]
## Chain file for R liftover package
mm39_mm10_chain <- snakemake@params[["mm39_mm10_chain"]]
## Reference genome version
if (snakemake@params[["species"]] == "homo_sapiens") {
  genome_version <- "hg38"
} else if (snakemake@params[["species"]] == "mus_musculus") {
  genome_version <- "mm10"
}  

if (nrow(read.csv(star_fusion_results, header=1, sep="\t")) != 0) {
  # Load the results file into a list of fusion objects
  fusions <- import_starfusion(star_fusion_results, genome_version)
  
  # Use R package liftover to convert coordinates from mm39 to mm10, 
  #  since chimeraviz only supports mouse genome version mm10.
  if (genome_version == "mm10") {
    library(rtracklayer)
    ch <- rtracklayer::import.chain(mm39_mm10_chain)
 
    convert_error_list <- c()
    for (i in 1:length(fusions)) {
      # gene_upstream
      chr <- fusions[[i]]@gene_upstream@chromosome
      bp <- fusions[[i]]@gene_upstream@breakpoint
      grObject <- GRanges(seqnames=c(chr), ranges=IRanges(start=bp, end=bp))
      grObject_mm10 <- rtracklayer::liftOver(grObject, ch)
      bp_mm10 <- grObject_mm10[[1]]@ranges@start
      # modify fusion breakpoint from mm39 to mm10
      fusions[[i]]@gene_upstream@breakpoint <- bp_mm10
      
      # if fusion has coordinates which cannot be converted to mm10
      if (length(bp_mm10) == 0) {
        convert_error_list <- c(convert_error_list, i)
      }
      
      # gene_downstream
      chr <- fusions[[i]]@gene_downstream@chromosome
      bp <- fusions[[i]]@gene_downstream@breakpoint
      grObject <- GRanges(seqnames=c(chr), ranges=IRanges(start=bp, end=bp))
      grObject_mm10 <- rtracklayer::liftOver(grObject, ch)
      bp_mm10 <- grObject_mm10[[1]]@ranges@start
      # modify fusion breakpoint from mm39 to mm10
      fusions[[i]]@gene_downstream@breakpoint <- bp_mm10
      
      # if fusion has coordinates which cannot be converted to mm10
      if (length(bp_mm10) == 0) {
        convert_error_list <- c(convert_error_list, i)
      }
    }
    
    # Remove gene fusions whose coordinates fail to be converted to mm10
    if (length(convert_error_list) != 0) {
      fusions_filtered <- fusions[-unique(convert_error_list)]
    } else {
      fusions_filtered <- fusions
    }
    
    
    # Remove gene fusions related to chrM since chimeraviz database doesn't have cytoband information for mm10 chrM.
    fusions_chrM_id <- c()
    for (id in c(1:length(fusions_filtered))) {
      if (length(get_fusion_by_chromosome(fusions_filtered[id], "chrM")) != 0) {
        fusions_chrM_id <- c(fusions_chrM_id, id)
      }
    }
    if (length(fusions_chrM_id) != 0) {
      fusions_filtered <- fusions_filtered[-fusions_chrM_id]
    } else {
      fusions_filtered <- fusions_filtered
    }
  } else {
    fusions_filtered <- fusions
  }
  
  # After removing chrM-related fusions, fusions_filtered might be empty
  if (length(fusions_filtered) != 0) {
    # Output fusion circos plot
    png(filename = circos_png,
        units = "in",
        width = 10,
        height = 10,
        res = 400)
    plot_circle(fusions_filtered)
    dev.off()
    pdf(file = circos_pdf,
        width = 10,
        height = 10)
    plot_circle(fusions_filtered)
    dev.off()
  } else {
    # Output an empty plot
    # Output an empty plot
    png(filename = circos_png,
        units = "in",
        width = 10,
        height = 10,
        res = 400)
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("未找到可以用于作图的融合基因！"), 
         cex = 1.6, col = "black")
    dev.off()
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    cairo_pdf(file = circos_pdf,
              width = 10,
              height = 10)
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("未找到可以用于作图的融合基因！"), 
         cex = 1.6, col = "black")
    dev.off()
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
} else {
  # Output an empty plot
  png(filename = circos_png,
      units = "in",
      width = 10,
      height = 10,
      res = 400)
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("未找到可以用于作图的融合基因！"), 
       cex = 1.6, col = "black")
  dev.off()
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  cairo_pdf(file = circos_pdf,
            width = 10,
            height = 10)
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("未找到可以用于作图的融合基因！"), 
       cex = 1.6, col = "black")
  dev.off()
  par(mar = c(5, 4, 4, 2) + 0.1)
}
