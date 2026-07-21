# pass parameters
args <- commandArgs(trailingOnly=TRUE)
rmd_script <- args[1]
output <- args[2]
project_name <- args[3]
project_number <- args[4]
customer_name <- args[5]
customer_unit <- args[6]
logo_pic <- args[7]
experiment_pipeline_pic <- args[8]
analysis_pipeline_pic <- args[9]
meta_data_file <- args[10]
fastq_qc_plot_dir <- args[11]
fastp_qc_stats <- args[12]
ref_genome <- args[13]
ref_genome_download <- args[14]
map_stats <- args[15]
hisat2_dir <- args[16]
library_assessment_dir <- args[17]
featureCounts_merged <- args[18]
fpkm_merged <- args[19]
fpkm_boxplot <- args[20]
fpkm_heatmap <- args[21]
fpkm_pca <- args[22]
diffexp_FC <- args[23]
diffexp_p_criteria <- args[24]
diffexp_dir <- args[25]
diffexp_stats <- args[26]
diffexp_stats_barplot <- args[27]
enrich_principle_pic <- args[28]
enrich_dir <- args[29]
pathview_example_pic <- args[30]
rmats_events_pic <- args[31]
rmats_fdr_cutoff <- args[32]
rmats_diff_cutoff <- args[33]
rmats_diff_stats <- args[34]
rmats2sashimiplot_example_plot <- args[35]
var_dir <- args[36]
gene_fusion_description_pic <- args[37]
fusion_results_dir <- args[38]
tf_dir <- args[39]
cuffcmp_gtf <- args[40]
cuffcmp_novel_gene <- args[41]
dexseq_padj_cutoff <- args[42]
dexseq_dir <- args[43]
dexseq_example_pic <- args[44]
recomend_params_table <- args[45]

# absolute path for these files otherwise Rmarkdown will raise an error
logo_pic <- paste0(normalizePath(dirname(logo_pic)), '/',  basename(logo_pic))
experiment_pipeline_pic <- paste0(normalizePath(dirname(experiment_pipeline_pic)), '/',  basename(experiment_pipeline_pic))
analysis_pipeline_pic <- paste0(normalizePath(dirname(analysis_pipeline_pic)), '/',  basename(analysis_pipeline_pic))
meta_data_file <- paste0(normalizePath(dirname(meta_data_file)), '/',  basename(meta_data_file))
fastq_qc_plot_dir <- paste0(normalizePath(dirname(fastq_qc_plot_dir)), '/',  basename(fastq_qc_plot_dir))
fastp_qc_stats <- paste0(normalizePath(dirname(fastp_qc_stats)), '/',  basename(fastp_qc_stats))
map_stats <- paste0(normalizePath(dirname(map_stats)), '/',  basename(map_stats))
hisat2_dir <- paste0(normalizePath(dirname(hisat2_dir)), '/',  basename(hisat2_dir))
library_assessment_dir <- paste0(normalizePath(dirname(library_assessment_dir)), '/',  basename(library_assessment_dir))
featureCounts_merged <- paste0(normalizePath(dirname(featureCounts_merged)), '/',  basename(featureCounts_merged))
fpkm_merged <- paste0(normalizePath(dirname(fpkm_merged)), '/',  basename(fpkm_merged))
fpkm_boxplot <- paste0(normalizePath(dirname(fpkm_boxplot)), '/',  basename(fpkm_boxplot))
fpkm_heatmap <- paste0(normalizePath(dirname(fpkm_heatmap)), '/',  basename(fpkm_heatmap))
fpkm_pca <- paste0(normalizePath(dirname(fpkm_pca)), '/',  basename(fpkm_pca))
diffexp_dir <- paste0(normalizePath(dirname(diffexp_dir)), '/',  basename(diffexp_dir))
diffexp_stats <- paste0(normalizePath(dirname(diffexp_stats)), '/',  basename(diffexp_stats))
diffexp_stats_barplot <- paste0(normalizePath(dirname(diffexp_stats_barplot)), '/',  basename(diffexp_stats_barplot))
enrich_principle_pic <- paste0(normalizePath(dirname(enrich_principle_pic)), '/',  basename(enrich_principle_pic))
enrich_dir <- paste0(normalizePath(dirname(enrich_dir)), '/',  basename(enrich_dir))
pathview_example_pic <- paste0(normalizePath(dirname(pathview_example_pic)), '/',  basename(pathview_example_pic))
rmats_events_pic <- paste0(normalizePath(dirname(rmats_events_pic)), '/',  basename(rmats_events_pic))
rmats_diff_stats <- paste0(normalizePath(dirname(rmats_diff_stats)), '/',  basename(rmats_diff_stats))
rmats2sashimiplot_example_plot <- paste0(normalizePath(dirname(rmats2sashimiplot_example_plot)), '/',  basename(rmats2sashimiplot_example_plot))
var_dir <- paste0(normalizePath(dirname(var_dir)), '/',  basename(var_dir))
gene_fusion_description_pic <- paste0(normalizePath(dirname(gene_fusion_description_pic)), '/',  basename(gene_fusion_description_pic))
fusion_results_dir <- paste0(normalizePath(dirname(fusion_results_dir)), '/',  basename(fusion_results_dir))
tf_dir <- paste0(normalizePath(dirname(tf_dir)), '/',  basename(tf_dir))
cuffcmp_gtf <- paste0(normalizePath(dirname(cuffcmp_gtf)), '/',  basename(cuffcmp_gtf))
cuffcmp_novel_gene <- paste0(normalizePath(dirname(cuffcmp_novel_gene)), '/',  basename(cuffcmp_novel_gene))
dexseq_dir <- paste0(normalizePath(dirname(dexseq_dir)), '/',  basename(dexseq_dir))
dexseq_example_pic <- paste0(normalizePath(dirname(dexseq_example_pic)), '/',  basename(dexseq_example_pic))
recomend_params_table <- paste0(normalizePath(dirname(recomend_params_table)), '/',  basename(recomend_params_table))

# run Rmarkdown script
rmarkdown::render(
  rmd_script,
  output_file = basename(output),
  output_dir = dirname(output),
  params = list(
    project_name = project_name,
    project_number = project_number,
    customer_name = customer_name,
    customer_unit = customer_unit,
    logo_pic = logo_pic,
    experiment_pipeline_pic = experiment_pipeline_pic,
    analysis_pipeline_pic = analysis_pipeline_pic,
    meta_data_file = meta_data_file,
    fastq_qc_plot_dir = fastq_qc_plot_dir,
    fastp_qc_stats = fastp_qc_stats,
    ref_genome = ref_genome,
    ref_genome_download = ref_genome_download,
    map_stats = map_stats,
    hisat2_dir = hisat2_dir,
    library_assessment_dir = library_assessment_dir,
    featureCounts_merged = featureCounts_merged,
    fpkm_merged = fpkm_merged,
    fpkm_boxplot = fpkm_boxplot,
    fpkm_heatmap = fpkm_heatmap,
    fpkm_pca = fpkm_pca,
    diffexp_FC = diffexp_FC,
    diffexp_p_criteria = diffexp_p_criteria,
    diffexp_dir = diffexp_dir,
    diffexp_stats = diffexp_stats,
    diffexp_stats_barplot = diffexp_stats_barplot,
    enrich_principle_pic = enrich_principle_pic,
    enrich_dir = enrich_dir,
    pathview_example_pic = pathview_example_pic,
    rmats_events_pic = rmats_events_pic,
    rmats_fdr_cutoff = rmats_fdr_cutoff,
    rmats_diff_cutoff = rmats_diff_cutoff,
    rmats_diff_stats = rmats_diff_stats,
    rmats2sashimiplot_example_plot = rmats2sashimiplot_example_plot,
    var_dir =var_dir,
    gene_fusion_description_pic = gene_fusion_description_pic,
    fusion_results_dir = fusion_results_dir,
    tf_dir = tf_dir,
    cuffcmp_gtf = cuffcmp_gtf,
    cuffcmp_novel_gene = cuffcmp_novel_gene,
    dexseq_padj_cutoff = dexseq_padj_cutoff,
    dexseq_dir = dexseq_dir,
    dexseq_example_pic = dexseq_example_pic,
    recomend_params_table = recomend_params_table
  )
)
