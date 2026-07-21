# function to get absolute path otherwise Rmarkdown will raise an error
get_absolute_path <- function(path) {
  absolute_path <- paste0(normalizePath(dirname(path)), '/',  basename(path))
  return(absolute_path)
}

# pass parameters
args <- commandArgs(trailingOnly=TRUE)
rmd_script <- args[1]
output <- args[2]
project_name <- args[3]
project_number <- args[4]
customer_name <- args[5]
customer_unit <- args[6]
amplicon_info <- args[7]
primer_seq_forward <- args[8]
primer_seq_reverse <- args[9]
tax_database <- args[10]
logo_pic <- get_absolute_path(args[11])
metadata <- get_absolute_path(args[12])
contrasts_all <- args[13]
experiment_pipeline_pic <- get_absolute_path(args[14])
analysis_pipeline_pic <- get_absolute_path(args[15])
fastq_qc_plot_dir <- get_absolute_path(args[16])
reads_qc_stats <- get_absolute_path(args[17])
reads_length_barplot_dir <- get_absolute_path(args[18])
all_ASVs_count_stats <- get_absolute_path(args[19])
ASVs_count_barplot_dir <- get_absolute_path(args[20])
ASVs_count_venn_dir <- get_absolute_path(args[21])
tax_barplot_dir <- get_absolute_path(args[22])
tax_heatmap_dir <- get_absolute_path(args[23])
krona_example_pic <- get_absolute_path(args[24])
alpha_stats <- get_absolute_path(args[25])
alpha_barplot_dir <- get_absolute_path(args[26])
alpha_feature_curve_dir <- get_absolute_path(args[27])
alpha_shannon_curve_dir <- get_absolute_path(args[28])
beta_pca_dir <- get_absolute_path(args[29])
beta_pcoa_dir <- get_absolute_path(args[30])
beta_nmds_dir <- get_absolute_path(args[31])
beta_upgma_dir <- get_absolute_path(args[32])
beta_heatmap_dir <- get_absolute_path(args[33])
beta_anosim_dir <- get_absolute_path(args[34])
anova_dir <- get_absolute_path(args[35])
rank_sum_dir <- get_absolute_path(args[36])
lefse_dir <- get_absolute_path(args[37])
lefse_lda_cutoff <- args[38]
lefse_p_cutoff <- args[39]
network_dir <- get_absolute_path(args[40])
spearman_r_cutoff <- args[41]
spearman_p_cutoff <- args[42]
picrust2_flowchart <- get_absolute_path(args[43])
func_kegg_dir <- get_absolute_path(args[44])
func_ec_dir <- get_absolute_path(args[45])
func_metacyc_dir <- get_absolute_path(args[46])

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
    amplicon_info = amplicon_info,
    primer_seq_forward = primer_seq_forward,
    primer_seq_reverse = primer_seq_reverse,
    tax_database = tax_database,
    logo_pic = logo_pic,
    metadata = metadata,
    contrasts_all = contrasts_all,
    experiment_pipeline_pic = experiment_pipeline_pic,
    analysis_pipeline_pic = analysis_pipeline_pic,
    fastq_qc_plot_dir = fastq_qc_plot_dir,
    reads_qc_stats = reads_qc_stats,
    reads_length_barplot_dir = reads_length_barplot_dir,
    all_ASVs_count_stats =  all_ASVs_count_stats,
    ASVs_count_barplot_dir = ASVs_count_barplot_dir,
    ASVs_count_venn_dir = ASVs_count_venn_dir,
    tax_barplot_dir = tax_barplot_dir,
    tax_heatmap_dir = tax_heatmap_dir,
    krona_example_pic = krona_example_pic,
    alpha_stats = alpha_stats,
    alpha_barplot_dir = alpha_barplot_dir,
    alpha_feature_curve_dir = alpha_feature_curve_dir,
    alpha_shannon_curve_dir = alpha_shannon_curve_dir,
    beta_pca_dir = beta_pca_dir,
    beta_pcoa_dir = beta_pcoa_dir,
    beta_nmds_dir = beta_nmds_dir,
    beta_upgma_dir = beta_upgma_dir,
    beta_heatmap_dir = beta_heatmap_dir,
    beta_anosim_dir = beta_anosim_dir,
    anova_dir = anova_dir,
    rank_sum_dir = rank_sum_dir,
    lefse_dir = lefse_dir,
    lefse_lda_cutoff = lefse_lda_cutoff,
    lefse_p_cutoff = lefse_p_cutoff,
    network_dir = network_dir,
    spearman_r_cutoff = spearman_r_cutoff,
    spearman_p_cutoff = spearman_p_cutoff,
    picrust2_flowchart = picrust2_flowchart,
    func_kegg_dir = func_kegg_dir,
    func_ec_dir = func_ec_dir,
    func_metacyc_dir = func_metacyc_dir
  )
)