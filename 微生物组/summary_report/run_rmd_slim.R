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
krona_example_pic <- get_absolute_path(args[23])
alpha_stats <- get_absolute_path(args[24])
alpha_feature_curve_dir <- get_absolute_path(args[25])
alpha_shannon_curve_dir <- get_absolute_path(args[26])

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
    krona_example_pic = krona_example_pic,
    alpha_stats = alpha_stats,
    alpha_feature_curve_dir = alpha_feature_curve_dir,
    alpha_shannon_curve_dir = alpha_shannon_curve_dir
  )
)