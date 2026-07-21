library(optparse)

option_list <- list(
  make_option(c("-t", "--template"), type = "character", default = NULL,
              help = "Path to report.Rmd template", metavar = "FILE"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Directory for rendered html output", metavar = "DIR"),
  make_option(c("-f", "--output_file"), type = "character", default = "report.html",
              help = "Rendered html filename (default: report.html)", metavar = "FILE"),
  make_option(c("-r", "--results_dir"), type = "character", default = NULL,
              help = "Plot/results directory passed into Rmd params", metavar = "DIR"),
  make_option(c("-p", "--rctd_log_png"), type = "character", default = NULL,
              help = "Path to RCTD_log.png passed into Rmd params", metavar = "FILE"),
  make_option(c("-n", "--report_title"), type = "character", default = "Xenium RCTD 反卷积分析报告",
              help = "Report title passed into Rmd params", metavar = "TEXT")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("template", "output_dir", "results_dir", "rctd_log_png")
missing_args <- sapply(required_args, function(x) is.null(opt[[x]]))
if (any(missing_args)) {
  missing_names <- paste(required_args[missing_args], collapse = ", ")
  stop(paste("Missing required arguments:", missing_names), call. = FALSE)
}

if (!file.exists(opt$template)) {
  stop(paste("Report template not found:", opt$template), call. = FALSE)
}

if (!dir.exists(opt$results_dir)) {
  stop(paste("Results directory not found:", opt$results_dir), call. = FALSE)
}

if (!file.exists(opt$rctd_log_png)) {
  stop(paste("RCTD_log.png not found:", opt$rctd_log_png), call. = FALSE)
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render the report.", call. = FALSE)
}

rmarkdown::render(
  input = opt$template,
  output_file = opt$output_file,
  output_dir = opt$output_dir,
  params = list(
    results_dir = opt$results_dir,
    rctd_log_png = opt$rctd_log_png,
    report_title = opt$report_title
  ),
  envir = new.env(parent = globalenv())
)

message("Report rendered successfully: ", file.path(opt$output_dir, opt$output_file))
