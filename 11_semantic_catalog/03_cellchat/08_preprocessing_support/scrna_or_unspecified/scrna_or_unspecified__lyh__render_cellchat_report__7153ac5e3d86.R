#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--input_dir"), type = "character", default = "."),
  make_option(c("--report_rmd"), type = "character", default = "cellchat_report.Rmd"),
  make_option(c("--output_file"), type = "character", default = "CellChat_Report.html")
)

opt <- parse_args(OptionParser(option_list = option_list))

input_dir <- normalizePath(opt$input_dir, winslash = "/", mustWork = TRUE)
report_rmd <- normalizePath(opt$report_rmd, winslash = "/", mustWork = TRUE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required to render the report. Please install it before running the final report step.")
}

rmarkdown::render(
  input = report_rmd,
  output_file = opt$output_file,
  output_dir = input_dir,
  params = list(input_dir = input_dir),
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

message(sprintf("Report generated: %s", file.path(input_dir, opt$output_file)))
