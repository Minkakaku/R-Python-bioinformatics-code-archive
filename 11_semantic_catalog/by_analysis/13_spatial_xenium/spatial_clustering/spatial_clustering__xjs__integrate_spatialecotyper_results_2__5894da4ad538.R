#!/usr/bin/env Rscript

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript spatial_ecotypes/integrate_spatialecotyper_results.R [options]",
      "",
      "Purpose:",
      "  Integrate finished single-sample SpatialEcoTyper results into conserved multi-sample SEs.",
      "  By default, sample directories are selected from the sample column in all_data_4slides_sample_filtered_reUMAP.qs2.",
      "",
      "Common options:",
      "  -h, --help                         Show this help message and exit.",
      "  --dry-run                          Only discover samples and write manifest; do not run integration.",
      "  --root=DIR                         Directory containing single-sample result folders.",
      "                                     Default: directory containing this script.",
      "  --outdir=DIR                       Integration output directory.",
      "                                     Default: <root>/integrated_SE.",
      "  --input-qs=FILE                    Seurat QS file used to read the sample whitelist.",
      "                                     Default: ../all_data_4slides_sample_filtered_reUMAP.qs2 relative to root.",
      "  --sample-col=COL                   Metadata column in --input-qs containing sample IDs.",
      "                                     Default: sample.",
      "  --samples=S1,S2,S3                 Optional comma-separated subset of samples to integrate.",
      "  --diagnose-subresolution=2:12      Test candidate subresolution values and exit before integration.",
      "                                     Writes subresolution_diagnosis_*.csv under --outdir.",
      "",
      "IntegrateSpatialEcoTyper options:",
      "  --normalization=METHOD             Default: None.",
      "  --nmf-ranks=10, --nmf-ranks=8,9,10, or --nmf-ranks=4:15",
      "                                     Underscore aliases like --nmf_ranks also work.",
      "                                     NMF rank(s) to test. Commas and inclusive ranges can be mixed.",
      "                                     Default: 10.",
      "  --nrun-per-rank=N                  NMF runs per rank. Default: 30.",
      "  --min-coph=FLOAT                   Minimum cophenetic coefficient. Default: 0.95.",
      "  --nfeatures=N                      Max top variable genes per cell type. Default: 300.",
      "  --min-features=N                   Minimum shared features across samples. Default: 10.",
      "  --region=COL                       Optional metadata region column. Default: NULL.",
      "  --downsample-by-region=true|false  Balance spatial clusters across regions. Default: true.",
      "  --subresolution=FLOAT              Within-sample clustering resolution. Default: 30.",
      "  --minibatch=N                      SNF minibatch size. Default: 5000.",
      "  --ncores=N                         Parallel cores. Default: 15.",
      "  --seed=N                           Random seed. Default: 1.",
      "",
      "Examples:",
      "  Rscript spatial_ecotypes/integrate_spatialecotyper_results.R --dry-run",
      "  Rscript spatial_ecotypes/integrate_spatialecotyper_results.R --ncores=15 --nmf-ranks=4:15",
      "  Rscript spatial_ecotypes/integrate_spatialecotyper_results.R --samples=SS01,SS02,WJ01 --dry-run",
      "",
      sep = "\n"
    )
  )
}

raw_args <- commandArgs(trailingOnly = TRUE)
if (any(raw_args %in% c("-h", "--help", "help"))) {
  print_usage()
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SpatialEcoTyper))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

parse_cli <- function(args) {
  opts <- list()
  flags <- character(0)

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }

    arg <- sub("^--", "", arg)
    if (grepl("=", arg, fixed = TRUE)) {
      key <- gsub("_", "-", sub("=.*$", "", arg))
      value <- sub("^[^=]+=", "", arg)
      opts[[key]] <- value
    } else {
      flags <- c(flags, gsub("_", "-", arg))
    }
  }

  opts$flags <- flags
  opts
}

as_bool <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

as_int_vector <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }

  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  values <- unlist(lapply(parts, function(part) {
    part <- trimws(part)
    if (!nzchar(part)) {
      return(integer(0))
    }

    if (grepl("^[0-9]+:[0-9]+$", part)) {
      bounds <- as.integer(strsplit(part, ":", fixed = TRUE)[[1]])
      return(seq.int(bounds[[1]], bounds[[2]]))
    }

    if (!grepl("^[0-9]+$", part)) {
      stop("Invalid integer/range value: ", part)
    }
    as.integer(part)
  }))

  unique(as.integer(values))
}

as_numeric_vector <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }

  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  values <- unlist(lapply(parts, function(part) {
    part <- trimws(part)
    if (!nzchar(part)) {
      return(numeric(0))
    }

    if (grepl("^[0-9.]+:[0-9.]+$", part)) {
      bounds <- as.numeric(strsplit(part, ":", fixed = TRUE)[[1]])
      if (any(is.na(bounds)) || any(bounds != as.integer(bounds))) {
        stop("Range syntax only supports integer endpoints: ", part)
      }
      return(seq.int(as.integer(bounds[[1]]), as.integer(bounds[[2]])))
    }

    value <- suppressWarnings(as.numeric(part))
    if (is.na(value)) {
      stop("Invalid numeric/range value: ", part)
    }
    value
  }))

  unique(as.numeric(values))
}

message_time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", ...)
}

make_sample_safe <- function(sample_id) {
  gsub("[^A-Za-z0-9_.-]+", "_", sample_id)
}

read_sample_map_from_qs <- function(input_qs, sample_col, sample_include = NULL) {
  adata <- qs2::qs_read(input_qs)
  if (!sample_col %in% colnames(adata@meta.data)) {
    stop("Missing sample column in input QS: ", sample_col)
  }

  sample_ids <- sort(unique(as.character(adata@meta.data[[sample_col]])))
  sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]
  if (!is.null(sample_include)) {
    sample_ids <- intersect(sample_include, sample_ids)
  }
  if (length(sample_ids) == 0) {
    stop("No samples found from ", input_qs, " column ", sample_col)
  }

  data.frame(
    sample = sample_ids,
    sample_safe = make_sample_safe(sample_ids),
    stringsAsFactors = FALSE
  )
}

discover_single_sample_results <- function(
    root_dir,
    out_dir = NULL,
    sample_include = NULL,
    sample_map = NULL,
    skip_dir_names = c("integrated_SE")) {
  if (!is.null(sample_map)) {
    sample_dirs <- file.path(root_dir, sample_map$sample_safe)
    sample_labels <- sample_map$sample
    sample_safe_labels <- sample_map$sample_safe
  } else {
    sample_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[file.info(sample_dirs)$isdir]
    sample_dirs <- sample_dirs[!basename(sample_dirs) %in% skip_dir_names]
    sample_labels <- basename(sample_dirs)
    sample_safe_labels <- basename(sample_dirs)
  }

  if (!is.null(out_dir)) {
    keep <- normalizePath(sample_dirs, mustWork = FALSE) != normalizePath(out_dir, mustWork = FALSE)
    sample_dirs <- sample_dirs[keep]
    sample_labels <- sample_labels[keep]
    sample_safe_labels <- sample_safe_labels[keep]
  }

  records <- lapply(seq_along(sample_dirs), function(i) {
    sample_dir <- sample_dirs[[i]]
    sample_id <- sample_labels[[i]]
    sample_safe <- sample_safe_labels[[i]]

    if (!is.null(sample_include) && !sample_id %in% sample_include) {
      return(NULL)
    }
    if (!dir.exists(sample_dir)) {
      return(data.frame(
        sample = sample_id,
        sample_safe = sample_safe,
        sample_dir = sample_dir,
        result_rds = NA_character_,
        scdata_qs = NA_character_,
        n_result_files = 0L,
        n_scdata_files = 0L,
        stringsAsFactors = FALSE
      ))
    }

    expected_result <- file.path(sample_dir, paste0(sample_safe, "_SpatialEcoTyper_results.rds"))
    expected_scdata <- file.path(sample_dir, paste0(sample_safe, "_scdata.qs2"))

    result_files <- if (file.exists(expected_result)) {
      expected_result
    } else {
      list.files(
        sample_dir,
        pattern = "_SpatialEcoTyper_results[.]rds$",
        full.names = TRUE
      )
    }

    scdata_files <- if (file.exists(expected_scdata)) {
      expected_scdata
    } else {
      list.files(sample_dir, pattern = "_scdata[.]qs2$", full.names = TRUE)
    }

    if (length(result_files) == 0 && length(scdata_files) == 0) {
      return(NULL)
    }

    data.frame(
      sample = sample_id,
      sample_safe = sample_safe,
      sample_dir = sample_dir,
      result_rds = if (length(result_files) == 1) result_files[[1]] else NA_character_,
      scdata_qs = if (length(scdata_files) == 1) scdata_files[[1]] else NA_character_,
      n_result_files = length(result_files),
      n_scdata_files = length(scdata_files),
      stringsAsFactors = FALSE
    )
  })

  manifest <- do.call(rbind, records)
  if (is.null(manifest) || nrow(manifest) == 0) {
    stop("No single-sample SpatialEcoTyper results were found under: ", root_dir)
  }

  manifest <- manifest[order(manifest$sample), , drop = FALSE]
  rownames(manifest) <- NULL
  manifest
}

validate_manifest <- function(manifest) {
  bad <- manifest[
    is.na(manifest$result_rds) |
      is.na(manifest$scdata_qs) |
      manifest$n_result_files != 1 |
      manifest$n_scdata_files != 1,
    ,
    drop = FALSE
  ]

  if (nrow(bad) == 0) {
    return(invisible(TRUE))
  }

  message("Invalid or ambiguous sample directories:")
  print(bad)
  stop(
    "Each sample directory must contain exactly one *_SpatialEcoTyper_results.rds ",
    "and exactly one *_scdata.qs2 file."
  )
}

read_spatialecotyper_result <- function(path, sample_id) {
  result <- readRDS(path)
  if (!is.list(result) || !all(c("obj", "metadata") %in% names(result))) {
    stop("Unexpected SpatialEcoTyper result structure for sample ", sample_id, ": ", path)
  }
  result
}

check_findclusters_eligible_samples <- function(spatialecotyper_list, manifest, subresolution) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    sample_id <- manifest$sample[[i]]
    result <- spatialecotyper_list[[sample_id]]

    tryCatch({
      obj <- Seurat::FindClusters(result$obj, resolution = subresolution, verbose = FALSE)
      clusters <- tryCatch(obj$seurat_clusters, error = function(e) NULL)

      data.frame(
        sample = sample_id,
        sample_safe = manifest$sample_safe[[i]],
        sample_dir = manifest$sample_dir[[i]],
        result_rds = manifest$result_rds[[i]],
        status = "ok",
        keep = TRUE,
        subresolution = subresolution,
        n_cells = ncol(obj),
        n_clusters = if (!is.null(clusters)) length(unique(clusters[!is.na(clusters)])) else NA_integer_,
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        sample = sample_id,
        sample_safe = manifest$sample_safe[[i]],
        sample_dir = manifest$sample_dir[[i]],
        result_rds = manifest$result_rds[[i]],
        status = "failed",
        keep = FALSE,
        subresolution = subresolution,
        n_cells = tryCatch(ncol(result$obj), error = function(e2) NA_integer_),
        n_clusters = NA_integer_,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
  })

  report <- do.call(rbind, rows)
  rownames(report) <- NULL
  report
}

build_combined_metadata <- function(spatialecotyper_list) {
  metadata_list <- lapply(names(spatialecotyper_list), function(sample_id) {
    metadata <- spatialecotyper_list[[sample_id]]$metadata
    if (!all(c("X", "Y", "CellType", "SE") %in% colnames(metadata))) {
      stop("Metadata for sample ", sample_id, " must contain X, Y, CellType and SE columns.")
    }

    metadata$sample <- sample_id
    metadata$local_SE <- metadata$SE
    metadata$sample_local_SE <- paste(sample_id, metadata$SE, sep = "__")
    metadata
  })

  combined <- do.call(rbind, metadata_list)
  rownames(combined) <- unlist(lapply(metadata_list, rownames), use.names = FALSE)
  combined
}

summarize_local_se <- function(combined_metadata) {
  tab <- as.data.frame(
    table(
      sample = combined_metadata$sample,
      local_SE = combined_metadata$local_SE,
      CellType = combined_metadata$CellType,
      useNA = "no"
    ),
    stringsAsFactors = FALSE
  )
  tab <- tab[tab$Freq > 0, , drop = FALSE]

  totals <- aggregate(Freq ~ sample + local_SE, tab, sum)
  colnames(totals)[colnames(totals) == "Freq"] <- "n_cells_in_local_SE"
  tab <- merge(tab, totals, by = c("sample", "local_SE"), all.x = TRUE)
  tab$celltype_fraction <- tab$Freq / tab$n_cells_in_local_SE
  tab[order(tab$sample, tab$local_SE, -tab$celltype_fraction), , drop = FALSE]
}

diagnose_subresolution_candidates <- function(spatialecotyper_list, candidates, out_dir) {
  all_summary <- list()
  all_celltypes <- list()
  all_failures <- list()

  for (res in candidates) {
    message_time("Diagnose subresolution: ", res)
    sample_rows <- list()
    celltype_rows <- list()

    for (sample_id in names(spatialecotyper_list)) {
      result <- spatialecotyper_list[[sample_id]]
      diagnosis <- tryCatch({
        obj <- Seurat::FindClusters(result$obj, resolution = res, verbose = FALSE)
        obj$SE <- paste0("InitSE", obj$seurat_clusters)
        metadata <- SpatialEcoTyper::AnnotateCells(result$metadata, obj)
        metadata <- metadata[!is.na(metadata$SE), , drop = FALSE]

        if (nrow(metadata) == 0) {
          stop("No cells remain after AnnotateCells().", call. = FALSE)
        }

        celltype_initse <- unique(metadata[, c("CellType", "SE"), drop = FALSE])
        celltype_counts <- aggregate(
          SE ~ CellType,
          celltype_initse,
          function(x) length(unique(x))
        )
        colnames(celltype_counts)[colnames(celltype_counts) == "SE"] <- "n_initse"
        celltype_counts$sample <- sample_id
        celltype_counts$subresolution <- res

        list(
          sample = data.frame(
            subresolution = res,
            sample = sample_id,
            status = "ok",
            n_cells = nrow(metadata),
            n_initse = length(unique(metadata$SE)),
            n_celltypes = length(unique(metadata$CellType)),
            error = NA_character_,
            stringsAsFactors = FALSE
          ),
          celltype = celltype_counts
        )
      }, error = function(e) {
        list(
          sample = data.frame(
            subresolution = res,
            sample = sample_id,
            status = "failed",
            n_cells = NA_integer_,
            n_initse = NA_integer_,
            n_celltypes = NA_integer_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          ),
          celltype = NULL
        )
      })

      sample_rows[[sample_id]] <- diagnosis$sample
      if (!is.null(diagnosis$celltype)) {
        celltype_rows[[sample_id]] <- diagnosis$celltype
      }
    }

    sample_df <- do.call(rbind, sample_rows)
    all_failures[[as.character(res)]] <- sample_df[sample_df$status != "ok", , drop = FALSE]

    celltype_df <- if (length(celltype_rows) > 0) {
      do.call(rbind, celltype_rows)
    } else {
      data.frame()
    }

    eligible_celltypes <- 0L
    if (nrow(celltype_df) > 0) {
      ct_split <- split(celltype_df, celltype_df$CellType)
      eligible_celltypes <- sum(vapply(ct_split, function(ct_df) {
        nrow(ct_df) >= 2 && min(ct_df$n_initse) >= 3
      }, logical(1)))
      all_celltypes[[as.character(res)]] <- celltype_df
    }

    ok_initse <- sample_df$n_initse[sample_df$status == "ok"]
    all_summary[[as.character(res)]] <- data.frame(
      subresolution = res,
      n_samples = nrow(sample_df),
      n_failed_samples = sum(sample_df$status != "ok"),
      min_initse = if (length(ok_initse) > 0) min(ok_initse) else NA_integer_,
      median_initse = if (length(ok_initse) > 0) median(ok_initse) else NA_real_,
      max_initse = if (length(ok_initse) > 0) max(ok_initse) else NA_integer_,
      eligible_celltypes_by_initse_count = eligible_celltypes,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, all_summary)
  sample_df <- do.call(rbind, lapply(candidates, function(res) {
    key <- as.character(res)
    failures <- all_failures[[key]]
    if (is.null(failures)) {
      return(NULL)
    }
    failures
  }))
  celltype_df <- if (length(all_celltypes) > 0) {
    do.call(rbind, all_celltypes)
  } else {
    data.frame()
  }

  write.csv(summary_df, file.path(out_dir, "subresolution_diagnosis_summary.csv"), row.names = FALSE, quote = TRUE)
  if (!is.null(sample_df) && nrow(sample_df) > 0) {
    write.csv(sample_df, file.path(out_dir, "subresolution_diagnosis_failures.csv"), row.names = FALSE, quote = TRUE)
  }
  if (!is.null(celltype_df) && nrow(celltype_df) > 0) {
    write.csv(celltype_df, file.path(out_dir, "subresolution_diagnosis_celltype_initse_counts.csv"), row.names = FALSE, quote = TRUE)
  }

  summary_df
}

write_run_config <- function(path, config) {
  lines <- unlist(lapply(names(config), function(key) {
    value <- config[[key]]
    if (length(value) > 1) {
      value <- paste(value, collapse = ",")
    }
    paste0(key, ": ", value)
  }))
  writeLines(lines, path)
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))

root_dir <- if (!is.null(opts$root)) opts$root else get_script_dir()
root_dir <- normalizePath(root_dir, mustWork = TRUE)

out_dir <- if (!is.null(opts$outdir)) {
  opts$outdir
} else {
  file.path(root_dir, "integrated_SE")
}
out_dir <- normalizePath(out_dir, mustWork = FALSE)

default_input_qs <- file.path(dirname(root_dir), "adata_cn.qs2")
input_qs <- if (!is.null(opts[["input-qs"]])) {
  opts[["input-qs"]]
} else if (file.exists(default_input_qs)) {
  default_input_qs
} else {
  NULL
}
if (!is.null(input_qs)) {
  input_qs <- normalizePath(input_qs, mustWork = TRUE)
}
sample_col <- if (!is.null(opts[["sample-col"]])) opts[["sample-col"]] else "sample"

sample_include <- NULL
if (!is.null(opts$samples) && nzchar(opts$samples)) {
  sample_include <- strsplit(opts$samples, ",", fixed = TRUE)[[1]]
}

normalization_method <- if (!is.null(opts$normalization)) opts$normalization else "None"
nmf_ranks <- as_int_vector(opts[["nmf-ranks"]], default = 10L)
nrun_per_rank <- if (!is.null(opts[["nrun-per-rank"]])) as.integer(opts[["nrun-per-rank"]]) else 30L
min_coph <- if (!is.null(opts[["min-coph"]])) as.numeric(opts[["min-coph"]]) else 0.95
nfeatures <- if (!is.null(opts$nfeatures)) as.integer(opts$nfeatures) else 300L
min_features <- if (!is.null(opts[["min-features"]])) as.integer(opts[["min-features"]]) else 10L
region <- if (!is.null(opts$region) && nzchar(opts$region)) opts$region else NULL
downsample_by_region <- if (!is.null(opts[["downsample-by-region"]])) {
  as_bool(opts[["downsample-by-region"]])
} else {
  TRUE
}
subresolution <- if (!is.null(opts$subresolution)) as.numeric(opts$subresolution) else 30
minibatch <- if (!is.null(opts$minibatch)) as.integer(opts$minibatch) else 5000L
ncores <- if (!is.null(opts$ncores)) as.integer(opts$ncores) else 15L
seed <- if (!is.null(opts$seed)) as.integer(opts$seed) else 1L
dry_run <- "dry-run" %in% opts$flags || as_bool(opts[["dry-run"]])
diagnose_subresolution <- if (!is.null(opts[["diagnose-subresolution"]])) {
  as_numeric_vector(opts[["diagnose-subresolution"]], default = numeric(0))
} else {
  numeric(0)
}

message_time("Root directory: ", root_dir)
message_time("Output directory: ", out_dir)

sample_map <- NULL
if (!is.null(input_qs)) {
  message_time("Read sample list from ", input_qs, " column ", sample_col)
  sample_map <- read_sample_map_from_qs(
    input_qs = input_qs,
    sample_col = sample_col,
    sample_include = sample_include
  )
  message_time("Samples from metadata: ", nrow(sample_map))
} else {
  message_time("No input QS found/provided; fall back to scanning sample directories.")
}

manifest <- discover_single_sample_results(
  root_dir = root_dir,
  out_dir = out_dir,
  sample_include = sample_include,
  sample_map = sample_map
)
validate_manifest(manifest)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(out_dir, "single_sample_manifest.csv")
discovered_manifest_path <- file.path(out_dir, "single_sample_manifest_discovered.csv")
write.csv(manifest, discovered_manifest_path, row.names = FALSE, quote = TRUE)

message_time("Discovered ", nrow(manifest), " samples.")
message_time("Discovered manifest written: ", discovered_manifest_path)

if (dry_run) {
  write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)
  message_time("Manifest written: ", manifest_path)
  print(manifest[, c("sample", "sample_safe", "result_rds", "scdata_qs"), drop = FALSE])
  message_time("Dry run finished. Remove --dry-run to run IntegrateSpatialEcoTyper().")
  quit(save = "no", status = 0)
}

message_time("Read single-sample SpatialEcoTyper results.")
spatialecotyper_list <- setNames(
  lapply(seq_len(nrow(manifest)), function(i) {
    read_spatialecotyper_result(manifest$result_rds[[i]], manifest$sample[[i]])
  }),
  manifest$sample
)

message_time("Precheck samples with Seurat::FindClusters().")
findclusters_report <- check_findclusters_eligible_samples(
  spatialecotyper_list = spatialecotyper_list,
  manifest = manifest,
  subresolution = subresolution
)
findclusters_report_path <- file.path(out_dir, "findclusters_sample_filter_report.csv")
write.csv(findclusters_report, findclusters_report_path, row.names = FALSE, quote = TRUE)
message_time("FindClusters filter report written: ", findclusters_report_path)

excluded_samples <- findclusters_report[!findclusters_report$keep, , drop = FALSE]
excluded_samples_path <- file.path(out_dir, "findclusters_excluded_samples.csv")
if (nrow(excluded_samples) > 0) {
  write.csv(excluded_samples, excluded_samples_path, row.names = FALSE, quote = TRUE)
  message_time("Excluded ", nrow(excluded_samples), " samples that failed FindClusters().")
  message_time("Excluded sample report written: ", excluded_samples_path)
} else if (file.exists(excluded_samples_path)) {
  unlink(excluded_samples_path)
}

keep_samples <- findclusters_report$sample[findclusters_report$keep]
manifest <- manifest[manifest$sample %in% keep_samples, , drop = FALSE]
manifest <- manifest[match(keep_samples, manifest$sample), , drop = FALSE]
rownames(manifest) <- NULL
spatialecotyper_list <- spatialecotyper_list[manifest$sample]

if (nrow(manifest) < 2) {
  stop(
    "Fewer than 2 samples remain after FindClusters precheck; ",
    "cannot run IntegrateSpatialEcoTyper(). Check: ",
    findclusters_report_path
  )
}

write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)
message_time("Manifest written after FindClusters filtering: ", manifest_path)
message_time("Samples kept for integration: ", nrow(manifest), " / ", nrow(findclusters_report))

if (length(diagnose_subresolution) > 0) {
  message_time("Run subresolution diagnosis and exit before integration.")
  diagnosis_summary <- diagnose_subresolution_candidates(
    spatialecotyper_list = spatialecotyper_list,
    candidates = diagnose_subresolution,
    out_dir = out_dir
  )
  print(diagnosis_summary)
  message_time("Subresolution diagnosis finished. Check outputs under: ", out_dir)
  quit(save = "no", status = 0)
}

message_time("Read single-sample expression matrices.")
data_list <- setNames(
  lapply(manifest$scdata_qs, qs2::qs_read),
  manifest$sample
)

combined_metadata <- build_combined_metadata(spatialecotyper_list)
combined_metadata_path <- file.path(out_dir, "single_sample_local_SE_metadata.qs2")
qs2::qs_save(combined_metadata, combined_metadata_path)

local_se_summary <- summarize_local_se(combined_metadata)
local_se_summary_path <- file.path(out_dir, "single_sample_local_SE_celltype_fraction.csv")
write.csv(local_se_summary, local_se_summary_path, row.names = FALSE, quote = TRUE)

config <- list(
  root_dir = root_dir,
  out_dir = out_dir,
  input_qs = if (is.null(input_qs)) "NULL" else input_qs,
  sample_col = sample_col,
  n_samples = nrow(manifest),
  samples = manifest$sample,
  normalization_method = normalization_method,
  nmf_ranks = nmf_ranks,
  nrun_per_rank = nrun_per_rank,
  min_coph = min_coph,
  nfeatures = nfeatures,
  min_features = min_features,
  region = if (is.null(region)) "NULL" else region,
  downsample_by_region = downsample_by_region,
  subresolution = subresolution,
  minibatch = minibatch,
  ncores = ncores,
  seed = seed
)
write_run_config(file.path(out_dir, "integration_run_config.txt"), config)

message_time("Start IntegrateSpatialEcoTyper().")
integration_result <- IntegrateSpatialEcoTyper(
  SpatialEcoTyper_list = spatialecotyper_list,
  data_list = data_list,
  outdir = out_dir,
  normalization.method = normalization_method,
  nmf_ranks = nmf_ranks,
  nrun.per.rank = nrun_per_rank,
  min.coph = min_coph,
  nfeatures = nfeatures,
  min.features = min_features,
  Region = region,
  downsample.by.region = downsample_by_region,
  subresolution = subresolution,
  minibatch = minibatch,
  ncores = ncores,
  seed = seed
)

if (!is.null(integration_result)) {
  saveRDS(integration_result, file.path(out_dir, "IntegrateSpatialEcoTyper_return.rds"))
}

message_time("Integration finished. Check outputs under: ", out_dir)
