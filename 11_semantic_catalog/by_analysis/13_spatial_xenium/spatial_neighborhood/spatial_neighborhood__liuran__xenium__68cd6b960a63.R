
# 1. Xenium function ------------------------------------------------------
# TODO: 算法缺陷在于赢者通吃
# 第一名的基因所在的细胞评分非常重，基础让其他细胞无法超过
# 最好可以达到如下效果：
# 多个老二可以干掉老一
# 多个老末永远也干不掉老一
annotate_cluster <- function(x, panel, adjust_length = FALSE){ # browser()
  # x是一组按照顺序排列的基因list
  # 一般是来源于差异表达之后，按照p值或者avg_log2fc排列的
  
  # panel是一个基因的注释数据list
  # name是细胞名，内容是基因向量
  
  # filter panel
  all_genes <- unlist(panel)
  inter_genes <- intersect(x, all_genes)
  
  panel_filter <- panel %>% map(~ .x[.x %in% inter_genes])
  
  # if(any(map_dbl(panel_filter, length) < 1)){
  #   stop("some celltype cannot to be annotation, causing none of marker genes found in x")
  # }
  
  x <- x[x %in% inter_genes]
  
  if(length(x) < 1) stop("none gene found in marker gene list")
  
  # score
  rank <- 1:length(x)
  # score_x <- 1/rank   # 用排名的倒数作为raw score
  # score_x <- 1/(rank + 1)^2 # 用排名的平方倒数作为raw score
  score_x <- 1/( (rank + 1) * (rank + 2) )# 用排名的倒数作为raw score
  
  # panel score
  panel_score <-
    panel_filter %>%
    map_dbl(function(g){
      idx  <- x %in% g
      res  <- sum(score_x[idx])
      
      if(adjust_length && length(g) > 1){
        # 修正因marker gene多则评分高的问题
        return(res/length(g))
      }
      
      res
    })
  
  return(panel_score)
  
}

# anno score
annotate_Q_val <- function(x){
  # x: annotation score
  # names(x): cell
  if(is.null(names(x))) stop("x must have names")
  
  # Q_val，来源于两部分
  # 第一是A/sum(x)
  # 第二是A/B - 1
  stopifnot(length(x) > 0)
  if(length(x) == 1) return(list(
    'A/sum(x)' = x,
    'A/B' = x, 
    warnning = "Single",
    first_cell = names(x)[1]
  ))
  
  x <- x %>% sort(decreasing = TRUE)
  
  c1 <-  x[1]/sum(x)
  c2 <-  x[1]/x[2] - 1
  
  q_val <- list(
    'A/sum(x)'   = c1,
    'delta(A/B)' = c2, 
    warnning = case_when(
      c1 < 0.1  || c2 < 0.1  ~ "Error",
      c1 < 0.2  || c2 < 0.25 ~ "Warnning",
      c1 < 0.3  || c2 < 0.5  ~ "Info", 
      TRUE ~ "Cool"
    ),
    first_cell  = names(x)[1],
    second_cell = names(x)[2],
    first_cell_score = x[1],
    second_cell_score = x[2]
  )
  
  return(q_val)
}

#' used to check Xenium project name 
#' 
check_xenium_pattern <-
  function(x, 
           pattern = "^output", 
           pattern2 = ".*?__.*?__.*?__.*?__.*?",
           use_basename = TRUE){
    
    stopifnot(length(x) > 0)
    
    if(use_basename){
      x <- basename(x)
    }
    
    stopifnot(length(pattern)  < 2)
    stopifnot(length(pattern2) < 2)
    
    if(length(pattern) == 1){
      idx <- x %>% str_detect(pattern = pattern)
    }else{
      idx <- TRUE
    }
    
    if(length(pattern2) == 1){
      idx2 <- x %>% str_detect(pattern = pattern2)
    }else{
      idx2 <- TRUE
    }
    
    # return
    return(
      idx & idx2
    )
    
  }

#' tool function used to extract character from xenium-output-file-name
extract_pattern <- function(
    x, 
    split_pattern = "__", 
    extract_postion =  c(3)){
  
  stopifnot(is.character(x))
  stopifnot(length(extract_postion) > 0)
  stopifnot(is.numeric(extract_postion))
  
  
  if(any(!str_detect(x, split_pattern))) stop(str_glue("Not found pattern {split_pattern} in one or more of {x}"))
  
  xx <- x %>% str_split(split_pattern)
  xx_max_len <- xx %>% map_dbl(length) %>% max
  
  
  if(max(extract_postion) > xx_max_len){
    stop("extract_postion must be wrong, causing maximum value of extract_postion bigger than split charater.")
  }
  
  res <- xx %>% map(~ .[extract_postion]) %>% map_chr(str_c, collapse = split_pattern)
  
  res
  
}


#' 解析Xenium project文件夹
#'
# parse project name
#' parse project name from 'data/*'
#' check unique proj name
#' <return> a data.frame
parseProjName_xenium <- 
  function(x,
           pattern = "^output.*",
           pattern2 = ".*?__.*?__.*?__.*?__.*?",
           split_pattern = "__", 
           extract_postion =  c(3),
           check.unique = TRUE) { 
    
    xenium_out_idx <- check_xenium_pattern(
      x, 
      pattern = pattern, 
      pattern2 = pattern2,
      use_basename = TRUE
    )
    
    if(sum(xenium_out_idx) < 1){
      stop(str_glue("No output diretory match the pattern: {pattern} and pattern: {pattern2}"))
    }
    
    x <- x[xenium_out_idx]
    
    
    stopifnot(length(x) > 0) 
    
    df <- tibble::tibble(.path = x)
    df <- df %>% 
      mutate(
        .dir = dirname(.path),
        .base = basename(.path),
        .proj = str_remove_all(.base, pattern)
      ) %>% 
      mutate(
        .proj = str_remove_all(.proj, "[._-]?$"),
        # .proj = .proj %>% str_split("__") %>% map_chr(3)
        .proj = .proj %>% extract_pattern(split_pattern = split_pattern, extract_postion = extract_postion)
      ) %>%
      dplyr::filter(
        str_detect(.base, pattern = pattern)
      )
    
    # if proj duplicated, send a message
    is_duplicated <- duplicated(df$.proj)
    if(sum(is_duplicated) > 0){
      message("Not all project names be uniquely identified...")
    }
    
    if(check.unique && sum(is_duplicated) > 0){
      
      message("Duplicated project name(s) set as unique label.")
      df$.proj <- make.names(df$.proj, unique = TRUE)
    }
    
    # return
    df <- as.data.frame(df)
    res_tryCatch <- 
      tryCatch({
        rownames(df) <- df$.proj
        "work"
      },
      error = function(x) return("error")
      )
    if(res_tryCatch == "error"){
      rownames(df) <- make.names(df$.path, unique = TRUE)
    }
    return(df)
  }

#' 解析ROI Selection文件夹
#' 

parseSelectionName <-
  function(x,
           filter_proj  = NULL,
           check.unique = TRUE) {
    # stopifnot(length(x) > 0) 
    if(length(x) < 1) return(data.frame())
    
    if(!is.null(filter_proj) && !is.character(filter_proj)){
      stop("filter_proj should only be character if not NULL")
    }
    
    df <- tibble::tibble(.path = x)
    df <- df %>%
      mutate(
        .dir = dirname(.path),
        .base = basename(.path),
        .proj = .base %>% str_remove_all("\\..*$") %>% str_remove_all("[._-]?$")
      ) 
    
    if(check.unique && any(duplicated(df$.proj))){
      stop(
        "Found dupicated project."
      )
    }
    
    if(!is.null(filter_proj)){
      df <- df %>%  dplyr::filter(.proj %in% filter_proj)
    }
    df
  }

# parseSelectionName <- 
#   function(x,
#            anno_pattern = "annotation.*",
#            check.unique = TRUE) {
#     
#     stopifnot(length(x) > 0) 
#     
#     df <- tibble::tibble(.path = x)
#     df_ls <- df %>% 
#       mutate(
#         .dir = dirname(.path),
#         .base = basename(.path),
#         .proj = .base %>% str_remove_all("\\..*$")
#       ) %>%
#       mutate(.grp = ifelse(str_detect(.proj, anno_pattern), "Annotation", "Selection")) %>%
#       
#       mutate(
#         .proj = .proj %>% str_remove_all(anno_pattern),
#         .proj = str_remove_all(.proj, "[._-]?$")
#       ) %>%
#       group_by(.grp) %>%
#       nest() %>%
#       deframe()
#     
#     # if proj duplicated, send a message
#     for(idx in seq_along(df_ls)){
#       
#       is_duplicated <- duplicated(df_ls[[idx]]$.proj)
#       if(sum(is_duplicated) > 0){
#         message("Not all project names be uniquely identified...")
#       }
#       
#       if(check.unique && sum(is_duplicated) > 0){
#         
#         message("Duplicated project name(s) set as unique label.")
#         df_ls[[idx]]$.proj <- make.names(df_ls[[idx]]$.proj, unique = TRUE)
#       }
#       
#       df_ls[[idx]] <- as.data.frame(df_ls[[idx]])
#       res_tryCatch <- 
#         tryCatch({
#           rownames(df_ls[[idx]]) <- df_ls[[idx]]$.proj
#           "work"
#         },
#         error = function(x) return("error")
#         )
#       if(res_tryCatch == "error"){
#         rownames(df_ls[[idx]]) <- make.names(df_ls[[idx]]$.path, unique = TRUE)
#       }
#     }
#     
#     # return
#     return(df_ls)
#   }
# 
# if(F){
#   parseSelectionName(list.files("data/ROI Selection/"))
# }

#' check all ROI Selection file have a annotation file
#' 


# 2. singlecell helping function  -----------------------------------------


# parse project name
#' parse project name from 'data/*'
#' check unique proj name
#' <return> a data.frame
parseProjName <- 
  function(x,
           pattern = "matrix.*",
           check.unique = TRUE) {
    
    stopifnot(length(x) > 0) 
    
    df <- tibble::tibble(.path = x)
    df <- df %>% 
      mutate(
        .dir = dirname(.path),
        .base = basename(.path),
        .proj = str_remove_all(.base, pattern)
      ) %>% 
      mutate(
        .proj = str_remove_all(.proj, "[._-]?$")
      ) %>%
      dplyr::filter(
        str_detect(.base, pattern = pattern)
      )
    
    # if proj duplicated, send a message
    is_duplicated <- duplicated(df$.proj)
    if(sum(is_duplicated) > 0){
      message("Not all project names be uniquely identified...")
    }
    
    if(check.unique && sum(is_duplicated) > 0){
      
      message("Duplicated project name(s) set as unique label.")
      df$.proj <- make.names(df$.proj, unique = TRUE)
    }
    
    # return
    df <- as.data.frame(df)
    res_tryCatch <- 
      tryCatch({
        rownames(df) <- df$.proj
        "work"
      },
      error = function(x) return("error")
      )
    if(res_tryCatch == "error"){
      rownames(df) <- make.names(df$.path, unique = TRUE)
    }
    return(df)
  }

# test 
if(F){
  dir("data") %>% parseProjName
  dir("data", full.names = TRUE) %>% parseProjName
  dir("data", full.names = TRUE) %>% c("data/C1A_matrix_") %>% parseProjName
  "data/I1A_mat" %>% parseProjName()
  rep("data/I1A_mat", 2) %>% parseProjName(pattern = "mat", check.unique = FALSE)
}


# merger seurat object
#' merge object 
#' add some condition validations
mergeSeurat <- function(x){
  if(!"list" %in% class(x)) stop("x can only be a list of SeuratObject.") 
  
  len <- length(x)
  if(len == 0) stop("No data found in x.")
  if(len == 1) return(x[[1]])
  
  merge(x[[1]], x[-1], add.cell.ids = names(x))
  
}

# findReference
#'
findReference <- function(obj = NULL, name){
  
  if(is.numeric(name)){
    return(name)
  }else{
    stopifnot(!is.null(obj))
    stopifnot(!is.null(names(obj)))
    name_idx <- match(name, names(obj), nomatch = NA)
    
    if(any(is.na(name_idx))){
      message(str_glue("{glue::glue_collapse(name[is.na(name_idx)])} not found..."))
      return(name_idx[!is.na(name_idx)])
    }else{
      return(name_idx)
    }
  }
}

if(F){
  findReference(list(), 1)
  findReference(list(a = 1), "a")
  findReference(list(a = 1, b = 2), c("a", "b"))
  findReference(list(a = 1, b = 2), c("a", "b", "c"))
  
}


# get assay name
getAssayName <- function(obj){
  return(
    names(obj@assays)
  )
}


# Find all markers using FindConservedMarkers
FindAllConservedMarkers <-
  function(obj, 
           ident.1 = NULL, 
           grouping.var,
           ...) {
    
    stopifnot(is(obj, "Seurat"))
    
    all_ident <- unique(Idents(obj))
    stopifnot(length(all_ident) > 1)
    
    
    if(is.null(ident.1)){
      message("ident.1 is null, will find all markers of each ident in all idents.")
    }else{
      all_ident <-  ident.1
    }
    
    all_ident %>%
      set_names(., .) %>%
      map(function(ident){
        FindConservedMarkers(
          object = obj,
          ident.1 = ident,
          grouping.var = grouping.var,
          ...
        )
      })
  }


# rewrite RenameIdents supported !!!
renameIdents <- function(object, 
                         name, 
                         old_ident = "old_idents_before_annotation",
                         force = FALSE,
                         ...){
  if(...length() > 0){
    "Only support object and name, other arguments will  be ignored..."
  }
  
  if(is.null(names(name))) stop("Please provide name like: c(old_ident = new_ident)")
  if(any(is.na(name))) {
    message("NA found in name will be ignored.")
    name <- na[!is.na(name)]
  }
  
  if(length(name) < 1) {
    warning("name you provided is blank, do nothing... ")
    return(object)
  }
  
  
  # stack old idents
  if(old_ident %in% colnames(object[[]])){
    if(!force){
      warning("old_ident name have been existed, re-write it using 'force = TRUE'")
    }
  }
  object <- StashIdent(object, save.name = old_ident)
  
  ex <- expr(
    RenameIdents(object, !!!name)
  )
  eval(ex)
}

# tool class for result of  parsed mito_qc criteria
# class: parsed_criteria
print.parsed_criteria <- function(x, ...){
  cat("parsed_criteria: ")
  x <- unclass(x)
  cat(str_c(x, collapse = " "))
  cat("\n")
}

# parse mito_qc criteria
# only support >, <, =
# x: string use as input
# values: a list start with a function,  then argument
# example:
# '> 200' would be parsed as: list('>', 200)


parse_qc_criteria <- function(x){
  
  # x first by trimed, also will be transformed a string if not
  orign_x <- x
  x <- str_trim(x)
  
  if(str_length(x) < 1){
    stop("x cannot be a  null character...")
  }
  
  # string split
  operator <- str_sub(x, end = 1)
  number   <- str_sub(x, start = 2) 
  if(str_trim(str_to_lower(number)) == "inf") {
    number <- Inf
  }else if(str_trim(str_to_lower(number)) %in% c("-inf", "- inf", "-  inf")) {
    number <- -Inf
  }else{
    number <- parse_number(number)
  }
  
  if(is.na(number)){
    stop(str_glue("Error number style after operator, your input: '{orign_x}'"))
  }
  
  # parse 
  res <- switch (operator,
                 ">" = list(">", number),
                 "<" = list("<", number),
                 "=" = list("==", number)
  )
  
  return(
    structure(res, class = "parsed_criteria")
  )
}

# example
if(F){
  local({
    print(test <- parse_qc_criteria("> 200"))
    
    (do.call(
      test[[1]], 
      c(1, test[-1])
    ))
  })
  
  parse_qc_criteria("> 200a")
  parse_qc_criteria("= a200a")
  parse_qc_criteria(">0")
  parse_qc_criteria("> ")
  parse_qc_criteria("> Inf")
  parse_qc_criteria("> inf")
  parse_qc_criteria("> --inf")
  parse_qc_criteria("> - inf")
}

# tools function to eval parsed_criteria 
# parsed_criteria record two things: (1) fucntion name; (2) the argument
# this fucntion also can add one more argument to eval the function

eval_parsed_criteria  <-
  function(x, 
           others = NULL,
           position = c("before", "after"),
           verbose = TRUE
  ){
    stopifnot("parsed_criteria" %in% class(x) )
    
    position <- match.arg(position)
    
    # extract data
    func <- x[[1]]
    argument <- x[-1]
    
    
    # add any argument?
    if(!is.null(others)){
      argument <- switch(
        position,
        before = c(list(others), argument),
        after = c(argument, list(others))
      )
    }
    # do.call
    if(verbose){
      cat("Eval expression: \n  ")
      print(rlang::call2(func, !!!argument))
    }
    
    do.call(func, argument)
    
  }

if(F){
  parse_qc_criteria("> 200a") %>%
    eval_parsed_criteria(others = 1)
  
  parse_qc_criteria("> 200a") %>%
    eval_parsed_criteria(others = 1:5)
}

# CellCycleScoring函数的nbin参数设置需要仔细筛选，否则会出错
# 使用一个递归的函数找到合适的nbin值，每次递归nbin减1，
# 如果nbin小于1，则直接返回原始对象，不再进行CellCycleScoring分析

reduceCellCycleScoring <- function(x, nbin = 30, verbose = TRUE, ...){
  if(nbin < 1){
    message("nbin < 1, return orginal object")
    return(x)
  }
  
  res <- tryCatch({
    message(str_glue("test nbin {nbin}" ))
    CellCycleScoring(x, nbin = nbin, ...)
  },
  error = function(e){
    message(str_glue("nbin {nbin} failed, test {nbin-1}" ))
    reduceCellCycleScoring(x, nbin = nbin - 1, ...)
  }
  )
  
  return(res)
}

# FindCluster函数的resolution决定了获得的分群的数量
# 使用递归函数获取一个相对合理的resolution值，要求如下：
# cluster number: 5 - 15之间
# resolution: 默认是0.8，其范围默认是0.1-2，0.1的间隔去尝试
# resolution_range: 默认0.1 - 2之间
# 会进行递归尝试，如果没有最恰当的结果：
# return_best_if_none_fit为TRUE代表返回最后一个最合适的结果，为FALSE值函数直接stop终止
reduceFindClusters <- function(
    x, 
    resolution = 0.8, 
    resolution_range = c(0.1, 2), 
    cluster_range = c(5, 15),
    verbose = TRUE, 
    return_best_if_none_fit = TRUE,
    ...
){
  if(!is.numeric(resolution_range)){
    warning("resolution_range must be numeric, use defatult range: 0.1 - 2")
    resolution_range <- c(0.1, 2)
  }
  stopifnot(all(resolution_range) > 0)
  stopifnot(length(resolution_range) <= 2)
  
  reso_seq <- seq(min(resolution_range), max(resolution_range), by = 0.1)
  if(resolution > max(reso_seq) || resolution < min(reso_seq)){
    warning("resolution not in resolution_range, use defatult value: 0.8")
    resolution <- 0.8
  }
  
  if(!is.numeric(cluster_range)){
    warning("cluster_range must be numeric, use defatult range: 5 - 15")
    resolution_range <- c(5, 15)
  }
  stopifnot(all(cluster_range > 0))
  stopifnot(length(cluster_range) <= 2)
  
  
  # 递归
  if(verbose){
    message(str_glue("\n### Test resolution {resolution} ###\n\n"))
    cat("\n")
  }
  res_cluster <-
    FindClusters(
      x,
      resolution = resolution,
      verbose = verbose,
      ...
    )
  res_cluster_number <- n_distinct(res_cluster[["seurat_clusters"]])
  
  # 
  
  # msg for output
  msg1 <- str_glue("None of the resolution in {paste(resolution_range, collapse ='-')}")
  msg2 <- str_glue("fit the cluster criteria: {paste(cluster_range, collapse ='-')}")
  msg3 <- str_glue("return the best cluster result, number of cluster: {res_cluster_number}")
  msg4 <- str_glue("Number of clusters: {res_cluster_number} not {msg2}")
  msg5 <- str_glue("Number of clusters: {res_cluster_number} {msg2}")
  
  # 核心代码，用于进行递归逻辑处理
  if(res_cluster_number >= min(cluster_range) &&
     res_cluster_number <= max(cluster_range)){
    #### 结果OK直接返回结果 ###
    message(msg5)
    cat("\n\n")
    return(res_cluster)
    
  }else if(res_cluster_number < min(cluster_range)){
    ### cluster太少，增加resolution ###
    reso_seq_final <- reso_seq[reso_seq > resolution]
    if(length(reso_seq_final) < 1){
      # 没有需要进行测试的resolution值，表明均不符合要求，但是最贴近cluster标准的结果就是res_cluster
      
      if(return_best_if_none_fit){
        if(verbose){
          message(paste(c(msg1, msg2, msg3), collapse = " "))
        }
        return(res_cluster)
      }else{
        stop(paste(c(msg1, msg2), collapse = " "))
      }
    }else{
      if(verbose){
        message(msg4)
        cat("\n\n")
      }
    }
    
    # 继续尝试
    reduceFindClusters(
      res_cluster, 
      resolution = reso_seq_final[1], 
      resolution_range = range(reso_seq_final), 
      cluster_range = cluster_range,
      verbose = verbose, 
      return_best_if_none_fit = return_best_if_none_fit,
      ...
    )
  }else{
    ### cluster太大，减少resolution ###
    
    # 降序查找
    reso_seq_final <- reso_seq[reso_seq < resolution] %>% rev()
    if(length(reso_seq_final) < 1){
      # 没有需要进行测试的resolution值，表明均不符合要求，但是最贴近cluster标准的结果就是上面的res_cluster
      if(return_best_if_none_fit){
        if(verbose){
          message(paste(c(msg1, msg2, msg3), collapse = " "))
        }
        return(res_cluster)
      }else{
        stop(paste(c(msg1, msg2), collapse = " "))
      }
    }else{
      if(verbose){
        message(msg4)
        cat("\n\n")
      }
    }
    
    # 继续尝试
    reduceFindClusters(
      res_cluster, 
      resolution = reso_seq_final[1], 
      resolution_range = range(reso_seq_final), 
      cluster_range = cluster_range,
      verbose = verbose, 
      return_best_if_none_fit = return_best_if_none_fit,
      ...
    )
  }
}

if(F){
  test_seurat <- pbmc_small %>% 
    SCTransform() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30) %>% 
    FindNeighbors(dims = 1:30)
  test_seurat %>% reduceFindClusters(verbose = T)
}


# sort monocle node
# 形如Y_1 Y_34
sort_node <- function(x, pattern = "Y_", decreasing = FALSE, ...){
  if(!all(str_detect(x, pattern))){
    stop("Not all input are legal node name: Y_1 Y_2 ...")
  }
  
  df <- data.frame(x = x)
  df$order <- as.numeric(str_remove(df$x, pattern = pattern))
  if(decreasing){
    df <- df %>% arrange(desc(order)) 
  }else{
    df <- df %>% arrange(order)
  }
  
  df$x
}

if(F){
  sort_node(c('Y_2', "Y_1", "Y_01", "Y_11"))
  sort_node(c('Y_2', "Y_a"))
  sort_node(c('cluster_2', "cluster_12"), pattern = "cluster_")
  
}



# import Seurat for monocle2
seurat2monocle <- function(
    x,
    assay = "SCT", 
    slot = c("counts", "data", "scale.data"), 
    expressionFamily = uninormal(), 
    setOrderingGene = TRUE,
    ...
){
  
  slot <- match.arg(slot)
  
  if(is.null(assay)) assay <- "RNA"
  
  assayName <- getAssayName(x)
  stopifnot(assay %in% assayName)
  
  # data
  data <- slot(x[[assay]], slot)
  
  pdata <- x[[]]
  pd <- new("AnnotatedDataFrame", data = pdata)
  
  fdata <- data.frame(
    gene_short_name = row.names(data),
    row.names = row.names(data)
  )
  fd <- new("AnnotatedDataFrame", data = fdata)
  
  if(is.null(expressionFamily)){
    if(slot == "data"){
      expressionFamily <- uninormal()
    }else{
      expressionFamily <- negbinomial.size()
    }
  }
  
  
  # return
  y <- newCellDataSet(
    data,
    phenoData = pd,
    featureData = fd,
    expressionFamily = expressionFamily
  )
  
  if(setOrderingGene){
    if(assay == "SCT"){
      orderingGene <- rownames(slot(x[["SCT"]], "scale.data"))
      if(length(orderingGene) < 1){
        stop("No scale.data in SCT assay.")
      }
      y <- setOrderingFilter(y, orderingGene)
    }
  }
  
  return(y)
}
