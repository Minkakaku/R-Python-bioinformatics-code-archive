library(tidyverse)
library(Seurat)

# 1. parse path ------------------------------------------------------------
parsePath <- function(path, obj_label = "ObjectData", roi_label = "Phenotype"){

  # 解析所有包含OjectData的文件夹

  # 情形1. obj_label在path里面，则Slide_ID默认为MIPCC-001，Tissue_ID默认为Tissue_001，ROI_ID则为ObjectData中每个文件的记录label
  # 形如data/ObjectData
  # 情形2. obj_label在path的子文件夹中，则Slide_ID为子文件夹的名字，Tissue_ID默认为Tissue_001，ROI_ID则为ObjectData中每个文件的记录label
  # 形如data/MIPCC-DSP-C-20220809C-01/ObjectData
  # 情形3. obj_label在path的子文件夹的子文件夹中，则Slide_ID为子文件夹的名字，Tissue_ID则为子文件夹的子文件夹名字，ROI_ID则为ObjectData中每个文件的记录label
  # 形如data/MIPCC-DSP-C-20220809C-01/Halo archive 2022-10-21 10-26 - v3.4.2986/ObjectData

  # 解析情形1
  path_1 <- list.files(path = path, full.names = TRUE)
  if(length(path_1) < 1){
    stop("No files found, please check you celldive data path.")
  }

  obj_1 <- path_1 %>% str_subset(obj_label)
  roi_1 <- obj_1  %>% list.roi(roi_label = roi_label)

  if(length(obj_1) > 0 && length(roi_1) > 0){
    obj_df_1 <-
      tibble(
        obj_path = obj_1
      ) %>%
      mutate(
        Slide_ID  = "MIPCC_001",
        Tissue_ID = "Tissue_001"
      ) %>%
      mutate(ROI_ID = "default") %>%
      mutate(path = list(roi_1)) %>%
      unnest(path)
  }else{
    obj_df_1 <- NULL
  }

  # 解析情形2
  path_2 <- path_1 %>%
    setdiff(obj_1) %>%
    map(function(x) list.files(path = x, full.names = TRUE)) %>%
    keep(function(x) length(x) > 0)

  obj_2 <- path_2 %>%
    map(~ str_subset(.x, obj_label)) %>%
    keep(function(x) length(x) > 0)

  roi_2 <- obj_2 %>% map(list.roi, roi_label = roi_label) %>% keep(function(x) length(x) > 0)
  if(length(obj_2) > 0 && length(roi_2) > 0){
    obj_df_2 <-
      tibble(
        obj_path = obj_2
      ) %>%
      unnest(obj_path) %>%
      mutate(
        Slide_ID  = obj_path %>% dirname() %>% basename(),
        Tissue_ID = "Tissue_001"
      ) %>%
      mutate(ROI_ID = "default") %>%
      mutate(path = roi_2) %>%
      unnest(path)
  }else{
    obj_df_2 <- NULL
  }

  # 解析情形3
  path_3 <- path_2 %>%
    setdiff(obj_2) %>%
    map(function(x) list.files(path = x, full.names = TRUE)) %>%
    keep(function(x) length(x) > 0)

  obj_3 <- path_3 %>%
    map(~ str_subset(.x, obj_label)) %>%
    map(~ keep(.x, function(x) str_detect(basename(x), obj_label))) %>%
    keep(function(x) length(x) > 0)

  roi_3 <- obj_3 %>% unlist %>% map(list.roi, roi_label = roi_label) %>% keep(function(x) length(x) > 0)
  if(length(obj_3) > 0 && length(roi_3) > 0){
    obj_df_3 <-
      tibble(
        obj_path = obj_3
      ) %>%
      unnest(obj_path) %>%
      mutate(
        Slide_ID  = obj_path %>% dirname() %>% dirname() %>% basename(),
        Tissue_ID = obj_path %>% dirname() %>% basename()
      ) %>%
      mutate(ROI_ID = "default") %>%
      mutate(path = roi_3) %>%
      unnest(path)
  }else{
    obj_df_3 <- NULL
  }

  # return
  obj_df <-
    list(obj_df_1, obj_df_2, obj_df_3) %>%
    keep(function(x) length(x) > 0)

  if(length(obj_df) < 1){
    stop("No files found, please check you celldive data path.")
  }

  return(bind_rows(obj_df))

}

if(F){
  parsePath("data", roi_label = "Phenotype|PhenoType")
}


# 2. parse intensity ------------------------------------------------------

parseIntensity <- function(x, intensity_label = "Intensity"){
  # 默认Intensity为Cell Intensity： KI67-AF488 Cell Intensity
  # Membrane Intensity: CD3-CY5 Membrane Intensity
  # Cytoplasm Intensity： CD68-CY5 Cytoplasm Intensity
  # Nucleus Intensity： FOXP3-CY5 Nucleus Intensity

  stopifnot(is.data.frame(x))

  all_colnames <- colnames(x)

  cell_inten_idx      <- all_colnames %>% str_subset(str_glue("Cell {intensity_label}$"))
  membrane_inten_idx  <- all_colnames %>% str_subset(str_glue("Membrane {intensity_label}$"))
  cytoplasm_inten_idx <- all_colnames %>% str_subset(str_glue("Cytoplasm {intensity_label}$"))
  nucleus_inten_idx   <- all_colnames %>% str_subset(str_glue("Nucleus {intensity_label}$"))

  # all markers
  all_markers <- cell_inten_idx %>% str_remove(str_glue("Cell {intensity_label}$")) %>% str_trim()
  dup_idx <- duplicated(all_markers)

  if(any(dup_idx)){
    warning("Found duplicated marker...", immediate. = TRUE)
  }

  cell_inten <- x[, cell_inten_idx, drop=FALSE]
  membrane_inten <- x[, membrane_inten_idx, drop=FALSE]
  cytoplasm_inten <- x[, cytoplasm_inten_idx, drop=FALSE]
  nucleus_inten <- x[, nucleus_inten_idx, drop=FALSE]

  res <- list(
    intensity = cell_inten,
    membrane  = membrane_inten,
    cytoplasm = cytoplasm_inten,
    nucleus   = nucleus_inten
  ) %>%
    purrr::keep(function(x) all(dim(x) > 0))

  if(length(res) < 1){
    stop("Not found any intensity...")
  }

  return(res)
}

if(F){
  undebug(parseIntensity)

  all_data$dat[[1]] %>%
    parseIntensity()
}


# 3. parse segment --------------------------------------------------------
parseSegmentLable <- function(x, seg_label = "Classifier Label", return_all = FALSE){
  stopifnot(is.data.frame(x))

  if(!seg_label %in% colnames(x)){
    stop(str_glue("column {seg_label} not found"))
  }

  levs <- unique(x[[seg_label]])

  if(length(levs) < 1 || (length(levs) == 1 && str_trim(levs) == "")){
    # 空字符
    warning(
      "Not found segment or tissue classification",
      immediate. = TRUE
    )
  }

  if(return_all){
    return(x[[seg_label]])
  }

  return(levs)
}

# 3+ parse channel --------------------------------------------------------

parseChannel <- function(x, intensity_label = "Intensity"){
  stopifnot(is.data.frame(x))

  all_colnames <- colnames(x)

  # all markers
  cell_inten_idx <- all_colnames %>% str_subset(str_glue("Cell {intensity_label}$"))
  all_markers  <- cell_inten_idx %>% str_remove(str_glue("Cell {intensity_label}$")) %>% str_trim()
  dup_idx <- duplicated(all_markers)

  if(any(dup_idx)){
    warning("Found duplicated marker...", immediate. = TRUE)
  }

  all_markers
}

# 4. parse cutoff ---------------------------------------------------------
parseCutoff <- function(x, intensity_label = "Intensity", classification_label = "Classification"){
  stopifnot(is.data.frame(x))

  all_colnames <- colnames(x)

  # cell intensity
  # cell_inten_idx <- all_colnames %>% str_subset(str_glue("Cell {intensity_label}$"))
  # cell_intensity <- x[, cell_inten_idx, drop=FALSE]
  all_dat <- x

  # all markers
  cell_inten_idx <- all_colnames %>% str_subset(str_glue("Cell {intensity_label}$"))
  all_markers  <- cell_inten_idx %>% str_remove(str_glue("Cell {intensity_label}$")) %>% str_trim()
  dup_idx <- duplicated(all_markers)

  if(any(dup_idx)){
    warning("Found duplicated marker...", immediate. = TRUE)
  }

  # giving classification and intensity, get cutoff
  getCutoff <- function(class, intensity, name){
    levs <- unique(class)
    if(length(levs) < 1){
      stop("No classification found.")
    }

    if(length(levs) == 1){

      if(unique(levs) == "0") return(max(intensity)) # 全是阴性
      if(unique(levs) == "1") return(min(intensity)) # 全是阳性

      warning(str_glue("classification levels: {unique(levs)}, not 0 or 1 ..."), immediate. = TRUE) # 非0 1
      return(NULL)
    }

    if(length(levs) > 2 ){

      if(all(c("0", "1") %in% levs)){

        lev1 <- '0'
        lev2 <- "1"

      }else{

        warning(str_glue("{name} found more than 2 levels ..."), immediate. = TRUE)
        return(NULL)
      }
    }

    lev1 <- levs[1]
    lev2 <- levs[2]

    range1 <- range(intensity[class == lev1])
    range2 <- range(intensity[class == lev2])

    if(range1[1] < range2[2] || range1[2] > range2[1]){
      message("Intensity of two classification overlay.")
    }

    if(range1[2] > range2[2]){
      return(
        (range1[1] + range2[2]) / 2
      )
    }else{
      return(
        (range2[1] + range1[2]) / 2
      )
    }

  }

  res_cutoff <-
    all_markers %>%
    set_names(., .) %>%
    map(function(marker){

      names_labels <-
        list("Cell", "Membrane", "Cytoplasm", "Nucleus") %>%
        set_names(., .) %>%
        map(function(y){
          inten_idx <- str_glue("{marker} {y} {intensity_label}")
          if(y == "Cell"){
            class_idx <- str_glue("{marker} Positive {classification_label}")
          }else{
            class_idx <- str_glue("{marker} Positive {y} {classification_label}")
          }
          list(
            inten_idx = inten_idx,
            class_idx = class_idx
          )
        })

      if( all(! map_chr(names_labels, "class_idx") %in% all_colnames) ){
        # 若"Cell", "Membrane", "Cytoplasm", "Nucleus"均没有cutoff，则默认代表全是阴性表达
        # 设置为intensity的最大值
        return(
          names_labels %>%
            map_chr("inten_idx" ) %>%
            map(function(x){
              if(x %in% colnames(all_dat) ) return( max(all_dat[[x]]) )
              return(NULL)
            })
        )
      }

      names_labels %>%
        map(function(x){
          inten_idx <- x$inten_idx
          class_idx <- x$class_idx

          if(!class_idx %in% all_colnames) return(NULL)
          getCutoff(
            class = all_dat[[class_idx]],
            intensity = all_dat[[inten_idx]],
            name = marker
          )
        })

    })
  res_cutoff
}

if(F){
  all_data$dat[[1]] %>%
    parseCutoff() %>%
    map(discard, is.null) %>%
    discard(function(x)length(x) < 1)
}


# 5. parse cell location --------------------------------------------------
parseLocation <- function(x){
  # x: result of parseCutoff
  x %>%
    map(discard, is.null) %>% # 每个通道下的亚定位为NULL则去除
    discard(function(x)length(x) < 1) %>% # 若一个通道均没有定位则去除此通道
    map(function(x) names(x))
}

if(F){
  all_data$dat[[1]] %>%
    parseCutoff() %>%
    parseLocation()
}



# 6. parse roi label ------------------------------------------------------
parseROILable <- function(x, roi_label = "Analysis Region"){
  stopifnot(is.data.frame(x))

  if(!roi_label %in% colnames(x)){
    stop(str_glue("column {roi_label} not found"))
  }

  levs <- unique(x[[roi_label]])

  if(length(levs) > 1){
    stop(str_glue("found more than 1 ROI label: {str_c(levs, collapse = ', ')}."))
  }
  if(length(levs) < 1 || str_trim(levs) == ""){
    # 空字符
    stop("ROI label cannot be null or white character.")
  }
  return(levs)
}


# 7. setROIID -------------------------------------------------------------

setROIId <- function(raw_meta, tissue_id = FALSE, sep = " | "){
  stopifnot(is.data.frame(raw_meta))

  id_cols <- c('Slide_ID', 'Tissue_ID', 'ROI_ID')
  if(!tissue_id) {
    id_cols <- id_cols %>% setdiff("Tissue_ID")
  }

  if(any(!id_cols %in% colnames(raw_meta))){
    stop(str_glue('glue::glue_collapse(id_cols, sep = ",", last = ",") should all in raw_meta.'))
  }

  # 设置ID
  if(tissue_id){

    res <- raw_meta %>% mutate(ID = str_c(Slide_ID, Tissue_ID, ROI_ID, sep = sep))

  }else{

    res <- raw_meta %>% mutate(ID = str_c(Slide_ID, ROI_ID, sep = sep))

  }

  # if("Segment_ID" %in% colnames(raw_meta)){
  #
  #   res <- raw_meta %>% mutate(ID = str_c(Slide_ID, Tissue_ID, ROI_ID, Segment_ID, sep = sep))
  #
  # }else {
  #
  #   res <- raw_meta %>% mutate(ID = str_c(Slide_ID, Tissue_ID, ROI_ID, sep = sep))
  #
  # }

  if("Segment_ID" %in% colnames(raw_meta)){

    res <- res %>% mutate(ID = str_c(ID, Segment_ID, sep = sep))

  }

  return(res)
}



# 8. splitSegment ---------------------------------------------------------

splitSegment <-
  function(x, tissue_id = FALSE, sep = " | "){
    stopifnot(is(x, "CellDive"))

    if(!all(x@raw_meta$ID == names(x@raw_data))){
      stop("CellDive: raw_meta and raw_data should be the same order.")
    }

    len <- length(x@raw_data)
    if(len < 1) stop("Length less than 1, please check the length of raw_data.")

    seg_ls <- list()
    raw_data <- list()

    for(i in 1:len){
      # parse segment
      seg <-
        # 没有segment，则置segment为TotalSegment
        tryCatch(
          parseSegmentLable(x@raw_data[[i]]),
          error = function(x) {
            if(str_detect(str_to_lower(x$message), "column .* not found") ){
              return("TotalSegment")
            }else{
              stop(x$message)
            }
          }
        )

      # update segment
      seg_ls[[i]] <- seg

      if(length(seg) == 1){
        temp <- list(x@raw_data[[i]])
        names(temp) <- seg
        raw_data[[i]] <- temp
        next
      }else if(length(seg) == 0){
        stop("Please check script, because segment should only be 'TotalSegment' when segment not found in data..")
      }

      # update data
      seg_all <- parseSegmentLable(x@raw_data[[i]], return_all = TRUE)
      seg_all <- factor(seg_all, levels = seg)

      dat_ls <- x@raw_data[[i]] %>% split(f = seg_all)
      raw_data[[i]] <- dat_ls
    }
    names(raw_data) <- names(x@raw_data)

    # split segment
    raw_meta <- x@raw_meta
    raw_meta$Segment_ID <- seg_ls

    if(! 'original_ID' %in% colnames(raw_meta)){
      warning("Not found original_ID in raw_meta, set as 1:n()", immediate. = TRUE)

      raw_meta <- raw_meta %>% mutate(original_ID = 1:n())

    }

    raw_meta <- raw_meta %>%
      ungroup() %>%
      # mutate(original_ID = 1:n()) %>% # move into createCellDive.
      unnest(Segment_ID) %>%
      # mutate(original_ID = str_c(Slide_ID, Tissue_ID, ROI_ID, sep = sep)) %>%
      # mutate(ID = str_c(Slide_ID, Tissue_ID, ROI_ID, Segment_ID, sep = sep)) %>%
      setROIId(tissue_id = tissue_id)

    x@raw_meta <- raw_meta

    # split data
    x@raw_data <- raw_data %>%  flatDoubleList(sep = sep)

    # 无Segment的ID与有Segment的ID对应表
    id2id <- raw_meta$ID %>% split( raw_meta$original_ID %>% factor(., levels = unique(.)) )


    # 更新cutoff
    cutoff <- id2id %>% map2(x@cutoff, function(xx, yy) xx %>% map(~ yy)) %>% flatDoubleList()
    names(cutoff) <- id2id %>% unlist(use.names = FALSE)
    x@cutoff <- cutoff

    # 更新CellSubLocation
    cell_subposition <- id2id %>% map2(x@cell_subposition, function(xx, yy) xx %>% map(~ yy)) %>% flatDoubleList()
    names(cell_subposition) <- id2id %>% unlist(use.names = FALSE)
    x@cell_subposition <- cell_subposition

    # 更新channel
    channel <- id2id %>% map2(x@channel, function(xx, yy) xx %>% map(~ yy)) %>% flatDoubleList()
    names(channel) <- id2id %>% unlist(use.names = FALSE)
    x@channel <- channel

    # verify
    stopifnot( all( x@raw_meta$ID == names(x@raw_data) ) )
    stopifnot( all( x@raw_meta$ID == names(x@cutoff) ) )
    stopifnot( all( x@raw_meta$ID == names(x@cell_subposition) ) )
    stopifnot( all( x@raw_meta$ID == names(x@channel) ) )

    x
  }

if(F){
  cd2 <- splitSegment(cd)
}



# 9. chechCutoff ----------------------------------------------------------

checkCutoff <- function(x, marker = NULL){
  stopifnot(is(x, "CellDive"))

  if(length(marker) < 1){
    stop("Marker is null")
  }

  if(length(marker) < 1) {
    warning("Marker provided more than 1, only check the first marker", immediate. = TRUE)
    marker <- marker[1]
  }

  all_markers <- getChannel(x)
  if(!marker %in% all_markers){
    stop(str_glue("Marker {marker} not found."))
  }

  cutoff <- x@cutoff %>% map(marker) %>% map(discard, .p = is.null)
  cutoff_subloc <- cutoff %>% map(names)

  cell_subloc <- getCellSubPosition(x) %>% .[[marker]]

  is_cutoff_exist <- map2_lgl(cutoff_subloc, cell_subloc, ~ .y %in% .x)

  all(is_cutoff_exist)
}


# 10. resetCutoff ---------------------------------------------------------
# if not all markers have cutoff, find other sub_location to find cutoff, then resetIntensity
resetCutoff <- function(x, marker = NULL){
  stopifnot(is(x, "CellDive"))

  if(length(marker) < 1){
    stop("Marker is null")
  }

  is_cutoff_exist <- marker %>% map_lgl(~checkCutoff(x, marker = .x))
  sub_loc <- x %>% getCellSubPosition()

  if(sum(!is_cutoff_exist) < 1){
    warning("All marker have cutoff, no changed.")
    return(x)
  }

  warning("Some marker have no cutoff, try changing cell sub_location...")
  for(m in marker[!is_cutoff_exist]){
    cutoff <- x@cutoff %>%
      map(m) %>%
      map(discard, .p = is.null) %>%
      map(names) %>%
      map_chr(function(x){
        if(length(x) < 1) return("NULL")
        if("Cell" %in% x){
          if(length(x) == 1) return("Cell")
          return(x %>% setdiff("Cell") %>% .[1])
        }
      })
    if(any(cutoff == "NULL")){
      stop(str_glue("Channel {m}: not found cutoff in all ROI"))
    }else{
      warning(str_glue("Channel {m}: reset sublocation successed."))
      sub_loc[[m]] <- cutoff
    }
  }

  # resetIntensity
  xx <-
    x %>%
    setCellSubLocation(sub_loc) %>%
    resetIntensity()

  return(xx)
}


# 11  singlePositiveStat --------------------------------------------------------

# not recommended to use
# use positiveStat first
singlePositiveStat <- function(x, marker = NULL, forceResetIntensity = FALSE){
  stopifnot(is(x, "CellDive"))

  all_markers <- getChannel(x)
  if (is.null(marker) || !marker %in% all_markers) {
    stop("Marker is null or not in Celldive.")
  }

  # cheack cutoff 存在
  cutoff <- x@cutoff %>% map(marker) %>% map(discard, .p = is.null)
  cutoff_subloc <- cutoff %>% map(names)
  cell_subloc <- getCellSubPosition(x) %>% .[[marker]]

  is_cutoff_exist <-
    map2_lgl(cutoff_subloc, cell_subloc, ~ .y %in% .x)
  if (!all(is_cutoff_exist)) {
    if (!forceResetIntensity) {
      stop("Some cutoff of intensity in current cell_sublocation not exist.")
    }
    df <- getCellSubPosition(x)
    df[[marker]] <- cutoff_subloc %>%
      map(setdiff, y = "Cell") %>%
      map_chr(function(x)
        if (length(x) < 1)
          "Cell"
        else
          x[1])
    if (any(length(df[[marker]])) < 1) {
      stop(str_glue("Not all {marker} channel have cutoff value."))
    }
    x <- x %>%
      setCellSubLocation(df) %>%
      resetIntensity()
    cell_subloc <- getCellSubPosition(x) %>% .[[marker]]
  }

  # 提取cutoff
  cutoff_value <- map2(cutoff, cell_subloc, ~ .x[[.y]])

  # 提取表达量
  all_dat <- x %>% GetAssayData()

  all_meta <- x[[]] %>%
    mutate(.g = str_c(original_ID, Segment_ID))

  roi2cell <-
    all_meta  %>%
    rownames_to_column(".cell") %>%
    dplyr::select(.g, .cell) %>%
    group_by(.g) %>%
    nest() %>%
    deframe() %>%
    map(unlist, use.names = TRUE)

  stopifnot(length(cutoff_value) == length(roi2cell))

  ratio_dat <-
    roi2cell %>%
    map2(cutoff_value, function(x, y) {
      dat <- all_dat[marker, x, drop = FALSE] %>% unlist
      ratio <- sum(dat > y) / length(dat)
      ratio
    }) %>%
    enframe() %>%
    unnest(value)

  res <-
    all_meta %>%
    left_join(ratio_dat, by = c(.g = "name")) %>%
    distinct(.g, .keep_all = TRUE) %>%
    dplyr::select(ID, ratio = value)
  res
}

# 12. positiveStat -------------------------------------------------------

positiveStat <- function(x,
                         marker = NULL,
                         forceResetIntensity = FALSE,
                         result = c("ratio", "density", "cell_table"),
                         sep = "-"){

  # browser()

  # 返回结果的类型
  # ratio，返回细胞比例
  # density，返回细胞的密度，每个ROI的阳性细胞数据除以Area（未实施，因原始数据缺少ROI的面积）
  # cell_table，返回每个细胞的阴阳性结果（前两种是每个ROI中的结果）
  result <- match.arg(result)

  if(result == "density"){
    stop(
      "Not implemented: cannot obtain ROI/AOI area."
    )
  }


  stopifnot(is(x, "CellDive"))

  all_markers <- getChannel(x)
  if(is.null(marker) || !marker %in% all_markers){
    stop("Marker is null or not in Celldive.")
  }

  # check cutoff 存在
  is_cutoff_exist <- marker %>% map_lgl(~checkCutoff(x, marker = .x))

  if(!all(is_cutoff_exist)){

    if(!forceResetIntensity){
      stop("Some cutoff not found. Try `resetCutoff` first to find another cutoff.")
    }
    warning("Some cutoff not found.")
    x <- x %>% resetCutoff(marker = marker)
  }

  # 提取cutoff
  cutoff_vals <-
    marker %>%
    set_names(., .) %>%
    map(function(m){
      cutoff <- x@cutoff %>% map(m)
      sub_loc <- x %>% getCellSubPosition() %>% .[[m]]

      map2(cutoff, sub_loc, ~ .x[[.y]])
    })
  # cutoff_value <- map2(cutoff, cell_subloc, ~ .x[[.y]])


  # 提取表达量
  all_dat <- x %>% GetAssayData()

  all_meta <- x[[]] %>%
    mutate(.g = str_c(original_ID, Segment_ID))

  roi2cell <-
    all_meta  %>%
    rownames_to_column(".cell") %>%
    dplyr::select(.g, .cell) %>%
    group_by(.g) %>%
    nest()%>%
    deframe() %>%
    map(unlist, use.names = TRUE)

  roi2cell_name <- all_meta$ID[match(names(roi2cell), all_meta$.g)]

  stopifnot(length(cutoff_vals[[1]]) == length(roi2cell))

  pos_ls <-
    # ROI -> Marker -> pos_idx
    marker %>%
    set_names(., .) %>%
    map(function(m){
      cutoff_value <- cutoff_vals[[m]]
      stopifnot(all(roi2cell_name == names(cutoff_value)))

      map2(roi2cell, cutoff_value, function(x, y){
        val <- all_dat[m, x, drop=FALSE] %>% as.matrix() %>% c
        res <- val > y

        names(res) <- x # 注明每个细胞的ID

        res
      })
    }) %>%
    purrr::transpose()

  marker_combn <- marker %>% map(~c("Pos", "Neg")) %>% expand.grid()

  # 返回所有的marker的阴阳组合
  pos_combn_dat <-
    pos_ls %>%
    map(function(x){
      marker_combn %>%
        apply(1, simplify = FALSE, function(y){
          idx_ls <- map2(x, y, function(idx, pos_neg) if(pos_neg == "Neg") !idx else idx)
          names(idx_ls) <- NULL
          idx <- idx_ls %>% purrr::reduce(`&`)

          if(result == "cell_table") return(idx)

          idx <- sum(idx)/length(idx)
          idx
        }) %>%
        do.call(cbind, .) %>%
        as.data.frame()
    })

  marker_combn_name <- marker_combn %>% apply(1, str_c, collapse = sep)

  # 如果需要返回每个细胞的阴性性结果
  if(result == "cell_table"){
    combn_dat <- pos_combn_dat %>%
      map(~ {names(.x) <- marker_combn_name; .x}) %>%
      map(~ rownames_to_column(.x, var = 'ID') ) %>%
      do.call(rbind, .)

    return(combn_dat)
  }

  # 默认返回比例

  ratio_dat <- pos_combn_dat %>%
    map(c) %>%
    map(~ {names(.x) <- marker_combn_name; .x}) %>%
    enframe() %>%
    mutate(value = bind_rows(value)) %>%
    unnest(value)

  res <-
    all_meta %>%
    distinct(.g, .keep_all = TRUE) %>%
    left_join(ratio_dat, by = c(.g = "name")) %>%
    dplyr::select(ID, all_of(marker_combn_name))

  res
}


# 13. panelPostiveStat ----------------------------------------------------
# 此函数主要是可以传入一组细胞及细胞的marker列表，然后返回一个行为细胞列为ROI的data.frame
# 这个marker叫做panel

panelPostiveStat <- function(x,
                             panel,
                             forceResetIntensity = FALSE,
                             sep = "-",
                             result = c("ratio", "density", "cell_table")){

  # browser()

  stopifnot(is(x, "CellDive"))
  stopifnot(is(panel, "list"))

  # 返回结果的类型
  # ratio，返回细胞比例
  # density，返回细胞的密度，每个ROI的阳性细胞数据除以Area（未实施，因原始数据缺少ROI的面积）
  # cell_table，返回每个细胞的阴阳性结果（前两种是每个ROI中的结果）
  result <- match.arg(result)

  # all_markers
  all_markers <- x %>% getChannel()
  stopifnot(length(all_markers) > 0)

  # panel markers should be in CellDive
  all_markers_panel <- panel %>% map(names) %>% unlist() %>% unique()

  if(all(!all_markers_panel %in% all_markers)) stop("Please check marker labels in your panel, all label not found.")
  if(any(!all_markers_panel %in% all_markers)) {
    warning("Some panel marker not found will be ignored.", immediate. = TRUE)
    panel <- panel[panel %>% map(names) %>% map_lgl(~any(! .x %in% all_markers))]
  }

  # computate cell proportion
  res <-
    panel %>%
    imap(function(marker, cell){
      channel <- names(marker)

      pos_neg <- marker %>%
        str_trim() %>%
        fct_recode(Pos = "+", Neg = "-") %>%
        str_c(collapse = "-")

      props_stat <- positiveStat(x,
                                 marker = channel,
                                 forceResetIntensity = forceResetIntensity,
                                 sep = sep,
                                 result = result)

      stat <- props_stat[, c("ID", pos_neg), drop=FALSE]
      colnames(stat)[2] <- cell
      stat
    })

  # 是否所有的细胞的ROI数量和顺序都是一样的
  all_same <- res %>% map('ID') %>% map_lgl(~ all(.x == res[[1]][["ID"]])) %>% all()
  if(all_same) {
    res_df <- res %>%
      map_at(-1, dplyr::select, -ID) %>%
      bind_cols() %>%
      remove_rownames() %>%
      column_to_rownames("ID") %>% t %>%
      as.data.frame()

    return(res_df)
  }

  res_df <- res %>%
    purrr::reduce(left_join, by = "ID") %>%
    column_to_rownames("ID") %>% t %>%
    as.data.frame()
  return(res_df)
}


# 14. panelPostiveStat_speed----------------------------------------------------------------
# 此函数等同于panelPostiveStat
# 但是速度可提升17倍
# 25s -- 1.4s

panelPostiveStat_speed <- function(x,
                                   panel,
                                   forceResetIntensity = FALSE,
                                   sep = "-",
                                   result = c("ratio", "density", "cell_table")) {

  # browser()

  stopifnot(is(x, "CellDive"))
  stopifnot(is(panel, "list"))

  # 返回结果的类型
  # ratio，返回细胞比例
  # density，返回细胞的密度，每个ROI的阳性细胞数据除以Area（未实施，因原始数据缺少ROI的面积）
  # cell_table，返回每个细胞的阴阳性结果（前两种是每个ROI中的结果）
  result <- match.arg(result)

  # all_markers
  all_markers <- x %>% getChannel()
  stopifnot(length(all_markers) > 0)

  # panel markers should be in CellDive
  all_markers_panel <- panel %>% map(names) %>% unlist() %>% unique()

  if(all(!all_markers_panel %in% all_markers)) stop("Please check marker labels in your panel, all label not found.")
  if(any(!all_markers_panel %in% all_markers)) {
    warning("Some panel marker not found will be ignored.", immediate. = TRUE)
    panel <- panel[panel %>% map(names) %>% map_lgl(~any(! .x %in% all_markers))]
  }


  # check cutoff existed
  is_cutoff_exist <- all_markers_panel %>% map_lgl(~checkCutoff(x, marker = .x))

  if(!all(is_cutoff_exist)){

    if(!forceResetIntensity){
      stop("Some cutoff not found. Try `resetCutoff` first to find another cutoff.")
    }
    warning("Some cutoff not found.")
    x <- x %>% resetCutoff(marker = marker)
  }

  # compute cutoff
  cutoff_vals <-
    all_markers_panel %>%
    set_names(., .) %>%
    map(function(m){
      cutoff <- x@cutoff %>% map(m)
      sub_loc <- x %>% getCellSubPosition() %>% .[[m]]

      map2(cutoff, sub_loc, ~ .x[[.y]])
    })

  cutoff_df <- cutoff_vals %>%
    lapply(unlist) %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("ID")

  # 汇总cutoff表，各个channel的cutoff、obj_id、ID
  roi_cell <- x[["ID"]] %>% rownames_to_column("obj_id")

  merge_cutoff_df <- merge(roi_cell, cutoff_df, by = "ID", sort = FALSE)
  rownames(merge_cutoff_df) <- merge_cutoff_df$obj_id

  # 保证cutoff和对象的细胞顺序一致
  stopifnot(nrow(merge_cutoff_df) == ncol(x))
  merge_cutoff_df <- merge_cutoff_df[colnames(x), ]

  # all panel positive data in all cells
  all_dat <- x %>% GetAssayData() %>% as.matrix() %>% t

  analysis_dat <- all_dat[, all_markers_panel, drop=FALSE]
  analysis_cutoff <- as.matrix( merge_cutoff_df[, all_markers_panel, drop=FALSE] )

  # res即为每个marker的阴阳性数据
  res <- analysis_dat > analysis_cutoff
  rownames(res) <- rownames(analysis_dat)

  # all_cell_table
  celltype_res_ls <-
    map(panel, function(marker){# browser()
      m <- names(marker)
      neg_idx <- which(marker == "-")

      analysis_res <- res[, m, drop = FALSE]
      if(length(neg_idx) > 0){
        analysis_res[, neg_idx] <- !analysis_res[, neg_idx, drop=FALSE]
      }

      cell_table <- rowSums(analysis_res) == length(m)

      names(cell_table) <- rownames(res)

      # cell_table or ratio
      if(result == "cell_table"){
        return(cell_table)
      }

      # sort by merge_cutoff_df
      cell_table  <- cell_table[rownames(merge_cutoff_df)]
      pos_count   <- tapply(cell_table, merge_cutoff_df$ID, FUN = sum)
      total_count <- tapply(cell_table, merge_cutoff_df$ID, FUN = length)

      roi_ratio <- pos_count/total_count

      roi_ratio
    })

  # browser()

  # cell_id/Segment为行, cell_type为列
  celltype_res <- as.data.frame(celltype_res_ls, check.names = FALSE)

  # return
  res_final <- as.data.frame(t(celltype_res), check.names = FALSE)
  res_final
}

# yyy. class definition ---------------------------------------------------

setClass(
  "CellDive",
  contains = "Seurat",
  slot = c(
    raw_meta = "data.frame",
    raw_data = "list",
    cutoff = "list",
    cell_subposition = "list",
    channel = "list"
  )
)

setMethod("show", "CellDive", function(object){
  class(object) <- "Seurat"
  cat("An object of class CellDive\n")
  cat(str_glue("{nrow(object@raw_meta)} ROIs total {n_distinct(unlist(object@channel))} channels\n\n"))
  cat("------------------------\n")

  # no seurat data
  if(length(object@assays) < 0){
    cat("An object of class Seurat\n")
    cat("None of assays found\n")
  }else{
    print(object)
  }

})

# [

`[.CellDive` <- function(x, i, j, ...){
  x <- as(x, "Seurat")

  if (missing(i) && missing(j)) {
    return(x)
  }
  if (missing(i)) {
    i <- NULL
  }
  else if (missing(j)) {
    j <- colnames(x)
  }

  x[i, j, ...]
}

# GetAssayData
GetAssayData.CellDive <- function(x, ...){
  x <- as(x, "Seurat")

  GetAssayData(x, ...)
}

# getCellSubPosition
getCellSubPosition <- function(x){
  stopifnot(is(x, "CellDive"))

  if(length(x@cell_subposition) < 1){
    stop("No cell_subposition info found.")
  }

  all_channels <- x@cell_subposition %>% map(names) %>% unlist() %>% unique()
  cell_subposition_df <-
    x@cell_subposition %>%
    map(function(x){
      others <- all_channels %>% setdiff(names(x))
      x[others] <- NA
      x
    }) %>%
    map(as.data.frame, check.names = FALSE) %>%
    bind_rows() %>%
    mutate(ROI = names(x@cell_subposition)) %>%
    dplyr::select(any_of("ROI"), everything())

  cell_subposition_df

}

if(F){
  getCellSubPosition(cd)

}

# setCellSubLocation
setCellSubLocation <- function(x, df){
  stopifnot(is(x, "CellDive"))
  stopifnot(is.data.frame(df) || is.matrix(df))

  df <- as.data.frame(df)

  roi_len <- nrow(x@raw_meta)
  if(roi_len != nrow(df)){
    stop("Wrong CellSubLocation number.")
  }

  if("ROI" != colnames(df)[1]){
    stop("The first colname of CellSubLocation should be 'ROI'")
  }

  cell_subloc <- df %>% apply(1, function(x)x[-1], simplify = FALSE)
  names(cell_subloc) <- df[["ROI"]]

  cell_subloc2 <- cell_subloc %>% map(~ .[!is.na(.)]) %>% map(as.list)

  x@cell_subposition <- cell_subloc2

  x
}

if(F){
  cell_sublocation <- getCellSubPosition(cd)
  cell_sublocation$CD20 <- "Cell"
  setCellSubLocation(cd, cell_sublocation) %>% getCellSubPosition()
}

# setCellSubLocationByList
# 手动指定一个sub_location表格，然后根据是否相应的channel有cutoff值，来决定是否采用此channel
setCellSubLocationByList <- function(x, sub_location, force = TRUE){ # browser()
  stopifnot(is(x, "CellDive"))
  stopifnot(is.list(sub_location))

  # provided channel
  provide_channel <- names(sub_location)
  if(is.null(provide_channel)) stop("sub_location should have names.")
  stopifnot(length(sub_location) > 0)

  all_channels <- x %>% getChannel()

  if(any(!provide_channel %in% all_channels)){
    stop("Check the provided sub_location, some channel not found.")
  }

  # sub_locatin_in_celldive
  sub_loc_in_cd <- x %>% getCellSubPosition()
  all_cutoff <- x %>% showCutoff()

  # 更新sub_loc_in_cd
  for(i in seq_along(sub_location)){
    sub_loc <- sub_location[[i]]
    chann <- names(sub_location)[i]

    # update sub_loc_in_cd，保证有cutoff，如果没有则置为cell，
    # 另外如果Membrane没有cutoff，则尝试Cytoplasm，若Cytoplasm也没有cutoff，则置为Cell

    sub_list <- c("Cell", "Membrane", "Cytoplasm", "Nucleus")

    # 每一个channel下面的需要尝试进行设置的sub_loc可能不止一个
    for(sub in sub_loc){
      if(! sub %in% sub_list) stop('sub_location should only be one of "Cell", "Membrane", "Cytoplasm", "Nucleus"')

      if('NULL' %in% all_cutoff[[chann]][[sub]]){# 不是所有的ROI都有cutoff

        if(sub == "Membrane"){
          # 改为Cytoplasm重新尝试
          if(! 'NULL' %in% all_cutoff[[chann]][["Cytoplasm"]]){
            sub_loc_in_cd[[chann]] <- "Cytoplasm"
            break # 只要设置成功，就不会再判断后续的备选sub_loc
          }
        }

        if(sub == sub_loc[length(sub_loc)]){
          # 如果是最后一个channel，依然没有，则置为Cell
          sub_loc_in_cd[[chann]] <- "Cell"
          break
        }

      }else{
        sub_loc_in_cd[[chann]] <- sub
        break
      }
    }
  }

  # 已更新，则设置sub_loc
  x <- setCellSubLocation(x, sub_loc_in_cd)

  # force resetIntensity
  if(force) x <- x %>% resetIntensity()

  return(x)
}

# getChannel
getChannel <- function(x){
  stopifnot(is(x, "CellDive"))

  if(length(x@channel) < 1){
    stop("No channel info found.")
  }

  # all channels
  all_channels <- x@channel %>% unlist() %>% unique()

  return(all_channels)
}

if(F){
  getChannel(cd)
}

# showCutoff
# 不同于parseCutoff，获得的结果是一个多级的列表
# showCutoff用于查看每个channel中，行为ROI列为sub_location的cutoff表

showCutoff <- function(x, marker = NULL, return_all = FALSE){
  stopifnot(is(x, "CellDive"))

  if(is.null(x@cutoff)){

    if(is.null(raw_data)) stop("No raw_data slot in CellDive.")

    x@cutoff <- x@raw_data %>% map(parseCutoff)

  }

  cutoff <- x@cutoff

  # all channel
  all_marker <- x %>% getChannel()

  if(is.null(marker)){
    warning(str_glue("marker not set , return all channel."), immediate. = TRUE)

  }else{

    if(!marker %in% all_marker) {
      warning("Wrong marker, return all channel.", immediate. = TRUE)
    }
    else{
      all_marker <- marker
    }

  }


  if(return_all){
    all_marker <- x %>% getChannel()
  }



  cutoff_by_marker <-
    all_marker %>%
    set_names(., .) %>%
    map(function(x){ # x -> channel
      cutoff %>%
        imap(function(y, name){ # y -> roi, name -> roi_name

          if(x %in% names(y)) return(y[[x]])

          # 没有此marker，则返回NULL
          warning(str_glue("ROI {name} not found channel {x}"))
          list(
            Cell = NULL,
            Membrane = NULL,
            Cytoplasm = NULL,
            Nucleus = NULL
          )
        })
    })

  res <-
    cutoff_by_marker %>%
    map(function(x){# x -> marker
      x_modified <-
        # 将NULL修改为字符'NULL'
        x %>%
        map(function(xx){
          xx %>%
            enframe %>%
            mutate(
              value = map_chr(
                value,
                ~ ifelse(is.null(.x), 'NULL', as.character(round(.x, 2)))
              )
            )
        })

      # NULL改为字符'NULL'后，就可以修改data.frame的样式了
      x_final <-
        x_modified %>%
        enframe(name = "ID") %>%
        unnest("value") %>%
        pivot_wider(
          names_from = "name",
          values_from = "value"
        )
      x_final
    })

  # return
  # if(return_all) return(res)

  # if(is.null(marker)){
  #
  #   warning(str_glue("marker not set , return all channel."), immediate. = TRUE)
  #   return(res)
  # }
  #
  # if(!marker %in% all_marker){
  #
  #   warning(str_glue("{marker} not found, return all channel."), immediate. = TRUE)
  #   return(res)
  # }

  return(res)

}

# fillCutoff
# 通过一些数理方法用于填补Cutoff的空白
# 注意，仅用于应急，可信度不明确
fillCutoff <- function(x, method = c("sd", "iqr", "others")){
  stopifnot(is(x, "CellDive"))

  method <- match.arg(method)

  # 获取cutoff
  if(is.null(x@cutoff)){

    if(is.null(raw_data)) stop("No raw_data slot in CellDive.")

    x@cutoff <- x@raw_data %>% map(parseCutoff)

  }
  cutoff <- x@cutoff

  # 获取raw_data
  if(is.null(x@raw_data)){
    stop("No raw_data slot in CellDive.")
  }
  raw_data <- x@raw_data

  # 三层循环
  # cutoff: ROI -> channel -> sub_loc
  for( roi_idx in seq_along(cutoff)){
    # roi cutoff
    roi_name <- names(cutoff)[roi_idx]
    roi_cutoff  <- cutoff[[roi_idx]]

    # roi intensity
    roi_intensity <- raw_data[[roi_idx]]

    # roi_cutoff: channel -> sub_loc
    for(marker_idx in seq_along(roi_cutoff)){
      marker_name <- names(roi_cutoff)[marker_idx]
      marker_cutoff <- roi_cutoff[[marker_idx]]

      # sub_loc: 一共四种
      for( sub_loc_idx in seq_along(marker_cutoff)){
        sub_loc <- names(marker_cutoff)[sub_loc_idx]
        cutoff_val <- marker_cutoff[[sub_loc_idx]]

        # 有cutoff，不做处理
        if(!is.null(cutoff_val)) next()

        # 没有Intensity，不做处理
        intensity_label <- createIntensityLabel(marker_name, sub_loc)
        if(! intensity_label %in% colnames(roi_intensity)) next()

        # 开始计算cutoff
        intensity <- roi_intensity[[intensity_label]]

        if(method == "quantile"){
        # mean + 2*sd

          cutoff_computed <- mean(intensity) + 2 * sd(intensity)

        }else if(method == "iqr"){
        # median + 1.5 * iqr

          cutoff_computed <- median(intensity) + 1.5 * IQR(intensity)

        }else{
          # others
          cutoff_computed <- quantile(intensity, 0.95)

        }
        cutoff[[roi_name]][[marker_name]][[sub_loc]] <- cutoff_computed
      }
    }
  }

  stopifnot(length(cutoff) == length(x@cutoff))
  stopifnot(length(flatDoubleList(cutoff)) == length(flatDoubleList(x@cutoff)))

  # return
  x@cutoff <- cutoff
  x
}

# resetIntensity
resetIntensity <- function(x,
                           sep = " | ",
                           merge_not_identical_channels = FALSE,
                           project = "celldive",
                           assay = "protein",
                           names.field = 1,
                           names.delim = "_",
                           meta.data = NULL,
                           ...) {

  stopifnot(is(x, "CellDive"))

  # channels
  stopifnot(!is.null(x@channel))
  stopifnot(!is.null(x@cell_subposition))
  stopifnot(!is.null(x@channel))

  celldive <- x

  # 根据sub_location, 设置counts
  all_channels <- unique(unlist(celldive@channel, use.names = FALSE) )
  omit_channels <- celldive@channel %>% map(~ setdiff(all_channels, .x)) %>% keep(function(x) length(x) > 0)
  if(length(omit_channels) > 0){
    if(!merge_not_identical_channels){
      stop("Not all ROIs have identity channels. If you still want to merge data please set merge_not_identical_channels as TRUE.")
    }
    warning("Not all ROIs have identity channels.", immediate. = TRUE)
  }

  # raw_counts
  raw_counts <- list()
  for(i in seq_along(celldive@raw_data)){
    .roi <- names(celldive@raw_data)[i]
    .dat <- celldive@raw_data[[i]]

    .sub_loc  <- celldive@cell_subposition[[i]]
    .channels <- celldive@channel[[i]]

    .sub_loc <- .sub_loc[.channels] # make sure identical channels order

    if(any(names(.sub_loc) != .channels)){
      stop("Channels recording in sub_location not equals to channels")
    }

    inten_name <- createIntensityLabel(.channels, .sub_loc)
    # print(inten_name)

    inten_dat <- as.data.frame(t(.dat[, inten_name, drop=FALSE]))
    colnames(inten_dat) <- str_c(.roi, .dat[["Object Id"]], sep = " | ")
    rownames(inten_dat) <- .channels

    absent_channel <- omit_channels[[.roi]]
    if(length(absent_channel) > 0){
      raw_counts[, absent_channel] <- NA
    }

    raw_counts[[i]] <- inten_dat
  }
  names(raw_counts) <- names(celldive@raw_data)


  # counts
  counts <-
    raw_counts %>%
    map(~ .x[all_channels, ,drop=FALSE]) %>%
    bind_cols()
  colnames(counts) <- make.names(colnames(counts))

  stopifnot(all(celldive@raw_meta$ID == names(celldive@raw_data)))

  # seurat meta_data
  # sc_meta <-
  #   celldive@raw_meta %>%
  #   mutate(Object_ID = celldive@raw_data %>% map("Object Id")) %>%
  #   unnest(Object_ID) %>%
  #   mutate(obj_id = str_c(ID, Object_ID, sep = sep)) %>%
  #   discard(is.list) %>%
  #   as.data.frame()
  # rownames(sc_meta) <- colnames(counts)
  sc_meta <-
    celldive@raw_meta %>%
    mutate(Object_ID = celldive@raw_data %>% map("Object Id")) %>%
    mutate(
      XMin = celldive@raw_data %>% map("XMin"),
      XMax = celldive@raw_data %>% map("XMax"),
      YMin = celldive@raw_data %>% map("YMin"),
      YMax = celldive@raw_data %>% map("YMax"),
    ) %>%
    unnest(c(Object_ID, XMin, XMax, YMin, YMax)) %>%
    mutate(obj_id = str_c(ID, Object_ID, sep = sep)) %>%
    discard(is.list) %>%
    as.data.frame()
  rownames(sc_meta) <- colnames(counts)

  # Seurat class
  sc <- CreateSeuratObject(
    counts      =  counts,
    project     =  project,
    assay       =  assay,
    names.field =  names.field,
    names.delim =  names.delim,
    meta.data   =  sc_meta,
    ...
  )

  # sc <- sc %>% AddMetaData(sc_meta)

  # copy Seurat data to celldive
  for(sl in slotNames(sc)){
    slot(celldive, sl) <- slot(sc, sl)
  }

  # warnning when duplicated ID foundn
  dup_idx <- duplicated(celldive@raw_meta$ID)
  if(sum(dup_idx) > 0){
    tab <- table(names(celldive@raw_meta$ID))
    warning("Some ROI ID(s) are duplicated: ", immediate. = TRUE)
    print(tab[celldive@raw_meta$ID[dup_idx]])
  }


  return(celldive)
}

# createCellDive
createCellDive <- function(counts = NULL,

                           raw_meta = NULL, # data.frame
                           raw_data = NULL, # list
                           cutoff = NULL, # list
                           cell_subposition = NULL, # list
                           channel = NULL, #

                           tissue_id = FALSE,
                           dat_col = 'dat',
                           sep = " | ",

                           verbose = TRUE,

                           force_dummy = FALSE,

                           parse_roi_label = TRUE,
                           merge_not_identical_channels = FALSE,

                           project = "celldive",
                           assay = "protein",
                           names.field = 1,
                           names.delim = "_",
                           meta.data = NULL,
                           ...){

  # new celldive
  celldive <- new(
    "CellDive",
    version = packageVersion("SeuratObject")
  )

  ### set ###
  # raw_meta
  if(verbose){
    warning("### Setting raw_meta ###", immediate. = TRUE)
  }
  if(!is.null(raw_meta)){

    if(!"ID" %in% colnames(raw_meta)){# 需要设置好ID
      if(is.null(raw_data)){

        # dat_col默认是'dat'，是保留列，如果存在于raw_meta中，则不再继续读取
        if(dat_col %in% colnames(raw_meta) && is.list(raw_meta[[dat_col]])){

          raw_data <- raw_meta[[dat_col]]

        }else{
          suppressMessages(
            suppressWarnings(
              raw_data <- raw_meta$path %>% map(read_csv)
            )
          )
        }

      }
      if(parse_roi_label){
        roi_id <- raw_data %>% map_chr(parseROILable)
        if(length(roi_id) == nrow(raw_meta)){
          raw_meta$ROI_label <- raw_meta$ROI_ID
          raw_meta$ROI_ID <- roi_id
        }else{
          warning("Length of parsed ROI_ID not equal to raw_meta")
        }
      }

      raw_meta <- raw_meta %>%
        setROIId(tissue_id = tissue_id, sep = sep) %>%
        mutate(original_ID = 1:n())

      if("ID" %in% colnames(raw_meta)) names(raw_data) <- raw_meta$ID

      if(any(duplicated(raw_meta$ID))) {
        warning("Some ID of ROI duplicated. Did you forget parse ROI label?", immediate. = TRUE)
      }
    }
    celldive@raw_meta <- raw_meta
  }

  # raw_data
  if(verbose){
    warning("### Setting raw_data ###", immediate. = TRUE)
  }
  if(!is.null(raw_data)){
    celldive@raw_data <- raw_data
  }
  # else if(!is.null(raw_meta)){
  #   celldive@raw_data <- raw_data
  # }

  # cutoff
  if(verbose){
    warning("### Setting cutoff ###", immediate. = TRUE)
  }
  if(!is.null(cutoff)){
    celldive@cutoff <- cutoff
  }else if(length(celldive@raw_data) > 0){
    cutoff <- celldive@raw_data %>% map(parseCutoff)
    celldive@cutoff <- cutoff
  }

  # cell_subposition
  if(verbose){
    warning("### Setting cell_subposition ###", immediate. = TRUE)
  }
  if(!is.null(cell_subposition)){
    celldive@cell_subposition <- cell_subposition
  }else if(!is.null(cutoff)){

    # set cell sub_loc
    cell_subposition <- cutoff %>% map(~map(.x, discard, is.null))
    cell_subposition <- cell_subposition %>%
      map(~map(.x, function(x){
        if(length(x) < 1) return("Cell")
        if(length(x) == 1) return(names(x))
        xx <- x %>% names() %>% setdiff("Cell") %>% sort %>%  .[1]
        return(xx)
      }))

    celldive@cell_subposition <- cell_subposition
  }

  # channel
  if(verbose){
    warning("### Setting channel ###", immediate. = TRUE)
  }
  if(!is.null(channel)){
    celldive@channel <- channel
  }else if(!is.null(cutoff)){
    channel <- cutoff %>% map(names)
    celldive@channel <- channel
  }


  # force_dummy: 返回假数据，不将CellDive数据更新到Seurat的slot中
  if(!force_dummy){
    # 数据充足，则直接resetIntensity
    # 否则仍然使用假数据
    if(!is.null(celldive@raw_data) &&
       !is.null(celldive@cell_subposition) &&
       !is.null(celldive@channel)){

        if(verbose){

        }
        return(resetIntensity(
          celldive,
          sep = sep,
          merge_not_identical_channels = merge_not_identical_channels,
          project     =  project,
          assay       =  assay,
          names.field =  names.field,
          names.delim =  names.delim,
          meta.data   =  meta.data,
          ...
        ))
    }
  }

  # 如果counts没有指定，则使用假数据
  if(is.null(counts)){
    counts <-
      matrix(
        1:12,
        nrow = 3,
        dimnames = list(
          paste0("gene", 1:3),
          paste0("cell", 1:4)
        )
      )
  }

  # Seurat class
  sc <- CreateSeuratObject(
    counts      =  counts,
    project     =  project,
    assay       =  assay,
    names.field =  names.field,
    names.delim =  names.delim,
    meta.data   =  meta.data,
    ...
  )

  sc <- sc %>% AddMetaData(sc_meta)

  # copy Seurat data to celldive
  for(sl in slotNames(sc)){
    slot(celldive, sl) <- slot(sc, sl)
  }

  return(celldive)

}

# SegmentLabel --> 1 -> 多
# raw_meta <--- all_data
# raw_data <--- all_data
if(F){
  # test
  undebug(createCellDive)
  cd <-
    createCellDive(
    raw_meta = all_data[1:10, ]
  )


  cd2 <-
    createCellDive(
      raw_meta = all_data[1:3, ]
    )

  # test
  createCellDive(
    iris[, -5] %>% t %>%
      as.data.frame() %>%
      set_names(., str_c("cell", 1:ncol(.)))
  ) %>%
    NormalizeData() %>%
    FindVariableFeatures()%>%
    ScaleData() %>%
    RunPCA(npcs = 3)  %>%
    AddMetaData(., iris$Species, col.name = "group") %>%
    DimPlot(group.by = "group")


  celldive <- new(
    "CellDive",
    version = packageVersion("SeuratObject")
  )
}




# zzz. tools ----------------------------------------------------------------

list.roi <- function(x, roi_label = "Phenotype"){

  if(length(x) < 1) return(character(0))
  if(length(x) > 1){
    stop(
      str_glue("please only provide 1 path to list roi:\n{str_c(x, collapse = ', ')}")
    )
  }

  all.files <- list.files(x, full = TRUE)

  roi.files <- all.files %>% str_subset(pattern = roi_label)

  return(roi.files)

}

createIntensityLabel <- function(channel, sub_loc){
  return(
    str_c(channel, " ", sub_loc, " Intensity")
  )
}

#  define func.
flatDoubleList <- function(x, sep = "---", ...){

  res <- list()

  name1 <- names(x)
  if(is.null(name1)) name1 <- seq_along(x)
  name2 <- x %>%
    map(function(x){
      name <- names(x)
      if(is.null(name)) name <- seq_along(x)
      name
    })
  name_ls <- list()
  id <- 1
  for(i in seq_along(x)){
    for(j in seq_along(x[[i]])){
      val <- x[[i]][[j]]
      name <- paste(name1[i], name2[[i]][j], sep = sep)

      name_ls[[id]] <- name
      res[[id]] <- val

      id <- id + 1
    }
  }
  names(res) <- name_ls
  res
}

# define nestDoubleList
nestDoubleList <- function(x, sep = "---", simple = TRUE, rev = FALSE, ...){
  res <- list()

  name <- names(x)
  if(any(!str_detect(name, sep))) stop(str_glue("All names should have separator '{sep}'"))

  name_pattern <- name %>% str_split(pattern = sep)

  if(rev){
    if(simple){
      warning("Simple will be 'FALSE' when rev is 'TRUE'")
      simple <- FALSE
    }
  }

  for(i in name_pattern){

    suppressWarnings(not_num <- is.na(parse_integer(i[2])))

    if(simple && !not_num){
      name_temp <- str_c(i, collapse = sep)
      res[[i[1]]] <- x[[name_temp]]
    }else{
      name_temp <- str_c(i, collapse = sep)
      if(rev){
        res[[i[2]]][[i[1]]] <- x[[name_temp]]
      }else{
        res[[i[1]]][[i[2]]] <- x[[name_temp]]
      }

    }

  }
  return(res)
}

## 验证anno的列是否填写内容
isNA <- function(df){

  if(is.data.frame(df) || is.matrix(df)){
    if(ncol(df) > 1) stop("Data Frame should only have one column!")
    dat <- unlist(df)
  }else{
    dat <- df
  }

  ratio <- length(which(is.na(dat)))/length(dat)
  ifelse(ratio > 0.9, TRUE, FALSE)
}


## 获取客户注释的列名
getUserColumnName <- function(anno_data, singleLevel = F){
  cols_name <- colnames(anno_data)
  cols_name <- setdiff(cols_name, c("No.",  "Scan_ID",  "ROI_ID",  "Segment_ID")) # 固定列，去除

  cols_name <- map_chr(cols_name,  # 不为NA的列返回名字，否则标记为"NA"
                       ~ifelse(!isNA(anno_data[, .]), ., "NA"))
  cols_name <- cols_name[!cols_name=="NA"]
  # 去除只有一个水平的列
  if(singleLevel){
    idx <- map_lgl(cols_name,
                   ~ifelse(length(unique(anno_data[, ., drop=T])) == 1, TRUE, FALSE))
    cols_name <- cols_name[!idx]
  }
  cols_name
}


# define a function used to add 'Others T'
add_others_t <- function(df, result = c("ratio", "density", "cell_table"), negative2zero = 0){
  # 算法：Others T = T - Treg - Tc - Th

  result <- match.arg(result)

  if("Others T" %in% rownames(df)) {
    warning("Others T existed.")
    return(df)
  }

  all_cells <- c("T cell", "T regulatory", "T cytotoxic", "T helper")

  if(any(!all_cells %in% rownames(df))){
    message_ <- glue::glue_collapse(all_cells, sep = ", ", last = " and ")
    stop(str_glue("One of more of {message_} not found in df, Others T cannot be compute ."))
  }


  t_total <- df[all_cells[1], ] %>% unlist(use.names = FALSE)
  t_reg   <- df[all_cells[2], ] %>% unlist(use.names = FALSE)
  t_c     <- df[all_cells[3], ] %>% unlist(use.names = FALSE)
  t_h     <- df[all_cells[4], ] %>% unlist(use.names = FALSE)

  if(result == "ratio"){

    others_t <-  t_total - t_reg - t_c - t_h

    # 'Others T' can be less than 0, because HALO cutoff may not be 100% accurate
    if(negative2zero){

      less_than_negative0.1 <- others_t < -0.1
      if(any(less_than_negative0.1)){
        stop("Others T less than -0.1, please check.")
      }

      less_than_0 <- others_t < 0
      if(any(less_than_0)){
        others_t[less_than_0] <- 0
      }
    }
    df["Others T", ] <- others_t

  }else if(result == "cell_table"){

    others_t <-  t_total & (!t_reg) & (!t_c) & (!t_h)
    df["Others T", ] <- others_t

  }else{

    stop(str_glue("result cannot be {result}"))

  }

  return(df)
}


# define a function used to add 'Others T'
add_others_t2 <- function(df, result = c("ratio", "density", "cell_table"), negative2zero = 0){
  # 算法：Others T = T - Treg - Tc - Th

  stopifnot(tibble::has_rownames(df))

  result <- match.arg(result)

  if("Others T" %in% rownames(df)) {
    warning("Others T existed.")
    return(df)
  }

  all_cells <- c(
    "T cell",
    "T regulatory",
    "T cytotoxic",
    "T helper",
    "T regulatory memory",
    "T cytotoxic memory",
    "T helper memory"
  )

  # all_cells > 2
  if(any(!all_cells %in% rownames(df))){
    message_ <- glue::glue_collapse(all_cells, sep = ", ", last = " and ")
    stop(str_glue("One of more of {message_} not found in df, Others T cannot be compute ."))
  }

  # COMPUTE
  t_total <- df[all_cells[1], ] %>% unlist(use.names = FALSE)
  t_named <- df[all_cells[-1], ,drop=FALSE]

  if(result == "ratio"){

    others_t <-  t_total - apply(t_named, 2, sum)

    # 'Others T' can be less than 0, because HALO cutoff may not be 100% accurate
    if(negative2zero){

      less_than_negative0.1 <- others_t < -0.1
      if(any(less_than_negative0.1)){
        stop("Others T less than -0.1, please check.")
      }

      less_than_0 <- others_t < 0
      if(any(less_than_0)){
        others_t[less_than_0] <- 0
      }
    }
    df["Others T", ] <- others_t

  }else if(result == "cell_table"){

    others_t <-  t_total & (!apply(t_named, 2, any))
    df["Others T", ] <- others_t

  }else{

    stop(str_glue("result cannot be {result}"))

  }

  return(df)
}
