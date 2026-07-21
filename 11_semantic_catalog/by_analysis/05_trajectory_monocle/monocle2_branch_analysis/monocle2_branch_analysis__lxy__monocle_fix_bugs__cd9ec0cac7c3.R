require(monocle)
# monocle fix bugs
BEAM <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
         reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", branch_states = NULL, 
         branch_point = 1, relative_expr = TRUE, branch_labels = NULL, 
         verbose = FALSE, cores = 1, ...){
  branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr, 
                               reducedModelFormulaStr = reducedModelFormulaStr, branch_states = branch_states, 
                               branch_point = branch_point, relative_expr = relative_expr, 
                               cores = cores, branch_labels = branch_labels, verbose = verbose, 
                               ...)
  cmbn_df <- branchTest_res[, 1:4]
  if (verbose) 
    message("pass branchTest")
  fd <- fData(cds)[row.names(cmbn_df), ]
  cmbn_df <- cbind(cmbn_df, fd)
  if (verbose) 
    message("return results")
  return(cmbn_df)
}


branchTest <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
    reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", branch_states = NULL, 
    branch_point = 1, relative_expr = TRUE, cores = 1, branch_labels = NULL, 
    verbose = FALSE, ...) {
  
    if ("Branch" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
        cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = branch_states, 
            branch_point = branch_point, branch_labels = branch_labels, 
            ...)
    }
    else cds_subset <- cds
    branchTest_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = fullModelFormulaStr, 
        reducedModelFormulaStr = reducedModelFormulaStr, cores = cores, 
        relative_expr = relative_expr, verbose = verbose)
    return(branchTest_res)
}

buildBranchCellDataSet <- 
  function (cds, progenitor_method = c("sequential_split", "duplicate"), 
    branch_states = NULL, branch_point = 1, branch_labels = NULL, 
    stretch = TRUE) {
    
    progenitor_method <- match.arg(progenitor_method)
    
    if (is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime)) 
        stop("Please first order the cells in pseudotime using orderCells()")
    if (is.null(branch_point) & is.null(branch_states)) 
        stop("Please either specify the branch_point or branch_states to select subset of cells")
    if (!is.null(branch_labels) & !is.null(branch_states)) {
        if (length(branch_labels) != length(branch_states)) 
            stop("length of branch_labels doesn't match with that of branch_states")
        branch_map <- setNames(branch_labels, as.character(branch_states))
    }
    if (cds@dim_reduce_type == "DDRTree") {
        pr_graph_cell_proj_mst <- minSpanningTree(cds)
    }
    else {
        pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
    }
    root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
    root_state <- pData(cds)[root_cell, ]$State
    pr_graph_root <- subset(pData(cds), State == root_state)
    if (cds@dim_reduce_type == "DDRTree") {
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root), 
            ]
    }
    else {
        root_cell_point_in_Y <- row.names(pr_graph_root)
    }
    root_cell <- names(which(igraph:::degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, 
        mode = "all") == 1, useNames = T))[1]
    paths_to_root <- list()
    if (is.null(branch_states) == FALSE) {
        for (leaf_state in branch_states) {
            curr_cell <- subset(pData(cds), State == leaf_state)
            if (cds@dim_reduce_type == "DDRTree") {
                closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell), 
                  ]
            }
            else {
                curr_cell_point_in_Y <- row.names(curr_cell)
            }
            curr_cell <- names(which(igraph:::degree(pr_graph_cell_proj_mst, 
                v = curr_cell_point_in_Y, mode = "all") == 1, 
                useNames = T))[1]
            path_to_ancestor <- igraph:::shortest_paths(pr_graph_cell_proj_mst, 
                curr_cell, root_cell)
            path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
            if (cds@dim_reduce_type == "DDRTree") {
                closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                ancestor_cells_for_branch <- row.names(closest_vertex)[which(igraph:::V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% 
                  path_to_ancestor)]
            }
            else if (cds@dim_reduce_type == "ICA") {
                ancestor_cells_for_branch <- path_to_ancestor
            }
            ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch, 
                colnames(cds))
            paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_branch
        }
    }
    else {
        if (cds@dim_reduce_type == "DDRTree") 
            pr_graph_cell_proj_mst <- minSpanningTree(cds)
        else pr_graph_cell_proj_mst <- cds@auxOrderingData$ICA$cell_ordering_tree
        mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
        branch_cell <- mst_branch_nodes[branch_point]
        mst_no_branch_point <- pr_graph_cell_proj_mst - igraph:::V(pr_graph_cell_proj_mst)[branch_cell]
        path_to_ancestor <- igraph:::shortest_paths(pr_graph_cell_proj_mst, 
            branch_cell, root_cell)
        path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
        for (backbone_nei in igraph:::V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name) {
            descendents <- igraph:::bfs(mst_no_branch_point, igraph:::V(mst_no_branch_point)[backbone_nei], 
                unreachable = FALSE)
            descendents <- descendents$order[!is.na(descendents$order)]
            descendents <- igraph:::V(mst_no_branch_point)[descendents]$name
            if (root_cell %in% descendents == FALSE) {
                path_to_root <- unique(c(path_to_ancestor, branch_cell, 
                  descendents))
                if (cds@dim_reduce_type == "DDRTree") {
                  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                  path_to_root <- row.names(closest_vertex)[which(igraph:::V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% 
                    path_to_root)]
                }
                else {
                  path_to_root <- path_to_root
                }
                closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
                path_to_root <- intersect(path_to_root, colnames(cds))
                paths_to_root[[backbone_nei]] <- path_to_root
            }
        }
    }
    all_cells_in_subset <- c()
    if (is.null(branch_labels) == FALSE) {
        if (length(branch_labels) != 2) 
            stop("Error: branch_labels must have exactly two entries")
        names(paths_to_root) <- branch_labels
    }
    for (path_to_ancestor in paths_to_root) {
        if (length(path_to_ancestor) == 0) {
            stop("Error: common ancestors between selected State values on path to root State")
        }
        all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
    }
    all_cells_in_subset <- unique(all_cells_in_subset)
    common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
    cds <- cds[, row.names(pData(cds[, all_cells_in_subset]))]
    Pseudotime <- pData(cds)$Pseudotime
    pData <- pData(cds)
    if (stretch) {
        max_pseudotime <- -1
        for (path_to_ancestor in paths_to_root) {
            max_pseudotime_on_path <- max(pData[path_to_ancestor, 
                ]$Pseudotime)
            if (max_pseudotime < max_pseudotime_on_path) {
                max_pseudotime <- max_pseudotime_on_path
            }
        }
        branch_pseudotime <- max(pData[common_ancestor_cells, 
            ]$Pseudotime)
        for (path_to_ancestor in paths_to_root) {
            max_pseudotime_on_path <- max(pData[path_to_ancestor, 
                ]$Pseudotime)
            path_scaling_factor <- (max_pseudotime - branch_pseudotime)/(max_pseudotime_on_path - 
                branch_pseudotime)
            if (is.finite(path_scaling_factor)) {
                branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
                pData[branch_cells, ]$Pseudotime <- ((pData[branch_cells, 
                  ]$Pseudotime - branch_pseudotime) * path_scaling_factor + 
                  branch_pseudotime)
            }
        }
        pData$Pseudotime <- 100 * pData$Pseudotime/max_pseudotime
    }
    pData$original_cell_id <- row.names(pData)
    pData$original_cell_id <- row.names(pData)
    if (length(paths_to_root) != 2) 
        stop("more than 2 branch states are used!")
    pData[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1]
    progenitor_pseudotime_order <- order(pData[common_ancestor_cells, 
        "Pseudotime"])
    if (progenitor_method == "duplicate") {
        ancestor_exprs <- exprs(cds)[, common_ancestor_cells]
        expr_blocks <- list()
        for (i in 1:length(paths_to_root)) {
            if (nrow(ancestor_exprs) == 1) 
                exprs_data <- t(as.matrix(ancestor_exprs))
            else exprs_data <- ancestor_exprs
            colnames(exprs_data) <- paste("duplicate", i, 1:length(common_ancestor_cells), 
                sep = "_")
            expr_lineage_data <- exprs(cds)[, setdiff(paths_to_root[[i]], 
                common_ancestor_cells)]
            exprs_data <- cbind(exprs_data, expr_lineage_data)
            expr_blocks[[i]] <- exprs_data
        }
        ancestor_pData_block <- pData[common_ancestor_cells, 
            ]
        pData_blocks <- list()
        weight_vec <- c()
        for (i in 1:length(paths_to_root)) {
            weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
            weight_vec_block <- rep(1, length(common_ancestor_cells))
            new_pData_block <- ancestor_pData_block
            row.names(new_pData_block) <- paste("duplicate", 
                i, 1:length(common_ancestor_cells), sep = "_")
            pData_lineage_cells <- pData[setdiff(paths_to_root[[i]], 
                common_ancestor_cells), ]
            weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
            weight_vec <- c(weight_vec, weight_vec_block)
            new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
            new_pData_block$Branch <- names(paths_to_root)[i]
            pData_blocks[[i]] <- new_pData_block
        }
        pData <- do.call(rbind, pData_blocks)
        exprs_data <- do.call(cbind, expr_blocks)
    }
    else if (progenitor_method == "sequential_split") {
        pData$Branch <- names(paths_to_root)[1]
        branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), 
            by = 2)]
        pData[common_ancestor_cells[branchA], "Branch"] <- names(paths_to_root)[1]
        branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), 
            by = 2)]
        pData[common_ancestor_cells[branchB], "Branch"] <- names(paths_to_root)[2]
        zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
        exprs_data <- cbind(exprs(cds), duplicate_root = exprs(cds)[, 
            zero_pseudotime_root_cell])
        pData <- rbind(pData, pData[zero_pseudotime_root_cell, 
            ])
        row.names(pData)[nrow(pData)] <- "duplicate_root"
        pData[nrow(pData), "Branch"] <- names(paths_to_root)[2]
        weight_vec <- rep(1, nrow(pData))
        for (i in 1:length(paths_to_root)) {
            path_to_ancestor <- paths_to_root[[i]]
            branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
            pData[branch_cells, ]$Branch <- names(paths_to_root)[i]
        }
    }
    pData$Branch <- as.factor(pData$Branch)
    pData$State <- factor(pData$State)
    Size_Factor <- pData$Size_Factor
    fData <- fData(cds)
    colnames(exprs_data) <- row.names(pData)
    cds_subset <- newCellDataSet(as.matrix(exprs_data), phenoData = new("AnnotatedDataFrame", 
        data = pData), featureData = new("AnnotatedDataFrame", 
        data = fData), expressionFamily = cds@expressionFamily, 
        lowerDetectionLimit = cds@lowerDetectionLimit)
    pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
    pData(cds_subset)$Size_Factor <- Size_Factor
    cds_subset@dispFitInfo <- cds@dispFitInfo
    return(cds_subset)
}