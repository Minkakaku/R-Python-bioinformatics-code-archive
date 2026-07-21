netAnalysis_signalingRole_network_mod<-function (object, signaling, slot.name = "netP", measure = c("outdeg",
    "indeg", "flowbet", "info"), measure.name = c("Sender", "Receiver",
    "Mediator", "Influencer"), color.use = NULL, color.heatmap = "BuGn",
    width = 6.5, height = 1.4, font.size = 8, font.size.title = 10,
    cluster.rows = FALSE, cluster.cols = FALSE)
{
    if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    centr <- slot(object, slot.name)$centr[signaling]
    for (i in 1:length(centr)) {
        centr0 <- centr[[i]]
        mat <- matrix(unlist(centr0), ncol = length(centr0),
            byrow = FALSE)
        mat <- t(mat)
        rownames(mat) <- names(centr0)
        colnames(mat) <- names(centr0$outdeg)
        if (!is.null(measure)) {
            mat <- mat[measure, , drop = FALSE]
            if (!is.null(measure.name)) {
                if (length(measure.name) != length(measure)) {
                  stop("The length of `measure.name` is not the same as that of `measure`! Please modify it! \n")
                }
                rownames(mat) <- measure.name
            }
        }
        mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
        colnames(mat)[15]<-'NK'
        colnames(mat)[18]<-'γδT'
        print(mat)
        if (is.null(color.use)) {
            color.use <- scPalette(length(colnames(mat)))
        }
        color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
            name = color.heatmap))))(100)
        df <- data.frame(group = colnames(mat))
        rownames(df) <- colnames(mat)
        cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
        print(cell.cols.assigned)
        col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),
            which = "column", show_legend = FALSE, show_annotation_name = FALSE,
            simple_anno_size = grid::unit(0.2, "cm"))
        ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white",
            name = "Importance", bottom_annotation = col_annotation,
            cluster_rows = cluster.rows, cluster_columns = cluster.cols,
            row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size),
            column_names_gp = gpar(fontsize = font.size), width = unit(width,
                "cm"), height = unit(height, "cm"), column_title = paste0(names(centr[i]),
                " signaling pathway network"), column_title_gp = gpar(fontsize = font.size.title),
            column_names_rot = 45, heatmap_legend_param = list(title = "Importance",
                title_gp = gpar(fontsize = 8, fontface = "plain"),
                title_position = "leftcenter-rot", border = NA,
                at = c(round(min(mat, na.rm = T), digits = 1),
                  round(max(mat, na.rm = T), digits = 1)), legend_height = unit(20,
                  "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2,
                  "mm")))
        draw(ht1)
    }
}

 netVisual_individual_mod<-function (object, signaling, signaling.name = NULL, pairLR.use = NULL,
    color.use = NULL, vertex.receiver = NULL, sources.use = NULL,
    targets.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = 1,
    vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex = 0.8,
    weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8,
    graphics.init = TRUE, layout = c("circle", "hierarchy", "chord",
        "spatial"), height = 5, thresh = 0.05, alpha.image = 0.15,
    point.size = 1.5, group = NULL, cell.order = NULL, small.gap = 1,
    big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE,
    legend.pos.x = 20, legend.pos.y = 20, nCol = NULL, ...)
{
    layout <- match.arg(layout)
    if (is.null(vertex.weight)) {
        vertex.weight <- as.numeric(table(object@idents))
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig,
        key = "pathway_name", matching.exact = T, pair.only = F)
    if (is.null(signaling.name)) {
        signaling.name <- signaling
    }
    net <- object@net
    pairLR.use.name <- dimnames(net$prob)[[3]]
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    if (!is.null(pairLR.use)) {
        if (is.data.frame(pairLR.use)) {
            pairLR.name <- intersect(pairLR.name, as.character(pairLR.use$interaction_name))
        }
        else {
            pairLR.name <- intersect(pairLR.name, as.character(pairLR.use))
        }
        if (length(pairLR.name) == 0) {
            stop("There is no significant communication for the input L-R pairs!")
        }
    }
    pairLR <- pairLR[pairLR.name, ]
    prob <- net$prob
    pval <- net$pval
    prob[pval > thresh] <- 0
    if (length(pairLR.name) > 1) {
        pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name],
            3, sum) != 0]
    }
    else {
        pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) !=
            0]
    }
    if (length(pairLR.name.use) == 0) {
        stop(paste0("There is no significant communication of ",
            signaling.name))
    }
    else {
        pairLR <- pairLR[pairLR.name.use, ]
    }
    nRow <- length(pairLR.name.use)
    prob <- prob[, , pairLR.name.use]
    pval <- pval[, , pairLR.name.use]
    if (is.null(nCol)) {
        nCol <- min(length(pairLR.name.use), 2)
    }
    if (length(dim(prob)) == 2) {
        prob <- replicate(1, prob, simplify = "array")
        pval <- replicate(1, pval, simplify = "array")
    }
    if (is.null(edge.weight.max)) {
        edge.weight.max = max(prob)
    }
    if (layout == "hierarchy") {
        if (graphics.init) {
            par(mfrow = c(nRow, 2), mar = c(5, 4, 4, 2) + 0.1)
        }
        for (i in 1:length(pairLR.name.use)) {
            signalName_i <- pairLR$interaction_name_2[i]
            prob.i <- prob[, , i]
            colnames(prob.i)[15]<-'NK'
						rownames(prob.i)[15]<-'NK'
						prob.i[prob.i < 0.01] <- 0
						prob.i[, c("Proliferative Tcells","Th2","Treg")] <- 0
						print(prob.i)
            netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver,
                sources.use = sources.use, targets.use = targets.use,
                remove.isolate = remove.isolate, top = top, color.use = color.use,
                vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max,
                vertex.size.max = vertex.size.max, weight.scale = weight.scale,
                edge.weight.max = edge.weight.max, edge.width.max = edge.width.max,
                title.name = signalName_i, ...)
            netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),
                vertex.receiver), sources.use = sources.use,
                targets.use = targets.use, remove.isolate = remove.isolate,
                top = top, color.use = color.use, vertex.weight = vertex.weight,
                vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max,
                weight.scale = weight.scale, edge.weight.max = edge.weight.max,
                edge.width.max = edge.width.max, title.name = signalName_i,
                ...)
        }
        gg <- recordPlot()
    }
    else if (layout == "circle") {
        if (graphics.init) {
            par(mfrow = c(ceiling(length(pairLR.name.use)/nCol),
                nCol), xpd = TRUE)
        }
        gg <- vector("list", length(pairLR.name.use))
        for (i in 1:length(pairLR.name.use)) {
            signalName_i <- pairLR$interaction_name_2[i]
            prob.i <- prob[, , i]
            gg[[i]] <- netVisual_circle(prob.i, sources.use = sources.use,
                targets.use = targets.use, remove.isolate = remove.isolate,
                top = top, color.use = color.use, vertex.weight = vertex.weight,
                vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max,
                weight.scale = weight.scale, edge.weight.max = edge.weight.max,
                edge.width.max = edge.width.max, title.name = signalName_i,
                vertex.label.cex = vertex.label.cex, ...)
        }
    }
    else if (layout == "spatial") {
        if (graphics.init) {
            par(mfrow = c(ceiling(length(pairLR.name.use)/nCol),
                nCol), xpd = TRUE)
        }
        coordinates <- object@images$coordinates
        labels <- object@idents
        meta.t <- object@meta
        meta.t$labels <- labels
        gg <- vector("list", length(pairLR.name.use))
        for (i in 1:length(pairLR.name.use)) {
            signalName_i <- pairLR$interaction_name_2[i]
            prob.i <- prob[, , i]
            gg[[i]] <- netVisual_spatial(prob.i, coordinates = coordinates,
                meta = meta.t, sample.use = sample.use, alpha.image = alpha.image,
                point.size = point.size, sources.use = sources.use,
                targets.use = targets.use, idents.use = idents.use,
                remove.isolate = remove.isolate, top = top, color.use = color.use,
                vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max,
                vertex.size.max = vertex.size.max, weight.scale = weight.scale,
                edge.weight.max = edge.weight.max, edge.width.max = edge.width.max,
                title.name = signalName_i, vertex.label.cex = vertex.label.cex,
                ...)
        }
    }
    else if (layout == "chord") {
        if (graphics.init) {
            par(mfrow = c(ceiling(length(pairLR.name.use)/nCol),
                nCol), xpd = TRUE)
        }
        gg <- vector("list", length(pairLR.name.use))
        for (i in 1:length(pairLR.name.use)) {
            title.name <- pairLR$interaction_name_2[i]
            net <- prob[, , i]
            gg[[i]] <- netVisual_chord_cell_internal(net, color.use = color.use,
                sources.use = sources.use, targets.use = targets.use,
                remove.isolate = remove.isolate, group = group,
                cell.order = cell.order, lab.cex = vertex.label.cex,
                small.gap = small.gap, big.gap = big.gap, scale = scale,
                reduce = reduce, title.name = title.name, show.legend = show.legend,
                legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
        }
    }
    return(gg)
}
